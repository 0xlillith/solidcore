//! Sphere geometry represented by its actual circular surface of revolution.
//!
//! Upstream: `source/SMLib/src/SmSphere.cpp`, `source/SMLib/inc/SmSphere.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`Sphere`], [`Sphere::new`], [`Sphere::empty`] | `SmSphere`, constructors |
//! | `Clone`, implicit `Drop` | copy constructor, destructor |
//! | [`Sphere::radius`], [`Sphere::max_analytic_domain`] | `GetRadius`, `GetMaxAnalyticDomain` |
//! | [`Sphere::is_kind_of`] | `IsKindOf` |
//! | [`Sphere::create_canonical`] | `CreateCanonical` |
//! | [`Sphere::create_iso_parametric_curve`] | `CreateIsoParametricCurve` |
//! | [`Sphere::evaluate_step`], [`Sphere::evaluate_step_point`] | `EvaluateSTEP`, `EvaluateSTEPPoint` |
//! | [`Sphere::make_nurb`] | inherited `MakeNurb` with sphere virtual iso calls |
//!
//! Deviations: owned generator instead of an owning pointer; explicit Kernel.
//! A failed canonical generator construction returns its error instead of
//! dereferencing the native NULL output. Unwritten derivative slots are zero. No recursion. This is a work in progress;
//! the remaining sphere functions and native comparisons belong to PR #79.

use crate::curve::curve::Curve;
use crate::curve::ellipse::{Ellipse, EllipseOptions};
use crate::error::{Result, SolidError};
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::vector2d::Point2;
use crate::kernel::Kernel;
use crate::math::{max2, min2, Point3, Vec3, EFF_ZERO, PI, UNDEF_DOUBLE};
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::nurbs::b_spline_surface::{
    BSplineSurface, BSplineSurfaceForm, IsoCurveOptions, SurfParamType,
};
use crate::surface::surf_of_revolution::{RevolutionIsoCurve, SurfOfRevolution};
use crate::surface::surface::p2::{scale_step_angular_derivatives, SPHERE_TYPE};

/// Optional native sphere constructor arguments, with their upstream defaults.
#[derive(Debug, Clone, Default)]
pub struct SphereOptions {
    /// Exchange the NURB sweep and generator parameter directions.
    pub swap_uv: bool,
    /// Reverse the circular generator's parameterization.
    pub inside_out: bool,
    /// An owned NURB generator; absent means an actual analytic circle.
    pub circle_nurb: Option<BSplineCurve>,
}

/// An analytic sphere (`SmSphere`) with its revolution and NURB state.
#[derive(Debug, Clone)]
pub struct Sphere {
    pub(crate) rev: SurfOfRevolution,
    pub(crate) radius: f64,
    pub(crate) inside_out: bool,
}

impl Sphere {
    /// Empty sphere for database input, preserving native unset radius.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            rev: SurfOfRevolution::empty(),
            radius: UNDEF_DOUBLE,
            inside_out: false,
        }
    }

    /// Construct over the angular domain in degrees, with a lazy surface NURB.
    ///
    /// # Errors
    /// Invalid generator geometry or a failed circular generator construction.
    pub fn new(
        kernel: &Kernel,
        origin: Point3,
        x_axis: Vec3,
        y_axis: Vec3,
        domain: Extent2d,
        radius: f64,
        options: SphereOptions,
    ) -> Result<Self> {
        let rev = if let Some(generator) = options.circle_nurb {
            SurfOfRevolution::new(
                kernel,
                Curve::from(generator),
                origin,
                x_axis,
                y_axis,
                domain,
                options.swap_uv,
            )?
        } else {
            // Native base construction with NULL skips the generator plane test.
            let mut rev = SurfOfRevolution::empty();
            rev.planar_generator = true;
            rev.position = Axis2Placement::new(origin, x_axis, y_axis);
            rev.anal_uv_domain = domain;
            rev.swap_uv = options.swap_uv;
            rev.base
                .set_b_spline_surface_form(BSplineSurfaceForm::SurfOfRevolution);
            let generator = Ellipse::create_canonical_circle(
                kernel,
                &Axis2Placement::new(origin, x_axis, x_axis.cross(y_axis)),
                radius,
                Some(&domain.v_interval()),
                Some(options.inside_out),
            )?;
            rev.gen_curve = Some(Curve::from(generator));
            rev
        };
        Ok(Self {
            rev,
            radius,
            inside_out: options.inside_out,
        })
    }

    /// Stored sphere radius (`GetRadius`).
    #[must_use]
    pub const fn radius(&self) -> f64 {
        self.radius
    }

    /// Whether the generating circle runs inside out (`GetInsideOut`).
    #[must_use]
    pub const fn inside_out(&self) -> bool {
        self.inside_out
    }

    /// The inherited geometric state, without changing the analytic type.
    #[must_use]
    pub const fn rev(&self) -> &SurfOfRevolution {
        &self.rev
    }

    /// The inherited NURB surface state.
    #[must_use]
    pub const fn base(&self) -> &BSplineSurface {
        &self.rev.base
    }

    /// Maximum analytic domain: longitude 0..360 and latitude -90..90 degrees.
    #[must_use]
    pub const fn max_analytic_domain() -> Extent2d {
        Extent2d::from_bounds(0.0, -90.0, 360.0, 90.0)
    }

    /// Native runtime type and base-class membership.
    ///
    /// # Errors
    /// Those of the inherited type query.
    pub fn is_kind_of(&self, type_code: i64) -> Result<bool> {
        if type_code == SPHERE_TYPE {
            Ok(true)
        } else {
            self.rev.is_kind_of(type_code)
        }
    }
}

fn sphere_iso_curve(
    rev: &SurfOfRevolution,
    radius: f64,
    inside_out: bool,
    kernel: &Kernel,
    param: SurfParamType,
    parameter: f64,
    domain: Option<&Extent2d>,
    uv_curve: bool,
) -> Result<RevolutionIsoCurve> {
    let iso = rev.base.create_iso_parametric_curve(
        kernel,
        param,
        parameter,
        IsoCurveOptions {
            domain: domain.copied(),
            uv_iso_curve: uv_curve,
        },
    )?;
    if iso.curve.is_degenerate(EFF_ZERO, None)? {
        return Ok(RevolutionIsoCurve {
            curve: Curve::from(iso.curve),
            max_gap_3d: iso.max_gap_3d,
            uv_iso_curve: iso.uv_iso_curve,
        });
    }
    let step_param = if rev.swap_uv {
        if param == SurfParamType::U {
            SurfParamType::V
        } else {
            SurfParamType::U
        }
    } else {
        param
    };
    let origin = rev.position.origin();
    let z = rev.position.z_axis();
    let interval = iso.curve.natural_interval();
    let min_point = iso.curve.evaluate_point(interval.min())?;
    let max_point = iso.curve.evaluate_point(interval.max())?;
    let min_vec = min_point - origin;
    let max_vec = max_point - origin;
    let (center, x, y, circle_radius, circle_inside, analytic) = if step_param == SurfParamType::U {
        let mid_point = iso.curve.evaluate_point(interval.evaluate(0.5))?;
        let mid_vec = mid_point - origin;
        (
            origin,
            z.cross(mid_vec.cross(z)).unitize_or_z(),
            z,
            radius,
            inside_out,
            rev.anal_uv_domain.v_interval(),
        )
    } else {
        let center = origin + z * min_vec.dot(z);
        (
            center,
            rev.position.x_axis(),
            rev.position.y_axis(),
            (min_point - center).length(),
            false,
            rev.anal_uv_domain.u_interval(),
        )
    };
    let circle_domain = if domain.is_none() {
        analytic
    } else {
        let circle_z = x.cross(y);
        let min_angle = analytic
            .contains_periodic_value(
                circle_z.ccw_angle_between(x, min_vec)? * 180.0 / PI,
                360.0,
                0.0,
            )
            .periodic_value;
        let max_angle = analytic
            .contains_periodic_value(
                circle_z.ccw_angle_between(x, max_vec)? * 180.0 / PI,
                360.0,
                0.0,
            )
            .periodic_value;
        let period = 360.0;
        let tol = EFF_ZERO * 100.0 * (1.0 + period);
        if analytic.is_closed(period)
            && analytic
                .is_value_on_periodic_boundary(min_angle, period, tol)
                .0
            && analytic
                .is_value_on_periodic_boundary(max_angle, period, tol)
                .0
        {
            analytic
        } else {
            Extent1d::from_min_max(min2(min_angle, max_angle), max2(min_angle, max_angle))
        }
    };
    let curve = Ellipse::new_circle(
        kernel,
        center,
        x,
        y,
        circle_domain,
        circle_radius,
        &EllipseOptions {
            dim: 3,
            nurb: Some(&iso.curve),
            inside_out: circle_inside,
        },
    )?;
    Ok(RevolutionIsoCurve {
        curve: Curve::from(curve),
        max_gap_3d: iso.max_gap_3d,
        uv_iso_curve: iso.uv_iso_curve,
    })
}

impl Sphere {
    /// Build the inherited NURB using the sphere's actual circular iso curves.
    ///
    /// # Errors
    /// Generator, NURB construction or polar conversion failure.
    pub fn make_nurb(&mut self, kernel: &Kernel) -> Result<()> {
        let radius = self.radius;
        let inside_out = self.inside_out;
        self.rev
            .make_nurb_with_sweep_params_using(kernel, None, &|rev, k, p, t| {
                Ok(sphere_iso_curve(rev, radius, inside_out, k, p, t, None, false)?.curve)
            })
    }

    /// Construct a full sphere; a negative radius selects inside-out orientation.
    ///
    /// # Errors
    /// Degenerate radius or failed NURB construction (`SM_ERR`).
    pub fn create_canonical(
        kernel: &Kernel,
        position: &Axis2Placement,
        radius: f64,
    ) -> Result<Self> {
        let inside_out = radius < 0.0;
        let abs_radius = radius.abs();
        let operation = "surface::sphere::Sphere::create_canonical";
        if abs_radius.is_nan() || abs_radius <= EFF_ZERO {
            return Err(SolidError::sm_err(operation));
        }
        let mut sphere = Self::new(
            kernel,
            position.origin(),
            position.x_axis(),
            position.y_axis(),
            Self::max_analytic_domain(),
            abs_radius,
            SphereOptions {
                inside_out,
                ..SphereOptions::default()
            },
        )?;
        match sphere.make_nurb(kernel) {
            Err(error) if error.is_unported() => return Err(error),
            Err(_) => return Err(SolidError::sm_err(operation)),
            Ok(()) => {}
        }
        if sphere.rev.base.gw_nurb().is_none() {
            return Err(SolidError::sm_err(operation));
        }
        Ok(sphere)
    }

    /// Extract a real circular iso curve, retaining degenerate B-spline results.
    ///
    /// # Errors
    /// Base iso extraction, evaluation, angular conversion or circle construction.
    pub fn create_iso_parametric_curve(
        &self,
        kernel: &Kernel,
        param: SurfParamType,
        parameter: f64,
        domain: Option<&Extent2d>,
        uv_curve: bool,
    ) -> Result<RevolutionIsoCurve> {
        sphere_iso_curve(
            &self.rev,
            self.radius,
            self.inside_out,
            kernel,
            param,
            parameter,
            domain,
            uv_curve,
        )
    }
}

impl Sphere {
    /// Evaluate the native triangular STEP derivative matrix through order three.
    /// Unwritten entries are zero, and angular derivatives are per degree.
    ///
    /// # Errors
    /// Orders above three (`SM_ERR`), or an invalid derivative layout.
    pub fn evaluate_step(&self, mut uv: Point2, hu: usize, hv: usize) -> Result<Vec<Vec3>> {
        let operation = "surface::sphere::Sphere::evaluate_step";
        if hu > 3 || hv > 3 {
            return Err(SolidError::sm_err(operation));
        }
        let origin = self.rev.position.origin();
        let x = self.rev.position.x_axis();
        let y = self.rev.position.y_axis();
        let z = self.rev.position.z_axis();
        let radius = self.radius;
        let min_highest = hu.min(hv);
        if !self.rev.is_out_of_bounds_enabled() && !self.rev.anal_uv_domain.contains_point(uv, 0.0)
        {
            uv = self.rev.anal_uv_domain.clamp_point(uv);
        }
        let cos_u = (uv.x * PI / 180.0).cos();
        let sin_u = (uv.x * PI / 180.0).sin();
        let cos_v = (uv.y * PI / 180.0).cos();
        let sin_v = (uv.y * PI / 180.0).sin();
        // Orders were bounded above: at most sixteen entries.
        let mut values = vec![Vec3::ZERO; (hu + 1) * (hv + 1)];
        let mut put = |index: usize, value: Vec3| -> Result<()> {
            let target = values
                .get_mut(index)
                .ok_or_else(|| SolidError::sm_err(operation))?;
            *target = value;
            Ok(())
        };
        put(
            0,
            origin + radius * (cos_v * (cos_u * x + sin_u * y) + sin_v * z),
        )?;
        if hu >= 1 {
            put(hv + 1, radius * cos_v * (-sin_u * x + cos_u * y))?;
        }
        if hv >= 1 {
            put(1, radius * (-sin_v * (cos_u * x + sin_u * y) + cos_v * z))?;
        }
        if hu >= 1 && hv >= 1 && min_highest >= 2 {
            put(hv + 2, radius * -sin_v * (-sin_u * x + cos_u * y))?;
        }
        if hu >= 1 && hv >= 2 && min_highest >= 3 {
            put(hv + 3, radius * -cos_v * (-sin_u * x + cos_u * y))?;
        }
        if hu >= 2 {
            put(2 * hv + 2, radius * cos_v * (-cos_u * x - sin_u * y))?;
        }
        if hv >= 2 {
            put(2, radius * (-cos_v * (cos_u * x + sin_u * y) - sin_v * z))?;
        }
        if hu >= 2 && hv >= 1 && min_highest >= 3 {
            put(2 * hv + 3, radius * -sin_v * (-cos_u * x - sin_u * y))?;
        }
        if hu >= 3 {
            put(3 * hv + 3, radius * cos_v * (sin_u * x - cos_u * y))?;
        }
        if hv >= 3 {
            put(3, radius * (sin_v * (cos_u * x + sin_u * y) - cos_v * z))?;
        }
        scale_step_angular_derivatives(&mut values, hu, hv, true, true)?;
        Ok(values)
    }

    /// Evaluate one STEP point, using the same clamping as the derivative call.
    ///
    /// # Errors
    /// Those of STEP evaluation.
    pub fn evaluate_step_point(&self, uv: Point2) -> Result<Point3> {
        self.evaluate_step(uv, 0, 0)?
            .first()
            .copied()
            .ok_or_else(|| SolidError::sm_err("surface::sphere::Sphere::evaluate_step_point"))
    }
}
