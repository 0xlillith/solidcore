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
//!
//! Deviations: owned generator instead of an owning pointer; explicit Kernel.
//! A failed canonical generator construction returns its error instead of
//! dereferencing the native NULL output. No recursion. This is a work in progress;
//! the remaining sphere functions and native comparisons belong to PR #79.

use crate::curve::curve::Curve;
use crate::curve::ellipse::Ellipse;
use crate::error::Result;
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::extent2d::Extent2d;
use crate::kernel::Kernel;
use crate::math::{Point3, Vec3, UNDEF_DOUBLE};
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::nurbs::b_spline_surface::{BSplineSurface, BSplineSurfaceForm};
use crate::surface::surf_of_revolution::SurfOfRevolution;
use crate::surface::surface::p2::SPHERE_TYPE;

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
