//! Cone offsets, NURB recognition and rebuilding, and ellipse trim intervals.
//!
//! Upstream: `source/SMLib/src/SmCone.cpp:5282`, `source/SMLib/inc/SmCone.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`ellipse_trim_intervals`] | `sm_GetEllipseTrimIntervals` |
//! | [`Cone::create_offset_surface`] | `SmCone::CreateOffsetSurface` |
//! | [`Cone::is_nurb_surface_cone`] | `SmCone::IsNurbSurfaceCone` |
//! | [`Cone::rebuild_step_from_nurb_parameters`] | `SmCone::RebuildSTEPFromNURBParameters` |
//!
//! Deviations: returns intervals instead of an output array. A failed primitive
//! intersection/evaluation is returned where upstream ignores its status and may
//! read unset output. A STEP classification actually written before an error is
//! retained because this helper reads only that classification. Circle-plane errors
//! retain the native initialized empty hit count. The six-transition limit is
//! checked; there is no recursion. Lazy NURB creation for offsets uses an owned
//! working copy instead of mutating the input cache. Native attribute notifications
//! have no Rust equivalent. Offset homogeneous conversion preserves the direct
//! divisions of `TO_EUCLID`, which differs from `N_CPtToPtEuclid`'s reciprocal multiply.
//! Cone STEP inversion deliberately uses the cone's own analytic domain, as its
//! native override ignores the domain argument passed by this helper.

use super::p1::Cone;
use crate::curve::ellipse::{Ellipse, EllipseOptions};
use crate::error::{Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::geom_utility::p1::{line_plane_intersect, shell_sort_double_array};
use crate::geom::geom_utility::p2::{circle_plane_intersect, CurveIntersection};
use crate::kernel::Kernel;
use crate::math::{max2, EFF_ZERO};
use crate::nurbs::b_spline_surface::LocationType;

const OP_TRIM: &str = "surface::cone::p2::ellipse_trim_intervals";

#[derive(Clone, Copy, PartialEq)]
enum Membership {
    Outside,
    Inside,
    Unknown,
}

// The native helper reads only eLoc and ignores STEPInversion's status.
// Preserve a classification actually written before a generator error, while
// keeping missing callees and failures before classification explicit.
// CirclePlaneIntersect sets its output count to zero before any failing
// operation. The native trim helper ignores the status and consumes that count.
fn trim_circle_hits(result: Result<CurveIntersection>) -> Result<CurveIntersection> {
    match result {
        Ok(hits) => Ok(hits),
        Err(error) if !error.is_unported() => Ok(CurveIntersection::Points(Vec::new())),
        Err(error) => Err(error),
    }
}

fn trim_point_location(
    kernel: &Kernel,
    cone: &Cone,
    point: crate::math::Point3,
    tolerance: f64,
) -> Result<LocationType> {
    let positive = cone.positive_point(kernel, point)?;
    let mut location = None;
    match cone.rev.step_inversion_tracking_location(
        kernel,
        &cone.rev.anal_uv_domain,
        positive.point,
        tolerance,
        None,
        &mut location,
    ) {
        Ok(inversion) => Ok(inversion.location),
        Err(error) if !error.is_unported() => location.ok_or(error),
        Err(error) => Err(error),
    }
}

fn push_angle(angles: &mut Vec<f64>, value: f64) -> Result<()> {
    if angles.len() == 6 {
        return Err(SolidError::invalid_input(
            OP_TRIM,
            "more than six trim transitions",
        ));
    }
    angles.push(value);
    Ok(())
}

/// Ellipse angular intervals (degrees) within the cone's NURB trim domain.
/// `ellipse` supplies the position and radii; `angular_interval` supplies the
/// requested trim. The native algorithm first intersects a full 360-degree ellipse.
///
/// # Errors
/// Primitive geometry failures and invalid transition data, rather than reading
/// an unset C++ output after an ignored failing status.
pub fn ellipse_trim_intervals(
    kernel: &Kernel,
    ellipse: &Ellipse,
    angular_interval: &Extent1d,
    cone: &Cone,
    uv_domain: &Extent2d,
    tolerance: f64,
) -> Result<Vec<Extent1d>> {
    let position = ellipse.position();
    let origin = position.origin();
    let x_axis = position.x_axis();
    let y_axis = position.y_axis();
    let ellipse_normal = x_axis.cross(y_axis);
    let full = Ellipse::new(
        kernel,
        origin,
        x_axis,
        y_axis,
        Extent1d::from_min_max(0.0, 360.0),
        ellipse.x_radius(),
        ellipse.y_radius(),
        &EllipseOptions {
            dim: 3,
            ..EllipseOptions::default()
        },
    )?;
    let height = cone.height(kernel)?;
    let cone_axis = cone.rev.position.z_axis();
    let bottom = cone.rev.position.origin() + cone.rev.step_uv_domain().min().y * cone_axis;
    let top = bottom + height * cone_axis;
    let domain =
        crate::surface::surface::p1::convert_domain_from_nurbs_to_step_with(uv_domain, |uv| {
            cone.convert_uv_from_nurbs_to_step(kernel, uv)
        })?;
    let corners = [
        cone.evaluate_step_point(domain.evaluate(0.0, 0.0))?,
        cone.evaluate_step_point(domain.evaluate(0.0, 1.0))?,
        cone.evaluate_step_point(domain.evaluate(1.0, 0.0))?,
        cone.evaluate_step_point(domain.evaluate(1.0, 1.0))?,
    ];
    let start_vector = corners[1] - corners[0];
    let end_vector = corners[3] - corners[2];
    let mut top_hits = trim_circle_hits(circle_plane_intersect(
        tolerance,
        cone.top_radius(),
        top,
        cone_axis,
        origin,
        ellipse_normal,
    ))?;
    let mut bottom_hits = trim_circle_hits(circle_plane_intersect(
        tolerance,
        cone.bot_radius(),
        bottom,
        cone_axis,
        origin,
        ellipse_normal,
    ))?;
    let start_parameter = line_plane_intersect(corners[0], start_vector, origin, ellipse_normal)?;
    let end_parameter = line_plane_intersect(corners[2], end_vector, origin, ellipse_normal)?;
    for hits in [&mut top_hits, &mut bottom_hits] {
        if let CurveIntersection::Points(points) = hits {
            if let [first, second] = points.as_slice() {
                if first.distance(*second) < tolerance {
                    points.truncate(1);
                }
            }
        }
    }
    let coincident = matches!(top_hits, CurveIntersection::Coincident)
        || matches!(bottom_hits, CurveIntersection::Coincident);
    let mut angles = Vec::new();
    // The native order is top-circle intersections, then bottom-circle intersections.
    for hits in [&top_hits, &bottom_hits] {
        if let CurveIntersection::Points(points) = hits {
            if points.len() == 2 {
                for point in points {
                    let location = trim_point_location(kernel, cone, *point, tolerance)?;
                    if location == LocationType::Interior {
                        push_angle(
                            &mut angles,
                            full.step_inversion(kernel, *point, false, 0.0).0,
                        )?;
                    }
                }
            }
        }
    }
    let scaled_zero =
        EFF_ZERO * (1.0 + max2(corners[0].max_dimension(), corners[3].max_dimension()));
    for (corner, vector, parameter) in [
        (corners[0], start_vector, start_parameter),
        (corners[2], end_vector, end_parameter),
    ] {
        if coincident || (parameter >= 0.0 - scaled_zero && parameter <= 1.0 + scaled_zero) {
            let point = corner + parameter * vector;
            push_angle(
                &mut angles,
                full.step_inversion(kernel, point, false, 0.0).0,
            )?;
        }
    }
    let count = angles.len();
    if count < 2 {
        return Ok(Vec::new());
    }
    shell_sort_double_array(&mut angles);
    let angular_zero = EFF_ZERO * 100.0 * 360.0;
    let mut membership = Vec::new();
    for (index, angle) in angles.iter().copied().enumerate() {
        let next = (index + 1) % count;
        let start = if next == 0 { angle - 360.0 } else { angle };
        let end = *angles.get(next).ok_or(SolidError::sm_err(OP_TRIM))?;
        if end - start < angular_zero {
            membership.push(Membership::Unknown);
            continue;
        }
        let midpoint = full.evaluate_step_point((end + start) / 2.0)?;
        let location = trim_point_location(kernel, cone, midpoint, tolerance)?;
        membership.push(if location == LocationType::Exterior {
            Membership::Outside
        } else {
            Membership::Inside
        });
    }
    let mut result = Vec::new();
    for (index, angle) in angles.iter().copied().enumerate() {
        let next = (index + 1) % count;
        let previous = (index + count - 1) % count;
        let state = *membership.get(index).ok_or(SolidError::sm_err(OP_TRIM))?;
        if state == Membership::Inside
            || (state == Membership::Unknown
                && membership.get(previous) == Some(&Membership::Outside)
                && membership.get(next) == Some(&Membership::Outside))
        {
            let start = if next == 0 { angle - 360.0 } else { angle };
            let end = *angles.get(next).ok_or(SolidError::sm_err(OP_TRIM))?;
            let mut interval = Extent1d::new();
            interval.add_value(start);
            interval.add_value(end);
            result.extend(angular_interval.intersect_periodic(&interval, 360.0, false)?);
        }
    }
    Ok(result)
}

// SmCone.cpp:4305 and :4695 share endpoint analysis, but only recognition
// checks the knot multiplicities and middle circle. Rebuild must not gain
// those extra filters, and uses the top frame's x-axis even when not inside-out.
use crate::geom::axis2_placement::Axis2Placement;
use crate::math::EFF_ZERO_SQRT;
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::nurbs::b_spline_surface::p2::IsoCurveOptions;
use crate::nurbs::b_spline_surface::{BSplineSurface, SurfParamType};

#[derive(Clone, Copy)]
struct ConeCircle {
    frame: Axis2Placement,
    radius: f64,
    start: f64,
    end: f64,
}

struct ConeDetection {
    top: ConeCircle,
    bottom: ConeCircle,
    height: f64,
    axis: crate::math::Vec3,
    inside_out: bool,
    swap_uv: bool,
    linear_interval: Extent1d,
}

fn circle_data(curve: &BSplineCurve, tolerance: f64) -> Result<Option<ConeCircle>> {
    Ok(crate::curve::curve::p5::is_arc_base(
        curve,
        |c, tol| c.is_degenerate(tol, None),
        5,
        tolerance,
    )?
    .map(|(frame, radius, start, end)| ConeCircle {
        frame,
        radius,
        start,
        end,
    }))
}

fn detect_cone(
    kernel: &Kernel,
    surface: &BSplineSurface,
    scale: f64,
    strict: bool,
) -> Result<Option<ConeDetection>> {
    let degree_u = surface.degree(SurfParamType::U);
    let degree_v = surface.degree(SurfParamType::V);
    let domain = surface.natural_uv_domain();
    let direction = if degree_u == 1 && surface.number_control_points(SurfParamType::U) == 2 {
        SurfParamType::V
    } else if degree_v == 1 && surface.number_control_points(SurfParamType::V) == 2 {
        SurfParamType::U
    } else {
        SurfParamType::Neither
    };
    if degree_u + degree_v != 3 || !surface.is_rational() || direction == SurfParamType::Neither {
        return Ok(None);
    }
    let linear_interval = if direction == SurfParamType::V {
        domain.u_interval()
    } else {
        domain.v_interval()
    };
    if strict {
        let (_, multiplicities) = surface.knots(direction, None)?;
        if multiplicities.first() != Some(&3)
            || multiplicities.last() != Some(&3)
            || multiplicities
                .iter()
                .skip(1)
                .take(multiplicities.len().saturating_sub(2))
                .any(|m| *m != 2)
        {
            return Ok(None);
        }
    }
    let min = surface.evaluate_point(domain.min())?;
    let max = surface.evaluate_point(domain.max())?;
    let angle_tolerance = scale * EFF_ZERO_SQRT;
    let zero = scale * 100.0 * EFF_ZERO * (1.0 + min.max_dimension() + max.max_dimension());
    let (swap_uv, sweep, bottom_parameter, top_parameter) = if degree_u == 1 {
        (true, SurfParamType::U, domain.min().x, domain.max().x)
    } else {
        (false, SurfParamType::V, domain.min().y, domain.max().y)
    };
    let top_curve = surface
        .create_iso_parametric_curve(kernel, sweep, top_parameter, IsoCurveOptions::default())?
        .curve;
    let bottom_curve = surface
        .create_iso_parametric_curve(kernel, sweep, bottom_parameter, IsoCurveOptions::default())?
        .curve;
    let top_degenerate = top_curve.is_degenerate(zero, None)?;
    let bottom_degenerate = bottom_curve.is_degenerate(zero, None)?;
    let mut top = if top_degenerate {
        None
    } else {
        circle_data(&top_curve, zero)?
    };
    let mut bottom = if bottom_degenerate {
        None
    } else {
        circle_data(&bottom_curve, zero)?
    };
    if let Some(circle) = top.filter(|_| bottom_degenerate) {
        let center = bottom_curve.evaluate_point(bottom_curve.natural_interval().mid())?;
        bottom = Some(ConeCircle {
            frame: Axis2Placement::new(center, circle.frame.x_axis(), circle.frame.y_axis()),
            radius: 0.0,
            ..circle
        });
    } else if let Some(circle) = bottom.filter(|_| top_degenerate) {
        let center = top_curve.evaluate_point(top_curve.natural_interval().mid())?;
        top = Some(ConeCircle {
            frame: Axis2Placement::new(center, circle.frame.x_axis(), circle.frame.y_axis()),
            radius: 0.0,
            ..circle
        });
    }
    let (Some(mut top), Some(bottom)) = (top, bottom) else {
        return Ok(None);
    };
    let twist = top.frame.x_axis().angle_between(bottom.frame.x_axis())?;
    if (top.start - bottom.start).abs() > angle_tolerance
        || (top.end - bottom.end).abs() > angle_tolerance
        || twist.abs() > angle_tolerance
    {
        return Ok(None);
    }
    let axis = top.frame.origin() - bottom.frame.origin();
    let height = axis.length();
    if height < EFF_ZERO_SQRT {
        return Ok(None);
    }
    let axis = axis.unitize()?;
    if !axis.is_parallel_to(bottom.frame.z_axis(), angle_tolerance) {
        return Ok(None);
    }
    if strict {
        let middle_curve = surface
            .create_iso_parametric_curve(
                kernel,
                sweep,
                (bottom_parameter + top_parameter) / 2.0,
                IsoCurveOptions::default(),
            )?
            .curve;
        let Some(middle) = circle_data(&middle_curve, zero)? else {
            return Ok(None);
        };
        if ((bottom.radius + top.radius) / 2.0 - middle.radius).abs() > zero
            || !middle.frame.z_axis().is_parallel_to(axis, angle_tolerance)
        {
            return Ok(None);
        }
        let (distance, _) = crate::geom::geom_utility::p1::line_point_distance(
            top.frame.origin(),
            axis,
            middle.frame.origin(),
        );
        if distance > zero {
            return Ok(None);
        }
    }
    let inside_out = axis.dot(bottom.frame.z_axis()) < 0.0;
    if (bottom.radius - top.radius).abs() < EFF_ZERO * (1.0 + bottom.radius) {
        top.radius = bottom.radius;
    }
    Ok(Some(ConeDetection {
        top,
        bottom,
        height,
        axis,
        inside_out,
        swap_uv,
        linear_interval,
    }))
}

impl Cone {
    /// Mutable inherited NURB base. Call `rebuild_step_from_nurb_parameters`
    /// after geometry edits to update the analytic and generator data.
    pub fn base_mut(&mut self) -> &mut BSplineSurface {
        &mut self.rev.base
    }

    /// Recognize a rational line/circle NURB as a cone (`IsNurbSurfaceCone`).
    /// The analytic object owns its regenerated NURB, with the input domain.
    ///
    /// # Errors
    /// Invalid NURB data, iso extraction/evaluation or reconstruction errors.
    pub fn is_nurb_surface_cone(
        kernel: &Kernel,
        surface: &BSplineSurface,
        tolerance_scale: f64,
    ) -> Result<Option<Self>> {
        let Some(detection) = detect_cone(kernel, surface, tolerance_scale, true)? else {
            return Ok(None);
        };
        let d = detection;
        let x_axis = if d.inside_out {
            d.top.frame.x_axis()
        } else {
            d.bottom.frame.x_axis()
        };
        let z_axis = if d.inside_out { -d.axis } else { d.axis };
        let mut cone = Self::new_bounded(
            if d.inside_out {
                d.top.frame.origin()
            } else {
                d.bottom.frame.origin()
            },
            x_axis,
            z_axis.cross(x_axis),
            if d.inside_out {
                d.top.radius
            } else {
                d.bottom.radius
            },
            if d.inside_out {
                d.bottom.radius
            } else {
                d.top.radius
            },
            d.top.start,
            d.top.end,
            d.height,
            d.swap_uv,
            d.inside_out,
            false,
        );
        cone.make_nurb(kernel)?;
        crate::error::ignore_status(cone.reparameterize(kernel, &surface.natural_uv_domain()))?;
        // Notify(SM_NO_COPY) copies attributes; the geometric base does not model them.
        Ok(Some(cone))
    }

    /// Reparameterize the inherited NURB with this cone's virtual iso curves.
    ///
    /// # Errors
    /// NURB or generator reparameterization and polar-converter errors.
    pub fn reparameterize(&mut self, kernel: &Kernel, domain: &Extent2d) -> Result<()> {
        if self.rev.base.gw_nurb().is_none() {
            self.make_nurb(kernel)?;
        }
        let inside_out = self.inside_out;
        self.rev
            .reparameterize_using(kernel, domain, &|rev, k, p, t| {
                Ok(super::p1::cone_iso_curve(rev, inside_out, k, p, t, None, false)?.curve)
            })
    }
}

impl Cone {
    /// Rebuild analytic and generator data from the current NURB
    /// (`RebuildSTEPFromNURBParameters`). Unlike recognition, this does not add
    /// a midpoint-circle or knot-multiplicity test.
    ///
    /// # Errors
    /// A NURB that no longer represents the required line/circle surface,
    /// extraction/evaluation, generator or polar-converter errors.
    pub fn rebuild_step_from_nurb_parameters(&mut self, kernel: &Kernel) -> Result<()> {
        use crate::curve::curve::Curve;
        use crate::curve::line::{BoundedLineOptions, Line};
        // GetDegree/GetNaturalUVDomain reach the native virtual NURB factory.
        if self.rev.base.gw_nurb().is_none() {
            self.make_nurb(kernel)?;
        }
        let Some(d) = detect_cone(kernel, &self.rev.base, 1000.0, false)? else {
            return Err(SolidError::sm_err(
                "surface::cone::p2::Cone::rebuild_step_from_nurb_parameters",
            ));
        };
        let radius_difference = if d.inside_out {
            d.bottom.radius - d.top.radius
        } else {
            d.top.radius - d.bottom.radius
        };
        // Upstream leaves the previous semi-angle untouched for equal radii.
        if radius_difference.abs() > EFF_ZERO {
            self.semi_angle_deg = radius_difference.atan2(d.height) * 180.0 / std::f64::consts::PI;
        }
        self.base_radius = if d.inside_out {
            d.top.radius
        } else {
            d.bottom.radius
        };
        self.inside_out = d.inside_out;
        self.rev.planar_generator = true;
        // In the native rebuild this is a frame copy, not a new frame made
        // from the independently computed cone axis and top x-axis.
        self.rev.position = if d.inside_out {
            d.top.frame
        } else {
            d.bottom.frame
        };
        self.rev.swap_uv = d.swap_uv;
        self.rev.anal_uv_domain = Extent2d::from_bounds(d.top.start, 0.0, d.top.end, d.height);
        let point = self.rev.position.origin() + self.base_radius * self.rev.position.x_axis();
        let vector = d.height * self.rev.position.z_axis()
            + (self.radius(d.height) - self.base_radius) * self.rev.position.x_axis();
        let length = vector.length();
        let mut line = Line::bounded(
            point,
            vector,
            &Extent1d::from_min_max(0.0, d.height),
            &BoundedLineOptions {
                scale: length / d.height,
                dim: 3,
                nurb: None,
                inside_out: d.inside_out,
            },
        );
        crate::error::ignore_status(line.edit_parameterization(&d.linear_interval))?;
        // The native order rebuilds the polar converter BEFORE replacing the
        // old generator. Keep this order, including the virtual cone iso call.
        let inside_out = self.inside_out;
        crate::error::ignore_status(self.rev.rebuild_polar_converter(kernel, &|rev, k, p, t| {
            Ok(super::p1::cone_iso_curve(rev, inside_out, k, p, t, None, false)?.curve)
        }))?;
        self.rev.gen_curve = Some(Curve::from(line));
        Ok(())
    }
}

// SmNurbsCrv.h:51 TO_EUCLID divides each coordinate directly by w.
// CPoint::to_euclid implements the DIFFERENT N_CPtToPtEuclid operation (multiply
// by reciprocal). Around a large cone apex that changes the computed height
// and radius enough to alter the following offset branch and movement.
fn offset_euclid(point: crate::nurbs::math::CPoint) -> crate::math::Point3 {
    use crate::nurbs::math::{NOW, NOZ};
    let z = if point.z == NOZ { 0.0 } else { point.z };
    if point.w == NOW {
        crate::math::Point3::new(point.x, point.y, z)
    } else {
        crate::math::Point3::new(
            point.x / point.w,
            point.y / point.w,
            if point.z == NOZ {
                0.0
            } else {
                point.z / point.w
            },
        )
    }
}

impl Cone {
    /// Analytic cone offset (`CreateOffsetSurface`), including apex truncation.
    /// `None` is the upstream successful empty result when both radii shrink
    /// below the approximation tolerance.
    ///
    /// # Errors
    /// Invalid geometry, canonical construction, missing generator, NURB access
    /// or unported callees. Ordinary statuses ignored by upstream stay ignored.
    pub fn create_offset_surface(
        &self,
        kernel: &Kernel,
        distance: f64,
        tolerance: f64,
    ) -> Result<Option<Self>> {
        let tolerance = crate::geom::tol::approx_tol3d(kernel, tolerance);
        // Native GetHeight precedes GetNaturalUVDomain's lazy MakeNurb.
        // Build a missing NURB on an owned working copy: this Rust API keeps
        // the input immutable while preserving the generated output geometry.
        let height = self.height(kernel)?;
        if self.rev.base.gw_nurb().is_none() {
            let mut ready = self.clone();
            ready.make_nurb(kernel)?;
            if ready.rev.base.gw_nurb().is_none() {
                return Err(SolidError::null_pointer(
                    "surface::cone::p2::Cone::create_offset_surface",
                ));
            }
            return ready.offset_with_height(kernel, distance, tolerance, height);
        }
        self.offset_with_height(kernel, distance, tolerance, height)
    }

    fn offset_with_height(
        &self,
        kernel: &Kernel,
        distance: f64,
        tolerance: f64,
        height: f64,
    ) -> Result<Option<Self>> {
        use crate::error::ignore_status;
        let sign = 1.0
            * (if self.rev.swap_uv { -1.0 } else { 1.0 })
            * (if self.inside_out { -1.0 } else { 1.0 });
        let origin = self.rev.position.origin();
        let x_axis = self.rev.position.x_axis();
        let y_axis = self.rev.position.y_axis();
        let z_axis = self.rev.position.z_axis();
        let domain = self.rev.base.natural_uv_domain();
        let mut sin = (self.semi_angle_deg * std::f64::consts::PI / 180.0).sin();
        let mut cos = (self.semi_angle_deg * std::f64::consts::PI / 180.0).cos();
        let normal = (cos * x_axis - sin * z_axis).unitize()?;
        let offset = sign * distance * normal;
        let offset_z = -sin * sign * distance;
        let offset_radius = cos * sign * distance;
        let old_bottom = self.bot_radius();
        let old_top = self.top_radius();
        let mut new_bottom = old_bottom + offset_radius;
        let mut new_top = old_top + offset_radius;
        if new_bottom < tolerance && new_top < tolerance {
            return Ok(None);
        }
        if new_bottom < tolerance || new_top < tolerance {
            let (new_height, bottom_origin) = if new_bottom < tolerance {
                let tangent = height / old_top;
                let new_height = new_top * tangent;
                new_bottom = 0.0;
                (
                    new_height,
                    origin + (height + offset_z - new_height) * z_axis,
                )
            } else {
                let tangent = height / old_bottom;
                let new_height = new_bottom * tangent;
                new_top = 0.0;
                (new_height, origin + offset_z * z_axis)
            };
            let placement = Axis2Placement::new(bottom_origin, x_axis, y_axis);
            let mut cone = Self::create_canonical(
                kernel,
                &placement,
                new_bottom,
                new_top,
                new_height,
                self.rev.swap_uv,
                self.inside_out,
            )?;
            // Preserve the native formula: it multiplies the input MAXIMUM,
            // rather than the input span, even for a nonzero domain minimum.
            let target = if self.rev.swap_uv {
                Extent2d::from_bounds(
                    domain.min().x,
                    domain.min().y,
                    domain.min().x + new_height / height * domain.max().x,
                    domain.max().y,
                )
            } else {
                Extent2d::from_bounds(
                    domain.min().x,
                    domain.min().y,
                    domain.max().x,
                    domain.min().y + new_height / height * domain.max().y,
                )
            };
            ignore_status(cone.rev.base.reparameterize(&target))?;
            let generator_interval = if self.rev.swap_uv {
                target.u_interval()
            } else {
                target.v_interval()
            };
            let generator = cone.rev.gen_curve.as_mut().ok_or(SolidError::null_pointer(
                "surface::cone::p2::Cone::create_offset_surface",
            ))?;
            ignore_status(generator.edit_parameterization(kernel, &generator_interval))?;
            let inside_out = cone.inside_out;
            ignore_status(cone.rev.rebuild_polar_converter(kernel, &|rev, k, p, t| {
                Ok(super::p1::cone_iso_curve(rev, inside_out, k, p, t, None, false)?.curve)
            }))?;
            return Ok(Some(cone));
        }
        let mut cone = self.clone();
        cone.base_radius = self.base_radius + offset_radius;
        ignore_status(
            cone.rev
                .position
                .set_canonical(origin + offset_z * z_axis, x_axis, y_axis),
        )?;
        if let Some(nurb) = cone.rev.base.nurb.as_mut() {
            let (outer, inner) = if self.rev.swap_uv {
                (nurb.m(), nurb.n())
            } else {
                (nurb.n(), nurb.m())
            };
            for i in 0..=outer {
                let (u, v) = if self.rev.swap_uv { (0, i) } else { (i, 0) };
                let first = offset_euclid(nurb.cpt(u, v)?);
                let z_distance = (first.x - origin.x) * z_axis.x
                    + (first.y - origin.y) * z_axis.y
                    + (first.z - origin.z) * z_axis.z;
                let center = origin + z_distance * z_axis;
                let radial = first - center;
                let radius = self.radius(z_distance);
                let movement = if crate::geom::tol::is_zero(radius) {
                    let (u, v) = if self.rev.swap_uv { (1, i) } else { (i, 1) };
                    let next = offset_euclid(nurb.cpt(u, v)?);
                    let next_vector = next - origin;
                    let next_height = next_vector.dot(z_axis);
                    let radial = next_vector - next_height * z_axis;
                    let actual_radius = radial.length();
                    let height_radius = self.radius(next_height);
                    let unit = radial / actual_radius;
                    cos = unit.dot(x_axis);
                    sin = unit.dot(y_axis);
                    actual_radius / height_radius * offset_radius * (cos * x_axis + sin * y_axis)
                        + offset_z * z_axis
                } else {
                    let scale = (radius + offset_radius) / radius;
                    radial * (scale - 1.0) + offset_z * z_axis
                };
                for j in 0..=inner {
                    let (u, v) = if self.rev.swap_uv { (j, i) } else { (i, j) };
                    let euclidean = offset_euclid(nurb.cpt(u, v)?);
                    let point = nurb.cpt_mut(u, v)?;
                    // The native code multiplies these coordinates by w.
                    point.x = (euclidean.x + movement.x) * point.w;
                    point.y = (euclidean.y + movement.y) * point.w;
                    point.z = (euclidean.z + movement.z) * point.w;
                }
            }
        }
        let mut translation = Axis2Placement::IDENTITY;
        translation.translate(offset);
        let generator = cone.rev.gen_curve.as_mut().ok_or(SolidError::null_pointer(
            "surface::cone::p2::Cone::create_offset_surface",
        ))?;
        ignore_status(generator.transform(kernel, &translation, None))?;
        Ok(Some(cone))
    }
}
