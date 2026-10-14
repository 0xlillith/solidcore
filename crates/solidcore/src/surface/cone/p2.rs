//! Cone offsets and recognition, part 2: ellipse intervals inside cone trims.
//!
//! Upstream: `source/SMLib/src/SmCone.cpp:5282`, `source/SMLib/inc/SmCone.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`ellipse_trim_intervals`] | `sm_GetEllipseTrimIntervals` |
//!
//! Deviations: returns intervals instead of an output array. A failed primitive
//! intersection/evaluation is returned where upstream ignores its status and may
//! read unset output. The six-transition limit is checked; there is no recursion.
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
    let mut top_hits = circle_plane_intersect(
        tolerance,
        cone.top_radius(),
        top,
        cone_axis,
        origin,
        ellipse_normal,
    )?;
    let mut bottom_hits = circle_plane_intersect(
        tolerance,
        cone.bot_radius(),
        bottom,
        cone_axis,
        origin,
        ellipse_normal,
    )?;
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
                    let location = cone
                        .step_inversion(kernel, *point, tolerance, None)?
                        .location;
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
        let location = cone
            .step_inversion(kernel, midpoint, tolerance, None)?
            .location;
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
