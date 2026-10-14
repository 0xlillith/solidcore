//! Differential coverage of real cone ellipse trims; setup failures fail cases.
use crate::{Case, HarnessError, Output, UnitDiff};
use serde_json::Value;
use solidcore::curve::ellipse::{Ellipse, EllipseOptions};
use solidcore::geom::{extent1d::Extent1d, extent2d::Extent2d};
use solidcore::math::{Point3, Vec3};
use solidcore::surface::cone::p1::Cone;
use solidcore::surface::cone::p2::ellipse_trim_intervals;
use solidcore_ffi::units::surface_cone_p2 as ffi;

/// Cone part-2 differential runner.
pub struct Runner;

fn input(case: &Case) -> Result<Vec<f64>, HarnessError> {
    let values: Vec<f64> =
        serde_json::from_value(case.inputs.get("trim").cloned().unwrap_or(Value::Null))
            .map_err(|e| HarnessError::new(e.to_string()))?;
    ffi::validate(&values).map_err(|e| HarnessError::new(e.to_string()))?;
    Ok(values)
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "surface.cone_p2"
    }
    fn cases(&self, ctx: &crate::case::CaseContext) -> Result<Vec<Case>, HarnessError> {
        crate::runner::default_cases(
            self,
            &crate::case::CaseContext {
                seeds: ctx.seeds.max(256),
                ..ctx.clone()
            },
        )
    }
    fn generate(&self, rng: &mut crate::Rng) -> Option<Value> {
        Some(generated_analytic(rng))
    }
    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        if case
            .inputs
            .get("op")
            .and_then(Value::as_str)
            .is_some_and(|op| op != "trim")
        {
            return native_analytic(case);
        }
        ffi::eval_trim(&input(case)?)
            .map(|v| Output::new().values("intervals", v))
            .map_err(|e| HarnessError::new(e.to_string()))
    }
    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        if case
            .inputs
            .get("op")
            .and_then(Value::as_str)
            .is_some_and(|op| op != "trim")
        {
            return rust_analytic(case);
        }
        let values = input(case)?;
        let [bottom, top, height, start, end, swap, inside, tilt, z, rx, ry, a, b, u0, v0, u1, v1, tol] =
            values.as_slice()
        else {
            return Err(HarnessError::new("invalid cone trim record"));
        };
        let kernel = solidcore::Kernel::default();
        let mut cone = Cone::new_bounded(
            Point3::ZERO,
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            *bottom,
            *top,
            *start,
            *end,
            *height,
            *swap != 0.0,
            *inside != 0.0,
            false,
        );
        cone.make_nurb(&kernel)
            .map_err(|e| HarnessError::new(e.to_string()))?;
        let ellipse = Ellipse::new(
            &kernel,
            Point3::new(0.0, 0.0, *z),
            Vec3::new(tilt.cos(), 0.0, tilt.sin()),
            Vec3::new(0.0, 1.0, 0.0),
            Extent1d::from_min_max(0.0, 360.0),
            *rx,
            *ry,
            &EllipseOptions {
                dim: 3,
                ..EllipseOptions::default()
            },
        )
        .map_err(|e| HarnessError::new(e.to_string()))?;
        let intervals = ellipse_trim_intervals(
            &kernel,
            &ellipse,
            &Extent1d::from_min_max(*a, *b),
            &cone,
            &Extent2d::from_bounds(*u0, *v0, *u1, *v1),
            *tol,
        )
        .map_err(|e| HarnessError::new(e.to_string()))?;
        let mut out = vec![intervals.len() as f64];
        for interval in intervals {
            out.push(interval.min());
            out.push(interval.max());
        }
        Ok(Output::new().values("intervals", out))
    }
}

#[derive(serde::Deserialize)]
struct AnalyticOptions {
    #[serde(default)]
    kind: i32,
    #[serde(default)]
    mutation: i32,
    #[serde(default = "one")]
    scale: f64,
    #[serde(default)]
    distance: f64,
    #[serde(default = "small_tolerance")]
    tolerance: f64,
}
fn one() -> f64 {
    1.0
}
fn small_tolerance() -> f64 {
    1e-7
}
fn analytic_options(case: &Case, values: &[f64]) -> Result<(i32, AnalyticOptions), HarnessError> {
    let op = match case.inputs.get("op").and_then(Value::as_str) {
        Some("recognize") => 1,
        Some("rebuild") => 2,
        Some("offset") => 3,
        _ => return Err(HarnessError::new("unknown cone analytic operation")),
    };
    let options: AnalyticOptions = serde_json::from_value(case.inputs.clone())
        .map_err(|e| HarnessError::new(e.to_string()))?;
    ffi::validate_analytic(
        values,
        op,
        options.kind,
        options.mutation,
        options.scale,
        options.distance,
        options.tolerance,
    )
    .map_err(|e| HarnessError::new(e.to_string()))?;
    Ok((op, options))
}
fn native_analytic(case: &Case) -> Result<Output, HarnessError> {
    let values = input(case)?;
    let (op, o) = analytic_options(case, &values)?;
    ffi::eval_analytic(
        &values,
        op,
        o.kind,
        o.mutation,
        o.scale,
        o.distance,
        o.tolerance,
    )
    .map(|v| Output::new().values("state", v))
    .map_err(|e| HarnessError::new(e.to_string()))
}
fn mutate_cone(cone: &mut Cone, mutation: i32) -> solidcore::Result<()> {
    if mutation == 0 {
        return Ok(());
    }
    let swap = cone.rev().swap_uv();
    let nurb = cone.base_mut().get_or_create_gw_nurb()?;
    if mutation == 3 {
        let knots = if swap { nurb.v_knots() } else { nurb.u_knots() };
        let first = knots
            .first()
            .copied()
            .ok_or(solidcore::SolidError::invalid_input(
                "cone fixture",
                "missing knots",
            ))?;
        let last = knots
            .last()
            .copied()
            .ok_or(solidcore::SolidError::invalid_input(
                "cone fixture",
                "missing knots",
            ))?;
        let index = knots.iter().position(|k| *k > first && *k < last);
        if let Some(index) = index {
            let value = (first
                + knots
                    .get(index)
                    .copied()
                    .ok_or(solidcore::SolidError::invalid_input(
                        "cone fixture",
                        "missing knot",
                    ))?)
                / 2.0;
            let target = if swap {
                nurb.v_knots_mut()
            } else {
                nurb.u_knots_mut()
            };
            *target
                .get_mut(index)
                .ok_or(solidcore::SolidError::invalid_input(
                    "cone fixture",
                    "missing knot",
                ))? = value;
        }
        return Ok(());
    }
    if mutation == 4 || mutation == 5 {
        for i in 0..=nurb.n() {
            for j in 0..=nurb.m() {
                let point = nurb.cpt_mut(i, j)?;
                if mutation == 4 {
                    point.z = 0.0;
                } else {
                    use solidcore::nurbs::math::{NOW, NOZ};
                    if point.w != NOW {
                        point.x /= point.w;
                        point.y /= point.w;
                        point.z = if point.z == NOZ {
                            0.0
                        } else {
                            point.z / point.w
                        };
                    }
                    point.w = NOW;
                }
            }
        }
        return Ok(());
    }
    let angle: f64 = 0.05;
    let cosine = angle.cos();
    let sine = angle.sin();
    let last_u = nurb.n();
    let last_v = nurb.m();
    for i in 0..=last_u {
        for j in 0..=last_v {
            if if swap { i != last_u } else { j != last_v } {
                continue;
            }
            let point = nurb.cpt_mut(i, j)?;
            if mutation == 1 {
                let x = point.x;
                let y = point.y;
                point.x = cosine * x - sine * y;
                point.y = sine * x + cosine * y;
            } else {
                point.x *= 2.0;
                point.y *= 2.0;
                point.z *= 2.0;
                point.w *= 2.0;
            }
        }
    }
    Ok(())
}
fn state_vec(out: &mut Vec<f64>, point: Vec3) {
    out.extend([point.x, point.y, point.z]);
}
fn analytic_state(
    kernel: &solidcore::Kernel,
    out: &mut Vec<f64>,
    cone: &Cone,
) -> solidcore::Result<()> {
    use solidcore::curve::curve::CurveKind;
    let rev = cone.rev();
    let position = rev.position();
    out.extend([
        cone.type_code() as f64,
        cone.base_radius(),
        cone.semi_angle_deg(),
        f64::from(u8::from(rev.swap_uv())),
        f64::from(u8::from(cone.inside_out())),
    ]);
    state_vec(out, position.origin());
    state_vec(out, position.x_axis());
    state_vec(out, position.y_axis());
    for domain in [rev.step_uv_domain(), rev.base().natural_uv_domain()] {
        out.extend([
            domain.min().x,
            domain.min().y,
            domain.max().x,
            domain.max().y,
        ]);
    }
    out.extend([cone.bot_radius(), cone.top_radius(), cone.height(kernel)?]);
    let generator = rev
        .gen_curve()
        .ok_or(solidcore::SolidError::null_pointer("cone analytic fixture"))?;
    let CurveKind::Line(line) = generator.kind() else {
        return Err(solidcore::SolidError::invalid_input(
            "cone analytic fixture",
            "unexpected generator kind",
        ));
    };
    let step = line.step_interval();
    let nurb = generator.natural_interval();
    out.extend([
        generator.type_code() as f64,
        step.min(),
        step.max(),
        nurb.min(),
        nurb.max(),
    ]);
    for fraction in [0.0, 0.4, 1.0] {
        state_vec(out, generator.evaluate_point(nurb.evaluate(fraction))?);
    }
    let angles = rev.polar_converter().angles();
    out.extend([
        f64::from(u8::from(rev.have_polar_conversion())),
        angles.len() as f64,
    ]);
    if let (Some(first), Some(last)) = (angles.first(), angles.last()) {
        out.extend([*first, *last]);
    }
    let nurb = rev.base().gw_nurb();
    out.push(f64::from(u8::from(nurb.is_some())));
    if let Some(nurb) = nurb {
        out.extend([
            nurb.p() as f64,
            nurb.q() as f64,
            (nurb.n() + 1) as f64,
            (nurb.m() + 1) as f64,
        ]);
        for row in nurb.rows() {
            for point in row {
                out.extend([point.x, point.y, point.z, point.w]);
            }
        }
        out.push(nurb.u_knots().len() as f64);
        out.extend_from_slice(nurb.u_knots());
        out.push(nurb.v_knots().len() as f64);
        out.extend_from_slice(nurb.v_knots());
    }
    Ok(())
}
fn operation_status(error: &solidcore::SolidError) -> f64 {
    match error {
        solidcore::SolidError::Upstream { code, .. } => f64::from(*code),
        solidcore::SolidError::InvalidInput { .. } => 1007.0,
        _ => 1001.0,
    }
}
fn rust_analytic(case: &Case) -> Result<Output, HarnessError> {
    let values = input(case)?;
    let (op, o) = analytic_options(case, &values)?;
    let [bottom, top, height, start, end, swap, inside, ..] = values.as_slice() else {
        return Err(HarnessError::new("invalid cone geometry"));
    };
    let kernel = solidcore::Kernel::default();
    let (origin, x_axis, y_axis) =
        if let Some(&[ox, oy, oz, xx, xy, xz, yx, yy, yz]) = values.get(18..) {
            (
                Point3::new(ox, oy, oz),
                Vec3::new(xx, xy, xz),
                Vec3::new(yx, yy, yz),
            )
        } else {
            (
                Point3::ZERO,
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            )
        };
    let mut cone = if o.kind & 4 != 0 {
        if o.kind & 1 == 0 {
            Cone::new_unbounded(
                &kernel,
                origin,
                x_axis,
                y_axis,
                *bottom,
                values
                    .get(7)
                    .copied()
                    .ok_or_else(|| HarnessError::new("missing semi-angle"))?
                    * 180.0
                    / std::f64::consts::PI,
                *swap != 0.0,
                *inside != 0.0,
            )
        } else {
            Cone::new_cylinder_unbounded(
                &kernel,
                origin,
                x_axis,
                y_axis,
                *bottom,
                *swap != 0.0,
                *inside != 0.0,
            )
        }
        .map_err(|e| HarnessError::new(e.to_string()))?
    } else if o.kind & 1 == 0 {
        Cone::new_bounded(
            origin,
            x_axis,
            y_axis,
            *bottom,
            *top,
            *start,
            *end,
            *height,
            *swap != 0.0,
            *inside != 0.0,
            false,
        )
    } else {
        Cone::new_cylinder_bounded(
            origin,
            x_axis,
            y_axis,
            *bottom,
            *start,
            *end,
            *height,
            *swap != 0.0,
            *inside != 0.0,
            false,
        )
    };
    // The native SmBSplineSurface copy constructor forces the source's
    // virtual NURB factory before producing recognition's pure B-spline input.
    if o.kind & 2 == 0 || op == 1 {
        cone.make_nurb(&kernel)
            .map_err(|e| HarnessError::new(e.to_string()))?;
    }
    mutate_cone(&mut cone, o.mutation).map_err(|e| HarnessError::new(e.to_string()))?;
    let mut status = 1000.0;
    let result = match op {
        1 => match Cone::is_nurb_surface_cone(&kernel, cone.rev().base(), o.scale) {
            Ok(value) => {
                status = f64::from(u8::from(value.is_some()));
                value
            }
            Err(error) => {
                if let Some(output) = Output::unported(&error) {
                    return Ok(output);
                }
                status = operation_status(&error);
                None
            }
        },
        2 => {
            if let Err(error) = cone.rebuild_step_from_nurb_parameters(&kernel) {
                if let Some(output) = Output::unported(&error) {
                    return Ok(output);
                }
                status = operation_status(&error);
            }
            Some(cone)
        }
        3 => match cone.create_offset_surface(&kernel, o.distance, o.tolerance) {
            Ok(value) => value,
            Err(error) => {
                if let Some(output) = Output::unported(&error) {
                    return Ok(output);
                }
                status = operation_status(&error);
                None
            }
        },
        _ => return Err(HarnessError::new("invalid cone operation")),
    };
    let mut out = vec![status, f64::from(u8::from(result.is_some()))];
    if let Some(cone) = result {
        analytic_state(&kernel, &mut out, &cone).map_err(|e| HarnessError::new(e.to_string()))?;
    }
    Ok(Output::new().values("state", out))
}

fn generated_analytic(rng: &mut crate::Rng) -> Value {
    let scale = match rng.below(3) {
        0 => 1e-7,
        1 => 1.0,
        _ => 1e6,
    };
    let kind = rng.below(4) as i32;
    let bottom = rng.uniform(0.25, 8.0) * scale;
    let top = if kind % 2 == 1 {
        bottom
    } else if rng.chance(0.15) {
        0.0
    } else {
        rng.uniform(0.1, 8.0) * scale
    };
    let height = rng.uniform(1.0, 30.0) * scale;
    let start = match rng.below(3) {
        0 => 0.0,
        1 => 20.0,
        _ => 45.0,
    };
    let end = match rng.below(3) {
        0 => 180.0,
        1 => 250.0,
        _ => 360.0,
    };
    let swap = f64::from(u8::from(rng.chance(0.5)));
    let inside = f64::from(u8::from(rng.chance(0.5)));
    let op = match rng.below(3) {
        0 => "recognize",
        1 => "rebuild",
        _ => "offset",
    };
    let mutation = if kind < 2 && op != "offset" && rng.chance(0.25) {
        1 + rng.below(2) as i32
    } else {
        0
    };
    let tolerance = match rng.below(4) {
        0 => 0.0,
        1 => 1e-12 * scale,
        2 => 1e-7 * scale,
        _ => 1e-3 * scale,
    };
    let distance = match rng.below(5) {
        0 => -0.5 * bottom,
        1 => -1.25 * bottom,
        2 => 0.0,
        3 => 0.25 * bottom,
        _ => bottom,
    };
    serde_json::json!({"op":op,"kind":kind,"mutation":mutation,"scale":1.0,
        "distance":distance,"tolerance":tolerance,
        "trim":[bottom,top,height,start,end,swap,inside,0.0,height/2.0,scale,scale,0.0,360.0,0.0,0.0,1.0,1.0,1e-7*scale]})
}
