//! query.brep: actual `SM_API` evaluation entry points, status and failure atomicity.
//! Initial surface API batch only; this does not certify the complete B-rep unit.
use super::surface_surface_cache_p1::{push3, status};
use super::surface_surface_cache_p3::{gen_surf, push_surf, surf_domain, surface_in, SurfIn};
use crate::case::CaseContext;
use crate::runner::default_cases;
use crate::{Case, HarnessError, Output, Rng, UnitDiff};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use solidcore::geom::vector2d::Vec2;
use solidcore::query::brep;
use solidcore_ffi::units::query_brep as ffi;

/// Native surface API comparison for `query.brep`.
pub struct Runner;
#[derive(Clone, Serialize, Deserialize)]
struct Inputs {
    surf: Option<SurfIn>,
    op: u8,
    uv: [f64; 2],
    /// 0 ordinary; 1/2 NaN u/v; 3/4 infinite u/v. JSON cannot encode NaN.
    #[serde(default)]
    special: u8,
    initial: [[f64; 3]; 3],
}
fn uv(i: &Inputs) -> [f64; 2] {
    let mut uv = i.uv;
    match i.special {
        1 => uv[0] = f64::NAN,
        2 => uv[1] = f64::NAN,
        3 => uv[0] = f64::INFINITY,
        4 => uv[1] = f64::NEG_INFINITY,
        _ => {}
    }
    uv
}
fn flat(i: &Inputs) -> Result<Vec<f64>, HarnessError> {
    let mut out = vec![if i.surf.is_some() { 1.0 } else { 0.0 }];
    if let Some(s) = &i.surf {
        push_surf(&mut out, s)?;
    }
    out.extend(uv(i));
    out.extend(i.initial.iter().flatten());
    Ok(out)
}
fn record(out: Vec<f64>) -> Output {
    Output::new().values("status_and_outputs", out)
}
impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "query.brep"
    }
    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let out = ffi::surface(i.op, &flat(&i)?).map_err(|e| HarnessError::new(e.to_string()))?;
        Ok(record(out))
    }
    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let s = i.surf.as_ref().map(surface_in).transpose()?;
        let [u, v] = uv(&i);
        let uv = Vec2::new(u, v);
        let mut result_status = 0.0;
        let mut outputs: Vec<f64> = i.initial.iter().flatten().copied().collect();
        let evaluated = match i.op {
            0 => brep::evaluate_surface_point(s.as_ref(), uv).map(|p| {
                let mut out = Vec::new();
                push3(&mut out, p);
                out
            }),
            1 => brep::evaluate_surface_normal(s.as_ref(), uv).map(|p| {
                let mut out = Vec::new();
                push3(&mut out, p);
                out
            }),
            2 => brep::evaluate_surface_derivatives(s.as_ref(), uv).map(|(p, du, dv)| {
                let mut out = Vec::new();
                push3(&mut out, p);
                push3(&mut out, du);
                push3(&mut out, dv);
                out
            }),
            _ => return Err(HarnessError::new("surface operation")),
        };
        match evaluated {
            Ok(v) => {
                for (target, value) in outputs.iter_mut().zip(v) {
                    *target = value;
                }
            }
            Err(e) => {
                if let Some(deferred) = Output::unported(&e) {
                    return Ok(deferred);
                }
                result_status = status(&e);
            }
        }
        let mut result = vec![result_status];
        result.extend(outputs);
        Ok(record(result))
    }
    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        for _ in 0..100 {
            let surf = gen_surf(rng);
            let [a, b, c, d] = surf_domain(&surf).ok()?;
            let fraction = match rng.below(5) {
                0 => 0.0,
                1 => 1.0,
                _ => rng.unit(),
            };
            let mut uv = [a + fraction * (c - a), b + rng.unit() * (d - b)];
            if rng.below(8) == 0 {
                uv[0] = a - (c - a);
            }
            let i = Inputs {
                surf: Some(surf),
                op: rng.below(3) as u8,
                uv,
                special: if rng.below(10) == 0 {
                    1 + rng.below(4) as u8
                } else {
                    0
                },
                initial: [
                    [17.0, -29.0, 41.0],
                    [53.0, -67.0, 79.0],
                    [83.0, -97.0, 101.0],
                ],
            };
            if ffi::validate(i.op, &flat(&i).ok()?).is_ok() {
                return serde_json::to_value(i).ok();
            }
        }
        None
    }
    fn cases(&self, ctx: &CaseContext) -> Result<Vec<Case>, HarnessError> {
        default_cases(
            self,
            &CaseContext {
                seeds: ctx.seeds.max(600),
                ..ctx.clone()
            },
        )
    }
}
