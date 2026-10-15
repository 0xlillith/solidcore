//! WIP native shell/region query comparisons. Other owned operations are pending.
//! Native fixtures and raw graph imports reuse the edge-use fixture protocol.
//! Setup failures fail the harness; topology and native integer flags compare exactly.

use super::topology_edgeuse_p1::{fixture_input, import};
use crate::{Case, HarnessError, Output, Rng, UnitDiff};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use solidcore::topology::shell::{shell_brep, shell_region};
use solidcore::topology::{Region, Shell};
use solidcore_ffi::units::{topology_edgeuse_p1 as graph, topology_shell as ffi};

/// The shell/region unit runner (coverage is still under development).
pub struct Runner;

#[derive(Serialize, Deserialize)]
struct Inputs {
    brep: Value,
    #[serde(default)]
    void_state: Option<i32>,
}

const TYPES: [i64; 5] = [16013, 16012, 16000, 16999, -1];

fn failure(e: impl std::fmt::Display) -> HarnessError {
    HarnessError::new(e.to_string())
}
fn flag(value: bool) -> i64 {
    i64::from(value)
}
fn list(out: &mut Vec<i64>, values: impl IntoIterator<Item = u32>) {
    let values: Vec<_> = values.into_iter().collect();
    out.push(values.len() as i64);
    out.extend(values.into_iter().map(i64::from));
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "topology.shell"
    }

    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let inputs: Inputs = case.parse()?;
        let (path, block) = fixture_input(inputs.brep)?;
        let native = ffi::query(&path, &block, inputs.void_state).map_err(failure)?;
        if native
            .iter()
            .any(|v| !v.is_finite() || v.fract() != 0.0 || v.abs() > (1_u64 << 53) as f64)
        {
            return Err(failure("native topology query returned an invalid integer"));
        }
        Ok(Output::new().ints(
            "records",
            native.into_iter().map(|v| v as i64).collect::<Vec<_>>(),
        ))
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let inputs: Inputs = case.parse()?;
        let (path, block) = fixture_input(inputs.brep)?;
        let export = graph::eval(graph::Op::Export, &path, &block).map_err(failure)?;
        let mut brep = import(&export)?;
        let mut out = vec![brep.shells().len() as i64, brep.regions().len() as i64];
        for (id, shell) in brep.shells().iter() {
            out.extend([
                Shell::type_code(),
                flag(shell.is_vertex_shell()),
                flag(shell.is_wire_shell()),
                flag(shell.is_faceuse_shell()),
                flag(shell.is_unknown_shell()),
                flag(shell.is_closed(&brep, id).map_err(failure)?),
                shell_region(&brep, id)
                    .map_err(failure)?
                    .map_or(-1, |v| i64::from(v.index())),
                flag(shell_brep(&brep, id).map_err(failure)?.is_some()),
                shell
                    .first_faceuse()
                    .map_err(failure)?
                    .map_or(-1, |v| i64::from(v.index())),
                shell
                    .wire_edge(&brep)
                    .map_err(failure)?
                    .map_or(-1, |v| i64::from(v.index())),
                shell
                    .vertex(&brep)
                    .map_err(failure)?
                    .map_or(-1, |v| i64::from(v.index())),
            ]);
            out.extend(TYPES.into_iter().map(|ty| flag(Shell::is_kind_of(ty))));
            list(
                &mut out,
                shell
                    .get_faceuses()
                    .map_err(failure)?
                    .into_iter()
                    .map(|v| v.index()),
            );
            list(
                &mut out,
                shell
                    .mated_shells(&brep)
                    .map_err(failure)?
                    .into_iter()
                    .map(|v| v.index()),
            );
        }
        let regions: Vec<_> = brep.regions().ids().collect();
        for id in regions {
            if let Some(value) = inputs.void_state {
                brep.region_mut(id).map_err(failure)?.set_is_void(value);
            }
            let region = brep.region(id).map_err(failure)?;
            out.extend([
                Region::type_code(),
                i64::from(region.is_void()),
                flag(region.is_infinite_region(&brep, id)),
                flag(region.get_brep(&brep).is_some()),
            ]);
            out.extend(TYPES.into_iter().map(|ty| flag(Region::is_kind_of(ty))));
            list(
                &mut out,
                region
                    .get_shells(None, |_, _| Ok(false))
                    .map_err(failure)?
                    .into_iter()
                    .map(|v| v.index()),
            );
        }
        Ok(Output::new().ints("records", out))
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        let scale = [1e-3, 1.0, 1e3].get(rng.below(3) as usize).copied()?;
        let origin = rng.point3(-2.0, 2.0).map(|v| v * scale);
        let normal = [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), 1.0];
        Some(
            json!({"brep": {"kind":"face", "surface":{"plane":{"o":origin,"n":normal}},
            "dom":[-scale,-scale,scale,scale]}, "void_state":rng.below(3) as i32}),
        )
    }
}
