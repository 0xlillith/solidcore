//! Native loop topology, orientation and ordered-list differential queries.
use super::topology_edgeuse_p1::{brep_block_of, gen_brep, import, path_of, BrepIn};
use crate::{Case, HarnessError, Output, Rng, UnitDiff};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType;
use solidcore::topology::loopuse as api;
use solidcore::topology::{ArenaId, Brep, Loop, Loopuse};
use solidcore_ffi::units::{topology_edgeuse_p1 as graph, topology_loopuse as ffi};

/// Differential runner for the in-progress loop topology port.
pub struct Runner;
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Inputs {
    brep: BrepIn,
    #[serde(default)]
    orientation: Option<u8>,
}
fn failure(e: impl std::fmt::Display) -> HarnessError {
    HarnessError::new(e.to_string())
}
fn fixture(i: &Inputs) -> Result<(String, Vec<f64>), HarnessError> {
    if i.orientation.is_some_and(|v| v > 2) {
        return Err(failure("invalid orientation"));
    }
    Ok((
        path_of(Some(&i.brep))?,
        brep_block_of(Some(&i.brep), false, false)?,
    ))
}
fn orient(code: u8) -> OrientType {
    match code {
        1 => OrientType::Same,
        2 => OrientType::Opposite,
        _ => OrientType::Unknown,
    }
}
fn code(value: OrientType) -> i64 {
    match value {
        OrientType::Unknown => 0,
        OrientType::Same => 1,
        OrientType::Opposite => 2,
    }
}
fn id<T: ArenaId>(value: Option<T>) -> i64 {
    value.map_or(-1, |v| i64::from(v.index()))
}
fn list<T: ArenaId>(out: &mut Vec<i64>, values: impl IntoIterator<Item = Option<T>>) {
    let values: Vec<_> = values.into_iter().map(id).collect();
    out.push(values.len() as i64);
    out.extend(values);
}
fn text(out: &mut Vec<i64>, value: &str) {
    out.push(value.len() as i64);
    out.extend(value.bytes().map(i64::from));
}
fn query(b: &mut Brep, orientation: Option<u8>) -> solidcore::Result<Vec<i64>> {
    let uses: Vec<_> = b.loopuses().ids().collect();
    let mut out = vec![uses.len() as i64, b.loops().len() as i64];
    for lu in uses {
        if let Some(value) = orientation {
            api::set_orientation(b, lu, orient(value))?;
        }
        let record = *b.loopuse(lu)?;
        let ty = api::loopuse_type(b, lu)?;
        out.extend([
            Loopuse::type_code(),
            ty,
            code(api::orientation(b, lu)?),
            i64::from(ty == 16005),
            i64::from(ty == 16006),
            i64::from(api::is_closed_ptrs(b, lu)?),
            id(record.faceuse),
            id(record.loop_),
            id(record.mate),
            id(solidcore::topology::edgeuse::p1::loopuse_shell(b, lu)?),
            id(api::vertexuse(b, lu)?),
        ]);
        for t in [16003, 16004, 16001, 16010, 12060, 11001, 16999, -1] {
            out.push(i64::from(Loopuse::is_kind_of(t)));
        }
        text(&mut out, Loopuse::class_string());
        text(&mut out, Loopuse::type_string());
        list(
            &mut out,
            api::edgeuses(b, lu, None, None, |_, _| Ok(true))?
                .into_iter()
                .map(Some),
        );
        list(&mut out, api::edges(b, lu, None, |_, _| Ok(true))?);
        list(&mut out, api::vertexuses(b, lu, None, |_, _| Ok(true))?);
        list(&mut out, api::vertices(b, lu, None, |_, _| Ok(true))?);
    }
    for (l, record) in b.loops().iter() {
        let (primary, other) = api::loopuses(b, l)?;
        out.extend([
            Loop::type_code(),
            id(record.loopuse),
            id(Some(primary)),
            id(other),
            id(api::loop_face(b, l)?),
            i64::from(api::is_outer_loop(b, l)?),
            i64::from(api::loop_is_closed_ptrs(b, l)?),
        ]);
        for t in [16003, 16004, 16001, 16010, 12060, 11001, 16999, -1] {
            out.push(i64::from(Loop::is_kind_of(t)));
        }
        text(&mut out, Loop::class_string());
        text(&mut out, Loop::type_string());
        list(
            &mut out,
            api::edgeuses(b, primary, None, None, |_, _| Ok(true))?
                .into_iter()
                .map(Some),
        );
        list(&mut out, api::edges(b, primary, None, |_, _| Ok(true))?);
        list(&mut out, api::vertices(b, primary, None, |_, _| Ok(true))?);
    }
    Ok(out)
}
impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "topology.loopuse"
    }
    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let (path, block) = fixture(&i)?;
        let data = ffi::query(&path, &block, i.orientation).map_err(failure)?;
        if data
            .iter()
            .any(|v| !v.is_finite() || v.fract() != 0.0 || v.abs() > (1_u64 << 53) as f64)
        {
            return Err(failure("invalid native integer record"));
        }
        Ok(Output::new().ints(
            "records",
            data.into_iter().map(|v| v as i64).collect::<Vec<_>>(),
        ))
    }
    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let (path, block) = fixture(&i)?;
        let ex = graph::eval(graph::Op::Export, &path, &block).map_err(failure)?;
        let mut b = import(&ex)?;
        Ok(Output::new().ints("records", query(&mut b, i.orientation).map_err(failure)?))
    }
    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        serde_json::to_value(Inputs {
            brep: gen_brep(rng),
            orientation: if rng.chance(0.5) {
                Some(rng.below(3) as u8)
            } else {
                None
            },
        })
        .ok()
    }
}
