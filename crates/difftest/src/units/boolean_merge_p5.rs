//! Wire/shell-vertex boolean stage on identical raw input topology. Construction
//! uses the shared real-B-rep fixture builder; no C++ merge results enter Rust.
use super::topology_vertexuse::{self as shared, BrepIn};
use crate::{Case, HarnessError, Output, Rng, UnitDiff};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use solidcore::boolean::merge::p5::{process_wires_and_shell_vertices, WireRelations};
use solidcore::geom::a_object::BooleanOperation;
use solidcore::geom::context::MarkLock;
use solidcore::kernel::Kernel;
use solidcore::topology::topology_traverser::{Entity, Marks};
use solidcore::topology::{edge, vertex, vertexuse, ArenaId, Brep, EdgeId, VertexId, VertexuseId};
use solidcore_ffi::units::boolean_merge_p5::{self as ffi, Operand, Request};
use std::collections::BTreeMap;

/// Differential runner for `SmMerge::ProcessWiresAndShellVertices`.
pub struct Runner;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum Spec {
    Empty,
    Vertex { point: [f64; 3] },
    Graph { brep: BrepIn },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Inputs {
    primary: Spec,
    other: Spec,
    operation: i32,
    /// Bit 1 marks edges, bit 2 marks vertices.
    #[serde(default)]
    primary_pre: u8,
    #[serde(default)]
    other_pre: u8,
    #[serde(default)]
    primary_unmark: Vec<usize>,
    #[serde(default)]
    other_unmark: Vec<usize>,
    #[serde(default)]
    mates: Option<Vec<i64>>,
    #[serde(default = "default_marks")]
    marks: [i32; 2],
    #[serde(default)]
    prefill: bool,
    #[serde(default = "one")]
    repeat: i32,
}
fn default_marks() -> [i32; 2] {
    [1, 2]
}
fn one() -> i32 {
    1
}
fn err(e: &(impl ToString + ?Sized)) -> HarnessError {
    HarnessError::new(e.to_string())
}

fn program(s: &Spec) -> Result<(String, Vec<f64>), HarnessError> {
    let mut block = Vec::new();
    let path = match s {
        Spec::Empty => {
            block.push(5.0);
            String::new()
        }
        Spec::Vertex { point } => {
            block.push(4.0);
            block.extend(point);
            String::new()
        }
        Spec::Graph { brep } => {
            match brep {
                BrepIn::File { rebuild, .. } => block.extend([0.0, f64::from(*rebuild)]),
                BrepIn::Wire { curve, t0, t1 } => {
                    block.push(1.0);
                    shared::push_curve(&mut block, curve)?;
                    block.extend([t0, t1]);
                }
                BrepIn::Face { surface, dom } => {
                    block.push(2.0);
                    shared::push_surface(&mut block, surface)?;
                    block.extend(dom);
                }
                BrepIn::Network => block.push(3.0),
            }
            block.extend([0.0, 0.0]);
            shared::path_of(brep)?
        }
    };
    Ok((path, block))
}
fn operand(p: &(String, Vec<f64>)) -> Operand<'_> {
    Operand {
        path: &p.0,
        block: &p.1,
    }
}
fn flags(c: [usize; 10], pre: u8, unmark: &[usize]) -> Vec<i64> {
    (0..c[7] + c[9])
        .map(|i| i64::from(pre & if i < c[7] { 1 } else { 2 } != 0 && !unmark.contains(&i)))
        .collect()
}
fn mapping(i: &Inputs, a: [usize; 10], b: [usize; 10]) -> Vec<i64> {
    i.mates.clone().unwrap_or_else(|| {
        (0..b[7])
            .map(|j| if j < a[7] { j as i64 } else { -1 })
            .collect()
    })
}
fn import(
    f: &solidcore_ffi::units::topology_topology_traverser::Fixture,
) -> Result<Brep, HarnessError> {
    let (mut b, _) = shared::import(&f.graph)?;
    for (i, list) in f.vertexuses.iter().enumerate() {
        b.vertex_mut(VertexId::from_index(i as u32))
            .map_err(|e| err(&e))?
            .vertexuses = list.iter().map(|&j| VertexuseId::from_index(j)).collect();
    }
    Ok(b)
}
fn seed(
    k: &Kernel,
    b: &Brep,
    marks: &mut Marks,
    mt: solidcore::geom::context::MarkType,
    f: &[i64],
) {
    for (i, e) in b.edges().ids().enumerate() {
        if f.get(i) == Some(&1) {
            marks.mark(k, Entity::Edge(e), mt);
        }
    }
    for (i, v) in b.vertices().ids().enumerate() {
        if f.get(b.edges().len() + i) == Some(&1) {
            marks.mark(k, Entity::Vertex(v), mt);
        }
    }
}
fn id<I: ArenaId>(value: Option<I>) -> i64 {
    value.map_or(-1, |i| i64::from(i.index()))
}
fn snapshot(
    k: &Kernel,
    b: &Brep,
    marks: &Marks,
    mt: solidcore::geom::context::MarkType,
    ids: &mut Vec<i64>,
    values: &mut Vec<f64>,
) -> Result<(), HarnessError> {
    ids.extend(
        [
            b.regions().len(),
            b.shells().len(),
            b.faceuses().len(),
            b.faces().len(),
            b.loopuses().len(),
            b.loops().len(),
            b.edgeuses().len(),
            b.edges().len(),
            b.vertexuses().len(),
            b.vertices().len(),
        ]
        .map(|n| n as i64),
    );
    for (e, r) in b.edges().iter() {
        let v = edge::vertex(b, e).map_err(|e| err(&e))?;
        ids.extend([
            id(v),
            id(edge::other_vertex(b, e, v).map_err(|e| err(&e))?),
            i64::from(vertexuse::edge_is_wire(b, e).map_err(|e| err(&e))?),
            i64::from(marks.is_marked(k, Entity::Edge(e), mt)),
            r.edgeuses.len() as i64,
        ]);
        ids.extend(r.edgeuses.iter().map(|e| i64::from(e.index())));
        values.extend([r.interval.min(), r.interval.max(), r.tolerance]);
    }
    for (v, r) in b.vertices().iter() {
        ids.extend([
            i64::from(vertex::is_shell_vertex(b, v).map_err(|e| err(&e))?),
            i64::from(marks.is_marked(k, Entity::Vertex(v), mt)),
            r.vertexuses.len() as i64,
        ]);
        ids.extend(r.vertexuses.iter().map(|v| i64::from(v.index())));
        values.extend([r.point.x, r.point.y, r.point.z, r.tolerance]);
    }
    for (_, r) in b.regions().iter() {
        ids.push(i64::from(r.in_brep));
    }
    for (_, sh) in b.shells().iter() {
        ids.push(id(sh.region));
    }
    for (_, f) in b.faces().iter() {
        ids.push(id(f.faceuse));
        values.extend([
            f.tolerance,
            f.uv_domain.min().x,
            f.uv_domain.min().y,
            f.uv_domain.max().x,
            f.uv_domain.max().y,
        ]);
    }
    for (_, fu) in b.faceuses().iter() {
        ids.extend([
            id(fu.shell),
            id(fu.face),
            id(fu.mate),
            i64::from(fu.orientation as i32),
            fu.loopuses.len() as i64,
        ]);
        ids.extend(fu.loopuses.iter().map(|l| i64::from(l.index())));
    }
    for (_, l) in b.loops().iter() {
        ids.push(id(l.loopuse));
    }
    for (_, lu) in b.loopuses().iter() {
        ids.extend([
            id(lu.faceuse),
            id(lu.loop_),
            id(lu.mate),
            i64::from(lu.orientation as i32),
        ]);
        ids.extend(match lu.start {
            solidcore::topology::LoopuseStart::Unknown => [0, -1],
            solidcore::topology::LoopuseStart::Edgeuse(e) => [1, id(e)],
            solidcore::topology::LoopuseStart::Vertexuse(v) => [2, id(v)],
        });
    }
    for (_, eu) in b.edgeuses().iter() {
        ids.extend([id(eu.edge), id(eu.vertexuse), id(eu.cw), id(eu.ccw)]);
        ids.extend(match eu.owner {
            solidcore::topology::EdgeuseOwner::Unknown => [0, -1],
            solidcore::topology::EdgeuseOwner::Loopuse(l) => [1, id(l)],
            solidcore::topology::EdgeuseOwner::Shell(sh) => [2, id(sh)],
        });
        ids.push(i64::from(eu.orientation as i32));
    }
    for (_, vu) in b.vertexuses().iter() {
        ids.push(id(vu.vertex));
        ids.extend(match vu.owner {
            solidcore::topology::VertexuseOwner::Unknown => [0, -1],
            solidcore::topology::VertexuseOwner::Edgeuse(e) => [1, id(e)],
            solidcore::topology::VertexuseOwner::Loopuse(l) => [2, id(l)],
            solidcore::topology::VertexuseOwner::Shell(sh) => [3, id(sh)],
        });
    }
    Ok(())
}
fn output(raw: &[f64]) -> Result<Output, HarnessError> {
    let n = *raw.first().ok_or_else(|| err("empty result"))? as usize;
    let end = n.checked_add(1).ok_or_else(|| err("output overflow"))?;
    let ints: Vec<i64> = raw
        .get(1..end)
        .ok_or_else(|| err("truncated topology"))?
        .iter()
        .map(|&x| x as i64)
        .collect();
    let values = raw
        .get(end..)
        .ok_or_else(|| err("truncated values"))?
        .to_vec();
    Ok(Output::new()
        .ints("topology", ints)
        .values("geometry", values))
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "boolean.merge_p5"
    }
    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let a = program(&i.primary)?;
        let b = program(&i.other)?;
        let aa = ffi::export(operand(&a)).map_err(|e| err(&e))?;
        let bb = ffi::export(operand(&b)).map_err(|e| err(&e))?;
        let ac = solidcore_ffi::units::topology_vertexuse::counts(&aa.graph)
            .ok_or_else(|| err("primary graph"))?;
        let bc = solidcore_ffi::units::topology_vertexuse::counts(&bb.graph)
            .ok_or_else(|| err("other graph"))?;
        let am = flags(ac, i.primary_pre, &i.primary_unmark);
        let bm = flags(bc, i.other_pre, &i.other_unmark);
        let mates = mapping(&i, ac, bc);
        let r = Request {
            operation: i.operation,
            marks: i.marks,
            primary_marked: &am,
            other_marked: &bm,
            mates: &mates,
            prefill: i.prefill,
            repeat: i.repeat,
        };
        output(&ffi::eval(operand(&a), operand(&b), r).map_err(|e| err(&e))?)
    }
    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let a = program(&i.primary)?;
        let b = program(&i.other)?;
        let aa = ffi::export(operand(&a)).map_err(|e| err(&e))?;
        let bb = ffi::export(operand(&b)).map_err(|e| err(&e))?;
        let ac = solidcore_ffi::units::topology_vertexuse::counts(&aa.graph)
            .ok_or_else(|| err("primary graph"))?;
        let bc = solidcore_ffi::units::topology_vertexuse::counts(&bb.graph)
            .ok_or_else(|| err("other graph"))?;
        let mut a = import(&aa)?;
        let mut b = import(&bb)?;
        let mut k = Kernel::default();
        let la = MarkLock::new(std::slice::from_mut(&mut k.context), 0, i.marks[0])
            .map_err(|e| err(&e))?;
        let lb = if i.marks[0] == i.marks[1] {
            None
        } else {
            Some(
                MarkLock::new(std::slice::from_mut(&mut k.context), 0, i.marks[1])
                    .map_err(|e| err(&e))?,
            )
        };
        let ma = la.mark_type();
        let mb = lb.as_ref().map_or(ma, MarkLock::mark_type);
        let mut am = Marks::default();
        let mut bm = Marks::default();
        seed(
            &k,
            &a,
            &mut am,
            ma,
            &flags(ac, i.primary_pre, &i.primary_unmark),
        );
        seed(
            &k,
            &b,
            &mut bm,
            mb,
            &flags(bc, i.other_pre, &i.other_unmark),
        );
        let mut mates = BTreeMap::new();
        // Fixture setup mirrors SmTopologyIntersector::Relate, including the
        // tolerance propagation it performs BEFORE the unit under test is called.
        for (j, m) in mapping(&i, ac, bc).into_iter().enumerate() {
            if m < 0 {
                continue;
            }
            let e = EdgeId::from_index(m as u32);
            let oe = EdgeId::from_index(j as u32);
            let at = a.edge(e).map_err(|e| err(&e))?.tolerance;
            let bt = b.edge(oe).map_err(|e| err(&e))?.tolerance;
            if at < bt {
                edge::set_tolerance(&mut a, e, bt, true, true).map_err(|e| err(&e))?;
            }
            if bt < at {
                edge::set_tolerance(&mut b, oe, at, true, true).map_err(|e| err(&e))?;
            }
            mates.insert(oe, e);
        }
        let vertices = BTreeMap::new();
        let relations = WireRelations {
            other_edges: &mates,
            vertices: &vertices,
        };
        let operation =
            BooleanOperation::from_code(i.operation).ok_or_else(|| err("boolean operation"))?;
        let mut keep = Vec::new();
        let mut ids = Vec::new();
        for _ in 0..i.repeat {
            if i.prefill && ac[7] > 0 {
                keep.push(EdgeId::from_index(0));
            }
            match process_wires_and_shell_vertices(
                &k, &mut a, &b, &mut am, &mut bm, relations, operation, ma, mb, &mut keep,
            ) {
                Ok(()) => {
                    ids.extend([0, keep.len() as i64]);
                    ids.extend(keep.iter().map(|e| i64::from(e.index())));
                }
                Err(e) => return Output::unported(&e).ok_or_else(|| err(&e)),
            }
        }
        let mut values = Vec::new();
        snapshot(&k, &a, &am, ma, &mut ids, &mut values)?;
        snapshot(&k, &b, &bm, mb, &mut ids, &mut values)?;
        Ok(Output::new()
            .ints("topology", ids)
            .values("geometry", values))
    }
    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        let scale = *rng.choose(&[1.0e-9, 1.0, 1.0e9])?;
        let graph = match rng.below(4) {
            0 => json!({"kind":"vertex","point":rng.point3(-scale,scale)}),
            1 => json!({"kind":"graph","brep":{"kind":"network"}}),
            2 => {
                // Shared plane fixtures support a bounded UV domain. Keep the
                // large world-space translation while using a valid sub-domain.
                let extent = scale.min(1.0e3);
                json!({"kind":"graph","brep":{"kind":"face","surface":{"plane":{"o":rng.point3(-scale,scale),"n":[0.0,0.0,1.0]}},"dom":[-extent,-extent,extent,extent]}})
            }
            _ => {
                let p = rng.point3(-scale, scale);
                let q = [p[0] + scale, p[1] + scale, p[2] + scale];
                json!({"kind":"graph","brep":{"kind":"wire","curve":{"dim":3,"cpts":[[p[0],p[1],p[2],-f64::MAX],[q[0],q[1],q[2],-f64::MAX]],"p":1,"u":[0.0,0.0,1.0,1.0]},"t0":0.2,"t1":0.8}})
            }
        };
        let mask_pairs = [[1, 2], [4, 8], [8, 16], [2, 4], [16, 1], [1, 1]];
        Some(
            json!({"primary":graph,"other":graph,"operation":rng.below(11),
            "primary_pre":3,"other_pre":3,"marks":rng.choose(&mask_pairs)?,"prefill":true,"repeat":2}),
        )
    }
}
