//! Runner for unit `blend.hcr`. C++: `solidcore_ffi::units::blend_hcr`. Rust:
//! `solidcore::blend::hcr`.
//!
//! The C++ side builds a real `SmBrep` (the B-rep blocks of unit `topology.vertexuse`:
//! upstream test files, a wire edge, a face from a B-spline or plane surface, the wire
//! network), constructs an `SmHCR` on it from a view vector or a view placement, seeds
//! silhouette edges (the silhouette branch of `ComputeSilhouettes`' edge loop) and
//! back-facing edges, and runs one member. It exports the scene B-rep
//! (`GetSceneBrep()`, the copy transformed by the unported `SmBrep` copy constructor and
//! `Transform`) with its `CalculateBoundingBox` (unit `query.brep`, unported); the Rust
//! side imports that scene record by record (the `topology.vertexuse` importer), builds
//! its HCR with [`Hcr::from_scene`], seeds the same edges and runs the same member.
//!
//! Compared: the construction state (transforms, projection box, ray offset tolerance,
//! view vector, options), `IsSilhouetteEdge` / `IsBackFacingEdge`, the intervals and
//! point classes of `IntersectCurveWithSilhouettes` / `ComputeCurveVisibility`, the
//! visible and hidden curve pieces of `ComputeEdgeVisibility` (natural interval and 5
//! points each), `GetSilhouetteEdges` (sorted: upstream iterates a pointer hash) and
//! `IsPointVisible`. Status `-2` marks an upstream path that is undefined (UB-170, null
//! topology): the shim does not run it and the port returns `InvalidInput`. Ops whose
//! edge curves (or, for op 5, face surfaces) the importer cannot hold exactly (an
//! `SmEllipse` edge curve is imported as its B-spline part) report `inexact` on both
//! sides.
//!
//! `ComputeSilhouettes` and `ComputeGlobalVisibility` are not run: their first statement
//! is a `topology.brep_p1` query (every case would defer), and their C++ path reaches
//! UB-170 on open silhouette chains.

use super::topology_vertexuse::{gen_brep, import, path_of, push_curve, push_surface, sel, BrepIn};
use crate::{Case, CaseContext, HarnessError, Output, Rng, UnitDiff};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use solidcore::blend::hcr::{view_transform_for_vector, Hcr, HcrObject, HcrOptions};
use solidcore::curve::curve::p2::{BSPLINE_CURVE_TYPE, LINE_TYPE};
use solidcore::curve::curve::Curve;
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::extent3d::Extent3d;
use solidcore::intersect::curve_class::p1::CurveClassification;
use solidcore::intersect::point_class::PointClassification;
use solidcore::kernel::Kernel;
use solidcore::math::Vec3;
use solidcore::topology::{edge, Brep, EdgeId, VertexId};
use solidcore::SolidError;
use solidcore_ffi::units::{blend_hcr as ffi, topology_vertexuse as fixture};
use std::collections::BTreeSet;

/// The unit runner.
pub struct Runner;

fn harness(msg: impl Into<String>) -> HarnessError {
    HarnessError::new(msg)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum OpIn {
    Construct,
    IsSilhouette { edge: i64 },
    IsBackFacing { edge: i64 },
    Intersect { edge: i64, fa: f64, fb: f64 },
    CurveVisibility { edge: i64, fa: f64, fb: f64, start_hint: u8, end_hint: u8 },
    EdgeVisibility { edge: i64, vertex: i64, visible: bool, premark: bool },
    SilhouetteEdges,
    PointVisible { p: [f64; 3] },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Inputs {
    brep: BrepIn,
    /// View vector (mode 0) or the placement origin (mode 1).
    view: [f64; 3],
    /// Placement `Z`, `Xref` (mode 1).
    #[serde(default)]
    placement: Option<[[f64; 3]; 2]>,
    approx_tol: f64,
    ang_tol: f64,
    #[serde(default = "yes")]
    seams: bool,
    #[serde(default = "yes")]
    smooth: bool,
    /// Seeded silhouette edges (selectors); `[-1]`: every edge.
    #[serde(default)]
    sil: Vec<i64>,
    /// Seeded back-facing edges (selectors).
    #[serde(default)]
    back: Vec<i64>,
    #[serde(flatten)]
    op: OpIn,
}

const fn yes() -> bool {
    true
}

fn flag(b: bool) -> f64 {
    f64::from(u8::from(b))
}

/// The B-rep block (as unit `topology.vertexuse`, without gap initialisation).
fn brep_block(b: &BrepIn) -> Result<Vec<f64>, HarnessError> {
    let mut v = Vec::new();
    match b {
        BrepIn::File { rebuild, .. } => v.extend([0.0, flag(*rebuild)]),
        BrepIn::Wire { curve, t0, t1 } => {
            v.push(1.0);
            push_curve(&mut v, curve)?;
            v.extend([*t0, *t1]);
        }
        BrepIn::Face { surface, dom } => {
            v.push(2.0);
            push_surface(&mut v, surface)?;
            v.extend(dom);
        }
        BrepIn::Network => v.push(3.0),
    }
    v.extend([0.0, 0.0]);
    Ok(v)
}

/// Selectors resolved against the fixture's edge and vertex counts.
struct Resolved {
    sil: Vec<i64>,
    back: Vec<i64>,
    op: OpIn,
}

fn resolve(i: &Inputs, n_e: usize, n_v: usize) -> Option<Resolved> {
    if n_e == 0 && !matches!(i.op, OpIn::Construct | OpIn::SilhouetteEdges | OpIn::PointVisible { .. }) {
        return None;
    }
    let e = |s: i64| sel(s, n_e).max(0);
    let sil = if i.sil == [-1] {
        (0..n_e.min(64) as i64).collect()
    } else if n_e == 0 {
        Vec::new()
    } else {
        i.sil.iter().map(|&s| e(s)).collect()
    };
    let back = if n_e == 0 { Vec::new() } else { i.back.iter().map(|&s| e(s)).collect() };
    let op = match i.op.clone() {
        OpIn::IsSilhouette { edge } => OpIn::IsSilhouette { edge: e(edge) },
        OpIn::IsBackFacing { edge } => OpIn::IsBackFacing { edge: e(edge) },
        OpIn::Intersect { edge, fa, fb } => OpIn::Intersect { edge: e(edge), fa, fb },
        OpIn::CurveVisibility { edge, fa, fb, start_hint, end_hint } => {
            OpIn::CurveVisibility { edge: e(edge), fa, fb, start_hint, end_hint }
        }
        OpIn::EdgeVisibility { edge, vertex, visible, premark } => OpIn::EdgeVisibility {
            edge: e(edge),
            vertex: if vertex < 0 { -1 } else { sel(vertex, n_v) },
            visible,
            premark,
        },
        other => other,
    };
    Some(Resolved { sil, back, op })
}

fn op_args(op: &OpIn) -> (u8, Vec<f64>) {
    match *op {
        OpIn::Construct => (0, vec![]),
        OpIn::IsSilhouette { edge } => (1, vec![edge as f64]),
        OpIn::IsBackFacing { edge } => (2, vec![edge as f64]),
        OpIn::Intersect { edge, fa, fb } => (3, vec![edge as f64, fa, fb, 0.0, 0.0]),
        OpIn::CurveVisibility { edge, fa, fb, start_hint, end_hint } => {
            (4, vec![edge as f64, fa, fb, f64::from(start_hint), f64::from(end_hint)])
        }
        OpIn::EdgeVisibility { edge, vertex, visible, premark } => {
            (5, vec![edge as f64, vertex as f64, flag(visible), flag(premark)])
        }
        OpIn::SilhouetteEdges => (6, vec![]),
        OpIn::PointVisible { p } => (7, p.to_vec()),
    }
}

fn hcr_block(i: &Inputs, r: &Resolved) -> Vec<f64> {
    let mut v = Vec::new();
    match i.placement {
        None => v.push(0.0),
        Some(_) => v.push(1.0),
    }
    v.extend(i.view);
    if let Some([z, x]) = i.placement {
        v.extend(z);
        v.extend(x);
    }
    v.extend([i.approx_tol, i.ang_tol, flag(i.seams), flag(i.smooth)]);
    v.push(r.sil.len() as f64);
    v.extend(r.sil.iter().map(|&x| x as f64));
    v.push(r.back.len() as f64);
    v.extend(r.back.iter().map(|&x| x as f64));
    v
}

/// Builds the fixture, resolves the selectors, and runs shim op `op` (`None`: the op of
/// the case).
fn call(i: &Inputs, scene_only: bool) -> Result<Result<(Resolved, Vec<f64>), String>, HarnessError> {
    let path = path_of(&i.brep)?;
    let block = brep_block(&i.brep)?;
    let ex = match fixture::eval(fixture::Op::Export, &path, &block) {
        Ok(ex) => ex,
        Err(e) => return Ok(Err(format!("rejected: {e}"))),
    };
    let Some(c) = fixture::counts(&ex) else {
        return Ok(Err("build_failed".into()));
    };
    let Some(r) = resolve(i, c[7], c[9]) else {
        return Ok(Err("no_edge".into()));
    };
    let (op, args) = op_args(&r.op);
    let mut rest = hcr_block(i, &r);
    let op = if scene_only {
        8
    } else {
        rest.extend(args);
        op
    };
    match ffi::eval(op, &path, &block, &rest) {
        Ok(v) => Ok(Ok((r, v))),
        Err(e) => Ok(Err(format!("rejected: {e}"))),
    }
}

/// `(status, record, edge curve types, scene box, scene export)` of a shim output.
fn split(v: &[f64]) -> Result<(f64, &[f64], &[f64], Extent3d, &[f64]), HarnessError> {
    let short = || harness("short output");
    let status = *v.first().ok_or_else(short)?;
    let len = *v.get(1).ok_or_else(short)? as usize;
    let record = v.get(2..2 + len).ok_or_else(short)?;
    let n_e = *v.get(2 + len).ok_or_else(short)? as usize;
    let at = 3 + len;
    let types = v.get(at..at + n_e).ok_or_else(short)?;
    let b = v.get(at + n_e..at + n_e + 6).ok_or_else(short)?;
    let (lo, hi) = (Vec3::new(b[0], b[1], b[2]), Vec3::new(b[3], b[4], b[5]));
    let scene_box = if lo.x <= hi.x && lo.y <= hi.y && lo.z <= hi.z {
        Extent3d::from_points(lo, hi)
    } else if Extent3d::new().min() == lo && Extent3d::new().max() == hi {
        Extent3d::new()
    } else {
        return Err(harness("scene box neither ordered nor empty"));
    };
    Ok((status, record, types, scene_box, v.get(at + n_e + 6..).ok_or_else(short)?))
}

fn exact_curve(types: &[f64], e: i64) -> bool {
    types
        .get(e as usize)
        .is_some_and(|&t| t == BSPLINE_CURVE_TYPE as f64 || t == LINE_TYPE as f64)
}

/// Whether the op depends on geometry the import does not hold exactly.
fn inexact(r: &Resolved, types: &[f64], exact_faces: &[bool]) -> bool {
    let curves = |e: i64| !exact_curve(types, e) || r.sil.iter().any(|&s| !exact_curve(types, s));
    match r.op {
        OpIn::Intersect { edge, .. } | OpIn::CurveVisibility { edge, .. } => curves(edge),
        OpIn::EdgeVisibility { edge, .. } => curves(edge) || exact_faces.iter().any(|x| !x),
        _ => false,
    }
}

fn output(status: f64, record: Vec<f64>) -> Output {
    Output::new().values("status", vec![status]).values("record", record)
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "blend.hcr"
    }

    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let (r, v) = match call(&i, false)? {
            Ok(x) => x,
            Err(msg) => return Ok(Output::err(msg)),
        };
        let (status, record, types, _, ex) = split(&v)?;
        let (_, exact_faces) = import(ex)?;
        if inexact(&r, types, &exact_faces) {
            return Ok(Output::err("inexact"));
        }
        Ok(output(status, record.to_vec()))
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let (r, v) = match call(&i, true)? {
            Ok(x) => x,
            Err(msg) => return Ok(Output::err(msg)),
        };
        let (_, _, types, scene_box, ex) = split(&v)?;
        let (scene, exact_faces) = import(ex)?;
        if inexact(&r, types, &exact_faces) {
            return Ok(Output::err("inexact"));
        }
        match run_rust(&i, &r, scene, scene_box) {
            Ok((status, record)) => Ok(output(status, record)),
            Err(e) => Ok(Output::unported(&e).unwrap_or_else(|| match e {
                SolidError::InvalidInput { operation: "blend::hcr", .. } => output(-2.0, vec![]),
                SolidError::Upstream { code, .. } => output(f64::from(code), vec![]),
                other => Output::err(format!("rust_error: {other}")),
            })),
        }
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        Some(gen_case(rng))
    }

    fn cases(&self, ctx: &CaseContext) -> Result<Vec<Case>, HarnessError> {
        let ctx = CaseContext {
            root: ctx.root.clone(),
            seeds: ctx.seeds.max(300),
        };
        crate::runner::default_cases(self, &ctx)
    }
}

fn eid(i: i64) -> EdgeId {
    EdgeId::from_index(i as u32)
}

fn put_axis(o: &mut Vec<f64>, a: &Axis2Placement) {
    for v in [a.origin(), a.x_axis(), a.y_axis(), a.z_axis()] {
        o.extend([v.x, v.y, v.z]);
    }
}

fn put_box(o: &mut Vec<f64>, b: &Extent3d) {
    let (lo, hi) = (b.min(), b.max());
    o.extend([lo.x, lo.y, lo.z, hi.x, hi.y, hi.z]);
}

fn put_pc(o: &mut Vec<f64>, pc: &PointClassification<HcrObject>) {
    o.push(f64::from(pc.point_class().code()));
    let (kind, idx) = match pc.object() {
        None => (0.0, -1.0),
        Some(HcrObject::Vertex(v)) => (1.0, f64::from(v.index())),
        Some(HcrObject::Edge(e)) => (2.0, f64::from(e.index())),
        Some(HcrObject::Region(r)) => (3.0, f64::from(r.index())),
        Some(HcrObject::Curve(_)) => (4.0, -1.0),
    };
    o.extend([kind, idx, pc.pre_snap_param(), pc.t_param()]);
}

fn put_cc(o: &mut Vec<f64>, cc: &CurveClassification<HcrObject>) {
    o.push(cc.size() as f64);
    for c in cc.intervals() {
        o.extend([c.interval.min(), c.interval.max()]);
        put_pc(o, &c.start);
        put_pc(o, &c.mid);
        put_pc(o, &c.end);
    }
}

fn put_piece(o: &mut Vec<f64>, c: &Curve) {
    let n = c.natural_interval();
    o.extend([n.min(), n.max()]);
    for i in 0..5 {
        match c.evaluate_point(n.evaluate(f64::from(i) / 4.0)) {
            Ok(p) => o.extend([0.0, p.x, p.y, p.z]),
            Err(e) => o.extend([
                match e {
                    SolidError::Upstream { code, .. } => f64::from(code),
                    _ => -1.0,
                },
                0.0,
                0.0,
                0.0,
            ]),
        }
    }
}

fn run_rust(i: &Inputs, r: &Resolved, scene: Brep, scene_box: Extent3d) -> Result<(f64, Vec<f64>), SolidError> {
    let k = Kernel::default();
    let transform = match i.placement {
        None => view_transform_for_vector(Vec3::new(i.view[0], i.view[1], i.view[2]))?,
        Some([z, x]) => {
            let mut a = Axis2Placement::default();
            let _ = a.set_step_canonical(
                Vec3::new(i.view[0], i.view[1], i.view[2]),
                Vec3::new(z[0], z[1], z[2]),
                Vec3::new(x[0], x[1], x[2]),
            );
            a
        }
    };
    let opts = HcrOptions {
        approx_tol_3d: i.approx_tol,
        ang_tol_rad: i.ang_tol,
        display_seams: i.seams,
        display_smooth_edges: i.smooth,
    };
    let mut h = Hcr::from_scene(&k, scene, &transform, scene_box, opts);
    for &e in &r.sil {
        h.add_silhouette_edge(&k, eid(e))?;
    }
    for &e in &r.back {
        h.add_back_facing_edge(eid(e));
    }
    let mut rec = Vec::new();
    let mut status = 0.0;
    match r.op {
        OpIn::Construct => {
            put_axis(&mut rec, h.view_transform());
            put_axis(&mut rec, h.inv_view_transform());
            put_box(&mut rec, h.proj_bbox());
            rec.push(h.ray_offset_tolerance());
            let v = h.view_vector();
            rec.extend([v.x, v.y, v.z, h.approx_tol_3d(), h.ang_tol_rad()]);
            rec.extend([flag(h.display_seams()), flag(h.display_smooth_edges())]);
        }
        OpIn::IsSilhouette { edge } => rec.push(flag(h.is_silhouette_edge(&k, eid(edge))?)),
        OpIn::IsBackFacing { edge } => rec.push(flag(h.is_back_facing_edge(&k, eid(edge))?)),
        OpIn::Intersect { edge, fa, fb } | OpIn::CurveVisibility { edge, fa, fb, .. } => {
            let e = eid(edge);
            let ivl = edge::interval(h.scene_brep(), e)?;
            let ivl = Extent1d::from_raw(ivl.evaluate(fa), ivl.evaluate(fb));
            let tol = edge::tolerance(h.scene_brep(), e)?;
            let mut cc = CurveClassification::with_curve_interval(&k, HcrObject::Curve(e), ivl, tol);
            let res = match r.op {
                OpIn::CurveVisibility { start_hint, end_hint, .. } => {
                    h.compute_curve_visibility(&k, &mut cc, u64::from(start_hint), u64::from(end_hint))
                }
                _ => h.intersect_curve_with_silhouettes(&k, &mut cc),
            };
            match res {
                Ok(()) => {}
                Err(e) if e.is_unported() || matches!(e, SolidError::InvalidInput { .. }) => return Err(e),
                Err(SolidError::Upstream { code, .. }) => status = f64::from(code),
                Err(e) => return Err(e),
            }
            put_cc(&mut rec, &cc);
        }
        OpIn::EdgeVisibility { edge, vertex, visible, premark } => {
            let mut marked = BTreeSet::new();
            if premark {
                marked.insert(eid(edge));
            }
            let v = (vertex >= 0).then(|| VertexId::from_index(vertex as u32));
            match h.compute_edge_visibility(&k, eid(edge), v, visible, &mut marked) {
                Ok(()) => {}
                Err(SolidError::Upstream { code, .. }) => status = f64::from(code),
                Err(e) => return Err(e),
            }
            rec.push(h.visible_curves().len() as f64);
            for c in h.visible_curves() {
                put_piece(&mut rec, c);
            }
            rec.push(h.invisible_curves().len() as f64);
            for c in h.invisible_curves() {
                put_piece(&mut rec, c);
            }
        }
        OpIn::SilhouetteEdges => {
            let edges = h.silhouette_edges();
            rec.push(edges.len() as f64);
            rec.extend(edges.iter().map(|e| f64::from(e.index())));
        }
        OpIn::PointVisible { p } => rec.push(flag(h.is_point_visible(&k, Vec3::new(p[0], p[1], p[2]))?)),
    }
    Ok((status, rec))
}

// ---- generation -------------------------------------------------------------------------

fn gen_case(rng: &mut Rng) -> Value {
    let brep = gen_brep(rng);
    let scale = *rng.choose(&[1e-4, 1.0, 1.0, 1.0, 1e4]).unwrap_or(&1.0);
    let view = match rng.below(6) {
        0 => [0.0, 0.0, 1.0],
        1 => [0.0, 0.0, -1.0],
        2 => [1.0, 0.0, 0.0],
        3 => [0.0, 1e-13, 0.0],
        _ => {
            let p = rng.point3(-1.0, 1.0);
            [p[0] * scale, p[1] * scale, p[2] * scale + 1e-3]
        }
    };
    let placement = rng.chance(0.25).then(|| {
        let z = rng.point3(-1.0, 1.0);
        let x = if rng.chance(0.2) { z } else { rng.point3(-1.0, 1.0) };
        [[z[0], z[1], z[2] + 0.01], [x[0], x[1], x[2]]]
    });
    let sil: Vec<i64> = match rng.below(4) {
        0 => vec![],
        1 => vec![-1],
        _ => (0..rng.below(4)).map(|_| rng.below(1 << 20) as i64).collect(),
    };
    let back: Vec<i64> = (0..rng.below(3)).map(|_| rng.below(1 << 20) as i64).collect();
    let edge = rng.below(1 << 20) as i64;
    let (fa, fb) = if rng.chance(0.5) {
        (0.0, 1.0)
    } else {
        let a = rng.uniform(0.0, 0.6);
        (a, a + rng.uniform(0.1, 0.4))
    };
    let hint = |rng: &mut Rng| *rng.choose(&[0u8, 1, 1, 2]).unwrap_or(&1);
    let op = match rng.below(16) {
        0 | 1 => json!({"op": "construct"}),
        2 | 3 => json!({"op": "is_silhouette", "edge": edge}),
        4 | 5 => json!({"op": "is_back_facing", "edge": edge}),
        6..=8 => json!({"op": "intersect", "edge": edge, "fa": fa, "fb": fb}),
        9 | 10 => json!({"op": "curve_visibility", "edge": edge, "fa": fa, "fb": fb,
                          "start_hint": hint(rng), "end_hint": hint(rng)}),
        11..=13 => json!({"op": "edge_visibility", "edge": edge,
                          "vertex": if rng.chance(0.5) { -1 } else { rng.below(1 << 20) as i64 },
                          "visible": rng.chance(0.5), "premark": rng.chance(0.2)}),
        14 => json!({"op": "silhouette_edges"}),
        _ => {
            let p = rng.point3(-2.0, 2.0);
            json!({"op": "point_visible", "p": p})
        }
    };
    let mut v = json!({
        "brep": brep,
        "view": view,
        "approx_tol": *rng.choose(&[1e-6, 1e-3, 0.1]).unwrap_or(&1e-3),
        "ang_tol": *rng.choose(&[0.01, 0.2]).unwrap_or(&0.2),
        "seams": rng.chance(0.7),
        "smooth": rng.chance(0.7),
        "sil": sil,
        "back": back,
    });
    if let Some(p) = placement {
        v["placement"] = json!(p);
    }
    if let (Some(o), Some(map)) = (op.as_object(), v.as_object_mut()) {
        for (key, val) in o {
            map.insert(key.clone(), val.clone());
        }
    }
    v
}
