//! Runner for unit `topology.brep_p7`. C++: `solidcore_ffi::units::topology_brep_p7`.
//! Rust: `solidcore::topology::brep::p7`.
//!
//! **Ops.** `read` (`SmBrep::ReadFromFile` of an upstream `.smb` file, as read or with
//! bytes truncated / overwritten), `make_face_with_curves` (planes and B-spline surfaces
//! with loops of B-spline curves, 3D only, UV only or both, and the input errors),
//! `make_face` (a null surface, or a square wire loop built in the shim), `glue`
//! (`GlueEdgesGeneral` of edge pairs of a B-rep of unit `topology.vertexuse`), `copy`
//! (`CopyFaces` of the first faces into a new B-rep) and `merge_self`
//! (`MergeSelfIntersectionCurves` without curves). Both sides produce one record
//! (`out`) laid out as in the ffi module; statuses are `SmStatus` values (a Rust
//! `InvalidInput` is `SM_ERR_INVALID_INPUT`, `1007`).
//!
//! **Guards.** A `read` stream goes through the stream guard of unit `topology.brep_data`
//! on both sides (a rejected stream is an `err("guard")` case). A modified stream whose
//! C++ `SmBrepData` read succeeds is not run in C++: `MakeTopologyFromData` uses its
//! cross-reference indices unchecked and crashes on some of them
//! (`plan/upstream-issues.md` UB-164). Both sides decide this with the C++ reader; the
//! C++ side reports `-3` and the Rust side defers to `cpp-guard:unvalidated-topology`.
//! **Revisit when `topology.brep_p10` lands** (e.g. run the C++ read out of process): a
//! guard-free sweep crashed in 13 of 15 files, while most such streams are safe. The other
//! shim guards answer `-2` (upstream would dereference null); both sides would only echo
//! that constant, so the Rust side defers them: `cpp-guard:ub-165`
//! (`MakeFaceWithCurves` with a curve on surface and no UV curves), `cpp-guard:null-curve`
//! (`GlueEdgesGeneral` with `SM_OT_SAME`, equal intervals and one edge without curve),
//! `cpp-guard:null-surface`
//! (`MergeSelfIntersectionCurves` on a face without surface).
//!
//! **Deferred.** Every read that gets past `SmBrepData::ReadFromFile` (it now runs
//! `MakeTopologyFromData` of `topology.brep_p10` on a B-rep with the infinite region
//! and defers later: `sm_CheckAndFixUVTrimCurves`, analytic copies, `HealBrep`'s
//! callees), every face made
//! (`MakeShellVertex`, `MakeFaceTopology`, `GetInfiniteRegion`), every glue past its
//! early returns, every copy (`SmBrep::SetTolerance`). The comparisons are the read
//! errors of damaged files, the input errors of the face makers, the null-surface status
//! of `MakeFace`, the glue no-ops and `MergeSelfIntersectionCurves` without curves.
//! Region nesting borrows the live B-rep through the shell query provider. Native
//! owner/list metadata uses the same build and traversal registry as each input;
//! topology limitations are checked before traversal, and geometry limitations
//! only when a query consumes that geometry.
//! Glue imports the native notification records, including all five stored marks,
//! flags and opaque user values, before invoking its real merge notification.
//! Attribute storage and callbacks have no decoder here: their measured presence
//! is a conservative fixture-level harness limit after the null/same-edge no-ops.
//! Native marks are restored exactly and never trigger that representation guard.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::rc::Rc;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use solidcore::curve::curve::Curve as SmCurve;
use solidcore::curve::line::Line;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::object::ObjectRef;
use solidcore::geom::tol::are_same;
use solidcore::geom::vector2d::Point2;
use solidcore::geom::vector3d::Vec3;
use solidcore::kernel::Kernel;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::b_spline_surface::{BSplineSurface, BSplineSurfaceForm};
use solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType;
use solidcore::nurbs::crv_basic::Curve as NlCurve;
use solidcore::nurbs::math::CPoint;
use solidcore::nurbs::srf_basic::Surface as NlSurface;
use solidcore::surface::plane::Plane;
use solidcore::surface::surface::Surface;
use solidcore::topology::brep::{p10, p7};
use solidcore::topology::brep_cache::{arena_identity, CacheObject};
use solidcore::topology::database_io::FileType;
use solidcore::topology::edge as ed;
use solidcore::topology::topology::Topology;
use solidcore::topology::{Arena, Brep, BrepId, EdgeId, FaceId, Region};
use solidcore::{EntityRef, SolidError};
use solidcore_ffi::units::topology_brep_data as bd;
use solidcore_ffi::units::topology_brep_p7::{self as ffi, Op};
use solidcore_ffi::units::topology_vertexuse as fixture;
use solidcore_ffi::FfiError;

use super::topology_shell::RayQueryContext;
use super::topology_vertexuse::{self as shared, apply_query_metadata, BrepIn, QueryMetadata};
use crate::case::CaseContext;
use crate::runner::default_cases;
use crate::{Case, HarnessError, Output, Rng, UnitDiff};

/// The difftest runner.
pub struct Runner;

fn harness(msg: impl Into<String>) -> HarnessError {
    HarnessError::new(msg)
}

/// `SM_ERR_INVALID_INPUT`.
const SM_ERR_INVALID_INPUT: f64 = 1007.0;
/// "Not run": the shim keeps the call away from upstream (`-3`).
const NOT_RUN: f64 = -3.0;
/// `NOW`: the weight of a non-rational control point.
const NOW: f64 = -f64::MAX;
/// `NOZ`: the z of a 2D control point.
const NOZ: f64 = f64::MAX;

/// The `.smb` files of the cases (the readable ASCII files of unit `topology.vertexuse`).
const FILES: [&str; 15] = [
    "Mouse.smb",
    "Endcap.smb",
    "box_with_large_fillet.smb",
    "cone.smb",
    "disk.smb",
    "ellips_ts2.smb",
    "plane_circ.smb",
    "plane_square.smb",
    "sphere.smb",
    "surf_extrude.smb",
    "surf_rev.smb",
    "surf_sphere150.smb",
    "torus.smb",
    "toruswvy.smb",
    "vase_rev.smb",
];

// ---- inputs -----------------------------------------------------------------------------

/// A B-spline curve (`w = NOW` non-rational).
#[derive(Debug, Clone, Serialize, Deserialize)]
struct CrvIn {
    dim: u64,
    cpts: Vec<[f64; 4]>,
    p: usize,
    u: Vec<f64>,
}

/// A surface: a plane or a B-spline surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum SrfIn {
    Plane {
        o: [f64; 3],
        n: [f64; 3],
    },
    Bspline {
        net: Vec<Vec<[f64; 4]>>,
        p: usize,
        q: usize,
        u: Vec<f64>,
        v: Vec<f64>,
    },
}

/// One curve of a face: its 3D and/or UV curve and its orientation code (`1` SAME, `2`
/// OPPOSITE).
#[derive(Debug, Clone, Serialize, Deserialize)]
struct FaceCurveIn {
    #[serde(default)]
    c3: Option<CrvIn>,
    #[serde(default)]
    uv: Option<CrvIn>,
    orient: u8,
}

/// A change made to a file's bytes before reading.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum Damage {
    /// Keep the first `n` bytes.
    Truncate(usize),
    /// Overwrite the byte at `at` (modulo the length) with `byte`.
    Overwrite { at: usize, byte: u8 },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum Inputs {
    Read {
        file: String,
        #[serde(default)]
        damage: Option<Damage>,
        #[serde(default)]
        rebuild: bool,
    },
    MakeFaceWithCurves {
        surface: SrfIn,
        dom: [f64; 4],
        orient: u8,
        region_null: bool,
        /// `pbOptCrvOnSurf`: 0 null, 1 FALSE, 2 TRUE.
        crv_on_surf: u8,
        loops: Vec<usize>,
        has3d: bool,
        has_uv: bool,
        curves: Vec<FaceCurveIn>,
        #[serde(default)]
        points: Vec<[f64; 3]>,
    },
    MakeFace {
        kind: u8,
        side: f64,
        orient: u8,
        region_null: bool,
    },
    Glue {
        brep: BrepIn,
        e1: i64,
        e2: i64,
        orient: u8,
        gap: f64,
        nesting: bool,
    },
    Copy {
        brep: BrepIn,
        n: u64,
        merge: bool,
    },
    MergeSelf {
        brep: BrepIn,
        face: i64,
        tol: f64,
    },
}

fn flag(b: bool) -> f64 {
    if b {
        1.0
    } else {
        0.0
    }
}

fn push_curve(v: &mut Vec<f64>, c: &CrvIn) {
    let n = c.cpts.len().saturating_sub(1);
    let m = c.u.len().saturating_sub(1);
    v.extend([c.dim as f64, n as f64, c.p as f64, m as f64]);
    for p in &c.cpts {
        v.extend(p);
    }
    v.extend(&c.u);
}

fn push_surface(v: &mut Vec<f64>, s: &SrfIn) {
    match s {
        SrfIn::Plane { o, n } => {
            v.push(1.0);
            v.extend(o);
            v.extend(n);
        }
        SrfIn::Bspline {
            net,
            p,
            q,
            u,
            v: vk,
        } => {
            let n = net.len().saturating_sub(1);
            let m = net.first().map_or(0, |r| r.len().saturating_sub(1));
            v.extend([0.0, 0.0, 0.0]);
            v.extend([n as f64, m as f64, *p as f64, *q as f64]);
            v.extend([
                (u.len().saturating_sub(1)) as f64,
                (vk.len().saturating_sub(1)) as f64,
            ]);
            for row in net {
                for pt in row {
                    v.extend(pt);
                }
            }
            v.extend(u);
            v.extend(vk);
        }
    }
}

/// The B-rep block (`validate_only` and `init_gaps` off).
fn brep_block(b: &BrepIn) -> Result<Vec<f64>, HarnessError> {
    let mut out = Vec::new();
    match b {
        BrepIn::File { rebuild, .. } => out.extend([0.0, flag(*rebuild)]),
        BrepIn::Wire { curve, t0, t1 } => {
            out.push(1.0);
            shared::push_curve(&mut out, curve)?;
            out.extend([*t0, *t1]);
        }
        BrepIn::Face { surface, dom } => {
            out.push(2.0);
            shared::push_surface(&mut out, surface)?;
            out.extend(dom);
        }
        BrepIn::Network => out.push(3.0),
    }
    out.extend([0.0, 0.0]);
    Ok(out)
}

/// Whether `bytes` are exactly the test file `file` (a damage that changed nothing): the
/// FFI runs the native read without the B-rep data check on exactly the shipped files.
fn is_test_file(file: &str, bytes: &[u8]) -> Result<bool, HarnessError> {
    Ok(read_bytes(file, None)? == bytes)
}

/// The bytes of a read case.
fn read_bytes(file: &str, damage: Option<Damage>) -> Result<Vec<u8>, HarnessError> {
    let p = crate::fixtures::upstream_test_file(&crate::case::repo_root(), file)?;
    let mut bytes = std::fs::read(&p).map_err(|e| harness(format!("{}: {e}", p.display())))?;
    match damage {
        None => {}
        Some(Damage::Truncate(n)) => bytes.truncate(n),
        Some(Damage::Overwrite { at, byte }) => {
            let len = bytes.len();
            if len > 0 {
                if let Some(b) = bytes.get_mut(at % len) {
                    *b = byte;
                }
            }
        }
    }
    Ok(bytes)
}

/// Whether the C++ `SmBrepData::ReadFromFile` (through the `topology.brep_data` shim)
/// reads the bytes successfully: the FFI's test for not running a modified stream.
fn cpp_brep_data_reads(bytes: &[u8]) -> Result<bool, Fail> {
    let req = bd::Request {
        mode: bd::Mode::ReadFile,
        // every read case is ASCII: the encoder sends `binary = 0` to the FFI
        binary: false,
        version: 0,
        dim: 0,
        persist: false,
        source: 0,
        pre: 0,
        data: bytes.to_vec(),
        upstream_file: false,
    };
    let items = bd::run(&req).map_err(|e| hfail(format!("c++ brep data: {e}")))?;
    Ok(items.first() == Some(&bd::Item::Int(0)))
}

/// The `topology.brep_data` stream guard of the bytes: `false` when it rejects them.
fn guard_accepts(bytes: &[u8], upstream_file: bool) -> Result<bool, HarnessError> {
    let req = bd::Request {
        mode: bd::Mode::ReadFile,
        binary: false,
        version: 0,
        dim: 0,
        persist: false,
        source: 0,
        pre: 0,
        data: bytes.to_vec(),
        upstream_file,
    };
    match bd::run(&req) {
        Ok(_) => Ok(true),
        Err(FfiError::InvalidArgument(_)) => Ok(false),
        Err(e) => Err(harness(format!("c++ guard: {e}"))),
    }
}

/// A cached export with its key (the `.smb` path and the block's bits).
type CachedExport = (String, Vec<u64>, Rc<Vec<f64>>);

thread_local! {
    /// The last few exports (both sides of a case export the same B-rep).
    static EXPORTS: RefCell<Vec<CachedExport>> = const { RefCell::new(Vec::new()) };
}

fn export_of(brep: &BrepIn) -> Result<Rc<Vec<f64>>, HarnessError> {
    let path = shared::path_of(brep)?;
    let block = brep_block(brep)?;
    let bits: Vec<u64> = block.iter().map(|x| x.to_bits()).collect();
    if let Some(e) = EXPORTS.with(|c| {
        c.borrow()
            .iter()
            .find(|(p, b, _)| *p == path && *b == bits)
            .map(|(_, _, e)| Rc::clone(e))
    }) {
        return Ok(e);
    }
    let e = Rc::new(
        fixture::eval(fixture::Op::Export, &path, &block)
            .map_err(|e| harness(format!("c++ export: {e}")))?,
    );
    EXPORTS.with(|c| {
        let mut c = c.borrow_mut();
        if c.len() >= 4 {
            c.remove(0);
        }
        c.push((path, bits, Rc::clone(&e)));
    });
    Ok(e)
}

/// The inputs with selectors resolved against the B-rep, and its export; `Err(output)`
/// when the case ends early (the B-rep did not build, or has nothing to select).
#[allow(clippy::type_complexity)]
fn resolve(i: &Inputs) -> Result<Result<(Inputs, Option<Rc<Vec<f64>>>), Output>, HarnessError> {
    let (Inputs::Glue { brep, .. } | Inputs::Copy { brep, .. } | Inputs::MergeSelf { brep, .. }) =
        i
    else {
        return Ok(Ok((i.clone(), None)));
    };
    let ex = export_of(brep)?;
    let Some(c) = fixture::counts(&ex) else {
        return Ok(Err(Output::new().values("out", ex.to_vec())));
    };
    let [_, _, _, n_f, _, _, _, n_e, ..] = c;
    let mut r = i.clone();
    match &mut r {
        Inputs::Glue { e1, e2, .. } => {
            *e1 = shared::sel(*e1, n_e);
            *e2 = shared::sel(*e2, n_e);
        }
        Inputs::MergeSelf { face, .. } => {
            if n_f == 0 {
                return Ok(Err(Output::err("no_face")));
            }
            *face = shared::sel((*face).max(0), n_f);
        }
        _ => {}
    }
    Ok(Ok((r, Some(ex))))
}

// ---- C++ side ---------------------------------------------------------------------------

fn encode(i: &Inputs) -> Result<(Op, String, Vec<f64>), HarnessError> {
    let mut v = Vec::new();
    let mut path = String::new();
    let op = match i {
        Inputs::Read {
            file,
            damage,
            rebuild,
        } => {
            let bytes = read_bytes(file, *damage)?;
            v.push(bytes.len() as f64);
            v.extend(bytes.iter().map(|&b| f64::from(b)));
            v.extend([0.0, flag(*rebuild)]);
            Op::Read
        }
        Inputs::MakeFaceWithCurves {
            surface,
            dom,
            orient,
            region_null,
            crv_on_surf,
            loops,
            has3d,
            has_uv,
            curves,
            points,
        } => {
            push_surface(&mut v, surface);
            v.extend(dom);
            v.extend([
                f64::from(*orient),
                flag(*region_null),
                f64::from(*crv_on_surf),
            ]);
            v.push(loops.len() as f64);
            v.extend(loops.iter().map(|&n| n as f64));
            v.extend([flag(*has3d), flag(*has_uv), curves.len() as f64]);
            for c in curves {
                if *has3d {
                    push_curve(
                        &mut v,
                        c.c3.as_ref().ok_or_else(|| harness("missing 3d curve"))?,
                    );
                }
                if *has_uv {
                    push_curve(
                        &mut v,
                        c.uv.as_ref().ok_or_else(|| harness("missing uv curve"))?,
                    );
                }
                v.push(f64::from(c.orient));
            }
            v.push(points.len() as f64);
            for p in points {
                v.extend(p);
            }
            Op::MakeFaceWithCurves
        }
        Inputs::MakeFace {
            kind,
            side,
            orient,
            region_null,
        } => {
            v.extend([
                f64::from(*kind),
                *side,
                f64::from(*orient),
                flag(*region_null),
            ]);
            Op::MakeFace
        }
        Inputs::Glue {
            brep,
            e1,
            e2,
            orient,
            gap,
            nesting,
        } => {
            v.extend([
                *e1 as f64,
                *e2 as f64,
                f64::from(*orient),
                *gap,
                flag(*nesting),
            ]);
            v.extend(brep_block(brep)?);
            path = shared::path_of(brep)?;
            Op::Glue
        }
        Inputs::Copy { brep, n, merge } => {
            v.extend([*n as f64, flag(*merge)]);
            v.extend(brep_block(brep)?);
            path = shared::path_of(brep)?;
            Op::Copy
        }
        Inputs::MergeSelf { brep, face, tol } => {
            v.extend([*face as f64, *tol]);
            v.extend(brep_block(brep)?);
            path = shared::path_of(brep)?;
            Op::MergeSelf
        }
    };
    Ok((op, path, v))
}

fn run_cpp(i: &Inputs) -> Result<Output, HarnessError> {
    if let Inputs::Read { file, damage, .. } = i {
        let bytes = read_bytes(file, *damage)?;
        if !guard_accepts(&bytes, is_test_file(file, &bytes)?)? {
            return Ok(Output::err("guard"));
        }
    }
    let (r, _) = match resolve(i)? {
        Ok(x) => x,
        Err(o) => return Ok(o),
    };
    let (op, path, input) = encode(&r)?;
    let out = ffi::eval(op, &path, &input).map_err(|e| harness(format!("c++: {e}")))?;
    Ok(Output::new().values("out", out))
}

// ---- Rust side --------------------------------------------------------------------------

/// What ends a Rust case early: a deferral or a harness error.
enum Fail {
    Deferred(Box<Output>),
    Harness(HarnessError),
}

impl From<HarnessError> for Fail {
    fn from(e: HarnessError) -> Self {
        Self::Harness(e)
    }
}

fn hfail(e: impl core::fmt::Display) -> Fail {
    Fail::Harness(harness(e.to_string()))
}

/// A case the shim answers with `-2` (upstream would dereference null): both sides would
/// only echo that constant, so it is deferred under the guard's `label`, never compared.
fn cpp_null_deref(label: &str, symbol: &str) -> Fail {
    Fail::Deferred(Box::new(Output::deferred(label, symbol)))
}

/// The `SmStatus` of a Rust result (`0` success, `SM_ERR_INVALID_INPUT` for
/// `InvalidInput`); a deferral or another error ends the case.
fn status<T>(r: &Result<T, SolidError>) -> Result<f64, Fail> {
    match r {
        Ok(_) => Ok(0.0),
        Err(SolidError::Upstream { code, .. }) => Ok(f64::from(*code)),
        Err(SolidError::InvalidInput { .. }) => Ok(SM_ERR_INVALID_INPUT),
        Err(e) => Err(match Output::unported(e) {
            Some(o) => Fail::Deferred(Box::new(o)),
            None => hfail(format!("rust: {e}")),
        }),
    }
}

fn rust_curve(c: &CrvIn) -> Result<SmCurve, HarnessError> {
    let pw = c
        .cpts
        .iter()
        .map(|p| CPoint::new(p[0], p[1], p[2], p[3]))
        .collect();
    let cur = NlCurve::new(pw, c.p, c.u.clone()).map_err(|e| harness(format!("curve: {e}")))?;
    Ok(SmCurve::from(BSplineCurve::from_nurb(cur, c.dim, false)))
}

fn rust_surface(s: &SrfIn) -> Result<Surface, HarnessError> {
    match s {
        SrfIn::Plane { o, n } => Ok(Surface::from(Plane::from_normal(
            Vec3::new(o[0], o[1], o[2]),
            Vec3::new(n[0], n[1], n[2]),
        ))),
        SrfIn::Bspline { net, p, q, u, v } => {
            let net = net
                .iter()
                .map(|r| {
                    r.iter()
                        .map(|x| CPoint::new(x[0], x[1], x[2], x[3]))
                        .collect()
                })
                .collect();
            let nl = NlSurface::new(net, *p, *q, u.clone(), v.clone())
                .map_err(|e| harness(format!("surface: {e}")))?;
            let mut s = BSplineSurface::from_nurb(nl, false);
            s.set_b_spline_surface_form(
                BSplineSurfaceForm::from_code(0).ok_or_else(|| harness("form"))?,
            );
            Ok(Surface::from(s))
        }
    }
}

/// `SmOrientType` codes `SM_OT_UNKNOWN` (0), `SM_OT_SAME` (1), `SM_OT_OPPOSITE` (2).
fn orient_of(code: u8) -> Result<OrientType, Fail> {
    match code {
        0 => Ok(OrientType::Unknown),
        1 => Ok(OrientType::Same),
        2 => Ok(OrientType::Opposite),
        _ => Err(hfail("SmOrientType code")),
    }
}

fn key_of(e: EntityRef) -> ObjectRef {
    let (type_code, key) = match e {
        EntityRef::Face(f) => (
            solidcore::topology::face::p2::TYPE,
            1 + u64::from(f.index()),
        ),
        EntityRef::Edge(x) => (ed::TYPE, 1 << 40 | u64::from(x.index())),
        EntityRef::Edgeuse(x) => (16005, 2 << 40 | u64::from(x.index())),
        _ => (0, 3 << 40),
    };
    ObjectRef {
        type_code,
        key: ObjectKey(key),
    }
}

/// A distinct identity for every topology object `ReadFromFile` builds (geometry copies
/// are keyed from `1 << 50` up).
fn read_key_of(e: EntityRef) -> ObjectRef {
    let (type_code, kind, i) = match e {
        EntityRef::Brep(x) => (16011, 0, x.index()),
        EntityRef::Region(x) => (16012, 1, x.index()),
        EntityRef::Shell(x) => (16013, 2, x.index()),
        EntityRef::Faceuse(x) => (16022, 3, x.index()),
        EntityRef::Face(x) => (16002, 4, x.index()),
        EntityRef::Loopuse(x) => (16004, 5, x.index()),
        EntityRef::Loop(x) => (16003, 6, x.index()),
        EntityRef::Edgeuse(x) => (16005, 7, x.index()),
        EntityRef::Edge(x) => (16016, 8, x.index()),
        EntityRef::Vertexuse(x) => (16006, 9, x.index()),
        EntityRef::Vertex(x) => (16021, 10, x.index()),
    };
    ObjectRef {
        type_code,
        key: ObjectKey((kind + 1) << 40 | u64::from(i)),
    }
}

fn read_query_key(object: CacheObject) -> ObjectRef {
    match object {
        CacheObject::Edge(edge) | CacheObject::Curve(edge) => read_key_of(EntityRef::Edge(edge)),
        CacheObject::Face(face) | CacheObject::Surface(face) => read_key_of(EntityRef::Face(face)),
        CacheObject::Vertex(vertex) => read_key_of(EntityRef::Vertex(vertex)),
        CacheObject::Brep => ObjectRef {
            type_code: 16011,
            key: ObjectKey(1u64 << 40),
        },
        other => arena_identity(other),
    }
}

fn this_face(f: FaceId) -> ObjectRef {
    key_of(EntityRef::Face(f))
}

fn query_key(object: CacheObject) -> ObjectRef {
    match object {
        CacheObject::Edge(edge) | CacheObject::Curve(edge) => key_of(EntityRef::Edge(edge)),
        CacheObject::Face(face) | CacheObject::Surface(face) => key_of(EntityRef::Face(face)),
        _ => arena_identity(object),
    }
}

/// Fixture identities include the class because notification records are keyed
/// only by `ObjectKey`. For example, edge 0 and vertex 0 must remain distinct.
fn glue_object_ref(type_code: i64, index: u64) -> ObjectRef {
    ObjectRef {
        type_code,
        key: ObjectKey(((type_code as u64) << 32) | (index + 1)),
    }
}

fn glue_key(entity: EntityRef) -> ObjectRef {
    let (type_code, index) = match entity {
        EntityRef::Region(id) => (16012, id.index()),
        EntityRef::Shell(id) => (16013, id.index()),
        EntityRef::Faceuse(id) => (16022, id.index()),
        EntityRef::Face(id) => (16002, id.index()),
        EntityRef::Loopuse(id) => (16004, id.index()),
        EntityRef::Loop(id) => (16003, id.index()),
        EntityRef::Edgeuse(id) => (16005, id.index()),
        EntityRef::Edge(id) => (16016, id.index()),
        EntityRef::Vertexuse(id) => (16006, id.index()),
        EntityRef::Vertex(id) => (16021, id.index()),
        // Every call contains exactly the B-rep imported for this fixture.
        EntityRef::Brep(_) => (16011, 0),
    };
    glue_object_ref(type_code, u64::from(index))
}

fn glue_query_key(object: CacheObject) -> ObjectRef {
    let identity = arena_identity(object);
    glue_object_ref(identity.type_code, identity.key.0)
}

struct GlueMetadataReader<'a> {
    raw: &'a [f64],
    at: usize,
}

impl<'a> GlueMetadataReader<'a> {
    fn word(&mut self) -> Result<u32, HarnessError> {
        let value = *self
            .raw
            .get(self.at)
            .ok_or_else(|| harness("glue metadata: truncated integer"))?;
        self.at += 1;
        if !value.is_finite()
            || value.fract() != 0.0
            || !(0.0..=f64::from(u32::MAX)).contains(&value)
        {
            return Err(harness("glue metadata: invalid integer"));
        }
        Ok(value as u32)
    }

    fn count(&mut self) -> Result<usize, HarnessError> {
        usize::try_from(self.word()?).map_err(|_| harness("glue metadata: count overflow"))
    }

    fn wide(&mut self) -> Result<u64, HarnessError> {
        let lo = u64::from(self.word()?);
        let hi = u64::from(self.word()?);
        Ok(lo | (hi << 32))
    }

    fn flag(&mut self) -> Result<bool, HarnessError> {
        match self.word()? {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(harness("glue metadata: invalid flag")),
        }
    }

    fn block(&mut self) -> Result<&'a [f64], HarnessError> {
        let count = self.count()?;
        let end = self
            .at
            .checked_add(count)
            .ok_or_else(|| harness("glue metadata: block length overflow"))?;
        let result = self
            .raw
            .get(self.at..end)
            .ok_or_else(|| harness("glue metadata: truncated block"))?;
        self.at = end;
        Ok(result)
    }
}

fn glue_topology_rows(brep: &Brep) -> Vec<(u32, u32, ObjectRef)> {
    let mut rows = vec![(10, 0, glue_query_key(CacheObject::Brep))];
    let mut push = |kind, index, entity| rows.push((kind, index, glue_key(entity)));
    for id in brep.regions().ids() {
        push(0, id.index(), EntityRef::Region(id));
    }
    for id in brep.shells().ids() {
        push(1, id.index(), EntityRef::Shell(id));
    }
    for id in brep.faceuses().ids() {
        push(2, id.index(), EntityRef::Faceuse(id));
    }
    for id in brep.faces().ids() {
        push(3, id.index(), EntityRef::Face(id));
    }
    for id in brep.loopuses().ids() {
        push(4, id.index(), EntityRef::Loopuse(id));
    }
    for id in brep.loops().ids() {
        push(5, id.index(), EntityRef::Loop(id));
    }
    for id in brep.edgeuses().ids() {
        push(6, id.index(), EntityRef::Edgeuse(id));
    }
    for id in brep.edges().ids() {
        push(7, id.index(), EntityRef::Edge(id));
    }
    for id in brep.vertexuses().ids() {
        push(8, id.index(), EntityRef::Vertexuse(id));
    }
    for id in brep.vertices().ids() {
        push(9, id.index(), EntityRef::Vertex(id));
    }
    rows
}

/// Actual notification inputs imported from the native fixture's graph.
struct GlueNotifications {
    records: BTreeMap<ObjectKey, Topology>,
    attributes: AttributeStore,
}

/// Decode all fields before the operation. The native getters serialize flags,
/// all stored marks and user data as exact u32 limbs. Attribute/callback presence
/// is checked only when `GlueEdgesGeneral` will reach notifications; its null and
/// same-edge returns do not consume those unsupported fields.
#[allow(clippy::too_many_lines)]
fn import_glue_metadata(
    brep: &mut Brep,
    raw: &[f64],
    reaches_notifications: bool,
) -> Result<(QueryMetadata, GlueNotifications), Fail> {
    let mut reader = GlueMetadataReader { raw, at: 0 };
    if reader.word()? != 1 {
        return Err(hfail("glue metadata: unsupported version"));
    }
    let query = reader.block()?;
    let callback_present = reader.flag()?;
    let expected = glue_topology_rows(brep);
    if reader.count()? != expected.len() {
        return Err(hfail("glue metadata: topology count differs from graph"));
    }
    let mut records = BTreeMap::new();
    let mut attributed_objects = 0_usize;
    for (kind, index, object) in expected {
        if reader.word()? != kind || reader.word()? != index {
            return Err(hfail("glue metadata: topology identity differs from graph"));
        }
        let mut record = Topology::new();
        record.set_all_flags(reader.wide()?);
        record.restore_marks([
            reader.wide()?,
            reader.wide()?,
            reader.wide()?,
            reader.wide()?,
            reader.wide()?,
        ]);
        record.set_user_index1(reader.wide()?);
        record.set_user_index2(reader.wide()?);
        record.set_user_ptr1(reader.wide()?);
        if reader.flag()? {
            attributed_objects += 1;
        }
        if records.insert(object.key, record).is_some() {
            return Err(hfail("glue metadata: colliding topology keys"));
        }
    }
    let geometry = brep
        .edges()
        .iter()
        .map(|(id, edge)| (0, id.index(), edge.curve.is_some()))
        .chain(
            brep.faces()
                .iter()
                .map(|(id, face)| (1, id.index(), face.surface.is_some())),
        )
        .chain(
            brep.edgeuses()
                .iter()
                .map(|(id, edgeuse)| (2, id.index(), edgeuse.uv_trim_curve.is_some())),
        )
        .collect::<Vec<_>>();
    if reader.count()? != geometry.len() {
        return Err(hfail("glue metadata: geometry count differs from graph"));
    }
    for (kind, index, present) in geometry {
        if reader.word()? != kind || reader.word()? != index || reader.flag()? != present {
            return Err(hfail("glue metadata: geometry identity differs from graph"));
        }
        let attributed = reader.flag()?;
        if attributed && !present {
            return Err(hfail("glue metadata: attributes on absent geometry"));
        }
        if attributed {
            attributed_objects += 1;
        }
    }
    if reader.at != raw.len() {
        return Err(hfail("glue metadata: trailing values"));
    }
    let metadata = apply_query_metadata(brep, query, &glue_query_key)?;
    // Keep the existing lazy geometry checks: a later curve query diagnoses its
    // actual dependency. Notification import does not consume all geometry.
    if let Some(output) = metadata.deferred_output_for(false, false) {
        return Err(Fail::Deferred(Box::new(output)));
    }
    if reaches_notifications && (callback_present || attributed_objects != 0) {
        return Err(Fail::Deferred(Box::new(Output::deferred(
            "harness:topology-notification-export",
            format!(
                "glue native input has callbacks={callback_present}, attributed_objects={attributed_objects}; those handlers are not decoded"
            ),
        ))));
    }
    Ok((
        metadata,
        GlueNotifications {
            records,
            // Every native attribute slot was checked. An empty store is used
            // only when exact, or on a no-op that never observes notification state.
            attributes: AttributeStore::default(),
        },
    ))
}

fn glue_metadata(
    input: &BrepIn,
    brep: &mut Brep,
    reaches_notifications: bool,
) -> Result<(QueryMetadata, GlueNotifications), Fail> {
    let raw = ffi::eval(
        Op::GlueMetadata,
        &shared::path_of(input)?,
        &brep_block(input)?,
    )
    .map_err(|error| hfail(format!("c++ glue metadata: {error}")))?;
    import_glue_metadata(brep, &raw, reaches_notifications)
}

fn query_metadata(input: &BrepIn, brep: &mut Brep) -> Result<QueryMetadata, Fail> {
    let raw = ffi::eval(
        Op::QueryMetadata,
        &shared::path_of(input)?,
        &brep_block(input)?,
    )
    .map_err(|error| hfail(format!("c++ query metadata: {error}")))?;
    let metadata = apply_query_metadata(brep, &raw, &query_key)?;
    if let Some(output) = metadata.deferred_output_for(false, false) {
        return Err(Fail::Deferred(Box::new(output)));
    }
    Ok(metadata)
}

fn query_status<T>(result: &Result<T, SolidError>, queries: &RayQueryContext) -> Result<f64, Fail> {
    if matches!(result, Err(SolidError::Unported { unit, .. }) if *unit == "harness:analytic-geometry-export")
    {
        if let Some(output) = queries.metadata().deferred_output() {
            return Err(Fail::Deferred(Box::new(output)));
        }
    }
    status(result)
}

/// The record of a run.
type Record = Result<Vec<f64>, Fail>;

/// The region, shell, face, loop, edge and vertex counts of a B-rep as the shim's walk
/// counts them: regions in the B-rep, shells in a region, faces with a face use, loops
/// with a loop use, edges with edge uses and vertices with vertex uses (arena orphans,
/// e.g. a glued-away edge, are not counted).
fn brep_counts(b: &Brep) -> [f64; 6] {
    [
        b.regions().iter().filter(|(_, r)| r.in_brep).count() as f64,
        b.shells()
            .iter()
            .filter(|(_, s)| s.region.is_some())
            .count() as f64,
        b.faces()
            .iter()
            .filter(|(_, f)| f.faceuse.is_some())
            .count() as f64,
        b.loops()
            .iter()
            .filter(|(_, l)| l.loopuse.is_some())
            .count() as f64,
        b.edges()
            .iter()
            .filter(|(_, e)| !e.edgeuses.is_empty())
            .count() as f64,
        b.vertices()
            .iter()
            .filter(|(_, v)| !v.vertexuses.is_empty())
            .count() as f64,
    ]
}

/// The tolerances of the edges and vertices [`brep_counts`] counts, each sorted ascending.
fn brep_tolerances(b: &Brep) -> Vec<f64> {
    let mut te: Vec<f64> = b
        .edges()
        .iter()
        .filter(|(_, e)| !e.edgeuses.is_empty())
        .map(|(_, e)| e.tolerance)
        .collect();
    let mut tv: Vec<f64> = b
        .vertices()
        .iter()
        .filter(|(_, v)| !v.vertexuses.is_empty())
        .map(|(_, v)| v.tolerance)
        .collect();
    te.sort_by(f64::total_cmp);
    tv.sort_by(f64::total_cmp);
    te.extend(tv);
    te
}

/// Whether `MakeFaceWithCurves` fails its input checks (`SmBrep.cpp:15722-15741`) on a
/// face with 3D curves: no loops, or a surface degenerate to a curve or a point (the
/// port's tests, which fail soft as upstream's flags).
fn mfwc_input_error(surface: &SrfIn, loops: &[usize]) -> Result<bool, Fail> {
    if loops.is_empty() {
        return Ok(true);
    }
    let srf = rust_surface(surface)?;
    let flag = |r: solidcore::Result<bool>| -> Result<bool, Fail> {
        match r {
            Ok(x) => Ok(x),
            Err(e) => match Output::unported(&e) {
                Some(o) => Err(Fail::Deferred(Box::new(o))),
                None => Ok(false),
            },
        }
    };
    let tol = solidcore::surface::surface::p1::DEGENERATE_TOL_3D;
    Ok(
        flag(srf.is_degenerate_curve(solidcore::math::EFF_ZERO_DEG, None, tol))?
            || flag(srf.is_degenerate_point(tol, None))?,
    )
}

/// The record of a face maker: `status [face shell region counts]`.
fn new_face_record(b: &Brep, r: &Result<p7::NewFace, SolidError>) -> Record {
    let mut o = vec![status(r)?];
    if let Ok(nf) = r {
        o.extend([
            flag(nf.face.is_some()),
            flag(nf.shell.is_some()),
            flag(nf.region.is_some()),
        ]);
        o.extend(brep_counts(b));
    }
    Ok(o)
}

#[allow(clippy::too_many_lines)] // one record per op
fn rust_record(k: &mut Kernel, i: &Inputs, ex: Option<&[f64]>) -> Record {
    match i {
        Inputs::Read {
            file,
            damage,
            rebuild,
        } => {
            let bytes = read_bytes(file, *damage)?;
            // the C++ side does not run a modified stream its SmBrepData read accepts
            // (MakeTopologyFromData does not validate the topology: UB-164); decided with
            // the same C++ reader as the FFI
            if !is_test_file(file, &bytes)? && cpp_brep_data_reads(&bytes)? {
                return Err(Fail::Deferred(Box::new(Output::deferred(
                    "cpp-guard:unvalidated-topology",
                    "SmBrep::MakeTopologyFromData",
                ))));
            }
            // a B-rep as `new SmBrep(ctx)` makes it (with the infinite region)
            let mut breps = Arena::<BrepId, Brep>::default();
            let id = breps
                .insert(Brep::new_with_infinite_region().map_err(hfail)?)
                .map_err(hfail)?;
            let b = breps
                .get_mut(id)
                .ok_or_else(|| hfail("missing registered B-rep"))?;
            // ReadFromFile builds real geometry and owners from BrepData, without
            // the legacy fixture importer's analytic stand-ins.
            let mut queries = RayQueryContext::new(b, id, QueryMetadata::default(), read_query_key);
            let mut records = BTreeMap::new();
            let mut attrs = AttributeStore::default();
            let mut env = p10::TopologyEnv {
                all_attributes: &[],
                attributes: &mut attrs,
                key_of: &read_key_of,
                records: &mut records,
                next_key: 1u64 << 50,
                region_queries: &mut queries,
            };
            let r = p7::read_from_file(k, b, bytes, FileType::Ascii, *rebuild, &mut env);
            let mut o = vec![status(&r)?];
            if r.is_ok() {
                o.extend(brep_counts(b));
            }
            Ok(o)
        }
        Inputs::MakeFaceWithCurves {
            surface,
            dom,
            orient,
            region_null,
            crv_on_surf,
            loops,
            has3d,
            has_uv,
            curves,
            points,
        } => {
            // the shim's guard: upstream dereferences the null UV curve array once MakeFace
            // succeeds (the port returns InvalidInput there); its own input checks (no
            // loops, a degenerate surface) return first
            if *crv_on_surf == 2 && !*has_uv && *has3d && !mfwc_input_error(surface, loops)? {
                return Err(cpp_null_deref(
                    "cpp-guard:ub-165",
                    "SmBrep::MakeFaceWithCurves",
                ));
            }
            let mut b = Brep::new();
            let region = if *region_null {
                None
            } else {
                Some(
                    b.add_region_record(Region {
                        in_brep: true,
                        ..Region::default()
                    })
                    .map_err(hfail)?,
                )
            };
            let mut c3 = Vec::new();
            let mut uv = Vec::new();
            for c in curves {
                if *has3d {
                    let x = c.c3.as_ref().ok_or_else(|| hfail("missing 3d curve"))?;
                    c3.push(Some(rust_curve(x)?));
                }
                if *has_uv {
                    let x = c.uv.as_ref().ok_or_else(|| hfail("missing uv curve"))?;
                    uv.push(rust_curve(x)?);
                }
            }
            let input = p7::FaceCurvesInput {
                loops: loops.clone(),
                curves_3d: has3d.then_some(c3),
                uv_curves: has_uv.then_some(uv),
                orientations: curves
                    .iter()
                    .map(|c| orient_of(c.orient))
                    .collect::<Result<_, _>>()?,
                loop_points: points.iter().map(|p| Vec3::new(p[0], p[1], p[2])).collect(),
            };
            let srf = rust_surface(surface)?;
            let domain =
                Extent2d::from_points(Point2::new(dom[0], dom[1]), Point2::new(dom[2], dom[3]));
            let cos = match crv_on_surf {
                0 => None,
                1 => Some(false),
                2 => Some(true),
                _ => return Err(hfail("crv_on_surf code")),
            };
            let mut records = BTreeMap::new();
            let mut attrs = AttributeStore::default();
            let mut env = ed::NotifyEnv {
                records: &mut records,
                attributes: &mut attrs,
                context: None,
                key_of: &key_of,
            };
            let r = p7::make_face_with_curves(
                k,
                &mut b,
                region,
                input,
                srf,
                &domain,
                orient_of(*orient)?,
                cos,
                &this_face,
                &mut env,
            );
            new_face_record(&b, &r)
        }
        Inputs::MakeFace {
            kind,
            side,
            orient,
            region_null,
        } => {
            if *kind == 1 {
                // the shim's square wire loop: its -3 when a side is too small for a line
                let p = [
                    Vec3::new(0.0, 0.0, 0.0),
                    Vec3::new(*side, 0.0, 0.0),
                    Vec3::new(*side, *side, 0.0),
                    Vec3::new(0.0, *side, 0.0),
                ];
                for (a, b) in p.iter().zip(p.iter().cycle().skip(1)) {
                    if Line::create_line_segment(3, *a, *b, None).is_err() {
                        return Ok(vec![NOT_RUN]);
                    }
                }
                // The shim's loop starts with MakeShellVertex in a fresh B-rep's region;
                // run it, then MakeWireEdgeVertex and MakeWireEdge (topology.brep_p3)
                // come before MakeFace sees the loop.
                let mut b = Brep::new();
                let region = b
                    .add_region_record(Region {
                        in_brep: true,
                        ..Region::default()
                    })
                    .map_err(hfail)?;
                let mut records = BTreeMap::new();
                let mut attrs = AttributeStore::default();
                let mut env = ed::NotifyEnv {
                    records: &mut records,
                    attributes: &mut attrs,
                    context: None,
                    key_of: &key_of,
                };
                if let Err(e) = solidcore::topology::brep::p9::make_shell_vertex(
                    &mut b,
                    Some(region),
                    p[0],
                    &mut env,
                ) {
                    return Err(match Output::unported(&e) {
                        Some(o) => Fail::Deferred(Box::new(o)),
                        None => hfail(e),
                    });
                }
                return Err(Fail::Deferred(Box::new(Output::deferred(
                    "topology.brep_p3",
                    "SmBrep::MakeWireEdgeVertex",
                ))));
            }
            let mut b = Brep::new();
            let region = if *region_null {
                None
            } else {
                Some(
                    b.add_region_record(Region {
                        in_brep: true,
                        ..Region::default()
                    })
                    .map_err(hfail)?,
                )
            };
            let surface = (*kind == 1).then(|| {
                Surface::from(Plane::from_normal(
                    Vec3::new(0.0, 0.0, 0.0),
                    Vec3::new(0.0, 0.0, 1.0),
                ))
            });
            let dom = Extent2d::from_points(Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));
            let r = p7::make_face(
                k,
                &mut b,
                region,
                &[],
                &[],
                &[],
                surface,
                &dom,
                orient_of(*orient)?,
                &p7::MakeFaceOptions::default(),
                &this_face,
            );
            new_face_record(&b, &r)
        }
        Inputs::Glue {
            brep,
            e1,
            e2,
            orient,
            gap,
            nesting,
        } => {
            let ex = ex.ok_or_else(|| hfail("no export"))?;
            let (mut b, _) = shared::import(ex)?;
            let reaches_notifications = *e1 >= 0 && *e2 >= 0 && e1 != e2;
            let (metadata, mut notifications) = glue_metadata(brep, &mut b, reaches_notifications)?;
            let mut breps = Arena::<BrepId, Brep>::default();
            let id = breps.insert(b).map_err(hfail)?;
            let b = breps
                .get_mut(id)
                .ok_or_else(|| hfail("missing registered B-rep"))?;
            let mut queries = RayQueryContext::new(b, id, metadata, glue_query_key);
            let eid = |x: i64| (x >= 0).then(|| EdgeId::from_index(x as u32));
            // the shim's guard: upstream dereferences a null edge curve
            // (CurveMaxDistanceBetween) when gluing two distinct non-wire edges with
            // SM_OT_SAME, the same intervals and exactly one curve (two null curves
            // compare equal)
            if let (Some(a), Some(d)) = (eid(*e1), eid(*e2)) {
                if a != d && !ed::is_wire(b, d).map_err(hfail)? && *orient == 1 {
                    let (ia, id) = (
                        ed::interval(b, a).map_err(hfail)?,
                        ed::interval(b, d).map_err(hfail)?,
                    );
                    let one_null = b.edge(a).map_err(hfail)?.curve.is_none()
                        != b.edge(d).map_err(hfail)?.curve.is_none();
                    if are_same(ia.min(), id.min()) && are_same(ia.max(), id.max()) && one_null {
                        return Err(cpp_null_deref(
                            "cpp-guard:null-curve",
                            "SmBrep::GlueEdgesGeneral",
                        ));
                    }
                }
            }
            let mut env = ed::NotifyEnv {
                records: &mut notifications.records,
                attributes: &mut notifications.attributes,
                context: None, // native callback absence was checked before notification
                key_of: &glue_key,
            };
            let r = p7::glue_edges_general(
                k,
                b,
                eid(*e1),
                orient_of(*orient)?,
                *gap,
                *nesting,
                eid(*e2),
                None,
                &mut env,
                &mut queries,
            );
            let s = query_status(&r, &queries)?;
            // the B-rep's counts and tolerances after the op (the shim walks the C++ B-rep
            // afterwards)
            let mut o = vec![0.0, s];
            o.extend(brep_counts(b));
            o.extend(brep_tolerances(b));
            Ok(o)
        }
        Inputs::Copy { n, merge, .. } => {
            let ex = ex.ok_or_else(|| hfail("no export"))?;
            let (mut b, _) = shared::import(ex)?;
            let faces: Vec<FaceId> = b.faces().ids().take(*n as usize).collect();
            let mut dest = Brep::new();
            let mut records = BTreeMap::new();
            let mut attrs = AttributeStore::default();
            let mut env = ed::NotifyEnv {
                records: &mut records,
                attributes: &mut attrs,
                context: None,
                key_of: &key_of,
            };
            let r = p7::copy_faces(k, &mut b, &faces, &mut dest, *merge, &this_face, &mut env);
            // (a face that could not be copied makes upstream's status SM_ERR; the copies
            // made are kept)
            let s = match &r {
                Ok((_, false)) => f64::from(solidcore::error::SM_ERR),
                _ => status(&r)?,
            };
            let mut o = vec![0.0, s, r.as_ref().map_or(0, |(f, _)| f.len()) as f64];
            o.extend(brep_counts(&dest));
            Ok(o)
        }
        Inputs::MergeSelf { brep, face, tol } => {
            let ex = ex.ok_or_else(|| hfail("no export"))?;
            let (mut b, _) = shared::import(ex)?;
            let metadata = query_metadata(brep, &mut b)?;
            let mut breps = Arena::<BrepId, Brep>::default();
            let id = breps.insert(b).map_err(hfail)?;
            let b = breps
                .get_mut(id)
                .ok_or_else(|| hfail("missing registered B-rep"))?;
            let mut queries = RayQueryContext::new(b, id, metadata, query_key);
            let f = FaceId::from_index(*face as u32);
            if b.face(f).map_err(hfail)?.surface.is_none() {
                return Err(cpp_null_deref(
                    "cpp-guard:null-surface",
                    "SmBrep::MergeSelfIntersectionCurves",
                ));
            }
            let mut records = BTreeMap::new();
            let mut attrs = AttributeStore::default();
            let mut env = ed::NotifyEnv {
                records: &mut records,
                attributes: &mut attrs,
                context: None,
                key_of: &key_of,
            };
            let r = p7::merge_self_intersection_curves(
                k,
                b,
                Some(f),
                *tol,
                &[],
                &mut env,
                &mut queries,
            );
            let s = query_status(&r, &queries)?;
            let n = r.map_or(0, |v| v.len());
            Ok(vec![0.0, s, n as f64])
        }
    }
}

fn run_rust(i: &Inputs) -> Result<Output, HarnessError> {
    if let Inputs::Read { file, damage, .. } = i {
        let bytes = read_bytes(file, *damage)?;
        if !guard_accepts(&bytes, is_test_file(file, &bytes)?)? {
            return Ok(Output::err("guard"));
        }
    }
    let (r, ex) = match resolve(i)? {
        Ok(x) => x,
        Err(o) => return Ok(o),
    };
    let mut k = Kernel::default();
    match rust_record(&mut k, &r, ex.as_deref().map(Vec::as_slice)) {
        Ok(v) => Ok(Output::new().values("out", v)),
        Err(Fail::Deferred(o)) => Ok(*o),
        Err(Fail::Harness(h)) => Err(h),
    }
}

// ---- cases ------------------------------------------------------------------------------

/// A degree-1 B-spline line of dimension `dim` from `a` to `b` (`z = NOZ` in 2D).
fn line(dim: u64, a: [f64; 3], b: [f64; 3]) -> CrvIn {
    let (za, zb) = if dim == 2 { (NOZ, NOZ) } else { (a[2], b[2]) };
    CrvIn {
        dim,
        cpts: vec![[a[0], a[1], za, NOW], [b[0], b[1], zb, NOW]],
        p: 1,
        u: vec![0.0, 0.0, 1.0, 1.0],
    }
}

/// A square loop of side `s` at `(x, y)` on `z = 0`, as 3D and UV lines, orientation
/// `orient` for every curve, run counter-clockwise or clockwise.
fn square(x: f64, y: f64, s: f64, ccw: bool) -> Vec<FaceCurveIn> {
    let mut pts = [[x, y], [x + s, y], [x + s, y + s], [x, y + s]];
    if !ccw {
        pts.reverse();
    }
    let at = |i: usize| pts.get(i % 4).copied().unwrap_or([0.0, 0.0]);
    (0..4)
        .map(|i| {
            let a = at(i);
            let b = at(i + 1);
            FaceCurveIn {
                c3: Some(line(3, [a[0], a[1], 0.0], [b[0], b[1], 0.0])),
                uv: Some(line(2, [a[0], a[1], 0.0], [b[0], b[1], 0.0])),
                orient: 1,
            }
        })
        .collect()
}

/// A bilinear B-spline patch over `[0, s]^2` on `z = 0`.
fn patch(s: f64) -> SrfIn {
    SrfIn::Bspline {
        net: vec![
            vec![[0.0, 0.0, 0.0, NOW], [0.0, s, 0.0, NOW]],
            vec![[s, 0.0, 0.0, NOW], [s, s, 0.0, NOW]],
        ],
        p: 1,
        q: 1,
        u: vec![0.0, 0.0, s, s],
        v: vec![0.0, 0.0, s, s],
    }
}

/// The face-from-curves cases: real faces (outer squares, with a hole, 3D / UV / both,
/// curve on surface), and the input errors.
fn face_cases() -> Vec<Inputs> {
    let plane = SrfIn::Plane {
        o: [0.0, 0.0, 0.0],
        n: [0.0, 0.0, 1.0],
    };
    let mut v = Vec::new();
    for (si, srf) in [plane.clone(), patch(10.0)].into_iter().enumerate() {
        for (has3d, has_uv) in [(true, false), (false, true), (true, true)] {
            for region_null in [false, true] {
                for cos in [0u8, 1, 2] {
                    if cos != 0 && has3d {
                        continue;
                    }
                    let mut curves = square(1.0, 1.0, 8.0, true);
                    let mut loops = vec![4];
                    if si == 0 {
                        curves.extend(square(3.0, 3.0, 2.0, false));
                        loops.push(4);
                    }
                    v.push(Inputs::MakeFaceWithCurves {
                        surface: srf.clone(),
                        dom: [0.0, 0.0, 10.0, 10.0],
                        orient: 1,
                        region_null,
                        crv_on_surf: cos,
                        loops,
                        has3d,
                        has_uv,
                        curves,
                        points: Vec::new(),
                    });
                }
            }
        }
    }
    // input errors: no curves, no loops, degenerate surfaces, all curves degenerate
    let base =
        |loops: Vec<usize>, has3d: bool, has_uv: bool, curves: Vec<FaceCurveIn>, srf: SrfIn| {
            Inputs::MakeFaceWithCurves {
                surface: srf,
                dom: [0.0, 0.0, 10.0, 10.0],
                orient: 1,
                region_null: false,
                crv_on_surf: 0,
                loops,
                has3d,
                has_uv,
                curves,
                points: Vec::new(),
            }
        };
    v.push(base(vec![], false, false, vec![], plane.clone()));
    v.push(base(vec![], true, false, vec![], plane.clone()));
    v.push(base(vec![4], false, false, Vec::new(), plane.clone()));
    let flat = SrfIn::Bspline {
        net: vec![
            vec![[0.0, 0.0, 0.0, NOW], [0.0, 0.0, 0.0, NOW]],
            vec![[5.0, 0.0, 0.0, NOW], [5.0, 0.0, 0.0, NOW]],
        ],
        p: 1,
        q: 1,
        u: vec![0.0, 0.0, 1.0, 1.0],
        v: vec![0.0, 0.0, 1.0, 1.0],
    };
    let point = SrfIn::Bspline {
        net: vec![
            vec![[1.0, 1.0, 1.0, NOW], [1.0, 1.0, 1.0, NOW]],
            vec![[1.0, 1.0, 1.0, NOW], [1.0, 1.0, 1.0, NOW]],
        ],
        p: 1,
        q: 1,
        u: vec![0.0, 0.0, 1.0, 1.0],
        v: vec![0.0, 0.0, 1.0, 1.0],
    };
    v.push(base(
        vec![4],
        true,
        false,
        square(1.0, 1.0, 8.0, true),
        flat,
    ));
    v.push(base(
        vec![4],
        true,
        false,
        square(1.0, 1.0, 8.0, true),
        point,
    ));
    v
}

/// A B-rep with faces (a wire or the network would make a case both sides reject).
fn gen_face_brep(rng: &mut Rng) -> BrepIn {
    loop {
        let b = shared::gen_brep(rng);
        if matches!(b, BrepIn::File { .. } | BrepIn::Face { .. }) {
            return b;
        }
    }
}

/// A generated op: a damaged read, a glue, a merge-self, a face from curves, a copy or a
/// `MakeFace`, each under the same guards as the built-in cases.
fn gen_inputs(rng: &mut Rng) -> Inputs {
    let file = (*rng.choose(&FILES).unwrap_or(&"disk.smb")).to_string();
    match rng.below(6) {
        0 => {
            // damaged streams the B-rep data guard accepts (a rejected one would be a
            // case both sides refuse; give up after a few tries)
            let mut damage = Damage::Truncate(0);
            for _ in 0..16 {
                let d = if rng.chance(0.5) {
                    Damage::Truncate(rng.below(4000) as usize)
                } else {
                    Damage::Overwrite {
                        at: rng.below(1 << 16) as usize,
                        byte: *rng.choose(b"0123456789 -.eX\n").unwrap_or(&b'9'),
                    }
                };
                damage = d;
                let ok = read_bytes(&file, Some(d))
                    .and_then(|bytes| guard_accepts(&bytes, false))
                    .unwrap_or(false);
                if ok {
                    break;
                }
            }
            Inputs::Read {
                file,
                damage: Some(damage),
                rebuild: rng.chance(0.5),
            }
        }
        1 => Inputs::Glue {
            brep: shared::gen_brep(rng),
            e1: rng.range_i64(-1, 1 << 10),
            e2: rng.range_i64(-1, 1 << 10),
            orient: rng.below(3) as u8,
            gap: rng.uniform(0.0, 1.0e-3),
            nesting: rng.chance(0.5),
        },
        2 => Inputs::MergeSelf {
            brep: gen_face_brep(rng),
            face: rng.below(1 << 10) as i64,
            tol: rng.uniform(1.0e-6, 1.0e-2),
        },
        3 => {
            // a square outer loop on a plane or a patch over [0, 10]^2, closed or open
            // (2 or 3 sides: a gap target), sometimes with a square hole; 3D, UV or both
            let s = rng.uniform(0.5, 8.0);
            let (x, y) = (rng.uniform(0.0, 10.0 - s), rng.uniform(0.0, 10.0 - s));
            let sides = 2 + rng.below(3) as usize;
            let mut curves = square(x, y, s, true);
            curves.truncate(sides);
            let mut loops = vec![sides];
            if sides == 4 && rng.chance(0.3) {
                curves.extend(square(x + s / 4.0, y + s / 4.0, s / 2.0, false));
                loops.push(4);
            }
            // sometimes every SmOrientType code on the curves
            if rng.chance(0.3) {
                for c in &mut curves {
                    c.orient = rng.below(3) as u8;
                }
            }
            let (has3d, has_uv) = *rng
                .choose(&[(true, false), (false, true), (true, true)])
                .unwrap_or(&(true, true));
            let surface = if rng.chance(0.5) {
                SrfIn::Plane {
                    o: [0.0, 0.0, 0.0],
                    n: [0.0, 0.0, 1.0],
                }
            } else {
                patch(10.0)
            };
            Inputs::MakeFaceWithCurves {
                surface,
                dom: [0.0, 0.0, 10.0, 10.0],
                orient: rng.below(3) as u8,
                region_null: rng.chance(0.5),
                crv_on_surf: rng.below(3) as u8,
                loops,
                has3d,
                has_uv,
                curves,
                points: Vec::new(),
            }
        }
        4 => Inputs::Copy {
            brep: shared::gen_brep(rng),
            n: rng.below(8),
            merge: rng.chance(0.5),
        },
        _ => Inputs::MakeFace {
            kind: rng.below(2) as u8,
            side: rng.uniform(0.1, 10.0),
            orient: rng.below(3) as u8,
            region_null: rng.chance(0.5),
        },
    }
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "topology.brep_p7"
    }

    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        run_cpp(&case.parse()?)
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        run_rust(&case.parse()?)
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        serde_json::to_value(gen_inputs(rng)).ok()
    }

    fn cases(&self, ctx: &CaseContext) -> Result<Vec<Case>, HarnessError> {
        let mut v = default_cases(self, ctx)?;
        let push = |v: &mut Vec<Case>, name: String, i: &Inputs| -> Result<(), HarnessError> {
            let index = v.len();
            v.push(Case {
                unit: "topology.brep_p7".to_string(),
                name,
                source: crate::CaseSource::Fixture {
                    path: "built-in: topology.brep_p7 cases".to_string(),
                    index,
                    upstream_ref: Some("tests/prog_test (pt_TestFiles B-reps)".to_string()),
                },
                inputs: serde_json::to_value(i).map_err(|e| harness(e.to_string()))?,
            });
            Ok(())
        };
        for file in FILES {
            for rebuild in [false, true] {
                push(
                    &mut v,
                    format!("read_{file}_{}", u8::from(rebuild)),
                    &Inputs::Read {
                        file: file.to_string(),
                        damage: None,
                        rebuild,
                    },
                )?;
            }
            // glue no-ops and errors, merge without curves, copies of small files
            let brep = BrepIn::File {
                file: file.to_string(),
                rebuild: false,
            };
            // ((0, 1), (0, 4) and (2, 5) include every distinct edge pair of the files that
            // shares both vertices: Mouse (0, 4) (2, 5); disk, surf_sphere150, torus, toruswvy (0, 1))
            for (j, (e1, e2)) in [
                (-1, 0),
                (0, -1),
                (-1, -1),
                (0, 0),
                (3, 3),
                (0, 1),
                (1, 5),
                (2, 7),
                (0, 4),
                (2, 5),
            ]
            .into_iter()
            .enumerate()
            {
                push(
                    &mut v,
                    format!("glue_{file}_{j}"),
                    &Inputs::Glue {
                        brep: brep.clone(),
                        e1,
                        e2,
                        orient: 1 + (j % 2) as u8,
                        gap: 1.0e-6,
                        nesting: j % 3 == 0,
                    },
                )?;
            }
            for f in 0..4 {
                push(
                    &mut v,
                    format!("merge_self_{file}_f{f}"),
                    &Inputs::MergeSelf {
                        brep: brep.clone(),
                        face: f,
                        tol: 1.0e-5,
                    },
                )?;
            }
        }
        for file in ["plane_square.smb", "disk.smb", "plane_circ.smb"] {
            for merge in [false, true] {
                push(
                    &mut v,
                    format!("copy_{file}_{}", u8::from(merge)),
                    &Inputs::Copy {
                        brep: BrepIn::File {
                            file: file.to_string(),
                            rebuild: false,
                        },
                        n: 2,
                        merge,
                    },
                )?;
            }
        }
        for (j, i) in face_cases().iter().enumerate() {
            push(&mut v, format!("face_curves_{j}"), i)?;
        }
        // pre-review reproducers: crv_on_surf TRUE with 3D but no UV curves (upstream
        // dereferences the null UV array, SmBrep.cpp:16600), a loop of no curves, and
        // more curves than the loops use
        let sq = square(1.0, 1.0, 8.0, true);
        let plane = SrfIn::Plane {
            o: [0.0, 0.0, 0.0],
            n: [0.0, 0.0, 1.0],
        };
        for (j, (loops, curves)) in [
            (vec![4usize], sq.clone()),
            (vec![0], Vec::new()),
            (vec![3], sq.clone()),
        ]
        .into_iter()
        .enumerate()
        {
            for cos in [1u8, 2] {
                push(
                    &mut v,
                    format!("face_curves_cos_no_uv_{j}_{cos}"),
                    &Inputs::MakeFaceWithCurves {
                        surface: plane.clone(),
                        dom: [0.0, 0.0, 10.0, 10.0],
                        orient: 1,
                        region_null: false,
                        crv_on_surf: cos,
                        loops: loops.clone(),
                        has3d: true,
                        has_uv: false,
                        curves: curves.clone(),
                        points: Vec::new(),
                    },
                )?;
            }
        }
        // an open UV-only loop (three sides of the square) with crv_on_surf TRUE: the gap
        // target takes upstream's pOptSurfCrvOnSurf branch of DeformToEndPointTargets
        push(
            &mut v,
            "face_curves_cos_open_loop".to_string(),
            &Inputs::MakeFaceWithCurves {
                // (a B-spline patch: a curve on surface over a plane defers earlier, to
                // curve.crv_on_surf)
                surface: patch(10.0),
                dom: [0.0, 0.0, 10.0, 10.0],
                orient: 1,
                region_null: false,
                crv_on_surf: 2,
                loops: vec![3],
                has3d: false,
                has_uv: true,
                curves: sq.iter().take(3).cloned().collect(),
                points: Vec::new(),
            },
        )?;
        // a square too small for its lines (CreateLineSegment fails: -3 on both sides) and
        // the smallest side that still makes a face in C++
        for (name, side) in [
            ("make_face_tiny_1e-300", 1.0e-300),
            ("make_face_small_1e-12", 1.0e-12),
        ] {
            push(
                &mut v,
                name.to_string(),
                &Inputs::MakeFace {
                    kind: 1,
                    side,
                    orient: 1,
                    region_null: false,
                },
            )?;
        }
        // damaged streams the B-rep data guard accepts and whose read fails: the compared
        // read errors (a digit of an integer token overwritten; deterministic search)
        let k = Kernel::default();
        for file in ["disk.smb", "plane_square.smb", "cone.smb", "Mouse.smb"] {
            let bytes = read_bytes(file, None)?;
            let mut found = 0;
            for at in (0..bytes.len().min(20_000)).step_by(3) {
                if found >= 3 {
                    break;
                }
                if !bytes.get(at).is_some_and(u8::is_ascii_digit) {
                    continue;
                }
                for byte in *b"90" {
                    let d = Damage::Overwrite { at, byte };
                    let b2 = read_bytes(file, Some(d))?;
                    if b2 == bytes || !guard_accepts(&b2, false)? {
                        continue;
                    }
                    let fails = solidcore::topology::brep_data::BrepData::read_from_file(
                        &k,
                        b2,
                        FileType::Ascii,
                    )
                    .is_err();
                    if fails {
                        push(
                            &mut v,
                            format!("read_{file}_error{found}"),
                            &Inputs::Read {
                                file: file.to_string(),
                                damage: Some(d),
                                rebuild: false,
                            },
                        )?;
                        found += 1;
                        break;
                    }
                }
            }
        }
        for kind in [0u8, 1] {
            for region_null in [false, true] {
                for orient in [1u8, 2] {
                    push(
                        &mut v,
                        format!("make_face_{kind}_{}_{orient}", u8::from(region_null)),
                        &Inputs::MakeFace {
                            kind,
                            side: 4.0,
                            orient,
                            region_null,
                        },
                    )?;
                }
            }
        }
        Ok(v)
    }
}

#[cfg(test)]
mod glue_notification_tests {
    use super::*;

    #[test]
    fn glue_with_native_records_reaches_native_vertex_rejection() {
        let input = Inputs::Glue {
            brep: BrepIn::File {
                file: "box_with_large_fillet.smb".to_owned(),
                rebuild: false,
            },
            e1: 0,
            e2: 1,
            orient: 2,
            gap: 1.0e-6,
            nesting: false,
        };
        let (resolved, export) = resolve(&input)
            .unwrap()
            .unwrap_or_else(|_| panic!("native regression fixture must have both edges"));
        let (op, path, encoded) = encode(&resolved).unwrap();
        let native = ffi::eval(op, &path, &encoded).unwrap();
        assert_eq!(native.get(..2), Some([0.0, 1001.0].as_slice()));
        let result = rust_record(
            &mut Kernel::default(),
            &resolved,
            export.as_deref().map(Vec::as_slice),
        );
        let actual = match result {
            Ok(values) => values,
            Err(Fail::Deferred(_)) => panic!("vertex rejection must be compared, not deferred"),
            Err(Fail::Harness(error)) => panic!("native glue fixture: {error}"),
        };
        assert_eq!(actual, native);
    }
}
