//! Runner for unit `blend.fillet_solver_p1`. C++: `solidcore_ffi::units::blend_fillet_solver_p1`.
//! Rust: `solidcore::blend::fillet_solver::p1`.
//!
//! Both sides build the same solver: the base `SmFilletSolver` (the shim derives a
//! class that only names itself, so every virtual is the base version) with the
//! same offset surfaces (offsets of `SmPlane`s or `SmBSplineSurface`s, blocks of unit
//! `topology.edgeuse_p1`). Ops: the constructor, getters and an offset-bookkeeping
//! script (`state`), `SetupOffsetExtension` over planes, open, pinched (singular side)
//! and periodic (rational circle x line) B-spline surfaces at extreme scales
//! (`extension`), the sphere-solve and rail-UV-curve function objects with and without
//! Jacobians (`sphere`, `rail_uv`; the base `LoadJacobian` fails, after the rail-UV
//! rows are computed), `ComputeSurfaceValues` (base `SetupOffsetValues` fails),
//! `ComputeUVVectors` (including zero and parallel tangents), `SplitFilletGeom` with
//! empty intervals and `FinishConstruction` (deferred with an offset surface to copy:
//! `SmOffsetSurface::Copy`, unit `surface.offset_surface`).
//!
//! Not reachable from the shim: the file-static helpers (`sm_TestAnalyticSurface`,
//! `sm_CheckEdgeCurveDomains`, ...) and the `.cpp`-local `SmRailRailIntersectENFO` /
//! `SmFindClippedRadiusEFO`; the B-rep queries need edge uses the base solver does not
//! have (their `NULL` paths are compared in `state`).

use serde::{Deserialize, Serialize};
use serde_json::Value;

use solidcore::blend::fillet_custom_solver::{JacobianLoad, OffsetSurface, SurfaceId};
use solidcore::nurbs::math::NOW;
use solidcore::blend::fillet_solver::p1::{
    compute_uv_vectors, FilStatus, FilletSolver, FilletSolverOffsetsData, FilletSolverOps,
    FilletSphereSolveEnfo, OffsetSurfaceSlot, RailUvCurveIntersectEnfo,
};
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::matrix::Matrix;
use solidcore::geom::vector2d::{Point2, Vec2};
use solidcore::geom::vector3d::Vec3;
use solidcore::intersect::curve_class::p1::CurveClassification;
use solidcore::intersect::point_class::{ClassObject, ClassObjectKind};
use solidcore::intersect::surface_intersector::p1::TsectPnt;
use solidcore::kernel::Kernel;
use solidcore::topology::Brep;
use solidcore::SolidError;
use solidcore_ffi::units::blend_fillet_solver_p1::{self as ffi, Op};

use super::topology_edgeuse_p1::{
    bspline_of, domain, gen_curve, gen_surface, push_curve, push_surface,
    surface_of, BssIn, CrvIn, SrfIn,
};
use crate::{HarnessError, Output, Rng, UnitDiff};

/// The runner.
pub struct Runner;

fn harness(msg: impl Into<String>) -> HarnessError {
    HarnessError::new(msg)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Arrays {
    x: Vec<f64>,
    f: Vec<f64>,
    jac: Option<[usize; 2]>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum Inputs {
    State {
        tols: [f64; 3],
        base: [SrfIn; 2],
        has: [bool; 2],
        off: [f64; 2],
        orient: [f64; 2],
        cmds: Vec<[f64; 3]>,
    },
    Extension {
        ext: f64,
        surface: SrfIn,
        dom: [f64; 4],
    },
    Sphere {
        c: [f64; 3],
        r: f64,
        adjusts: bool,
        arrays: Arrays,
    },
    RailUv {
        rail: usize,
        #[serde(default)]
        load: bool,
        #[serde(default)]
        found: bool,
        curve: CrvIn,
        ivl: [f64; 2],
        arrays: Arrays,
    },
    SurfValues {
        uv: [f64; 4],
    },
    UvVectors {
        vecs: [[f64; 3]; 5],
        uvd: [f64; 4],
    },
    Split,
    Finish {
        surface: SrfIn,
        has: [bool; 2],
        off: f64,
    },
}

impl Inputs {
    fn ffi(&self) -> Op {
        match self {
            Inputs::State { .. } => Op::State,
            Inputs::Extension { .. } => Op::Extension,
            Inputs::Sphere { .. } => Op::Sphere,
            Inputs::RailUv { .. } => Op::RailUv,
            Inputs::SurfValues { .. } => Op::SurfValues,
            Inputs::UvVectors { .. } => Op::UvVectors,
            Inputs::Split => Op::Split,
            Inputs::Finish { .. } => Op::Finish,
        }
    }
}

fn push_arrays(o: &mut Vec<f64>, a: &Arrays) {
    o.push(a.x.len() as f64);
    o.extend_from_slice(&a.x);
    o.push(a.f.len() as f64);
    o.extend_from_slice(&a.f);
    match a.jac {
        Some([r, c]) => o.extend_from_slice(&[1.0, r as f64, c as f64]),
        None => o.extend_from_slice(&[0.0, 0.0, 0.0]),
    }
}

fn cpp_input(i: &Inputs) -> Result<Vec<f64>, HarnessError> {
    let mut o = Vec::new();
    match i {
        Inputs::State {
            tols,
            base,
            has,
            off,
            orient,
            cmds,
        } => {
            o.extend_from_slice(tols);
            push_surface(&mut o, &base[0])?;
            push_surface(&mut o, &base[1])?;
            for (h, d) in has.iter().zip(off) {
                o.push(if *h { 1.0 } else { 0.0 });
                o.push(*d);
            }
            o.extend_from_slice(orient);
            o.push(cmds.len() as f64);
            for c in cmds {
                o.extend_from_slice(c);
            }
        }
        Inputs::Extension { ext, surface, dom } => {
            o.push(*ext);
            push_surface(&mut o, surface)?;
            o.extend_from_slice(dom);
        }
        Inputs::Sphere {
            c,
            r,
            adjusts,
            arrays,
        } => {
            o.extend_from_slice(c);
            o.push(*r);
            o.push(if *adjusts { 1.0 } else { 0.0 });
            push_arrays(&mut o, arrays);
        }
        Inputs::RailUv {
            rail,
            load,
            found,
            curve,
            ivl,
            arrays,
        } => {
            o.push(*rail as f64);
            o.push(b(*load));
            o.push(b(*found));
            push_curve(&mut o, curve)?;
            o.extend_from_slice(ivl);
            push_arrays(&mut o, arrays);
        }
        Inputs::SurfValues { uv } => o.extend_from_slice(uv),
        Inputs::UvVectors { vecs, uvd } => {
            for v in vecs {
                o.extend_from_slice(v);
            }
            o.extend_from_slice(uvd);
        }
        Inputs::Split => o.push(0.0),
        Inputs::Finish { surface, has, off } => {
            push_surface(&mut o, surface)?;
            o.push(if has[0] { 1.0 } else { 0.0 });
            o.push(if has[1] { 1.0 } else { 0.0 });
            o.push(*off);
        }
    }
    Ok(o)
}

fn cpp(i: &Inputs) -> Result<Output, HarnessError> {
    let input = cpp_input(i)?;
    let flat = ffi::eval(i.ffi(), &input).map_err(|e| harness(e.to_string()))?;
    Ok(Output::new().values("flat", flat))
}

// ---- Rust side ------------------------------------------------------------------------

fn status(e: &SolidError) -> f64 {
    match e {
        SolidError::InvalidInput { .. } => 1007.0,
        SolidError::Upstream { code, .. } => f64::from(*code),
        _ => -1.0,
    }
}

fn b(v: bool) -> f64 {
    if v {
        1.0
    } else {
        0.0
    }
}

fn fil_status(v: f64) -> FilStatus {
    match v as i64 {
        1 => FilStatus::Processed,
        2 => FilStatus::Failure,
        3 => FilStatus::SurfIntFailure,
        4 => FilStatus::EndPntDropFailure,
        5 => FilStatus::ToBeSqueezed,
        6 => FilStatus::CubicRailRailInterpolate,
        7 => FilStatus::NoIntRailEu,
        8 => FilStatus::SolverNotConverge,
        9 => FilStatus::ToBeDeleted,
        10 => FilStatus::SurfTracingFailure,
        11 => FilStatus::SelfIntHandlingFailure,
        12 => FilStatus::SurfSkinningFailure,
        _ => FilStatus::Unprocessed,
    }
}

fn put_state(o: &mut Vec<f64>, fs: &FilletSolver) -> Result<(), SolidError> {
    for i in 0..2 {
        let s = fs.surface(i)?;
        o.push(b(s.is_some()));
        o.push(s.map_or(-1.0, |s| s.key as f64));
        o.push(s.map_or(0.0, |s| s.surface.offset_distance));
        o.push(fs.orientation(i)?);
        o.push(fs.saved_offsets.get(i).copied().unwrap_or(0.0));
    }
    o.push(fs.extension_factor());
    o.push(b(fs.reverse_trim()));
    o.push(f64::from(fs.status() as i32));
    Ok(())
}

fn rust_state(
    tols: &[f64; 3],
    base: &[SrfIn; 2],
    has: [bool; 2],
    off: &[f64; 2],
    orient: &[f64; 2],
    cmds: &[[f64; 3]],
    o: &mut Vec<f64>,
) -> Result<Result<(), SolidError>, HarnessError> {
    let bases = [surface_of(&base[0])?, surface_of(&base[1])?];
    let mut fs = FilletSolver::new(tols[0], tols[1], tols[2]);
    let mut next_key = 0u64;
    let mut make = |rail: usize, d: f64| {
        let key = next_key;
        next_key += 1;
        let base = if rail == 0 { &bases[0] } else { &bases[1] };
        OffsetSurfaceSlot {
            key,
            surface: OffsetSurface::new(d, base.clone(), SurfaceId::Free(rail as u32)),
        }
    };
    let r = (|| -> Result<(), SolidError> {
        for (k, (h, d)) in has.iter().zip(off).enumerate() {
            if *h {
                fs.set_surface(k, Some(make(k, *d)))?;
            }
        }
        fs.orientations = *orient;
        o.extend_from_slice(&[
            fs.this_approx_tol_3d(),
            fs.this_ang_tol_rad(),
            fs.tangency_tolerance(),
            fs.conversion_tol(),
            fs.extension_factor(),
            b(fs.reverse_trim()),
            f64::from(fs.status() as i32),
            fs.jacobian_size() as f64,
            b(fs.can_make_analytic()),
            b(fs.offset_radii_can_change()),
            fs.fillet_geoms().len() as f64,
            b(fs.first_fillet_geom()? == fs.last_fillet_geom()?),
            FilletSolver::TYPE as f64,
            b(FilletSolver::is_kind_of(FilletSolver::TYPE)),
            b(FilletSolver::is_kind_of(solidcore::geom::object::OBJECT_TYPE)),
            b(FilletSolver::is_kind_of(16002)),
            b(fs.edgeuse(0)?.is_none()),
            b(fs.vertex(&Brep::default(), 0)?.is_none()),
            b(fs.fillet_surface_generator().is_none()),
            b(fs.fillet_geom_for_surfaces().is_none()),
            b(fs.fillet_executive().is_none()),
        ]);
        let mut saved = FilletSolverOffsetsData::default();
        for c in cmds {
            match c[0] as i64 {
                0 => fs.set_surfaces_to_zero_offset()?,
                1 => fs.reload_surface_offsets()?,
                2 => saved = fs.offsets_data()?,
                3 => fs.set_offsets_data(&saved)?,
                4 => {
                    let rail = c[1] as usize;
                    fs.set_surface(rail, Some(make(rail, c[2])))?;
                }
                5 => fs.set_extension_factor(c[1]),
                6 => fs.set_reverse_trim(c[1] != 0.0),
                7 => fs.set_status(fil_status(c[1])),
                _ => {}
            }
            put_state(o, &fs)?;
        }
        Ok(())
    })();
    Ok(r)
}

fn put_arrays(o: &mut Vec<f64>, f: &[f64], j: Option<&Matrix>) {
    o.extend_from_slice(f);
    if let Some(j) = j {
        for r in 0..j.num_rows() {
            for c in 0..j.num_columns() {
                o.push(j.at(r, c).unwrap_or(f64::NAN));
            }
        }
    }
}

fn rust(i: &Inputs) -> Result<Output, HarnessError> {
    let mut o = Vec::new();
    let res: Result<(), SolidError> = match i {
        Inputs::State {
            tols,
            base,
            has,
            off,
            orient,
            cmds,
        } => rust_state(tols, base, *has, off, orient, cmds, &mut o)?,
        Inputs::Extension { ext, surface, dom } => {
            let k = Kernel::default();
            let s = surface_of(surface)?;
            let mut fs = FilletSolver::new(1.0e-6, 0.1, 0.01);
            fs.set_extension_factor(*ext);
            let mut d = Extent2d::from_points(Point2::new(dom[0], dom[1]), Point2::new(dom[2], dom[3]));
            let r = fs.setup_offset_extension(&k, &s, &mut d);
            if r.is_ok() {
                let (u, v) = (d.u_interval(), d.v_interval());
                o.extend_from_slice(&[0.0, u.min(), v.min(), u.max(), v.max()]);
            }
            r
        }
        Inputs::Sphere {
            c,
            r,
            adjusts,
            arrays,
        } => {
            let mut fs = FilletSolver::new(1.0e-6, 0.1, 0.01);
            let mut f = arrays.f.clone();
            let mut j = arrays.jac.map(|[r, c]| Matrix::new(r, c));
            let mut e = FilletSphereSolveEnfo::new(Vec3::new(c[0], c[1], c[2]), *r, &mut fs);
            e.set_plane_adjusts(*adjusts);
            let res = e.evaluate(&arrays.x, &mut f, j.as_mut());
            match &res {
                Ok(found) => o.extend_from_slice(&[0.0, b(*found)]),
                Err(e) => {
                    if let Some(out) = Output::unported(e) {
                        return Ok(out);
                    }
                    o.extend_from_slice(&[status(e), 0.0]);
                }
            }
            put_arrays(&mut o, &f, j.as_ref());
            Ok(())
        }
        Inputs::RailUv {
            rail,
            load,
            found,
            curve,
            ivl,
            arrays,
        } => {
            let crv = bspline_of(curve)?;
            let mut fs = ShimSolver {
                base: FilletSolver::new(1.0e-6, 0.1, 0.01),
                load: *load,
                found: *found,
            };
            let mut f = arrays.f.clone();
            let mut j = arrays.jac.map(|[r, c]| Matrix::new(r, c));
            let mut e = RailUvCurveIntersectEnfo::new(
                &crv,
                Extent1d::from_min_max(ivl[0], ivl[1]),
                *rail,
                &mut fs,
            );
            let res = e.evaluate(&arrays.x, &mut f, j.as_mut());
            match &res {
                Ok(found) => o.extend_from_slice(&[0.0, b(*found)]),
                Err(e) => {
                    if let Some(out) = Output::unported(e) {
                        return Ok(out);
                    }
                    o.extend_from_slice(&[status(e), 0.0]);
                }
            }
            put_arrays(&mut o, &f, j.as_ref());
            Ok(())
        }
        Inputs::SurfValues { uv } => {
            let mut fs = FilletSolver::new(1.0e-6, 0.1, 0.01);
            let mut t = TsectPnt::new();
            let r = FilletSolver::compute_surface_values(
                &mut fs as &mut dyn FilletSolverOps,
                [Point2::new(uv[0], uv[1]), Point2::new(uv[2], uv[3])],
                &mut t,
            );
            if r.is_ok() {
                o.push(0.0);
            }
            r
        }
        Inputs::UvVectors { vecs, uvd } => {
            let v = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
            let mut t = TsectPnt::new();
            let [du0, dv0, du1, dv1, cd] = *vecs;
            let [a0, b0, a1, b1] = *uvd;
            t.surface_pv[0][1][0] = v(du0);
            t.surface_pv[0][0][1] = v(dv0);
            t.surface_pv[1][1][0] = v(du1);
            t.surface_pv[1][0][1] = v(dv1);
            t.uv_curve_pv[0][1] = Vec2::new(a0, b0);
            t.uv_curve_pv[1][1] = Vec2::new(a1, b1);
            t.curve_pv[1] = v(cd);
            let r = compute_uv_vectors(&mut t);
            o.push(match &r {
                Ok(()) => 0.0,
                Err(e) => status(e),
            });
            for uvp in &t.uv_curve_pv {
                o.push(uvp[1].x);
                o.push(uvp[1].y);
            }
            match r {
                Err(e @ SolidError::Unported { .. }) => Err(e),
                _ => Ok(()),
            }
        }
        Inputs::Split => {
            let k = Kernel::default();
            let mut fs = FilletSolver::new(1.0e-6, 0.1, 0.01);
            let cc: CurveClassification<H> = CurveClassification::new(&k);
            let r = fs.split_fillet_geom(0, 0, &[], &[], &cc, &cc);
            if r.is_ok() {
                o.extend_from_slice(&[0.0, fs.fillet_geoms().len() as f64]);
            }
            r
        }
        Inputs::Finish { surface, has, off } => {
            let base = surface_of(surface)?;
            let mut fs = FilletSolver::new(1.0e-6, 0.1, 0.01);
            let r = (|| -> Result<(), SolidError> {
                if has[0] {
                    let s = OffsetSurface::new(*off, base.clone(), SurfaceId::Free(0));
                    fs.set_surface(0, Some(OffsetSurfaceSlot { key: 0, surface: s }))?;
                }
                if has[1] {
                    let s = OffsetSurface::new(-*off, base.clone(), SurfaceId::Free(0));
                    fs.set_surface(1, Some(OffsetSurfaceSlot { key: 1, surface: s }))?;
                }
                fs.finish_construction()
            })();
            if r.is_ok() {
                o.push(0.0);
            }
            r
        }
    };
    match res {
        Ok(()) => Ok(Output::new().values("flat", o)),
        Err(e) => {
            if let Some(out) = Output::unported(&e) {
                return Ok(out);
            }
            if o.is_empty() || matches!(i, Inputs::State { .. }) {
                o.push(status(&e));
            }
            Ok(Output::new().values("flat", o))
        }
    }
}

/// The shim's solver: the base, or with `load` the deterministic `LoadJacobian` of
/// `ShimFS` (row `ne` = `x[surf1_offset] - 0.5 x[0] + rail` and its partials, written
/// only inside `f` and the Jacobian; `surf2_offset = np`; `ne += 1`, `np += 2`).
struct ShimSolver {
    base: FilletSolver,
    load: bool,
    found: bool,
}

impl FilletSolverOps for ShimSolver {
    fn solver(&self) -> &FilletSolver {
        &self.base
    }

    fn solver_mut(&mut self) -> &mut FilletSolver {
        &mut self.base
    }

    fn load_jacobian(
        &mut self,
        num_equations: &mut usize,
        num_parameters: &mut usize,
        load: &mut JacobianLoad<'_>,
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
    ) -> Result<bool, SolidError> {
        if !self.load {
            return self.base.load_jacobian(num_equations, num_parameters, load, f, jacobian);
        }
        let ne = *num_equations;
        let o1 = load.surf1_offset;
        if let (Some(slot), Some(&xo), Some(&x0)) = (f.get_mut(ne), load.x.get(o1), load.x.first()) {
            *slot = xo - 0.5 * x0 + load.rail_index as f64;
        }
        if let Some(j) = jacobian {
            if ne < j.num_rows() {
                if o1 < j.num_columns() {
                    j.set_at(ne, o1, 1.0)?;
                }
                if j.num_columns() > 0 {
                    j.set_at(ne, 0, -0.5)?;
                }
            }
        }
        load.surf2_offset = *num_parameters;
        *num_equations += 1;
        *num_parameters += 2;
        Ok(self.found)
    }
}

/// A classification handle for the (empty) `SplitFilletGeom` classifications.
#[derive(Debug, Clone, Copy, PartialEq)]
struct H;

impl ClassObject for H {
    fn class_kind(&self) -> ClassObjectKind {
        ClassObjectKind::Other
    }
}

// ---- generators ------------------------------------------------------------------------

/// A degree-2 rational full circle of radius `r` (u) swept along z (v, degree 1): a
/// cylinder periodic in u.
fn cylinder(r: f64, h: f64) -> SrfIn {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let ring: [[f64; 3]; 9] = [
        [1.0, 0.0, 1.0],
        [1.0, 1.0, w],
        [0.0, 1.0, 1.0],
        [-1.0, 1.0, w],
        [-1.0, 0.0, 1.0],
        [-1.0, -1.0, w],
        [0.0, -1.0, 1.0],
        [1.0, -1.0, w],
        [1.0, 0.0, 1.0],
    ];
    let net = ring
        .iter()
        .map(|&[x, y, wt]| {
            vec![
                [x * r * wt, y * r * wt, 0.0, wt],
                [x * r * wt, y * r * wt, h * wt, wt],
            ]
        })
        .collect();
    SrfIn::Bspline(BssIn {
        oob: false,
        form: 0,
        net,
        p: 2,
        q: 1,
        u: vec![0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0],
        v: vec![0.0, 0.0, 1.0, 1.0],
    })
}

fn gen_srf(rng: &mut Rng) -> SrfIn {
    if rng.chance(0.25) {
        let scale = *rng.choose(&[1e-4, 1.0, 1.0, 1e4]).unwrap_or(&1.0);
        return cylinder(scale * rng.uniform(0.5, 2.0), scale * rng.uniform(0.5, 3.0));
    }
    gen_surface(rng)
}

fn gen_arrays(rng: &mut Rng, min_x: usize, min_f: usize, min_rows: usize, min_cols: usize) -> Arrays {
    let nx = min_x + rng.below(4) as usize;
    let x = (0..nx).map(|_| rng.uniform(-1.0, 2.0)).collect();
    let jac = rng.chance(0.6).then(|| {
        let rows = min_rows + rng.below(4) as usize;
        let cols = min_cols + rng.below(4) as usize;
        [rows, cols]
    });
    let nf = jac.map_or(min_f, |[r, _]| r.max(min_f)) + rng.below(3) as usize;
    let f = (0..nf).map(|_| rng.uniform(-5.0, 5.0)).collect();
    Arrays { x, f, jac }
}

fn gen_state(rng: &mut Rng) -> Inputs {
    let n = rng.below(12) as usize;
    let cmds = (0..n)
        .map(|_| match rng.below(10) {
            0 | 1 => [0.0, 0.0, 0.0],
            2 => [1.0, 0.0, 0.0],
            3 | 4 => [2.0, 0.0, 0.0],
            5 | 6 => [3.0, 0.0, 0.0],
            7 => [4.0, rng.below(2) as f64, rng.uniform(-3.0, 3.0)],
            8 => [5.0 + rng.below(2) as f64, rng.uniform(0.0, 2.0).round(), 0.0],
            _ => [7.0, rng.below(13) as f64, 0.0],
        })
        .collect();
    let o = |rng: &mut Rng| if rng.chance(0.5) { 1.0 } else { -1.0 };
    Inputs::State {
        tols: [rng.uniform(1e-8, 1e-3), rng.uniform(0.01, 0.5), rng.uniform(0.001, 0.1)],
        base: [gen_surface(rng), gen_surface(rng)],
        has: [rng.chance(0.8), rng.chance(0.8)],
        off: [rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0)],
        orient: [o(rng), o(rng)],
        cmds,
    }
}

fn gen_extension(rng: &mut Rng) -> Inputs {
    let surface = gen_srf(rng);
    let [u0, v0, u1, v1] = match &surface {
        SrfIn::Plane { .. } => {
            let s = *rng.choose(&[1e-6, 1.0, 1e6]).unwrap_or(&1.0);
            let a = rng.uniform(-1.0, 1.0) * s;
            let c = rng.uniform(-1.0, 1.0) * s;
            [a, c, a + rng.uniform(0.0, 2.0) * s, c + rng.uniform(0.0, 2.0) * s]
        }
        s @ SrfIn::Bspline(_) => {
            let d = domain(s);
            if rng.chance(0.7) {
                d
            } else {
                // a sub-domain
                let t = |rng: &mut Rng, a: f64, b: f64| a + (b - a) * rng.uniform(0.0, 1.0);
                let (ua, ub) = (t(rng, d[0], d[2]), t(rng, d[0], d[2]));
                let (va, vb) = (t(rng, d[1], d[3]), t(rng, d[1], d[3]));
                [ua.min(ub), va.min(vb), ua.max(ub), va.max(vb)]
            }
        }
    };
    let ext = *rng.choose(&[0.5, 0.5, 0.0, 0.1, 2.0, -0.3]).unwrap_or(&0.5);
    Inputs::Extension {
        ext,
        surface,
        dom: [u0, v0, u1, v1],
    }
}

fn gen_case(rng: &mut Rng) -> Inputs {
    match rng.below(20) {
        0..=4 => gen_state(rng),
        5..=10 => gen_extension(rng),
        11 => Inputs::Sphere {
            c: rng.point3(-2.0, 2.0),
            r: rng.uniform(-2.0, 2.0),
            adjusts: rng.chance(0.5),
            arrays: gen_arrays(rng, 0, 0, 0, 0),
        },
        12..=14 => {
            let scale = *rng.choose(&[1e-3, 1.0, 1e3]).unwrap_or(&1.0);
            let dim = 2 + rng.below(2);
            let curve = gen_curve(rng, dim, scale);
            let a = curve.u.first().copied().unwrap_or(0.0);
            let z = curve.u.last().copied().unwrap_or(1.0);
            let mut arrays = gen_arrays(rng, 3, 2, 2, 3);
            let load = rng.chance(0.7);
            if let Some(x0) = arrays.x.first_mut() {
                *x0 = a + (z - a) * rng.uniform(-0.2, 1.2);
            }
            // half the loading cases start on the curve's first point (a clamped
            // non-rational curve passes through it): the converged branch
            if let (true, Some([px, py, _, w])) = (load && rng.chance(0.5), curve.cpts.first().copied()) {
                if w == NOW {
                    if let [x0, x1, x2, ..] = arrays.x.as_mut_slice() {
                        *x0 = a;
                        *x1 = px;
                        *x2 = py;
                    }
                }
            }
            Inputs::RailUv {
                load,
                found: rng.chance(0.6),
                rail: rng.below(2) as usize,
                curve,
                ivl: [a, z],
                arrays,
            }
        }
        15 => Inputs::SurfValues {
            uv: [rng.uniform(0.0, 1.0), rng.uniform(0.0, 1.0), rng.uniform(0.0, 1.0), rng.uniform(0.0, 1.0)],
        },
        16 | 17 => {
            let mut vecs = [[0.0; 3]; 5];
            for v in &mut vecs {
                *v = rng.point3(-2.0, 2.0);
            }
            match rng.below(5) {
                0 => vecs[0] = [0.0, 0.0, 0.0],
                1 => vecs[1] = vecs[0].map(|c| c * 2.0),
                2 => vecs[2] = vecs[2].map(|c| c * 1e-5),
                _ => {}
            }
            Inputs::UvVectors {
                vecs,
                uvd: [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0)],
            }
        }
        18 => Inputs::Split,
        _ => Inputs::Finish {
            surface: gen_surface(rng),
            has: [rng.chance(0.5), rng.chance(0.5)],
            off: rng.uniform(-1.0, 1.0),
        },
    }
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "blend.fillet_solver_p1"
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        serde_json::to_value(gen_case(rng)).ok()
    }

    fn cases(&self, ctx: &crate::case::CaseContext) -> Result<Vec<crate::Case>, HarnessError> {
        let ctx = crate::case::CaseContext {
            seeds: ctx.seeds.max(200),
            ..ctx.clone()
        };
        crate::runner::default_cases(self, &ctx)
    }

    fn run_cpp(&self, case: &crate::Case) -> Result<Output, HarnessError> {
        cpp(&case.parse()?)
    }

    fn run_rust(&self, case: &crate::Case) -> Result<Output, HarnessError> {
        rust(&case.parse()?)
    }
}
