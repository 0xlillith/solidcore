//! FFI wrappers for unit `surface.offset_surface`.
//!
//! Shim: `harness/cpp/shim/surface_offset_surface.cpp` (`sc_surface_offset_surface_eval`).
//! Upstream: `source/SMLib/src/SmOffsetSurface.cpp`, `source/SMLib/inc/SmOffsetSurface.h`.
//!
//! An offset block is `d owns sing allow base`: `SmOffsetSurface(d, *base, owns, sing)`
//! then `SetAllowExtension(allow)`. A base block starts with its kind:
//!
//! | kind | base                                                    | data                                  |
//! |------|---------------------------------------------------------|---------------------------------------|
//! | `0`  | `SmBSplineSurface(gw_SURFACE*, FALSE, NULL)` (+ oob, form) | `oob form n m p q r s cpts(x y z w) u v` |
//! | `1`  | `SmPlane(O, N)`                                         | `O[3] N[3]`                           |
//! | `2`  | `SmCurveBoundedSurface` without boundary curves         | basis block (kind 0 or 1), `implicit` |
//! | `3`  | `SmSTEPSurface(*base, TRUE)`                            | base block (kind 0, 1 or 4)           |
//! | `4`  | `SmSurfOfRevolution` (unit `surface.surf_of_revolution_p1`'s block) | that block            |
//! | `5`  | a nested offset                                         | an offset block                       |
//!
//! [`Op::Script`] builds the offset (the record starts with the construction status:
//! a surface of revolution can fail to build) and runs its ops on it in order; every op
//! writes its status (`0` = `SM_SUCCESS`, else the `SmStatus`) then its outputs:
//!
//! | [`Step`]       | upstream                              | outputs on success                     |
//! |----------------|---------------------------------------|----------------------------------------|
//! | `State`        | the members, `GetNaturalUVDomain`, a 3x3 `EvaluatePoint` grid | see [`Step::State`] |
//! | `Eval`         | `Evaluate`                            | `(hu+1)(hv+1)` vectors, untouched ones `SM_UNDEF_DOUBLE` |
//! | `EvalSimple`   | `EvaluateSimple`                      | as `Eval`                              |
//! | `EvalPt`       | `EvaluatePoint`                       | `P`                                    |
//! | `Normal`       | `EvaluateNormal`                      | `N`                                    |
//! | `Knots`        | `GetKnots`                            | `n knots mults`                        |
//! | `MaxAnal`      | `GetMaxAnalyticDomain`                | `D`                                    |
//! | `BBox`         | `CalculateBoundingBox`                | the wanted boxes                       |
//! | `StepInv`      | `STEPInversion`                       | `u v loc`                              |
//! | `Reverse`, `Swap`, `Transform`, `SetDist`, `SetAllow` | `Reverse`, `SwapUV`, `Transform`, `SetOffsetDistance`, `SetAllowExtension` | - |
//! | `Trim`         | `TrimWithDomain`                      | the (snapped) domain                   |
//! | `Split`        | `SplitAt`                             | both halves' states                    |
//! | `CreateOffset` | `CreateOffsetSurface`                 | the result's state                     |
//! | `Extend`, `ExtendDir` | `CreateExtendedSurface`        | the result's state                     |
//! | `Copy`         | `Copy`                                | the copy's state                       |
//! | `Eq`           | `operator==` with another offset and itself | `eq self_eq`                     |
//! | `IsPlanar`     | `IsPlanar` (`SmSurface` body, virtual box) | `planar P N dev`                  |
//! | `IsKindOf`     | `IsKindOf`                            | `flag`                                 |
//!
//! [`script`] accepts only what upstream handles without undefined behaviour or runaway
//! work: well-formed clamped B-splines, planes with a non-zero normal, offsets nested at
//! most [`MAX_NESTING`] deep, finite values below [`MAX_COORD`], parameters within
//! [`MAX_REACH`] domain lengths of the base domain, derivative orders of at most 6 and
//! never written past the caller's array (`bOnlyUpperHalf` with `hv > hu` writes
//! `Duv` past it), no bounding boxes or planarity tests over a curve-bounded or STEP
//! base (their polar box recurses without bound upstream, UB-84), no box whose domain
//! may lie outside the base domain (the base-class box of the offset asks the offset
//! again for the outside part, without bound once an extension exists, and its surface
//! cache subdivision runs away for small trimmed bases; the Rust side defers that path
//! on `surface.surface_p6`), and no planarity test once an extension may exist.

use super::surface_surf_of_revolution_p1 as rev;
use crate::FfiError;

/// Unchanged C++ did not finish within the five-second oracle deadline.
/// This is a harness error, not a kernel status or a deferred comparison.
pub const REFERENCE_TIMEOUT: i32 = crate::offset_reference::TIMEOUT;
/// The isolated reference process crashed or its transport failed.
pub const REFERENCE_PROCESS_ERROR: i32 = crate::offset_reference::PROCESS_ERROR;

/// Largest control point index.
pub const MAX_INDEX: usize = 24;
/// Largest degree.
pub const MAX_DEGREE: usize = 5;
/// Largest coordinate magnitude.
pub const MAX_COORD: f64 = 1.0e6;
/// Deepest offset nesting in a block.
pub const MAX_NESTING: usize = 3;
/// Largest number of ops in a script.
pub const MAX_OPS: usize = 16;
/// How many domain lengths outside the base domain a parameter may be.
pub const MAX_REACH: f64 = 3.0;
/// Largest `WriteToDB` record [`read_db`] passes on.
pub const MAX_DB_BYTES: usize = 1 << 20;

const NOW: f64 = -f64::MAX;

const OP_SCRIPT: i32 = 0;
const OP_WRITE_DB: i32 = 1;
const OP_READ_DB: i32 = 2;

/// An `SmBSplineSurface` over a NURB net.
#[derive(Debug, Clone, PartialEq)]
pub struct Bss {
    /// `SetOutOfBoundsEnabled`.
    pub oob: bool,
    /// `SetBSplineSurfaceForm` code.
    pub form: u8,
    /// Control points `x y z w` (`w` = `NL_NOW` for a polynomial surface), `net[i][j]`.
    pub net: Vec<Vec<[f64; 4]>>,
    /// Degree in u.
    pub p: usize,
    /// Degree in v.
    pub q: usize,
    /// Knots in u.
    pub u: Vec<f64>,
    /// Knots in v.
    pub v: Vec<f64>,
}

/// A base surface.
#[derive(Debug, Clone, PartialEq)]
pub enum Srf {
    /// Kind 0.
    BSpline(Bss),
    /// Kind 1: `SmPlane(o, n)`.
    Plane {
        /// Origin.
        o: [f64; 3],
        /// Normal.
        n: [f64; 3],
    },
    /// Kind 2: a curve-bounded surface without boundary curves.
    Bounded {
        /// Basis (a B-spline or a plane).
        basis: Box<Srf>,
        /// `bImplicitOuter`.
        implicit: bool,
    },
    /// Kind 3: an owning `SmSTEPSurface`.
    Step {
        /// Base (a B-spline, a plane or a surface of revolution).
        base: Box<Srf>,
    },
    /// Kind 4: a surface of revolution.
    Revolution(rev::Surf),
    /// Kind 5: a nested offset.
    Offset(Box<Off>),
}

/// An offset surface to build.
#[derive(Debug, Clone, PartialEq)]
pub struct Off {
    /// `dOffsetDistance`.
    pub d: f64,
    /// `bOwnsSurface`.
    pub owns: bool,
    /// `lSingularities` (`SM_SS_*` flags, `16` = `SM_SS_UNKNOWN`).
    pub sing: u32,
    /// `SetAllowExtension`.
    pub allow: bool,
    /// The base surface.
    pub base: Srf,
}

/// One op of a script.
#[derive(Debug, Clone, PartialEq)]
pub enum Step {
    /// `dom[4] d allow owns per_u per_v sing has_ext [ext_dom[4] ext_type] base_type
    /// base_dom[4] type`, then for a 3x3 grid over the natural domain each
    /// `EvaluatePoint`'s status and point.
    State,
    /// `Evaluate(uv, hu, hv, ul, vl, upper, .., nzt)`.
    Eval {
        /// Parameter.
        uv: [f64; 2],
        /// Highest u derivative.
        hu: usize,
        /// Highest v derivative.
        hv: usize,
        /// `bUFromLeft`.
        ul: bool,
        /// `bVFromLeft`.
        vl: bool,
        /// `bOnlyUpperHalf`.
        upper: bool,
        /// `bNonZeroTangents`.
        nzt: bool,
    },
    /// `EvaluateSimple(uv, hu, hv, ul, vl, upper, ..)`.
    EvalSimple {
        /// Parameter.
        uv: [f64; 2],
        /// Highest u derivative.
        hu: usize,
        /// Highest v derivative.
        hv: usize,
        /// `bUFromLeft`.
        ul: bool,
        /// `bVFromLeft`.
        vl: bool,
        /// `bOnlyUpperHalf`.
        upper: bool,
    },
    /// `EvaluatePoint(uv)`.
    EvalPt {
        /// Parameter.
        uv: [f64; 2],
    },
    /// `EvaluateNormal(uv, ul, vl)`.
    Normal {
        /// Parameter.
        uv: [f64; 2],
        /// `bUFromLeft`.
        ul: bool,
        /// `bVFromLeft`.
        vl: bool,
    },
    /// `GetKnots(param, .., ivl)`.
    Knots {
        /// `0` = `SM_SP_U`, `1` = `SM_SP_V`.
        param: u8,
        /// Interval of interest.
        ivl: Option<[f64; 2]>,
    },
    /// `GetMaxAnalyticDomain`.
    MaxAnal,
    /// `CalculateBoundingBox(dom, normal?, pseudo?, polar?, NULL, NULL, expand)`.
    BBox {
        /// Domain.
        dom: [f64; 4],
        /// Normal, pseudo and polar box wanted.
        want: [bool; 3],
        /// `bExpandPosBoxesByZoneTol3d`.
        expand: bool,
    },
    /// `STEPInversion(dom, p, tol, uv, loc, guess)`.
    StepInv {
        /// Analytic domain.
        dom: [f64; 4],
        /// Point.
        p: [f64; 3],
        /// Distance tolerance.
        tol: f64,
        /// Guess.
        guess: Option<[f64; 2]>,
    },
    /// `Reverse(param)`.
    Reverse {
        /// `0` = `SM_SP_U`, `1` = `SM_SP_V`.
        param: u8,
    },
    /// `SwapUV()`.
    Swap,
    /// `Transform(SmAxis2Placement(origin, x, y), scale)`.
    Transform {
        /// Origin.
        origin: [f64; 3],
        /// X axis.
        x: [f64; 3],
        /// Y axis.
        y: [f64; 3],
        /// Optional scale.
        scale: Option<[f64; 3]>,
    },
    /// `TrimWithDomain(dom)`.
    Trim {
        /// Domain.
        dom: [f64; 4],
    },
    /// `SplitAt(t, param)`.
    Split {
        /// Split parameter.
        t: f64,
        /// `0` = `SM_SP_U`, `1` = `SM_SP_V`.
        param: u8,
    },
    /// `CreateOffsetSurface(d, 0.0)`.
    CreateOffset {
        /// Signed distance.
        d: f64,
    },
    /// `CreateExtendedSurface(dist, cont)`.
    Extend {
        /// Distance.
        dist: f64,
        /// `SmContinuityType` code.
        cont: u8,
    },
    /// `CreateExtendedSurface(param, dist, cont)`.
    ExtendDir {
        /// `SmSurfParamType` code.
        param: u8,
        /// Distance.
        dist: f64,
        /// `SmContinuityType` code.
        cont: u8,
    },
    /// `Copy`.
    Copy,
    /// `operator==` with another offset, then with itself.
    Eq {
        /// The other offset.
        other: Box<Off>,
    },
    /// `SetOffsetDistance(d)`.
    SetDist {
        /// Distance.
        d: f64,
    },
    /// `SetAllowExtension(b)`.
    SetAllow {
        /// Flag.
        b: bool,
    },
    /// `IsPlanar(tol, ..)`.
    IsPlanar {
        /// Tolerance.
        tol: f64,
    },
    /// `IsKindOf(t)`.
    IsKindOf {
        /// Type code.
        t: i64,
    },
}

impl Step {
    fn code(&self) -> f64 {
        f64::from(match self {
            Step::State => 0,
            Step::Eval { .. } => 1,
            Step::EvalSimple { .. } => 2,
            Step::EvalPt { .. } => 3,
            Step::Normal { .. } => 4,
            Step::Knots { .. } => 5,
            Step::MaxAnal => 6,
            Step::BBox { .. } => 7,
            Step::StepInv { .. } => 8,
            Step::Reverse { .. } => 9,
            Step::Swap => 10,
            Step::Transform { .. } => 11,
            Step::Trim { .. } => 12,
            Step::Split { .. } => 13,
            Step::CreateOffset { .. } => 14,
            Step::Extend { .. } => 15,
            Step::ExtendDir { .. } => 16,
            Step::Copy => 17,
            Step::Eq { .. } => 18,
            Step::SetDist { .. } => 19,
            Step::SetAllow { .. } => 20,
            Step::IsPlanar { .. } => 21,
            Step::IsKindOf { .. } => 22,
        })
    }
}

fn bad(msg: impl Into<String>) -> FfiError {
    FfiError::InvalidArgument(msg.into())
}

fn b(x: bool) -> f64 {
    f64::from(u8::from(x))
}

fn finite(v: &[f64], what: &str) -> Result<(), FfiError> {
    if v.iter().any(|c| !c.is_finite() || c.abs() > MAX_COORD) {
        return Err(bad(format!("{what}: non-finite or huge value")));
    }
    Ok(())
}

fn nonzero(v: [f64; 3], what: &str) -> Result<(), FfiError> {
    finite(&v, what)?;
    if v.iter().map(|c| c * c).sum::<f64>() < 1e-12 {
        return Err(bad(format!("{what}: zero vector")));
    }
    Ok(())
}

fn domain(d: [f64; 4], what: &str) -> Result<(), FfiError> {
    finite(&d, what)?;
    let [u0, v0, u1, v1] = d;
    if u0 > u1 || v0 > v1 {
        return Err(bad(format!("{what}: inverted domain")));
    }
    Ok(())
}

/// Checks a clamped knot vector of degree `p`; returns its range.
fn check_knots(knots: &[f64], p: usize, what: &str) -> Result<(f64, f64), FfiError> {
    finite(knots, what)?;
    let m = knots
        .len()
        .checked_sub(1)
        .ok_or_else(|| bad(format!("{what}: no knots")))?;
    if knots.windows(2).any(|w| w.first() > w.last()) {
        return Err(bad(format!("{what}: knots not non-decreasing")));
    }
    let (Some(&a), Some(&b)) = (knots.first(), knots.last()) else {
        return Err(bad(format!("{what}: empty knots")));
    };
    let short = || bad(format!("{what}: knots too short"));
    if m < 2 * p + 1 {
        return Err(short());
    }
    let head = knots.get(..=p).ok_or_else(short)?;
    let tail = knots.get(m - p..).ok_or_else(short)?;
    let inner = knots.get(p + 1..m - p).ok_or_else(short)?;
    if a >= b || head.iter().any(|&k| k != a) || tail.iter().any(|&k| k != b) {
        return Err(bad(format!("{what}: knots not clamped")));
    }
    if b - a < 1e-6 * a.abs().max(b.abs()).max(1.0) {
        return Err(bad(format!("{what}: parameter range too small")));
    }
    if inner.iter().any(|&k| k == a || k == b) {
        return Err(bad(format!("{what}: end multiplicity above degree + 1")));
    }
    let mut run = 1;
    for w in inner.windows(2) {
        run = if w.first() == w.last() { run + 1 } else { 1 };
        if run > p {
            return Err(bad(format!(
                "{what}: interior multiplicity above the degree"
            )));
        }
    }
    Ok((a, b))
}

/// Checks a B-spline block; returns its natural domain.
fn check_bss(s: &Bss) -> Result<[f64; 4], FfiError> {
    if s.form > 11 {
        return Err(bad("surface: bad form"));
    }
    let n = s
        .net
        .len()
        .checked_sub(1)
        .ok_or_else(|| bad("surface: empty net"))?;
    let cols = s.net.first().map_or(0, Vec::len);
    let m = cols
        .checked_sub(1)
        .ok_or_else(|| bad("surface: empty net row"))?;
    if s.net.iter().any(|r| r.len() != cols) {
        return Err(bad("surface: ragged net"));
    }
    if n > MAX_INDEX || m > MAX_INDEX || s.p > MAX_DEGREE || s.q > MAX_DEGREE {
        return Err(bad("surface: too large"));
    }
    if s.p == 0 || s.p > n || s.q == 0 || s.q > m {
        return Err(bad("surface: need 1 <= p <= n, 1 <= q <= m"));
    }
    if s.u.len() != n + s.p + 2 || s.v.len() != m + s.q + 2 {
        return Err(bad("surface: knot counts"));
    }
    let mut rational = None;
    for &[x, y, z, w] in s.net.iter().flatten() {
        finite(&[x, y, z], "surface point")?;
        let r = w != NOW;
        if *rational.get_or_insert(r) != r {
            return Err(bad("surface: mixed rationality"));
        }
        if r && !(1e-3..=1e3).contains(&w) {
            return Err(bad("surface: weight out of range"));
        }
    }
    let (u0, u1) = check_knots(&s.u, s.p, "surface u")?;
    let (v0, v1) = check_knots(&s.v, s.q, "surface v")?;
    Ok([u0, v0, u1, v1])
}

/// What the validator knows about a base: its natural domain (`None` when not known
/// in advance), whether a bounding box of it recurses upstream, and its distance.
#[derive(Debug, Clone, Copy)]
struct BaseInfo {
    dom: Option<[f64; 4]>,
    unboxable: bool,
}

fn push_bss(o: &mut Vec<f64>, s: &Bss) -> Result<[f64; 4], FfiError> {
    let dom = check_bss(s)?;
    let n = s.net.len() - 1;
    let m = s.net.first().map_or(0, Vec::len) - 1;
    o.extend([0.0, b(s.oob), f64::from(s.form)]);
    o.extend([n, m, s.p, s.q, s.u.len() - 1, s.v.len() - 1].map(|x| x as f64));
    for p in s.net.iter().flatten() {
        o.extend_from_slice(p);
    }
    o.extend_from_slice(&s.u);
    o.extend_from_slice(&s.v);
    Ok(dom)
}

fn push_srf(o: &mut Vec<f64>, s: &Srf, depth: usize) -> Result<BaseInfo, FfiError> {
    Ok(match s {
        Srf::BSpline(bs) => BaseInfo {
            dom: Some(push_bss(o, bs)?),
            unboxable: false,
        },
        Srf::Plane { o: org, n } => {
            finite(org, "plane origin")?;
            nonzero(*n, "plane normal")?;
            o.push(1.0);
            o.extend_from_slice(org);
            o.extend_from_slice(n);
            BaseInfo {
                dom: None,
                unboxable: false,
            }
        }
        Srf::Bounded { basis, implicit } => {
            if !matches!(**basis, Srf::BSpline(_) | Srf::Plane { .. }) {
                return Err(bad("curve-bounded basis: a B-spline or a plane"));
            }
            o.push(2.0);
            let info = push_srf(o, basis, depth)?;
            o.push(b(*implicit));
            BaseInfo {
                dom: info.dom,
                unboxable: true,
            }
        }
        Srf::Step { base } => {
            if !matches!(
                **base,
                Srf::BSpline(_) | Srf::Plane { .. } | Srf::Revolution(_)
            ) {
                return Err(bad("STEP base: a B-spline, plane or surface of revolution"));
            }
            o.push(3.0);
            let info = push_srf(o, base, depth)?;
            BaseInfo {
                // (the STEP domain of a surface of revolution differs)
                dom: if matches!(**base, Srf::BSpline(_)) {
                    info.dom
                } else {
                    None
                },
                unboxable: true,
            }
        }
        Srf::Revolution(r) => {
            o.push(4.0);
            rev::push_surf(o, r)?;
            BaseInfo {
                dom: None,
                unboxable: false,
            }
        }
        Srf::Offset(off) => {
            if depth >= MAX_NESTING {
                return Err(bad("offsets nested too deep"));
            }
            o.push(5.0);
            push_off(o, off, depth + 1)?
        }
    })
}

fn push_off(o: &mut Vec<f64>, off: &Off, depth: usize) -> Result<BaseInfo, FfiError> {
    finite(&[off.d], "offset distance")?;
    if off.sing > 16 {
        return Err(bad("singularities: SM_SS_* flags or SM_SS_UNKNOWN"));
    }
    o.extend([off.d, b(off.owns), f64::from(off.sing), b(off.allow)]);
    push_srf(o, &off.base, depth)
}

/// The largest index `EvaluatePlainOrSimple` writes for a non-zero offset.
fn max_write_index(hu: usize, hv: usize, upper: bool) -> usize {
    let hod = if upper { hu.max(hv) } else { hu + hv };
    let mut idx = vec![0];
    if hu >= 1 {
        idx.push(hv + 1);
    }
    if hv >= 1 {
        idx.push(1);
    }
    if hu >= 2 {
        idx.push(2 * (hv + 1));
    }
    if (upper && hod >= 2) || (!upper && hu >= 1 && hv >= 1) {
        idx.push(hv + 2);
    }
    if hv >= 2 {
        idx.push(2);
    }
    if hu >= 3 {
        idx.push(3 * (hv + 1));
    }
    if (upper && hod >= 3) || (!upper && hu >= 2 && hv >= 1) {
        idx.push(2 * (hv + 1) + 1);
    }
    if (upper && hod >= 3) || (!upper && hu >= 1 && hv >= 2) {
        idx.push((hv + 1) + 2);
    }
    if hv >= 3 {
        idx.push(3);
    }
    idx.into_iter().max().unwrap_or(0)
}

fn check_orders(hu: usize, hv: usize, upper: bool) -> Result<(), FfiError> {
    if hu > 6 || hv > 6 {
        return Err(bad("eval: order above 6"));
    }
    // The shim checks zero-offset delegation to SmPlane against actual current
    // distances (including nesting and preceding SetDist/Transform operations).
    let hod = if upper { hu.max(hv) } else { hu + hv };
    // (more than 3 fails before anything is written)
    if hod <= 3 && max_write_index(hu, hv, upper) >= (hu + 1) * (hv + 1) {
        return Err(bad("eval: upstream writes past the derivative array"));
    }
    Ok(())
}

/// A parameter within [`MAX_REACH`] domain lengths of `dom` (any finite one when the
/// domain is not known in advance, bounded by a fixed reach).
fn check_uv(uv: [f64; 2], dom: Option<[f64; 4]>, what: &str) -> Result<(), FfiError> {
    finite(&uv, what)?;
    let [u, v] = uv;
    let ok = match dom {
        Some([u0, v0, u1, v1]) => {
            let (du, dv) = ((u1 - u0) * MAX_REACH, (v1 - v0) * MAX_REACH);
            u >= u0 - du && u <= u1 + du && v >= v0 - dv && v <= v1 + dv
        }
        None => u.abs() <= 1.0e4 && v.abs() <= 1.0e4,
    };
    if !ok {
        return Err(bad(format!("{what}: parameter too far outside the domain")));
    }
    Ok(())
}

fn contained(inner: [f64; 4], outer: Option<[f64; 4]>) -> bool {
    match outer {
        Some([u0, v0, u1, v1]) => {
            let tol = 1e-12;
            inner[0] >= u0 - tol
                && inner[1] >= v0 - tol
                && inner[2] <= u1 + tol
                && inner[3] <= v1 + tol
        }
        None => false,
    }
}

/// The base a box asks for its boxes, through nested offsets.
fn box_base(off: &Off) -> &Srf {
    match &off.base {
        Srf::Offset(inner) => box_base(inner),
        s => s,
    }
}

fn push_steps(o: &mut Vec<f64>, off: &Off, info: BaseInfo, steps: &[Step]) -> Result<(), FfiError> {
    if steps.len() > MAX_OPS {
        return Err(bad("too many ops"));
    }
    // the base domain as the script goes (unknown once an edit may have changed it)
    let mut dom = info.dom;
    let unboxable =
        info.unboxable || matches!(box_base(off), Srf::Bounded { .. } | Srf::Step { .. });
    // an extension may exist once an evaluation ran with extension allowed
    let mut allow = off.allow;
    let mut ext_possible = false;
    o.push(steps.len() as f64);
    for s in steps {
        o.push(s.code());
        match s {
            Step::State | Step::MaxAnal | Step::Copy => {}
            Step::Eval {
                uv,
                hu,
                hv,
                ul,
                vl,
                upper,
                nzt,
            } => {
                check_uv(*uv, dom, "eval")?;
                check_orders(*hu, *hv, *upper)?;
                ext_possible |= allow;
                o.extend_from_slice(uv);
                o.extend([*hu as f64, *hv as f64, b(*ul), b(*vl), b(*upper), b(*nzt)]);
            }
            Step::EvalSimple {
                uv,
                hu,
                hv,
                ul,
                vl,
                upper,
            } => {
                check_uv(*uv, dom, "eval simple")?;
                check_orders(*hu, *hv, *upper)?;
                ext_possible |= allow;
                o.extend_from_slice(uv);
                o.extend([*hu as f64, *hv as f64, b(*ul), b(*vl), b(*upper)]);
            }
            Step::EvalPt { uv } => {
                check_uv(*uv, dom, "eval point")?;
                ext_possible |= allow;
                o.extend_from_slice(uv);
            }
            Step::Normal { uv, ul, vl } => {
                check_uv(*uv, dom, "normal")?;
                ext_possible |= allow;
                o.extend_from_slice(uv);
                o.extend([b(*ul), b(*vl)]);
            }
            Step::Knots { param, ivl } => {
                if *param > 1 {
                    return Err(bad("knots: param"));
                }
                let iv = ivl.unwrap_or([0.0, 0.0]);
                finite(&iv, "knots interval")?;
                o.extend([f64::from(*param), b(ivl.is_some()), iv[0], iv[1]]);
            }
            Step::BBox {
                dom: d,
                want,
                expand,
            } => {
                domain(*d, "box domain")?;
                if unboxable {
                    return Err(bad("box: curve-bounded or STEP base (UB-84)"));
                }
                // Unknown/edited domains are checked against the actual base by the
                // isolated shim. Unknown is not evidence of containment (UB-101).
                if dom.is_some() && !contained(*d, dom) {
                    return Err(bad("box: domain may lie outside the base domain"));
                }
                o.extend_from_slice(d);
                o.extend([b(want[0]), b(want[1]), b(want[2]), b(*expand)]);
            }
            Step::StepInv {
                dom: d,
                p,
                tol,
                guess,
            } => {
                domain(*d, "inversion domain")?;
                finite(p, "inversion point")?;
                if !(tol.is_finite() && *tol > 0.0 && *tol <= 1.0) {
                    return Err(bad("inversion tolerance"));
                }
                let g = guess.unwrap_or([0.0, 0.0]);
                check_uv(g, dom, "inversion guess")?;
                ext_possible |= allow;
                o.extend_from_slice(d);
                o.extend_from_slice(p);
                o.extend([*tol, b(guess.is_some()), g[0], g[1]]);
            }
            Step::Reverse { param } => {
                if *param > 1 {
                    return Err(bad("reverse: param"));
                }
                dom = None;
                o.push(f64::from(*param));
            }
            Step::Swap => {
                dom = None;
            }
            Step::Transform {
                origin,
                x,
                y,
                scale,
            } => {
                finite(origin, "transform origin")?;
                nonzero(*x, "transform x")?;
                nonzero(*y, "transform y")?;
                let c = [
                    x[1] * y[2] - x[2] * y[1],
                    x[2] * y[0] - x[0] * y[2],
                    x[0] * y[1] - x[1] * y[0],
                ];
                if c.iter().map(|v| v * v).sum::<f64>() < 1e-6 {
                    return Err(bad("transform: axes parallel"));
                }
                let sc = scale.unwrap_or([1.0, 1.0, 1.0]);
                finite(&sc, "transform scale")?;
                if sc.iter().any(|&v| v.abs() < 1e-3) {
                    return Err(bad("transform: scale too small"));
                }
                o.extend_from_slice(origin);
                o.extend_from_slice(x);
                o.extend_from_slice(y);
                o.push(b(scale.is_some()));
                o.extend_from_slice(&sc);
            }
            Step::Trim { dom: d } => {
                domain(*d, "trim domain")?;
                if d[2] - d[0] < 1e-6 || d[3] - d[1] < 1e-6 {
                    return Err(bad("trim: domain too small"));
                }
                check_uv([d[0], d[1]], dom, "trim min")?;
                check_uv([d[2], d[3]], dom, "trim max")?;
                ext_possible |= allow;
                dom = None;
                o.extend_from_slice(d);
            }
            Step::Split { t, param } => {
                finite(&[*t], "split")?;
                if *param > 1 {
                    return Err(bad("split: param"));
                }
                // (the base must be split strictly inside its current domain; a plane or a
                // surface of revolution, or a base an earlier op edited, has no domain
                // known here)
                let Some(d) = dom else {
                    return Err(bad("split: base domain unknown"));
                };
                let (a, z) = if *param == 0 {
                    (d[0], d[2])
                } else {
                    (d[1], d[3])
                };
                if *t <= a || *t >= z {
                    return Err(bad("split: parameter not inside the domain"));
                }
                o.extend([*t, f64::from(*param)]);
            }
            Step::CreateOffset { d } | Step::SetDist { d } => {
                finite(&[*d], "offset distance")?;
                o.push(*d);
            }
            Step::Extend { dist, cont } => {
                if !(dist.is_finite() && *dist > 0.0 && *dist <= 1.0e3) || *cont > 13 {
                    return Err(bad("extend: distance or continuity"));
                }
                o.extend([*dist, f64::from(*cont)]);
            }
            Step::ExtendDir { param, dist, cont } => {
                if *param > 8 || !(dist.is_finite() && *dist > 0.0 && *dist <= 1.0e3) || *cont > 13
                {
                    return Err(bad("extend: direction, distance or continuity"));
                }
                o.extend([f64::from(*param), *dist, f64::from(*cont)]);
            }
            Step::Eq { other } => {
                push_off(o, other, 1)?;
            }
            Step::SetAllow { b: flag } => {
                allow = *flag;
                o.push(b(*flag));
            }
            Step::IsPlanar { tol } => {
                if unboxable || ext_possible {
                    return Err(bad("planarity: box not computable without UB"));
                }
                if !(tol.is_finite() && *tol >= 0.0) {
                    return Err(bad("planarity tolerance"));
                }
                o.push(*tol);
            }
            Step::IsKindOf { t } => {
                if t.unsigned_abs() > 1 << 20 {
                    return Err(bad("type code"));
                }
                o.push(*t as f64);
            }
        }
    }
    Ok(())
}

fn call(op: i32, input: &[f64]) -> Result<Vec<f64>, FfiError> {
    crate::offset_reference::call(op, input)
}

/// Static validation of `off` and `steps` (without calling upstream).
///
/// State-dependent guards run inside [`script`]'s isolated reference worker.
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for inputs upstream cannot handle (module docs).
pub fn check(off: &Off, steps: &[Step]) -> Result<(), FfiError> {
    let mut input = Vec::new();
    let info = push_off(&mut input, off, 0)?;
    push_steps(&mut input, off, info, steps)
}

/// Builds `off` and runs `steps` on it; returns the shim's record (the construction
/// status, then each op's status and outputs).
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for inputs upstream cannot handle (module docs);
/// [`FfiError::Status`] for a shim failure, [`REFERENCE_TIMEOUT`], or
/// [`REFERENCE_PROCESS_ERROR`]. Neither deadline nor crash is a successful comparison.
pub fn script(off: &Off, steps: &[Step]) -> Result<Vec<f64>, FfiError> {
    let mut input = Vec::new();
    let info = push_off(&mut input, off, 0)?;
    push_steps(&mut input, off, info, steps)?;
    call(OP_SCRIPT, &input)
}

fn db_header(binary: bool, version: u32) -> Result<[f64; 2], FfiError> {
    if !(30..=1000).contains(&version) {
        return Err(bad("db version"));
    }
    Ok([b(binary), f64::from(version)])
}

/// The bytes upstream's `WriteToDB` writes for `off` (ASCII with 17 digits, or binary,
/// at `version`); `None` when the write fails.
///
/// # Errors
///
/// As [`script`].
pub fn write_db_bytes(off: &Off, binary: bool, version: u32) -> Result<Option<Vec<u8>>, FfiError> {
    let mut input = Vec::new();
    push_off(&mut input, off, 0)?;
    input.extend(db_header(binary, version)?);
    let rec = call(OP_WRITE_DB, &input)?;
    match rec.as_slice() {
        [build, w, bytes @ ..] if *build == 0.0 && *w == 0.0 => bytes
            .iter()
            .map(|&x| {
                if (0.0..=255.0).contains(&x) && x.fract() == 0.0 {
                    Ok(Some(x as u8))
                } else {
                    Err(bad("write db: byte out of range"))
                }
            })
            .collect::<Result<Option<Vec<u8>>, FfiError>>(),
        _ => Ok(None),
    }
}

/// Reads `bytes` with `SmOffsetSurface::ReadFromDB`; returns the status, then the
/// state of [`Step::State`] on success.
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] unless `bytes` are exactly what [`write_db_bytes`]
/// writes for `off` with the same file type and version.
pub fn read_db(off: &Off, binary: bool, version: u32, bytes: &[u8]) -> Result<Vec<f64>, FfiError> {
    if bytes.len() > MAX_DB_BYTES {
        return Err(bad("read db: record too large"));
    }
    if write_db_bytes(off, binary, version)?.as_deref() != Some(bytes) {
        return Err(bad("read db: bytes are not the offset's WriteToDB record"));
    }
    let mut input = Vec::with_capacity(bytes.len() + 3);
    input.extend(db_header(binary, version)?);
    input.push(bytes.len() as f64);
    input.extend(bytes.iter().map(|&x| f64::from(x)));
    call(OP_READ_DB, &input)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plane(d: f64) -> Off {
        Off {
            d,
            owns: true,
            sing: 0,
            allow: false,
            base: Srf::Plane {
                o: [0.0; 3],
                n: [0.0, 0.0, 1.0],
            },
        }
    }

    fn upper() -> Step {
        Step::EvalSimple {
            uv: [0.5, 0.5],
            hu: 2,
            hv: 0,
            ul: true,
            vl: true,
            upper: true,
        }
    }

    #[test]
    fn offset_plane_upper_orders_depend_on_actual_distance() {
        assert!(script(&plane(1.0), &[upper()]).is_ok());
        assert!(script(&plane(1e-20), &[upper()]).is_ok());
        assert!(matches!(
            script(&plane(0.0), &[upper()]),
            Err(FfiError::InvalidArgument(_))
        ));
        assert!(matches!(
            script(&plane(1.0), &[Step::SetDist { d: 0.0 }, upper()]),
            Err(FfiError::InvalidArgument(_))
        ));
        assert!(script(&plane(0.0), &[Step::SetDist { d: 1.0 }, upper()]).is_ok());
        let nested = Off {
            base: Srf::Offset(Box::new(plane(-1.0))),
            ..plane(1.0)
        };
        assert!(matches!(
            script(&nested, &[upper()]),
            Err(FfiError::InvalidArgument(_))
        ));
    }

    #[test]
    fn unknown_revolution_domain_cannot_bypass_polar_guard() {
        let off = Off {
            d: 0.25,
            owns: true,
            sing: 0,
            allow: true,
            base: Srf::Revolution(rev::Surf {
                gen: rev::Gen::Line {
                    a: [1.0, 0.0, 0.0],
                    b: [1.0, 0.0, 1.0],
                },
                build: rev::Build::Canonical {
                    point: [0.0; 3],
                    axis: [0.0, 0.0, 1.0],
                },
                oob: true,
                swap: false,
            }),
        };
        let steps = [
            Step::EvalPt { uv: [0.5, 2.0] },
            Step::BBox {
                dom: [0.0, 0.0, 1.0, 2.0],
                want: [false, false, true],
                expand: true,
            },
        ];
        assert!(check(&off, &steps).is_ok()); // actual domain is a runtime fact
        assert!(matches!(
            script(&off, &steps),
            Err(FfiError::InvalidArgument(_))
        ));
        assert!(script(
            &off,
            &[Step::BBox {
                dom: [0.0, 0.0, 1.0, 1.0],
                want: [false, false, true],
                expand: true
            }]
        )
        .is_ok());
    }

    #[test]
    fn tiny_g2_extension_is_rejected_before_a_nonprogressing_search() {
        assert!(matches!(
            script(
                &plane(1.0),
                &[Step::Extend {
                    dist: 1e-300,
                    cont: 5
                }]
            ),
            Err(FfiError::InvalidArgument(_))
        ));
        // The rejected child is reaped and does not poison the next call.
        assert!(script(&plane(1.0), &[upper()]).is_ok());
    }

    #[test]
    fn gen5909_is_bounded_and_keeps_exact_input_bits() {
        let input: Vec<f64> = include_str!("../../tests/data/offset-gen5909.txt")
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();
        assert_eq!(input.len(), 60);
        assert_eq!(
            input[22].to_bits(),
            (-0.916_536_277_098_507_5_f64).to_bits()
        );
        assert!(matches!(
            call(OP_SCRIPT, &input),
            Err(FfiError::InvalidArgument(_))
        ));
    }
}
