//! Safe wrapper for `blend.hcr`; shim: `harness/cpp/shim/blend_hcr.cpp`.
//!
//! Input: the B-rep block of unit `topology.vertexuse` (validated by building it through
//! [`fixture::eval`]), then the HCR block `view_mode V[3] (O[3] Z[3] Xref[3] for mode 1:
//! V is then the origin) approx_tol ang_tol seams smooth n_sil sil[..] n_back back[..]`
//! and the op arguments:
//! - 0 construction state, 6 `GetSilhouetteEdges`, 8 scene only: none;
//! - 1 `IsSilhouetteEdge`, 2 `IsBackFacingEdge`: `edge`;
//! - 3 `IntersectCurveWithSilhouettes`, 4 `ComputeCurveVisibility`: `edge fa fb
//!   start_hint end_hint` (the classified interval is the edge interval evaluated at
//!   `0 <= fa < fb <= 1`; hints 0..=2);
//! - 5 `ComputeEdgeVisibility`: `edge vertex(-1 none) visible premark`;
//! - 7 `IsPointVisible`: `P[3]`.
//!
//! Output: `status record_len record scene_box(6) scene_export` (status `-2`: the
//! upstream path is undefined and was not run). Rejected here: a view vector that cannot
//! be unitized (UB-169: upstream builds the view transform from uninitialised vectors),
//! entity indices out of range, seeded silhouette edges without curve, non-finite
//! numbers, more than 64 seeds.
use super::topology_vertexuse::{self as fixture, Op};
use crate::{call_with_buffer, FfiError};
use std::ffi::CString;

mod sys {
    extern "C" {
        pub fn sc_blend_hcr_eval(
            op: i32,
            path: *const std::ffi::c_char,
            input: *const f64,
            n: i64,
            out: *mut f64,
            cap: i64,
            len: *mut i64,
        ) -> i32;
    }
}

fn bad(msg: &str) -> FfiError {
    FfiError::InvalidArgument(format!("blend.hcr: {msg}"))
}

struct Reader<'a> {
    v: &'a [f64],
    at: usize,
}

impl Reader<'_> {
    fn num(&mut self) -> Result<f64, FfiError> {
        let x = *self.v.get(self.at).ok_or_else(|| bad("short input"))?;
        self.at += 1;
        if !x.is_finite() {
            return Err(bad("non-finite number"));
        }
        Ok(x)
    }
    fn index(&mut self, n: usize, allow_none: bool) -> Result<i64, FfiError> {
        let x = self.num()?;
        if x.fract() != 0.0 || x < if allow_none { -1.0 } else { 0.0 } || x >= n as f64 {
            return Err(bad("index out of range"));
        }
        Ok(x as i64)
    }
    fn flag(&mut self) -> Result<(), FfiError> {
        let x = self.num()?;
        if x != 0.0 && x != 1.0 {
            return Err(bad("flag"));
        }
        Ok(())
    }
}

/// Runs op `op` (0..=8) on the B-rep block `brep` (file `path` for a file block) and the
/// HCR block `rest`.
pub fn eval(op: u8, path: &str, brep: &[f64], rest: &[f64]) -> Result<Vec<f64>, FfiError> {
    if op > 8 {
        return Err(bad("op"));
    }
    let export = fixture::eval(Op::Export, path, brep)?;
    let c = fixture::counts(&export).ok_or_else(|| bad("fixture rejected"))?;
    let (n_e, n_v) = (c[7], c[9]);
    let mut rd = Reader { v: rest, at: 0 };
    let mode = rd.num()?;
    let (x, y, z) = (rd.num()?, rd.num()?, rd.num()?);
    if mode == 0.0 {
        // UB-169: SmVector3d::Unitize fails below SM_EFF_ZERO_SQ
        if x * x + y * y + z * z < 1.0e-24 {
            return Err(bad("UB-169: view vector cannot be unitized"));
        }
    } else if mode == 1.0 {
        for _ in 0..6 {
            rd.num()?;
        }
    } else {
        return Err(bad("view mode"));
    }
    rd.num()?;
    rd.num()?;
    rd.flag()?;
    rd.flag()?;
    for _ in 0..2 {
        let n = rd.num()?;
        if n.fract() != 0.0 || !(0.0..=64.0).contains(&n) {
            return Err(bad("seed count"));
        }
        for _ in 0..n as usize {
            rd.index(n_e, false)?;
        }
    }
    match op {
        1 | 2 => {
            rd.index(n_e, false)?;
        }
        3 | 4 => {
            rd.index(n_e, false)?;
            let (fa, fb) = (rd.num()?, rd.num()?);
            if !(0.0 <= fa && fa < fb && fb <= 1.0) {
                return Err(bad("interval fractions"));
            }
            for _ in 0..2 {
                let h = rd.num()?;
                if h != 0.0 && h != 1.0 && h != 2.0 {
                    return Err(bad("hint"));
                }
            }
        }
        5 => {
            rd.index(n_e, false)?;
            rd.index(n_v, true)?;
            rd.flag()?;
            rd.flag()?;
        }
        7 => {
            for _ in 0..3 {
                rd.num()?;
            }
        }
        _ => {}
    }
    if rd.at != rest.len() {
        return Err(bad("trailing input"));
    }
    let mut input = brep.to_vec();
    input.extend_from_slice(rest);
    let cpath = CString::new(path).map_err(|_| bad("NUL in path"))?;
    call_with_buffer(1 << 16, |out, len| {
        // SAFETY: the B-rep block and every index are validated above against the built
        // fixture (the scene is its copy); the shim rejects seeded edges without curve
        // and skips (status -2) the undefined upstream paths. The buffers and path live
        // for the call and no C++ pointer escapes.
        unsafe {
            sys::sc_blend_hcr_eval(
                i32::from(op),
                cpath.as_ptr(),
                input.as_ptr(),
                input.len() as i64,
                out.as_mut_ptr(),
                out.len() as i64,
                &raw mut *len,
            )
        }
    })
}
