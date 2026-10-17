//! FFI wrappers for unit `blend.fillet_solver_p1`.
//!
//! Shim: `harness/cpp/shim/blend_fillet_solver_p1.cpp` (`sc_blend_fillet_solver_p1_eval`).
//! Upstream: `source/SMLib/src/SmFilletSolver.cpp`.
//!
//! Surface and curve blocks are those of unit `topology.edgeuse_p1` (`1 O[3] N[3]` an
//! `SmPlane`, `0 oob form n m p q r s cpts u v` an `SmBSplineSurface`; a curve
//! `dim n p m cpts knots`). The solver is a minimal `SmFilletSolver` subclass (every
//! virtual is the base version), no edge uses or executive, orientations `1.0`, saved
//! offsets `0.0`. Arrays are `nx x[nx] nf f[nf] has_j rows cols` (a zero-filled
//! `rows x cols` Jacobian when `has_j`).
//!
//! | [`Op`]       | upstream | parameters |
//! |--------------|----------|------------|
//! | `State`      | constructor, getters, `SetSurface`, `SetSurfacesToZeroOffset`, `ReloadSurfaceOffsets`, `GetOffsetsData`, `SetOffsetsData`, `SetExtensionFactor`, `SetReverseTrim`, `SetStatus` | `tol ang tang surfA surfB has0 off0 has1 off1 o0 o1 k`, then `k` x `cmd a b` |
//! | `Extension`  | `SetupOffsetExtension` | `ext surf u0 v0 u1 v1` |
//! | `Sphere`     | `SmFilletSphereSolveENFO::Evaluate` | `c[3] r adjusts arrays` |
//! | `RailUv`     | `SmRailUVCurveIntersectENFO::Evaluate` | `rail curve t0 t1 arrays` |
//! | `SurfValues` | `ComputeSurfaceValues` | `uv0[2] uv1[2]` |
//! | `UvVectors`  | `ComputeUVVectors` | `du0 dv0 du1 dv1 crv_deriv` (each `[3]`), `uvd0[2] uvd1[2]` |
//! | `Split`      | `SplitFilletGeom` with empty interval arrays | `index` |
//! | `Finish`     | `FinishConstruction` | `surf has0 has1 off` |
//!
//! `State` commands: `0` zero offsets, `1` reload offsets, `2` save offsets data, `3`
//! restore it, `4` `SetSurface(a, new SmOffsetSurface(b, base[a]))`, `5`
//! `SetExtensionFactor(a)`, `6` `SetReverseTrim(a != 0)`, `7` `SetStatus(a)`.
//!
//! Output: `State` the 21 constructor values, then per command the state (per rail
//! `present key offset orientation saved`, then `ext reverse status`); `Extension` `0`
//! and the domain; `Sphere`, `RailUv` `status found f[nf]` and the Jacobian; others the
//! status (`Split` plus the geometry count, `UvVectors` plus both UV derivatives).
//!
//! [`eval`] accepts only what upstream handles without undefined behaviour: valid
//! blocks (see `topology.edgeuse_p1`), finite numbers, at most [`MAX_ARRAY`] array
//! entries and [`MAX_COMMANDS`] commands, rails `0`/`1`, statuses of `SmFilStatus`,
//! `Sphere` with `nf >= rows`, `RailUv` with `nx >= 3`, `nf >= 2` and a Jacobian of at
//! least 2 rows, 3 columns and `rows <= nf`; `Split` only for geometry index 0 (the
//! solver has one geometry: `m_vFilletGeoms[i]` past it is out of bounds; non-empty
//! intervals reach `GetRailVertex` on the constructor's rail-less geometry, a `NULL`
//! dereference).

use crate::units::topology_edgeuse_p1::{bad, curve, surface, Reader};
use crate::{call_with_buffer, FfiError};

/// Largest array length accepted by [`eval`].
pub const MAX_ARRAY: usize = 32;
/// Largest number of `State` commands accepted by [`eval`].
pub const MAX_COMMANDS: usize = 32;
/// Largest `SmFilStatus` value.
const MAX_STATUS: usize = 12;

/// The upstream routine run by the shim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// Constructor, getters and the offset bookkeeping script.
    State,
    /// `SetupOffsetExtension`.
    Extension,
    /// `SmFilletSphereSolveENFO::Evaluate`.
    Sphere,
    /// `SmRailUVCurveIntersectENFO::Evaluate`.
    RailUv,
    /// `ComputeSurfaceValues`.
    SurfValues,
    /// `ComputeUVVectors`.
    UvVectors,
    /// `SplitFilletGeom` with empty interval arrays.
    Split,
    /// `FinishConstruction`.
    Finish,
}

impl Op {
    fn code(self) -> i32 {
        match self {
            Op::State => 0,
            Op::Extension => 1,
            Op::Sphere => 2,
            Op::RailUv => 3,
            Op::SurfValues => 4,
            Op::UvVectors => 5,
            Op::Split => 6,
            Op::Finish => 7,
        }
    }
}

mod sys {
    extern "C" {
        pub fn sc_blend_fillet_solver_p1_eval(
            op: i32,
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
    }
}

/// Reads an array block; returns `(nx, nf, jacobian rows/cols)`.
fn arrays(rd: &mut Reader<'_>) -> Result<(usize, usize, Option<(usize, usize)>), FfiError> {
    let nx = rd.index(MAX_ARRAY)?;
    rd.finites(nx)?;
    let nf = rd.index(MAX_ARRAY)?;
    rd.finites(nf)?;
    let has_j = rd.flag()?;
    let rows = rd.index(MAX_ARRAY)?;
    let cols = rd.index(MAX_ARRAY)?;
    Ok((nx, nf, has_j.then_some((rows, cols))))
}

fn validate_state(rd: &mut Reader<'_>) -> Result<(), FfiError> {
    rd.finites(3)?;
    surface(rd)?;
    surface(rd)?;
    for _ in 0..2 {
        rd.flag()?;
        rd.finite()?;
    }
    rd.finites(2)?;
    let k = rd.index(MAX_COMMANDS)?;
    for _ in 0..k {
        let cmd = rd.index(7)?;
        match cmd {
            4 => {
                rd.index(1)?;
                rd.finite()?;
            }
            7 => {
                rd.index(MAX_STATUS)?;
                rd.finite()?;
            }
            _ => {
                rd.finites(2)?;
            }
        }
    }
    Ok(())
}

fn validate(op: Op, input: &[f64]) -> Result<(), FfiError> {
    let mut rd = Reader { input, at: 0 };
    match op {
        Op::State => validate_state(&mut rd)?,
        Op::Extension => {
            rd.finite()?;
            surface(&mut rd)?;
            rd.finites(4)?;
        }
        Op::Sphere => {
            rd.finites(4)?;
            rd.flag()?;
            let (_, nf, j) = arrays(&mut rd)?;
            if let Some((rows, _)) = j {
                if nf < rows {
                    return Err(bad("sphere: rF shorter than the Jacobian"));
                }
            }
        }
        Op::RailUv => {
            rd.index(1)?;
            curve(&mut rd)?;
            rd.finites(2)?;
            let (nx, nf, j) = arrays(&mut rd)?;
            if nx < 3 || nf < 2 {
                return Err(bad("rail-uv: need nx >= 3 and nf >= 2"));
            }
            if let Some((rows, cols)) = j {
                if rows < 2 || cols < 3 || nf < rows {
                    return Err(bad("rail-uv: Jacobian needs rows >= 2, cols >= 3, rows <= nf"));
                }
            }
        }
        Op::SurfValues => {
            rd.finites(4)?;
        }
        Op::UvVectors => {
            rd.finites(15 + 4)?;
        }
        Op::Split => {
            rd.index(0)?;
        }
        Op::Finish => {
            surface(&mut rd)?;
            rd.flag()?;
            rd.flag()?;
            rd.finite()?;
        }
    }
    rd.done()
}

/// Runs `op` on `input` (layouts in the module docs); returns the shim's output record.
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for input the shim must not see; [`FfiError::Status`]
/// for a shim failure.
pub fn eval(op: Op, input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(op, input)?;
    call_with_buffer(4096, |buf, len| {
        // SAFETY: input is valid for input.len() reads and holds exactly what the shim
        // reads for `op` (validated above: well-formed surface and curve blocks, every
        // index, count and flag the shim passes to upstream, the array sizes upstream
        // indexes); buf is valid for buf.len() writes and the shim honours the capacity;
        // len is valid.
        unsafe {
            sys::sc_blend_fillet_solver_p1_eval(
                op.code(),
                input.as_ptr(),
                input.len() as i64,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })
}
