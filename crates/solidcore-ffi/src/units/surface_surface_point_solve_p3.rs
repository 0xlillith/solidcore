//! Native wrappers for the five `surface.surface_point_solve_p3` functions.
//!
//! Common input is p1's validated problem `W S D solver test(3) loose`.
//! Global adds `has_target target request`; local adds `guess(2)`; stepper
//! operations add `uv(2) tight next_uv(2) iteration hits_u hits_v`.
//! Outputs retain upstream statuses, termination reasons, iteration counts,
//! solution metadata and parameter values. The final global/local flag verifies
//! restoration of the surface's out-of-bounds evaluation state.

use super::surface_surface_cache_p1::{bad, Reader};
use super::surface_surface_point_solve_p1::problem;
use crate::{call_with_buffer, FfiError};

/// Which assigned upstream function to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum Op {
    /// `SmSurface::GlobalPointSolve_1`.
    Global,
    /// `SmSurface::LocalPointSolve_1`.
    Local,
    /// `SmSurfPtSolveStepper::FindBoundarySolution`.
    Boundary,
    /// `SmSurfPtSolveStepper::FindNextStep`.
    Next,
    /// `SmSurfPtSolveStepper::CheckCrease` (disabled upstream).
    Crease,
}

mod sys {
    extern "C" {
        pub fn sc_surface_surface_point_solve_p3_eval(
            op: i32,
            input: *const f64,
            n: i64,
            output: *mut f64,
            capacity: i64,
            len: *mut i64,
        ) -> i32;
    }
}

/// Validate the complete flat input before it reaches the C++ reader.
///
/// Surface counts, degrees, knot vectors and domains use the p1 guards. Stepper
/// counters model states reached within a local solve, with at most 35 iterations
/// and 64 boundary hits. Finite local guesses may lie outside the domain: the
/// upstream entry point clamps them. No caller-declared count is allocated here.
pub fn validate(op: Op, input: &[f64]) -> Result<(), FfiError> {
    let mut rd = Reader { input, at: 0 };
    problem(&mut rd)?;
    validate_tail(op, &mut rd)?;
    rd.done()
}

fn validate_tail(op: Op, rd: &mut Reader<'_>) -> Result<(), FfiError> {
    match op {
        Op::Global => {
            rd.flag()?;
            rd.finite()?;
            rd.index(3)?;
        }
        Op::Local => {
            rd.finite()?;
            rd.finite()?;
        }
        Op::Boundary | Op::Next | Op::Crease => {
            for _ in 0..5 {
                rd.finite()?;
            }
            rd.index(35)?;
            for _ in 0..2 {
                let count = rd.finite()?;
                if count.fract() != 0.0 || count.abs() > 64.0 {
                    return Err(bad("invalid stepper boundary counter"));
                }
            }
        }
    }
    Ok(())
}

/// Run a native case after checking its full wire representation.
pub fn eval(op: Op, input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(op, input)?;
    call_with_buffer(64, |output, len| {
        // SAFETY: validation proves the complete input layout and all C++ allocation
        // counts. Output/len belong to the caller, and the shim checks capacity and
        // keeps exceptions and every allocated geometry object inside C++.
        unsafe {
            sys::sc_surface_surface_point_solve_p3_eval(
                op as i32,
                input.as_ptr(),
                input.len() as i64,
                output.as_mut_ptr(),
                output.len() as i64,
                &raw mut *len,
            )
        }
    })
}

/// Run a case on a real surface of revolution using p2's native construction guard.
pub fn eval_rev(
    op: Op,
    wrap: bool,
    surface: &super::surface_surf_of_revolution_p1::Surf,
    rest: &[f64],
) -> Result<Vec<f64>, FfiError> {
    use super::surface_surface_cache_p1::Srf;
    use super::surface_surface_point_solve_p2::{rev_head, rev_info};
    let (status, info) = rev_info(wrap, surface)?;
    let Some((dom, _, _)) = info else {
        return Err(bad(format!(
            "surface of revolution construction failed: {status}"
        )));
    };
    if !dom.iter().all(|v| v.is_finite()) || dom[0] >= dom[2] || dom[1] >= dom[3] {
        return Err(bad("invalid surface of revolution domain"));
    }
    let mut rd = Reader { input: rest, at: 0 };
    rd.domain_in(
        &Srf {
            dom,
            u_distinct: 0,
            v_distinct: 0,
        },
        true,
    )?;
    if !super::surface_surface_point_solve_p1::OPS.contains(&rd.index(10)?) {
        return Err(bad("unsupported solver operation"));
    }
    rd.vec3()?;
    rd.finite()?;
    validate_tail(op, &mut rd)?;
    rd.done()?;
    let mut input = rev_head(wrap, surface)?;
    input.extend_from_slice(rest);
    call_with_buffer(64, |output, len| {
        // SAFETY: the surface constructor and full problem/tail are validated above;
        // buffers and exception containment are the same as in eval.
        unsafe {
            sys::sc_surface_surface_point_solve_p3_eval(
                op as i32,
                input.as_ptr(),
                input.len() as i64,
                output.as_mut_ptr(),
                output.len() as i64,
                &raw mut *len,
            )
        }
    })
}
