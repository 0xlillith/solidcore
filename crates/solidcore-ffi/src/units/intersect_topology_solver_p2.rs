//! FFI for the topology p2 traversal using the validated p1 fixture protocol.
//!
//! Upstream: `source/SMLib/src/SmTopologySolver.cpp:3024` (`LocalSolve`).
//! The C++ entry runs real `SolveTrees`/`SolveTwoTrees`, which invoke `LocalSolve`.
use super::intersect_topology_solver_p1::{validate, Op};
use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_intersect_topology_solver_p2_local(
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_intersect_topology_solver_p2_point_storage(counts: *mut i64) -> i32;
        pub fn sc_intersect_topology_solver_p2_trees(
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
    }
}

/// Runs topology-tree traversal on a p1-protocol tree record.
///
/// # Errors
/// Invalid counts, flags, objects or geometry are rejected by the shared validator;
/// C++ setup errors remain harness errors.
pub fn eval_trees(input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(Op::Trees, input)?;
    call_with_buffer(1024, |out, len| {
        // SAFETY: the complete record is validated by the same parser as the
        // delegated C++ reference; output capacity is enforced by the shim.
        unsafe {
            sys::sc_intersect_topology_solver_p2_trees(
                input.as_ptr(),
                input.len() as i64,
                out.as_mut_ptr(),
                out.len() as i64,
                &raw mut *len,
            )
        }
    })
}

/// C++ object-list counts for external and managed point nodes, respectively.
/// Used to verify the fixture storage convention independently of the port.
///
/// # Errors
/// C++ setup failures.
pub fn point_storage_counts() -> Result<[i64; 2], FfiError> {
    let mut counts = [0; 2];
    // SAFETY: the fixed-size shim writes exactly two i64 entries.
    let status = unsafe { sys::sc_intersect_topology_solver_p2_point_storage(counts.as_mut_ptr()) };
    crate::check(status)?;
    Ok(counts)
}

/// Runs a validated local-solver fixture (including p2's `CurveVertex` class).
///
/// # Errors
/// Invalid records or C++ setup failures.
pub fn eval_local(input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(Op::Local, input)?;
    call_with_buffer(1024, |out, len| {
        // SAFETY: the complete p1-protocol local record is validated above;
        // the delegated reference enforces output capacity.
        unsafe {
            sys::sc_intersect_topology_solver_p2_local(
                input.as_ptr(),
                input.len() as i64,
                out.as_mut_ptr(),
                out.len() as i64,
                &raw mut *len,
            )
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cpp_distinguishes_external_and_managed_point_nodes() {
        assert_eq!(point_storage_counts().unwrap(), [0, 1]);
    }
}
