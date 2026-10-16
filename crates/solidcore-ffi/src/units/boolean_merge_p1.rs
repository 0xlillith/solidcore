//! Safe wrappers for `boolean.merge_p1`; shim: `harness/cpp/shim/boolean_merge_p1.cpp`.
//!
//! Operations 1..7 use the validated B-rep block of `topology_vertexuse`. Selectors
//! are nonnegative integers below 32 and resolved modulo actual entity counts by
//! the shim. Operations requiring absent entities reject the fixture before call.
use super::topology_vertexuse::{self as fixture, Op};
use crate::{call_with_buffer, FfiError};
use std::ffi::CString;
mod sys {
    extern "C" {
        pub fn sc_boolean_merge_p1_eval(
            op: i32,
            a: i32,
            b: i32,
            c: i32,
            path: *const std::ffi::c_char,
            input: *const f64,
            n: i64,
            out: *mut f64,
            cap: i64,
            len: *mut i64,
        ) -> i32;
    }
}
/// Call a merge helper. Result: status, record length, record, then a raw graph
/// export in the `topology_vertexuse` format. Op 0 returns only control flags.
pub fn eval(op: u8, args: [u8; 3], path: &str, block: &[f64]) -> Result<Vec<f64>, FfiError> {
    if op > 7 || args.iter().any(|&a| a > 31) {
        return Err(FfiError::InvalidArgument("merge op/selectors".into()));
    }
    if op != 0 {
        let export = fixture::eval(Op::Export, path, block)?;
        let counts = fixture::counts(&export)
            .ok_or_else(|| FfiError::InvalidArgument("fixture rejected".into()))?;
        let valid = match op {
            2 => counts[3] > 0,
            3 => counts[9] > 0 && counts[0] > 0,
            4 => counts[0] > 0,
            5 => counts[2] > 0,
            6 => counts[7] > 0,
            _ => true,
        };
        if !valid {
            return Err(FfiError::InvalidArgument(
                "fixture lacks selected entity".into(),
            ));
        }
    } else if !block.is_empty() {
        return Err(FfiError::InvalidArgument("control takes no fixture".into()));
    }
    let path = CString::new(path).map_err(|_| FfiError::InvalidArgument("NUL in path".into()))?;
    let [a, b, c] = args;
    call_with_buffer(4096, |out, len| {
        // SAFETY: the fixture parser validates every input read and construction;
        // selectors are bounded and checked against actual graph counts in C++.
        // All buffers and the NUL-terminated path live for this call; no pointer escapes.
        unsafe {
            sys::sc_boolean_merge_p1_eval(
                i32::from(op),
                i32::from(a),
                i32::from(b),
                i32::from(c),
                path.as_ptr(),
                block.as_ptr(),
                block.len() as i64,
                out.as_mut_ptr(),
                out.len() as i64,
                &raw mut *len,
            )
        }
    })
}
