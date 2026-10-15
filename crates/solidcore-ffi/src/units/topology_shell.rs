//! Native shell/region queries on genuine SMLib B-rep fixtures.
//!
//! Upstream: `SmShell.cpp`, `SmRegion.cpp` and their headers.
//! Fixture blocks use the validated `topology.edgeuse_p1` protocol. Results are
//! flat integer-valued records; pointers remain inside the native callback.

use super::topology_edgeuse_p1::{self as fixture, Op};
use crate::{call_with_buffer, FfiError};

fn invalid(reason: impl Into<String>) -> FfiError {
    FfiError::InvalidArgument(reason.into())
}

/// Compare-ready shell/region records, optionally setting each region's native
/// void state. `SmBoolean` is an integer; its indeterminate state is preserved.
///
/// # Errors
/// Malformed fixture inputs, failed native construction, or unsafe native
/// mate/owner references are harness errors rather than comparison results.
pub fn query(path: &str, block: &[f64], void_state: Option<i32>) -> Result<Vec<f64>, FfiError> {
    let export = fixture::eval(Op::Export, path, block)?;
    if fixture::counts(&export).is_none() {
        return Err(invalid("native topology fixture construction failed"));
    }
    let path = std::ffi::CString::new(path).map_err(|_| invalid("path contains NUL"))?;
    let count = i64::try_from(block.len()).map_err(|_| invalid("fixture too long"))?;
    call_with_buffer(1024, |buffer, len| {
        // SAFETY: the existing fixture export validated the entire input block.
        // Strings, arrays and length pointers remain live throughout this call;
        // the shim guards native pointers and honours the output capacity.
        unsafe {
            sys::sc_topology_shell_query(
                path.as_ptr(),
                block.as_ptr(),
                count,
                i32::from(void_state.is_some()),
                void_state.unwrap_or(0),
                buffer.as_mut_ptr(),
                buffer.len() as i64,
                &raw mut *len,
            )
        }
    })
}

mod sys {
    unsafe extern "C" {
        pub fn sc_topology_shell_query(
            path: *const std::ffi::c_char,
            input: *const f64,
            count: i64,
            set_void: i32,
            void_state: i32,
            output: *mut f64,
            capacity: i64,
            output_len: *mut i64,
        ) -> i32;
    }
}
