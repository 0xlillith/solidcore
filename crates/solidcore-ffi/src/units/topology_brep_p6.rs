//! Checked native `SmBrep::{CopyEdges, MakeFaceTopology, TestClosureSplitRegion}`.
//! Uses the topology.vertexuse fixture block and the traverser's raw input export.

use super::topology_topology_traverser::{self as graph, Fixture};
use super::topology_vertexuse::counts;
use crate::{call_with_buffer, FfiError};

/// One request. `op`: closure=0, copy=1, edgeless face=2. Flag bits are
/// infinite-region testing / calculate trims (1), optional copy list (2), null
/// surface (4), opposite orientation (8), and radial verification (16).
#[derive(Debug, Clone)]
pub struct Request {
    /// Operation, 0..=2.
    pub op: i32,
    /// Existing face-use index for closure testing.
    pub seed: i64,
    /// Edge indices in copy order (duplicates are valid).
    pub indices: Vec<i64>,
    /// Repeat with the same topology/context once or twice.
    pub repeat: i32,
    /// Real mark to hold locked across the operation, or zero.
    pub lock_mask: i32,
    /// Operation flag bits, 0..=31.
    pub flags: i32,
}

mod sys {
    extern "C" {
        pub fn sc_topology_brep_p6_eval(
            op: i32, path: *const std::ffi::c_char, input: *const f64, n: i64,
            seed: i64, indices: *const i64, count: i64, repeat: i32,
            lock_mask: i32, flags: i32, out: *mut f64, cap: i64, len: *mut i64,
        ) -> i32;
    }
}

/// Validate construction and all operation-dependent selectors before native code.
pub fn export(path: &str, input: &[f64]) -> Result<Fixture, FfiError> {
    graph::export(path, input)
}

/// Execute a checked native request. `[-1, status]` means fixture construction
/// failed, `[-2]` means a trimmed curve constructor returned null; both happen
/// before the operation. Other statuses/exceptions are operation/harness errors.
pub fn eval(path: &str, input: &[f64], request: &Request) -> Result<Vec<f64>, FfiError> {
    let invalid = || FfiError::InvalidArgument("brep_p6 request".into());
    if !(0..=2).contains(&request.op)
        || !(1..=2).contains(&request.repeat)
        || ![0, 1, 2, 4, 8, 16].contains(&request.lock_mask)
        || !(0..=31).contains(&request.flags)
        || request.indices.len() > 64
        || (request.op != 1 && !request.indices.is_empty())
    { return Err(invalid()); }
    let fixture = export(path, input)?;
    if let Some(c) = counts(&fixture.graph) {
        if request.op == 0 && (request.seed < 0 || request.seed as usize >= c[2]) {
            return Err(invalid());
        }
        if request.indices.iter().any(|i| *i < 0 || *i as usize >= c[7]) {
            return Err(invalid());
        }
    }
    let path = std::ffi::CString::new(path).map_err(|_| invalid())?;
    call_with_buffer(512, |buffer, len| {
        // SAFETY: the shared fixture parser validated every byte of input; all
        // indices, flags and repeat bounds were checked above. Arrays and output
        // buffers remain live for the call; native objects never cross FFI.
        unsafe {
            sys::sc_topology_brep_p6_eval(
                request.op, path.as_ptr(), input.as_ptr(), input.len() as i64,
                request.seed, request.indices.as_ptr(), request.indices.len() as i64,
                request.repeat, request.lock_mask, request.flags,
                buffer.as_mut_ptr(), buffer.len() as i64, &raw mut *len,
            )
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_selectors_flags_and_counts_before_native_execution() {
        let base = Request { op: 1, seed: 0, indices: vec![], repeat: 1, lock_mask: 0, flags: 2 };
        for request in [
            Request { indices: vec![-1], ..base.clone() },
            Request { indices: vec![i64::MAX], ..base.clone() },
            Request { indices: vec![0; 65], ..base.clone() },
            Request { op: 0, seed: -1, ..base.clone() },
            Request { op: 0, seed: i64::MAX, ..base.clone() },
            Request { repeat: 0, ..base.clone() },
            Request { repeat: 3, ..base.clone() },
            Request { lock_mask: 31, ..base.clone() },
            Request { flags: 32, ..base },
        ] {
            assert!(matches!(eval("", &[3.0, 0.0, 0.0], &request), Err(FfiError::InvalidArgument(_))));
        }
    }
}
