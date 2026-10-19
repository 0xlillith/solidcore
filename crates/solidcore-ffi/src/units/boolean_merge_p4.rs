//! Validated entry points to the real `SmMerge::ManifoldBoolean`.
//! Fixed-size operands avoid native file readers and input-sized allocations.
use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_boolean_merge_p4_eval(
            mode: i32,
            a: *const f64,
            b: *const f64,
            operation: i32,
            flags: *const i32,
            manifold: i32,
            distinct_contexts: i32,
            caller_lock: i32,
            out: *mut f64,
            cap: i64,
            len: *mut i64,
        ) -> i32;
    }
}

/// Operand kind (-1 null, 0 empty, 1 isolated vertex, 2 B-spline line wire,
/// 3 analytic plane sheet), origin xyz and dimension/direction xyz.
pub type Operand = [f64; 7];

/// Native operation parameters. Explicit tolerances in the shim are 1e-9 and
/// 0.3 radians; the constructor's automatic tolerance calculation is not tested.
#[derive(Debug, Clone, Copy)]
pub struct Request {
    /// One of the eleven real `SmBooleanOperationType` values, 0..=10.
    pub operation: i32,
    /// Cookie cutter, imprint, postprocess, classify faces, keep other.
    pub flags: [bool; 5],
    /// Try closing single-vertex gaps (`m_bManifoldBoolean`).
    pub manifold: bool,
    /// Give the other B-rep its own native context.
    pub distinct_contexts: bool,
    /// Supply an initially unbound caller-owned mark lock.
    pub caller_lock: bool,
}

fn validate(operand: &Operand) -> Result<(), FfiError> {
    if !matches!(operand[0], -1.0 | 0.0 | 1.0 | 2.0 | 3.0) || operand.iter().any(|v| !v.is_finite())
    {
        return Err(FfiError::InvalidArgument("ManifoldBoolean operand".into()));
    }
    Ok(())
}

fn call(mode: i32, a: &Operand, b: &Operand, request: Request) -> Result<Vec<f64>, FfiError> {
    validate(a)?;
    validate(b)?;
    if !(0..=10).contains(&request.operation) {
        return Err(FfiError::InvalidArgument(
            "ManifoldBoolean operation".into(),
        ));
    }
    let flags = request.flags.map(i32::from);
    call_with_buffer(32768, |out, len| {
        // SAFETY: operand arrays have exactly seven validated finite elements;
        // flags and output buffers are live, and the shim owns every native object.
        unsafe {
            sys::sc_boolean_merge_p4_eval(
                mode,
                a.as_ptr(),
                b.as_ptr(),
                request.operation,
                flags.as_ptr(),
                i32::from(request.manifold),
                i32::from(request.distinct_contexts),
                i32::from(request.caller_lock),
                out.as_mut_ptr(),
                out.len() as i64,
                &raw mut *len,
            )
        }
    })
}

/// Execute the native Boolean and export status, state, surviving operands and
/// every extracted component graph (each component has a record-length prefix).
/// A sole -1 means native fixture construction failed before the operation ran.
pub fn eval(a: &Operand, b: &Operand, request: Request) -> Result<Vec<f64>, FfiError> {
    call(0, a, b, request)
}

/// Export one input graph before any Boolean work; no native result enters Rust.
pub fn export(a: &Operand) -> Result<Vec<f64>, FfiError> {
    call(
        -1,
        a,
        a,
        Request {
            operation: 0,
            flags: [false; 5],
            manifold: true,
            distinct_contexts: false,
            caller_lock: false,
        },
    )
}
