//! Real native PRE-state and ClassifyFaceInterval dispatch probes.
use crate::{call_with_buffer, FfiError};
mod sys {
    extern "C" {
        pub fn sc_intersect_curve_class_p7_dispatch(
            branch: i32,
            export_only: i32,
            out: *mut f64,
            cap: i64,
            len: *mut i64,
        ) -> i32;
    }
}
/// Branch 0: edge endpoint; 1: 3D midpoint; 2: UV midpoint.
/// Native premises require a fresh finite real plane with no retained trim cache.
pub fn probe(branch: u8, export_only: bool) -> Result<Vec<f64>, FfiError> {
    if branch > 2 {
        return Err(FfiError::InvalidArgument("p7 dispatch branch".into()));
    }
    call_with_buffer(4096, |out, len| {
        // SAFETY: scalar branch is validated, and output storage belongs to this call.
        unsafe {
            sys::sc_intersect_curve_class_p7_dispatch(
                i32::from(branch),
                i32::from(export_only),
                out.as_mut_ptr(),
                out.len() as i64,
                &raw mut *len,
            )
        }
    })
}
