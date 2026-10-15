//! Guarded native SmLoop/SmLoopuse fixture queries.
use crate::{call_with_buffer,FfiError};
use super::topology_edgeuse_p1 as graph;

/// Compare native loop getters, optionally setting the actual loop-use orientation.
/// # Errors
/// Invalid fixture records, native setup failure or unsupported orientation.
pub fn query(path:&str,block:&[f64],orientation:Option<u8>)->Result<Vec<f64>,FfiError> {
    if orientation.is_some_and(|v| v>2) {return Err(FfiError::InvalidArgument("orientation outside 0..=2".into()));}
    let snapshot=graph::eval(graph::Op::Export,path,block)?;
    if graph::counts(&snapshot).is_none() {return Err(FfiError::InvalidArgument("native loop fixture setup failed".into()));}
    let path=std::ffi::CString::new(path).map_err(|_|FfiError::InvalidArgument("path contains NUL".into()))?;
    call_with_buffer(1024,|out,len|{
        // SAFETY: the complete geometry record was validated by the shared
        // graph wrapper; the orientation is bounded and output buffers stay live.
        unsafe {sys::sc_topology_loopuse_query(path.as_ptr(),block.as_ptr(),block.len() as i64,
            orientation.map_or(-1,i32::from),out.as_mut_ptr(),out.len() as i64,&raw mut *len)}
    })
}
mod sys {
    unsafe extern "C" {
        pub fn sc_topology_loopuse_query(path:*const std::ffi::c_char,input:*const f64,count:i64,
            orientation:i32,output:*mut f64,capacity:i64,output_len:*mut i64)->i32;
    }
}
