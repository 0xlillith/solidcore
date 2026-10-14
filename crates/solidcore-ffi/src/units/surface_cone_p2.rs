//! Safe fixture protocol for actual cone ellipse trim intervals (SmCone.cpp:5282).
use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_surface_cone_p2_trim(
            input: *const f64,
            n: i64,
            out: *mut f64,
            cap: i64,
            len: *mut i64,
        ) -> i32;
    }
}

/// Validate the fixed record before either implementation consumes it.
///
/// # Errors
/// Malformed, nonfinite or unbounded fixture data.
pub fn validate(input: &[f64]) -> Result<(), FfiError> {
    let Some(
        &[bottom, top, height, start, end, swap, inside, tilt, z, rx, ry, a, b, u0, v0, u1, v1, tol],
    ) = input.get(..)
    else {
        return Err(FfiError::InvalidArgument(
            "expected 18 cone trim scalars".into(),
        ));
    };
    if input.iter().any(|v| !v.is_finite() || v.abs() > 1e9)
        || bottom <= 0.0
        || top < 0.0
        || height <= 0.0
        || rx <= 0.0
        || ry <= 0.0
        || start >= end
        || end - start > 360.0
        || start.abs() > 360.0
        || end.abs() > 360.0
        || a > b
        || b - a > 360.0
        || a.abs() > 720.0
        || b.abs() > 720.0
        || tilt.abs() > std::f64::consts::PI
        || tol < 0.0
        || [swap, inside].iter().any(|v| *v != 0.0 && *v != 1.0)
        || [u0, v0, u1, v1].iter().any(|v| !(0.0..=1.0).contains(v))
        || u0 >= u1
        || v0 >= v1
    {
        return Err(FfiError::InvalidArgument(format!(
            "invalid bounded cone trim at z={z}"
        )));
    }
    Ok(())
}

/// Run the native helper after successful geometry and primitive preflight.
///
/// # Errors
/// Invalid fixture data, primitive/setup failure, or reference failure.
pub fn eval_trim(input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(input)?;
    call_with_buffer(32, |out, len| {
        // SAFETY: the complete fixed-size input is validated; the shim bounds
        // output writes and checks the primitive statuses ignored by upstream.
        unsafe {
            sys::sc_surface_cone_p2_trim(
                input.as_ptr(),
                input.len() as i64,
                out.as_mut_ptr(),
                out.len() as i64,
                &raw mut *len,
            )
        }
    })
}
