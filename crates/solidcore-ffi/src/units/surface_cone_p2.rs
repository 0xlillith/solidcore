//! Guarded native cone trim, recognition, rebuilding and offset fixture protocol.
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
    ) = input.get(..18)
    else {
        return Err(FfiError::InvalidArgument(
            "expected 18 cone trim scalars".into(),
        ));
    };
    if !matches!(input.len(), 18 | 27)
        || input.iter().any(|v| !v.is_finite() || v.abs() > 1e9)
        || bottom < 0.0
        || (bottom == 0.0 && top == 0.0)
        || top < 0.0
        || height <= 0.0
        || rx < 0.0
        || ry < 0.0
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
        || u0 > u1
        || v0 > v1
        // The linear NURB direction must span a nonzero seam line. The angular
        // direction may collapse; swap exchanges which direction is linear.
        || (swap == 0.0 && v0 == v1)
        || (swap == 1.0 && u0 == u1)
    {
        return Err(FfiError::InvalidArgument(format!(
            "invalid bounded cone trim at z={z}"
        )));
    }
    if let Some(&[_, _, _, x0, x1, x2, y0, y1, y2]) = input.get(18..) {
        let x_norm = x0 * x0 + x1 * x1 + x2 * x2;
        let y_norm = y0 * y0 + y1 * y1 + y2 * y2;
        let dot = x0 * y0 + x1 * y1 + x2 * y2;
        if (x_norm - 1.0).abs() > 1e-8 || (y_norm - 1.0).abs() > 1e-8 || dot.abs() > 1e-8 {
            return Err(FfiError::InvalidArgument(
                "pose is not a canonical frame".into(),
            ));
        }
    }
    Ok(())
}

/// Run the native helper after successful geometry and primitive preflight.
///
/// # Errors
/// Invalid fixture data, primitive/setup failure, or reference failure.
pub fn eval_trim(input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(input)?;
    if input.len() != 18 {
        return Err(FfiError::InvalidArgument(
            "trim protocol requires 18 scalars".into(),
        ));
    }
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

mod analytic_sys {
    extern "C" {
        pub fn sc_surface_cone_p2_analytic(
            op: i32,
            kind: i32,
            mutation: i32,
            scale: f64,
            distance: f64,
            tolerance: f64,
            input: *const f64,
            n: i64,
            out: *mut f64,
            cap: i64,
            len: *mut i64,
        ) -> i32;
    }
}

/// Validate analytic fixture geometry and all native options.
///
/// Records contain 18 geometry scalars and optionally nine position/axis scalars.
/// Kind bits select a real cylinder (1), an initially uncached surface (2), and
/// an unbounded surface (4); zero selects a bounded, cached real cone.
/// Mutations select none (0), end twist (1), homogeneous end rescaling (2),
/// noncanonical interior knot (3), flattened control net (4), or nonrational net (5).
/// Mutations require an initially cached NURB.
///
/// # Errors
/// Malformed or out-of-bounds fixture data.
pub fn validate_analytic(
    input: &[f64],
    op: i32,
    kind: i32,
    mutation: i32,
    scale: f64,
    distance: f64,
    tolerance: f64,
) -> Result<(), FfiError> {
    // Shared validation keeps the native and Rust fixture domains identical.
    validate(input)?;
    if !(1..=3).contains(&op)
        || !(0..=7).contains(&kind)
        || !(0..=5).contains(&mutation)
        || !scale.is_finite()
        || !(0.0..=1e6).contains(&scale)
        || !distance.is_finite()
        || distance.abs() > 1e9
        || !tolerance.is_finite()
        || tolerance.abs() > 1e9
        || (kind % 2 == 1 && input.first() != input.get(1))
        || (kind & 2 != 0 && mutation != 0)
        || (kind & 4 != 0
            && input
                .get(7)
                .is_none_or(|angle| angle.abs() >= std::f64::consts::FRAC_PI_2))
    {
        return Err(FfiError::InvalidArgument(
            "invalid cone analytic fixture options".into(),
        ));
    }
    Ok(())
}

/// Run independent native recognition, rebuilding or offsetting.
///
/// # Errors
/// Invalid data, reference construction or state extraction failure.
pub fn eval_analytic(
    input: &[f64],
    op: i32,
    kind: i32,
    mutation: i32,
    scale: f64,
    distance: f64,
    tolerance: f64,
) -> Result<Vec<f64>, FfiError> {
    validate_analytic(input, op, kind, mutation, scale, distance, tolerance)?;
    call_with_buffer(512, |out, len| {
        // SAFETY: the fixed record, operation, kind, mutation and finite options
        // are validated; the shim constructs its own objects and bounds output.
        unsafe {
            analytic_sys::sc_surface_cone_p2_analytic(
                op,
                kind,
                mutation,
                scale,
                distance,
                tolerance,
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
    fn record() -> [f64; 18] {
        [
            1.0, 2.0, 4.0, 20.0, 250.0, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 0.0, 360.0, 0.0, 0.0, 1.0,
            1.0, 1e-7,
        ]
    }
    #[test]
    fn rejects_malformed_records_and_uninitialized_mutation_before_cpp() {
        let input = record();
        assert!(validate(&input[..17]).is_err());
        assert!(validate_analytic(&input, 3, 2, 1, 1.0, 0.0, 1e-7).is_err());
        let mut pose = input.to_vec();
        pose.extend([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0]);
        assert!(validate(&pose).is_err());
        let mut nonfinite = input;
        nonfinite[0] = f64::INFINITY;
        assert!(validate(&nonfinite).is_err());
    }
    #[test]
    fn validates_the_linear_axis_after_uv_swap() {
        let mut input = record();
        input[13] = 0.3;
        input[15] = 0.3;
        assert!(validate(&input).is_ok());
        input[5] = 1.0;
        assert!(matches!(
            validate(&input),
            Err(FfiError::InvalidArgument(_))
        ));
        input[13] = 0.2;
        input[15] = 0.9;
        input[14] = 0.3;
        input[16] = 0.3;
        assert!(validate(&input).is_ok());
        assert!(eval_trim(&input).is_ok());
        input[5] = 0.0;
        assert!(matches!(
            validate(&input),
            Err(FfiError::InvalidArgument(_))
        ));
    }

    #[test]
    fn native_truncation_discards_the_partial_sweep() {
        let mut input = record();
        input[0] = 2.0;
        input[1] = 1.0;
        let output = eval_analytic(&input, 3, 0, 0, 1.0, -1.5, 1e-7).unwrap();
        assert_eq!(output.get(..2), Some([1000.0, 1.0].as_slice()));
        // Source input sweep [20,250] becomes CreateCanonical's [0,360].
        assert_eq!(output.get(16), Some(&0.0));
        assert_eq!(output.get(18), Some(&360.0));
    }
    #[test]
    fn rejects_null_truncation_factory_before_the_unsafe_offset_call() {
        let input = record();
        let baseline = eval_analytic(&input, 3, 0, 0, 1.0, 0.0, 0.0).unwrap();
        assert_eq!(baseline.first(), Some(&1000.0));
        let semi = 1.0_f64.atan2(4.0) * 180.0 / std::f64::consts::PI;
        let cosine = (semi * std::f64::consts::PI / 180.0).cos();
        let mut distance = -2.0 / cosine;
        for _ in 0..4 {
            distance = distance.next_up();
        }
        assert!(matches!(
            eval_analytic(&input, 3, 0, 0, 1.0, distance, 0.0),
            Err(FfiError::Status(-1002))
        ));
    }
}
