//! Native cone/plane and cone/cone intersection adapters.
//!
//! Cone records: origin, x axis, y axis, bottom/top radii, height, start/end
//! angles, swap UV, inside out, cylinder class (17 doubles). Plane records:
//! origin, x axis, y axis, STEP domain (13 doubles). Domain records contain
//! two rectangles as fractions of each natural NURB domain. The result stores
//! status, needs-general, full curve geometry and NURBS, UV nullness, types and
//! deviations. No native object or allocation crosses the boundary.

use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_surface_cone_p3_intersect(
            op: i32,
            cone: *const f64,
            other: *const f64,
            domains: *const f64,
            tolerance: f64,
            has_tolerance: i32,
            same: i32,
            out: *mut f64,
            cap: i64,
            len: *mut i64,
        ) -> i32;
    }
}

fn cone_valid(v: &[f64]) -> bool {
    v.len() == 17
        && v.iter().all(|x| x.is_finite())
        && v.get(14..)
            .unwrap_or(&[])
            .iter()
            .all(|x| *x == 0.0 || *x == 1.0)
}

/// Calls the native analytic intersection (`op` 0 cone/cone, 1 cone/plane,
/// 2 global cone/plane, 3 global cone/cone). `same` uses the exact same native cone pointer.
///
/// # Errors
/// Invalid records, malformed domains, or native adapter failures.
pub fn intersect(
    op: i32,
    cone: &[f64],
    other: &[f64],
    domains: &[f64; 8],
    tolerance: Option<f64>,
    same: bool,
) -> Result<Vec<f64>, FfiError> {
    let other_valid = if op == 0 || op == 3 {
        cone_valid(other)
    } else {
        other.len() == 13 && other.iter().all(|x| x.is_finite())
    };
    if !(0..=3).contains(&op)
        || !cone_valid(cone)
        || !other_valid
        || !domains.iter().all(|x| x.is_finite())
        || tolerance.is_some_and(|x| !x.is_finite())
        || (same && op != 0 && op != 3)
    {
        return Err(FfiError::InvalidArgument(
            "cone intersection records".into(),
        ));
    }
    // Coincident native cylinders dereference the optional tolerance unconditionally.
    // The shim guards that exact path after computing its geometry, rather than
    // rejecting all calls without a tolerance here.
    call_with_buffer(4096, |out, len| {
        // SAFETY: fixed record lengths and values checked above; output and length
        // are writable, and native-owned curves are destroyed before returning.
        unsafe {
            sys::sc_surface_cone_p3_intersect(
                op,
                cone.as_ptr(),
                other.as_ptr(),
                domains.as_ptr(),
                tolerance.unwrap_or(0.0),
                i32::from(tolerance.is_some()),
                i32::from(same),
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

    fn cone() -> Vec<f64> {
        vec![
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 3.0, 0.0, 360.0, 0.0, 0.0, 0.0,
        ]
    }
    const DOMAINS: [f64; 8] = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0];

    #[test]
    fn rejects_malformed_record_before_entering_cpp() {
        assert!(matches!(
            intersect(0, &[], &cone(), &DOMAINS, None, false),
            Err(FfiError::InvalidArgument(_))
        ));
    }

    #[test]
    fn guards_only_default_tolerance_path_that_dereferences_null() {
        let a = cone();
        assert_eq!(
            intersect(0, &a, &a, &DOMAINS, None, false),
            Err(FfiError::Status(-1002))
        );
        // Identity returns before reading the optional tolerance.
        assert!(intersect(0, &a, &a, &DOMAINS, None, true).is_ok());
        let mut other = cone();
        // Same axis but disjoint heights: the native routine returns before the
        // coincident branch, so omitting the tolerance remains valid.
        if let Some(z) = other.get_mut(2) {
            *z = 10.0;
        }
        assert!(intersect(0, &a, &other, &DOMAINS, None, false).is_ok());
        // Nested radii likewise return before the null dereference.
        let mut nested = cone();
        for i in [9, 10] {
            if let Some(r) = nested.get_mut(i) {
                *r = 2.0;
            }
        }
        assert!(intersect(0, &a, &nested, &DOMAINS, None, false).is_ok());
    }
}
