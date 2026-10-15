//! Native public-query calls. Shim: `harness/cpp/shim/query_api_queries.cpp`.
//!
//! Fixed layout: origin[3], direction[3], scale, domain[4], query[3], parameter,
//! tight, null, `remove_curve`, `inside_out`, `edge_tolerance`, selector, kind.
//! Kinds: 0 line, 1 circle, 2 ellipse, 3 cubic B-spline, 4 plane, 5 quadratic
//! B-spline surface, 6 box, 7 empty B-rep, 8 rational cylinder, 9 rational sphere.
//! Graph exports append `(upward_faceuse, rectangular_trim)` per face to the
//! shared `topology_edgeuse_p1` format.
//! Op 0 exports the native graph; 1..25
//! call `SmApi` queries (see the runner). Op 26 exports the graph with the selected
//! face boundary's last orientation set to Unknown when `inside_out` is true, as
//! op 11 does before querying. Results start with the exact `SmStatus`.
use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_query_api_queries_eval(
            op: i32,
            input: *const f64,
            n: i64,
            out: *mut f64,
            cap: i64,
            len: *mut i64,
        ) -> i32;
    }
}

/// Calls a native query after checking the complete fixed layout and the
/// construction ranges used by the kernel. Op 0 returns a raw topology export.
pub fn eval(op: i32, input: &[f64]) -> Result<Vec<f64>, FfiError> {
    let bad = || FfiError::InvalidArgument("query.api_queries layout or construction range".into());
    if input.len() != 22
        || !(0..=26).contains(&op)
        || input.iter().any(|x| !x.is_finite() || x.abs() > 1e12)
    {
        return Err(bad());
    }
    let a: &[f64; 22] = input.try_into().map_err(|_| bad())?;
    if !(0.0..=9.0).contains(&a[21])
        || a[21].fract() != 0.0
        || a[6] <= 1e-10
        || a[6] > 1e6
        || a[20] < 0.0
        || a[20].fract() != 0.0
        || a[20] > 1000.0
        || [a[15], a[16], a[17], a[18]]
            .iter()
            .any(|&flag| flag != 0.0 && flag != 1.0)
        || a[19] < 0.0
        || a[19] > 1.0
    {
        return Err(bad());
    }
    let kind = a[21] as i32;
    if (kind == 0 || kind == 4) && a[3] * a[3] + a[4] * a[4] + a[5] * a[5] < 1e-20 {
        return Err(bad());
    }
    if kind <= 3 && (a[7] >= a[8] || (kind != 0 && (a[7] < 0.0 || a[8] > 1.0))) {
        return Err(bad());
    }
    if (kind == 4 || kind == 5 || kind >= 8)
        && (a[7] >= a[9]
            || a[8] >= a[10]
            || (kind == 5 && (a[7] < 0.0 || a[8] < 0.0 || a[9] > 1.0 || a[10] > 1.0)))
    {
        return Err(bad());
    }
    if (op == 1 || op == 2) && kind > 3 || op == 3 && !matches!(kind, 4 | 5 | 8 | 9) {
        return Err(bad());
    }
    // SmEdge::ClosestPoint dereferences its curve; all other public calls guard it.
    if op == 21 && a[17] != 0.0 && a[16] == 0.0 {
        return Err(bad());
    }
    call_with_buffer(2048, |out, len| {
        // SAFETY: layout and all construction/selector inputs checked above;
        // all output buffers are supplied by call_with_buffer.
        unsafe {
            sys::sc_query_api_queries_eval(
                op,
                input.as_ptr(),
                22,
                out.as_mut_ptr(),
                out.len() as i64,
                len,
            )
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_malformed_layout_before_ffi() {
        assert!(eval(1, &[]).is_err());
        assert!(eval(1, &[f64::NAN; 22]).is_err());
        assert!(eval(99, &[0.0; 22]).is_err());
    }
}
