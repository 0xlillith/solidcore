//! Safe native oracle for curve-classification part 6.
//!
//! Each input interval is `lo hi` followed by three eight-number point records:
//! `kind gap t pre_snap stored_point x y z`. Kinds: 0 unknown, 1 curve, 2 face,
//! 3 edge, 4 vertex, 5 region, 6 nothing. Topology comes from a real unit box.
//! Flags bit 0 selects vertices-only merging; bit 1 enables paired fillet handling.
//! Geometry: 0 line, 1 circle, 2 ellipse, 3 quadratic B-spline, 4 constant B-spline.
//! The classified range is [0,1]; supplied intervals form a partition of that range.
//! Native constructor status/null checks precede topology access. All objects are scoped.

use crate::{call_with_buffer, FfiError};

/// Numbers per interval (two endpoints, three eight-number classifications).
pub const WIDTH: usize = 26;
/// Largest fixture partition (native algorithms may create additional intervals).
pub const MAX_INTERVALS: usize = 32;

mod sys {
    extern "C" {
        pub fn sc_intersect_curve_class_p6_run(
            op: i32,
            kind: i32,
            scale: f64,
            tol: f64,
            flags: i32,
            index: i64,
            a: *const f64,
            na: i64,
            b: *const f64,
            nb: i64,
            input: *const f64,
            out: *mut f64,
            cap: i64,
            len: *mut i64,
        ) -> i32;
    }
}

fn bad() -> FfiError {
    FfiError::InvalidArgument("invalid curve-classification fixture".into())
}

fn interval(v: &[f64]) -> Result<(), FfiError> {
    let v: &[f64; WIDTH] = v.try_into().map_err(|_| bad())?;
    if v.len() != WIDTH
        || v.iter().any(|x| !x.is_finite())
        || v[0] < 0.0
        || v[1] > 1.0
        || v[0] > v[1]
    {
        return Err(bad());
    }
    for pc in v[2..].as_chunks::<8>().0 {
        if pc[0].fract() != 0.0
            || !(0.0..=6.0).contains(&pc[0])
            || pc[1] < 0.0
            || pc[1] > 1e20
            || pc[2].abs() > 1e20
            || pc[3].abs() > 1e20
            || ![0.0, 1.0].contains(&pc[4])
            || pc[5..].iter().any(|x| x.abs() > 1e100)
        {
            return Err(bad());
        }
    }
    Ok(())
}

fn partition(v: &[f64]) -> Result<(), FfiError> {
    if !v.len().is_multiple_of(WIDTH) || v.len() / WIDTH > MAX_INTERVALS {
        return Err(bad());
    }
    let mut previous = 0.0;
    for ivl in v.as_chunks::<WIDTH>().0 {
        interval(ivl)?;
        if ivl[0] != previous || ivl[0] >= ivl[1] {
            return Err(bad());
        }
        previous = ivl[1];
    }
    if !v.is_empty() && previous != 1.0 {
        return Err(bad());
    }
    Ok(())
}

/// Validates the complete native fixture before passing it to C++.
/// Operations: 0 insert, 1 homogenize, 2 fix endpoint, 3 merge interval.
///
/// # Errors
/// Rejects malformed layouts, flags, partitions and invalid operation indices.
#[allow(clippy::too_many_arguments)]
pub fn validate(
    op: i32,
    kind: i32,
    scale: f64,
    tol: f64,
    flags: i32,
    index: usize,
    a: &[f64],
    b: &[f64],
    input: &[f64],
) -> Result<(), FfiError> {
    if !(0..=3).contains(&op)
        || !(0..=4).contains(&kind)
        || !scale.is_finite()
        || !(0.0..=1e100).contains(&scale)
        || scale == 0.0
        || !tol.is_finite()
        || !(0.0..=1e10).contains(&tol)
        || !(0..=3).contains(&flags)
    {
        return Err(bad());
    }
    partition(a)?;
    partition(b)?;
    interval(input)?;
    let count = a.len() / WIDTH;
    if (op != 1 && count == 0)
        || (op == 1 && a.is_empty() != b.is_empty())
        || (op == 2 && index > count)
        || (op == 3 && index >= count)
    {
        return Err(bad());
    }
    if op == 3 && flags & 1 == 0 {
        let record: &[f64; WIDTH] = a
            .get(index * WIDTH..(index + 1) * WIDTH)
            .and_then(|x| x.try_into().ok())
            .ok_or_else(bad)?;
        let classes = [record[2], record[10], record[18]];
        let early = (classes[1] == 0.0 && flags & 1 == 0) || classes == [4.0, 3.0, 4.0];
        // MergeInterval casts classification objects to B-rep topology. A geometry
        // or NOTHING record cannot supply that topology; UNKNOWN is allowed.
        if !early && classes.iter().any(|x| *x == 1.0 || *x == 6.0) {
            return Err(bad());
        }
    }
    // Index is ignored by insert/homogenize, but must fit the native integer type.
    i64::try_from(index).map_err(|_| bad())?;
    Ok(())
}

/// Calls the real native operation and returns status plus complete classification state.
///
/// # Errors
/// See [`validate`]; construction failures are reported as shim errors.
#[allow(clippy::too_many_arguments)]
pub fn run(
    op: i32,
    kind: i32,
    scale: f64,
    tol: f64,
    flags: i32,
    index: usize,
    a: &[f64],
    b: &[f64],
    input: &[f64],
) -> Result<Vec<f64>, FfiError> {
    validate(op, kind, scale, tol, flags, index, a, b, input)?;
    call_with_buffer(32768, |out, len| {
        // SAFETY: full records, finite values, enum/flag ranges, partitions and indices
        // are validated above. The shim checks constructor results and output capacity.
        unsafe {
            sys::sc_intersect_curve_class_p6_run(
                op,
                kind,
                scale,
                tol,
                flags,
                index as i64,
                a.as_ptr(),
                (a.len() / WIDTH) as i64,
                b.as_ptr(),
                (b.len() / WIDTH) as i64,
                input.as_ptr(),
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
    fn ivl() -> Vec<f64> {
        let mut v = vec![0.0, 1.0];
        for _ in 0..3 {
            v.extend([0.0, 0.0, 0.0, 1e20, 0.0, 0.0, 0.0, 0.0]);
        }
        v
    }
    #[test]
    fn rejects_malformed_records_and_invalid_indices() {
        let v = ivl();
        for n in 0..WIDTH {
            assert!(run(0, 0, 1.0, 1e-6, 0, 0, &v[..n], &v, &v).is_err());
        }
        for at in 0..WIDTH {
            for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut bad = v.clone();
                bad[at] = x;
                assert!(run(0, 0, 1.0, 1e-6, 0, 0, &bad, &v, &v).is_err());
            }
        }
        assert!(run(2, 0, 1.0, 1e-6, 0, 2, &v, &v, &v).is_err());
        assert!(run(3, 0, 1.0, 1e-6, 0, 1, &v, &v, &v).is_err());
        assert!(run(0, 0, 1.0, 1e-6, 4, 0, &v, &v, &v).is_err());
        assert!(run(1, 0, 1.0, 1e-6, 0, 0, &[], &v, &v).is_err());
    }
    #[test]
    fn native_degenerate_and_tiny_geometry_constructor_paths() {
        let v = ivl();
        for kind in 0..=4 {
            for scale in [1e-100, 1e-15, 1e-9, 1.0, 1e9, 1e100] {
                // Merge's unknown-interior exit must still invoke the real native API.
                let out = run(3, kind, scale, 1e-6, 0, 0, &v, &v, &v).unwrap();
                assert_eq!(out[0], 0.0);
            }
        }
    }
}

#[cfg(test)]
mod isolated_tests {
    use super::*;
    use std::process::{Command, Stdio};
    use std::time::{Duration, Instant};
    const CHILD: &str = "units::intersect_curve_class_p6::isolated_tests::native_probe_child";
    #[test]
    #[ignore = "only run by the subprocess safety matrix"]
    fn native_probe_child() {
        let s = std::env::var("SC_CC_P6_PROBE").unwrap();
        let args: Vec<f64> = s.split(',').map(|x| x.parse().unwrap()).collect();
        let (op, kind, flags, scale) = (args[0] as i32, args[1] as i32, args[2] as i32, args[3]);
        let mut a = vec![0.0, 1.0];
        for c in &args[4..7] {
            a.extend([*c, 0.0, 0.5, 1e20, 0.0, 0.0, 0.0, 0.0]);
        }
        let input = a.clone();
        // Malformed classification kinds for merge are rejected before native code.
        if validate(op, kind, scale, 1e-6, flags, 0, &a, &a, &input).is_err() {
            return;
        }
        let result = run(op, kind, scale, 1e-6, flags, 0, &a, &a, &input);
        assert!(result.is_ok(), "{s}: {result:?}");
    }
    #[test]
    fn native_operation_geometry_flag_matrix_in_subprocesses() {
        let exe = std::env::current_exe().unwrap();
        let listing = Command::new(&exe).args(["--ignored", "--list"]).output().unwrap();
        assert!(String::from_utf8(listing.stdout).unwrap().contains(CHILD));
        for kind in 0..=4 {
            for flags in 0..=3 {
                for op in 0..=3 {
                    for classes in [
                        [0, 0, 0],
                        [1, 1, 1],
                        [2, 2, 2],
                        [3, 3, 3],
                        [4, 3, 4],
                        [5, 5, 5],
                        [6, 6, 6],
                        [0, 5, 0],
                        [4, 2, 4],
                    ] {
                        let probe = format!(
                            "{op},{kind},{flags},1,{},{},{}",
                            classes[0], classes[1], classes[2]
                        );
                        let mut child = Command::new(&exe)
                            .args(["--ignored", "--exact", CHILD])
                            .env("SC_CC_P6_PROBE", &probe)
                            .stdout(Stdio::null())
                            .stderr(Stdio::null())
                            .spawn()
                            .unwrap();
                        let start = Instant::now();
                        loop {
                            if let Some(status) = child.try_wait().unwrap() {
                                assert!(
                                    status.success(),
                                    "native crash/failure for {probe}: {status}"
                                );
                                break;
                            }
                            if start.elapsed() > Duration::from_secs(5) {
                                child.kill().unwrap();
                                child.wait().unwrap();
                                panic!("native timeout for {probe}");
                            }
                            std::thread::sleep(Duration::from_millis(5));
                        }
                    }
                }
            }
        }
    }
}
