//! FFI for the topology p2 traversal using the validated p1 fixture protocol.
//!
//! Upstream: `source/SMLib/src/SmTopologySolver.cpp:3024` (`LocalSolve`).
//! The C++ entry runs real `SolveTrees`/`SolveTwoTrees`, which invoke `LocalSolve`.
use super::intersect_topology_solver_p1::{validate, Op};
use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_intersect_topology_solver_p2_empty_lists(results: *mut i64) -> i32;
        pub fn sc_intersect_topology_solver_p2_edge_vertex_operand() -> i32;
        pub fn sc_intersect_topology_solver_p2_static_fixture(
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_intersect_topology_solver_p2_static(
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_intersect_topology_solver_p2_local(
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_intersect_topology_solver_p2_point_storage(counts: *mut i64) -> i32;
        pub fn sc_intersect_topology_solver_p2_trees(
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
    }
}

/// Runs topology-tree traversal on a p1-protocol tree record.
///
/// # Errors
/// Invalid counts, flags, objects or geometry are rejected by the shared validator;
/// C++ setup errors remain harness errors.
pub fn eval_trees(input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(Op::Trees, input)?;
    call_with_buffer(1024, |out, len| {
        // SAFETY: the complete record is validated by the same parser as the
        // delegated C++ reference; output capacity is enforced by the shim.
        unsafe {
            sys::sc_intersect_topology_solver_p2_trees(
                input.as_ptr(),
                input.len() as i64,
                out.as_mut_ptr(),
                out.len() as i64,
                &raw mut *len,
            )
        }
    })
}

/// C++ object-list counts for external and managed point nodes, respectively.
/// Used to verify the fixture storage convention independently of the port.
///
/// # Errors
/// C++ setup failures.
pub fn point_storage_counts() -> Result<[i64; 2], FfiError> {
    let mut counts = [0; 2];
    // SAFETY: the fixed-size shim writes exactly two i64 entries.
    let status = unsafe { sys::sc_intersect_topology_solver_p2_point_storage(counts.as_mut_ptr()) };
    crate::check(status)?;
    Ok(counts)
}

/// Native solve status and solution count for each order of an empty-list tree.
///
/// # Errors
/// C++ setup failures.
pub fn empty_list_results() -> Result<[i64; 4], FfiError> {
    let mut results = [0; 4];
    // SAFETY: the shim writes exactly four i64 entries to the provided array.
    let status = unsafe { sys::sc_intersect_topology_solver_p2_empty_lists(results.as_mut_ptr()) };
    crate::check(status)?;
    Ok(results)
}

/// Runs a validated local-solver fixture (including p2's `CurveVertex` class).
///
/// # Errors
/// Invalid records or C++ setup failures.
pub fn eval_local(input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(Op::Local, input)?;
    call_with_buffer(1024, |out, len| {
        // SAFETY: the complete p1-protocol local record is validated above;
        // the delegated reference enforces output capacity.
        unsafe {
            sys::sc_intersect_topology_solver_p2_local(
                input.as_ptr(),
                input.len() as i64,
                out.as_mut_ptr(),
                out.len() as i64,
                &raw mut *len,
            )
        }
    })
}

fn validate_static(input: &[f64]) -> Result<(), FfiError> {
    let invalid = || FfiError::InvalidArgument("unsafe topology p2 static fixture".into());
    let Some(&[entry, operation, request, tolerance, best, n1, n2, _, _, _]) = input.get(..10)
    else {
        return Err(invalid());
    };
    if input.iter().any(|v| !v.is_finite())
        || entry.fract() != 0.0
        || !(0.0..=8.0).contains(&entry)
        || operation.fract() != 0.0
        || !(0.0..=18.0).contains(&operation)
        || request.fract() != 0.0
        || !(0.0..=3.0).contains(&request)
        || !(0.0..=1e9).contains(&tolerance)
        || best.abs() > 1e30
        || [n1, n2]
            .iter()
            .any(|n| n.fract() != 0.0 || !(1.0..=64.0).contains(n))
        || input
            .get(7..)
            .is_none_or(|v| v.iter().any(|x| x.abs() > 1e9))
    {
        return Err(invalid());
    }
    let count = n1 as usize + n2 as usize;
    if input.len() != 12 + 3 * count
        || input
            .get(input.len().saturating_sub(2)..)
            .is_none_or(|flags| flags.iter().any(|v| *v != 0.0 && *v != 1.0))
    {
        return Err(invalid());
    }
    Ok(())
}

/// Runs bounded static topology entrypoints on shell-vertex B-reps or empty shapes.
/// Header: entry, operation, request, tolerance, best, n1, n2, query xyz, vertices xyz;
/// the final two flags select swapped arguments and inclusion without a shape limit.
///
/// # Errors
/// Malformed/nonfinite/unbounded input, or a reference construction failure.
pub fn eval_static(input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate_static(input)?;
    call_with_buffer(1024, |out, len| {
        // SAFETY: the complete fixed header, bounded counts and finite payload are
        // validated above; the reference copies only within the supplied capacity.
        unsafe {
            sys::sc_intersect_topology_solver_p2_static(
                input.as_ptr(),
                input.len() as i64,
                out.as_mut_ptr(),
                out.len() as i64,
                &raw mut *len,
            )
        }
    })
}

/// Input-only graph and query metadata of one native static fixture B-rep.
#[derive(Debug)]
pub struct StaticBrepFixture {
    /// Shared topology graph, including its success prefix.
    pub graph: Vec<f64>,
    /// Native lists, tolerances, geometry kinds, owner links and face state.
    pub query_metadata: Vec<f64>,
    /// The exact face created for this operand, absent for entries 0–4.
    pub face: Option<u32>,
    /// Original query vertices in input order, as graph registry ids.
    pub vertices: Vec<u32>,
    /// Edge reached by the query curve's real owner chain (entry 7 only).
    pub curve_edge: Option<u32>,
}

fn invalid_fixture() -> FfiError {
    FfiError::InvalidArgument("malformed topology p2 native fixture export".into())
}

struct FixtureReader<'a> {
    values: &'a [f64],
    at: usize,
}
impl FixtureReader<'_> {
    fn number(&mut self) -> Result<f64, FfiError> {
        let value = *self.values.get(self.at).ok_or_else(invalid_fixture)?;
        self.at += 1;
        if !value.is_finite() {
            return Err(invalid_fixture());
        }
        Ok(value)
    }
    fn count(&mut self) -> Result<usize, FfiError> {
        let value = self.number()?;
        if value < 0.0
            || value.fract() != 0.0
            || value > self.values.len().saturating_sub(self.at) as f64
        {
            return Err(invalid_fixture());
        }
        Ok(value as usize)
    }
    fn block(&mut self) -> Result<Vec<f64>, FfiError> {
        let count = self.count()?;
        let end = self.at.checked_add(count).ok_or_else(invalid_fixture)?;
        let block = self
            .values
            .get(self.at..end)
            .ok_or_else(invalid_fixture)?
            .to_vec();
        self.at = end;
        Ok(block)
    }
    fn id(&mut self, count: usize) -> Result<Option<u32>, FfiError> {
        let value = self.number()?;
        if value == -1.0 {
            return Ok(None);
        }
        if value < 0.0
            || value.fract() != 0.0
            || value >= count as f64
            || value > f64::from(u32::MAX)
        {
            return Err(invalid_fixture());
        }
        Ok(Some(value as u32))
    }
    fn brep(&mut self, vertex_count: usize) -> Result<StaticBrepFixture, FfiError> {
        let graph = self.block()?;
        let counts = super::topology_edgeuse_p1::counts(&graph).ok_or_else(invalid_fixture)?;
        let query_metadata = self.block()?;
        let face = self.id(counts[3])?;
        let count = self.count()?;
        if count != vertex_count {
            return Err(invalid_fixture());
        }
        let mut vertices = Vec::with_capacity(count);
        for _ in 0..count {
            vertices.push(self.id(counts[9])?.ok_or_else(invalid_fixture)?);
        }
        let curve_edge = self.id(counts[7])?;
        Ok(StaticBrepFixture {
            graph,
            query_metadata,
            face,
            vertices,
            curve_edge,
        })
    }
}

/// Export the exact native input B-reps used by [`eval_static`], before running
/// any topology solve. The versioned envelope contains two length-prefixed graph
/// and query-metadata pairs, followed by the selected face/vertex/curve-owner ids.
/// It includes no reference solutions, classifications or expected statuses.
///
/// # Errors
/// The same input/setup errors as [`eval_static`], or a malformed export envelope.
pub fn static_fixture(input: &[f64]) -> Result<[StaticBrepFixture; 2], FfiError> {
    validate_static(input)?;
    let values = call_with_buffer(4096, |out, len| {
        // SAFETY: the common validator bounds every constructor input; the shim
        // uses the same constructors as eval_static and checks output capacity.
        unsafe {
            sys::sc_intersect_topology_solver_p2_static_fixture(
                input.as_ptr(),
                input.len() as i64,
                out.as_mut_ptr(),
                out.len() as i64,
                &raw mut *len,
            )
        }
    })?;
    let mut reader = FixtureReader {
        values: &values,
        at: 0,
    };
    if reader.number()? != 1.0 {
        return Err(invalid_fixture());
    }
    let n1 = *input.get(5).ok_or_else(invalid_fixture)? as usize;
    let n2 = *input.get(6).ok_or_else(invalid_fixture)? as usize;
    let first = reader.brep(n1)?;
    let second = reader.brep(n2)?;
    if reader.at != values.len() {
        return Err(invalid_fixture());
    }
    Ok([first, second])
}

/// Whether the native vertex selected as a curve by `EdgeVertex` is a curve.
/// This runtime-class probe avoids invoking the upstream undefined cast.
///
/// # Errors
/// Reference fixture construction failure.
pub fn edge_vertex_operand_is_curve() -> Result<bool, FfiError> {
    // SAFETY: the native helper constructs and frees its own fixed-size fixture
    // and performs a runtime type query, never the invalid cast.
    match unsafe { sys::sc_intersect_topology_solver_p2_edge_vertex_operand() } {
        0 => Ok(false),
        1 => Ok(true),
        code => Err(FfiError::Status(code)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cpp_edge_vertex_selected_operand_is_not_a_curve() {
        assert!(!edge_vertex_operand_is_curve().unwrap());
    }

    #[test]
    fn static_protocol_rejects_unbounded_or_malformed_records() {
        let valid = [
            0.0, 3.0, 1.0, 1e-5, 1e30, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0,
        ];
        for (index, value) in [
            (0, 9.0),
            (1, 19.0),
            (2, 4.0),
            (5, 1e9),
            (5, 1.5),
            (6, -1.0),
            (7, f64::NAN),
            (16, 0.5),
        ] {
            let mut malformed = valid;
            malformed[index] = value;
            assert!(matches!(
                eval_static(&malformed),
                Err(FfiError::InvalidArgument(_))
            ));
        }
        assert!(matches!(
            eval_static(&valid[..17]),
            Err(FfiError::InvalidArgument(_))
        ));
    }
    #[test]
    fn cpp_face_fixtures_have_valid_setup_and_complete() {
        for entry in 5..=8 {
            for included in [0.0, 1.0] {
                let input = [
                    f64::from(entry),
                    3.0,
                    1.0,
                    1e-5,
                    1e30,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    0.0,
                    included,
                ];
                // A normal operation error is distinct from failure to construct
                // the reference fixture, which must fail this test.
                let output = eval_static(&input).unwrap();
                assert!(!output.is_empty());
                if included == 0.0 {
                    assert_eq!(output[0], 0.0);
                }
            }
        }
    }
    #[test]
    fn native_static_export_keeps_created_face_topology_and_query_owner() {
        for entry in 5..=8 {
            let input = [
                f64::from(entry),
                3.0,
                1.0,
                1e-5,
                1e30,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                1.0,
            ];
            let [first, second] = static_fixture(&input).unwrap();
            for brep in [&first, &second] {
                let counts = super::super::topology_edgeuse_p1::counts(&brep.graph).unwrap();
                assert_eq!(brep.vertices.len(), 1);
                assert!(brep.face.is_some());
                assert!(counts[6] >= 4 && counts[7] >= 4 && counts[9] >= 4);
                assert!(brep.query_metadata.len() > 10);
            }
            assert_eq!(first.curve_edge.is_some(), entry == 7);
            assert!(second.curve_edge.is_none());
        }
    }

    #[test]
    fn cpp_empty_object_lists_succeed_in_both_orders() {
        assert_eq!(empty_list_results().unwrap(), [0, 0, 0, 0]);
    }

    #[test]
    fn cpp_distinguishes_external_and_managed_point_nodes() {
        assert_eq!(point_storage_counts().unwrap(), [0, 1]);
    }
}
