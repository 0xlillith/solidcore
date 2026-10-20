//! Validated entry points to the real `SmMerge::ManifoldBoolean`.
//! Fixed-size operands avoid native file readers and input-sized allocations.
//! Coordinates/dimensions are bounded by 1e100 through the shared reader, leaving
//! headroom for endpoint addition and squared-distance products in constructors.
//! This is the oracle's supported numeric domain, not a kernel input restriction.
//! UB-175 [PSIRT] in `IntersectInsertRelate` / `~SmProtoTopologyManager` is not
//! exercised by the current empty/vertex/open-line/plane operand families (no
//! closed edges). Adding circles or closed wires requires rechecking that guard.
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

const MAX_COORD: f64 = 1e100;
crate::reader::unit_reader!(MAX_COORD);

fn validate(operand: &Operand) -> Result<(), FfiError> {
    let mut reader = Reader::new(operand);
    reader.choice(&[-1.0, 0.0, 1.0, 2.0, 3.0])?;
    reader.vec3()?;
    reader.vec3()?;
    reader.done()
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


/// The same native input operand before any Boolean work, including its actual
/// concrete geometry types, list ordering and direct owner identities.
#[derive(Debug, Clone)]
pub struct Fixture {
    /// Shared raw topology/geometry export, including its initial status word.
    pub graph: Vec<f64>,
    /// Native query metadata, decoded by the shared differential-runner adapter.
    pub query_metadata: Vec<f64>,
    /// Actual pre-operation topology flags, stored marks and user fields.
    pub notification_records: Vec<NotificationRecord>,
    /// The native operand context has user or system notification callbacks.
    pub callback_present: bool,
    /// Native topology or directly owned geometry has attached attributes.
    pub attributes_present: bool,
    /// Existing trimmed-surface cache state for every face in graph order.
    /// Observed without building or validating a native cache.
    pub surface_caches: Vec<SurfaceCacheState>,
}

/// Observed native trimmed-surface cache entry before the Boolean operation.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceCacheState {
    /// A cache is already attached to this face's surface.
    pub present: bool,
    /// The attached cache records an edited face; false when no cache exists.
    pub face_was_modified: bool,
}

/// Original operand inputs and their live context states before Boolean work.
#[derive(Debug, Clone)]
pub struct InputPair {
    /// Null operands are absent; each present operand retains its native metadata.
    pub operands: [Option<Fixture>; 2],
    /// Primary and other context states after both constructors. For shared
    /// contexts these are equal, including when an operand itself is null.
    pub contexts: [ContextState; 2],
}

/// Native context metadata captured before Access construction or mark locking.
#[derive(Debug, Clone, PartialEq)]
pub struct ContextState {
    /// Mark, Mark2, Mark3, IO and assertion current-mark counters, in that order.
    pub current_marks: [u64; 5],
    /// Lock flags in the same order as `current_marks`.
    pub mark_locks: [bool; 5],
    /// The context's original Boolean-in-progress flag.
    pub doing_boolean: bool,
    /// Raw model size estimate, including the native default sentinel.
    pub model_size_estimate: f64,
    /// Raw large/small ratio, including the native default sentinel.
    pub large_small_size_ratio: f64,
    /// Stored zone tolerance, without recomputing it from another field.
    pub zone_tol_3d: f64,
}

fn context_input(raw: &[f64]) -> Result<ContextState, FfiError> {
    let invalid = || FfiError::InvalidArgument("ManifoldBoolean input context".into());
    if raw.len() != 20 || raw.first() != Some(&1.0) { return Err(invalid()); }
    let word = |value: f64| -> Result<u32, FfiError> {
        if !value.is_finite() || value < 0.0 || value > f64::from(u32::MAX) || value.fract() != 0.0 { return Err(invalid()); }
        Ok(value as u32)
    };
    let mut current_marks = [0; 5];
    for (mark, pair) in current_marks.iter_mut().zip(raw.get(1..11).ok_or_else(invalid)?.as_chunks::<2>().0) {
        *mark = u64::from(word(pair[0])?) | (u64::from(word(pair[1])?) << 32);
    }
    let mut mark_locks = [false; 5];
    for (locked, value) in mark_locks.iter_mut().zip(raw.get(11..16).ok_or_else(invalid)?) {
        if !matches!(*value, 0.0 | 1.0) { return Err(invalid()); }
        *locked = *value == 1.0;
    }
    if !matches!(*raw.get(16).ok_or_else(invalid)?, 0.0 | 1.0) || raw.get(17..).ok_or_else(invalid)?.iter().any(|v| !v.is_finite()) { return Err(invalid()); }
    Ok(ContextState { current_marks, mark_locks, doing_boolean: *raw.get(16).ok_or_else(invalid)? == 1.0,
        model_size_estimate: *raw.get(17).ok_or_else(invalid)?, large_small_size_ratio: *raw.get(18).ok_or_else(invalid)?, zone_tol_3d: *raw.get(19).ok_or_else(invalid)? })
}

/// Native input metadata keyed by the same graph traversal used for topology.
#[derive(Debug, Clone)]
pub struct NotificationRecord {
    /// 0 B-rep, 1 region, 2 shell, 3 faceuse, 4 face, 5 loopuse, 6 loop,
    /// 7 edgeuse, 8 edge, 9 vertexuse, 10 vertex.
    pub kind: u32,
    /// Index within that graph traversal; the B-rep has index zero.
    pub index: u32,
    /// Flags, five stored marks, user index 1, user index 2 and user pointer 1.
    /// Each native word is encoded through two exact unsigned 32-bit limbs.
    pub values: [u64; 9],
}

fn notification_inputs(raw: &[f64]) -> Result<(Vec<NotificationRecord>, bool, bool), FfiError> {
    let invalid = || FfiError::InvalidArgument("ManifoldBoolean notification input record".into());
    let word = |v: f64| -> Result<u32, FfiError> {
        if !v.is_finite() || v < 0.0 || v > f64::from(u32::MAX) || v.fract() != 0.0 {
            return Err(invalid());
        }
        Ok(v as u32)
    };
    let [version, callback, attributes, count] = raw.get(..4).ok_or_else(invalid)? else { return Err(invalid()); };
    if *version != 1.0 || !matches!(*callback, 0.0 | 1.0) || !matches!(*attributes, 0.0 | 1.0) { return Err(invalid()); }
    let count = usize::try_from(word(*count)?).map_err(|_| invalid())?;
    let expected = count.checked_mul(20).and_then(|n| n.checked_add(4)).ok_or_else(invalid)?;
    if raw.len() != expected { return Err(invalid()); }
    let mut records = Vec::new();
    for row in raw.get(4..).ok_or_else(invalid)?.as_chunks::<20>().0 {
        let kind=word(row[0])?; let index=word(row[1])?;
        if kind > 10 || (kind == 0 && index != 0) { return Err(invalid()); }
        let mut values=[0;9];
        for (value, limbs) in values.iter_mut().zip(row[2..].as_chunks::<2>().0) {
            *value=u64::from(word(limbs[0])?) | (u64::from(word(limbs[1])?) << 32);
        }
        records.push(NotificationRecord { kind,index,values });
    }
    Ok((records,*callback==1.0,*attributes==1.0))
}

fn surface_cache_inputs(raw: &[f64]) -> Result<Vec<SurfaceCacheState>, FfiError> {
    let invalid = || FfiError::InvalidArgument("ManifoldBoolean input surface caches".into());
    let [version,count] = raw.get(..2).ok_or_else(invalid)? else { return Err(invalid()); };
    if *version != 1.0 || !count.is_finite() || *count < 0.0 || count.fract() != 0.0 || *count > raw.len() as f64 { return Err(invalid()); }
    let expected = (*count as usize).checked_mul(2).and_then(|n| n.checked_add(2)).ok_or_else(invalid)?;
    if raw.len() != expected { return Err(invalid()); }
    raw.get(2..).ok_or_else(invalid)?.as_chunks::<2>().0.iter().map(|[present,modified]| {
        if !matches!(*present,0.0 | 1.0) || !matches!(*modified,0.0 | 1.0) || (*present == 0.0 && *modified != 0.0) { return Err(invalid()); }
        Ok(SurfaceCacheState { present:*present == 1.0, face_was_modified:*modified == 1.0 })
    }).collect()
}

fn input_fixture(raw: &[f64]) -> Result<Fixture, FfiError> {
    let invalid = || FfiError::InvalidArgument("ManifoldBoolean input metadata record".into());
    let take = |at: &mut usize| -> Result<Vec<f64>, FfiError> {
        let value = *raw.get(*at).ok_or_else(invalid)?;
        if !value.is_finite() || value < 0.0 || value.fract() != 0.0 || value > raw.len() as f64 {
            return Err(invalid());
        }
        *at = at.checked_add(1).ok_or_else(invalid)?;
        let end = at.checked_add(value as usize).ok_or_else(invalid)?;
        let values = raw.get(*at..end).ok_or_else(invalid)?.to_vec();
        *at = end;
        Ok(values)
    };
    let mut at = 0;
    let graph = take(&mut at)?;
    let query_metadata = take(&mut at)?;
    let notifications = take(&mut at)?;
    let surface_caches = surface_cache_inputs(&take(&mut at)?)?;
    let (notification_records,callback_present,attributes_present)=notification_inputs(&notifications)?;
    if at != raw.len() || graph.is_empty() || query_metadata.is_empty() {
        return Err(invalid());
    }
    Ok(Fixture { graph, query_metadata, notification_records, callback_present, attributes_present, surface_caches })
}

fn input_pair(raw: &[f64]) -> Result<InputPair, FfiError> {
    let invalid = || FfiError::InvalidArgument("ManifoldBoolean input metadata pair".into());
    if raw.first() != Some(&2.0) { return Err(invalid()); }
    let mut at: usize = 1;
    let mut next_context = || -> Result<ContextState, FfiError> {
        if raw.get(at) != Some(&20.0) { return Err(invalid()); }
        at += 1;
        let end = at.checked_add(20).ok_or_else(invalid)?;
        let result = context_input(raw.get(at..end).ok_or_else(invalid)?)?;
        at = end;
        Ok(result)
    };
    let contexts = [next_context()?, next_context()?];
    let mut result = [None, None];
    for slot in &mut result {
        let present = *raw.get(at).ok_or_else(invalid)?;
        at = at.checked_add(1).ok_or_else(invalid)?;
        match present {
            0.0 => {},
            1.0 => {
                let count = *raw.get(at).ok_or_else(invalid)?;
                if !count.is_finite() || count < 0.0 || count.fract() != 0.0 || count > raw.len() as f64 {
                    return Err(invalid());
                }
                at = at.checked_add(1).ok_or_else(invalid)?;
                let end = at.checked_add(count as usize).ok_or_else(invalid)?;
                *slot = Some(input_fixture(raw.get(at..end).ok_or_else(invalid)?)?);
                at = end;
            },
            _ => return Err(invalid()),
        }
    }
    if at != raw.len() { return Err(invalid()); }
    Ok(InputPair { operands: result, contexts })
}

/// Export both original operands in the same construction order and contexts as
/// `eval`, before Boolean work. Null inputs remain absent. A native constructor
/// rejection is an error here; `eval` classifies rejection before the Rust run.
pub fn export_fixtures(a: &Operand, b: &Operand, distinct_contexts: bool)
    -> Result<InputPair, FfiError> {
    input_pair(&call(-2, a, b, Request {
        operation: 0, flags: [false; 5], manifold: true,
        distinct_contexts, caller_lock: false,
    })?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn input_metadata_lengths_are_checked_before_slicing() {
        for raw in [vec![-1.0], vec![f64::NAN], vec![1.5], vec![1e300],
                    vec![2.0, 0.0], vec![1.0, 0.0, 3.0, 1.0],
                    vec![1.0, 0.0, 1.0, 1.0, 9.0]] {
            assert!(input_fixture(&raw).is_err());
        }
    }

    #[test]
    fn input_pair_rejects_bad_presence_lengths_and_trailing_values() {
        let mut context = vec![1.0];
        for _ in 0..5 { context.extend([f64::from(u32::MAX), f64::from(u32::MAX)]); }
        context.extend([1.0,0.0,1.0,0.0,1.0,1.0,-1.0,-1.0,1e-5]);
        let mut pair = vec![2.0];
        for _ in 0..2 { pair.push(20.0); pair.extend(&context); }
        pair.extend([0.0,0.0]);
        let imported = input_pair(&pair).unwrap();
        assert!(imported.operands.iter().all(Option::is_none));
        assert_eq!(imported.contexts[0].current_marks, [u64::MAX;5]);
        assert_eq!(imported.contexts[0].mark_locks, [true,false,true,false,true]);
        assert!(imported.contexts[0].doing_boolean);
        assert_eq!(imported.contexts[0], imported.contexts[1]);
        pair.push(0.0);
        assert!(input_pair(&pair).is_err());
        context[11] = 2.0;
        assert!(context_input(&context).is_err());
        context[11] = 1.0;
        context[1] = f64::from(u32::MAX) + 1.0;
        assert!(context_input(&context).is_err());
        context[1] = 0.0;
        context[19] = f64::NAN;
        assert!(context_input(&context).is_err());
        for raw in [vec![-1.0], vec![2.0, 0.0, 0.0], vec![1.0, 2.0],
                    vec![1.0, 1.0, f64::NAN], vec![1.0, 1.0, 1e300],
                    vec![1.0, 0.0, 0.0, 9.0]] {
            assert!(input_pair(&raw).is_err());
        }
    }

    #[test]
    fn notification_input_preserves_full_words_and_rejects_invalid_limbs() {
        let mut raw=vec![1.0,0.0,0.0,1.0,10.0,0.0];
        raw.extend([f64::from(u32::MAX);18]);
        let (records,callback,attributes)=notification_inputs(&raw).unwrap();
        assert_eq!(records[0].values,[u64::MAX;9]);
        assert!(!callback && !attributes);
        raw[6]=f64::from(u32::MAX)+1.0;
        assert!(notification_inputs(&raw).is_err());
        raw[6]=f64::NAN;
        assert!(notification_inputs(&raw).is_err());
        assert!(notification_inputs(&raw[..raw.len()-1]).is_err());
    }

    #[test]
    fn operand_coordinate_bound_precedes_native_construction() {
        for value in [1e308, -1e308, f64::INFINITY, f64::NAN] {
            let mut operand = [2.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
            operand[1] = value;
            assert!(validate(&operand).is_err());
            operand[1] = 0.0;
            operand[4] = value;
            assert!(validate(&operand).is_err());
        }
        assert!(validate(&[2.0, MAX_COORD, -MAX_COORD, 0.0, MAX_COORD, 0.0, 0.0]).is_ok());
    }
}
