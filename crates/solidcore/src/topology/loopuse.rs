//! Loop and loop-use topology in the typed B-rep arenas.
//!
//! Upstream: `source/SMLib/inc/SmLoop.h`, `source/SMLib/inc/SmLoopuse.h`,
//! `source/SMLib/src/SmLoop.cpp`, `source/SMLib/src/SmLoopuse.cpp`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`orientation`], [`set_orientation`] | `SmLoopuse::GetOrientation`, `SetOrientation` |
//! | [`loopuse_type`], [`vertexuse`], [`set_loop`] | `GetLoopuseType`, `GetVertexuse`, `SetLoop` |
//! | [`set_edgeuse`], [`collect_edgeuses`] | `SetEdgeuse`, `CollectEdgeuses` |
//! | [`loop_face`], [`loopuses`], [`is_outer_loop`] | `SmLoop::GetFace`, `GetLoopuses`, `IsOuterLoop` |
//! | [`edgeuses`], [`edges`], [`vertexuses`], [`vertices`] | `SmLoopuse` array getters |
//!
//! Deviations: typed ids replace pointers; nullable array entries remain optional
//! ids. Invalid native dereferences return errors. Attribute lookup is supplied by
//! the caller. The existing arena-bounded CCW traversal preserves native order;
//! `CollectEdgeuses` retains the upstream 100000-step limit and partial writes.
//! No recursion. This is an incomplete unit port pending native certification.

use crate::error::{EntityRef, Result, SolidError};
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::topology::edgeuse::p1;
use crate::topology::{
    Brep, EdgeId, EdgeuseId, EdgeuseOwner, FaceId, LoopId, LoopuseId, LoopuseStart, VertexId,
    VertexuseId,
};

const OP: &str = "topology::loopuse";

/// Native loop-use orientation (`GetOrientation`).
/// # Errors
/// A dangling loop-use id.
pub fn orientation(b: &Brep, lu: LoopuseId) -> Result<OrientType> {
    Ok(b.loopuse(lu)?.orientation)
}

/// Store the orientation independently of the owning face-use orientation.
/// # Errors
/// A dangling loop-use id.
pub fn set_orientation(b: &mut Brep, lu: LoopuseId, value: OrientType) -> Result<()> {
    b.loopuse_mut(lu)?.orientation = value;
    Ok(())
}

/// Native member type (`GetLoopuseType`).
/// # Errors
/// A dangling loop-use id.
pub fn loopuse_type(b: &Brep, lu: LoopuseId) -> Result<i64> {
    Ok(match b.loopuse(lu)?.start {
        LoopuseStart::Unknown => 16999,
        LoopuseStart::Edgeuse(_) => 16005,
        LoopuseStart::Vertexuse(_) => 16006,
    })
}

/// The vertex-use head of a vertex loop-use, or native NULL for other types.
/// # Errors
/// A dangling loop-use id.
pub fn vertexuse(b: &Brep, lu: LoopuseId) -> Result<Option<VertexuseId>> {
    Ok(match b.loopuse(lu)?.start {
        LoopuseStart::Vertexuse(v) => v,
        _ => None,
    })
}

/// Friend write of the loop pointer (`SetLoop`); native assertions are diagnostic.
/// # Errors
/// A dangling loop-use id.
pub fn set_loop(b: &mut Brep, lu: LoopuseId, value: Option<LoopId>) -> Result<()> {
    b.loopuse_mut(lu)?.loop_ = value;
    Ok(())
}

/// Rotate an edge loop-use's head to one of its owned edge uses (`SetEdgeuse`).
/// # Errors
/// Native `SM_ERR` for a wrong loop type, NULL candidate or different owner; dangling ids.
pub fn set_edgeuse(b: &mut Brep, lu: LoopuseId, value: Option<EdgeuseId>) -> Result<()> {
    if matches!(b.loopuse(lu)?.start, LoopuseStart::Edgeuse(_)) {
        if let Some(eu) = value {
            if b.edgeuse(eu)?.owner == EdgeuseOwner::Loopuse(Some(lu)) {
                b.loopuse_mut(lu)?.start = LoopuseStart::Edgeuse(Some(eu));
                return Ok(());
            }
        }
    }
    Err(SolidError::sm_err(OP))
}

/// Collect a CCW edge-use chain, retaining writes before an error.
/// # Errors
/// NULL/dangling links or the native 100000-step limit.
pub fn collect_edgeuses(b: &mut Brep, lu: LoopuseId, start: Option<EdgeuseId>) -> Result<()> {
    let first = start.ok_or_else(|| SolidError::null_pointer(OP))?;
    b.loopuse_mut(lu)?.start = LoopuseStart::Edgeuse(Some(first));
    let mut current = Some(first);
    let mut count = 0;
    loop {
        let eu = current.ok_or_else(|| SolidError::null_pointer(OP))?;
        let record = b.edgeuse_mut(eu)?;
        record.owner = EdgeuseOwner::Loopuse(Some(lu));
        current = record.ccw;
        if current == Some(first) {
            return Ok(());
        }
        count += 1;
        if count > 100000 {
            return Err(SolidError::sm_err(OP));
        }
    }
}

/// Face of the actual primary loop-use (`SmLoop::GetFace`).
/// # Errors
/// Native undefined missing primary/face-use links or dangling ids.
pub fn loop_face(b: &Brep, l: LoopId) -> Result<Option<FaceId>> {
    let lu = b
        .loop_(l)?
        .loopuse
        .ok_or_else(|| SolidError::null_pointer(OP))?;
    let fu = b
        .loopuse(lu)?
        .faceuse
        .ok_or_else(|| SolidError::null_pointer(OP))?;
    Ok(b.faceuse(fu)?.face)
}

/// Primary and opposite loop uses (`SmLoop::GetLoopuses`).
/// # Errors
/// Missing primary loop use or dangling ids.
pub fn loopuses(b: &Brep, l: LoopId) -> Result<(LoopuseId, Option<LoopuseId>)> {
    let first = b
        .loop_(l)?
        .loopuse
        .ok_or_else(|| SolidError::null_pointer(OP))?;
    Ok((first, b.loopuse(first)?.mate))
}

/// Whether the primary loop-use is the first in its face-use list.
/// # Errors
/// Missing required links, empty native list or dangling ids.
pub fn is_outer_loop(b: &Brep, l: LoopId) -> Result<bool> {
    let (lu, _) = loopuses(b, l)?;
    let fu = b
        .loopuse(lu)?
        .faceuse
        .ok_or_else(|| SolidError::null_pointer(OP))?;
    let first = b
        .faceuse(fu)?
        .loopuses
        .first()
        .ok_or_else(|| SolidError::null_pointer(OP))?;
    Ok(*first == lu)
}

/// Ordered edge uses with optional native attribute and connection filters.
/// # Errors
/// Dangling/malformed links or those of the caller's attribute lookup.
pub fn edgeuses<F>(
    b: &Brep,
    lu: LoopuseId,
    attribute: Option<u64>,
    connected: Option<EntityRef>,
    mut has: F,
) -> Result<Vec<EdgeuseId>>
where
    F: FnMut(EntityRef, u64) -> Result<bool>,
{
    let mut out = Vec::new();
    for eu in p1::loopuse_edgeuses(b, lu)? {
        if let Some(id) = attribute {
            if !has(EntityRef::Edgeuse(eu), id)? {
                continue;
            }
        }
        if connected.is_some() && !p1::is_connected_to(b, eu, connected)? {
            continue;
        }
        out.push(eu);
    }
    Ok(out)
}

/// Unique edges in first-seen CCW order, retaining defined NULL array entries.
/// # Errors
/// Dangling links, undefined NULL attribute lookup, or caller lookup errors.
pub fn edges<F>(
    b: &Brep,
    lu: LoopuseId,
    attribute: Option<u64>,
    mut has: F,
) -> Result<Vec<Option<EdgeId>>>
where
    F: FnMut(EntityRef, u64) -> Result<bool>,
{
    let mut out = Vec::new();
    for eu in p1::loopuse_edgeuses(b, lu)? {
        let edge = b.edgeuse(eu)?.edge;
        if let Some(id) = attribute {
            let e = edge.ok_or_else(|| SolidError::null_pointer(OP))?;
            if !has(EntityRef::Edge(e), id)? {
                continue;
            }
        }
        if !out.contains(&edge) {
            out.push(edge);
        }
    }
    Ok(out)
}

/// Vertex uses in CCW order, retaining repetitions and defined NULL entries.
/// # Errors
/// Dangling links, undefined NULL attribute lookup, or caller lookup errors.
pub fn vertexuses<F>(
    b: &Brep,
    lu: LoopuseId,
    attribute: Option<u64>,
    mut has: F,
) -> Result<Vec<Option<VertexuseId>>>
where
    F: FnMut(EntityRef, u64) -> Result<bool>,
{
    let members = match b.loopuse(lu)?.start {
        LoopuseStart::Vertexuse(v) => vec![v],
        _ => p1::loopuse_edgeuses(b, lu)?
            .into_iter()
            .map(|eu| Ok(b.edgeuse(eu)?.vertexuse))
            .collect::<Result<Vec<_>>>()?,
    };
    let mut out = Vec::new();
    for vu in members {
        if let Some(id) = attribute {
            let v = vu.ok_or_else(|| SolidError::null_pointer(OP))?;
            if !has(EntityRef::Vertexuse(v), id)? {
                continue;
            }
        }
        out.push(vu);
    }
    Ok(out)
}

/// Unique vertices in first-seen CCW order with the actual vertex attribute filter.
/// # Errors
/// NULL vertex-use links, dangling ids or caller lookup errors.
pub fn vertices<F>(
    b: &Brep,
    lu: LoopuseId,
    attribute: Option<u64>,
    mut has: F,
) -> Result<Vec<Option<VertexId>>>
where
    F: FnMut(EntityRef, u64) -> Result<bool>,
{
    let mut out = Vec::new();
    for vu in vertexuses(b, lu, None, |_, _| Ok(true))? {
        let vu = vu.ok_or_else(|| SolidError::null_pointer(OP))?;
        let vertex = b.vertexuse(vu)?.vertex;
        if let Some(id) = attribute {
            let v = vertex.ok_or_else(|| SolidError::null_pointer(OP))?;
            if !has(EntityRef::Vertex(v), id)? {
                continue;
            }
        }
        if !out.contains(&vertex) {
            out.push(vertex);
        }
    }
    Ok(out)
}

impl crate::topology::Loopuse {
    /// Native `SmLoopuse_TYPE`.
    #[must_use]
    pub const fn type_code() -> i64 {
        16004
    }
    /// Native class name (`GetClassString`).
    #[must_use]
    pub const fn class_string() -> &'static str {
        "SmLoopuse"
    }
    /// Native type name (`GetTypeString`).
    #[must_use]
    pub const fn type_string() -> &'static str {
        "SmLoopuse_TYPE"
    }
    /// Native class and base-class membership (`IsKindOf`).
    #[must_use]
    pub const fn is_kind_of(code: i64) -> bool {
        code == Self::type_code() || crate::topology::topology::is_kind_of(code)
    }
}
impl crate::topology::Loop {
    /// Native `SmLoop_TYPE` (`GetType`, `GetClassType`).
    #[must_use]
    pub const fn type_code() -> i64 {
        16003
    }
    /// Native class name (`GetClassString`).
    #[must_use]
    pub const fn class_string() -> &'static str {
        "SmLoop"
    }
    /// Native type name (`GetTypeString`).
    #[must_use]
    pub const fn type_string() -> &'static str {
        "SmLoop_TYPE"
    }
    /// Native class and base-class membership (`IsKindOf`).
    #[must_use]
    pub const fn is_kind_of(code: i64) -> bool {
        code == Self::type_code() || crate::topology::topology::is_kind_of(code)
    }
}

/// Whether the loop-use CCW chain closes (`IsClosedPtrs`). Vertex loops are closed.
/// # Errors
/// Dangling links or a non-closing cycle exceeding the arena's edge-use count.
pub fn is_closed_ptrs(b: &Brep, lu: LoopuseId) -> Result<bool> {
    let first = match b.loopuse(lu)?.start {
        LoopuseStart::Vertexuse(_) => return Ok(true),
        LoopuseStart::Edgeuse(value) => value,
        LoopuseStart::Unknown => None,
    };
    let Some(first) = first else {
        return Ok(false);
    };
    let mut current = Some(first);
    for _ in 0..=b.edgeuses().len() {
        let Some(eu) = current else {
            return Ok(false);
        };
        current = b.edgeuse(eu)?.ccw;
        if current == Some(first) {
            return Ok(true);
        }
    }
    Err(SolidError::InvalidInput {
        operation: OP,
        reason: "CCW cycle does not close at its head within the arena bound",
    })
}

/// Loop pointer closure, true for a loop without a primary use (`SmLoop::IsClosedPtrs`).
/// # Errors
/// Dangling ids or an invalid CCW cycle.
pub fn loop_is_closed_ptrs(b: &Brep, l: LoopId) -> Result<bool> {
    match b.loop_(l)?.loopuse {
        Some(lu) => is_closed_ptrs(b, lu),
        None => Ok(true),
    }
}
