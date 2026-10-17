//! B-rep editing, part 5: coincident-vertex targets and dependency-limited edits.
//!
//! Upstream: `source/SMLib/src/SmBrep.cpp`, `source/SMLib/inc/SmBrep.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`coincident_vertex_target`], [`combine_coincident_vertices`] | `SmBrep::CombineCoincidentVertices` |
//! | [`create_wire_edge_from_curve`] | `SmBrep::CreateWireEdgeFromCurve` |
//! | [`delete_topological_vertex`] | `SmBrep::DeleteTopologicalVertex` |
//! | [`find_and_set_infinite_region`] | `SmBrep::FindAndSetInfiniteRegion` |
//! | [`glue_vertices`] | `SmBrep::GlueVertices` |
//! | [`make_edge_in_face`] | `SmBrep::MakeEdgeInFace` |
//! | [`make_manifold`] | `SmBrep::MakeManifold` |
//! | [`remove_topological_edges_and_vertices`] | `SmBrep::RemoveTopologicalEdgesAndVertices` |
//! | [`remove_faces`] | `SmBrep::RemoveFaces` |
//! | [`stitch_and_orient`] | `SmBrep::StitchAndOrient` |
//! | [`squeeze_face`] | `SmBrep::SqueezeFace` |
//! | [`transform`] | `SmBrep::Transform` |
//!
//! Deviations: **review WIP, plan status `todo`**. Complete control flow is present
//! for coincident-vertex combining, wire creation, infinite-region reset, face-edge
//! insertion, manifold/face cleanup and stitching/orientation, with per-call
//! dependency barriers. The remaining omitted bodies are identified at their
//! early returns below; those continuations must be written before back-edges
//! can be wired or this unit can be counted as ported.
//! `SmBrep::Notify` is owned by `topology.brep_p2`, not this unit. Owning-list
//! traversal and editing state belong to `topology.brep_p1`: arena allocation order
//! is not substituted for list order. Callers supply the actual editing flag where
//! it is tested before any missing call. Pointer inputs are optional typed ids;
//! out-parameters are mutable outputs where failure must clear them. The old-tol
//! target calculation is factored out for direct numerical-stage verification;
//! that verification does not certify vertex deletion or gluing. The internal
//! `SmPendingTopologyDelete` RAII helper belongs to `MakeTopologyFromData`'s pending
//! allocations; it needs no stand-alone arena API (Rust drops owned inputs on error).
//! Attribute detachment and topology insertion remain that reader's responsibility.
//! Aliased `CombineCoincidentVertices` is refused only for an isolated vertex,
//! where native deletion would invalidate its survivor. With incident edges it
//! retains the native path. Explicit editing arguments are temporary: part 1's
//! owning record/accessor will replace them without assuming that editing is on.
//! Debug drawing/assertions are omitted. No recursion is introduced.

use std::collections::BTreeMap;

use crate::curve::curve::Curve;
use crate::error::{Result, SolidError};
use crate::geom::attribute::{AttributeStore, ObjectKey};
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::context::{MarkLock, ALL_MARKS};
use crate::geom::extent1d::Extent1d;
use crate::geom::extent3d::Extent3d;
use crate::geom::object::ObjectRef;
use crate::geom::tol::ZONE_TOL_3D;
use crate::kernel::Kernel;
use crate::math::{max2, Point3, Vec3};
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::topology::topology::Topology;
use crate::topology::VertexuseId;
use crate::topology::{vertex, Brep, EdgeId, FaceId, LoopId, RegionId, ShellId, VertexId};

const OP: &str = "topology::brep::p5";

fn missing<T>(symbol: &'static str, unit: &'static str) -> Result<T> {
    Err(SolidError::unported(symbol, unit))
}

fn require_editing(enabled: bool) -> Result<()> {
    if enabled {
        Ok(())
    } else {
        Err(SolidError::sm_err(OP))
    }
}

/// Reads the owning B-rep's editing member once its record support lands. Unlike
/// the explicit-state entry points, opaque existing callers cannot assume TRUE.
pub(crate) fn editing_enabled(_b: &Brep) -> Result<bool> {
    // UNPORTED(topology.brep_p1): SmBrep::m_bEditingEnabled
    missing("SmBrep::m_bEditingEnabled", "topology.brep_p1")
}

/// The position, tolerance and temporary-edge length threshold calculated before
/// `CombineCoincidentVertices` modifies topology (`SM_USE_OLDTOL`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoincidentVertexTarget {
    /// Position of the vertex with more incident edges; ties keep the survivor.
    pub point: Point3,
    /// The tolerance assigned after combining the vertices.
    pub tolerance: f64,
    /// `dMinDist`, used when deciding whether a common edge is short enough to squeeze.
    pub minimum_edge_length: f64,
}

/// Computes the initial old-tolerance stage of `CombineCoincidentVertices`.
/// Does not edit either vertex. Distances and the strict edge-count/tolerance
/// comparisons retain the native operation order.
pub fn coincident_vertex_target(
    b: &Brep,
    survivor: VertexId,
    deleted: VertexId,
) -> Result<CoincidentVertexTarget> {
    let edges = vertex::edges(b, survivor, None)?;
    let edges2 = vertex::edges(b, deleted, None)?;
    let s = b.vertex(survivor)?;
    let d = b.vertex(deleted)?;
    let mut point = s.point;
    let mut min_dist = point.distance(d.point);
    if min_dist > d.tolerance {
        let doubled = 2.0 * min_dist;
        min_dist = doubled;
    } else {
        min_dist += d.tolerance;
    }
    let mut tolerance = max2(s.tolerance, min_dist);
    if edges2.len() > edges.len() {
        min_dist = point.distance(d.point);
        if min_dist > s.tolerance {
            let doubled = 2.0 * min_dist;
            min_dist = doubled;
        } else {
            min_dist += s.tolerance;
        }
        point = d.point;
        tolerance = max2(d.tolerance, min_dist);
    }
    Ok(CoincidentVertexTarget {
        point,
        tolerance,
        minimum_edge_length: min_dist,
    })
}

/// Combines coincident vertices. Clears the deleted-edge output before reading
/// either vertex, including failure paths. The target stage is available, but
/// successful deletion and gluing remain dependency limited.
pub fn combine_coincident_vertices(
    k: &mut Kernel,
    b: &mut Brep,
    survivor: VertexId,
    deleted: VertexId,
    deleted_edge: &mut Option<EdgeId>,
) -> Result<()> {
    combine_with_notifications(k, b, survivor, deleted, deleted_edge, None)
}

/// The caller's real notification records and attribute store for vertex edits.
/// No replacement records or context are synthesized by the editing operation.
pub struct VertexNotifyEnv<'a> {
    /// Topology records addressed by the caller's object keys.
    pub records: &'a mut BTreeMap<ObjectKey, Topology>,
    /// Attributes attached to those objects.
    pub attributes: &'a mut AttributeStore,
    /// Maps a vertex identity to its notified object.
    pub vertex_ref: &'a dyn Fn(VertexId) -> ObjectRef,
    /// Maps a vertex-use identity to its notified object.
    pub vertexuse_ref: &'a dyn Fn(VertexuseId) -> ObjectRef,
}

/// `CombineCoincidentVertices` with the actual notification environment needed
/// by `SmVertex::SetPoint`. The arena-only entry defers when it reaches this
/// environment boundary; supplying it never bypasses an unported native callee.
pub fn combine_with_notifications(
    k: &mut Kernel,
    b: &mut Brep,
    survivor: VertexId,
    deleted: VertexId,
    deleted_edge: &mut Option<EdgeId>,
    mut notifications: Option<&mut VertexNotifyEnv<'_>>,
) -> Result<()> {
    *deleted_edge = None;
    let target = coincident_vertex_target(b, survivor, deleted)?;
    let edges = vertex::edges(b, survivor, None)?;
    let edges2 = vertex::edges(b, deleted, None)?;
    if survivor == deleted && edges.is_empty() {
        return Err(SolidError::invalid_input(
            OP,
            "isolated aliased merge would delete its survivor",
        ));
    }
    if edges2.is_empty() {
        notify_merge(k, b, survivor, deleted, notifications.as_deref_mut())?;
        set_vertex_point(k, b, survivor, target.point, notifications.as_deref_mut())?;
        vertex::set_tolerance(b, survivor, target.tolerance, true, true)?;
        crate::error::ignore_status(delete_vertex(b, deleted))?;
        return Ok(());
    } else if edges.is_empty() {
        let first_uses = vertex::vertexuses(b, survivor)?;
        let last_uses = vertex::vertexuses(b, deleted)?;
        for vu in first_uses {
            crate::error::ignore_status(remove_vertexuse(b, survivor, vu))?;
            crate::error::ignore_status(post_insert_vertexuse(b, deleted, vu))?;
        }
        for vu in last_uses {
            crate::error::ignore_status(remove_vertexuse(b, deleted, vu))?;
            crate::error::ignore_status(post_insert_vertexuse(b, survivor, vu))?;
        }
        notify_merge(k, b, survivor, deleted, notifications.as_deref_mut())?;
        set_vertex_point(k, b, survivor, target.point, notifications.as_deref_mut())?;
        vertex::set_tolerance(b, survivor, target.tolerance, true, true)?;
        crate::error::ignore_status(delete_vertex(b, deleted))?;
        return Ok(());
    }
    let mut common = None;
    for e in edges {
        if crate::topology::edge::other_vertex(b, e, Some(survivor))? != Some(deleted) {
            continue;
        }
        let eus = &b.edge(e)?.edgeuses;
        let lamina = if let [head, _] = eus.as_slice() {
            crate::topology::edgeuse::p1::is_loop_edgeuse(b, *head)?
        } else {
            false
        };
        if crate::topology::edge::is_manifold(b, e)? || lamina {
            common = Some(e);
            break;
        }
    }
    if let Some(e) = common {
        let record = b.edge(e)?;
        let curve = record
            .curve
            .as_ref()
            .ok_or_else(|| SolidError::invalid_input(OP, "common edge without curve"))?;
        let length = crate::error::ignored(curve.length(k, &record.interval, target.tolerance))?
            .unwrap_or(0.0);
        if length > target.minimum_edge_length {
            // Retain a long common edge; GlueVertices will close it.
            common = None;
        } else {
            *deleted_edge = Some(e);
        }
        // A long common edge is retained as a closed edge by GlueVertices.
    } else {
        let faces = vertex::faces(b, survivor)?;
        let faces2 = vertex::faces(b, deleted)?;
        if let Some(face) = faces.into_iter().find(|f| faces2.contains(f)) {
            for vu in vertex::vertexuses(b, survivor)? {
                let faceuse = crate::topology::vertexuse::faceuse(b, vu)?;
                let owning_face = match faceuse {
                    Some(fu) => b.faceuse(fu)?.face,
                    None => None,
                };
                if owning_face == Some(face)
                    && matches!(
                        b.vertexuse(vu)?.owner,
                        crate::topology::VertexuseOwner::Loopuse(_)
                    )
                {
                    return Err(SolidError::sm_err(OP));
                }
            }
            let start = b.vertex(survivor)?.point;
            let mut end = b.vertex(deleted)?.point;
            if crate::geom::tol::in_tol(
                start.distance(end),
                crate::geom::tol::scaled_zero2(&start, &end),
            ) {
                end = start + Vec3::new(target.tolerance, target.tolerance, target.tolerance);
            }
            // Preserve failures from the temporary line before the missing edit
            // metadata is consulted. The owned temporary drops on deferral.
            let line = BSplineCurve::create_line_segment(3, start, end)?;
            let interval = line.natural_interval();
            let editing = editing_enabled(b)?;
            let mut output = EdgeInFace::default();
            make_edge_in_face(
                b,
                editing,
                Some(face),
                Some(survivor),
                Some(deleted),
                Curve::from(line),
                None,
                interval,
                OrientType::Same,
                crate::math::EFF_ZERO,
                &mut output,
            )?;
            common = output.edge;
        }
    }
    if let Some(edge) = common {
        squeeze_edge(b, edge, survivor)?;
    } else {
        crate::error::ignore_status(glue_vertices(b, Some(survivor), Some(deleted), None))?;
    }
    set_vertex_point(k, b, survivor, target.point, notifications.as_deref_mut())?;
    vertex::set_tolerance(b, survivor, target.tolerance, true, true)?;
    update_vertex_edge_intervals(b, survivor)
}

fn notify_merge(
    _k: &mut Kernel,
    _b: &mut Brep,
    _survivor: VertexId,
    _deleted: VertexId,
    _notifications: Option<&mut VertexNotifyEnv<'_>>,
) -> Result<()> {
    // UNPORTED(topology.brep_p2): SmBrep::Notify
    missing("SmBrep::Notify", "topology.brep_p2")
}

fn delete_vertex(_b: &mut Brep, _vertex: VertexId) -> Result<()> {
    // UNPORTED(topology.brep_p4): SmBrep::DeleteVertex
    missing("SmBrep::DeleteVertex", "topology.brep_p4")
}

fn squeeze_edge(_b: &mut Brep, _edge: EdgeId, _survivor: VertexId) -> Result<()> {
    // UNPORTED(topology.brep_p3): SmBrep::SqueezeEdge
    missing("SmBrep::SqueezeEdge", "topology.brep_p3")
}

fn remove_vertexuse(b: &mut Brep, vertex: VertexId, vu: VertexuseId) -> Result<()> {
    if b.vertexuse(vu)?.vertex != Some(vertex) {
        return Err(SolidError::sm_err(OP));
    }
    let index = b
        .vertex(vertex)?
        .vertexuses
        .iter()
        .position(|&id| id == vu)
        .ok_or_else(|| SolidError::sm_err(OP))?;
    b.vertex_mut(vertex)?.vertexuses.remove(index);
    b.vertexuse_mut(vu)?.vertex = None;
    Ok(())
}

fn post_insert_vertexuse(b: &mut Brep, vertex: VertexId, vu: VertexuseId) -> Result<()> {
    b.vertex(vertex)?;
    if b.vertexuse(vu)?.vertex.is_some() {
        return Err(SolidError::sm_err(OP));
    }
    b.vertex_mut(vertex)?.vertexuses.push(vu);
    b.vertexuse_mut(vu)?.vertex = Some(vertex);
    Ok(())
}

fn set_vertex_point(
    k: &mut Kernel,
    b: &mut Brep,
    vertex: VertexId,
    point: Point3,
    notifications: Option<&mut VertexNotifyEnv<'_>>,
) -> Result<()> {
    let Some(env) = notifications else {
        // UNPORTED(topology.brep_p2): SmBrep::Notify (arena-only notification environment)
        return missing("SmBrep::Notify", "topology.brep_p2");
    };
    vertex::set_point(
        b,
        vertex,
        point,
        env.records,
        env.attributes,
        Some(&mut k.context),
        (env.vertex_ref)(vertex),
        env.vertexuse_ref,
    )
}

fn update_vertex_edge_intervals(_b: &mut Brep, _vertex: VertexId) -> Result<()> {
    // UNPORTED(topology.brep_p1): SmBrep::UpdateVertexEdgeIntervals
    missing("SmBrep::UpdateVertexEdgeIntervals", "topology.brep_p1")
}

/// Makes a wire edge over the caller's curve interval. Evaluates start before end,
/// then needs the owning B-rep's infinite region; the owned curve drops on error.
pub fn create_wire_edge_from_curve(
    _k: &Kernel,
    b: &mut Brep,
    curve: Curve,
    interval: Extent1d,
) -> Result<EdgeId> {
    let start = curve.evaluate_point(interval.min())?;
    let end = curve.evaluate_point(interval.max())?;
    let region = infinite_region(b)?.ok_or_else(|| SolidError::null_pointer(OP))?;
    let (_, first) = make_shell_vertex(b, region, start)?;
    let (edge, last) =
        make_wire_edge_vertex(b, region, first, curve, interval, OrientType::Same, end)?;
    if start.distance(end) < b.tolerance() {
        crate::error::ignore_status(glue_vertices(b, Some(first), Some(last), None))?;
    }
    Ok(edge)
}

fn infinite_region(_b: &Brep) -> Result<Option<RegionId>> {
    // UNPORTED(topology.brep_p1): SmBrep::GetInfiniteRegion
    missing("SmBrep::GetInfiniteRegion", "topology.brep_p1")
}

fn make_shell_vertex(_b: &mut Brep, _r: RegionId, _point: Point3) -> Result<(ShellId, VertexId)> {
    // UNPORTED(topology.brep_p9): SmBrep::MakeShellVertex
    missing("SmBrep::MakeShellVertex", "topology.brep_p9")
}

fn make_wire_edge_vertex(
    _b: &mut Brep,
    _r: RegionId,
    _first: VertexId,
    _curve: Curve,
    _interval: Extent1d,
    _orientation: OrientType,
    _end: Point3,
) -> Result<(EdgeId, VertexId)> {
    // UNPORTED(topology.brep_p3): SmBrep::MakeWireEdgeVertex
    missing("SmBrep::MakeWireEdgeVertex", "topology.brep_p3")
}

/// Deletes a degree-two topological vertex, optionally selecting the shorter edge.
/// The native editing-disabled refusal occurs before Notify and input checks.
pub fn delete_topological_vertex(
    _k: &Kernel,
    _b: &mut Brep,
    editing_enabled: bool,
    _vertex: Option<VertexId>,
    _edge: Option<EdgeId>,
) -> Result<()> {
    require_editing(editing_enabled)?;
    // OMITTED WIP: SmBrep.cpp:7616-7859 (degree checks, geometry joining,
    // UV-curve cleanup, merge notification, squeeze and UV reconstruction).
    // This is not a completed port: replacing Notify alone is insufficient.
    // UNPORTED(topology.brep_p2): SmBrep::Notify
    missing("SmBrep::Notify", "topology.brep_p2")
}

/// Finds the infinite region, then resets its void state. The ray-selection
/// operation belongs to part 9; its optional shell and bounding box are retained
/// in the interface for that integration. `void_setting_code` is the native
/// `SmSetRegionVoidType` code (0 nested-solid, 1 infinite-only, 2 swap). Its common
/// enum is being ported in part 1's independent PR; this deferred boundary does
/// not duplicate that type.
pub fn find_and_set_infinite_region(
    b: &mut Brep,
    shell: Option<ShellId>,
    mut void_setting_code: i32,
    bbox: Option<&mut Extent3d>,
) -> Result<()> {
    let region =
        find_infinite_region(b, shell, bbox)?.ok_or_else(|| SolidError::null_pointer(OP))?;
    if let Some(shell) = shell {
        if shell_is_faceuse_shell(b, shell)? {
            void_setting_code = 2; // SM_VS_INFINITEREGION_SWAP
        }
    }
    reset_infinite_region(b, region, void_setting_code)
}

fn find_infinite_region(
    _b: &mut Brep,
    _shell: Option<ShellId>,
    _bbox: Option<&mut Extent3d>,
) -> Result<Option<RegionId>> {
    // UNPORTED(topology.brep_p9): SmBrep::FindInfiniteRegion
    missing("SmBrep::FindInfiniteRegion", "topology.brep_p9")
}

fn shell_is_faceuse_shell(_b: &Brep, _shell: ShellId) -> Result<bool> {
    // UNPORTED(topology.shell): SmShell::IsFaceuseShell
    missing("SmShell::IsFaceuseShell", "topology.shell")
}

fn reset_infinite_region(_b: &mut Brep, _region: RegionId, _mode: i32) -> Result<()> {
    // UNPORTED(topology.brep_p1): SmBrep::ResetInfiniteRegion
    missing("SmBrep::ResetInfiniteRegion", "topology.brep_p1")
}

/// Glues vertices, leaving the survivor in place. `wire_edges` is the optional
/// caller-supplied wire list, rather than an empty-list substitute for no list.
pub fn glue_vertices(
    _b: &mut Brep,
    survivor: Option<VertexId>,
    deleted: Option<VertexId>,
    _wire_edges: Option<&[EdgeId]>,
) -> Result<()> {
    // OMITTED WIP: SmBrep.cpp:9599-9797 (notifications, tolerances, vertex-use
    // transfer, list removal and shell merging). The equal-input exit is complete;
    // nonaliases must not become successful merely by wiring the entry notification.
    glue_vertices_entry(survivor == deleted)
}

/// Common entry barrier for arena and opaque-handle callers. The first native
/// action is the B-rep notification, before dereferencing either vertex; no
/// temporary or substitute B-rep is constructed for handle-only adapters.
pub(crate) fn glue_vertices_entry(aliases: bool) -> Result<()> {
    if aliases {
        return Ok(());
    }
    // UNPORTED(topology.brep_p2): SmBrep::Notify
    missing("SmBrep::Notify", "topology.brep_p2")
}

/// Outputs of `MakeEdgeInFace`, initialized to null after its editing and face checks.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct EdgeInFace {
    /// The inserted edge.
    pub edge: Option<EdgeId>,
    /// The new loop when insertion splits an existing loop.
    pub loop_: Option<LoopId>,
    /// The new face when insertion splits the face.
    pub face: Option<FaceId>,
}

/// Finds both vertex sectors then inserts an edge into a face. The native output
/// lifecycle is preserved: editing-disabled/null-face failures leave outputs alone;
/// a valid face clears all three before sector classification.
pub fn make_edge_in_face(
    b: &mut Brep,
    editing_enabled: bool,
    face: Option<FaceId>,
    start: Option<VertexId>,
    end: Option<VertexId>,
    curve: Curve,
    uv_curve: Option<BSplineCurve>,
    interval: Extent1d,
    orientation: OrientType,
    max_gap_3d: f64,
    output: &mut EdgeInFace,
) -> Result<()> {
    require_editing(editing_enabled)?;
    let f = face.ok_or_else(|| SolidError::null_pointer(OP))?;
    b.face(f)?;
    *output = EdgeInFace::default();
    let (t_start, t_end, end_orientation) = if orientation == OrientType::Opposite {
        (interval.max(), interval.min(), OrientType::Same)
    } else {
        (interval.min(), interval.max(), OrientType::Opposite)
    };
    let first = find_vertex_sector(
        b,
        f,
        start,
        &curve,
        uv_curve.as_ref(),
        t_start,
        t_end,
        orientation,
        b.tolerance(),
        false,
    )?
    .ok_or_else(|| SolidError::null_pointer(OP))?;
    let last = find_vertex_sector(
        b,
        f,
        end,
        &curve,
        uv_curve.as_ref(),
        t_end,
        t_start,
        end_orientation,
        b.tolerance(),
        false,
    )?
    .ok_or_else(|| SolidError::null_pointer(OP))?;
    make_edge_in_face_with_vu(
        b,
        f,
        first,
        last,
        curve,
        uv_curve,
        interval,
        orientation,
        max_gap_3d,
        output,
    )
}

fn find_vertex_sector(
    _b: &Brep,
    _f: FaceId,
    _v: Option<VertexId>,
    _curve: &Curve,
    _uv: Option<&BSplineCurve>,
    _start: f64,
    _end: f64,
    _orientation: OrientType,
    _tolerance: f64,
    _all_sectors: bool,
) -> Result<Option<VertexuseId>> {
    // UNPORTED(topology.face_p4): SmFace::FindVertexSector
    missing("SmFace::FindVertexSector", "topology.face_p4")
}

fn make_edge_in_face_with_vu(
    _b: &mut Brep,
    _f: FaceId,
    _first: VertexuseId,
    _last: VertexuseId,
    _curve: Curve,
    _uv: Option<BSplineCurve>,
    _interval: Extent1d,
    _orientation: OrientType,
    _max_gap_3d: f64,
    _output: &mut EdgeInFace,
) -> Result<()> {
    // UNPORTED(topology.brep_p9): SmBrep::MakeEdgeInFaceWithVU
    missing("SmBrep::MakeEdgeInFaceWithVU", "topology.brep_p9")
}

/// Makes the selected regions manifold. `None` means all non-infinite regions;
/// `Some(&[])` means no regions and returns the native `SM_ERR`.
pub fn make_manifold(
    b: &mut Brep,
    regions_to_keep: Option<&[RegionId]>,
    keep_interior_faces: bool,
) -> Result<()> {
    // The temporary editing change has no lasting effect on this early return.
    if regions_to_keep.is_some_and(<[RegionId]>::is_empty) {
        return Err(SolidError::sm_err(OP));
    }
    let previous = set_editing_enabled(b, true)?;
    let result = (|| {
        let mut regions = if let Some(given) = regions_to_keep {
            given.to_vec()
        } else {
            let mut all = brep_regions(b)?;
            let infinite = infinite_region(b)?;
            if let Some(index) = all.iter().position(|r| Some(*r) == infinite) {
                all.remove(index);
            }
            all
        };
        if regions.is_empty() {
            return Err(SolidError::sm_err(OP));
        }
        for face in brep_faces(b)? {
            let (first, last) = face_shells(b, face)?;
            let r1 = b
                .shell(first)?
                .region
                .ok_or_else(|| SolidError::null_pointer(OP))?;
            let r2 = b
                .shell(last)?
                .region
                .ok_or_else(|| SolidError::null_pointer(OP))?;
            let keep1 = regions.contains(&r1);
            let keep2 = regions.contains(&r2);
            if r1 != r2 && (keep1 != keep2 || (keep_interior_faces && keep1)) {
                continue;
            }
            delete_face(b, face, false, true, None)?;
        }
        for edge in brep_edges(b)? {
            if crate::topology::edge::is_wire(b, edge)? {
                delete_edge(b, edge)?;
            }
        }
        for vertex in brep_vertices(b)? {
            if vertex::is_shell_vertex(b, vertex)? {
                delete_vertex(b, vertex)?;
            }
        }
        regions = brep_regions(b)?;
        for region in regions {
            let is_infinite = Some(region) == infinite_region(b)?;
            set_region_is_void(b, region, is_infinite)?;
        }
        Ok(())
    })();
    set_editing_enabled(b, previous)?;
    result
}

fn set_editing_enabled(_b: &mut Brep, _enabled: bool) -> Result<bool> {
    // UNPORTED(topology.brep_p1): SmBrep::m_bEditingEnabled
    missing("SmBrep::m_bEditingEnabled", "topology.brep_p1")
}

fn brep_regions(_b: &Brep) -> Result<Vec<RegionId>> {
    // UNPORTED(topology.brep_p1): SmBrep::GetRegions
    missing("SmBrep::GetRegions", "topology.brep_p1")
}
fn brep_faces(_b: &Brep) -> Result<Vec<FaceId>> {
    // UNPORTED(topology.brep_p1): SmBrep::GetFaces
    missing("SmBrep::GetFaces", "topology.brep_p1")
}
fn brep_edges(_b: &Brep) -> Result<Vec<EdgeId>> {
    // UNPORTED(topology.brep_p1): SmBrep::GetEdges
    missing("SmBrep::GetEdges", "topology.brep_p1")
}
fn brep_vertices(_b: &Brep) -> Result<Vec<VertexId>> {
    // UNPORTED(topology.brep_p1): SmBrep::GetVertices
    missing("SmBrep::GetVertices", "topology.brep_p1")
}
fn face_shells(b: &Brep, face: FaceId) -> Result<(ShellId, ShellId)> {
    let (a, z) = crate::topology::face::p2::faceuses(b, face)?;
    let a = a.ok_or_else(|| SolidError::null_pointer(OP))?;
    let z = z.ok_or_else(|| SolidError::null_pointer(OP))?;
    Ok((
        b.faceuse(a)?
            .shell
            .ok_or_else(|| SolidError::null_pointer(OP))?,
        b.faceuse(z)?
            .shell
            .ok_or_else(|| SolidError::null_pointer(OP))?,
    ))
}
fn delete_face(
    _b: &mut Brep,
    _face: FaceId,
    _connected: bool,
    _nest: bool,
    _keep_region: Option<RegionId>,
) -> Result<()> {
    // UNPORTED(topology.brep_p3): SmBrep::DeleteFace
    missing("SmBrep::DeleteFace", "topology.brep_p3")
}
fn delete_edge(_b: &mut Brep, _edge: EdgeId) -> Result<()> {
    // UNPORTED(topology.brep_p4): SmBrep::DeleteEdge
    missing("SmBrep::DeleteEdge", "topology.brep_p4")
}
fn set_region_is_void(_b: &mut Brep, _region: RegionId, _is_void: bool) -> Result<()> {
    // UNPORTED(topology.shell): SmRegion::SetIsVoid
    missing("SmRegion::SetIsVoid", "topology.shell")
}

/// Removes redundant edges/vertices, preserving optional keep lists. The native
/// owning-list traversal is required even for apparently empty arena contents.
pub fn remove_topological_edges_and_vertices(
    _k: &Kernel,
    _b: &mut Brep,
    _keep_edges: Option<&[EdgeId]>,
    _keep_vertices: Option<&[VertexId]>,
) -> Result<()> {
    // OMITTED WIP: SmBrep.cpp:25663-26158 (ordered edge/vertex cleanup and
    // temporary tolerance, orientation, geometry and removal bookkeeping).
    // UNPORTED(topology.brep_p1): SmBrep::GetEdges
    missing("SmBrep::GetEdges", "topology.brep_p1")
}

/// Removes the requested faces then their wire edges and shell vertices.
/// Even an empty list advances and releases one unlocked context traversal mark.
/// Editing-disabled refusal occurs before the mark is acquired.
pub fn remove_faces(
    k: &mut Kernel,
    b: &mut Brep,
    editing_enabled: bool,
    faces: &[FaceId],
) -> Result<()> {
    require_editing(editing_enabled)?;
    let mut mark = MarkLock::new(core::slice::from_mut(&mut k.context), 0, ALL_MARKS)?;
    let result = (|| {
        let mut all_edges = Vec::new();
        let mut all_vertices = Vec::new();
        for &face in faces {
            b.face(face)?;
            for edge in face_edges(b, face)? {
                if !all_edges.contains(&edge) {
                    all_edges.push(edge);
                }
            }
            for vertex in face_vertices(b, face)? {
                if !all_vertices.contains(&vertex) {
                    all_vertices.push(vertex);
                }
            }
            delete_face(b, face, false, true, None)?;
        }
        for edge in all_edges {
            if crate::topology::edge::is_wire(b, edge)? {
                delete_edge(b, edge)?;
            }
        }
        for vertex in all_vertices {
            if vertex::is_shell_vertex(b, vertex)? {
                delete_vertex(b, vertex)?;
            }
        }
        Ok(())
    })();
    mark.release(core::slice::from_mut(&mut k.context))?;
    result
}

fn face_edges(_b: &Brep, _face: FaceId) -> Result<Vec<EdgeId>> {
    // UNPORTED(topology.face_p4): SmFace::GetEdges
    missing("SmFace::GetEdges", "topology.face_p4")
}
fn face_vertices(_b: &Brep, _face: FaceId) -> Result<Vec<VertexId>> {
    // UNPORTED(topology.face_p4): SmFace::GetVertices
    missing("SmFace::GetVertices", "topology.face_p4")
}

/// Stitching and orientation options with the native defaults.
#[derive(Debug, Clone, PartialEq)]
pub struct StitchAndOrientOptions {
    /// The raw stitching tolerance; nonpositive values select twice the B-rep tolerance.
    pub stitch_tolerance: f64,
    /// Allows formation of spine edges.
    pub allow_spine_edges: bool,
    /// Optional faces to swap UV after import.
    pub faces_to_swap: Option<Vec<FaceId>>,
    /// Overrides the lamina-count-based decision when present, including `Some(false)`.
    pub just_unify: Option<bool>,
}

impl Default for StitchAndOrientOptions {
    fn default() -> Self {
        Self {
            stitch_tolerance: -0.1,
            allow_spine_edges: false,
            faces_to_swap: None,
            just_unify: None,
        }
    }
}

/// Computes the first stitching tolerance exactly as `StitchAndOrient` does.
/// Only a nonpositive requested tolerance is capped at `100 * SM_ZONE_TOL_3D`.
pub fn initial_stitch_tolerance(brep_tolerance: f64, requested: f64) -> f64 {
    let local = 2.0 * brep_tolerance;
    if requested <= 0.0 {
        if local > 100.0 * ZONE_TOL_3D {
            100.0 * ZONE_TOL_3D
        } else {
            local
        }
    } else {
        requested
    }
}

/// Stitches twice if needed, then orients surfaces. No success is manufactured
/// from the statuses ignored upstream when a required callee is missing.
pub fn stitch_and_orient(b: &mut Brep, options: &StitchAndOrientOptions) -> Result<()> {
    let mut lamina = 0;
    let mut tolerance = initial_stitch_tolerance(b.tolerance(), options.stitch_tolerance);
    crate::error::ignore_status(stitch_faces(
        b,
        tolerance,
        options.allow_spine_edges,
        &mut lamina,
    ))?;
    if lamina > 1 {
        tolerance = if options.stitch_tolerance <= 0.0 {
            tolerance * 20.0
        } else {
            tolerance * 10.0
        };
        crate::error::ignore_status(stitch_faces(
            b,
            tolerance,
            options.allow_spine_edges,
            &mut lamina,
        ))?;
    }
    crate::error::ignore_status(orient_trimmed_surfaces(
        b,
        true,
        options.just_unify.unwrap_or(lamina > 0),
        options.faces_to_swap.as_deref(),
    ))?;
    Ok(())
}

fn stitch_faces(
    _b: &mut Brep,
    _tolerance: f64,
    _allow_spine: bool,
    _lamina: &mut usize,
) -> Result<()> {
    // UNPORTED(topology.brep_p4): SmBrep::StitchFaces
    missing("SmBrep::StitchFaces", "topology.brep_p4")
}

fn orient_trimmed_surfaces(
    _b: &mut Brep,
    _normals_outward: bool,
    _just_unify: bool,
    _faces_to_swap: Option<&[FaceId]>,
) -> Result<()> {
    // UNPORTED(topology.brep_p3): SmBrep::OrientTrimmedSurfaces
    missing("SmBrep::OrientTrimmedSurfaces", "topology.brep_p3")
}

/// Output lists of `SqueezeFace`. All are reset even when editing is disabled.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct SqueezedFace {
    /// Deleted vertices, as stale identities.
    pub deleted_vertices: Vec<VertexId>,
    /// New edges.
    pub new_edges: Vec<EdgeId>,
    /// Deleted edges, as stale identities.
    pub deleted_edges: Vec<EdgeId>,
}

/// Deletes a degenerate face and stitches its former boundary topology.
pub fn squeeze_face(
    b: &mut Brep,
    editing_enabled: bool,
    face: FaceId,
    _stitch_tolerance: Option<f64>,
    output: &mut SqueezedFace,
) -> Result<()> {
    *output = SqueezedFace::default();
    require_editing(editing_enabled)?;
    b.face(face)?;
    // OMITTED WIP: SmBrep.cpp:27348-27469 (box lock and stitch configuration,
    // region choice, face deletion and all three cleanup passes/output updates).
    // The box lock precedes GetVertices in native code and must be restored too.
    // UNPORTED(topology.face_p4): SmFace::GetVertices
    missing("SmFace::GetVertices", "topology.face_p4")
}

/// Scales about the current origin then applies the placement. The identity exit
/// precedes notifications, geometry conversion and tolerance scaling, even with
/// `scale_tolerances = true`. Nonidentity transformation requires part 2's Notify
/// before any point or geometry is changed; those continuations are uncertified.
pub fn transform(
    _k: &Kernel,
    _b: &mut Brep,
    placement: &Axis2Placement,
    scale: Option<Vec3>,
    _scale_tolerances: bool,
) -> Result<()> {
    if placement.is_identity(crate::math::EFF_ZERO)
        && scale.is_none_or(|s| s.approx_eq(Vec3::new(1.0, 1.0, 1.0)))
    {
        return Ok(());
    }
    // OMITTED WIP: SmBrep.cpp:28799-28975 (NURB conversion, geometry edits,
    // trim flags, reflection correction, model-size and tolerance scaling).
    // UNPORTED(topology.brep_p2): SmBrep::Notify
    missing("SmBrep::Notify", "topology.brep_p2")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::Vertex;

    #[test]
    fn failure_outputs_follow_native_lifecycle() {
        let mut b = Brep::new();
        let v = b.add_vertex_record(Vertex::default()).unwrap();
        let mut edge = Some(EdgeId::from_index(3));
        assert!(
            combine_coincident_vertices(&mut Kernel::default(), &mut b, v, v, &mut edge).is_err()
        );
        assert!(edge.is_none());
        let mut out = SqueezedFace {
            deleted_vertices: vec![v],
            new_edges: vec![EdgeId::from_index(1)],
            deleted_edges: vec![EdgeId::from_index(2)],
        };
        assert!(squeeze_face(&mut b, false, FaceId::from_index(0), None, &mut out).is_err());
        assert_eq!(out, SqueezedFace::default());
    }

    #[test]
    fn identity_precedes_tolerance_update_and_notifications() {
        let mut b = Brep::new();
        let mut k = Kernel::default();
        let v = b
            .add_vertex_record(Vertex {
                point: Vec3::new(1., 2., 3.),
                tolerance: 42.,
                ..Vertex::default()
            })
            .unwrap();
        let original = b.vertex(v).unwrap().clone();
        transform(&k, &mut b, &Axis2Placement::IDENTITY, None, true).unwrap();
        assert_eq!(b.vertex(v).unwrap(), &original);
        assert_eq!(b.vertices().len(), 1);
        assert!(remove_faces(&mut k, &mut b, true, &[]).is_ok());
        assert!(remove_faces(&mut k, &mut b, false, &[]).is_err());
    }

    #[test]
    fn empty_face_removal_preserves_an_existing_mark_lock() {
        use crate::geom::context::MarkType;
        let mut k = Kernel::default();
        let mut b = Brep::new();
        let mut outer = MarkLock::new(core::slice::from_mut(&mut k.context), 0, ALL_MARKS).unwrap();
        let held = outer.mark_type();
        let marks = [
            MarkType::Mark,
            MarkType::Mark2,
            MarkType::Mark3,
            MarkType::MarkIo,
            MarkType::MarkAssert,
        ];
        let before = marks.map(|m| k.context.current_mark(m));
        remove_faces(&mut k, &mut b, true, &[]).unwrap();
        assert!(k.context.is_mark_locked(held));
        let advanced: Vec<_> = marks
            .into_iter()
            .zip(before)
            .filter(|(m, n)| k.context.current_mark(*m) != *n)
            .collect();
        assert_eq!(advanced.len(), 1);
        assert_ne!(advanced[0].0, held);
        assert_eq!(k.context.current_mark(advanced[0].0), advanced[0].1 + 2);
        assert!(!k.context.is_mark_locked(advanced[0].0));
        let before_refusal = marks.map(|m| k.context.current_mark(m));
        assert!(remove_faces(&mut k, &mut b, false, &[]).is_err());
        assert_eq!(marks.map(|m| k.context.current_mark(m)), before_refusal);
        outer
            .release(core::slice::from_mut(&mut k.context))
            .unwrap();
        assert!(!k.context.is_mark_locked(held));
    }
}
