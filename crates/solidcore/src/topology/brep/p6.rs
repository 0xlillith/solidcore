//! Edge copying, face topology construction and region-closure testing.
//!
//! Upstream: `source/SMLib/src/SmBrep.cpp`, `source/SMLib/inc/SmBrep.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`copy_edges`] | `SmBrep::CopyEdges` |
//! | [`make_face_topology`] | `SmBrep::MakeFaceTopology` |
//! | [`test_closure_split_region`] | `SmBrep::TestClosureSplitRegion(SmFaceuse*)` |
//! | [`test_closure_split_region_attempt`] | same, retaining native status and partial outputs |
//!
//! Deviations: pointers become arena ids; optional outputs become return values.
//! Edits preserve the upstream mutations preceding an error. Missing callees stop
//! at their first call: face construction at `SetTolerance` (`face_p3`),
//! Infinite-region reassignment uses the caller's shell and B-rep ray provider.
//! Closed-shell splits calculate their
//! real rough box after two `SmBrep::Notify` calls
//! ([`crate::topology::brep::p2::notify`] with the caller's [`NotifyEnv`], statuses
//! ignored; the new region gets its default `SmTopology` record when it is made, as
//! `new SmRegion`). Residual distribution uses the caller's shell ray context;
//! inner-shell classification reaches the traverser's real face-intersection
//! dependency. Empty shell/region deletion unlinks records while retaining their
//! arena slots, so no live id is renumbered. Native output ids remain published.
//! The shell implementation comes from PR #154. Surface owner identities are
//! supplied by the caller, independently of arena indices. A null trimmed copy
//! returns a null-pointer error instead of the upstream dereference. There is no recursion;
//! the traverser uses its bounded, iterative implementation.
//! The copy notifications (`SM_NO_COPY` on each source edge and its vertices) run on a
//! `topology.brep_p2` arena state seeded from both B-reps: the callers keep no attribute
//! state, and both B-reps' objects are named by their arena identities, which overlap
//! across the two B-reps (an overlapping pair propagates nothing, as for an object
//! notified about itself). `SmEdge::Notify` and `SmVertex::Notify` pass `SM_NO_COPY`
//! straight to `SmTopology::Notify`, which is called directly (the source B-rep is only
//! read). They notify through the destination kernel's context, the only one given.
//! An attempt exposes native status separately for callers that ignore it;
//! missing dependencies and Rust representation failures remain outer errors.

use crate::error::{EntityRef, Result, SolidError};
use crate::geom::attribute::ObjectKey;
use crate::geom::context::{MarkLock, ALL_MARKS};
use crate::geom::object::{NotifyOperation, ObjectRef};
use crate::intersect::point_class::PointClassification;
use crate::kernel::Kernel;
use crate::math::{Point3, Vec3, EFF_ZERO_SQRT};
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::surface::surface::Surface;
use crate::topology::brep::p2::{self as brep_p2, NotifyArg};
use crate::topology::edge::{self, NotifyEnv};
use crate::topology::edgeuse::p1 as eu;
use crate::topology::edgeuse::{p2 as eu_p2, p3 as eu_p3};
use crate::topology::shell::{self, ShellType};
use crate::topology::topology::Topology;
use crate::topology::topology_traverser::{self as tr, Entity, Marks};
use crate::topology::{
    Brep, EdgeId, EdgeuseId, EdgeuseOwner, Face, FaceId, Faceuse, FaceuseId, LoopuseStart, Region,
    RegionId, Shell, ShellId, VertexId,
};

const OP: &str = "topology::brep::p6";

fn required<T>(value: Option<T>) -> Result<T> {
    value.ok_or_else(|| SolidError::null_pointer(OP))
}

// The split's bare pointer dereferences have no native error status to ignore.
// Refuse missing arena links as representation errors, not synthetic SmStatus.
fn split_link<T>(value: Option<T>) -> Result<T> {
    value.ok_or_else(|| SolidError::invalid_input(OP, "unlinked topology in region split"))
}

/// Copy edges to independent wires in `destination`, in input order, returning the new
/// edges (`pNewEdges`). An empty list succeeds without inspecting the destination. Each
/// edge's trimmed NURBS copy is made into a wire edge by part 5's
/// `CreateWireEdgeFromCurve` over the copy's natural interval; then the source edge and
/// its start vertex (and its end vertex, when both edges have two) report the copy
/// (`Notify(SM_NO_COPY, pNew, pDestinationBrep, this)`, statuses ignored; see the
/// module's Deviations).
///
/// # Errors
///
/// Those of `CreateTrimmedNURBSCurve` and `CreateWireEdgeFromCurve` (`SER`); a
/// null-pointer error for a null curve copy or a missing vertex, which upstream
/// dereferences.
pub fn copy_edges(
    destination_kernel: &mut Kernel,
    source: &Brep,
    edges: &[EdgeId],
    destination: &mut Brep,
) -> Result<Vec<EdgeId>> {
    let mut result = Vec::new();
    for &edge in edges {
        let curve = required(edge::create_trimmed_nurbs_curve(
            destination_kernel,
            source,
            edge,
        )?)?;
        let interval = curve.natural_interval();
        let new_edge = crate::topology::brep::p5::create_wire_edge_from_curve(
            destination_kernel,
            destination,
            curve,
            interval,
        )?;

        // Notify Copy Edge - propogates attributes and reports copy
        let mut state = brep_p2::ArenaNotifyState::new(source);
        state.seed(destination);
        let data = [brep_p2::brep_ref(destination), brep_p2::brep_ref(source)];
        let mut notify_copy = |old: EntityRef, new: EntityRef| {
            use crate::topology::loopuse::notify_ref;
            crate::error::ignore_status(crate::topology::topology::notify(
                &mut state.records,
                &mut state.attributes,
                Some(&mut destination_kernel.context),
                NotifyOperation::Copy,
                notify_ref(old),
                Some(notify_ref(new)),
                Some(data[0]),
                Some(data[1]),
            ))
        };
        notify_copy(EntityRef::Edge(edge), EntityRef::Edge(new_edge))?;

        // Notify Copy Vertices - propogates attributes and reports copy
        let new_start = required(edge::start_vertex(destination, new_edge)?)?;
        let new_end = required(edge::other_vertex(destination, new_edge, Some(new_start))?)?;
        let old_start = required(edge::start_vertex(source, edge)?)?;
        let old_end = required(edge::other_vertex(source, edge, Some(old_start))?)?;
        notify_copy(EntityRef::Vertex(old_start), EntityRef::Vertex(new_start))?;
        if old_start != old_end && new_start != new_end {
            notify_copy(EntityRef::Vertex(old_end), EntityRef::Vertex(new_end))?;
        }

        // set optional output
        result.push(new_edge);
    }
    Ok(result)
}

/// Topology created by `MakeFaceTopology` when construction succeeds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaceTopology {
    /// Newly created face.
    pub face: FaceId,
    /// New shell, including the shell created for an edgeless face.
    pub shell: Option<ShellId>,
    /// New region when the face splits an existing region.
    pub region: Option<RegionId>,
}

/// Input loops, isolated vertices, and the upstream face-construction flags.
#[derive(Debug, Clone, Default)]
pub struct FaceTopologyOptions<'a> {
    /// Explicit target region; `None` selects the B-rep's infinite region.
    pub region: Option<RegionId>,
    /// Number of input edge uses in each loop, outer loop first.
    pub edge_loops: &'a [usize],
    /// Concatenated, ordered edge uses, one closed sequence per loop.
    pub edgeuses: &'a [EdgeuseId],
    /// Isolated loop vertices.
    pub loop_vertices: &'a [VertexId],
    /// Verify radial sectors; upstream default is false.
    pub verify_radial_sectors: bool,
    /// Calculate UV trims; upstream default is false.
    pub calculate_trims: bool,
    /// Optional UV curves (including null entries) for the trim stage.
    pub uv_curves: Option<&'a [Option<crate::nurbs::b_spline_curve::BSplineCurve>]>,
}

/// Begin constructing a face from ordered edge-use loops and isolated vertices.
/// `face_key` is the caller's identity for the new face. A missing region or
/// surface fails before allocation. Malformed loop counts return `InvalidInput`
/// instead of reading past the input. The active `SM_USE_OLDTOL` call to
/// `SmFace::SetTolerance` currently defers every otherwise valid construction.
/// On error, allocated records remain in `brep`, as upstream retains its partial
/// edits; callers must not treat an error as an atomic rollback.
pub fn make_face_topology(
    brep: &mut Brep,
    surface: Option<Surface>,
    _orientation: OrientType,
    face_key: ObjectKey,
    options: &FaceTopologyOptions<'_>,
) -> Result<FaceTopology> {
    let region = required(options.region.or(brep.infinite_region()))?;
    brep.region(region)?;
    let mut surface = required(surface)?;
    let shell = if let Some(&first) = options.edgeuses.first() {
        required(eu::shell(brep, first)?)?
    } else {
        let shell = brep.add_shell_record(Shell {
            shell_type: ShellType::Faceuse,
            region: Some(region),
            ..Shell::default()
        })?;
        brep.region_mut(region)?.shells.push(shell);
        shell
    };

    let mut used = 0usize;
    for &count in options.edge_loops {
        // Upstream skips zero-sized loops after its release assertion.
        if count == 0 {
            continue;
        }
        let end = used
            .checked_add(count)
            .ok_or_else(|| SolidError::invalid_input(OP, "edge loop count overflow"))?;
        let loop_edges = options
            .edgeuses
            .get(used..end)
            .ok_or_else(|| SolidError::invalid_input(OP, "edge loop exceeds input edge uses"))?;
        let mut last = *required(loop_edges.last())?;
        for &current in loop_edges {
            let mate = required(eu::mate(brep, last)?)?;
            let current_vertex = required(brep.edgeuse(current)?.vertexuse)?;
            let last_vertex = required(brep.edgeuse(mate)?.vertexuse)?;
            // SE(INVALID_INPUT), not SER: inconsistent ordering only reports a
            // diagnostic upstream. Do not turn it into an early return here.
            let _connected =
                brep.vertexuse(current_vertex)?.vertex == brep.vertexuse(last_vertex)?.vertex;
            last = current;
        }
        used = end;
    }

    let face = brep.add_face_record(Face::new(shell))?;
    let first = brep.add_faceuse_record(Faceuse {
        face: Some(face),
        ..Faceuse::default()
    })?;
    let second = brep.add_faceuse_record(Faceuse {
        face: Some(face),
        mate: Some(first),
        ..Faceuse::default()
    })?;
    brep.faceuse_mut(first)?.mate = Some(second);
    surface.set_owner(Some(ObjectRef {
        type_code: crate::topology::face::p2::TYPE,
        key: face_key,
    }));
    let domain = surface.natural_uv_domain();
    let record = brep.face_mut(face)?;
    record.surface = Some(surface);
    record.uv_domain = domain;
    // This call precedes primary-faceuse assignment, list insertion, loop
    // allocation, radial sectors, trims and TestClosureSplitRegion upstream.
    // UNPORTED(topology.face_p3): SmFace::SetTolerance -- justified: findings.md topology.face_p3 "Stubs left for this unit"
    // (`topology::face::p3::set_tolerance` is ported; this is the stop point of the
    // unported remainder of MakeFaceTopology, owned by topology.brep_p6)
    Err(SolidError::unported(
        "SmFace::SetTolerance",
        "topology.face_p3",
    ))
}

/// Caller-owned shell and B-rep queries used while splitting a region.
///
/// The closure implementation receives `Some(shell)` for a shell ray and `None`
/// for a complete B-rep ray; callers must implement both selections.
///
/// The provider operates on the actual, mutating B-rep and preserves the ray's
/// partial output state even when its native status is an error. A harness may
/// use `check_geometry` to reject an input representation only when geometry is
/// read; native topology-only fast paths do not call that guard.
pub trait RegionSplitQueries {
    /// Verify that the geometry kinds about to be read are represented faithfully.
    ///
    /// # Errors
    /// An explicit representation limitation; the default provider has none.
    fn check_geometry(
        &mut self,
        _brep: &Brep,
        _read_curves: bool,
        _read_surfaces: bool,
    ) -> Result<()> {
        Ok(())
    }

    /// Verify the representation of one edge curve immediately before reading it.
    ///
    /// The default retains the provider's existing curve check. Providers with
    /// metadata for individual edges can avoid rejecting unvisited geometry.
    /// Topological early returns do not call this check.
    ///
    /// # Errors
    /// An explicit representation limitation from [`Self::check_geometry`].
    fn check_edge_curve_geometry(&mut self, brep: &Brep, _edge: EdgeId) -> Result<()> {
        self.check_geometry(brep, true, false)
    }

    /// Fire against the selected shell using the caller's real query context.
    fn ray_fire(
        &mut self,
        k: &mut Kernel,
        brep: &mut Brep,
        shell: ShellId,
        input: shell::ShellRayInput<EntityRef>,
    ) -> shell::ShellRayAttempt<EntityRef>;

    /// Fire against the complete B-rep using its actual identity and live cache.
    /// Preserve native status and all partial output fields. A supplied shell
    /// can fall back to this operation while finding the infinite region.
    fn brep_ray_fire(
        &mut self,
        k: &mut Kernel,
        brep: &mut Brep,
        input: shell::ShellRayInput<EntityRef>,
    ) -> shell::ShellRayAttempt<EntityRef>;
}

impl<F> RegionSplitQueries for F
where
    F: FnMut(
        &mut Kernel,
        &mut Brep,
        Option<ShellId>,
        shell::ShellRayInput<EntityRef>,
    ) -> shell::ShellRayAttempt<EntityRef>,
{
    fn ray_fire(
        &mut self,
        k: &mut Kernel,
        brep: &mut Brep,
        shell: ShellId,
        input: shell::ShellRayInput<EntityRef>,
    ) -> shell::ShellRayAttempt<EntityRef> {
        self(k, brep, Some(shell), input)
    }

    fn brep_ray_fire(
        &mut self,
        k: &mut Kernel,
        brep: &mut Brep,
        input: shell::ShellRayInput<EntityRef>,
    ) -> shell::ShellRayAttempt<EntityRef> {
        self(k, brep, None, input)
    }
}

/// Split a region when the seed face use bounds a closed set of neighbors.
/// Returns `(new_shell, new_region)`. Already separated mates, lamina faces and
/// non-closed traversals succeed without creating either object. Closed sets move
/// to a new solid region/shell, both notified through `env`, before calculating
/// the rough B-rep box and distributing residual faces, wires and nested shells.
/// The caller's marks and kernel context persist across calls. Every acquired
/// mark lock is released on success and on error.
pub fn test_closure_split_region(
    k: &mut Kernel,
    brep: &mut Brep,
    marks: &mut Marks,
    env: &mut NotifyEnv<'_>,
    do_infinite_region_testing: bool,
    start: FaceuseId,
    queries: &mut dyn RegionSplitQueries,
    identity: crate::intersect::curve_class::p7::FaceQueryIdentity<'_>,
) -> Result<(Option<ShellId>, Option<RegionId>)> {
    let attempt = test_closure_split_region_attempt(
        k,
        brep,
        marks,
        env,
        do_infinite_region_testing,
        start,
        queries,
        identity,
    )?;
    attempt.status?;
    Ok((attempt.shell, attempt.region))
}

/// Native status and outputs written by `TestClosureSplitRegion` before it exits.
#[derive(Debug)]
pub struct RegionSplitAttempt {
    /// The native `SmStatus`, represented as success or [`SolidError::Upstream`].
    /// A caller whose C++ counterpart ignores the status may ignore this field.
    pub status: Result<()>,
    /// New shell once the upstream output assignment is reached, even on error.
    pub shell: Option<ShellId>,
    /// New region once the upstream output assignment is reached, even on error.
    pub region: Option<RegionId>,
}

fn finish_split_attempt(
    mut attempt: RegionSplitAttempt,
    status: Result<()>,
) -> Result<RegionSplitAttempt> {
    match status {
        Ok(()) => Ok(attempt),
        Err(error @ SolidError::Upstream { .. }) => {
            attempt.status = Err(error);
            Ok(attempt)
        }
        // These failures do not correspond to a native status. In particular,
        // an unavailable callee or missing typed resolver must never become a
        // successful empty split just because native callers ignore SmStatus.
        Err(error) => Err(error),
    }
}

/// Test region closure, retaining native error status and partial output writes.
/// The outer result propagates unported dependencies and Rust representation
/// failures. Native `SmStatus` errors are in [`RegionSplitAttempt::status`], so
/// unchecked C++ callers can continue with the same partially written outputs.
pub fn test_closure_split_region_attempt(
    k: &mut Kernel,
    brep: &mut Brep,
    marks: &mut Marks,
    env: &mut NotifyEnv<'_>,
    do_infinite_region_testing: bool,
    start: FaceuseId,
    queries: &mut dyn RegionSplitQueries,
    identity: crate::intersect::curve_class::p7::FaceQueryIdentity<'_>,
) -> Result<RegionSplitAttempt> {
    let mut attempt = RegionSplitAttempt {
        status: Ok(()),
        shell: None,
        region: None,
    };
    let status = split_region_into(
        k,
        brep,
        marks,
        env,
        start,
        &mut SplitRun {
            attempt: &mut attempt,
            queries,
            do_infinite_region_testing,
            identity,
        },
    );
    finish_split_attempt(attempt, status)
}

struct SplitRun<'a> {
    identity: crate::intersect::curve_class::p7::FaceQueryIdentity<'a>,
    attempt: &'a mut RegionSplitAttempt,
    queries: &'a mut dyn RegionSplitQueries,
    do_infinite_region_testing: bool,
}

fn split_region_into(
    k: &mut Kernel,
    brep: &mut Brep,
    marks: &mut Marks,
    env: &mut NotifyEnv<'_>,
    start: FaceuseId,
    run: &mut SplitRun<'_>,
) -> Result<()> {
    let seed = brep.faceuse(start)?;
    if let Some(mate) = seed.mate {
        if brep.faceuse(mate)?.shell != seed.shell {
            return Ok(());
        }
    }
    if tr::has_lamina_edge(brep, split_link(seed.face)?)? {
        return Ok(());
    }
    let mut lock = MarkLock::new(std::slice::from_mut(&mut k.context), 0, ALL_MARKS)?;
    let result = split_closed(k, brep, marks, env, start, &lock, run);
    lock.release(std::slice::from_mut(&mut k.context))?;
    result
}

#[allow(clippy::too_many_lines)]
fn split_closed(
    k: &mut Kernel,
    brep: &mut Brep,
    marks: &mut Marks,
    env: &mut NotifyEnv<'_>,
    start: FaceuseId,
    lock: &MarkLock,
    run: &mut SplitRun<'_>,
) -> Result<()> {
    if !tr::closure_traversal(k, brep, marks, start, lock, None)? {
        return Ok(());
    }
    let new_shell = brep.add_shell_record(Shell {
        shell_type: ShellType::Faceuse,
        ..Shell::default()
    })?;
    let new_region = brep.add_region_record(Region {
        in_brep: true,
        shells: vec![new_shell],
        is_void: 0,
    })?;
    brep.shell_mut(new_shell)?.region = Some(new_region);
    // the new region's SmTopology members (default-constructed)
    let region_key = (env.key_of)(EntityRef::Region(new_region)).key;
    env.records.insert(region_key, Topology::default());

    let old_shell = split_link(brep.faceuse(start)?.shell)?;
    let old_region = split_link(brep.shell(old_shell)?.region)?;
    // Keep these native list queries in their upstream order; arena iteration
    // would include detached objects and change the traversal.
    let wires = shell::wire_edges(k, brep, old_shell)?;
    let faceuses = brep.shell(old_shell)?.get_faceuses()?;
    let shells = brep.region(old_region)?.shells.clone();
    let mut replaced = Vec::new();
    tr::faceuse_shell_replacement(k, brep, marks, start, new_shell, &mut replaced, lock)?;
    let mut other = Vec::new();
    tr::faceuse_traversal(
        k,
        brep,
        marks,
        split_link(brep.faceuse(start)?.mate)?,
        lock.mark_type(),
        Some(&mut other),
    )?;
    let mut all_replaced = replaced.clone();
    all_replaced.extend(other);
    // Upstream publishes these pointers only after both traversals. They stay
    // available to unchecked callers when a later native operation fails.
    run.attempt.shell = Some(new_shell);
    run.attempt.region = Some(new_region);
    // notify this Brep and the OldRegion for attribute propagation (SmBrep.cpp:27954-5;
    // both statuses ignored)
    let this = Some(NotifyArg::Object(brep_p2::brep_ref(brep)));
    ignore_status(brep_p2::notify(
        brep,
        env,
        None,
        NotifyOperation::PreEdit,
        [this, None, None],
    ))?;
    let region = |r| Some(NotifyArg::Entity(EntityRef::Region(r)));
    ignore_status(brep_p2::notify(
        brep,
        env,
        None,
        NotifyOperation::SplitInBrep,
        [region(old_region), region(old_region), region(new_region)],
    ))?;
    // SmBrep::TestClosureSplitRegion: SER(CalculateBoundingBox(sBrepBBox, TRUE, FALSE));
    // the omitted expansion argument retains its TRUE header default.
    if brep.locked_bounding_box.is_none() {
        let vertices = super::p1::vertices(brep, None)?;
        let edges = super::p1::edges(brep, None)?;
        run.queries
            .check_geometry(brep, !(vertices.len() > 500 && edges.len() > 100), false)?;
    }
    let mut brep_bbox = crate::query::brep::calculate_bounding_box(
        k,
        brep,
        crate::query::brep::BoundingBoxOptions {
            approximate: true,
            tight: false,
            expand_by_zone_tolerance: true,
        },
    )?;
    let radius = brep_bbox.sphere_bound().radius;
    reassign_split_infinite_region(
        k,
        brep,
        old_region,
        new_shell,
        &mut brep_bbox,
        run.do_infinite_region_testing,
        run.queries,
    )?;

    let faceuses_remain = faceuses.len() != all_replaced.len();
    if faceuses_remain || !wires.is_empty() || shells.len() > 1 {
        // This status is not checked natively. IsInnerShell first writes UNSURE
        // and leaves that nonzero value on every native error, so the caller's
        // boolean branch still selects the new shell. Surface sampling may also
        // read boundary curves; its harness guard does not replace the operation.
        run.queries.check_geometry(brep, true, true)?;
        let new_is_inner =
            ignored_native(tr::is_inner_shell(k, brep, &replaced, run.identity))?.unwrap_or(true);
        let (inner_shell, other_shell) = if new_is_inner {
            (new_shell, old_shell)
        } else {
            (old_shell, new_shell)
        };
        let mut distribution = Distribution {
            identity: run.identity,
            old_shell,
            new_shell,
            old_region,
            new_region,
            inner_shell,
            other_shell,
            radius,
            queries: &mut *run.queries,
        };
        if faceuses_remain {
            distribution.faces(k, brep, marks, lock, &faceuses, &mut all_replaced)?;
        }
        if !wires.is_empty() {
            distribution.wires(k, brep, marks, lock, &wires)?;
        }
        if shells.len() > 1 {
            distribution.shells(k, brep, &shells)?;
        }
    }
    unlink_empty_split_shell(brep, old_shell)?;
    unlink_empty_split_shell(brep, new_shell)?;
    Ok(())
}

struct Distribution<'a> {
    identity: crate::intersect::curve_class::p7::FaceQueryIdentity<'a>,
    old_shell: ShellId,
    new_shell: ShellId,
    old_region: RegionId,
    new_region: RegionId,
    inner_shell: ShellId,
    other_shell: ShellId,
    radius: f64,
    queries: &'a mut dyn RegionSplitQueries,
}

// Scope-equivalent to SmTemporaryChangeValue<SmBoolean>(m_bEditingEnabled,FALSE).
// Capturing the result before restoring also handles every fallible exit.
fn without_editing<T>(brep: &mut Brep, operation: impl FnOnce(&mut Brep) -> T) -> T {
    let previous = brep.editing_enabled();
    brep.set_editing_enabled(false);
    let result = operation(brep);
    brep.set_editing_enabled(previous);
    result
}

impl Distribution<'_> {
    fn recipient(&self, region: Option<RegionId>) -> ShellId {
        match region {
            None => self.other_shell,
            Some(region) if region == self.new_region => self.new_shell,
            Some(_) => self.old_shell,
        }
    }

    fn face_region(
        &mut self,
        k: &mut Kernel,
        brep: &mut Brep,
        faceuse: FaceuseId,
    ) -> Result<Option<RegionId>> {
        let edgeuses = eu::faceuse_edgeuses(brep, faceuse, false, None)?;
        let edgeuse = *edgeuses.first().ok_or_else(|| SolidError::sm_err(OP))?;
        let edge = split_link(eu::edge(brep, edgeuse)?)?;
        let tolerance = edge::tolerance(brep, edge)?;
        let mut step = 20.0 * tolerance;
        let mut state = shell::ShellRayState {
            region: None,
            classification: PointClassification::new(tolerance),
        };
        for iteration in 0..100 {
            self.queries.check_geometry(brep, true, true)?;
            let parameter = edge::interval(brep, edge)?.evaluate(0.4567);
            let stepped = ignored_native(eu_p2::evaluate_binormal_step_off(
                k,
                brep,
                edgeuse,
                parameter,
                step,
                &eu_p2::StepOffOptions::default(),
            ))?;
            let last_try = stepped.is_none();
            let point = if let Some(stepped) = stepped {
                stepped.point
            } else {
                if iteration > 0 {
                    return Err(SolidError::sm_err(OP));
                }
                let parameter = edge::interval(brep, edge)?.evaluate(0.4567);
                eu_p3::evaluate_binormal(
                    k,
                    brep,
                    edgeuse,
                    parameter,
                    false,
                    &eu_p3::BinormalOptions::default(),
                )?
                .point
            };
            let attempt = without_editing(brep, |brep| {
                self.queries.ray_fire(
                    k,
                    brep,
                    self.inner_shell,
                    shell::ShellRayInput {
                        start: point,
                        direction: Vec3::new(1.0, 0.0, 0.0),
                        zone_tolerance: tolerance,
                        sphere_bound: Some(self.radius),
                        state,
                    },
                )
            });
            // PC and region are in/out values across native retries, including
            // writes made by a failed ray. Non-native refusals must still escape.
            state = attempt.state;
            if ignored_native(attempt.status)?.is_some() {
                break;
            }
            if last_try {
                return Err(SolidError::sm_err(OP));
            }
            step *= 2.0;
        }
        // Native release code uses the last region even after all 100 failures.
        Ok(state.region)
    }

    fn faces(
        &mut self,
        k: &mut Kernel,
        brep: &mut Brep,
        marks: &mut Marks,
        lock: &MarkLock,
        faceuses: &[FaceuseId],
        all_replaced: &mut Vec<FaceuseId>,
    ) -> Result<()> {
        let mut index = 0;
        let mut more = Vec::new();
        while let Some(&faceuse) = faceuses.get(index) {
            if !marks.is_marked(k, Entity::Faceuse(faceuse), lock.mark_type()) {
                let region = self.face_region(k, brep, faceuse)?;
                tr::faceuse_shell_replacement(
                    k,
                    brep,
                    marks,
                    faceuse,
                    self.recipient(region),
                    &mut more,
                    lock,
                )?;
                all_replaced.extend_from_slice(&more);
                // SmBrep.cpp reuses the outer `ii` in this remark loop. Retain
                // that traversal quirk instead of silently fixing native order.
                index = 0;
                for &replaced in all_replaced.iter() {
                    marks.mark(k, Entity::Faceuse(replaced), lock.mark_type());
                    index += 1;
                }
            }
            index += 1;
        }
        Ok(())
    }

    fn wires(
        &mut self,
        k: &mut Kernel,
        brep: &mut Brep,
        marks: &mut Marks,
        lock: &MarkLock,
        wires: &[EdgeId],
    ) -> Result<()> {
        let mut current = wires.first().copied();
        let tolerance = edge::tolerance(brep, split_link(current)?)?;
        let mut pending: Vec<_> = wires.iter().copied().map(Some).collect();
        let mut connected = Vec::new();
        while let Some(wire) = current {
            lock.new_mark(std::slice::from_mut(&mut k.context))?;
            tr::collect_wire_edges(
                k,
                brep,
                marks,
                Some(wire),
                lock.mark_type(),
                &mut connected,
                None,
            )?;
            // A valid wire is marked by CollectWireEdges. Otherwise the native
            // rescan has no progress and repeats forever on malformed topology.
            if !marks.is_marked(k, Entity::Edge(wire), lock.mark_type()) {
                return Err(SolidError::invalid_input(
                    OP,
                    "wire traversal made no progress",
                ));
            }
            let primary = split_link(edge::primary_edgeuse(brep, wire)?)?;
            self.queries.check_geometry(brep, true, false)?;
            let point =
                eu_p3::normalized_evaluate(k, brep, primary, 0.4567, false, false, true)?.point;
            current = without_editing(brep, |brep| {
                let attempt = self.queries.ray_fire(
                    k,
                    brep,
                    self.inner_shell,
                    shell::ShellRayInput {
                        start: point,
                        direction: Vec3::new(1.0, 0.0, 0.0),
                        zone_tolerance: tolerance,
                        sphere_bound: Some(self.radius),
                        state: shell::ShellRayState {
                            region: None,
                            classification: PointClassification::new(tolerance),
                        },
                    },
                );
                attempt.status?;
                let recipient = self.recipient(attempt.state.region);
                for &wire in &connected {
                    let first = split_link(edge::primary_edgeuse(brep, wire)?)?;
                    let second = split_link(eu::mate(brep, first)?)?;
                    brep.edgeuse_mut(first)?.owner = EdgeuseOwner::Shell(Some(recipient));
                    brep.edgeuse_mut(second)?.owner = EdgeuseOwner::Shell(Some(recipient));
                }
                let mut next = None;
                for entry in &mut pending {
                    if let Some(wire) = *entry {
                        if marks.is_marked(k, Entity::Edge(wire), lock.mark_type()) {
                            *entry = None;
                        } else {
                            // Native deliberately selects the LAST unmarked wire.
                            next = Some(wire);
                        }
                    }
                }
                Ok::<_, SolidError>(next)
            })?;
        }
        Ok(())
    }

    fn shell_point(&mut self, k: &mut Kernel, brep: &mut Brep, id: ShellId) -> Result<Point3> {
        let shell = brep.shell(id)?;
        if shell.is_faceuse_shell() {
            let faceuse = required(shell.first_faceuse()?)?;
            let edgeuses = eu::faceuse_edgeuses(brep, faceuse, false, None)?;
            if let Some(&edgeuse) = edgeuses.first() {
                let edge = split_link(eu::edge(brep, edgeuse)?)?;
                let step = 20.0 * edge::tolerance(brep, edge)?;
                let parameter = edge::interval(brep, edge)?.evaluate(0.4567);
                self.queries.check_geometry(brep, true, true)?;
                return Ok(eu_p2::evaluate_binormal_step_off(
                    k,
                    brep,
                    edgeuse,
                    parameter,
                    step,
                    &eu_p2::StepOffOptions::default(),
                )?
                .point);
            }
            let loopuse = *brep
                .faceuse(faceuse)?
                .loopuses
                .first()
                .ok_or_else(|| SolidError::sm_err(OP))?;
            let LoopuseStart::Vertexuse(vertexuse) = brep.loopuse(loopuse)?.start else {
                return Err(SolidError::sm_err(OP));
            };
            let vertex = split_link(brep.vertexuse(split_link(vertexuse)?)?.vertex)?;
            Ok(brep.vertex(vertex)?.point)
        } else if shell.is_wire_shell() {
            let edge = required(shell.wire_edge(brep)?)?;
            let edgeuse = split_link(edge::primary_edgeuse(brep, edge)?)?;
            self.queries.check_geometry(brep, true, false)?;
            Ok(eu_p3::normalized_evaluate(k, brep, edgeuse, 0.4567, false, false, true)?.point)
        } else if shell.is_vertex_shell() {
            let vertex = required(shell.vertex(brep)?)?;
            Ok(brep.vertex(vertex)?.point)
        } else {
            Ok(Point3::UNDEF)
        }
    }

    fn shells(&mut self, k: &mut Kernel, brep: &mut Brep, shells: &[ShellId]) -> Result<()> {
        let mut moved = false;
        for &shell in shells {
            if shell == self.old_shell || shell == self.new_shell {
                continue;
            }
            let point = self.shell_point(k, brep, shell)?;
            let tolerance = brep.tolerance;
            let step = EFF_ZERO_SQRT * (1.0 + point.max_dimension());
            let direction = Vec3::new(1.0, 1.0, 1.0).unitize()?;
            let point = point + direction * step;
            without_editing(brep, |brep| {
                let attempt = self.queries.ray_fire(
                    k,
                    brep,
                    self.inner_shell,
                    shell::ShellRayInput {
                        start: point,
                        direction,
                        zone_tolerance: tolerance,
                        sphere_bound: Some(self.radius),
                        state: shell::ShellRayState {
                            region: None,
                            classification: PointClassification::new(tolerance),
                        },
                    },
                );
                attempt.status?;
                let found = match attempt.state.region {
                    Some(region) => region,
                    None => split_link(brep.shell(self.other_shell)?.region)?,
                };
                let old = split_link(brep.shell(shell)?.region)?;
                if old != found {
                    remove_shell_ignored(brep, old, shell)?;
                    brep.region_mut(found)?.shells.push(shell);
                    brep.shell_mut(shell)?.region = Some(found);
                    moved = true;
                }
                Ok::<_, SolidError>(())
            })?;
        }
        if moved {
            ignore_status(shell::order_region_shells(
                k,
                brep,
                self.old_region,
                self.identity,
            ))?;
        }
        Ok(())
    }
}

// Keep the native infinite-region test between the rough bounding box and
// residual distribution. The same live provider can select shell or B-rep rays.
fn reassign_split_infinite_region(
    k: &mut Kernel,
    brep: &mut Brep,
    old_region: RegionId,
    new_shell: ShellId,
    bbox: &mut crate::geom::extent3d::Extent3d,
    enabled: bool,
    queries: &mut dyn RegionSplitQueries,
) -> Result<()> {
    if enabled && brep.infinite_region() == Some(old_region) {
        super::p5::find_and_set_infinite_region_with_queries(
            k,
            brep,
            Some(new_shell),
            2, // SM_VS_INFINITEREGION_SWAP
            Some(bbox),
            queries,
        )?;
    }
    Ok(())
}

fn remove_shell_ignored(brep: &mut Brep, region: RegionId, shell: ShellId) -> Result<()> {
    // This existing helper represents native Remove's SM_ERR as InvalidInput.
    // Its status is explicitly unchecked here, as in brep::p2's same operation.
    match shell::region_remove_shell(brep, region, shell) {
        Ok(()) | Err(SolidError::InvalidInput { .. }) => Ok(()),
        Err(error) => Err(error),
    }
}

fn unlink_empty_split_shell(brep: &mut Brep, shell: ShellId) -> Result<()> {
    let region = split_link(brep.shell(shell)?.region)?;
    if brep.shell(shell)?.get_faceuses()?.is_empty() {
        remove_shell_ignored(brep, region, shell)?;
        // Dense arenas retain disconnected records; deleting a native object
        // must not renumber the handles of the surviving topology.
    }
    if brep.region(region)?.shells.is_empty() {
        ignore_status(brep.remove_region_from_list(region))?;
        let head = brep.region_list_head();
        brep.set_infinite_region_record(head);
    }
    Ok(())
}

/// A status upstream ignores: only [`SolidError::Unported`] and the port's
/// representation failures ([`SolidError::InvalidInput`], [`SolidError::NotFound`])
/// propagate.
fn ignore_status(r: Result<()>) -> Result<()> {
    ignored_native(r).map(|_| ())
}

// Preserve the native ordinary-error fallback without turning a missing arena
// link, typed resolver or unported callee into a successful classification.
fn ignored_native<T>(result: Result<T>) -> Result<Option<T>> {
    match result {
        Err(
            e @ (SolidError::Unported { .. }
            | SolidError::InvalidInput { .. }
            | SolidError::NotFound { .. }),
        ) => Err(e),
        Ok(value) => Ok(Some(value)),
        Err(_) => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::context::MarkType;
    use crate::math::Vec3;
    use crate::surface::plane::Plane;
    use crate::topology::shell::ShellMembers;

    fn no_ray(
        _: &mut Kernel,
        _: &mut Brep,
        _: Option<ShellId>,
        _: shell::ShellRayInput<EntityRef>,
    ) -> shell::ShellRayAttempt<EntityRef> {
        panic!("this topology-only split must not cast a ray")
    }

    fn closed_pair() -> (Brep, RegionId, ShellId, FaceuseId, FaceuseId) {
        let mut b = Brep::new();
        let region = b
            .add_region_record(Region {
                in_brep: true,
                is_void: 1,
                ..Region::default()
            })
            .unwrap();
        let shell = b
            .add_shell_record(Shell {
                region: Some(region),
                shell_type: ShellType::Faceuse,
                ..Shell::default()
            })
            .unwrap();
        b.region_mut(region).unwrap().shells.push(shell);
        let face = b.add_face_record(Face::new(shell)).unwrap();
        let first = b
            .add_faceuse_record(Faceuse {
                shell: Some(shell),
                face: Some(face),
                ..Faceuse::default()
            })
            .unwrap();
        let second = b
            .add_faceuse_record(Faceuse {
                shell: Some(shell),
                face: Some(face),
                mate: Some(first),
                ..Faceuse::default()
            })
            .unwrap();
        b.faceuse_mut(first).unwrap().mate = Some(second);
        b.face_mut(face).unwrap().faceuse = Some(first);
        b.shell_mut(shell).unwrap().members = ShellMembers::Faceuses(vec![first, second]);
        b.set_infinite_region_record(Some(region));
        b.set_editing_enabled(true);
        (b, region, shell, first, second)
    }

    #[test]
    fn split_releases_its_lock_and_preserves_partial_outputs_on_native_infinite_error() {
        let mut k = Kernel::default();
        let (mut b, region, shell, first, second) = closed_pair();
        let mut outer = MarkLock::new(
            std::slice::from_mut(&mut k.context),
            0,
            MarkType::Mark.code(),
        )
        .unwrap();
        let mut state = crate::topology::brep::p2::RegionNotifyState::new(&b);
        let result = test_closure_split_region_attempt(
            &mut k,
            &mut b,
            &mut Marks::default(),
            &mut state.env(),
            true,
            first,
            &mut no_ray,
            crate::intersect::curve_class::p7::FaceQueryIdentity {
                keys: &crate::topology::brep_cache::arena_identity,
                owner_of: &|_| None,
            },
        );
        let attempt = result.unwrap();
        assert_eq!(
            attempt.status,
            Err(SolidError::sm_err("topology::brep::p9"))
        );
        assert_eq!(attempt.shell, Some(ShellId::from_index(1)));
        assert_eq!(attempt.region, Some(RegionId::from_index(1)));
        // the new region got its (default) SmTopology record
        assert_eq!(state.records.len(), 2);
        assert!(k.context.is_mark_locked(MarkType::Mark));
        assert!(!k.context.is_mark_locked(MarkType::Mark2));
        assert_eq!(b.regions().len(), 2);
        assert_eq!(b.region(RegionId::from_index(1)).unwrap().is_void(), 0);
        assert_eq!(b.region(region).unwrap().is_void(), 1);
        assert_eq!(
            b.faceuse(first).unwrap().shell,
            Some(ShellId::from_index(1))
        );
        assert_eq!(b.faceuse(second).unwrap().shell, Some(shell));
        outer.release(std::slice::from_mut(&mut k.context)).unwrap();
    }

    #[test]
    fn closed_pair_completes_without_infinite_testing_or_residual_rays() {
        let mut k = Kernel::default();
        let (mut b, region, old_shell, first, second) = closed_pair();
        let mut notifications = crate::topology::brep::p2::RegionNotifyState::new(&b);
        let mut marks = Marks::default();
        let attempt = test_closure_split_region_attempt(
            &mut k,
            &mut b,
            &mut marks,
            &mut notifications.env(),
            false,
            first,
            &mut no_ray,
            crate::intersect::curve_class::p7::FaceQueryIdentity {
                keys: &crate::topology::brep_cache::arena_identity,
                owner_of: &|_| None,
            },
        )
        .unwrap();
        attempt.status.unwrap();
        let new_shell = attempt.shell.unwrap();
        let new_region = attempt.region.unwrap();
        assert_eq!(
            b.shell(old_shell).unwrap().get_faceuses().unwrap(),
            [second]
        );
        assert_eq!(b.shell(new_shell).unwrap().get_faceuses().unwrap(), [first]);
        assert_eq!(b.region(region).unwrap().shells, [old_shell]);
        assert_eq!(b.region(new_region).unwrap().shells, [new_shell]);
        assert_eq!(b.infinite_region(), Some(region));
        assert!(b.editing_enabled());
        assert!(!k.context.is_mark_locked(MarkType::Mark));
    }

    #[test]
    fn infinite_reassignment_uses_live_brep_rays_and_native_swap_or_error_state() {
        use crate::curve::curve::Curve;
        use crate::curve::line::Line;
        use crate::geom::extent1d::Extent1d;
        use crate::geom::extent3d::Extent3d;
        use crate::topology::shell::{ShellRayAttempt, ShellRayInput};

        for fail in [false, true] {
            let mut k = Kernel::default();
            let mut b = Brep::new_with_infinite_region().unwrap();
            let old = b.infinite_region().unwrap();
            let new = b
                .add_region_record(Region {
                    in_brep: true,
                    is_void: 0,
                    ..Region::default()
                })
                .unwrap();
            let shell = b
                .add_shell_record(Shell {
                    region: Some(new),
                    shell_type: ShellType::Wire,
                    ..Shell::default()
                })
                .unwrap();
            b.region_mut(new).unwrap().shells.push(shell);
            let start = b.add_vertex(Point3::ZERO).unwrap();
            let end = b.add_vertex(Point3::new(1.0, 0.0, 0.0)).unwrap();
            let edge = b.add_edge(start, end).unwrap();
            let record = b.edge_mut(edge).unwrap();
            record.curve = Some(Curve::from(Line::from_points(
                Point3::ZERO,
                Point3::new(1.0, 0.0, 0.0),
                3,
            )));
            record.interval = Extent1d::from_min_max(0.0, 1.0);
            let mut bbox =
                Extent3d::from_points(Point3::new(-2.0, -3.0, -4.0), Point3::new(2.0, 3.0, 4.0));
            let expected_first_start = bbox.max() + (bbox.max() - bbox.min()) * 0.5;
            let mut seen = Vec::new();
            let mut queries = |_: &mut Kernel,
                               brep: &mut Brep,
                               selected: Option<ShellId>,
                               mut input: ShellRayInput<EntityRef>| {
                // Wire-shell selection must use the complete B-rep branch.
                assert_eq!(selected, None);
                assert_eq!(brep.regions().len(), 2);
                assert_eq!(input.sphere_bound, None);
                seen.push(input.start);
                brep.set_editing_enabled(false);
                input.state.region = Some(new);
                ShellRayAttempt {
                    status: if fail {
                        Err(SolidError::sm_err("test::native_ray"))
                    } else {
                        Ok(())
                    },
                    state: input.state,
                }
            };
            b.set_editing_enabled(true);
            let result = reassign_split_infinite_region(
                &mut k,
                &mut b,
                old,
                shell,
                &mut bbox,
                true,
                &mut queries,
            );
            assert_eq!(seen.first(), Some(&expected_first_start));
            assert!(!b.editing_enabled());
            if fail {
                assert_eq!(seen.len(), 3);
                assert_eq!(result, Err(SolidError::sm_err("topology::brep::p9")));
                assert_eq!(b.infinite_region(), Some(old));
                assert_eq!(b.region(old).unwrap().is_void(), 1);
                assert_eq!(b.region(new).unwrap().is_void(), 0);
            } else {
                result.unwrap();
                assert_eq!(seen.len(), 1);
                assert_eq!(b.infinite_region(), Some(new));
                assert_eq!(super::super::p1::regions(&b, None).unwrap(), [new, old]);
                assert_eq!(b.region(old).unwrap().is_void(), 0);
                assert_eq!(b.region(new).unwrap().is_void(), 1);
            }
        }
    }

    #[test]
    fn deleting_empty_split_shell_preserves_surviving_ids_and_updates_infinite_head() {
        let mut b = Brep::new();
        let empty = b
            .add_region_record(Region {
                in_brep: true,
                ..Region::default()
            })
            .unwrap();
        let remaining = b
            .add_region_record(Region {
                in_brep: true,
                ..Region::default()
            })
            .unwrap();
        let shell = b
            .add_shell_record(Shell {
                region: Some(empty),
                shell_type: ShellType::Faceuse,
                ..Shell::default()
            })
            .unwrap();
        b.region_mut(empty).unwrap().shells.push(shell);
        b.set_region_list_head(Some(empty));
        b.set_infinite_region_record(Some(empty));
        unlink_empty_split_shell(&mut b, shell).unwrap();
        assert_eq!(super::super::p1::regions(&b, None).unwrap(), [remaining]);
        assert_eq!(b.infinite_region(), Some(remaining));
        assert_eq!(b.shell(shell).unwrap().region, None);
        assert_eq!(b.regions().len(), 2);
        assert_eq!(b.shells().len(), 1);
    }

    #[test]
    fn ignored_native_status_preserves_outputs_but_never_hides_rust_refusals() {
        let partial = || RegionSplitAttempt {
            status: Ok(()),
            shell: Some(ShellId::from_index(7)),
            region: Some(RegionId::from_index(11)),
        };
        let error = SolidError::sm_err("native::split");
        let attempt = finish_split_attempt(partial(), Err(error.clone())).unwrap();
        assert_eq!(attempt.status, Err(error));
        assert_eq!(attempt.shell, Some(ShellId::from_index(7)));
        assert_eq!(attempt.region, Some(RegionId::from_index(11)));
        for error in [
            SolidError::invalid_input(OP, "typed resolver required"),
            SolidError::unported("SmBrep::CalculateBoundingBox", "query.brep"),
            SolidError::NotFound {
                operation: OP,
                entity: EntityRef::Face(FaceId::from_index(23)),
            },
        ] {
            assert_eq!(
                finish_split_attempt(partial(), Err(error.clone())).unwrap_err(),
                error
            );
            assert_eq!(ignored_native::<()>(Err(error.clone())).unwrap_err(), error);
        }
        assert_eq!(
            ignored_native::<()>(Err(SolidError::sm_err(OP))).unwrap(),
            None
        );
    }

    #[test]
    fn missing_split_links_are_outer_errors_even_for_unchecked_callers() {
        let mut b = Brep::new();
        let seed = b.add_faceuse_record(Faceuse::default()).unwrap();
        let mut state = crate::topology::brep::p2::RegionNotifyState::default();
        let result = test_closure_split_region_attempt(
            &mut Kernel::default(),
            &mut b,
            &mut Marks::default(),
            &mut state.env(),
            true,
            seed,
            &mut no_ray,
            crate::intersect::curve_class::p7::FaceQueryIdentity {
                keys: &crate::topology::brep_cache::arena_identity,
                owner_of: &|_| None,
            },
        );
        assert!(matches!(result, Err(SolidError::InvalidInput { .. })));
        assert_eq!(b.shells().len(), 0);
        assert_eq!(b.regions().len(), 0);
    }

    #[test]
    fn extreme_loop_counts_fail_before_allocation_from_the_count() {
        let mut b = Brep::new();
        let region = b.add_region_record(Region::default()).unwrap();
        let surface = Surface::from(Plane::from_normal(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ));
        let options = FaceTopologyOptions {
            region: Some(region),
            edge_loops: &[usize::MAX],
            ..FaceTopologyOptions::default()
        };
        assert!(matches!(
            make_face_topology(
                &mut b,
                Some(surface),
                OrientType::Same,
                ObjectKey(999),
                &options
            ),
            Err(SolidError::InvalidInput { .. })
        ));
        assert_eq!(b.faces().len(), 0);
        assert_eq!(b.loopuses().len(), 0);
    }
}
