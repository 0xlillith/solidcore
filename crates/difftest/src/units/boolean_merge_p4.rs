//! Real native `ManifoldBoolean` versus the Rust orchestration entry point.
//! Native input graphs and query metadata are imported from the same paired
//! constructors. The real kernel adapter delegates the ported later stages;
//! genuine missing dependencies still defer at their exact call sites.
use super::topology_shell::RayQueryContext;
use super::topology_vertexuse::{apply_query_metadata, QueryMetadata};
use crate::{Case, HarnessError, Output, Rng, UnitDiff};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use solidcore::boolean::merge::p1::MergeOptions;
use solidcore::boolean::merge::p4::kernel_adapter::{query_key, OperandGeometryCheck};
use solidcore::boolean::merge::p4::{BooleanKernels, KernelManifoldIntersector, ManifoldState};
use solidcore::geom::a_object::BooleanOperation;
use solidcore::geom::attribute::AttributeStore;
use solidcore::geom::context::{ContextMetadata, MarkLock};
use solidcore::geom::object::ObjectRef;
use solidcore::kernel::Kernel;
use solidcore::topology::brep::{p2 as brep2, p6::RegionSplitQueries};
use solidcore::topology::brep_cache::{self, BrepCache, CacheObject, CacheSource};
use solidcore::topology::shell::{ShellRayAttempt, ShellRayInput, ShellRayState};
use solidcore::topology::topology::Topology;
use solidcore::topology::topology_traverser::{Entity, Marks};
use solidcore::topology::Brep;
use solidcore::topology::{BrepId, ShellId};
use solidcore::EntityRef;
use solidcore::SolidError;
use solidcore_ffi::units::boolean_merge_p4::{self as ffi, Operand};
use std::collections::BTreeMap;

/// Differential runner for the manifold Boolean unit.
pub struct Runner;

#[derive(Debug, Serialize, Deserialize)]
struct Inputs {
    primary: Operand,
    other: Operand,
    operation: i32,
    #[serde(default = "defaults")]
    flags: [bool; 5],
    #[serde(default = "yes")]
    manifold: bool,
    #[serde(default)]
    distinct_contexts: bool,
    #[serde(default)]
    caller_lock: bool,
}
fn defaults() -> [bool; 5] {
    [false, false, true, false, false]
}
fn yes() -> bool {
    true
}
fn error(e: &(impl ToString + ?Sized)) -> HarnessError {
    HarnessError::new(e.to_string())
}
fn request(i: &Inputs) -> ffi::Request {
    ffi::Request {
        operation: i.operation,
        flags: i.flags,
        manifold: i.manifold,
        distinct_contexts: i.distinct_contexts,
        caller_lock: i.caller_lock,
    }
}
fn context_metadata(native: &ffi::ContextState) -> ContextMetadata {
    ContextMetadata {
        current_marks: native.current_marks,
        mark_locks: native.mark_locks,
        doing_boolean: native.doing_boolean,
        model_size_estimate: native.model_size_estimate,
        large_small_size_ratio: native.large_small_size_ratio,
        zone_tol_3d: native.zone_tol_3d,
    }
}

fn primary_key(object: CacheObject) -> ObjectRef {
    query_key(0, false, object)
}
fn other_key(object: CacheObject) -> ObjectRef {
    query_key(0, true, object)
}

struct ImportedOperand {
    brep: Brep,
    metadata: QueryMetadata,
    notification_records: Vec<ffi::NotificationRecord>,
    callback_present: bool,
    attributes_present: bool,
    surface_caches: Vec<ffi::SurfaceCacheState>,
}

fn import_fixture(
    fixture: Option<ffi::Fixture>,
    keys: fn(CacheObject) -> ObjectRef,
) -> Result<Option<ImportedOperand>, HarnessError> {
    let Some(fixture) = fixture else {
        return Ok(None);
    };
    let (mut brep, _) = super::topology_vertexuse::import(&fixture.graph)?;
    let metadata = apply_query_metadata(&mut brep, &fixture.query_metadata, &keys)?;
    if fixture.surface_caches.len() != brep.faces().len() {
        return Err(error(
            "p4 surface-cache coverage differs from imported faces",
        ));
    }
    if fixture
        .surface_caches
        .iter()
        .any(|cache| !cache.present && cache.face_was_modified)
    {
        return Err(error("p4 absent surface cache has modified-face state"));
    }
    Ok(Some(ImportedOperand {
        brep,
        metadata,
        notification_records: fixture.notification_records,
        callback_present: fixture.callback_present,
        attributes_present: fixture.attributes_present,
        surface_caches: fixture.surface_caches,
    }))
}

/// Graph positions are checked before assigning native records to typed ids.
fn import_notifications(
    adapter: &mut KernelManifoldIntersector,
    other: bool,
    input: &ImportedOperand,
    marks: &mut Marks,
) -> Result<(), HarnessError> {
    let b = &input.brep;
    let mut expected = vec![(0, 0, EntityRef::Brep(BrepId::from_index(u32::from(other))))];
    macro_rules! append {
        ($kind:expr, $arena:expr, $variant:ident) => {
            for (index, id) in $arena.ids().enumerate() {
                expected.push((
                    $kind,
                    u32::try_from(index).map_err(|e| error(&e))?,
                    EntityRef::$variant(id),
                ));
            }
        };
    }
    append!(1, b.regions(), Region);
    append!(2, b.shells(), Shell);
    append!(3, b.faceuses(), Faceuse);
    append!(4, b.faces(), Face);
    append!(5, b.loopuses(), Loopuse);
    append!(6, b.loops(), Loop);
    append!(7, b.edgeuses(), Edgeuse);
    append!(8, b.edges(), Edge);
    append!(9, b.vertexuses(), Vertexuse);
    append!(10, b.vertices(), Vertex);
    if expected.len() != input.notification_records.len() {
        return Err(error("p4 notification count differs from imported graph"));
    }
    for ((kind, index, entity), native) in expected.into_iter().zip(&input.notification_records) {
        if (native.kind, native.index) != (kind, index) {
            return Err(error(
                "p4 notification identity differs from imported graph",
            ));
        }
        let [flags, m0, m1, m2, m3, m4, user1, user2, userptr] = native.values;
        let mut record = Topology::new();
        record.set_all_flags(flags);
        record.restore_marks([m0, m1, m2, m3, m4]);
        let traversal_entity = match entity {
            EntityRef::Face(id) => Some(Entity::Face(id)),
            EntityRef::Edge(id) => Some(Entity::Edge(id)),
            EntityRef::Vertex(id) => Some(Entity::Vertex(id)),
            EntityRef::Faceuse(id) => Some(Entity::Faceuse(id)),
            EntityRef::Vertexuse(id) => Some(Entity::Vertexuse(id)),
            EntityRef::Shell(id) => Some(Entity::Shell(id)),
            _ => None,
        };
        if let Some(entity) = traversal_entity {
            marks.restore_marks(entity, [m0, m1, m2, m3, m4]);
        }
        record.set_user_index1(user1);
        record.set_user_index2(user2);
        record.set_user_ptr1(userptr);
        let object = adapter.notification_object(other, entity);
        if adapter.records.insert(object.key, record).is_some() {
            return Err(error("p4 notification namespace collision"));
        }
    }
    Ok(())
}

/// Read owner links installed from the native export, never synthesize an owner.
fn owner_links(brep: &Brep, keys: fn(CacheObject) -> ObjectRef) -> BTreeMap<ObjectRef, ObjectRef> {
    let mut links = BTreeMap::new();
    for (id, edge) in brep.edges().iter() {
        if let Some(owner) = edge
            .curve
            .as_ref()
            .and_then(solidcore::curve::curve::Curve::owner)
        {
            links.insert(keys(CacheObject::Curve(id)), owner);
        }
    }
    for (id, face) in brep.faces().iter() {
        if let Some(owner) = face
            .surface
            .as_ref()
            .and_then(solidcore::surface::surface::Surface::owner)
        {
            links.insert(keys(CacheObject::Surface(id)), owner);
        }
    }
    links
}

#[derive(Debug)]
struct ImportedGeometryCheck(QueryMetadata);
impl OperandGeometryCheck for ImportedGeometryCheck {
    fn check(&self) -> solidcore::Result<()> {
        self.0.require_geometry()
    }
    fn check_edge_curve(&self, edge: solidcore::topology::EdgeId) -> solidcore::Result<()> {
        self.0.require_curve_geometry(edge)
    }
}

#[derive(Debug)]
struct ImportedNotificationCheck {
    callback_present: bool,
    attributes_present: bool,
}
impl OperandGeometryCheck for ImportedNotificationCheck {
    fn check(&self) -> solidcore::Result<()> {
        if self.callback_present || self.attributes_present {
            Err(SolidError::unported(
                "native p4 notification attributes or callbacks are not represented",
                "harness:topology-notification-export",
            ))
        } else {
            Ok(())
        }
    }
}

/// A warm native cache is observed but its retained tree is not imported.
/// The Rust adapter's actual flow retains no trim caches between query checkouts.
#[derive(Debug)]
struct NativeSurfaceCacheCheck {
    // Preserve native observations in exact graph-face order.
    states: Vec<(bool, bool)>,
}
impl OperandGeometryCheck for NativeSurfaceCacheCheck {
    fn check(&self) -> solidcore::Result<()> {
        if self.states.iter().any(|(present, _)| *present) {
            Err(SolidError::unported(
                "native retained trim-surface cache contents are not represented",
                "harness:trim-surface-cache-export",
            ))
        } else {
            Ok(())
        }
    }
}

/// Explicitly bound to the primary operand, including its actual absence.
/// There is no empty-B-rep substitute and no operand inference from addresses.
struct PrimaryQueries {
    shell: Option<RayQueryContext>,
    cache: Option<BrepCache>,
    owners: BTreeMap<ObjectRef, ObjectRef>,
}
fn make_ray_context(input: &ImportedOperand, other: bool) -> Result<RayQueryContext, HarnessError> {
    let keys: fn(CacheObject) -> ObjectRef = if other { other_key } else { primary_key };
    let mut context = RayQueryContext::new(
        &input.brep,
        BrepId::from_index(u32::from(other)),
        input.metadata.clone(),
        keys,
    );
    for native in &input.notification_records {
        let index = usize::try_from(native.index).map_err(|e| error(&e))?;
        let member = match native.kind {
            4 => Some(CacheObject::Face(
                input
                    .brep
                    .faces()
                    .ids()
                    .nth(index)
                    .ok_or_else(|| error("ray mark face outside imported graph"))?,
            )),
            8 => Some(CacheObject::Edge(
                input
                    .brep
                    .edges()
                    .ids()
                    .nth(index)
                    .ok_or_else(|| error("ray mark edge outside imported graph"))?,
            )),
            10 => Some(CacheObject::Vertex(
                input
                    .brep
                    .vertices()
                    .ids()
                    .nth(index)
                    .ok_or_else(|| error("ray mark vertex outside imported graph"))?,
            )),
            _ => None,
        };
        if let Some(member) = member {
            let [_, m0, m1, m2, m3, m4, _, _, _] = native.values;
            context
                .restore_member_marks(member, [m0, m1, m2, m3, m4])
                .map_err(|e| error(&e))?;
        }
    }
    Ok(context)
}

impl PrimaryQueries {
    // One fresh provider per manifold call. Core consumes it only at the final
    // editing-stage FindAndSetInfiniteRegion; no 3D topology mutation occurs
    // between its internal retries. Never reuse this cache across manifold calls.
    fn new(primary: Option<&ImportedOperand>) -> Result<Self, HarnessError> {
        let queries = Self {
            shell: primary
                .map(|input| make_ray_context(input, false))
                .transpose()?,
            cache: None,
            owners: primary
                .map_or_else(BTreeMap::new, |input| owner_links(&input.brep, primary_key)),
        };
        Ok(queries)
    }
}
impl RegionSplitQueries for PrimaryQueries {
    fn check_geometry(
        &mut self,
        brep: &Brep,
        curves: bool,
        surfaces: bool,
    ) -> solidcore::Result<()> {
        self.shell
            .as_mut()
            .ok_or_else(|| SolidError::null_pointer("p4 primary ray operand"))?
            .check_geometry(brep, curves, surfaces)
    }
    fn ray_fire(
        &mut self,
        kernel: &mut Kernel,
        brep: &mut Brep,
        shell: ShellId,
        input: ShellRayInput<EntityRef>,
    ) -> ShellRayAttempt<EntityRef> {
        match self.shell.as_mut() {
            Some(context) => context.ray_fire(kernel, brep, shell, input),
            None => ShellRayAttempt {
                status: Err(SolidError::null_pointer("p4 primary ray operand")),
                state: input.state,
            },
        }
    }

    fn brep_ray_fire(
        &mut self,
        kernel: &mut Kernel,
        brep: &mut Brep,
        input: ShellRayInput<EntityRef>,
    ) -> ShellRayAttempt<EntityRef> {
        let mut output = brep2::RayFireResult {
            region: input.state.region,
            intersection: input.state.classification,
        };
        let status = (|| {
            self.shell
                .as_ref()
                .ok_or_else(|| SolidError::null_pointer("p4 primary ray operand"))?
                .metadata()
                .require_geometry()?;
            let cache = &mut self.cache;
            let owners = &self.owners;
            brep2::ray_fire_with_cache_into(
                kernel,
                brep,
                BrepId::from_index(0),
                |brep| {
                    if cache.is_none() {
                        *cache = Some(brep_cache::make_brep_cache(
                            kernel,
                            CacheSource::Brep(brep),
                            &primary_key,
                        )?);
                    }
                    cache.as_ref().ok_or_else(|| {
                        SolidError::invalid_input(
                            "p4 primary ray cache",
                            "missing constructed cache",
                        )
                    })
                },
                &primary_key,
                &|object| owners.get(&object).copied(),
                input.start,
                input.direction,
                input.zone_tolerance,
                None,
                &mut output,
            )
        })();
        ShellRayAttempt {
            status,
            state: ShellRayState {
                region: output.region,
                classification: output.intersection,
            },
        }
    }
}

fn append_graph(output: &mut Output, prefix: &str, brep: &Brep) -> Result<(), HarnessError> {
    let graph = super::topology_brep_p10::snapshot(brep)?;
    for (key, values) in graph.ints {
        output.ints.insert(format!("{prefix}.{key}"), values);
    }
    for (key, values) in graph.values {
        output.values.insert(format!("{prefix}.{key}"), values);
    }
    Ok(())
}

struct Record<'a> {
    remaining: &'a [f64],
}
impl<'a> Record<'a> {
    fn num(&mut self) -> Result<f64, HarnessError> {
        let (v, rest) = self
            .remaining
            .split_first()
            .ok_or_else(|| error("short native record"))?;
        self.remaining = rest;
        Ok(*v)
    }
    fn take(&mut self, n: usize) -> Result<&'a [f64], HarnessError> {
        if n > self.remaining.len() {
            return Err(error("short native graph"));
        }
        let (v, rest) = self.remaining.split_at(n);
        self.remaining = rest;
        Ok(v)
    }
    fn done(&self) -> Result<(), HarnessError> {
        if self.remaining.is_empty() {
            Ok(())
        } else {
            Err(error("trailing native state"))
        }
    }
}

fn native_snapshot(raw: &[f64]) -> Result<Output, HarnessError> {
    let mut rd = Record { remaining: raw };
    let mut header = Vec::new();
    for _ in 0..9 {
        header.push(rd.num().map_err(|e| error(&e))? as i64);
    }
    let component_count = header
        .get(7)
        .zip(header.get(8))
        .and_then(|(&a, &b)| a.checked_add(b))
        .and_then(|n| usize::try_from(n).ok())
        .ok_or_else(|| error("invalid component count"))?;
    let mut output = Output::new().ints("status_and_state", header);
    let mut graph = |rd: &mut Record<'_>, prefix: &str| -> Result<(), HarnessError> {
        let n = rd.num().map_err(|e| error(&e))? as usize;
        let record = rd.take(n).map_err(|e| error(&e))?;
        let (editing, body) = record
            .split_first()
            .ok_or_else(|| error("empty native graph"))?;
        let mut exported = vec![0.0];
        exported.extend_from_slice(body);
        let (mut b, exact) = super::topology_vertexuse::import(&exported)?;
        if exact.iter().any(|v| !v) {
            return Err(error("inexact post-state surface import"));
        }
        b.set_editing_enabled(*editing != 0.0);
        append_graph(&mut output, prefix, &b)
    };
    let mut presence = Vec::new();
    for prefix in ["primary", "other"] {
        let present = rd.num().map_err(|e| error(&e))? != 0.0;
        presence.push(i64::from(present));
        if present {
            graph(&mut rd, prefix)?;
        }
    }
    for i in 0..component_count {
        graph(&mut rd, &format!("component.{i}"))?;
    }
    rd.done().map_err(|e| error(&e))?;
    Ok(output.ints("survivors", presence))
}

fn rust_snapshot(
    status: i32,
    returned: bool,
    kernels: &BooleanKernels<'_>,
    primary: Option<&Brep>,
    other: Option<&Brep>,
    state: &ManifoldState,
) -> Result<Output, HarnessError> {
    let other_kernel = kernels.other.as_deref().unwrap_or(kernels.primary);
    let mut output = Output::new()
        .ints(
            "status_and_state",
            vec![
                i64::from(status),
                i64::from(returned),
                i64::from(kernels.primary.context.doing_boolean()),
                i64::from(other_kernel.context.doing_boolean()),
                i64::from(state.intersection_loop_closed),
                state.intersection_edges.len() as i64,
                state.intersection_vertices.len() as i64,
                state.primary_components.len() as i64,
                state.other_components.len() as i64,
            ],
        )
        .ints(
            "survivors",
            vec![i64::from(primary.is_some()), i64::from(other.is_some())],
        );
    for (prefix, b) in [("primary", primary), ("other", other)] {
        if let Some(b) = b {
            append_graph(&mut output, prefix, b)?;
        }
    }
    for (i, b) in state
        .primary_components
        .iter()
        .chain(&state.other_components)
        .enumerate()
    {
        append_graph(&mut output, &format!("component.{i}"), b)?;
    }
    Ok(output)
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "boolean.merge_p4"
    }
    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let raw = ffi::eval(&i.primary, &i.other, request(&i)).map_err(|e| error(&e))?;
        if raw == [-1.0] {
            return Ok(Output::rejected("native primitive constructor"));
        }
        // Status is data: Output::err would suppress all state comparison.
        native_snapshot(&raw)
    }
    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let operation = BooleanOperation::from_code(i.operation)
            .ok_or_else(|| error("unknown Boolean operation"))?;
        let input_pair = ffi::export_fixtures(&i.primary, &i.other, i.distinct_contexts)
            .map_err(|e| error(&e))?;
        let [primary_fixture, other_fixture] = input_pair.operands;
        let [primary_context, other_context] = input_pair.contexts;
        let primary_input = import_fixture(primary_fixture, primary_key)?;
        let other_input = import_fixture(other_fixture, other_key)?;
        for input in [primary_input.as_ref(), other_input.as_ref()]
            .into_iter()
            .flatten()
        {
            if let Some(output) = input.metadata.deferred_output_for(false, false) {
                return Ok(output);
            }
        }
        let mut primary_queries = PrimaryQueries::new(primary_input.as_ref())?;
        let primary_metadata = primary_input.as_ref().map(|input| input.metadata.clone());
        let other_metadata = other_input.as_ref().map(|input| input.metadata.clone());
        let mut k = Kernel::default();
        let mut other_k = Kernel::default();
        k.context
            .restore_metadata(context_metadata(&primary_context));
        other_k
            .context
            .restore_metadata(context_metadata(&other_context));
        let mut kernels = BooleanKernels {
            primary: &mut k,
            other: if i.distinct_contexts {
                Some(&mut other_k)
            } else {
                None
            },
        };
        let mut state = ManifoldState::default();
        let mut attributes = AttributeStore::new();
        let mut adapter = KernelManifoldIntersector::new(0).map_err(|e| error(&e.to_string()))?;
        adapter.set_intersection_tolerances(1e-9, 0.3);
        if let Some(input) = primary_input.as_ref() {
            import_notifications(&mut adapter, false, input, &mut state.primary_marks)?;
            adapter.set_region_queries(false, Box::new(make_ray_context(input, false)?));
            adapter.set_notification_check(
                false,
                Box::new(ImportedNotificationCheck {
                    callback_present: input.callback_present,
                    attributes_present: input.attributes_present,
                }),
            );
        }
        if let Some(input) = other_input.as_ref() {
            import_notifications(&mut adapter, true, input, &mut state.other_marks)?;
            adapter.set_region_queries(true, Box::new(make_ray_context(input, true)?));
            adapter.set_notification_check(
                true,
                Box::new(ImportedNotificationCheck {
                    callback_present: input.callback_present,
                    attributes_present: input.attributes_present,
                }),
            );
        }
        // This is the actual Rust flow's empty retained-cache registry, not a
        // claim that native queries cannot create caches later. A reached warm
        // native input cache stays an explicit representation limitation.
        if let Some(input) = primary_input.as_ref() {
            adapter.set_surface_cache_registry(
                BTreeMap::new(),
                Some(Box::new(NativeSurfaceCacheCheck {
                    states: input
                        .surface_caches
                        .iter()
                        .map(|cache| (cache.present, cache.face_was_modified))
                        .collect(),
                })),
            );
        }
        let mut primary = primary_input.map(|input| input.brep);
        let mut other = other_input.map(|input| input.brep);
        if let Some(brep) = primary.as_ref() {
            adapter.owner_links.extend(owner_links(brep, primary_key));
        }
        if let Some(brep) = other.as_ref() {
            adapter.owner_links.extend(owner_links(brep, other_key));
        }
        if let Some(metadata) = primary_metadata {
            adapter.set_query_geometry_check(false, Box::new(ImportedGeometryCheck(metadata)));
        }
        if let Some(metadata) = other_metadata {
            adapter.set_query_geometry_check(true, Box::new(ImportedGeometryCheck(metadata)));
        }
        let mut lock = MarkLock::unbound();
        let [cookie_cutter, imprinting, boolean_post_process, imprint_and_classify_faces, keep_other_brep] =
            i.flags;
        let options = MergeOptions {
            cookie_cutter,
            imprinting,
            boolean_post_process,
            imprint_and_classify_faces,
            keep_other_brep,
        };
        adapter.merge.set_merge_options(Some(&options));
        let result = adapter.manifold_boolean_with_queries(
            &mut kernels,
            &mut primary,
            &mut other,
            &mut state,
            &mut attributes,
            i.manifold,
            operation,
            if i.caller_lock { Some(&mut lock) } else { None },
            &mut primary_queries,
        );
        match result {
            Err(e) if e.is_unported() => {
                Ok(Output::unported(&e).ok_or_else(|| error("missing deferral"))?)
            }
            Err(SolidError::Upstream { code, .. }) => rust_snapshot(
                code,
                false,
                &kernels,
                primary.as_ref(),
                other.as_ref(),
                &state,
            ),
            Err(e) => Err(error(&e)),
            Ok(returned) => rust_snapshot(
                0,
                returned,
                &kernels,
                primary.as_ref(),
                other.as_ref(),
                &state,
            ),
        }
    }
    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        let scale = match rng.below(5) {
            0 => 1e-9,
            1 => 1e-3,
            2 => 1.0,
            3 => 1e3,
            _ => 1e6,
        };
        let mut operand = || {
            let kind = rng.below(4) as f64;
            let [x, y, z] = rng.point3(-2.0 * scale, 2.0 * scale);
            [
                kind,
                x,
                y,
                z,
                rng.uniform(0.5, 2.0) * scale,
                rng.uniform(0.5, 2.0) * scale,
                rng.uniform(0.5, 2.0) * scale,
            ]
        };
        let a = operand();
        let b = operand();
        Some(json!({"primary":a,"other":b,"operation":rng.below(11),
            "flags":[rng.chance(0.2),rng.chance(0.2),rng.chance(0.5),rng.chance(0.2),rng.chance(0.5)],
            "manifold":rng.chance(0.8),"distinct_contexts":rng.chance(0.5),"caller_lock":rng.chance(0.5)}))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use solidcore::topology::VertexId;

    #[test]
    fn error_post_state_detects_mutations() {
        let mut k = Kernel::default();
        let mut other_k = Kernel::default();
        let pair = ffi::export_fixtures(
            &[1.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0],
            &[-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            true,
        )
        .unwrap();
        let [primary_fixture, _] = pair.operands;
        let [primary_context, other_context] = pair.contexts;
        k.context
            .restore_metadata(context_metadata(&primary_context));
        other_k
            .context
            .restore_metadata(context_metadata(&other_context));
        let mut primary = import_fixture(primary_fixture, primary_key)
            .unwrap()
            .map(|input| input.brep);
        let other: Option<Brep> = None;
        let mut state = ManifoldState::default();
        let mut snapshot = |k: &mut Kernel, primary: Option<&Brep>, state: &ManifoldState| {
            rust_snapshot(
                1008,
                false,
                &BooleanKernels {
                    primary: k,
                    other: Some(&mut other_k),
                },
                primary,
                other.as_ref(),
                state,
            )
            .unwrap()
        };
        let baseline = snapshot(&mut k, primary.as_ref(), &state);
        k.context.set_doing_boolean(true);
        assert_ne!(baseline, snapshot(&mut k, primary.as_ref(), &state));
        k.context.set_doing_boolean(false);
        primary.as_mut().unwrap().set_editing_enabled(true);
        assert_ne!(baseline, snapshot(&mut k, primary.as_ref(), &state));
        primary.as_mut().unwrap().set_editing_enabled(false);
        state.intersection_vertices.push(VertexId::from_index(0));
        assert_ne!(baseline, snapshot(&mut k, primary.as_ref(), &state));
        state.intersection_vertices.clear();
        let b = primary.as_mut().unwrap();
        let id = b.vertices().ids().next().unwrap();
        b.vertex_mut(id).unwrap().point.x += 1.0;
        assert_ne!(baseline, snapshot(&mut k, primary.as_ref(), &state));
        primary = None;
        assert_ne!(baseline, snapshot(&mut k, primary.as_ref(), &state));
    }
}

#[cfg(test)]
mod imported_notification_tests {
    use super::*;

    #[test]
    fn unrepresented_native_handlers_fail_only_when_checked() {
        for (callback_present, attributes_present) in [(false, true), (true, false), (true, true)] {
            let check = ImportedNotificationCheck {
                callback_present,
                attributes_present,
            };
            assert!(matches!(
                check.check(),
                Err(SolidError::Unported {
                    unit: "harness:topology-notification-export",
                    ..
                })
            ));
        }
        assert!(ImportedNotificationCheck {
            callback_present: false,
            attributes_present: false,
        }
        .check()
        .is_ok());
    }

    #[test]
    fn native_records_keep_wide_fields_and_operand_identity() {
        let make = || ImportedOperand {
            brep: Brep::new(),
            metadata: QueryMetadata::default(), // No query runs in this record-decoder test.
            notification_records: vec![ffi::NotificationRecord {
                kind: 0,
                index: 0,
                values: [u64::MAX, 1, 2, 3, 4, 5, 1_u64 << 60, 17, u64::MAX - 1],
            }],
            callback_present: false,
            attributes_present: false,
            surface_caches: Vec::new(), // No faces in this decoder-only fixture.
        };
        let mut adapter = KernelManifoldIntersector::new(0).unwrap();
        let first = make();
        import_notifications(&mut adapter, false, &first, &mut Marks::default()).unwrap();
        import_notifications(&mut adapter, true, &make(), &mut Marks::default()).unwrap();
        let a = adapter.notification_object(false, EntityRef::Brep(BrepId::from_index(0)));
        let b = adapter.notification_object(true, EntityRef::Brep(BrepId::from_index(1)));
        assert_ne!(a.key, b.key);
        let record = adapter.records.get(&a.key).unwrap();
        assert_eq!(record.all_flags(), u64::MAX);
        assert_eq!(record.user_index1(), 1_u64 << 60);
        assert_eq!(record.user_ptr1(), u64::MAX - 1);
        let mut invalid = make();
        invalid.notification_records[0].kind = 10;
        assert!(import_notifications(
            &mut KernelManifoldIntersector::new(0).unwrap(),
            false,
            &invalid,
            &mut Marks::default(),
        )
        .is_err());
        let mut invalid = make();
        invalid.notification_records.clear();
        assert!(import_notifications(
            &mut KernelManifoldIntersector::new(0).unwrap(),
            false,
            &invalid,
            &mut Marks::default(),
        )
        .is_err());
    }
}

#[cfg(test)]
mod imported_surface_cache_tests {
    use super::*;

    #[test]
    fn observed_warm_cache_is_not_replaced_with_a_cold_assumption() {
        for modified in [false, true] {
            assert!(matches!(
                NativeSurfaceCacheCheck {
                    states: vec![(false, false), (true, modified)],
                }
                .check(),
                Err(SolidError::Unported {
                    unit: "harness:trim-surface-cache-export",
                    ..
                })
            ));
        }
        assert!(NativeSurfaceCacheCheck {
            states: vec![(false, false)]
        }
        .check()
        .is_ok());
    }
}

#[cfg(test)]
mod imported_traversal_mark_tests {
    use super::*;
    use solidcore::geom::context::MarkType;
    use solidcore::math::Vec3;
    use solidcore::topology::{Faceuse, Vertexuse};

    #[test]
    fn six_supported_classes_restore_marks_in_both_operand_stores() {
        let mut brep = Brep::new();
        let shell = brep.add_shell().unwrap();
        let face = brep.add_face(shell).unwrap();
        let vertex = brep.add_vertex(Vec3::ZERO).unwrap();
        let edge = brep.add_edge(vertex, vertex).unwrap();
        let faceuse = brep.add_faceuse_record(Faceuse::default()).unwrap();
        let vertexuse = brep.add_vertexuse_record(Vertexuse::default()).unwrap();
        let records = [0, 2, 3, 4, 8, 9, 10].map(|kind| ffi::NotificationRecord {
            kind,
            index: 0,
            values: [0, 7, 11, 13, 17, 19, 0, 0, 0],
        });
        let input = ImportedOperand {
            brep,
            metadata: QueryMetadata::default(), // Record decoding only; no geometry query.
            notification_records: records.into_iter().collect(),
            callback_present: false,
            attributes_present: false,
            surface_caches: Vec::new(), // No cache decoder is invoked in this test.
        };
        let mut adapter = KernelManifoldIntersector::new(0).unwrap();
        let mut state = ManifoldState::default();
        import_notifications(&mut adapter, false, &input, &mut state.primary_marks).unwrap();
        import_notifications(&mut adapter, true, &input, &mut state.other_marks).unwrap();
        let mut kernel = Kernel::default();
        let mut native_context = kernel.context.metadata();
        native_context.current_marks = [7, 11, 13, 17, 19];
        kernel.context.restore_metadata(native_context);
        for entity in [
            Entity::Face(face),
            Entity::Edge(edge),
            Entity::Vertex(vertex),
            Entity::Faceuse(faceuse),
            Entity::Vertexuse(vertexuse),
            Entity::Shell(shell),
        ] {
            for mark in [
                MarkType::Mark,
                MarkType::Mark2,
                MarkType::Mark3,
                MarkType::MarkIo,
                MarkType::MarkAssert,
            ] {
                assert!(state.primary_marks.is_marked(&kernel, entity, mark));
                assert!(state.other_marks.is_marked(&kernel, entity, mark));
            }
        }
        assert_eq!(kernel.context.metadata(), native_context);
    }
}
