//! Real native `ManifoldBoolean` versus the Rust orchestration entry point.
//! All non-null Rust inputs currently stop at `merge_p1`'s `query.brep` transform
//! dependency. The later adapter deliberately fails if reached: wiring the
//! prerequisite must not silently turn these cases into fake parity passes.
use crate::{Case, HarnessError, Output, Rng, UnitDiff};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use solidcore::boolean::merge::p1::MergeOptions;
use solidcore::boolean::merge::p4::{
    self, BooleanKernels, ManifoldIntersector, ManifoldState, MergeEntity,
};
use solidcore::geom::a_object::BooleanOperation;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::context::{MarkLock, MarkType};
use solidcore::kernel::Kernel;
use solidcore::topology::{Brep, EdgeId, FaceId, FaceuseId, RegionId, VertexId};
use solidcore::{Result, SolidError};
use solidcore_ffi::units::boolean_merge_p4::{self as ffi, Operand};

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
fn import(input: &Operand) -> std::result::Result<Option<Brep>, HarnessError> {
    if input[0] == -1.0 {
        return Ok(None);
    }
    let raw = ffi::export(input).map_err(|e| error(&e))?;
    if raw == [-1.0] {
        return Err(error("native input constructor rejected"));
    }
    super::topology_vertexuse::import(&raw).map(|(b, _)| Some(b))
}

struct DeferredAdapter;
fn needs_adapter<T>() -> Result<T> {
    Err(SolidError::invalid_input("harness::boolean_merge_p4", "query.brep prerequisite is now wired: integrate the canonical intersector adapter before enabling success comparisons"))
}
#[allow(clippy::panic)] // Harness sentinel: reaching an unwired adapter must fail loudly.
impl ManifoldIntersector for DeferredAdapter {
    fn mapped_edges(&self) -> Result<Vec<EdgeId>> {
        needs_adapter()
    }
    fn mapped_vertices(&self) -> Result<Vec<VertexId>> {
        needs_adapter()
    }
    fn common_edges(&self, _: &Brep, _: &Brep) -> Result<(Vec<EdgeId>, Vec<EdgeId>)> {
        needs_adapter()
    }
    fn common_vertices(&self, _: &Brep, _: &Brep) -> Result<(Vec<VertexId>, Vec<VertexId>)> {
        needs_adapter()
    }
    fn edge_mate(&self, _: EdgeId) -> Option<EdgeId> {
        panic!("integrate canonical intersector adapter")
    }
    fn face_mate(&self, _: FaceId) -> Option<FaceId> {
        panic!("integrate canonical intersector adapter")
    }
    fn relate_faces(&mut self, _: &mut Brep, _: FaceId, _: &mut Brep, _: FaceId) -> Result<()> {
        needs_adapter()
    }
    fn fix_face_coincidence(
        &mut self,
        _: &mut Brep,
        _: &mut Brep,
        _: FaceId,
        _: FaceId,
    ) -> Result<(Vec<FaceId>, Vec<FaceId>, bool)> {
        needs_adapter()
    }
    fn process_wires(
        &mut self,
        _: &mut BooleanKernels<'_>,
        _: &mut Brep,
        _: &mut Brep,
        _: &mut ManifoldState,
        _: BooleanOperation,
        _: MarkType,
        _: MarkType,
        _: &mut Vec<EdgeId>,
    ) -> Result<()> {
        needs_adapter()
    }
    fn find_region_of_other_face(
        &mut self,
        _: &Kernel,
        _: &mut Brep,
        _: &mut Brep,
        _: FaceuseId,
    ) -> Result<Option<RegionId>> {
        needs_adapter()
    }
    fn propagate_regions(
        &mut self,
        _: &mut BooleanKernels<'_>,
        _: &Brep,
        _: &Brep,
        _: &mut AttributeStore,
        _: BooleanOperation,
    ) -> Result<()> {
        needs_adapter()
    }
    fn post_process(
        &mut self,
        _: &mut BooleanKernels<'_>,
        _: &mut Brep,
        _: &ManifoldState,
        _: bool,
    ) -> Result<()> {
        needs_adapter()
    }
    fn object_key(&self, _: bool, _: MergeEntity) -> ObjectKey {
        panic!("integrate canonical intersector adapter")
    }
    fn excluded_surface(&self, _: &Brep, _: FaceId) -> Result<bool> {
        needs_adapter()
    }
    fn clear_topology_deleted(&mut self) {
        panic!("integrate canonical intersector adapter")
    }
    fn topology_deleted(&self) -> bool {
        panic!("integrate canonical intersector adapter")
    }
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "boolean.merge_p4"
    }
    fn run_cpp(&self, case: &Case) -> std::result::Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let raw = ffi::eval(&i.primary, &i.other, request(&i)).map_err(|e| error(&e))?;
        if raw == [-1.0] {
            return Ok(Output::rejected("native primitive constructor"));
        }
        let status = *raw.first().ok_or_else(|| error("missing native status"))?;
        if status != 0.0 {
            return Ok(Output::err(format!("status:{}", status as i32)));
        }
        // No Rust success case reaches this record yet; it is retained as native
        // baseline evidence, never imported as Rust classification or topology.
        Ok(Output::new().values("native_result", raw))
    }
    fn run_rust(&self, case: &Case) -> std::result::Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let operation = BooleanOperation::from_code(i.operation)
            .ok_or_else(|| error("unknown Boolean operation"))?;
        let mut primary = import(&i.primary)?;
        let mut other = import(&i.other)?;
        let mut k = Kernel::default();
        let mut other_k = Kernel::default();
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
        let mut adapter = DeferredAdapter;
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
        let result = p4::manifold_boolean(
            &mut kernels,
            &mut primary,
            &mut other,
            &mut adapter,
            &mut state,
            &mut attributes,
            options,
            i.manifold,
            operation,
            if i.caller_lock { Some(&mut lock) } else { None },
        );
        match result {
            Err(e) if e.is_unported() => {
                Ok(Output::unported(&e).ok_or_else(|| error("missing deferral"))?)
            }
            Err(SolidError::Upstream { code, .. }) => Ok(Output::err(format!("status:{code}"))),
            Err(e) => Err(error(&e)),
            Ok(_) => Err(error(
                "unexpected success before canonical intersector adapter integration",
            )),
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
