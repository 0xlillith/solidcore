//! Differential coverage of topology p2's local tree traversal.
//!
//! The shared p1 fixture adapter builds the same objects and trees on both sides;
//! the p2 C++ entry reaches the pinned `SmTopologySolver::LocalSolve`.

use super::intersect_topology_solver_p1 as shared;
use super::topology_vertexuse::{apply_query_metadata, QueryMetadata};
use crate::{Case, CaseContext, HarnessError, Output, Rng, UnitDiff};
use serde::Deserialize;
use serde_json::Value;
use solidcore::curve::curve::Curve;
use solidcore::error::{Result as SolidResult, SolidError};
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::extent3d::Extent3d;
use solidcore::geom::global_solver::{SolutionRequested, SolverOperation};
use solidcore::geom::object::{Shape, ShapeMember};
use solidcore::geom::solution_array::{ObjectKind, SolutionArray, SolutionObject, SolutionType};
use solidcore::geom::tree::{NodeGeomType, Tree};
use solidcore::intersect::point_class::PointClassification;
use solidcore::intersect::topology_solver::p1::{
    BrepTreeSolve, SolverCurve, TopologyObjects, TopologySolver,
};
use solidcore::intersect::topology_solver::p2::{TopologySolveOptions, TopologyTrees};
use solidcore::kernel::Kernel;
use solidcore::math::{Point3, Vec2, Vec3};
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::crv_basic::Curve as NlCurve;
use solidcore::nurbs::math::{CPoint, NOW};
use solidcore::surface::surface::Surface;
use solidcore::topology::brep_cache::{self, arena_identity, CacheObject};
use solidcore::topology::{Arena, Brep, BrepId, EdgeId, FaceId, VertexId};
use solidcore::EntityRef;
use solidcore_ffi::units::intersect_topology_solver_p2 as ffi;
use std::cell::Cell;

/// Topology solver p2 differential runner.
pub struct Runner;

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "intersect.topology_solver_p2"
    }
    fn cases(&self, ctx: &CaseContext) -> Result<Vec<Case>, HarnessError> {
        crate::runner::default_cases(
            self,
            &CaseContext {
                seeds: ctx.seeds.max(600),
                ..ctx.clone()
            },
        )
    }
    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        if case.inputs.get("op").and_then(Value::as_str) == Some("static") {
            return ffi::eval_static(&static_input(case)?)
                .map(|v| static_output(&v))
                .map_err(|e| HarnessError::new(e.to_string()));
        }
        let output = if case.inputs.get("op").and_then(Value::as_str) == Some("local") {
            ffi::eval_local(&shared::p2_local_input(case)?)
        } else {
            ffi::eval_trees(&shared::p2_tree_input(case)?)
        };
        output
            .map(|v| shared::p2_tree_output(&v))
            .map_err(|e| HarnessError::new(e.to_string()))
    }
    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        if case.inputs.get("op").and_then(Value::as_str) == Some("static") {
            return static_rust(case);
        }
        shared::Runner.run_rust(case)
    }
    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        shared::p2_tree_generate(rng)
    }
}

// The shell-vertex static fixtures retain their independent, equivalent trees.
// Face fixtures import the topology, concrete geometry and owners of the same
// native input before any solve. They invoke local solvers without constructing
// a B-rep cache, just as the reference does.

#[derive(Deserialize)]
struct StaticInputs {
    #[serde(default)]
    swap: bool,
    #[serde(default)]
    included: bool,
    entry: i32,
    operation: i32,
    request: i32,
    tolerance: f64,
    point: [f64; 3],
    first: Vec<[f64; 3]>,
    second: Vec<[f64; 3]>,
}
fn static_best(operation: i32) -> f64 {
    if matches!(operation, 4 | 13 | 15) {
        0.0
    } else {
        1e30
    }
}
fn static_input(case: &Case) -> Result<Vec<f64>, HarnessError> {
    let i: StaticInputs = case.parse()?;
    let mut out = vec![
        f64::from(i.entry),
        f64::from(i.operation),
        f64::from(i.request),
        i.tolerance,
        static_best(i.operation),
        i.first.len() as f64,
        i.second.len() as f64,
    ];
    out.extend(i.point);
    for p in i.first.iter().chain(&i.second) {
        out.extend(p);
    }
    out.push(f64::from(u8::from(i.swap)));
    out.push(f64::from(u8::from(i.included)));
    Ok(out)
}
fn static_output(flat: &[f64]) -> Output {
    Output::new()
        .scalar("status", flat.first().copied().unwrap_or(-1.0))
        .values("static", flat.get(1..).unwrap_or_default().to_vec())
}
/// The shim's `id(object)` of a solution object: the vertex index, `999` for the
/// curve, `-2` for anything else (the faces and surfaces of face/face solves).
fn static_object_id(h: StaticHandle) -> f64 {
    match h {
        StaticHandle::Fixture(id @ 0..=999) => f64::from(id),
        _ => -2.0,
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
enum StaticHandle {
    /// The identities explicitly listed by the native static fixture protocol.
    Fixture(i32),
    /// Other live topology reached through a face classification. The operand
    /// and typed arena record preserve identity without aliasing the query edge.
    Member(usize, CacheObject),
}
impl StaticHandle {
    fn fixture_id(self) -> SolidResult<i32> {
        match self {
            Self::Fixture(id) => Ok(id),
            Self::Member(_, _) => Err(static_error()),
        }
    }
}
impl SolutionObject for StaticHandle {
    fn kind(&self) -> ObjectKind {
        match *self {
            Self::Member(_, object) => object.kind(),
            Self::Fixture(id) => match id {
                999 | -2 => ObjectKind::Curve,
                1000 | 1001 => ObjectKind::Other,
                2000 | 2001 => ObjectKind::Face,
                2002 | 2003 => ObjectKind::Surface,
                2004 => ObjectKind::Edge,
                _ => ObjectKind::Vertex,
            },
        }
    }
    fn natural_interval(&self) -> Option<Extent1d> {
        matches!(self, Self::Fixture(999 | -2)).then(|| Extent1d::from_min_max(0.0, 1.0))
    }
}
struct NativeOperand {
    id: BrepId,
    face: FaceId,
    vertices: Vec<VertexId>,
    curve_edge: Option<EdgeId>,
    metadata: QueryMetadata,
}

struct NativeScene {
    breps: Arena<BrepId, Brep>,
    operands: [NativeOperand; 2],
}

impl NativeScene {
    fn import(input: &[f64]) -> Result<Self, HarnessError> {
        fn operand(
            breps: &mut Arena<BrepId, Brep>,
            fixture: ffi::StaticBrepFixture,
        ) -> Result<NativeOperand, HarnessError> {
            let (mut brep, _) = super::topology_edgeuse_p1::import_with_types(&fixture.graph)?;
            let metadata =
                apply_query_metadata(&mut brep, &fixture.query_metadata, &arena_identity)?;
            let face = fixture
                .face
                .map(FaceId::from_index)
                .ok_or_else(|| HarnessError::new("native static fixture has no selected face"))?;
            brep.face(face)
                .map_err(|e| HarnessError::new(e.to_string()))?;
            let vertices = fixture
                .vertices
                .into_iter()
                .map(VertexId::from_index)
                .collect();
            let curve_edge = fixture.curve_edge.map(EdgeId::from_index);
            let id = breps
                .insert(brep)
                .map_err(|e| HarnessError::new(e.to_string()))?;
            Ok(NativeOperand {
                id,
                face,
                vertices,
                curve_edge,
                metadata,
            })
        }
        let [first, second] =
            ffi::static_fixture(input).map_err(|e| HarnessError::new(e.to_string()))?;
        let mut breps = Arena::default();
        let first = operand(&mut breps, first)?;
        let second = operand(&mut breps, second)?;
        Ok(Self {
            breps,
            operands: [first, second],
        })
    }

    fn operand(&self, slot: usize) -> SolidResult<&NativeOperand> {
        self.operands.get(slot).ok_or_else(static_error)
    }

    fn brep(&self, slot: usize) -> SolidResult<&Brep> {
        self.breps
            .get(self.operand(slot)?.id)
            .ok_or_else(static_error)
    }
}

struct StaticScene<'a> {
    points: Vec<Point3>,
    trees: Option<&'a [Tree<StaticHandle>; 4]>,
    curve: Curve,
    native: Option<NativeScene>,
    geometry_failure: Cell<Option<(usize, bool, bool)>>,
}
fn static_error() -> SolidError {
    SolidError::invalid_input(
        "difftest topology p2 static adapter",
        "object absent from the native static fixture",
    )
}
fn static_tree(
    points: &[Point3],
    offset: usize,
    tolerance: f64,
) -> SolidResult<Tree<StaticHandle>> {
    let mut bbox = Extent3d::new();
    for point in points {
        let mut b = Extent3d::from_point(*point);
        b.expand_absolute(tolerance);
        bbox = bbox.union(&b);
    }
    bbox.expand_absolute(tolerance * 100.0);
    let size = bbox.size();
    bbox.add_point(bbox.max() + size * 0.15);
    bbox.add_point(bbox.min() - size * 0.15);
    let mut tree = Tree::new_spatial(&bbox);
    for (index, point) in points.iter().enumerate() {
        let mut b = Extent3d::from_point(*point);
        b.expand_absolute(tolerance);
        tree.add_to_spatial_tree(
            &b,
            Some(StaticHandle::Fixture((offset + index) as i32)),
            NodeGeomType::Default,
        )?;
    }
    Ok(tree)
}
impl<'a> StaticScene<'a> {
    fn tree_view(&self, brep: StaticHandle) -> SolidResult<TopologyTrees<'a, StaticHandle>> {
        // These trees describe only the shell-vertex input of entries 0–4. The
        // face entries use explicit local solvers and never request B-rep trees.
        let [first, second, empty_first, empty_second] = self.trees.ok_or_else(static_error)?;
        let (vertices, empty) = match brep.fixture_id()? {
            1000 => (first, empty_first),
            1001 => (second, empty_second),
            _ => return Err(static_error()),
        };
        Ok(TopologyTrees {
            vertices,
            curves: empty,
            surfaces: empty,
        })
    }

    fn native(&self) -> SolidResult<&NativeScene> {
        self.native.as_ref().ok_or_else(static_error)
    }

    fn face_scope(&self, face: StaticHandle) -> SolidResult<(usize, FaceId)> {
        let (slot, id) = match face {
            StaticHandle::Fixture(2000) => (0, self.native()?.operand(0)?.face),
            StaticHandle::Fixture(2001) => (1, self.native()?.operand(1)?.face),
            StaticHandle::Member(slot, CacheObject::Face(id)) => (slot, id),
            _ => return Err(static_error()),
        };
        self.native()?.brep(slot)?.face(id)?;
        Ok((slot, id))
    }

    fn surface_scope(&self, surface: StaticHandle) -> SolidResult<(usize, FaceId)> {
        match surface {
            StaticHandle::Fixture(2002) => self.face_scope(StaticHandle::Fixture(2000)),
            StaticHandle::Fixture(2003) => self.face_scope(StaticHandle::Fixture(2001)),
            StaticHandle::Member(slot, CacheObject::Surface(id)) => {
                self.face_scope(StaticHandle::Member(slot, CacheObject::Face(id)))
            }
            _ => Err(static_error()),
        }
    }

    fn edge_scope(&self, edge: StaticHandle) -> SolidResult<(usize, EdgeId)> {
        let (slot, id) = match edge {
            StaticHandle::Fixture(2004) => (0, self.query_edge()?.ok_or_else(static_error)?),
            StaticHandle::Member(slot, CacheObject::Edge(id)) => (slot, id),
            _ => return Err(static_error()),
        };
        self.native()?.brep(slot)?.edge(id)?;
        Ok((slot, id))
    }

    fn curve_scope(&self, curve: StaticHandle) -> SolidResult<Option<(usize, EdgeId)>> {
        match curve {
            StaticHandle::Fixture(999) => self
                .query_edge()?
                .map(|id| self.edge_scope(StaticHandle::Member(0, CacheObject::Edge(id))))
                .transpose(),
            StaticHandle::Member(slot, CacheObject::Curve(id)) => self
                .edge_scope(StaticHandle::Member(slot, CacheObject::Edge(id)))
                .map(Some),
            _ => Err(static_error()),
        }
    }

    fn vertex_scope(&self, vertex: StaticHandle) -> SolidResult<(usize, VertexId)> {
        if let StaticHandle::Member(slot, CacheObject::Vertex(id)) = vertex {
            self.native()?.brep(slot)?.vertex(id)?;
            return Ok((slot, id));
        }
        let mut index = usize::try_from(vertex.fixture_id()?).map_err(|_| static_error())?;
        for (slot, operand) in self.native()?.operands.iter().enumerate() {
            if let Some(id) = operand.vertices.get(index) {
                self.native()?.brep(slot)?.vertex(*id)?;
                return Ok((slot, *id));
            }
            index = index
                .checked_sub(operand.vertices.len())
                .ok_or_else(static_error)?;
        }
        Err(static_error())
    }

    fn query_edge(&self) -> SolidResult<Option<EdgeId>> {
        self.native
            .as_ref()
            .map_or(Ok(None), |native| Ok(native.operand(0)?.curve_edge))
    }

    /// Prefer the native protocol's listed identity when this is that same
    /// record. Unlisted boundary topology retains its typed, operand-local id.
    fn member_handle(&self, slot: usize, object: CacheObject) -> SolidResult<StaticHandle> {
        let native = self.native()?;
        let operand = native.operand(slot)?;
        let brep = native.brep(slot)?;
        match object {
            CacheObject::Face(id) | CacheObject::Surface(id) => {
                let face = brep.face(id)?;
                if matches!(object, CacheObject::Surface(_)) {
                    face.surface.as_ref().ok_or_else(static_error)?;
                }
                if id == operand.face {
                    let slot = i32::try_from(slot).map_err(|_| static_error())?;
                    let base = if matches!(object, CacheObject::Surface(_)) {
                        2002
                    } else {
                        2000
                    };
                    return Ok(StaticHandle::Fixture(base + slot));
                }
            }
            CacheObject::Edge(id) | CacheObject::Curve(id) => {
                let edge = brep.edge(id)?;
                if matches!(object, CacheObject::Curve(_)) {
                    edge.curve.as_ref().ok_or_else(static_error)?;
                }
                if slot == 0 && Some(id) == operand.curve_edge {
                    return Ok(StaticHandle::Fixture(
                        if matches!(object, CacheObject::Curve(_)) {
                            999
                        } else {
                            2004
                        },
                    ));
                }
            }
            CacheObject::Vertex(id) => {
                brep.vertex(id)?;
                if let Some(index) = operand
                    .vertices
                    .iter()
                    .position(|candidate| *candidate == id)
                {
                    let offset: usize = native
                        .operands
                        .iter()
                        .take(slot)
                        .map(|operand| operand.vertices.len())
                        .sum();
                    let id = i32::try_from(offset + index).map_err(|_| static_error())?;
                    return Ok(StaticHandle::Fixture(id));
                }
            }
            _ => return Err(static_error()),
        }
        Ok(StaticHandle::Member(slot, object))
    }

    fn check_geometry(&self, slot: usize, curves: bool, surfaces: bool) -> SolidResult<()> {
        if self
            .native()?
            .operand(slot)?
            .metadata
            .deferred_output_for(curves, surfaces)
            .is_some()
        {
            self.geometry_failure.set(Some((slot, curves, surfaces)));
            return Err(SolidError::unported(
                "native topology solver fixture geometry export",
                "harness:analytic-geometry-export",
            ));
        }
        Ok(())
    }

    fn live_face_mut(&mut self, face: StaticHandle) -> SolidResult<(&mut Brep, FaceId)> {
        let (slot, face) = self.face_scope(face)?;
        self.check_geometry(slot, true, true)?;
        let native = self.native.as_mut().ok_or_else(static_error)?;
        let id = native.operand(slot)?.id;
        Ok((native.breps.get_mut(id).ok_or_else(static_error)?, face))
    }

    fn topology_limitation(&self) -> Option<Output> {
        self.native
            .as_ref()?
            .operands
            .iter()
            .find_map(|operand| operand.metadata.deferred_output_for(false, false))
    }

    fn geometry_limitation(&self) -> Option<Output> {
        let (slot, curves, surfaces) = self.geometry_failure.get()?;
        self.native
            .as_ref()?
            .operands
            .get(slot)?
            .metadata
            .deferred_output_for(curves, surfaces)
    }
}

impl TopologyObjects<StaticHandle> for StaticScene<'_> {
    fn with_brep_trees(
        &mut self,
        _kernel: &Kernel,
        first: StaticHandle,
        second: Option<StaticHandle>,
        solve: &mut BrepTreeSolve<'_, StaticHandle>,
    ) -> SolidResult<SolutionArray<StaticHandle>> {
        // External immutable tree storage survives the independent mutable face
        // borrow. No tree or B-rep is copied, rebuilt or unpublished by the solve.
        let first = self.tree_view(first)?;
        let second = second.map(|object| self.tree_view(object)).transpose()?;
        solve(self, first, second)
    }

    fn brep_trees(
        &self,
        _: &Kernel,
        brep: StaticHandle,
    ) -> SolidResult<TopologyTrees<'_, StaticHandle>> {
        self.tree_view(brep)
    }

    fn object_zone_tol_3d(
        &self,
        kernel: &Kernel,
        object: Option<StaticHandle>,
    ) -> SolidResult<f64> {
        let system = || {
            solidcore::geom::tol::zone_tol_3d(
                kernel,
                solidcore::geom::tol::SmallTopoType::NormalSize,
            )
        };
        let Some(object) = object else {
            return Ok(system());
        };
        if object == StaticHandle::Fixture(999) && self.query_edge()?.is_none() {
            self.curve(object)?;
            // The external query curve has no owner in the native constructor.
            return Ok(system());
        }
        if self.native.is_none() {
            return match object {
                StaticHandle::Fixture(1000 | 1001) => Ok(1e-5),
                _ => self.vertex_tolerance(object),
            };
        }
        let (slot, source) = match object {
            StaticHandle::Member(slot, source) => {
                self.member_handle(slot, source)?;
                match source {
                    CacheObject::Curve(_) => self.check_geometry(slot, true, false)?,
                    CacheObject::Surface(_) => self.check_geometry(slot, false, true)?,
                    _ => {}
                }
                (slot, source)
            }
            StaticHandle::Fixture(999 | 2004) => {
                self.check_geometry(0, true, false)?;
                let edge = self.query_edge()?.ok_or_else(static_error)?;
                (
                    0,
                    if object == StaticHandle::Fixture(999) {
                        CacheObject::Curve(edge)
                    } else {
                        CacheObject::Edge(edge)
                    },
                )
            }
            StaticHandle::Fixture(1000) => (0, CacheObject::Brep),
            StaticHandle::Fixture(1001) => (1, CacheObject::Brep),
            StaticHandle::Fixture(2000 | 2001) => {
                let (slot, face) = self.face_scope(object)?;
                (slot, CacheObject::Face(face))
            }
            StaticHandle::Fixture(2002 | 2003) => {
                let (slot, face) = self.surface_scope(object)?;
                self.check_geometry(slot, false, true)?;
                (slot, CacheObject::Surface(face))
            }
            StaticHandle::Fixture(_) => {
                let (slot, vertex) = self.vertex_scope(object)?;
                (slot, CacheObject::Vertex(vertex))
            }
        };
        brep_cache::object_zone_tol_3d(
            kernel,
            self.native()?.brep(slot)?,
            Some(source),
            &arena_identity,
            &|_| None,
        )
    }

    fn face_point_classify(
        &mut self,
        kernel: &Kernel,
        face: StaticHandle,
        uv: Vec2,
        source_zone_tolerance: f64,
        do_boundaries: bool,
        use_uv: bool,
        classification: &mut PointClassification<EntityRef>,
    ) -> SolidResult<()> {
        let (brep, face) = self.live_face_mut(face)?;
        solidcore::query::face::p2::point_classify(
            kernel,
            brep,
            face,
            uv,
            classification,
            solidcore::query::face::p2::PointClassifyOptions {
                src_zone_tol_3d: source_zone_tolerance,
                do_boundary_intersections: do_boundaries,
                use_uv_classification: use_uv,
                ..solidcore::query::face::p2::PointClassifyOptions::default()
            },
        )
    }

    fn face_point_3d_classify(
        &mut self,
        kernel: &Kernel,
        face: StaticHandle,
        point: Point3,
        source_zone_tolerance: f64,
        do_boundaries: bool,
        classification: &mut PointClassification<EntityRef>,
    ) -> SolidResult<()> {
        let (_, face_id) = self.face_scope(face)?;
        let Some(uv) = solidcore::intersect::topology_solver::p1::face_point_3d_uv_with_context(
            kernel,
            self,
            face,
            face_id,
            point,
            source_zone_tolerance,
            classification,
            None,
        )?
        else {
            return Ok(());
        };
        let (brep, face_id) = self.live_face_mut(face)?;
        solidcore::query::face::p1::point_3d_classify_finish(
            kernel,
            brep,
            face_id,
            uv,
            point,
            source_zone_tolerance,
            do_boundaries,
            classification,
            solidcore::query::face::p1::Point3dClassifyOptions::default(),
        )
    }

    fn face_brep_editing_enabled(&self, face: StaticHandle) -> SolidResult<bool> {
        let (slot, _) = self.face_scope(face)?;
        Ok(self.native()?.brep(slot)?.editing_enabled())
    }

    fn set_face_brep_editing_enabled(
        &mut self,
        face: StaticHandle,
        editing: bool,
    ) -> SolidResult<()> {
        let (slot, _) = self.face_scope(face)?;
        let native = self.native.as_mut().ok_or_else(static_error)?;
        let id = native.operand(slot)?.id;
        native
            .breps
            .get_mut(id)
            .ok_or_else(static_error)?
            .set_editing_enabled(editing);
        Ok(())
    }

    fn face_create_uv_trim_curves(
        &mut self,
        kernel: &Kernel,
        face: StaticHandle,
    ) -> SolidResult<()> {
        let (brep, face_id) = self.live_face_mut(face)?;
        solidcore::topology::face::p3::create_uv_trim_curves(kernel, brep, face_id, true, None)
            .map(|_| ())
    }

    fn face_classification_object(
        &self,
        face: StaticHandle,
        entity: EntityRef,
    ) -> SolidResult<StaticHandle> {
        let (slot, _) = self.face_scope(face)?;
        let object = match entity {
            EntityRef::Vertex(id) => CacheObject::Vertex(id),
            EntityRef::Edge(id) => CacheObject::Edge(id),
            EntityRef::Face(id) => CacheObject::Face(id),
            _ => return Err(static_error()),
        };
        self.member_handle(slot, object)
    }

    fn curve(&self, object: StaticHandle) -> SolidResult<&Curve> {
        if let Some((slot, edge)) = self.curve_scope(object)? {
            self.check_geometry(slot, true, false)?;
            self.native()?
                .brep(slot)?
                .edge(edge)?
                .curve
                .as_ref()
                .ok_or_else(static_error)
        } else {
            Ok(&self.curve)
        }
    }
    fn analytic_curve_handle(&self, curve: StaticHandle) -> Option<StaticHandle> {
        match curve {
            StaticHandle::Fixture(_) => Some(StaticHandle::Fixture(-2)),
            StaticHandle::Member(_, CacheObject::Curve(_)) => Some(curve),
            StaticHandle::Member(..) => None,
        }
    }
    fn vertex_point(&self, vertex: StaticHandle) -> SolidResult<Point3> {
        if let Some(native) = &self.native {
            let (slot, vertex) = self.vertex_scope(vertex)?;
            return Ok(native.brep(slot)?.vertex(vertex)?.point);
        }
        usize::try_from(vertex.fixture_id()?)
            .ok()
            .and_then(|i| self.points.get(i))
            .copied()
            .ok_or_else(static_error)
    }
    fn vertex_tolerance(&self, vertex: StaticHandle) -> SolidResult<f64> {
        if let Some(native) = &self.native {
            let (slot, vertex) = self.vertex_scope(vertex)?;
            return Ok(native.brep(slot)?.vertex(vertex)?.tolerance);
        }
        self.vertex_point(vertex)?;
        Ok(1e-5)
    }
    fn curve_edge(&self, curve: StaticHandle) -> SolidResult<Option<StaticHandle>> {
        let Some(owner) = self.curve(curve)?.edge(&|_| None)? else {
            return Ok(None);
        };
        let (slot, _) = self.curve_scope(curve)?.ok_or_else(static_error)?;
        let edge = self
            .native()?
            .brep(slot)?
            .edges()
            .ids()
            .find(|id| arena_identity(CacheObject::Edge(*id)) == owner)
            .ok_or_else(static_error)?;
        self.member_handle(slot, CacheObject::Edge(edge)).map(Some)
    }
    fn edge_curve(&self, edge: StaticHandle) -> SolidResult<Option<StaticHandle>> {
        let (slot, id) = self.edge_scope(edge)?;
        self.check_geometry(slot, true, false)?;
        self.native()?
            .brep(slot)?
            .edge(id)?
            .curve
            .as_ref()
            .map(|_| self.member_handle(slot, CacheObject::Curve(id)))
            .transpose()
    }
    fn edge_interval(&self, edge: StaticHandle) -> SolidResult<Extent1d> {
        let (slot, id) = self.edge_scope(edge)?;
        Ok(self.native()?.brep(slot)?.edge(id)?.interval)
    }
    fn surface_face(&self, surface: StaticHandle) -> SolidResult<Option<StaticHandle>> {
        let Some(owner) = self.surface(surface)?.face()? else {
            return Ok(None);
        };
        let (slot, _) = self.surface_scope(surface)?;
        let face = self
            .native()?
            .brep(slot)?
            .faces()
            .ids()
            .find(|id| arena_identity(CacheObject::Face(*id)) == owner)
            .ok_or_else(static_error)?;
        self.member_handle(slot, CacheObject::Face(face)).map(Some)
    }
    fn face_surface(&self, face: StaticHandle) -> SolidResult<StaticHandle> {
        let (slot, id) = self.face_scope(face)?;
        self.member_handle(slot, CacheObject::Surface(id))
    }
    fn surface(&self, surface: StaticHandle) -> SolidResult<&Surface> {
        let (slot, face) = self.surface_scope(surface)?;
        self.check_geometry(slot, false, true)?;
        self.native()?
            .brep(slot)?
            .face(face)?
            .surface
            .as_ref()
            .ok_or_else(static_error)
    }
    fn face_uv_domain(&self, face: StaticHandle) -> SolidResult<Extent2d> {
        let (slot, face) = self.face_scope(face)?;
        Ok(self.native()?.brep(slot)?.face(face)?.uv_domain)
    }
    fn face_tolerance(&self, face: StaticHandle) -> SolidResult<f64> {
        let (slot, face) = self.face_scope(face)?;
        Ok(self.native()?.brep(slot)?.face(face)?.tolerance)
    }
    fn surface_owner(&self, surface: StaticHandle) -> SolidResult<Option<StaticHandle>> {
        let Some(owner) = self.surface(surface)?.owner() else {
            return Ok(None);
        };
        let (slot, _) = self.surface_scope(surface)?;
        let face = self
            .native()?
            .brep(slot)?
            .faces()
            .ids()
            .find(|id| arena_identity(CacheObject::Face(*id)) == owner)
            .ok_or_else(static_error)?;
        self.member_handle(slot, CacheObject::Face(face)).map(Some)
    }
    fn shape_member(&self, object: StaticHandle) -> ShapeMember {
        let Some(native) = &self.native else {
            return ShapeMember::Other;
        };
        let scoped = match object {
            StaticHandle::Fixture(2000 | 2001) | StaticHandle::Member(_, CacheObject::Face(_)) => {
                self.face_scope(object)
                    .ok()
                    .map(|(slot, id)| (slot, CacheObject::Face(id)))
            }
            StaticHandle::Fixture(2004) | StaticHandle::Member(_, CacheObject::Edge(_)) => self
                .edge_scope(object)
                .ok()
                .map(|(slot, id)| (slot, CacheObject::Edge(id))),
            StaticHandle::Fixture(0..=127) | StaticHandle::Member(_, CacheObject::Vertex(_)) => {
                self.vertex_scope(object)
                    .ok()
                    .map(|(slot, id)| (slot, CacheObject::Vertex(id)))
            }
            _ => None,
        };
        scoped
            .and_then(|(slot, object)| {
                let brep = native.operands.get(slot)?.id;
                match object {
                    CacheObject::Face(id) => Some(ShapeMember::Face((brep, id))),
                    CacheObject::Edge(id) => Some(ShapeMember::Edge((brep, id))),
                    CacheObject::Vertex(id) => Some(ShapeMember::Vertex((brep, id))),
                    _ => None,
                }
            })
            .unwrap_or(ShapeMember::Other)
    }
}

fn static_rust(case: &Case) -> Result<Output, HarnessError> {
    let input = static_input(case)?;
    // These bounds match the native protocol. Face inputs additionally use its
    // validated graph-only export; no native solver result is read by this side.
    let i: StaticInputs = case.parse()?;
    if input.len() > 396
        || i.first.is_empty()
        || i.second.is_empty()
        || i.first.len() > 64
        || i.second.len() > 64
    {
        return Err(HarnessError::new("unbounded static fixture"));
    }
    let kernel = Kernel::default();
    let first: Vec<Point3> = i
        .first
        .iter()
        .map(|p| Point3::new(p[0], p[1], p[2]))
        .collect();
    let second: Vec<Point3> = i
        .second
        .iter()
        .map(|p| Point3::new(p[0], p[1], p[2]))
        .collect();
    let mut all = first.clone();
    all.extend(&second);
    let trees = if i.entry < 5 {
        let tree1 = static_tree(&first, 0, 1e-5).map_err(|e| HarnessError::new(e.to_string()))?;
        let tree2 = static_tree(&second, first.len(), 1e-5)
            .map_err(|e| HarnessError::new(e.to_string()))?;
        let empty = |tree: &Tree<StaticHandle>| -> Result<Tree<StaticHandle>, HarnessError> {
            let root = tree
                .top_node()
                .ok_or_else(|| HarnessError::new("missing fixture root"))?;
            let bbox = &tree
                .node(root)
                .map_err(|e| HarnessError::new(e.to_string()))?
                .bbox;
            Ok(Tree::new_spatial(bbox))
        };
        let empty_first = empty(&tree1)?;
        let empty_second = empty(&tree2)?;
        Some([tree1, tree2, empty_first, empty_second])
    } else {
        None
    };
    let nurb = NlCurve::new(
        vec![
            CPoint::new(-1.0, 0.0, 0.0, NOW),
            CPoint::new(1.0, 0.0, 0.0, NOW),
        ],
        1,
        vec![0.0, 0.0, 1.0, 1.0],
    )
    .map_err(|e| HarnessError::new(e.to_string()))?;
    let mut scene = StaticScene {
        points: all,
        trees: trees.as_ref(),
        curve: Curve::from(BSplineCurve::from_nurb(nurb, 3, false)),
        native: if i.entry >= 5 {
            Some(NativeScene::import(&input)?)
        } else {
            None
        },
        geometry_failure: Cell::new(None),
    };
    if let Some(deferred) = scene.topology_limitation() {
        return Ok(deferred);
    }
    let vectors = [Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)];
    let options = TopologySolveOptions {
        operation: SolverOperation::from_code(i.operation)
            .ok_or_else(|| HarnessError::new("operation"))?,
        requested: SolutionRequested::from_code(i.request)
            .ok_or_else(|| HarnessError::new("request"))?,
        tolerance_3d: i.tolerance,
        best_answer_so_far_sq: static_best(i.operation),
        opt_vectors: &vectors,
    };
    let point = Point3::new(i.point[0], i.point[1], i.point[2]);
    let curve = SolverCurve::Object(StaticHandle::Fixture(999));
    let interval = Extent1d::from_min_max(0.0, 1.0);
    let shape = Shape::empty().map_err(|e| HarnessError::new(e.to_string()))?;
    let result = match i.entry {
        0 => TopologySolver::brep_point_solve(
            &kernel,
            &mut scene,
            Some(StaticHandle::Fixture(1000)),
            point,
            options,
        ),
        1 => TopologySolver::brep_curve_solve(
            &kernel,
            &mut scene,
            Some(StaticHandle::Fixture(1000)),
            curve,
            &interval,
            options,
        ),
        2 => TopologySolver::brep_brep_solve(
            &kernel,
            &mut scene,
            Some(StaticHandle::Fixture(1000)),
            Some(StaticHandle::Fixture(1001)),
            options,
        ),
        3 => TopologySolver::shape_point_solve(&kernel, &mut scene, &shape, point, options),
        4 => TopologySolver::shape_curve_solve(
            &kernel, &mut scene, &shape, curve, &interval, options,
        ),
        5..=8 => {
            use solidcore::intersect::topology_solver::p1::{LocalSolver, LocalSolverKind};
            let mut ts = TopologySolver::new(
                &kernel,
                &mut scene,
                options.operation,
                options.requested,
                options.tolerance_3d,
                options.best_answer_so_far_sq,
                options.opt_vectors,
            );
            let kind = match i.entry {
                5 => LocalSolverKind::VertexFace,
                6 => LocalSolverKind::CurveFace,
                7 => LocalSolverKind::EdgeFace,
                _ => LocalSolverKind::FaceFace,
            };
            let mut local =
                LocalSolver::new(kind, i.swap, if i.included { None } else { Some(&shape) });
            if i.entry == 6 {
                local.set_curve_object(Some(StaticHandle::Fixture(999)), interval);
            }
            let (first, second) = match i.entry {
                5 => (StaticHandle::Fixture(0), StaticHandle::Fixture(2002)),
                6 | 7 => (StaticHandle::Fixture(999), StaticHandle::Fixture(2002)),
                _ => (StaticHandle::Fixture(2002), StaticHandle::Fixture(2003)),
            };
            let (first, second) = if i.swap {
                (second, first)
            } else {
                (first, second)
            };
            local
                .solve_it(&mut ts, Some(first), Some(second))
                .map(|()| ts.base.solutions)
        }
        _ => return Err(HarnessError::new("static entry")),
    };
    let mut flat = vec![0.0];
    for slot in 0..2 {
        let tolerance = scene
            .native
            .as_ref()
            .map_or(Ok(1e-5), |native| native.brep(slot).map(Brep::tolerance))
            .map_err(|e| HarnessError::new(e.to_string()))?;
        flat.push(tolerance);
    }
    for index in 0..scene.points.len() {
        flat.push(
            scene
                .vertex_tolerance(StaticHandle::Fixture(index as i32))
                .map_err(|e| HarnessError::new(e.to_string()))?,
        );
    }
    match result {
        Ok(solutions) => {
            flat.push(solutions.len() as f64);
            for s in solutions.as_slice() {
                flat.extend([
                    f64::from(s.solution_type.code()),
                    s.num_variables as f64,
                    s.num_objects as f64,
                    s.objects[0].map_or(-1.0, static_object_id),
                    s.objects[1].map_or(-1.0, static_object_id),
                    s.start.solution_value,
                ]);
                flat.extend(s.start.parameters.iter().take(s.num_variables).copied());
                if s.solution_type == SolutionType::RangeOfValues {
                    flat.push(s.end.solution_value);
                    flat.extend(s.end.parameters.iter().take(s.num_variables).copied());
                }
                flat.extend(
                    s.nodes
                        .iter()
                        .map(|node| f64::from(u8::from(node.is_some()))),
                );
            }
        }
        Err(e) => {
            if let Some(deferred) = Output::unported(&e) {
                return Ok(scene.geometry_limitation().unwrap_or(deferred));
            }
            let status = match e {
                SolidError::InvalidInput {
                    operation: "difftest topology p2 static adapter",
                    ..
                } => {
                    return Err(HarnessError::new(e.to_string()));
                }
                SolidError::InvalidInput { .. } => 1007.0,
                SolidError::Upstream { code, .. } => f64::from(code),
                _ => return Err(HarnessError::new(e.to_string())),
            };
            if let Some(slot) = flat.first_mut() {
                *slot = status;
            }
            flat.push(0.0);
        }
    }
    Ok(static_output(&flat))
}

#[cfg(test)]
mod context_tests {
    use super::*;

    #[test]
    fn classification_handles_keep_unlisted_topology_in_its_live_operand() {
        fn operand(
            breps: &mut Arena<BrepId, Brep>,
            tolerance: f64,
        ) -> (NativeOperand, EdgeId, VertexId) {
            let mut brep = Brep::new();
            let shell = brep.add_shell().unwrap();
            let face = brep.add_face(shell).unwrap();
            let listed = brep.add_vertex(Point3::new(0.0, 0.0, 0.0)).unwrap();
            let boundary = brep.add_vertex(Point3::new(1.0, 0.0, 0.0)).unwrap();
            let edge = brep.add_edge(listed, boundary).unwrap();
            brep.vertex_mut(boundary).unwrap().tolerance = tolerance;
            brep.edge_mut(edge).unwrap().tolerance = tolerance * 2.0;
            let id = breps.insert(brep).unwrap();
            (
                NativeOperand {
                    id,
                    face,
                    vertices: vec![listed],
                    curve_edge: None,
                    metadata: QueryMetadata::default(),
                },
                edge,
                boundary,
            )
        }

        let mut breps = Arena::default();
        let (first, edge0, vertex0) = operand(&mut breps, 0.125);
        let (second, edge1, vertex1) = operand(&mut breps, 0.5);
        // Both arena-local ids coincide; they denote distinct live objects.
        assert_eq!(edge0, edge1);
        assert_eq!(vertex0, vertex1);
        let second_id = second.id;
        let curve = Curve::from(BSplineCurve::from_nurb(
            NlCurve::new(
                vec![
                    CPoint::new(0.0, 0.0, 0.0, NOW),
                    CPoint::new(1.0, 0.0, 0.0, NOW),
                ],
                1,
                vec![0.0, 0.0, 1.0, 1.0],
            )
            .unwrap(),
            3,
            false,
        ));
        let scene = StaticScene {
            points: Vec::new(),
            trees: None,
            curve,
            native: Some(NativeScene {
                breps,
                operands: [first, second],
            }),
            geometry_failure: Cell::new(None),
        };
        let kernel = Kernel::default();
        let face0 = StaticHandle::Fixture(2000);
        let face1 = StaticHandle::Fixture(2001);
        let boundary0 = scene
            .face_classification_object(face0, EntityRef::Edge(edge0))
            .unwrap();
        let boundary1 = scene
            .face_classification_object(face1, EntityRef::Edge(edge1))
            .unwrap();
        assert_ne!(boundary0, boundary1);
        assert_eq!(scene.object_zone_tol_3d(&kernel, Some(boundary0)), Ok(0.25));
        assert_eq!(scene.object_zone_tol_3d(&kernel, Some(boundary1)), Ok(1.0));
        assert_eq!(
            scene.shape_member(boundary1),
            ShapeMember::Edge((second_id, edge1))
        );
        assert_eq!(static_object_id(boundary0), -2.0);

        let vertex = scene
            .face_classification_object(face1, EntityRef::Vertex(vertex1))
            .unwrap();
        assert_eq!(scene.vertex_tolerance(vertex), Ok(0.5));
        assert_eq!(static_object_id(vertex), -2.0);
        let listed = scene.native().unwrap().operand(1).unwrap().vertices[0];
        let listed_handle = scene
            .face_classification_object(face1, EntityRef::Vertex(listed))
            .unwrap();
        assert_eq!(static_object_id(listed_handle), 1.0);

        let missing = VertexId::from_index(200);
        assert!(matches!(
            scene.face_classification_object(face1, EntityRef::Vertex(missing)),
            Err(SolidError::NotFound { .. }),
        ));
    }
}
