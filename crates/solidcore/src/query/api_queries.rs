//! Public closest-point, boundary, bounding-box and classification queries.
//!
//! Upstream: `source/SM_API/src/SmApiQueries.cpp`, `source/SM_API/inc/SmApiQueries.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`curve_closest_point`], [`curve_closest_point_all`] | `SmApiCurveClosestPoint`, `SmApiCurveClosestPointAll` |
//! | [`surface_closest_point`] | `SmApiSurfaceClosestPoint` |
//! | [`get_closest_point`], [`brep_closest_point`], [`brep_closest_point_all`] | `SmApiGetClosestPoint`, `SmApiBrepClosestPoint`, `SmApiBrepClosestPointAll` |
//! | [`brep_relationship`], [`brep_distance`] | `SmApiBrepRelationship`, `SmApiBrepDistance` |
//! | [`get_edges`], [`get_faces`], [`find_edge`], [`find_faces`] | `SmApiGetEdges`, `SmApiGetFaces`, `SmApiFindEdge`, `SmApiFindFaces` |
//! | [`face_get_boundary_loops`] | `SmApiFaceGetBoundaryLoops` |
//! | [`brep_is_manifold_solid`], [`brep_bounding_box`] | `SmApiBrepIsManifoldSolid`, `SmApiBrepBoundingBox` |
//! | [`face_bounding_box`], [`face_internal_point`], [`face_classify_uv`], [`face_closest_point`] | `SmApiFaceBoundingBox`, `SmApiFaceInternalPoint`, `SmApiFaceClassifyUV`, `SmApiFaceClosestPoint` |
//! | [`edge_bounding_box`], [`edge_closest_point`], [`edge_parameter_range`], [`edge_classify_parameter`], [`edge_tangent`] | `SmApiEdgeBoundingBox`, `SmApiEdgeClosestPoint`, `SmApiEdgeParameterRange`, `SmApiEdgeClassifyParameter`, `SmApiEdgeTangent` |
//! | [`vertex_get_point`] | `SmApiVertexGetPoint` |
//! | `count_contained_samples`, `point_in_solid` (private) | `sm_CountContainedSamples`, `sm_PointInSolid` |
//!
//! Deviations: typed references/ids replace nullable pointers; results replace
//! output parameters (including append-only arrays), so errors expose no partial
//! arrays. Callers supply the Kernel and topology-solver object adapter/cache;
//! handles retain the caller's identity. Missing topology queries return Unported
//! at the call site. No new recursion or allocation from serialized counts.

use crate::curve::curve::p1::GlobalPointSolveOptions;
use crate::curve::curve::Curve;
use crate::error::{ignored, Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::extent3d::Extent3d;
use crate::geom::global_solver::{SolutionRequested, SolverOperation};
use crate::geom::object::{FaceKey, ObjectRef, Shape, ShapeMember};
use crate::geom::solution_array::{ObjectKind, Solution, SolutionArray, SolutionObject};
use crate::geom::tools::{self, EdgePick};
use crate::intersect::point_class::PointClassType;
use crate::intersect::topology_solver::p1::{TopologyObjects, TopologySolver};
use crate::intersect::topology_solver::p2::TopologySolveOptions;
use crate::kernel::Kernel;
use crate::math::{Point2, Point3, Vec3};
use crate::nurbs::b_spline_curve::BoundingBoxRequest;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::surface::surface::p1::zone_tol_3d;
use crate::surface::surface::p4::NoObject;
use crate::surface::surface::p5::DropPointOptions;
use crate::surface::surface::Surface;
use crate::topology::brep_cache::CacheObject;
use crate::topology::{
    edge, edgeuse::p1 as eu, face::p2 as face, Brep, BrepId, EdgeId, FaceId, VertexId,
};

const OP: &str = "query::api_queries";

/// A closest spatial point and its upstream (unsquared) solution distance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClosestPoint {
    /// Point on the queried entity.
    pub point: Point3,
    /// Distance supplied by the solver.
    pub distance: f64,
}

/// A closest curve point with its parameter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveClosestPoint {
    /// Spatial result.
    pub closest: ClosestPoint,
    /// Curve parameter.
    pub parameter: f64,
}

/// One spatial point and all returned UV representations of that point.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceClosestPoint {
    /// Spatial result.
    pub closest: ClosestPoint,
    /// UVs in solver order, including seam representations.
    pub uvs: Vec<Point2>,
}

#[derive(Clone, Copy, PartialEq)]
struct CurveHandle(u8, Option<Extent1d>);
impl SolutionObject for CurveHandle {
    fn kind(&self) -> ObjectKind {
        ObjectKind::Curve
    }
    fn natural_interval(&self) -> Option<Extent1d> {
        self.1
    }
}

fn curve_solutions(
    k: &Kernel,
    curve: &Curve,
    point: Point3,
    requested: SolutionRequested,
) -> Result<SolutionArray<CurveHandle>> {
    curve.global_point_solve(
        k,
        &curve.natural_interval(),
        SolverOperation::Minimize,
        point,
        crate::geom::tol::ZONE_TOL_3D,
        requested,
        &GlobalPointSolveOptions {
            handle: Some(CurveHandle(0, Some(curve.natural_interval()))),
            analytic_handle: Some(CurveHandle(1, None)),
            ..GlobalPointSolveOptions::default()
        },
    )
}

fn curve_result<O: Copy>(curve: &Curve, solution: &Solution<O>) -> Result<CurveClosestPoint> {
    let parameter = solution.start.param(0)?;
    let point = curve.evaluate_point(parameter)?;
    Ok(CurveClosestPoint {
        closest: ClosestPoint {
            point,
            distance: solution.start.solution_value,
        },
        parameter,
    })
}

/// Closest point on the curve's natural interval (`SmApiCurveClosestPoint`).
/// Solver and evaluation errors propagate; an empty solution set returns `SM_ERR`.
pub fn curve_closest_point(k: &Kernel, curve: &Curve, point: Point3) -> Result<CurveClosestPoint> {
    let solutions = curve_solutions(k, curve, point, SolutionRequested::Single)?;
    curve_result(
        curve,
        solutions.as_slice().first().ok_or(SolidError::sm_err(OP))?,
    )
}

/// All closest curve solutions, retaining solver order (`SmApiCurveClosestPointAll`).
/// An empty solution set or an evaluation failure returns an error.
pub fn curve_closest_point_all(
    k: &Kernel,
    curve: &Curve,
    point: Point3,
) -> Result<Vec<CurveClosestPoint>> {
    let solutions = curve_solutions(k, curve, point, SolutionRequested::All)?;
    if solutions.is_empty() {
        return Err(SolidError::sm_err(OP));
    }
    solutions
        .as_slice()
        .iter()
        .map(|s| curve_result(curve, s))
        .collect()
}

/// Closest surface point and UV representations (`SmApiSurfaceClosestPoint`).
/// Strictly smaller distance replaces the selected solution; failed secondary
/// evaluations are ignored. UVs must evaluate within the surface zone tolerance.
pub fn surface_closest_point(
    k: &Kernel,
    surface: &Surface,
    point: Point3,
) -> Result<SurfaceClosestPoint> {
    let drop = surface.drop_point_with_kernel::<NoObject>(
        k,
        point,
        &surface.natural_uv_domain(),
        &DropPointOptions {
            want_solutions: true,
            ..DropPointOptions::default()
        },
    )?;
    if !drop.success || drop.solutions.is_empty() {
        return Err(SolidError::sm_err(OP));
    }
    let tolerance = zone_tol_3d(k, surface)?;
    let mut best = drop
        .solutions
        .as_slice()
        .first()
        .ok_or(SolidError::sm_err(OP))?;
    for solution in drop.solutions.as_slice().iter().skip(1) {
        if solution.start.solution_value < best.start.solution_value {
            best = solution;
        }
    }
    let selected = Point2::new(best.start.param(0)?, best.start.param(1)?);
    let point = surface.evaluate_point(selected)?;
    let mut uvs = Vec::new();
    for solution in drop.solutions.as_slice() {
        let uv = Point2::new(solution.start.param(0)?, solution.start.param(1)?);
        if let Some(candidate) = ignored(surface.evaluate_point(uv))? {
            if candidate.closer_than(tolerance, point) {
                uvs.push(uv);
            }
        }
    }
    if uvs.is_empty() {
        return Err(SolidError::sm_err(OP));
    }
    Ok(SurfaceClosestPoint {
        closest: ClosestPoint {
            point,
            distance: best.start.solution_value,
        },
        uvs,
    })
}

fn solve_options(_k: &Kernel, requested: SolutionRequested) -> TopologySolveOptions<'static> {
    TopologySolveOptions {
        operation: SolverOperation::Minimize,
        requested,
        tolerance_3d: crate::geom::tol::ZONE_TOL_3D,
        best_answer_so_far_sq: 1.0e20 * 1.0e20,
        opt_vectors: &[],
    }
}

// SmSolution::GetPoint using the same caller-owned geometry adapter as the solver.
fn solution_point<O: SolutionObject>(
    objects: &dyn TopologyObjects<O>,
    solution: &Solution<O>,
    index: usize,
) -> Result<Point3> {
    let mut parameter = 0;
    for (i, object) in solution
        .objects
        .iter()
        .take(solution.num_objects)
        .enumerate()
    {
        let object = object.ok_or(SolidError::null_pointer(OP))?;
        let kind = object.kind();
        if i == index {
            if kind == ObjectKind::Vertex {
                return objects.vertex_point(object);
            }
            if matches!(kind, ObjectKind::Curve | ObjectKind::Edge) {
                let curve = if kind == ObjectKind::Edge {
                    objects
                        .edge_curve(object)?
                        .ok_or(SolidError::null_pointer(OP))?
                } else {
                    object
                };
                return objects
                    .curve(curve)?
                    .evaluate_point(solution.start.param(parameter)?);
            }
            if matches!(kind, ObjectKind::Surface | ObjectKind::Face) {
                let surface = if kind == ObjectKind::Face {
                    objects.face_surface(object)?
                } else {
                    object
                };
                return objects.surface(surface)?.evaluate_point(Point2::new(
                    solution.start.param(parameter)?,
                    solution.start.param(parameter + 1)?,
                ));
            }
            return Err(SolidError::sm_err(OP));
        }
        parameter += match kind {
            ObjectKind::Curve | ObjectKind::Edge => 1,
            ObjectKind::Surface | ObjectKind::Face => 2,
            _ => 0,
        };
    }
    Err(SolidError::sm_err(OP))
}

/// Closest point on a B-rep (`SmApiGetClosestPoint`); the adapter supplies its cache.
/// Solver errors, point-extraction failures and empty solution sets propagate.
pub fn get_closest_point<O: SolutionObject>(
    k: &Kernel,
    objects: &dyn TopologyObjects<O>,
    brep: O,
    point: Point3,
) -> Result<ClosestPoint> {
    brep_closest_point(k, objects, brep, point)
}

/// Single closest B-rep point (`SmApiBrepClosestPoint`).
/// Errors are those of `BrepPointSolve` and solution point extraction.
pub fn brep_closest_point<O: SolutionObject>(
    k: &Kernel,
    objects: &dyn TopologyObjects<O>,
    brep: O,
    point: Point3,
) -> Result<ClosestPoint> {
    let solutions = TopologySolver::brep_point_solve(
        k,
        objects,
        Some(brep),
        point,
        solve_options(k, SolutionRequested::Single),
    )?;
    let solution = solutions.as_slice().first().ok_or(SolidError::sm_err(OP))?;
    Ok(ClosestPoint {
        point: solution_point(objects, solution, 0)?,
        distance: solution.start.solution_value,
    })
}

/// All closest B-rep points in solver order (`SmApiBrepClosestPointAll`).
/// An empty solution set is `SM_ERR`; extraction failures propagate.
pub fn brep_closest_point_all<O: SolutionObject>(
    k: &Kernel,
    objects: &dyn TopologyObjects<O>,
    brep: O,
    point: Point3,
) -> Result<Vec<ClosestPoint>> {
    let solutions = TopologySolver::brep_point_solve(
        k,
        objects,
        Some(brep),
        point,
        solve_options(k, SolutionRequested::All),
    )?;
    if solutions.is_empty() {
        return Err(SolidError::sm_err(OP));
    }
    solutions
        .as_slice()
        .iter()
        .map(|s| {
            Ok(ClosestPoint {
                point: solution_point(objects, s, 0)?,
                distance: s.start.solution_value,
            })
        })
        .collect()
}

/// B-rep arena and its handle in the caller's topology-solver adapter.
#[derive(Clone, Copy)]
pub struct BrepRef<'a, O> {
    /// Exact topology.
    pub brep: &'a Brep,
    /// Adapter handle.
    pub handle: O,
}

/// Relationship code and boundary distance (`SmApiBrepRelationship`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepRelationship {
    /// 0 separate, 1 touching, 2 interpenetrating, 3 A contains B, 4 B contains A.
    pub relationship: i32,
    /// Boundary gap, or zero when within the sum of the B-rep tolerances.
    pub distance: f64,
}

// Returns the point class and whether its region exists and is material.
fn classify_point(
    _target: &Brep,
    _point: Point3,
    _tolerance: f64,
) -> Result<(PointClassType, bool)> {
    // UNPORTED(query.brep): SmBrep::Point3DClassify
    Err(SolidError::unported(
        "SmBrep::Point3DClassify",
        "query.brep",
    ))
}

fn point_in_solid(target: &Brep, point: Point3, tolerance: f64) -> Result<bool> {
    Ok(ignored(classify_point(target, point, tolerance))?
        .is_some_and(|(class, material)| class == PointClassType::Region && material))
}

fn brep_vertices(_brep: &Brep) -> Result<Vec<VertexId>> {
    // UNPORTED(topology.brep_p1): SmBrep::GetVertices
    Err(SolidError::unported(
        "SmBrep::GetVertices",
        "topology.brep_p1",
    ))
}

fn count_contained_samples(
    k: &Kernel,
    sample: &Brep,
    target: &Brep,
    tolerance: f64,
) -> Result<(usize, usize)> {
    let mut points = Vec::new();
    for vertex in brep_vertices(sample)? {
        points.push(sample.vertex(vertex)?.point);
    }
    for face in brep_faces(sample)? {
        if let Some(point) = ignored(face_internal_point(k, sample, face))? {
            points.push(point);
        }
    }
    let (mut inside, mut outside) = (0, 0);
    for point in points {
        if let Some((class, material)) = ignored(classify_point(target, point, tolerance))? {
            match class {
                PointClassType::Region => {
                    if material {
                        inside += 1;
                    } else {
                        outside += 1;
                    }
                }
                PointClassType::Unknown => outside += 1,
                _ => {}
            }
        }
    }
    Ok((inside, outside))
}

/// Boundary distance followed by containment/touching classification.
/// Missing point-classification/list dependencies return `Unported` rather than
/// being swallowed by upstream's ignored ordinary classification statuses.
pub fn brep_relationship<O: SolutionObject>(
    k: &Kernel,
    objects: &dyn TopologyObjects<O>,
    a: BrepRef<'_, O>,
    b: BrepRef<'_, O>,
) -> Result<BrepRelationship> {
    let solutions = TopologySolver::brep_brep_solve(
        k,
        objects,
        Some(a.handle),
        Some(b.handle),
        solve_options(k, SolutionRequested::Single),
    )?;
    let solution = solutions.as_slice().first().ok_or(SolidError::sm_err(OP))?;
    let distance = solution.start.solution_value;
    let tolerance =
        crate::geom::tol::context_zone_tol_3d(k, None, crate::geom::tol::SmallTopoType::NormalSize);
    if distance > a.brep.tolerance() + b.brep.tolerance() {
        let pa = ignored(solution_point(objects, solution, 0))?;
        let pb = if pa.is_some() {
            ignored(solution_point(objects, solution, 1))?
        } else {
            None
        };
        let relationship = if let (Some(pa), Some(pb)) = (pa, pb) {
            let bin_a = point_in_solid(a.brep, pb, tolerance)?;
            let ain_b = point_in_solid(b.brep, pa, tolerance)?;
            if bin_a {
                3
            } else if ain_b {
                4
            } else {
                0
            }
        } else {
            let (bin, bout) = count_contained_samples(k, b.brep, a.brep, tolerance)?;
            let (ain, aout) = count_contained_samples(k, a.brep, b.brep, tolerance)?;
            if bin > 0 && bout == 0 {
                3
            } else if ain > 0 && aout == 0 {
                4
            } else {
                0
            }
        };
        Ok(BrepRelationship {
            relationship,
            distance,
        })
    } else {
        let (bin, _) = count_contained_samples(k, b.brep, a.brep, tolerance)?;
        let (ain, _) = count_contained_samples(k, a.brep, b.brep, tolerance)?;
        Ok(BrepRelationship {
            relationship: if bin > 0 || ain > 0 { 2 } else { 1 },
            distance: 0.0,
        })
    }
}

/// Boundary distance through the same relationship query (`SmApiBrepDistance`).
/// Containment/classification errors propagate, even though only distance is returned.
pub fn brep_distance<O: SolutionObject>(
    k: &Kernel,
    objects: &dyn TopologyObjects<O>,
    a: BrepRef<'_, O>,
    b: BrepRef<'_, O>,
) -> Result<f64> {
    Ok(brep_relationship(k, objects, a, b)?.distance)
}

fn brep_edges(_brep: &Brep) -> Result<Vec<EdgeId>> {
    // UNPORTED(topology.brep_p1): SmBrep::GetEdges
    Err(SolidError::unported("SmBrep::GetEdges", "topology.brep_p1"))
}
fn brep_faces(_brep: &Brep) -> Result<Vec<FaceId>> {
    // UNPORTED(topology.brep_p1): SmBrep::GetFaces
    Err(SolidError::unported("SmBrep::GetFaces", "topology.brep_p1"))
}

/// Exact B-rep edge list (`SmApiGetEdges`). Empty lists return `SM_ERR`.
/// Missing list traversal is `Unported`.
pub fn get_edges(brep: &Brep) -> Result<Vec<EdgeId>> {
    let edges = brep_edges(brep)?;
    if edges.is_empty() {
        return Err(SolidError::sm_err(OP));
    }
    Ok(edges)
}

/// Exact B-rep face list (`SmApiGetFaces`). Empty lists return `SM_ERR`.
/// Missing list traversal is `Unported`.
pub fn get_faces(brep: &Brep) -> Result<Vec<FaceId>> {
    let faces = brep_faces(brep)?;
    if faces.is_empty() {
        return Err(SolidError::sm_err(OP));
    }
    Ok(faces)
}

/// Flattened edges and orientations in upward face-use loop order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryLoops {
    /// Number of edge uses in each loop (including zero for a vertex loop).
    pub edge_counts: Vec<usize>,
    /// Edge of each edge use, retaining repeated edges.
    pub edges: Vec<EdgeId>,
    /// Same or opposite orientation of each edge use.
    pub orientations: Vec<OrientType>,
}

/// Upward boundary loops (`SmApiFaceGetBoundaryLoops`).
/// Missing face use/loops and invalid edge/orientation return `SM_ERR`.
pub fn face_get_boundary_loops(brep: &Brep, id: FaceId) -> Result<BoundaryLoops> {
    let mut result = BoundaryLoops {
        edge_counts: Vec::new(),
        edges: Vec::new(),
        orientations: Vec::new(),
    };
    let loops = if let Some(fu) = face::upward_faceuse(brep, id)? {
        eu::faceuse_loopuses(brep, fu, false, None)?
    } else {
        Vec::new()
    };
    if loops.is_empty() {
        return Err(SolidError::sm_err(OP));
    }
    for lu in loops {
        let edges = eu::loopuse_edgeuses(brep, lu)?;
        result.edge_counts.push(edges.len());
        for edgeuse in edges {
            let edgeuse = brep.edgeuse(edgeuse)?;
            let edge = edgeuse.edge.ok_or(SolidError::sm_err(OP))?;
            if !matches!(edgeuse.orientation, OrientType::Same | OrientType::Opposite) {
                return Err(SolidError::sm_err(OP));
            }
            result.edges.push(edge);
            result.orientations.push(edgeuse.orientation);
        }
    }
    Ok(result)
}

/// Finds an edge index and parameter (`SmApiFindEdge`); missing match is `SM_ERR`.
pub fn find_edge(brep: &Brep, point: Point3) -> Result<EdgePick> {
    tools::find_edge(brep, point)?.ok_or(SolidError::sm_err(OP))
}

/// Finds a face (`SmApiFindFaces`); missing match is `SM_ERR`.
pub fn find_faces(
    k: &Kernel,
    brep: &Brep,
    id: BrepId,
    point: Point3,
    keys: &dyn Fn(CacheObject) -> ObjectRef,
    owners: &dyn Fn(ObjectRef) -> Option<ObjectRef>,
) -> Result<FaceId> {
    tools::find_face(k, brep, id, point, keys, owners)?.ok_or(SolidError::sm_err(OP))
}

/// Manifold solid predicate (`SmApiBrepIsManifoldSolid`); waits for B-rep traversal.
pub fn brep_is_manifold_solid(_brep: &Brep) -> Result<bool> {
    // UNPORTED(topology.brep_p2): SmBrep::IsManifoldSolid
    Err(SolidError::unported(
        "SmBrep::IsManifoldSolid",
        "topology.brep_p2",
    ))
}

/// Loose or tight B-rep box (`SmApiBrepBoundingBox`); missing queries are `Unported`.
pub fn brep_bounding_box(_k: &Kernel, _brep: &Brep, tight: bool) -> Result<Extent3d> {
    // UNPORTED(query.brep): SmBrep::CalculateBoundingBox / CalculateTightBoundingBox
    Err(SolidError::unported(
        if tight {
            "SmBrep::CalculateTightBoundingBox"
        } else {
            "SmBrep::CalculateBoundingBox"
        },
        "query.brep",
    ))
}

/// Loose or tight face box (`SmApiFaceBoundingBox`); missing queries are `Unported`.
pub fn face_bounding_box(_k: &Kernel, brep: &Brep, face: FaceId, tight: bool) -> Result<Extent3d> {
    brep.face(face)?;
    // UNPORTED(query.face_p1): SmFace::CalculateBoundingBox
    // UNPORTED(query.face_p3): SmFace::CalculateTightBoundingBox
    Err(SolidError::unported(
        if tight {
            "SmFace::CalculateTightBoundingBox"
        } else {
            "SmFace::CalculateBoundingBox"
        },
        if tight {
            "query.face_p3"
        } else {
            "query.face_p1"
        },
    ))
}

/// Interior face point (`SmApiFaceInternalPoint`). Errors from the face query
/// propagate, including explicit deferrals of its unfinished dependencies.
pub fn face_internal_point(k: &Kernel, brep: &Brep, face: FaceId) -> Result<Point3> {
    face::calculate_an_internal_point(k, brep, face)
}

/// UV classification (`SmApiFaceClassifyUV`), using UV containment rather than a
/// 3D boundary snap. The missing face classifier is `Unported`.
pub fn face_classify_uv(
    _k: &Kernel,
    brep: &Brep,
    face: FaceId,
    _uv: Point2,
) -> Result<PointClassType> {
    brep.face(face)?;
    // UNPORTED(query.face_p2): SmFace::PointClassify
    Err(SolidError::unported(
        "SmFace::PointClassify",
        "query.face_p2",
    ))
}

/// Trim-aware closest point on one face (`SmApiFaceClosestPoint`).
/// Shape construction and solver errors propagate.
pub fn face_closest_point<O: SolutionObject>(
    k: &Kernel,
    objects: &dyn TopologyObjects<O>,
    face: FaceKey,
    point: Point3,
) -> Result<ClosestPoint> {
    let shape = Shape::new(&[], &[face], &[], &[])?;
    let mut options = solve_options(k, SolutionRequested::Single);
    options.tolerance_3d =
        objects.face_tolerance(objects.shape_member_object(ShapeMember::Face(face))?)?;
    let solutions = TopologySolver::shape_point_solve(k, objects, &shape, point, options)?;
    let solution = solutions.as_slice().first().ok_or(SolidError::sm_err(OP))?;
    Ok(ClosestPoint {
        point: solution_point(objects, solution, 0)?,
        distance: solution.start.solution_value,
    })
}

/// Axis-aligned edge box (`SmApiEdgeBoundingBox`), preserving the distinction
/// between the direct tight-box call and loose-box calculation.
pub fn edge_bounding_box(k: &Kernel, brep: &Brep, edge: EdgeId, tight: bool) -> Result<Extent3d> {
    if tight {
        return Ok(
            crate::query::edge::calculate_tight_bounding_box(k, brep, edge)?.unwrap_or_default(),
        );
    }
    crate::query::edge::calculate_bounding_box(
        k,
        brep,
        edge,
        BoundingBoxRequest {
            normal: true,
            pseudo: false,
            polar: false,
        },
        false,
    )?
    .normal
    .ok_or(SolidError::sm_err(OP))
}

/// Closest point within an edge interval (`SmApiEdgeClosestPoint`).
/// A failed drop or absent curve returns `SM_ERR`; evaluation errors propagate.
pub fn edge_closest_point(
    k: &Kernel,
    brep: &Brep,
    id: EdgeId,
    point: Point3,
) -> Result<CurveClosestPoint> {
    let result =
        crate::query::edge::closest_point(k, brep, id, point)?.ok_or(SolidError::sm_err(OP))?;
    let curve = brep
        .edge(id)?
        .curve
        .as_ref()
        .ok_or(SolidError::sm_err(OP))?;
    Ok(CurveClosestPoint {
        closest: ClosestPoint {
            point: curve.evaluate_point(result.parameter)?,
            distance: result.distance,
        },
        parameter: result.parameter,
    })
}

/// Parameter interval (`SmApiEdgeParameterRange`); dangling ids are errors.
pub fn edge_parameter_range(brep: &Brep, id: EdgeId) -> Result<Extent1d> {
    Ok(brep.edge(id)?.interval)
}

/// Classifies a parameter (`SmApiEdgeClassifyParameter`). Vertices are tested
/// first with strict distance comparisons, including parameters outside the edge
/// interval. An ordinary evaluation failure falls through to interval containment.
pub fn edge_classify_parameter(
    _k: &Kernel,
    brep: &Brep,
    id: EdgeId,
    parameter: f64,
) -> Result<PointClassType> {
    let e = brep.edge(id)?;
    let curve = e.curve.as_ref().ok_or(SolidError::sm_err(OP))?;
    if let Some(point) = ignored(curve.evaluate_point(parameter))? {
        for v in [edge::start_vertex(brep, id)?, edge::end_vertex(brep, id)?]
            .into_iter()
            .flatten()
        {
            let tolerance = crate::topology::vertex::tolerance(brep, v)?;
            if point.distance(brep.vertex(v)?.point) < tolerance {
                return Ok(PointClassType::Vertex);
            }
        }
    }
    Ok(if e.interval.contains_value(parameter, 0.0) {
        PointClassType::Edge
    } else {
        PointClassType::Unknown
    })
}

/// Unit tangent (`SmApiEdgeTangent`), evaluating from the left as upstream.
/// A missing curve, evaluation failure or zero derivative returns an error.
pub fn edge_tangent(brep: &Brep, id: EdgeId, parameter: f64) -> Result<Vec3> {
    let curve = brep
        .edge(id)?
        .curve
        .as_ref()
        .ok_or(SolidError::sm_err(OP))?;
    let values = curve.evaluate(parameter, 1, true, true)?;
    values
        .get(1)
        .copied()
        .ok_or(SolidError::sm_err(OP))?
        .unitize()
}

/// Vertex position (`SmApiVertexGetPoint`); dangling ids are errors.
pub fn vertex_get_point(brep: &Brep, id: VertexId) -> Result<Point3> {
    Ok(brep.vertex(id)?.point)
}
