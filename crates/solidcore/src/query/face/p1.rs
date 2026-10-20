//! Face, shell and region bounding boxes, face point and edge-crossing classification,
//! and the checked curve evaluation API entry point.
//!
//! Upstream: `source/SMLib/src/SmFace.cpp`, `source/SMLib/src/SmShell.cpp`,
//! `source/SMLib/src/SmRegion.cpp`, `source/SM_API/src/SmApiCurves.cpp`,
//! `source/SMLib/inc/SmFace.h`, `source/SMLib/inc/SmCurveClass.h`.
//!
//! | Rust                                  | Upstream                         |
//! |---------------------------------------|----------------------------------|
//! | [`evaluate_curve`]                    | `SmApiEvaluateCurve`             |
//! | [`calculate_bounding_box_into`], [`calculate_bounding_box`], [`calculate_loose_bounding_box`] | `SmFace::CalculateBoundingBox` |
//! | [`calculate_shell_bounding_box_into`], [`calculate_shell_bounding_box`] | `SmShell::CalculateBoundingBox` |
//! | [`calculate_region_bounding_box_into`], [`calculate_region_bounding_box`] | `SmRegion::CalculateBoundingBox` |
//! | [`point_3d_classify`], [`Point3dClassifyOptions`] | `SmFace::Point3DClassify` |
//! | [`edge_curve_on_classify`], [`EdgeCurveOnOptions`], [`CurveOnClassification`] | `SmFace::EdgeCurveOnClassify` |
//! | [`ClassifyProblem`]                   | `SmClassifyProblem` (also re-exported by `query.face_p3`) |
//!
//! Deviations: external tolerance-owner walks are iterative and bounded to 64 links
//! (as `Curve::edge`), refusing cycles instead of unbounded native traversal.
//! A nullable curve is `Option<&Curve>`; `SmApiEvaluateCurve`'s three
//! optional output pointers become one returned `[position, d1, d2]` (all three are
//! evaluated upstream too). The face box's in/out form retains the boundary box
//! published before fallible surface work in the tight path; the loose path writes
//! only on success. The shell and region boxes initialise `output` first and keep
//! the union of the children that succeeded before a failing child, as upstream's
//! in/out box does; [`calculate_shell_bounding_box`] / [`calculate_region_bounding_box`]
//! return the box instead. Box queries that expose the tight option take `&mut Brep`
//! because the tight face calculation can create UV trim curves during classification;
//! [`calculate_loose_bounding_box`] keeps the read-only face query. Shell and region
//! boxes take `&mut Kernel` because the wire edge traversal claims a context mark. The surface box of a face expands by the
//! face's own tolerance (what `SmBSplineSurface::CalculateBoundingBox` reads through
//! `GetFace()`): arena surfaces hold no owning face, so the B-spline/plane/revolution
//! box is computed unexpanded and expanded here, as the override's last step does.
//! `Point3DClassify` requests the surface cache before its local drop, allowing
//! ordinary cache failure as upstream does, and reuses a successful checkout for
//! the global drop with [`FacePointTest::Disabled`] (SmFace.cpp:12796-12801).
//! This expresses the native `m_bPointTestEnabled` guard without changing shared
//! cache state or the surface owner. The local drop does not read that flag.
//! The subsequent explicit UV `PointClassify` still runs on the original B-rep;
//! this entry point supplies it with no trimmed surface cache.
//! The drop and classification phases are shared with providers that can resolve
//! a surface's actual owner and mutably access that owner's B-rep. The direct
//! [`point_3d_classify`] entry has no object-key mapping; a cache needing such a
//! mapping reports the unresolved context instead of treating an owned surface
//! as a standalone copy.
//! `EdgeCurveOnClassify` lifts a 2D curve to a curve on the
//! face's surface only when an upward edge use of the face is found (upstream lifts it
//! first and ignores `SetAll`'s status); an `SmCrvOnSurf` over a surface other than a
//! B-spline is a dispatch gap (`curve.crv_on_surf`). No recursion.

use crate::curve::crv_on_surf::CrvOnSurf;
use crate::curve::curve::{Curve, CurveKind};
use crate::error::ignored;
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::extent3d::Extent3d;
use crate::geom::global_solver::{SolutionRequested, SolverOperation};
use crate::geom::object::ObjectRef;
use crate::geom::solution_array::SolutionObject;
use crate::geom::tol::xsect_tol_3d_from_zones;
use crate::geom::vector2d::Point2;
use crate::geom::vector3d::{Point3, Vec3};
use crate::intersect::point_class::{ClassGeometry, PointClassType, PointClassification};
use crate::kernel::Kernel;
use crate::nurbs::b_spline_curve::BoundingBoxRequest;
use crate::nurbs::b_spline_surface::SurfParamType;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::query::face::p3;
use crate::surface::surface::p1::get_surface_cache_with_context;
use crate::surface::surface::p4::{FacePointTest, NoObject};
use crate::surface::surface::{Surface, SurfaceKind};
use crate::surface::surface_cache::p1::{OpaqueSurfaceCacheContext, SurfaceCacheContext};
use crate::surface::volume::Volume;
use crate::topology::brep_cache::{object_zone_tol_3d, CacheObject};
use crate::topology::face::p2 as fp2;
use crate::topology::shell::{self, ShellType};
use crate::topology::{
    edge, edgeuse::p1 as eu1, topology_traverser, vertex, Brep, EdgeId, EdgeuseId, FaceId, LoopId,
    RegionId, ShellId,
};
use crate::{EntityRef, Result, SolidError};

const OP_CURVE: &str = "query::face::p1::evaluate_curve";
const OP_FACE_BOX: &str = "query::face::p1::calculate_bounding_box";
const OP_SHELL_BOX: &str = "query::face::p1::calculate_shell_bounding_box";
const OP_REGION_BOX: &str = "query::face::p1::calculate_region_bounding_box";
const OP_POINT: &str = "query::face::p1::point_3d_classify";
const OP_EDGE_CURVE: &str = "query::face::p1::edge_curve_on_classify";

// ---- SmApiEvaluateCurve ------------------------------------------------------------------

/// `SmApiEvaluateCurve`: position, first and second derivative of `curve` at
/// `parameter`, which must lie in the curve's natural interval (inclusive). Evaluates
/// from the left with nonzero tangents (upstream's defaults).
///
/// # Errors
///
/// [`SolidError::InvalidInput`] (`SM_ERR_INVALID_INPUT`) for a null curve or a parameter
/// outside the natural interval (NaN included); those of the curve's evaluator.
pub fn evaluate_curve(curve: Option<&Curve>, parameter: f64) -> Result<[Vec3; 3]> {
    let curve = curve.ok_or_else(|| SolidError::invalid_input(OP_CURVE, "null curve"))?;

    // Negated in-range test so a NaN parameter is rejected too.
    let domain = curve.natural_interval();
    if !(parameter >= domain.min() && parameter <= domain.max()) {
        return Err(SolidError::invalid_input(
            OP_CURVE,
            "parameter outside natural domain",
        ));
    }

    let values = curve.evaluate(parameter, 2, true, true)?;
    let mut iter = values.into_iter();
    let missing = || SolidError::invalid_input(OP_CURVE, "curve returned too few derivatives");
    Ok([
        iter.next().ok_or_else(missing)?,
        iter.next().ok_or_else(missing)?,
        iter.next().ok_or_else(missing)?,
    ])
}

// ---- bounding boxes ------------------------------------------------------------------------

fn tight_face_box(k: &Kernel, b: &mut Brep, f: FaceId, output: &mut Extent3d) -> Result<()> {
    p3::calculate_tight_bounding_box_into(k, b, f, output)
}

const NORMAL_BOX: BoundingBoxRequest = BoundingBoxRequest {
    normal: true,
    pseudo: false,
    polar: false,
};

/// `SmFace::CalculateBoundingBox(rFaceBBox, bTight)` (`tight` defaults to `false`
/// upstream): with `tight`, `CalculateTightBoundingBox`. Otherwise, when the face's
/// surface is linear in U or V, the union of its vertex points (each expanded by the
/// vertex tolerance) and its edge curve boxes over the edge intervals (each expanded by
/// the edge tolerance); else the surface's box over the face's UV domain, expanded by
/// the face tolerance.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for a face without surface (upstream dereferences null)
/// or a vertex use without vertex; those of the curve and surface boxes and of the
/// face's vertex and edge lists; those of [`p3::calculate_tight_bounding_box`] for
/// `tight`.
pub fn calculate_bounding_box(
    k: &Kernel,
    b: &mut Brep,
    f: FaceId,
    tight: bool,
) -> Result<Extent3d> {
    let mut output = Extent3d::new();
    calculate_bounding_box_into(k, b, f, tight, &mut output)?;
    Ok(output)
}

/// [`calculate_bounding_box`] with the native in/out box contract.
///
/// Loose calculations replace `output` only on success. Tight calculations retain
/// their partial box as documented by [`p3::calculate_tight_bounding_box_into`].
///
/// # Errors
///
/// Those of [`calculate_bounding_box`].
pub fn calculate_bounding_box_into(
    k: &Kernel,
    b: &mut Brep,
    f: FaceId,
    tight: bool,
    output: &mut Extent3d,
) -> Result<()> {
    // when asked - call expensive TightBoundingBox
    if tight {
        return tight_face_box(k, b, f, output);
    }

    *output = calculate_loose_bounding_box(k, b, f)?;
    Ok(())
}

/// The read-only, non-tight `SmFace::CalculateBoundingBox` path.
///
/// The returned box includes the native vertex, edge and surface tolerance
/// expansions. No point classification or trim-curve creation is needed.
///
/// # Errors
///
/// Those of the face's surface, vertex and edge lists and their box calculations.
pub fn calculate_loose_bounding_box(k: &Kernel, b: &Brep, f: FaceId) -> Result<Extent3d> {
    let mut face_box = Extent3d::new();

    let surface = fp2::surface(b, f)?
        .ok_or_else(|| SolidError::invalid_input(OP_FACE_BOX, "face without surface"))?;
    let u_deg = surface.degree(SurfParamType::U);
    let v_deg = surface.degree(SurfParamType::V);

    // If one of the directions of the surface is linear we can compute the bounding box
    // using the edges and vertices.
    if u_deg == 1 || v_deg == 1 {
        for v in super::p2::face_vertices(b, f, false, None)? {
            let mut v_box = Extent3d::from_point(vertex::point(b, v)?);
            v_box.expand_absolute(vertex::tolerance(b, v)?);
            face_box = face_box.union(&v_box);
        }

        // Now load the Edge Tree
        for e in topology_traverser::face_edges(b, f)? {
            if let Some(curve) = edge::curve(b, e)? {
                let mut e_box = curve
                    .calculate_bounding_box(k, &edge::interval(b, e)?, NORMAL_BOX, true)?
                    .normal
                    .ok_or_else(|| SolidError::invalid_input(OP_FACE_BOX, "missing edge box"))?;
                e_box.expand_absolute(edge::tolerance(b, e)?);
                face_box = face_box.union(&e_box);
            }
        }
    } else {
        // Get the bounding box of the surface
        face_box = surface_box(k, surface, &fp2::uv_domain(b, f)?, fp2::tolerance(b, f)?)?;
    }

    // Set output and return
    Ok(face_box)
}

/// `pSurface->CalculateBoundingBox(GetUVDomain(), &box)` of a face's surface (expanded by
/// the zone tolerance, the default): the B-spline override (also that of planes and
/// surfaces of revolution) expands its normal box by the face tolerance as its last
/// step; the other kinds dispatch as is.
fn surface_box(
    k: &Kernel,
    surface: &Surface,
    uv_domain: &Extent2d,
    face_tolerance: f64,
) -> Result<Extent3d> {
    let missing = || SolidError::invalid_input(OP_FACE_BOX, "missing surface box");
    match surface.kind() {
        SurfaceKind::BSpline(_) | SurfaceKind::Plane(_) | SurfaceKind::SurfOfRevolution(_) => {
            let mut normal = surface
                .calculate_bounding_box(k, uv_domain, [true, false, false], false)?
                .normal_box
                .ok_or_else(missing)?;
            normal.expand_absolute(face_tolerance);
            Ok(normal)
        }
        _ => surface
            .calculate_bounding_box(k, uv_domain, [true, false, false], true)?
            .normal_box
            .ok_or_else(missing),
    }
}

/// `SmShell::CalculateBoundingBox(rBBox, bTight)` (`tight` defaults to `false`): `output`
/// is initialised, then holds a vertex shell's vertex point, the union of a wire
/// shell's edge boxes ([`crate::query::edge::calculate_bounding_box`] over
/// [`shell::wire_edges`]) or the union of a face shell's face boxes
/// ([`calculate_bounding_box`] of each face use's face, in face-use list order); it
/// stays empty for any other shell type.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for a vertex shell without vertex or a face use without
/// face (upstream dereferences null); those of the child boxes and lists. `output` then
/// holds the union of the children before the failing one.
pub fn calculate_shell_bounding_box_into(
    k: &mut Kernel,
    b: &mut Brep,
    id: ShellId,
    tight: bool,
    output: &mut Extent3d,
) -> Result<()> {
    output.init();

    let s = b.shell(id)?;
    match s.shell_type {
        ShellType::Vertex => {
            let v = s.vertex(b)?.ok_or_else(|| {
                SolidError::invalid_input(OP_SHELL_BOX, "vertex shell without vertex")
            })?;
            output.add_point(vertex::point(b, v)?);
        }
        ShellType::Wire => {
            for e in shell::wire_edges(k, b, id)? {
                let this_box =
                    crate::query::edge::calculate_bounding_box(k, b, e, NORMAL_BOX, tight)?
                        .normal
                        .ok_or_else(|| {
                            SolidError::invalid_input(OP_SHELL_BOX, "missing edge box")
                        })?;
                *output = output.union(&this_box);
            }
        }
        ShellType::Faceuse => {
            for fu in s.get_faceuses()? {
                let f = b.faceuse(fu)?.face.ok_or_else(|| {
                    SolidError::invalid_input(OP_SHELL_BOX, "face use without face")
                })?;
                let this_box = calculate_bounding_box(k, b, f, tight)?;
                *output = output.union(&this_box);
            }
        }
        ShellType::Unknown | ShellType::Other(_) => {}
    }

    Ok(())
}

/// [`calculate_shell_bounding_box_into`], returning the box.
///
/// # Errors
///
/// Those of [`calculate_shell_bounding_box_into`].
pub fn calculate_shell_bounding_box(
    k: &mut Kernel,
    b: &mut Brep,
    id: ShellId,
    tight: bool,
) -> Result<Extent3d> {
    let mut output = Extent3d::new();
    calculate_shell_bounding_box_into(k, b, id, tight, &mut output)?;
    Ok(output)
}

/// `SmRegion::CalculateBoundingBox(rBBox, bTight)` (`tight` defaults to `false`): `output`
/// is initialised, then holds the union of the region's shell boxes
/// ([`calculate_shell_bounding_box_into`]) in shell list order.
///
/// # Errors
///
/// `SM_ERR` ([`SolidError::Upstream`]) for a region without shells; those of the shell
/// boxes. `output` then holds the union of the shells before the failing one.
pub fn calculate_region_bounding_box_into(
    k: &mut Kernel,
    b: &mut Brep,
    id: RegionId,
    tight: bool,
    output: &mut Extent3d,
) -> Result<()> {
    output.init();

    let shells = b.region(id)?.get_shells(None, |_, _| Ok(false))?;
    if shells.is_empty() {
        return Err(SolidError::sm_err(OP_REGION_BOX));
    }
    let mut this_box = Extent3d::new();
    for s in shells {
        calculate_shell_bounding_box_into(k, b, s, tight, &mut this_box)?;
        *output = output.union(&this_box);
    }

    Ok(())
}

/// [`calculate_region_bounding_box_into`], returning the box.
///
/// # Errors
///
/// Those of [`calculate_region_bounding_box_into`].
pub fn calculate_region_bounding_box(
    k: &mut Kernel,
    b: &mut Brep,
    id: RegionId,
    tight: bool,
) -> Result<Extent3d> {
    let mut output = Extent3d::new();
    calculate_region_bounding_box_into(k, b, id, tight, &mut output)?;
    Ok(output)
}

// ---- SmFace::Point3DClassify ---------------------------------------------------------------

/// The optional arguments of [`point_3d_classify`]; [`Default`] gives upstream's
/// defaults (`pOptGuessUV = NULL`, `bTgtLoopsOnly = FALSE`, `pOptTgtLoops = NULL`,
/// `bUVSpaceOkay = TRUE`).
#[derive(Debug, Clone, Copy)]
pub struct Point3dClassifyOptions<'a> {
    /// `pOptGuessUV`: a guess for a local point drop before the global one.
    pub guess_uv: Option<Point2>,
    /// `bTgtLoopsOnly`: only classify against [`Point3dClassifyOptions::tgt_loops`].
    pub tgt_loops_only: bool,
    /// `pOptTgtLoops`: with `tgt_loops_only`, the loops to classify against (`None` or
    /// empty: the outer loop).
    pub tgt_loops: Option<&'a [LoopId]>,
    /// `bUVSpaceOkay`: UV-space classification may be tried when the 3D one is doubtful.
    pub uv_space_okay: bool,
}

impl Default for Point3dClassifyOptions<'_> {
    fn default() -> Self {
        Self {
            guess_uv: None,
            tgt_loops_only: false,
            tgt_loops: None,
            uv_space_okay: true,
        }
    }
}

/// The geometry and stored topology tolerances of one B-rep, for a classification of
/// its entities.
struct BrepGeometry<'a> {
    b: &'a Brep,
}

impl ClassGeometry<EntityRef> for BrepGeometry<'_> {
    fn curve(&self, o: EntityRef) -> Option<&Curve> {
        match o {
            EntityRef::Edge(e) => edge::curve(self.b, e).ok().flatten(),
            _ => None,
        }
    }

    fn surface(&self, o: EntityRef) -> Option<&Surface> {
        match o {
            EntityRef::Face(f) => fp2::surface(self.b, f).ok().flatten(),
            _ => None,
        }
    }

    fn volume(&self, _o: EntityRef) -> Option<&Volume> {
        None
    }

    fn point(&self, o: EntityRef) -> Option<Point3> {
        match o {
            EntityRef::Vertex(v) => vertex::point(self.b, v).ok(),
            _ => None,
        }
    }

    fn classification_vertex_point(&self, v: EntityRef) -> Option<Point3> {
        self.point(v)
    }

    fn topology_zone_tol_3d(&self, o: EntityRef) -> Option<f64> {
        match o {
            EntityRef::Vertex(v) => vertex::tolerance(self.b, v).ok(),
            EntityRef::Edge(e) => edge::tolerance(self.b, e).ok(),
            EntityRef::Face(f) => fp2::tolerance(self.b, f).ok(),
            _ => None,
        }
    }
}

/// The UV start of the first solution of a point solve.
fn solution_uv(parameters: &[f64]) -> Result<Point2> {
    let missing = || SolidError::invalid_input(OP_POINT, "point solve without parameters");
    Ok(Point2::new(
        *parameters.first().ok_or_else(missing)?,
        *parameters.get(1).ok_or_else(missing)?,
    ))
}

/// `SmFace::Point3DClassify`: classifies the 3D point `point` against face `f` into `pc`:
/// the point is dropped onto the face's surface within `src_zone_tol_3d` (a local drop
/// from [`Point3dClassifyOptions::guess_uv`] first, else the global drop over the face's
/// UV domain); a drop outside the UV domain (by more than a thousandth of its diagonal)
/// leaves `pc` unknown; else the UV point is classified by
/// [`super::p2::point_classify`], `pc`'s gap becomes the distance from `point` to the
/// surface point, and a classification to the face, an edge or a vertex whose gap
/// exceeds the intersection tolerance of `src_zone_tol_3d` and the classified object's
/// zone tolerance is reset to unknown. `pc` starts unknown on this face with source
/// tolerance `src_zone_tol_3d`; its other members are kept unless written.
///
/// `b` is mutable because [`super::p2::point_classify`] is.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for a face without surface (upstream dereferences null);
/// those of the point solves and of [`super::p2::point_classify`]. `pc` keeps what was
/// written before the failure.
pub fn point_3d_classify(
    k: &Kernel,
    b: &mut Brep,
    f: FaceId,
    point: Point3,
    src_zone_tol_3d: f64,
    do_boundary_intersections: bool,
    pc: &mut PointClassification<EntityRef>,
    opts: Point3dClassifyOptions<'_>,
) -> Result<()> {
    let uv_domain = fp2::uv_domain(b, f)?;
    point_3d_classify_initialize(f, src_zone_tol_3d, pc);
    let uv = {
        let surface = fp2::surface(b, f)?
            .ok_or_else(|| SolidError::invalid_input(OP_POINT, "face without surface"))?;
        point_3d_classify_uv_with_context(
            k,
            surface,
            None::<NoObject>,
            &uv_domain,
            point,
            src_zone_tol_3d,
            opts.guess_uv,
            &mut OpaqueSurfaceCacheContext,
        )?
    };
    let Some(uv) = uv else {
        return Ok(());
    };
    point_3d_classify_finish(
        k,
        b,
        f,
        uv,
        point,
        src_zone_tol_3d,
        do_boundary_intersections,
        pc,
        opts,
    )
}

const POINT_CONTEXT: &str = "query::face::point_3d_classify_with_context";

/// Live caller topology, borrowed only for this face query. A surface
/// handle is distinct from its face owner; no B-rep cache is built in advance.
pub(crate) struct FacePointContext<'a> {
    brep: &'a mut Brep,
    keys: &'a dyn Fn(CacheObject) -> ObjectRef,
    owner_of: &'a dyn Fn(ObjectRef) -> Option<ObjectRef>,
}

/// Resolve native `GetFace` through the original B-rep and caller owner registry.
/// This read-only helper is shared by surface-cache contexts and their entry gate.
/// Surface-owner walks are iterative and bounded to 64 inspected links.
///
/// # Errors
/// Refuses unknown face/surface identities and cyclic or overlong owner chains.
pub(crate) fn surface_face_with_context(
    brep: &Brep,
    surface: &Surface,
    keys: &dyn Fn(CacheObject) -> ObjectRef,
    owner_of: &dyn Fn(ObjectRef) -> Option<ObjectRef>,
) -> Result<Option<FaceId>> {
    use crate::surface::surface::p1::is_surface_type;
    let mut owner = surface.owner();
    let mut inherited_face = surface.face()?;
    // Native GetFace walks only Surface owners and stops at Curve owners.
    // Keep the actual caller links authoritative; an inline owned base can
    // instead carry its enclosing surface's maintained owner_chain_face.
    for _ in 0..64 {
        let Some(current) = owner else {
            return Ok(None);
        };
        if current.type_code == fp2::TYPE {
            return brep
                .faces()
                .ids()
                .find(|&face| keys(CacheObject::Face(face)) == current)
                .map(Some)
                .ok_or_else(|| {
                    SolidError::invalid_input(
                        POINT_CONTEXT,
                        "face owner is outside the caller B-rep",
                    )
                });
        }
        if !is_surface_type(current.type_code) {
            return Ok(None);
        }
        if let Some(next) = owner_of(current) {
            owner = Some(next);
            inherited_face = None;
        } else if let Some(face) = brep.faces().ids().find(|&face| {
            keys(CacheObject::Surface(face)) == current
                && brep.face(face).is_ok_and(|record| record.surface.is_some())
        }) {
            let parent =
                fp2::surface(brep, face)?.ok_or(SolidError::null_pointer(POINT_CONTEXT))?;
            owner = parent.owner();
            inherited_face = parent.face()?;
        } else if let Some(face) = inherited_face {
            return brep
                .faces()
                .ids()
                .find(|&candidate| keys(CacheObject::Face(candidate)) == face)
                .map(Some)
                .ok_or_else(|| {
                    SolidError::invalid_input(
                        POINT_CONTEXT,
                        "face owner is outside the caller B-rep",
                    )
                });
        } else {
            return Err(SolidError::invalid_input(
                POINT_CONTEXT,
                "unregistered surface owner chain",
            ));
        }
    }
    Err(SolidError::invalid_input(
        POINT_CONTEXT,
        "cyclic surface owner chain",
    ))
}

impl<'a> FacePointContext<'a> {
    /// Borrow the original B-rep and its caller-owned identity/owner registry.
    /// Contextual consumers share this provider instead of opaque cache defaults.
    pub(crate) fn new(
        brep: &'a mut Brep,
        keys: &'a dyn Fn(CacheObject) -> ObjectRef,
        owner_of: &'a dyn Fn(ObjectRef) -> Option<ObjectRef>,
    ) -> Self {
        Self {
            brep,
            keys,
            owner_of,
        }
    }

    fn face(&self, object: CacheObject) -> Result<FaceId> {
        let CacheObject::Face(face) = object else {
            return Err(SolidError::invalid_input(
                POINT_CONTEXT,
                "expected face handle",
            ));
        };
        self.brep.face(face)?;
        Ok(face)
    }

    fn surface(&self, handle: Option<CacheObject>) -> Result<&Surface> {
        let Some(CacheObject::Surface(face)) = handle else {
            return Err(SolidError::invalid_input(
                POINT_CONTEXT,
                "expected surface handle",
            ));
        };
        fp2::surface(self.brep, face)?.ok_or(SolidError::null_pointer(POINT_CONTEXT))
    }

    fn resolved_surface_face(&self, surface: &Surface) -> Result<Option<CacheObject>> {
        surface_face_with_context(self.brep, surface, self.keys, self.owner_of)
            .map(|face| face.map(CacheObject::Face))
    }

    fn registered_object(&self, owner: ObjectRef) -> Option<CacheObject> {
        // Resolve complete identities, never derive an arena id from a key.
        let mut objects = self
            .brep
            .faces()
            .ids()
            .map(CacheObject::Face)
            .chain(self.brep.edges().ids().map(CacheObject::Edge))
            .chain(self.brep.vertices().ids().map(CacheObject::Vertex))
            .chain(std::iter::once(CacheObject::Brep))
            .chain(
                self.brep
                    .faces()
                    .ids()
                    .filter(|&face| {
                        self.brep
                            .face(face)
                            .is_ok_and(|record| record.surface.is_some())
                    })
                    .map(CacheObject::Surface),
            )
            .chain(
                self.brep
                    .edges()
                    .ids()
                    .filter(|&edge| {
                        self.brep
                            .edge(edge)
                            .is_ok_and(|record| record.curve.is_some())
                    })
                    .map(CacheObject::Curve),
            );
        objects.find(|&object| (self.keys)(object) == owner)
    }

    fn registered_zone_tolerance(&self, k: &Kernel, object: CacheObject) -> Result<f64> {
        if let CacheObject::Curve(edge) = object {
            let curve = self
                .brep
                .edge(edge)?
                .curve
                .as_ref()
                .ok_or(SolidError::null_pointer(POINT_CONTEXT))?;
            let mut unresolved = false;
            let owner = curve.edge(|intermediate| {
                let owner = (self.owner_of)(intermediate);
                unresolved |= owner.is_none();
                owner
            })?;
            if unresolved {
                return Err(SolidError::invalid_input(
                    POINT_CONTEXT,
                    "unregistered curve tolerance owner",
                ));
            }
            let edge = owner
                .map(|owner| {
                    self.registered_object(owner).ok_or_else(|| {
                        SolidError::invalid_input(
                            POINT_CONTEXT,
                            "curve owner leaves the caller B-rep",
                        )
                    })
                })
                .transpose()?;
            return self.object_zone_tol_3d(k, edge);
        }
        if let CacheObject::Surface(face) = object {
            let surface = self.surface(Some(CacheObject::Surface(face)))?;
            return self.object_zone_tol_3d(k, self.resolved_surface_face(surface)?);
        }
        self.object_zone_tol_3d(k, Some(object))
    }

    fn owner_zone_tolerance(&self, k: &Kernel, owner: Option<ObjectRef>) -> Result<f64> {
        use crate::curve::curve::p2::is_curve_class;
        use crate::surface::surface::p1::is_surface_type;
        let Some(mut current) = owner else {
            return self.object_zone_tol_3d(k, None);
        };
        if let Some(object) = self.registered_object(current) {
            return self.registered_zone_tolerance(k, object);
        }
        // An intermediate geometry object need not be a top-level arena slot.
        // Its caller-owned links still have to resolve exactly; an unknown owner
        // is not an ownerless object and must never become a system tolerance.
        let surface_chain = is_surface_type(current.type_code);
        let curve_chain = is_curve_class(current.type_code)
            || current.type_code == crate::topology::edgeuse::p1::TYPE;
        if !surface_chain && !curve_chain {
            return Err(SolidError::invalid_input(
                POINT_CONTEXT,
                "unregistered tolerance owner",
            ));
        }
        for _ in 0..64 {
            let next = (self.owner_of)(current).ok_or_else(|| {
                SolidError::invalid_input(
                    POINT_CONTEXT,
                    "unregistered intermediate tolerance owner",
                )
            })?;
            if (surface_chain && next.type_code == fp2::TYPE)
                || (curve_chain && next.type_code == edge::TYPE)
            {
                let object = self.registered_object(next).ok_or_else(|| {
                    SolidError::invalid_input(POINT_CONTEXT, "owner chain leaves the caller B-rep")
                })?;
                return self.registered_zone_tolerance(k, object);
            }
            let continues = if surface_chain {
                is_surface_type(next.type_code)
            } else {
                is_curve_class(next.type_code)
                    || next.type_code == crate::topology::edgeuse::p1::TYPE
            };
            if !continues {
                // GetFace stops at curve owners; GetEdge stops at non-curve and
                // non-edgeuse owners. Validate that terminal identity first.
                self.registered_object(next).ok_or_else(|| {
                    SolidError::invalid_input(
                        POINT_CONTEXT,
                        "unregistered terminal tolerance owner",
                    )
                })?;
                return self.object_zone_tol_3d(k, None);
            }
            // A represented geometry has its own exact GetFace/GetEdge state.
            if let Some(object) = self.registered_object(next) {
                return self.registered_zone_tolerance(k, object);
            }
            current = next;
        }
        Err(SolidError::invalid_input(
            POINT_CONTEXT,
            "cyclic tolerance owner chain",
        ))
    }

    fn entity(&self, object: CacheObject) -> Result<EntityRef> {
        match object {
            CacheObject::Face(face) => {
                self.brep.face(face)?;
                Ok(EntityRef::Face(face))
            }
            CacheObject::Edge(edge) => {
                self.brep.edge(edge)?;
                Ok(EntityRef::Edge(edge))
            }
            CacheObject::Vertex(vertex) => {
                self.brep.vertex(vertex)?;
                Ok(EntityRef::Vertex(vertex))
            }
            _ => Err(SolidError::invalid_input(
                POINT_CONTEXT,
                "expected classified topology",
            )),
        }
    }

    fn classification_object(&self, entity: EntityRef) -> Result<CacheObject> {
        let object = match entity {
            EntityRef::Face(face) => CacheObject::Face(face),
            EntityRef::Edge(edge) => CacheObject::Edge(edge),
            EntityRef::Vertex(vertex) => CacheObject::Vertex(vertex),
            _ => {
                return Err(SolidError::invalid_input(
                    POINT_CONTEXT,
                    "unexpected classification entity",
                ));
            }
        };
        self.entity(object)?;
        Ok(object)
    }
}

impl SurfaceCacheContext<CacheObject> for FacePointContext<'_> {
    fn surface_face(
        &self,
        surface: &Surface,
        handle: Option<CacheObject>,
    ) -> Result<Option<CacheObject>> {
        self.surface(handle)?;
        self.resolved_surface_face(surface)
    }

    fn surface_has_direct_face_owner(
        &self,
        surface: &Surface,
        handle: Option<CacheObject>,
    ) -> Result<bool> {
        self.surface(handle)?;
        Ok(surface
            .owner()
            .is_some_and(|owner| owner.type_code == crate::topology::face::p2::TYPE))
    }

    fn surface_zone_tol_3d(
        &self,
        k: &Kernel,
        surface: &Surface,
        handle: Option<CacheObject>,
    ) -> Result<f64> {
        self.surface(handle)?;
        // The evaluated surface may be an owned base, not the top-level handle.
        let face = self.resolved_surface_face(surface)?;
        self.object_zone_tol_3d(k, face)
    }

    fn surface_owner_zone_tol_3d(
        &self,
        k: &Kernel,
        surface: &Surface,
        handle: Option<CacheObject>,
    ) -> Result<f64> {
        self.surface(handle)?;
        self.owner_zone_tolerance(k, surface.owner())
    }

    fn object_zone_tol_3d(&self, k: &Kernel, object: Option<CacheObject>) -> Result<f64> {
        object_zone_tol_3d(k, self.brep, object, self.keys, self.owner_of)
    }

    fn face_brep_editing_enabled(&self, face: CacheObject) -> Result<bool> {
        self.face(face)?;
        Ok(self.brep.editing_enabled())
    }

    fn set_face_brep_editing_enabled(&mut self, face: CacheObject, editing: bool) -> Result<()> {
        self.face(face)?;
        self.brep.set_editing_enabled(editing);
        Ok(())
    }

    fn face_create_uv_trim_curves(&mut self, k: &Kernel, face: CacheObject) -> Result<()> {
        let face = self.face(face)?;
        crate::topology::face::p3::create_uv_trim_curves(k, self.brep, face, true, None).map(|_| ())
    }

    fn face_point_classify(
        &mut self,
        k: &Kernel,
        face: CacheObject,
        uv: Point2,
        tolerance: f64,
        do_boundaries: bool,
        use_uv: bool,
        classification: &mut PointClassification<CacheObject>,
    ) -> Result<()> {
        let face = self.face(face)?;
        let mut local = classification.try_map_objects(|object| self.entity(object))?;
        let result = crate::query::face::p2::point_classify_with_context(
            k,
            self.brep,
            face,
            uv,
            &mut local,
            crate::query::face::p2::PointClassifyOptions {
                src_zone_tol_3d: tolerance,
                do_boundary_intersections: do_boundaries,
                use_uv_classification: use_uv,
                ..crate::query::face::p2::PointClassifyOptions::default()
            },
            self.keys,
            self.owner_of,
        );
        // The native status does not discard writes to any classification field.
        *classification = local.try_map_objects(|entity| self.classification_object(entity))?;
        result
    }
}

/// Classify on the caller's live mutable B-rep with its actual identity registry.
/// Cache acquisition and lazy UV writes operate on `b`; only the bounded numerical
/// surface view is copied for evaluation. The same `pc` retains all partial writes.
/// `keys` and `owner_of` must resolve the actual owners of this B-rep, including
/// intermediate owner chains; an unresolved owner is not treated as absent.
///
/// # Errors
/// Face/geometry/owner lookup, cache construction, point solve and classification
/// errors propagate without discarding earlier B-rep or classification changes.
pub fn point_3d_classify_with_context(
    k: &Kernel,
    b: &mut Brep,
    face: FaceId,
    point: Point3,
    source_tolerance: f64,
    do_boundaries: bool,
    pc: &mut PointClassification<EntityRef>,
    opts: Point3dClassifyOptions<'_>,
    keys: &dyn Fn(CacheObject) -> ObjectRef,
    owner_of: &dyn Fn(ObjectRef) -> Option<ObjectRef>,
) -> Result<()> {
    point_3d_classify_with_cache_context(
        k,
        b,
        face,
        point,
        source_tolerance,
        do_boundaries,
        pc,
        opts,
        keys,
        owner_of,
        None,
    )
}

/// Classify on the original B-rep with its actual retained trimmed-surface cache.
/// The cache is borrowed through the final UV classification; it is never copied
/// or substituted for the surface-drop cache. All partial writes survive errors.
///
/// # Errors
/// Those of [`point_3d_classify_with_context`] and the retained cache's real queries.
pub fn point_3d_classify_with_cache_context<'cache>(
    k: &Kernel,
    b: &mut Brep,
    face: FaceId,
    point: Point3,
    source_tolerance: f64,
    do_boundaries: bool,
    pc: &mut PointClassification<EntityRef>,
    opts: Point3dClassifyOptions<'cache>,
    keys: &dyn Fn(CacheObject) -> ObjectRef,
    owner_of: &dyn Fn(ObjectRef) -> Option<ObjectRef>,
    trim_cache: Option<&'cache mut dyn super::p2::FacePointCache>,
) -> Result<()> {
    let domain = fp2::uv_domain(b, face)?;
    point_3d_classify_initialize(face, source_tolerance, pc);
    let uv = with_face_surface_evaluation(b, face, |surface, brep| {
        let mut context = FacePointContext::new(brep, keys, owner_of);
        point_3d_classify_uv_with_context(
            k,
            surface,
            Some(CacheObject::Surface(face)),
            &domain,
            point,
            source_tolerance,
            opts.guess_uv,
            &mut context,
        )
    })?;
    let Some(uv) = uv else {
        return Ok(());
    };
    point_3d_classify_finish_impl(
        k,
        b,
        face,
        uv,
        point,
        source_tolerance,
        do_boundaries,
        pc,
        opts,
        Some(PointFinishContext {
            keys,
            owner_of,
            trim_cache,
        }),
    )
}

/// Initialise the fields written by `SmFace::Point3DClassify`, retaining all
/// other fields of the caller's classification. Call this after reading the
/// face's UV domain and before resolving its surface, as the native method does.
pub fn point_3d_classify_initialize(
    face: FaceId,
    source_zone_tolerance: f64,
    classification: &mut PointClassification<EntityRef>,
) {
    // SM_PC_UNKNOWN is also used for points outside the face's trim boundary.
    classification.set_class_object(PointClassType::Unknown, None);
    classification.set_face(Some(EntityRef::Face(face)));
    classification.set_src_zone_tol_3d(source_zone_tolerance);
}

/// Runs a face query with an exact numerical surface view and the original mutable
/// B-rep. The view retains its surface's owner identity and NURBS data; resolve
/// topology, owner state and trim writes through `brep` while using the view for
/// numerical evaluation. The view lives only for the callback, so the returned
/// value cannot borrow it. Changes made to `brep` remain visible even on error.
///
/// For `Point3DClassify`, read the UV domain and call
/// [`point_3d_classify_initialize`] before entering this bridge, preserving the
/// native order of classification writes and surface lookup.
///
/// # Errors
///
/// Propagates face lookup errors, [`SolidError::InvalidInput`] for a face without
/// a surface or geometry exceeding the evaluation-copy depth bound, and the
/// callback's error without replacing it.
pub fn with_face_surface_evaluation<T>(
    brep: &mut Brep,
    face: FaceId,
    query: impl FnOnce(&Surface, &mut Brep) -> Result<T>,
) -> Result<T> {
    let surface = fp2::surface(brep, face)?
        .ok_or_else(|| SolidError::invalid_input(OP_POINT, "face without surface"))?
        .clone_for_evaluation()?;
    query(&surface, brep)
}

/// The surface-drop phase of `SmFace::Point3DClassify`, with the live owner and
/// trim-creation context supplied separately from the numerical surface view.
/// Cache acquisition precedes even a hinted local drop. The global drop reuses
/// that checkout and disables its point filter; explicit UV classification is
/// performed by [`point_3d_classify_finish`].
///
/// Returns `None` when there is no drop or it falls outside the face's UV domain.
/// The caller must first use [`point_3d_classify_initialize`] so these successful
/// early exits retain the native partial output.
///
/// # Errors
///
/// Propagates point-solve errors and unresolved owner/context errors. Ordinary
/// native cache acquisition failures permit a null cache, as upstream does.
pub fn point_3d_classify_uv_with_context<O: SolutionObject>(
    k: &Kernel,
    surface: &Surface,
    surface_handle: Option<O>,
    uv_domain: &Extent2d,
    point: Point3,
    src_zone_tol_3d: f64,
    guess_uv: Option<Point2>,
    context: &mut dyn SurfaceCacheContext<O>,
) -> Result<Option<Point2>> {
    let mut cache = match get_surface_cache_with_context(k, surface, surface_handle, context) {
        Ok(cache) => Some(cache),
        Err(SolidError::Upstream { .. } | SolidError::Nlib { .. }) => None,
        Err(error) => return Err(error),
    };

    // When given a hint, use it in LocalSolve before trying a global solve.
    let mut found = None;
    if let Some(guess) = guess_uv {
        let result = surface.local_point_solve_with_tolerance_context::<NoObject>(
            k,
            uv_domain,
            SolverOperation::Intersect,
            point,
            guess,
            &mut |surface| context.surface_owner_zone_tol_3d(k, surface, surface_handle),
            &mut |surface| {
                context
                    .surface_zone_tol_3d(k, surface, surface_handle)
                    .map(Some)
            },
        )?;
        if result.found {
            found = Some(solution_uv(&result.solution.start.parameters)?);
        }
    }

    let uv = if let Some(uv) = found {
        uv
    } else {
        let solutions = surface.global_point_intersect_with_cache_context(
            k,
            surface_handle,
            uv_domain,
            point,
            src_zone_tol_3d,
            SolutionRequested::Single,
            FacePointTest::Disabled,
            &mut cache,
            context,
        )?;
        let Some(first) = solutions.as_slice().first() else {
            return Ok(None);
        };
        solution_uv(&first.start.parameters)?
    };

    // low work - UVPoint not contained by Face's domain
    let size = uv_domain.size().length();
    if !uv_domain.contains_point(uv, size / 1000.0) {
        return Ok(None);
    }
    Ok(Some(uv))
}

/// Complete `SmFace::Point3DClassify` after a successful in-domain surface drop.
/// This is the shared explicit UV classification, 3D-gap update and final object
/// tolerance check; `b` is the original B-rep, including trims created during
/// cache acquisition. The same classification record is used throughout.
///
/// # Errors
///
/// Propagates UV classification, geometry and tolerance errors, retaining every
/// field written to `pc` before a failure.
pub fn point_3d_classify_finish(
    k: &Kernel,
    b: &mut Brep,
    f: FaceId,
    uv: Point2,
    point: Point3,
    src_zone_tol_3d: f64,
    do_boundary_intersections: bool,
    pc: &mut PointClassification<EntityRef>,
    opts: Point3dClassifyOptions<'_>,
) -> Result<()> {
    point_3d_classify_finish_impl(
        k,
        b,
        f,
        uv,
        point,
        src_zone_tol_3d,
        do_boundary_intersections,
        pc,
        opts,
        None,
    )
}

struct PointFinishContext<'identity, 'cache> {
    keys: &'identity dyn Fn(CacheObject) -> ObjectRef,
    owner_of: &'identity dyn Fn(ObjectRef) -> Option<ObjectRef>,
    trim_cache: Option<&'cache mut dyn super::p2::FacePointCache>,
}

fn point_3d_classify_finish_impl<'cache>(
    k: &Kernel,
    b: &mut Brep,
    f: FaceId,
    uv: Point2,
    point: Point3,
    src_zone_tol_3d: f64,
    do_boundary_intersections: bool,
    pc: &mut PointClassification<EntityRef>,
    opts: Point3dClassifyOptions<'cache>,
    context: Option<PointFinishContext<'_, 'cache>>,
) -> Result<()> {
    // Preserve the actual retained cache through the native final PointClassify.
    let (identity, trim_cache) = match context {
        Some(context) => (Some((context.keys, context.owner_of)), context.trim_cache),
        None => (None, None),
    };
    let options = super::p2::PointClassifyOptions {
        src_zone_tol_3d,
        do_boundary_intersections,
        use_uv_classification: false,
        tgt_loops_only: opts.tgt_loops_only,
        tgt_loops: opts.tgt_loops,
        uv_space_okay: opts.uv_space_okay,
        trim_cache,
    };
    if let Some((keys, owner_of)) = identity {
        super::p2::point_classify_with_context(k, b, f, uv, pc, options, keys, owner_of)?;
    } else {
        super::p2::point_classify(k, b, f, uv, pc, options)?;
    }

    // PointClassify set the gap to that of the surface point; use the distance between
    // the 3D point and its surface point instead.
    let surface = fp2::surface(b, f)?
        .ok_or_else(|| SolidError::invalid_input(OP_POINT, "face without surface"))?;
    if let Some(surface_point) = ignored(surface.evaluate_point(uv))? {
        pc.set_gap_3d(point.distance(surface_point));
    }

    // If we classified to something but not within the specified tolerance, return no
    // classification.
    if pc.point_class() != PointClassType::Unknown && pc.point_class() != PointClassType::Region {
        let geometry = BrepGeometry { b };
        let xsect_tol_3d =
            xsect_tol_3d_from_zones(k, src_zone_tol_3d, pc.obj_zone_tol_3d(k, &geometry)?);
        if pc.gap_3d(k, &geometry, None)? > xsect_tol_3d {
            pc.set_class_object(PointClassType::Unknown, None);
        }
    }

    Ok(())
}

// ---- SmFace::EdgeCurveOnClassify -----------------------------------------------------------

/// Why a sector classification may be wrong (`SmClassifyProblem`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ClassifyProblem {
    /// `SM_CP_NONE`: the classification can be trusted.
    #[default]
    None = 0,
    /// `SM_CP_VERTEX_ON_SINGULARITY`.
    VertexOnSingularity = 1,
    /// `SM_CP_TOGGLED_CLASSIFICATION`: growing step-offs changed it once.
    ToggledClassification = 2,
    /// `SM_CP_SHORT_SECTOR_EDGEUSE1`: the first edge use is shorter than the step-off.
    ShortSectorEdgeuse1 = 3,
    /// `SM_CP_SHORT_SECTOR_EDGEUSE2`: the second edge use is shorter than the step-off.
    ShortSectorEdgeuse2 = 4,
    /// `SM_CP_SHORT_INTERVAL`: the curve interval is shorter than the step-off.
    ShortInterval = 5,
    /// `SM_CP_SEESAW_CLASSIFICATION`: growing step-offs changed it more than once.
    SeesawClassification = 6,
    /// `SM_CP_UNRESOLVED_CLASSIFICATION`: the iteration limit was reached.
    UnresolvedClassification = 7,
}

/// The answers of [`edge_curve_on_classify`] (`rbIsInsideFace`, `rbPossibleProblem`,
/// `reClassifyProblem`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CurveOnClassification {
    /// Whether the classified curve segment is inside the face.
    pub inside: bool,
    /// Whether a classification reported a possible problem.
    pub possible_problem: bool,
    /// The problem of the last classification that reported one.
    pub problem: ClassifyProblem,
}

/// The arguments of [`edge_curve_on_classify`] after the edge and curve.
#[derive(Debug, Clone, Copy)]
pub struct EdgeCurveOnOptions {
    /// `dEdgeParam`: the edge parameter of the edge/curve intersection.
    pub edge_param: f64,
    /// `crCurveInterval`: the range of interest of the curve.
    pub curve_interval: Extent1d,
    /// `eCurveOrientation`: [`OrientType::Same`] for the curve segment after the
    /// intersection, [`OrientType::Opposite`] for the one before.
    pub curve_orientation: OrientType,
    /// `d3dCurveTolerance`: the zone tolerance of the curve's object.
    pub curve_tolerance: f64,
    /// `bUVSpaceOkay` (upstream default `TRUE`).
    pub uv_space_okay: bool,
    /// `lFaceIndx`, for debugging a target face (upstream default `SM_UNDEF_ULONG`).
    pub face_index: Option<u64>,
}

impl EdgeCurveOnOptions {
    /// The required arguments with upstream's defaults for the others.
    #[must_use]
    pub const fn new(
        edge_param: f64,
        curve_interval: Extent1d,
        curve_orientation: OrientType,
        curve_tolerance: f64,
    ) -> Self {
        Self {
            edge_param,
            curve_interval,
            curve_orientation,
            curve_tolerance,
            uv_space_okay: true,
            face_index: None,
        }
    }
}

/// `SmFace::SectorClassifyCurveOn` with the same edge use as both sector bounds (in
/// opposite directions), the edge's tolerance for both and no singularity.
fn sector_classify_curve_on(
    k: &Kernel,
    b: &mut Brep,
    f: FaceId,
    edgeuse: EdgeuseId,
    edge_tolerance: f64,
    curve: &Curve,
    opts: &EdgeCurveOnOptions,
) -> Result<CurveOnClassification> {
    let sector_edge = |orientation| p3::SectorEdgeuse {
        edgeuse,
        param: opts.edge_param,
        zone_tol_3d: edge_tolerance,
        orientation,
    };
    let result = p3::sector_classify_curve_on(
        k,
        b,
        f,
        &p3::SectorClassifyInput {
            eu1: sector_edge(OrientType::Same),
            eu2: sector_edge(OrientType::Opposite),
            curve,
            curve_interval: opts.curve_interval,
            curve_orientation: opts.curve_orientation,
            curve_zone_tol_3d: opts.curve_tolerance,
            on_singularity: false,
            singularity_uv: Point2::new(0.0, 0.0),
            uv_space_okay: opts.uv_space_okay,
        },
    )?;
    Ok(CurveOnClassification {
        inside: result.inside,
        possible_problem: result.possible_problem,
        problem: result.problem,
    })
}

/// The 3D curve of a 2D `curve` on face `f`'s surface (`SmCrvOnSurf::SetAll`).
fn lift_to_face(k: &Kernel, b: &Brep, f: FaceId, curve: &Curve) -> Result<Curve> {
    let surface = fp2::surface(b, f)?
        .ok_or_else(|| SolidError::invalid_input(OP_EDGE_CURVE, "face without surface"))?;
    let SurfaceKind::BSpline(base) = surface.kind() else {
        // DISPATCH GAP: Rust's SmCrvOnSurf holds a B-spline surface only (curve.crv_on_surf,
        // partial; no UNPORTED marker)
        return Err(SolidError::unported(
            "SmSurface (non-B-spline) in SmCrvOnSurf",
            "curve.crv_on_surf",
        ));
    };
    Ok(Curve::new(CurveKind::CrvOnSurf(Box::new(CrvOnSurf::new(
        k,
        curve.clone(),
        base.clone(),
        None,
        0,
    )?))))
}

/// `SmFace::EdgeCurveOnClassify`: classifies the segment of `curve` that crosses the
/// face's boundary at `edge` (at [`EdgeCurveOnOptions::edge_param`]) as inside or outside
/// the face, against every upward edge use of the edge on this face (a seam edge has
/// two; an outside answer does not stop the search). A possible problem is passed on;
/// the first inside answer ends the search. A 2D `curve` is taken as a curve on the
/// face's surface.
///
/// `b` is mutable because the sector classification is.
///
/// # Errors
///
/// `SM_ERR_NULL_POINTER` for a null `edge`; those of the edge use list, the curve lift
/// and the sector classification ([`p3::sector_classify_curve_on`]).
pub fn edge_curve_on_classify(
    k: &Kernel,
    b: &mut Brep,
    f: FaceId,
    edge: Option<EdgeId>,
    curve: &Curve,
    opts: &EdgeCurveOnOptions,
) -> Result<CurveOnClassification> {
    let edge = edge.ok_or_else(|| SolidError::null_pointer(OP_EDGE_CURVE))?;

    // init outputs
    let mut out = CurveOnClassification::default();

    let edgeuses = edge::edgeuses(b, edge, None, None)?;
    let upward = fp2::upward_faceuse(b, f)?;

    // work in 3d: a 2D curve is taken as an SmCrvOnSurf on this face's surface
    let mut lifted = None;

    // for every edgeuse
    for eu in edgeuses {
        // skip all but the upward edgeuses (closed curves get 2) connected to this face
        if eu1::face(b, eu)? != Some(f) || eu1::faceuse(b, eu)? != upward {
            continue;
        }

        // Have found upward edgeuse connected to this face - classify curve relative to it
        let eu_tol = edge::tolerance(b, edge)?;
        let this_curve = if curve.dim() == 2 {
            if lifted.is_none() {
                lifted = Some(lift_to_face(k, b, f, curve)?);
            }
            lifted
                .as_ref()
                .ok_or_else(|| SolidError::invalid_input(OP_EDGE_CURVE, "missing lifted curve"))?
        } else {
            curve
        };
        let mine = sector_classify_curve_on(k, b, f, eu, eu_tol, this_curve, opts)?;

        // pass on any detected problems
        if mine.possible_problem {
            out.possible_problem = true;
            out.problem = mine.problem;
        }

        // pass on the answer
        if mine.inside {
            out.inside = true;
            return Ok(out);
        }

        // not all done: a seam edge has two upward edge uses connected to this face
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::line::Line;

    #[test]
    fn evaluate_curve_rejects_null_nan_and_outside() {
        assert!(evaluate_curve(None, 0.0).is_err());
        let c: Curve = Line::from_points(Vec3::new(0., 0., 0.), Vec3::new(1., 2., 3.), 3).into();
        assert!(evaluate_curve(Some(&c), f64::NAN).is_err());
        assert!(evaluate_curve(Some(&c), 2.0).is_err());
        let [p, d1, _] = evaluate_curve(Some(&c), 0.5).unwrap();
        assert_eq!(p, Vec3::new(0.5, 1.0, 1.5));
        assert_eq!(d1, Vec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn region_without_shells_is_sm_err() {
        let mut k = Kernel::default();
        let mut b = Brep::default();
        let r = b
            .add_region_record(crate::topology::Region::default())
            .unwrap();
        let mut out = Extent3d::from_point(Vec3::new(1., 2., 3.));
        let e = calculate_region_bounding_box_into(&mut k, &mut b, r, false, &mut out).unwrap_err();
        assert!(matches!(e, SolidError::Upstream { .. }));
        assert_eq!(out, Extent3d::new());
    }

    #[test]
    fn null_edge_is_null_pointer() {
        let k = Kernel::default();
        let mut b = Brep::default();
        let c: Curve = Line::from_points(Vec3::new(0., 0., 0.), Vec3::new(1., 0., 0.), 3).into();
        let opts = EdgeCurveOnOptions::new(
            0.0,
            Extent1d::from_min_max(0.0, 1.0),
            OrientType::Same,
            1e-6,
        );
        let e = edge_curve_on_classify(&k, &mut b, FaceId::from_index(0), None, &c, &opts);
        assert!(e.is_err());
    }
}

#[cfg(test)]
mod live_context_tests {
    use super::*;
    use crate::geom::attribute::ObjectKey;
    use crate::surface::plane::Plane;

    #[test]
    fn context_resolves_actual_owner_keys_and_mutates_the_original_editing_flag() {
        let mut brep = Brep::new();
        let shell = brep.add_shell().unwrap();
        let face = brep.add_face(shell).unwrap();
        let key = ObjectRef {
            type_code: fp2::TYPE,
            key: ObjectKey(0x1234_5678),
        };
        let mut surface = Surface::from(Plane::from_normal(Point3::ZERO, Vec3::new(0.0, 0.0, 1.0)));
        surface.set_owner(Some(key));
        brep.face_mut(face).unwrap().surface = Some(surface);
        let keys = |object| match object {
            CacheObject::Face(id) if id == face => key,
            _ => panic!("only the real face owner is queried"),
        };
        {
            let mut context = FacePointContext {
                brep: &mut brep,
                keys: &keys,
                owner_of: &|_| None,
            };
            let handle = Some(CacheObject::Surface(face));
            assert_eq!(
                context
                    .surface_face(context.surface(handle).unwrap(), handle)
                    .unwrap(),
                Some(CacheObject::Face(face))
            );
            let unknown = ObjectRef {
                key: ObjectKey(key.key.0 + 1),
                ..key
            };
            let mut evaluated =
                Surface::from(Plane::from_normal(Point3::ZERO, Vec3::new(0.0, 0.0, 1.0)));
            evaluated.set_owner(Some(unknown));
            assert!(context.surface_face(&evaluated, handle).is_err());
            context
                .set_face_brep_editing_enabled(CacheObject::Face(face), true)
                .unwrap();
        }
        assert!(brep.editing_enabled());
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn context_uses_the_evaluated_surface_and_resolves_non_face_owner_tolerances() {
        use crate::curve::curve::Curve;
        use crate::curve::line::Line;
        let kernel = Kernel::default();
        let mut brep = Brep::new();
        let shell = brep.add_shell().unwrap();
        let face = brep.add_face(shell).unwrap();
        let start = brep.add_vertex(Point3::ZERO).unwrap();
        let end = brep.add_vertex(Point3::new(1.0, 0.0, 0.0)).unwrap();
        let edge = brep.add_edge(start, end).unwrap();
        let face_ref = ObjectRef {
            type_code: fp2::TYPE,
            key: ObjectKey(1000),
        };
        let edge_ref = ObjectRef {
            type_code: edge::TYPE,
            key: ObjectKey(2000),
        };
        let vertex_ref = ObjectRef {
            type_code: vertex::TYPE,
            key: ObjectKey(3000),
        };
        let brep_ref = ObjectRef {
            type_code: crate::geom::a_object::BREP_TYPE,
            key: ObjectKey(4000),
        };
        let mut curve = Curve::from(Line::from_points(
            Point3::ZERO,
            Point3::new(1.0, 0.0, 0.0),
            3,
        ));
        let curve_ref = ObjectRef {
            type_code: curve.type_code(),
            key: ObjectKey(5000),
        };
        let _ = curve.set_owner(Some(edge_ref), Some(brep_ref));
        brep.edge_mut(edge).unwrap().curve = Some(curve);
        brep.edge_mut(edge).unwrap().tolerance = 0.125;
        brep.vertex_mut(start).unwrap().tolerance = 0.25;
        let mut surface = Surface::from(Plane::from_normal(Point3::ZERO, Vec3::new(0.0, 0.0, 1.0)));
        let surface_ref = ObjectRef {
            type_code: surface.type_code(),
            key: ObjectKey(6000),
        };
        surface.set_owner(Some(face_ref));
        brep.face_mut(face).unwrap().surface = Some(surface);
        brep.face_mut(face).unwrap().tolerance = 0.5;
        let keys = |object| match object {
            CacheObject::Face(id) if id == face => face_ref,
            CacheObject::Edge(id) if id == edge => edge_ref,
            CacheObject::Vertex(id) if id == start => vertex_ref,
            CacheObject::Brep => brep_ref,
            CacheObject::Curve(id) if id == edge => curve_ref,
            CacheObject::Surface(id) if id == face => surface_ref,
            CacheObject::Other(owner) => owner,
            _ => panic!("unexpected registration lookup"),
        };
        // Register all vertices because resolution scans the complete identity set.
        let keys = |object| match object {
            CacheObject::Vertex(id) if id == end => ObjectRef {
                type_code: vertex::TYPE,
                key: ObjectKey(3001),
            },
            object => keys(object),
        };
        let intermediate = ObjectRef {
            type_code: surface_ref.type_code,
            key: ObjectKey(7000),
        };
        let intermediate_curve = ObjectRef {
            type_code: curve_ref.type_code,
            key: ObjectKey(7001),
        };
        let cycle = ObjectRef {
            type_code: surface_ref.type_code,
            key: ObjectKey(7002),
        };
        let owner_of = |owner| {
            if owner == intermediate {
                Some(face_ref)
            } else if owner == intermediate_curve {
                Some(edge_ref)
            } else if owner == cycle {
                Some(cycle)
            } else {
                None
            }
        };
        let context = FacePointContext {
            brep: &mut brep,
            keys: &keys,
            owner_of: &owner_of,
        };
        let mut evaluated =
            Surface::from(Plane::from_normal(Point3::ZERO, Vec3::new(0.0, 0.0, 1.0)));
        let handle = Some(CacheObject::Surface(face));
        for (owner, expected) in [
            (edge_ref, 0.125),
            (vertex_ref, 0.25),
            (curve_ref, 0.125),
            (surface_ref, 0.5),
            (brep_ref, context.brep.tolerance()),
            (intermediate, 0.5),
            (intermediate_curve, 0.125),
        ] {
            evaluated.set_owner(Some(owner));
            assert_eq!(
                context
                    .surface_owner_zone_tol_3d(&kernel, &evaluated, handle)
                    .unwrap(),
                expected
            );
        }
        evaluated.set_surface_owner(intermediate, Some(face_ref));
        assert_eq!(
            context.surface_face(&evaluated, handle).unwrap(),
            Some(CacheObject::Face(face))
        );
        assert!(!context
            .surface_has_direct_face_owner(&evaluated, handle)
            .unwrap());
        evaluated.set_owner(Some(edge_ref));
        assert!(!context
            .surface_has_direct_face_owner(&evaluated, handle)
            .unwrap());
        assert_eq!(context.surface_face(&evaluated, handle).unwrap(), None);
        assert_eq!(
            context
                .surface_zone_tol_3d(&kernel, &evaluated, handle)
                .unwrap(),
            object_zone_tol_3d(&kernel, context.brep, None, &keys, &owner_of).unwrap()
        );
        evaluated.set_owner(Some(cycle));
        assert!(matches!(
            context.surface_owner_zone_tol_3d(&kernel, &evaluated, handle),
            Err(SolidError::InvalidInput { .. })
        ));
        for type_code in [
            edge_ref.type_code,
            curve_ref.type_code,
            surface_ref.type_code,
        ] {
            evaluated.set_owner(Some(ObjectRef {
                type_code,
                key: ObjectKey(9999),
            }));
            assert!(matches!(
                context.surface_owner_zone_tol_3d(&kernel, &evaluated, handle),
                Err(SolidError::InvalidInput { .. })
            ));
        }
    }

    #[test]
    fn evaluated_surface_resolves_live_surface_chains_without_inherited_face_metadata() {
        let kernel = Kernel::default();
        let mut brep = Brep::new();
        let shell = brep.add_shell().unwrap();
        let face = brep.add_face(shell).unwrap();
        let face_ref = ObjectRef {
            type_code: fp2::TYPE,
            key: ObjectKey(42_001),
        };
        let mut evaluated =
            Surface::from(Plane::from_normal(Point3::ZERO, Vec3::new(0.0, 0.0, 1.0)));
        let first = ObjectRef {
            type_code: evaluated.type_code(),
            key: ObjectKey(42_002),
        };
        let second = ObjectRef {
            type_code: evaluated.type_code(),
            key: ObjectKey(42_003),
        };
        let cycle = ObjectRef {
            type_code: evaluated.type_code(),
            key: ObjectKey(42_004),
        };
        let stop = ObjectRef {
            type_code: evaluated.type_code(),
            key: ObjectKey(42_005),
        };
        let curve = ObjectRef {
            type_code: crate::curve::curve::p2::TYPE,
            key: ObjectKey(42_006),
        };
        let unknown = ObjectRef {
            type_code: evaluated.type_code(),
            key: ObjectKey(42_007),
        };
        let mut stored = Surface::from(Plane::from_normal(Point3::ZERO, Vec3::new(0.0, 0.0, 1.0)));
        stored.set_owner(Some(face_ref));
        brep.face_mut(face).unwrap().surface = Some(stored);
        brep.face_mut(face).unwrap().tolerance = 0.125;
        let keys = |object| match object {
            CacheObject::Face(id) if id == face => face_ref,
            _ => crate::topology::brep_cache::arena_identity(object),
        };
        let owner_of = |owner| {
            if owner == first {
                Some(second)
            } else if owner == second {
                Some(face_ref)
            } else if owner == cycle {
                Some(cycle)
            } else if owner == stop {
                Some(curve)
            } else {
                None
            }
        };
        let context = FacePointContext::new(&mut brep, &keys, &owner_of);
        let handle = Some(CacheObject::Surface(face));
        let system = context.object_zone_tol_3d(&kernel, None).unwrap();
        assert_ne!(system, 0.125);
        evaluated.set_owner(Some(first));
        assert_eq!(evaluated.face().unwrap(), None);
        assert_eq!(
            context.surface_face(&evaluated, handle).unwrap(),
            Some(CacheObject::Face(face))
        );
        assert_eq!(
            context
                .surface_zone_tol_3d(&kernel, &evaluated, handle)
                .unwrap(),
            0.125
        );
        assert!(!context
            .surface_has_direct_face_owner(&evaluated, handle)
            .unwrap());
        for owner in [unknown, cycle] {
            evaluated.set_owner(Some(owner));
            assert!(matches!(
                context.surface_face(&evaluated, handle),
                Err(SolidError::InvalidInput { .. })
            ));
            assert!(matches!(
                context.surface_zone_tol_3d(&kernel, &evaluated, handle),
                Err(SolidError::InvalidInput { .. })
            ));
        }
        // A Curve owner stops GetFace, even when that curve could have an Edge.
        evaluated.set_owner(Some(stop));
        assert_eq!(context.surface_face(&evaluated, handle).unwrap(), None);
        assert_eq!(
            context
                .surface_zone_tol_3d(&kernel, &evaluated, handle)
                .unwrap(),
            system
        );
    }

    #[test]
    fn missing_surface_keeps_native_partial_classification_on_error() {
        let kernel = Kernel::default();
        let mut brep = Brep::new();
        let shell = brep.add_shell().unwrap();
        let face = brep.add_face(shell).unwrap();
        let mut pc = PointClassification::new(0.25);
        pc.set_gap_3d(17.0);
        let result = point_3d_classify_with_context(
            &kernel,
            &mut brep,
            face,
            Point3::ZERO,
            0.125,
            true,
            &mut pc,
            Point3dClassifyOptions::default(),
            &crate::topology::brep_cache::arena_identity,
            &|_| None,
        );
        assert!(matches!(result, Err(SolidError::InvalidInput { .. })));
        assert_eq!(pc.face(), Some(EntityRef::Face(face)));
        assert_eq!(pc.point_class(), PointClassType::Unknown);
        assert_eq!(
            pc.gap_3d(&kernel, &BrepGeometry { b: &brep }, None)
                .unwrap(),
            17.0
        );
    }
}

#[cfg(test)]
mod retained_point3d_cache_tests {
    use super::*;
    use crate::geom::tree::{NodeClass, NodeData, NodeId};
    use crate::query::face::p2::{FacePointCache, FacePointCacheDecision, MappedFacePointCache};
    use crate::surface::plane::Plane;
    use crate::surface::surface_cache::p1::{FaceContainment, SurfaceCacheOptions};
    use crate::surface::trim_srf_cache::TrimSrfCache;

    #[derive(Clone, Copy, Debug, PartialEq)]
    enum Stored {
        Surface,
        Topology(EntityRef),
    }
    fn decode(handle: Stored) -> Result<EntityRef> {
        match handle {
            Stored::Topology(entity) => Ok(entity),
            Stored::Surface => Err(SolidError::invalid_input(
                OP_POINT,
                "surface is not topology",
            )),
        }
    }

    #[test]
    fn point3d_finish_consumes_the_original_warm_tree_instead_of_loop_fallback() {
        let kernel = Kernel::default();
        let mut b = Brep::new();
        let shell = b.add_shell().unwrap();
        let face = b.add_face(shell).unwrap();
        let domain = Extent2d::from_bounds(-1.0, -1.0, 1.0, 1.0);
        let mut surface = Surface::from(Plane::new(
            Point3::ZERO,
            Vec3::X,
            Vec3::Y,
            crate::math::Vec2::new(1.0, 1.0),
            &domain,
        ));
        let keys = |object| {
            let mut key = crate::topology::brep_cache::arena_identity(object);
            key.key.0 += 0x70000;
            key
        };
        surface.set_owner(Some(keys(CacheObject::Face(face))));
        let cache_surface = surface.clone_for_evaluation().unwrap();
        b.face_mut(face).unwrap().surface = Some(surface);
        b.face_mut(face).unwrap().uv_domain = domain;
        b.face_mut(face).unwrap().tolerance = 0.125;
        b.set_editing_enabled(true);
        let mut cache = TrimSrfCache::new(cache_surface, SurfaceCacheOptions::default());
        cache.base_mut().set_surface_handle(Some(Stored::Surface));
        cache.base_mut().have_t_surface_cache = true;
        cache.base_mut().face_containment_done = FaceContainment::Done;
        let encode = |entity| Ok(Stored::Topology(entity));
        {
            let mut view = MappedFacePointCache::new(&mut cache, &decode, &encode);
            assert_eq!(
                view.classify_ready(&kernel, &mut b, face, Point2::ZERO, &keys, &|_| None)
                    .unwrap(),
                FacePointCacheDecision::NoDecision
            );
        }
        // Real built cache storage, preloaded with a retained inside classification.
        // The fallback has no loops and therefore cannot independently return Face.
        let tree = cache.base_mut().tree_mut().unwrap();
        for index in 0..tree.nodes().len() {
            let NodeData::Surface(aux) = &mut tree.node_mut(NodeId(index)).unwrap().data else {
                panic!("surface cache node");
            };
            aux.node_class = NodeClass::Inside;
            aux.face = Some(Stored::Topology(EntityRef::Face(face)));
        }
        let node_count = cache.base().tree().unwrap().nodes().len();
        cache.base_mut().have_t_surface_cache = true;
        cache.base_mut().face_containment_done = FaceContainment::Done;
        cache.base_mut().face_was_modified = false;
        let mut cached = PointClassification::new(0.125);
        cached.set_face(Some(EntityRef::Face(face)));
        let mut view = MappedFacePointCache::new(&mut cache, &decode, &encode);
        point_3d_classify_finish_impl(
            &kernel,
            &mut b,
            face,
            Point2::ZERO,
            Point3::ZERO,
            0.125,
            false,
            &mut cached,
            Point3dClassifyOptions::default(),
            Some(PointFinishContext {
                keys: &keys,
                owner_of: &|_| None,
                trim_cache: Some(&mut view),
            }),
        )
        .unwrap();
        assert_eq!(cached.point_class(), PointClassType::Face);
        assert_eq!(cached.object(), Some(EntityRef::Face(face)));
        assert_eq!(cached.face(), Some(EntityRef::Face(face)));
        assert_eq!(cached.uv_param(), Point2::ZERO);
        assert_eq!(cache.base().tree().unwrap().nodes().len(), node_count);
        assert_eq!(cache.base().surface_handle(), Some(Stored::Surface));
        assert!(!cache.base().face_was_modified);
        let mut no_cache = PointClassification::new(0.125);
        point_3d_classify_finish(
            &kernel,
            &mut b,
            face,
            Point2::ZERO,
            Point3::ZERO,
            0.125,
            false,
            &mut no_cache,
            Point3dClassifyOptions::default(),
        )
        .unwrap();
        assert_eq!(no_cache.point_class(), PointClassType::Unknown);
        assert!(b.editing_enabled());
    }
}
