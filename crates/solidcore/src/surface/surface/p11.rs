//! `SmSurface`, part 11: the local surface/surface intersection from a pair of guess
//! points.
//!
//! Upstream: `source/SMLib/src/SmSurface.cpp`, `source/SMLib/inc/SmSurface.h`.
//!
//! | Rust                                              | Upstream                                         |
//! |---------------------------------------------------|--------------------------------------------------|
//! | [`Surface::local_surface_intersect`], [`SurfaceIntersection`], [`LocalSurfaceIntersectOptions`] | `SmSurface::LocalSurfaceIntersect` (not virtual) |
//!
//! Deviations:
//! - **Out-params** become return values: `LocalSurfaceIntersect` returns `None` for
//!   `rbFoundAnswer == FALSE` and otherwise the solution it writes (`m_vStart[0..4]`
//!   and its value, the distance between the two surface points) as a
//!   [`SurfaceIntersection`]; the solution's objects (this surface and the other one,
//!   upstream pointers) are the caller's to attach.
//! - **Recursion.** `LocalSurfaceIntersect` calls itself once, with the iteration
//!   solution as guess, when that solution reached a domain boundary and no boundary
//!   was tried (the second call then starts on that boundary). The port carries a
//!   depth and returns `SM_ERR` past [`MAX_RECURSION`] calls (upstream reaches depth 1).
//! - **Statuses.** Every status upstream checks (`SER`) is propagated; the two
//!   `Unitize` calls of the plane normal ignore theirs (a zero normal becomes
//!   `(0, 0, 1)`, [`Vec3::unitize_or_z`]).
//! - **Unported callees** (return [`SolidError::Unported`]): the file-static
//!   `sm_IntersectSrfBndryCrv` (boundary iso curve against the other surface, only with
//!   `do_boundary_testing`; unit `surface.surface_p8`) and `SmSurface::LocalSurfaceSolve`
//!   (the tangent-touching refinements when the normals at the solution are within 3
//!   degrees of parallel; unit `surface.surface_p3`).
//! - Debug-only work (`SM_ASSERT`) is not ported.

use crate::error::{Result, SolidError};
use crate::geom::extent1d::{Extent1d, ExtentNd};
use crate::geom::extent2d::Extent2d;
use crate::geom::global_solver::SolverOperation;
use crate::geom::local_solve_nd::{BoundaryHandler, LocalSolveNd, TerminationReason};
use crate::geom::tol::are_same;
use crate::geom::vector2d::Point2;
use crate::geom::vector3d::{Point3, Vec3};
use crate::math::{EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI};
use crate::nurbs::b_spline_surface::SurfParamType;
use crate::nurbs::knots::Side;
use crate::surface::plane::Plane;

use super::p1::{SsiNormalizeEnfo, SspIntersectEnfo};
use super::Surface;

const OP_LSI: &str = "surface::surface::p11::Surface::local_surface_intersect";

/// The deepest `LocalSurfaceIntersect` self call accepted (upstream reaches 1).
pub const MAX_RECURSION: usize = 16;

/// Three degrees, in radians (`3.0 * SM_PI / 180.0`).
fn three_degrees() -> f64 {
    3.0 * PI / 180.0
}

// ---- calls into units that are not ported yet ------------------------------------------------

/// What `sm_IntersectSrfBndryCrv` finds: the surface point and the iso curve
/// parameter (the deviation and tangency flag are not read by this unit).
#[derive(Debug, Clone, Copy)]
struct BoundaryHit {
    uv: Point2,
    t: f64,
}

/// `sm_IntersectSrfBndryCrv(surface, domain, curve_surface, .., param, iso, tol,
/// guess_uv, guess_t, ..)`: the intersection of `surface` with the iso parameter curve
/// of `curve_surface`; `None` when not found.
#[allow(clippy::too_many_arguments)]
fn intersect_srf_bndry_crv(
    surface: &Surface,
    uv_domain: &Extent2d,
    curve_surface: &Surface,
    param: SurfParamType,
    iso_parameter: f64,
    tol_3d: f64,
    guess_uv: Point2,
    guess_t: f64,
) -> Result<Option<BoundaryHit>> {
    let _ = (
        surface,
        uv_domain,
        curve_surface,
        param,
        iso_parameter,
        tol_3d,
        guess_uv,
        guess_t,
    );
    // UNPORTED(surface.surface_p8): sm_IntersectSrfBndryCrv
    Err(SolidError::unported(
        "sm_IntersectSrfBndryCrv",
        "surface.surface_p8",
    ))
}

/// `SmSurface::LocalSurfaceSolve(domain, other, other_domain, op, tol, NULL, NULL,
/// uv, other_uv, found, solution)`: `[u, v, other u, other v]` and the solution value;
/// `None` when not found.
#[allow(clippy::too_many_arguments)]
fn local_surface_solve(
    surface: &Surface,
    uv_domain: &Extent2d,
    other: &Surface,
    other_uv_domain: &Extent2d,
    op: SolverOperation,
    tol: f64,
    uv: Point2,
    other_uv: Point2,
) -> Result<Option<([f64; 4], f64)>> {
    let _ = (surface, uv_domain, other, other_uv_domain, op, tol, uv, other_uv);
    // UNPORTED(surface.surface_p3): SmSurface::LocalSurfaceSolve
    Err(SolidError::unported(
        "SmSurface::LocalSurfaceSolve",
        "surface.surface_p3",
    ))
}

// ---- LocalSurfaceIntersect -------------------------------------------------------------------

/// A surface/surface intersection point found by [`Surface::local_surface_intersect`]
/// (upstream's `rSolution`: four variables, single valued).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceIntersection {
    /// `m_vStart[0..2]`: the point on this surface.
    pub uv: Point2,
    /// `m_vStart[2..4]`: the point on the other surface.
    pub other_uv: Point2,
    /// `m_vStart.m_dSolutionValue`: the distance between the two surface points.
    pub deviation: f64,
}

/// The optional arguments of [`Surface::local_surface_intersect`].
#[derive(Debug, Clone, Copy, Default)]
pub struct LocalSurfaceIntersectOptions {
    /// `cpOptPlaneNormal`: the solution lies on the plane through the guess points
    /// with this normal (else the cross product of the surface normals, recomputed
    /// while iterating: steepest descent).
    pub plane_normal: Option<Vec3>,
    /// `cpOptPreviousUV`: the previous point on this surface when stepping, so a
    /// boundary solution right next to it is not preferred.
    pub previous_uv: Option<Point2>,
    /// `bDoBoundaryTesting`: look for intersections with the surfaces' boundaries.
    pub do_boundary_testing: bool,
}

/// `sSolutionVector[..]` as `Point2`s.
fn split(sv: &[f64; 4]) -> (Point2, Point2) {
    (Point2::new(sv[0], sv[1]), Point2::new(sv[2], sv[3]))
}

/// The four parameters of a local solve's solution vector.
fn four(x: &[f64]) -> Result<[f64; 4]> {
    match x {
        [a, b, c, d, ..] => Ok([*a, *b, *c, *d]),
        _ => Err(SolidError::sm_err(OP_LSI)),
    }
}

/// The intervals of the two UV domains (`SmExtentNd(4)`).
fn intervals(uv_domain: &Extent2d, other_uv_domain: &Extent2d) -> ExtentNd {
    let mut ivls = ExtentNd::new(4);
    let values = [
        Extent1d::from_min_max(uv_domain.min().x, uv_domain.max().x),
        Extent1d::from_min_max(uv_domain.min().y, uv_domain.max().y),
        Extent1d::from_min_max(other_uv_domain.min().x, other_uv_domain.max().x),
        Extent1d::from_min_max(other_uv_domain.min().y, other_uv_domain.max().y),
    ];
    for (i, ivl) in values.into_iter().enumerate() {
        if let Some(slot) = ivls.get_mut(i) {
            *slot = ivl;
        }
    }
    ivls
}

impl Surface {
    /// The surface/surface intersection point near the guesses `uv_guess` (on this
    /// surface over `uv_domain`) and `other_uv_guess` (on `other` over
    /// `other_uv_domain`) (`LocalSurfaceIntersect`): a Newton solve of the two
    /// surfaces and a plane (through the guess points, normal to the given normal or to
    /// the cross product of the surface normals), with boundary iso curve solves when
    /// testing boundaries, tangent-touching refinements when the normals are nearly
    /// parallel, and a normalized surface/surface/plane solve when the solution leaves
    /// a gap. `None` when no point within `distance_tolerance` was found (also when the
    /// guess normals are parallel and no plane normal is given).
    ///
    /// Dispatch: `EvaluatePoint`/`EvaluateNormal` are the virtuals of each surface's
    /// own kind ([`Surface::evaluate_point`], [`Surface::evaluate_normal`]).
    ///
    /// # Errors
    ///
    /// Those of the evaluations, the angle between zero normals and the solves (all
    /// checked upstream); `SM_ERR` past [`MAX_RECURSION`] self calls;
    /// [`SolidError::Unported`]: `sm_IntersectSrfBndryCrv` (`surface.surface_p8`) and
    /// `SmSurface::LocalSurfaceSolve` (`surface.surface_p3`).
    pub fn local_surface_intersect(
        &self,
        uv_domain: &Extent2d,
        other: &Surface,
        other_uv_domain: &Extent2d,
        distance_tolerance: f64,
        uv_guess: Point2,
        other_uv_guess: Point2,
        opts: &LocalSurfaceIntersectOptions,
    ) -> Result<Option<SurfaceIntersection>> {
        self.local_surface_intersect_at(
            uv_domain,
            other,
            other_uv_domain,
            distance_tolerance,
            uv_guess,
            other_uv_guess,
            opts,
            0,
        )
    }

    #[allow(clippy::too_many_lines)]
    fn local_surface_intersect_at(
        &self,
        uv_domain: &Extent2d,
        other: &Surface,
        other_uv_domain: &Extent2d,
        distance_tolerance: f64,
        uv_guess: Point2,
        other_uv_guess: Point2,
        opts: &LocalSurfaceIntersectOptions,
        depth: usize,
    ) -> Result<Option<SurfaceIntersection>> {
        if depth > MAX_RECURSION {
            return Err(SolidError::sm_err(OP_LSI));
        }
        let do_boundary = opts.do_boundary_testing;
        let three_deg = three_degrees();

        // locals for NewtonRaphson Solutions
        let mut found_solution = false;

        // locals for boundary solutions
        let mut tried_boundary = false;
        let mut this_on_boundary = false;
        let mut other_on_boundary = false;
        let mut found_boundary = false;
        let mut uv_boundary = Point2::UNDEF;
        let mut uv_other_boundary = Point2::UNDEF;

        // 1st - When testing boundaries: when guess point is on a natural boundary look
        //       for 1st boundaryCurve/Surface intersection
        if do_boundary {
            // when guessPoint is on Surface1 U Boundary
            if are_same(uv_guess.x, uv_domain.min().x) || are_same(uv_guess.x, uv_domain.max().x)
            {
                tried_boundary = true;
                this_on_boundary = true;

                // intersect OtherSurface with const u Surface boundary isoParameterCurve
                // (SER: checked; sm_IntersectSrfBndryCrv over the other surface's kind)
                let hit = intersect_srf_bndry_crv(
                    other,
                    other_uv_domain,
                    self,
                    SurfParamType::U,
                    uv_guess.x,
                    distance_tolerance,
                    other_uv_guess,
                    uv_guess.y,
                )?;
                found_solution = hit.is_some();

                // save found solutions
                if let Some(h) = hit {
                    found_boundary = true;
                    uv_other_boundary = h.uv;
                    uv_boundary = Point2::new(uv_guess.x, h.t);
                }
            }

            // when guessPoint is on Surface1 V Boundary
            if !found_solution
                && (are_same(uv_guess.y, uv_domain.min().y)
                    || are_same(uv_guess.y, uv_domain.max().y))
            {
                tried_boundary = true;
                this_on_boundary = true;

                // intersect OtherSurface with const v Surface boundary isoParameterCurve
                // (SER: checked)
                let hit = intersect_srf_bndry_crv(
                    other,
                    other_uv_domain,
                    self,
                    SurfParamType::V,
                    uv_guess.y,
                    distance_tolerance,
                    other_uv_guess,
                    uv_guess.x,
                )?;
                found_solution = hit.is_some();
                if let Some(h) = hit {
                    found_boundary = true;
                    uv_other_boundary = h.uv;
                    uv_boundary = Point2::new(h.t, uv_guess.y);
                }
            }

            // when guessPoint is on Surface2 U Boundary
            if !found_solution
                && (are_same(other_uv_guess.x, other_uv_domain.min().x)
                    || are_same(other_uv_guess.x, other_uv_domain.max().x))
            {
                other_on_boundary = true;
                tried_boundary = true;

                // intersect Surface with const u OtherSurface boundary isoParameterCurve
                // (SER: checked)
                let hit = intersect_srf_bndry_crv(
                    self,
                    uv_domain,
                    other,
                    SurfParamType::U,
                    other_uv_guess.x,
                    distance_tolerance,
                    uv_guess,
                    other_uv_guess.y,
                )?;
                found_solution = hit.is_some();
                if let Some(h) = hit {
                    found_boundary = true;
                    uv_boundary = h.uv;
                    uv_other_boundary = Point2::new(other_uv_guess.x, h.t);
                }
            }

            // when guessPoint is on Surface2 V Boundary
            if !found_solution
                && (are_same(other_uv_guess.y, other_uv_domain.min().y)
                    || are_same(other_uv_guess.y, other_uv_domain.max().y))
            {
                other_on_boundary = true;
                tried_boundary = true;

                // intersect Surface with const v OtherSurface boundary isoParameterCurve
                // (SER: checked)
                let hit = intersect_srf_bndry_crv(
                    self,
                    uv_domain,
                    other,
                    SurfParamType::V,
                    other_uv_guess.y,
                    distance_tolerance,
                    uv_guess,
                    other_uv_guess.x,
                )?;
                // (its found flag is overwritten by the solve below before any read)
                if let Some(h) = hit {
                    found_boundary = true;
                    uv_boundary = h.uv;
                    uv_other_boundary = Point2::new(h.t, other_uv_guess.y);
                }
            }
        }

        // 2nd - iterate to Surf/Surf intersection Point from GuessPoint
        let ivls = intervals(uv_domain, other_uv_domain);
        // Do not do periodicity for this function.
        let periodicities = [false; 4];
        let guess = [uv_guess.x, uv_guess.y, other_uv_guess.x, other_uv_guess.y];

        // get Surface GuessPoints (SER: checked; each surface's own EvaluatePoint)
        let pnt1 = self.evaluate_point(uv_guess)?;
        let pnt2 = other.evaluate_point(other_uv_guess)?;

        // Select through point for plane containing answer
        let plane_orig: Point3 = if this_on_boundary {
            pnt1
        } else if other_on_boundary {
            pnt2
        } else {
            (pnt1 + pnt2) / 2.0
        };

        // Select normal for plane containing answer
        let steepest_descent;
        let mut plane_norm = if let Some(n) = opts.plane_normal {
            // In some cases like tracing tangency curves we need to get a normal from
            // outside because of the problems with the surface normals being the same.
            steepest_descent = false;
            n
        } else {
            // (SER: checked; EvaluateNormal(uv, TRUE, TRUE) of each surface's kind)
            let norm1 = self.evaluate_normal(uv_guess, Side::Left, Side::Left)?;
            let norm2 = other.evaluate_normal(other_uv_guess, Side::Left, Side::Left)?;
            steepest_descent = true;
            let n = norm1.cross(norm2);

            // check for bad starting points
            if n.length_squared() < EFF_ZERO_SQ {
                // Bad starting point - just quit and try a different one
                return Ok(None);
            }
            n
        };

        // build plane containing answers (SE: Unitize status ignored)
        plane_norm = plane_norm.unitize_or_z();
        let plane_surface = Plane::from_normal(plane_orig, plane_norm);

        // instance SmSSPIntersectENFO surf/surf/plane NewtonRaphson Solver
        let ssi_node_tol = distance_tolerance / 100.0;
        let mut eval_fun = SspIntersectEnfo::new(
            self,
            other,
            steepest_descent,
            Some(&plane_surface),
            ssi_node_tol * ssi_node_tol,
        );
        let (ls_found, ls_x, termination) = {
            let mut ls = LocalSolveNd::new(&mut eval_fun, Some(&ivls), Some(&periodicities));
            ls.set_desired_accuracy(EFF_ZERO_SQ);
            ls.set_boundary_handler(BoundaryHandler::TotalBoundaryHits, 3);
            ls.set_maximum_iterations(20);

            // Solve at 4 times tolerance for now and cull solutions to tolerance later
            // (SER: checked)
            let r = ls.solve_it(&guess, distance_tolerance * 4.0)?;
            (r.found, r.x, ls.termination_reason())
        };
        found_solution = ls_found;
        let mut sv = four(&ls_x)?;

        // set solution points
        let (mut uv, mut uv_other) = split(&sv);

        // When Solution iterated to a boundary from guessPoints not on the boundary and
        // we are doing BoundaryTesting
        if do_boundary
            && found_solution
            && !tried_boundary
            && (uv_domain.is_point_on_boundary(uv, EFF_ZERO).0
                || other_uv_domain.is_point_on_boundary(uv_other, EFF_ZERO).0)
        {
            // do intersection over with guessPoint = iteration Solution so that
            // BoundarySolutions can be found and selected.
            return self.local_surface_intersect_at(
                uv_domain,
                other,
                other_uv_domain,
                distance_tolerance,
                uv,
                uv_other,
                opts,
                depth + 1,
            );
        }

        // 3rd - Select Good Boundary Solutions over Iteration Solutions when we have them
        if found_boundary {
            let mut use_boundary = true;

            // skip boundary solutions that are very close to the last solution in favor
            // of iterationSolutions that aren't.
            if let (true, Some(prev)) = (found_solution, opts.previous_uv) {
                let uv_pt = Point2::new(sv[0], sv[1]);
                let dist_interior = prev.distance(uv_pt);
                let dist_boundary = prev.distance(uv_boundary);
                if dist_interior > 100.0 * dist_boundary {
                    use_boundary = false;
                }
            }

            // save the boundarySolution
            if use_boundary {
                found_solution = true;
                sv = [
                    uv_boundary.x,
                    uv_boundary.y,
                    uv_other_boundary.x,
                    uv_other_boundary.y,
                ];
                (uv, uv_other) = split(&sv);
            }
        }

        // 4th - When no Iteration solution was found because boundaries were hit - look
        //       for a surf/BoundaryCurve solution at the iteration stopping point
        if !found_solution && termination == TerminationReason::OutOfBounds && do_boundary {
            // when solution search stopped on Surface1 U Boundary
            if are_same(sv[0], uv_domain.min().x) || are_same(sv[0], uv_domain.max().x) {
                // (SER: checked)
                let hit = intersect_srf_bndry_crv(
                    other,
                    other_uv_domain,
                    self,
                    SurfParamType::U,
                    sv[0],
                    distance_tolerance,
                    Point2::new(sv[2], sv[3]),
                    sv[1],
                )?;
                found_solution = hit.is_some();
                if let Some(h) = hit {
                    sv = [sv[0], h.t, h.uv.x, h.uv.y];
                    (uv, uv_other) = split(&sv);
                }
            }

            // when solution search stopped on Surface1 V Boundary (upstream's
            // precedence: `(!found && min) || max`)
            if (!found_solution && are_same(sv[1], uv_domain.min().y))
                || are_same(sv[1], uv_domain.max().y)
            {
                // (SER: checked)
                let hit = intersect_srf_bndry_crv(
                    other,
                    other_uv_domain,
                    self,
                    SurfParamType::V,
                    sv[1],
                    distance_tolerance,
                    Point2::new(sv[2], sv[3]),
                    sv[0],
                )?;
                found_solution = hit.is_some();
                if let Some(h) = hit {
                    sv = [h.t, sv[1], h.uv.x, h.uv.y];
                    (uv, uv_other) = split(&sv);
                }
            }

            // when solution search stopped on Surface2 U Boundary
            if (!found_solution && are_same(sv[2], other_uv_domain.min().x))
                || are_same(sv[2], other_uv_domain.max().x)
            {
                // (SER: checked)
                let hit = intersect_srf_bndry_crv(
                    self,
                    uv_domain,
                    other,
                    SurfParamType::U,
                    sv[2],
                    distance_tolerance,
                    Point2::new(sv[0], sv[1]),
                    sv[3],
                )?;
                found_solution = hit.is_some();
                if let Some(h) = hit {
                    sv = [h.uv.x, h.uv.y, sv[2], h.t];
                    (uv, uv_other) = split(&sv);
                }
            }

            // when solution search stopped on Surface2 V Boundary
            if (!found_solution && are_same(sv[3], other_uv_domain.min().y))
                || are_same(sv[3], other_uv_domain.max().y)
            {
                // (SER: checked)
                let hit = intersect_srf_bndry_crv(
                    self,
                    uv_domain,
                    other,
                    SurfParamType::V,
                    sv[3],
                    distance_tolerance,
                    Point2::new(sv[0], sv[1]),
                    sv[2],
                )?;
                found_solution = hit.is_some();
                if let Some(h) = hit {
                    sv = [h.uv.x, h.uv.y, h.t, sv[3]];
                    (uv, uv_other) = split(&sv);
                }
            }
        }

        // 5th - when a solution was found and bSteepestDescent is TRUE and the surface
        //       normals are parallel to one another - look for a solution where the
        //       surfaces are a minDist from one another. A TangentTouching case
        if found_solution && steepest_descent {
            // get angle between surface normals at the solution (SER: checked)
            let norm1 = self.evaluate_normal(uv, Side::Left, Side::Left)?;
            let norm2 = other.evaluate_normal(uv_other, Side::Left, Side::Left)?;
            let angle = norm1.angle_between(norm2)?;

            // If the angle is within 3 degrees of being parallel then try refining
            // solution with a normalized convergence.
            if angle < three_deg || angle > PI - three_deg {
                // find point of closest approach between the two surfaces in current
                // solution neighborhood (SER: checked)
                let sol_norm = local_surface_solve(
                    self,
                    uv_domain,
                    other,
                    other_uv_domain,
                    SolverOperation::Minimize,
                    distance_tolerance,
                    uv,
                    uv_other,
                )?;
                // (SER: checked)
                let p1 = self.evaluate_point(uv)?;
                let p2 = other.evaluate_point(uv_other)?;
                let distance_between_points = p1.distance(p2);

                // use new solution when deviation is better
                if let Some((x, value)) = sol_norm {
                    if value < distance_between_points {
                        sv = x;
                        (uv, uv_other) = split(&sv);
                    }
                }
            }
        }

        // 6th - when a Non-Zero Deviation solution has been found
        if found_solution {
            // (SER: checked)
            let p1 = self.evaluate_point(uv)?;
            let p2 = other.evaluate_point(uv_other)?;
            let mut dist = p1.distance(p2);

            // If we don't converge to an exact answer then try converging to
            // surf/surf/plane intersection with a variation on the set of equations
            if dist > EFF_ZERO {
                // (SER: checked)
                let norm1 = self.evaluate_normal(uv, Side::Left, Side::Left)?;
                let norm2 = other.evaluate_normal(uv_other, Side::Left, Side::Left)?;
                let angle = norm1.angle_between(norm2)?;

                // If the angle is within 3 degrees of being perpendicular than we should
                // try to do a normalized convergence for minimization.
                if angle < three_deg || angle > PI - three_deg {
                    // Recompute plane normal if needed
                    if opts.plane_normal.is_none() {
                        let new_plane_norm = norm1.cross(norm2);
                        if new_plane_norm.length_squared() > EFF_ZERO_SQ {
                            plane_norm = new_plane_norm;
                        }
                    }
                    // (SE: Unitize status ignored)
                    plane_norm = plane_norm.unitize_or_z();
                    let plane_d = -p1.dot(plane_norm);

                    let ivls2 = intervals(uv_domain, other_uv_domain);
                    let periodicities2 = [false; 4];
                    let guess2 = [uv.x, uv.y, uv_other.x, uv_other.y];

                    // set up surf/surf/plane intersection NewtonRaphson solver using
                    // different eqns
                    let mut eval_fun2 =
                        SsiNormalizeEnfo::new(self, other, EFF_ZERO * EFF_ZERO, plane_norm, plane_d);
                    let r2 = {
                        let mut ls2 =
                            LocalSolveNd::new(&mut eval_fun2, Some(&ivls2), Some(&periodicities2));
                        ls2.set_boundary_handler(BoundaryHandler::TotalBoundaryHits, 8);

                        // solve it (SER: checked)
                        ls2.solve_it(&guess2, EFF_ZERO_SQRT)?
                    };

                    // when an answer was found
                    if r2.found {
                        let sv2 = four(&r2.x)?;
                        let (uv2, uv_other2) = split(&sv2);
                        // (SER: checked)
                        let pt1 = self.evaluate_point(uv2)?;
                        let pt2 = other.evaluate_point(uv_other2)?;
                        let dist2 = pt1.distance(pt2);

                        // save solution when deviation is better than previous solution
                        // deviation
                        if dist2 < dist {
                            dist = dist2;
                            sv = sv2;
                            (uv, uv_other) = split(&sv);
                        }
                    }

                    // when distance is not exactly zero
                    if dist != 0.0 {
                        // find point of closest approach between the two surfaces in
                        // current solution neighborhood (SER: checked; the found flag
                        // overwrites bFoundSolution upstream)
                        let sol = local_surface_solve(
                            self,
                            uv_domain,
                            other,
                            other_uv_domain,
                            SolverOperation::Normalize,
                            distance_tolerance,
                            uv,
                            uv_other,
                        )?;
                        if let Some((x, _)) = sol {
                            // get the Surf1Point/Surf2Point gap size (SER: checked)
                            let (uv2, uv_other2) = split(&x);
                            let pt1 = self.evaluate_point(uv2)?;
                            let pt2 = other.evaluate_point(uv_other2)?;
                            let dist2 = pt1.distance(pt2);

                            // save solution when this is a better answer
                            if dist2 < dist {
                                dist = dist2;
                                sv = x;
                            }
                        }
                    }
                }
            }

            // when deviation is within tolerance - save output
            if dist < distance_tolerance {
                let (uv, other_uv) = split(&sv);
                return Ok(Some(SurfaceIntersection {
                    uv,
                    other_uv,
                    deviation: dist,
                }));
            }
        }

        // all done
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plane(origin: Point3, normal: Vec3) -> Surface {
        Surface::from(Plane::from_normal(origin, normal))
    }

    #[test]
    fn two_perpendicular_planes_meet_on_their_line() {
        let a = plane(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
        let b = plane(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let dom = Extent2d::from_points(Point2::new(-5.0, -5.0), Point2::new(5.0, 5.0));
        let r = a
            .local_surface_intersect(
                &dom,
                &b,
                &dom,
                1e-6,
                Point2::new(0.3, 0.2),
                Point2::new(0.1, 0.4),
                &LocalSurfaceIntersectOptions::default(),
            )
            .unwrap()
            .unwrap();
        let p = a.evaluate_point(r.uv).unwrap();
        let q = b.evaluate_point(r.other_uv).unwrap();
        assert!(p.distance(q) < 1e-6);
        assert!(p.x.abs() < 1e-9 && p.z.abs() < 1e-9);
    }

    #[test]
    fn parallel_guess_normals_give_no_answer() {
        let a = plane(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
        let b = plane(Point3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 0.0, 1.0));
        let dom = Extent2d::from_points(Point2::new(-5.0, -5.0), Point2::new(5.0, 5.0));
        let r = a
            .local_surface_intersect(
                &dom,
                &b,
                &dom,
                1e-6,
                Point2::new(0.3, 0.2),
                Point2::new(0.1, 0.4),
                &LocalSurfaceIntersectOptions::default(),
            )
            .unwrap();
        assert!(r.is_none());
    }
}
