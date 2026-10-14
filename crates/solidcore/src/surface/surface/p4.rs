//! `SmSurface`, part 4 (work in progress).

use crate::curve::curve::Curve;
use crate::error::{Result, SolidError};
use crate::geom::global_solver::SolverOperation;
use crate::geom::local_solve_nd::{compute_jacobian, EvalNFunctions};
use crate::geom::matrix::Matrix;
use crate::math::{Point2, Vec3, EFF_ZERO, EFF_ZERO_SQ};
use crate::nurbs::knots::Side;

use super::p7::{evaluate_1st_derivatives, evaluate_2nd_derivatives};
use crate::surface::plane::Plane;

use super::Surface;

const OP_CS_INTERSECT: &str = "surface::surface::p4::CsIntersect::evaluate";
const OP_CS_NORMALIZE: &str = "surface::surface::p4::CsNormalize::evaluate";
const OP_SSP_INTERSECT: &str = "surface::surface::p4::SspIntersect::evaluate";
const OP_SSI_NORMALIZE: &str = "surface::surface::p4::SsiNormalize::evaluate";
const OP_SS_NORMALIZE: &str = "surface::surface::p4::SsNormalize::evaluate";

/// Entry `i` of an evaluation result.
fn nth(v: &[Vec3], i: usize, op: &'static str) -> Result<Vec3> {
    v.get(i)
        .copied()
        .ok_or_else(|| SolidError::invalid_input(op, "evaluation result too short"))
}

/// The `x` and `f` slices must hold `n` entries.
fn check_len(x: &[f64], f: &[f64], n: usize, op: &'static str) -> Result<()> {
    if x.len() == n && f.len() == n {
        Ok(())
    } else {
        Err(SolidError::invalid_input(
            op,
            "parameter or function array has the wrong length",
        ))
    }
}

/// Writes `values` into the leading entries of `f`.
fn set_all(f: &mut [f64], values: &[f64]) {
    for (slot, value) in f.iter_mut().zip(values) {
        *slot = *value;
    }
}

/// Newton equations of a curve/surface intersection (`SmCSIntersectENFO`):
/// `F = S(u, v) - C(t)`; after slow convergence near a tangency it switches to the
/// normalize equations. Holds the call count, last gap and the switch flag.
#[derive(Debug)]
pub struct CsIntersect<'a> {
    curve: &'a Curve,
    surface: &'a Surface,
    count: u64,
    last_sum: f64,
    do_minimize: bool,
}

impl<'a> CsIntersect<'a> {
    /// The equations for `curve` against `surface` (the constructor; the squared
    /// tolerance argument is unused upstream).
    #[must_use]
    pub const fn new(curve: &'a Curve, surface: &'a Surface) -> Self {
        Self {
            curve,
            surface,
            count: 0,
            last_sum: 0.0,
            do_minimize: false,
        }
    }
}

impl EvalNFunctions for CsIntersect<'_> {
    /// `Evaluate`: `x = [t u v]`; never reports an answer found.
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        mut jacobian: Option<&mut Matrix>,
    ) -> Result<bool> {
        check_len(x, f, 3, OP_CS_INTERSECT)?;
        let (Some(&t), Some(&u), Some(&v)) = (x.first(), x.get(1), x.get(2)) else {
            return Err(SolidError::sm_err(OP_CS_INTERSECT));
        };

        // count calls
        self.count += 1;

        // Evaluate CurvePoint
        let crv = self.curve.evaluate(t, 1, true, true)?;
        let c = nth(&crv, 0, OP_CS_INTERSECT)?;
        let ct = nth(&crv, 1, OP_CS_INTERSECT)?;

        // Evaluate SurfacePoint
        let uv = Point2::new(u, v);
        let (mut s, mut du, mut dv) =
            evaluate_1st_derivatives(self.surface, uv, Side::Left, Side::Left)?;

        // Fail at singular Curve and Surface Points - zero tangent lengths
        if du.length_squared() < EFF_ZERO_SQ
            || dv.length_squared() < EFF_ZERO_SQ
            || ct.length_squared() < EFF_ZERO_SQ
        {
            return Err(SolidError::sm_err(OP_CS_INTERSECT));
        }

        // get Gap Vector = S - C
        let mut gap: Vec3 = s - c;

        if let Some(jac) = jacobian.as_deref_mut() {
            jac.set_at(0, 0, -ct.x)?;
            jac.set_at(0, 1, du.x)?;
            jac.set_at(0, 2, dv.x)?;

            jac.set_at(1, 0, -ct.y)?;
            jac.set_at(1, 1, du.y)?;
            jac.set_at(1, 2, dv.y)?;

            jac.set_at(2, 0, -ct.z)?;
            jac.set_at(2, 1, du.z)?;
            jac.set_at(2, 2, dv.z)?;
        }

        // Compute function values
        set_all(f, &[gap.x, gap.y, gap.z]);

        // Check for solver not converging: too many iterations near a tangency.
        let iter_lim: u64 = 10;
        let tangent_tol_deg = 0.5;
        let scaled_zero = EFF_ZERO * (1.0 + s.max_dimension());

        if self.do_minimize || (self.count > iter_lim && gap.length_squared() > self.last_sum) {
            // if we are near a tangent point
            let crv2 = self.curve.evaluate(t, 2, true, true)?;
            let norm = self.surface.evaluate_normal(uv, Side::Left, Side::Left)?;
            let c0 = nth(&crv2, 0, OP_CS_INTERSECT)?;
            let c1 = nth(&crv2, 1, OP_CS_INTERSECT)?;
            if norm.is_perpendicular_to(c1, tangent_tol_deg) {
                // Slow convergence problem.  Switch to a normalize operation.
                self.do_minimize = true;
                let d2 = evaluate_2nd_derivatives(self.surface, uv, Side::Left, Side::Left)?;
                s = d2.point;
                du = d2.du;
                dv = d2.dv;
                let (duv, duu, dvv) = (d2.duv, d2.duu, d2.dvv);
                let c2 = nth(&crv2, 2, OP_CS_INTERSECT)?;

                gap = s - c0;

                // translate the surface point by its normal when the gap is tiny
                if gap.length_squared() < 1000.0 * scaled_zero {
                    if norm.dot(gap) >= 0.0 {
                        gap = (s + norm) - c0;
                    } else {
                        gap = (s - norm) - c0;
                    }
                }

                if let Some(jac) = jacobian.as_deref_mut() {
                    jac.set_at(0, 0, c2.dot(gap) - c1.dot(c1))?;
                    jac.set_at(0, 1, c1.dot(du))?;
                    jac.set_at(0, 2, c1.dot(dv))?;

                    jac.set_at(1, 0, -du.dot(c1))?;
                    jac.set_at(1, 1, duu.dot(gap) + du.dot(du))?;
                    let j12 = duv.dot(gap) + du.dot(dv);
                    jac.set_at(1, 2, j12)?;

                    jac.set_at(2, 0, -dv.dot(c1))?;
                    jac.set_at(2, 1, j12)?;
                    jac.set_at(2, 2, dvv.dot(gap) + dv.dot(dv))?;
                }

                set_all(f, &[c1.dot(gap), du.dot(gap), dv.dot(gap)]);
            }
        }

        if self.count > 9 {
            self.last_sum = gap.length_squared();
        }

        Ok(false)
    }
}

/// Newton equations of the curve/surface normalize, minimize and maximize solves, and
/// of the 3D signed directed minimize (`SmCSNormalizeENFO`).
#[derive(Debug)]
pub struct CsNormalize<'a> {
    curve: &'a Curve,
    surface: &'a Surface,
    operation: SolverOperation,
    opt_vectors: &'a [Vec3],
}

impl<'a> CsNormalize<'a> {
    /// The equations for `curve` against `surface` (the constructor; the distance
    /// and squared tolerances are unused upstream). `opt_vectors[0]` is the
    /// direction of a signed directed minimize.
    #[must_use]
    pub const fn new(
        curve: &'a Curve,
        surface: &'a Surface,
        operation: SolverOperation,
        opt_vectors: &'a [Vec3],
    ) -> Self {
        Self {
            curve,
            surface,
            operation,
            opt_vectors,
        }
    }
}

impl EvalNFunctions for CsNormalize<'_> {
    /// `Evaluate`: `x = [t u v]` (and `w` for the signed directed minimize); never
    /// reports an answer found.
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        mut jacobian: Option<&mut Matrix>,
    ) -> Result<bool> {
        let directed = self.operation == SolverOperation::Signed3dDirectedMinimize;
        let count = if directed { 4 } else { 3 };
        check_len(x, f, count, OP_CS_NORMALIZE)?;
        let (Some(&t), Some(&u), Some(&v)) = (x.first(), x.get(1), x.get(2)) else {
            return Err(SolidError::sm_err(OP_CS_NORMALIZE));
        };

        // evaluate Curve Point and derivs
        let crv = self.curve.evaluate(t, 2, true, true)?;
        let c = nth(&crv, 0, OP_CS_NORMALIZE)?;
        let ct = nth(&crv, 1, OP_CS_NORMALIZE)?;
        let ctt = nth(&crv, 2, OP_CS_NORMALIZE)?;

        // evaluate Surface Point and derivs
        let uv = Point2::new(u, v);
        let d2 = evaluate_2nd_derivatives(self.surface, uv, Side::Left, Side::Left)?;
        let (s, su, sv, suu, suv, svv) = (d2.point, d2.du, d2.dv, d2.duu, d2.duv, d2.dvv);

        // Fail at singular curve or surface points - zero length tangents
        if ct.length_squared() < EFF_ZERO_SQ
            || su.length_squared() < EFF_ZERO_SQ
            || sv.length_squared() < EFF_ZERO_SQ
        {
            return Err(SolidError::sm_err(OP_CS_NORMALIZE));
        }

        // get gap vector
        let mut gap: Vec3 = s - c;

        let scaled_zero = EFF_ZERO * 1000.0 * (1.0 + s.max_dimension());

        // Translate the surface point by its normal when the gap gets very small.
        if gap.length_squared() < scaled_zero {
            let norm = if directed {
                nth(self.opt_vectors, 0, OP_CS_NORMALIZE)?
            } else {
                self.surface.evaluate_normal(uv, Side::Left, Side::Left)?
            };
            if norm.dot(gap) >= 0.0 {
                gap = (s + norm) - c;
            } else {
                gap = (s - norm) - c;
            }
        }

        // Unitize Ct and get its derivative.
        let ct_unit = ct.unitize_or_z();
        let ct_unit_t = ct.unitized_derivative(ctt);

        if directed {
            if self.opt_vectors.is_empty() {
                return Err(SolidError::null_pointer(OP_CS_NORMALIZE));
            }
            let dir = nth(self.opt_vectors, 0, OP_CS_NORMALIZE)?;

            // Surface Normal, unit, and derivatives
            let cross = su.cross(sv);
            let cross_u = su.derivative_of_cross_product(sv, suu, suv);
            let cross_v = su.derivative_of_cross_product(sv, suv, svv);

            let n = cross.unitize_or_z();
            let n_u = cross.unitized_derivative(cross_u);
            let n_v = cross.unitized_derivative(cross_v);

            if let Some(jac) = jacobian.as_deref_mut() {
                jac.set_at(0, 0, ct_unit_t.dot(n))?;
                jac.set_at(0, 1, ct_unit.dot(n_u))?;
                jac.set_at(0, 2, ct_unit.dot(n_v))?;
                jac.set_at(0, 3, 0.0)?;

                jac.set_at(1, 0, -ct.x)?;
                jac.set_at(1, 1, su.x)?;
                jac.set_at(1, 2, sv.x)?;
                jac.set_at(1, 3, -dir.x)?;

                jac.set_at(2, 0, -ct.y)?;
                jac.set_at(2, 1, su.y)?;
                jac.set_at(2, 2, sv.y)?;
                jac.set_at(2, 3, -dir.y)?;

                jac.set_at(3, 0, -ct.z)?;
                jac.set_at(3, 1, su.z)?;
                jac.set_at(3, 2, sv.z)?;
                jac.set_at(3, 3, -dir.z)?;
            }

            // Compute function values
            let w = x.get(3).copied().unwrap_or(0.0);
            let vec: Vec3 = (s - c) - w * dir;

            set_all(f, &[ct_unit.dot(n), vec.x, vec.y, vec.z]);
        } else {
            // Normalization, Minimization, Maximization: unitized Ct, Su and Sv.
            let su_unit = su.unitize_or_z();
            let su_unit_u = su.unitized_derivative(suu);
            let su_unit_v = su.unitized_derivative(suv);
            let sv_unit = sv.unitize_or_z();
            let sv_unit_u = sv.unitized_derivative(suv);
            let sv_unit_v = sv.unitized_derivative(svv);

            if let Some(jac) = jacobian.as_deref_mut() {
                jac.set_at(0, 0, ct_unit_t.dot(gap) - ct_unit.dot(ct))?;
                jac.set_at(0, 1, ct_unit.dot(su))?;
                jac.set_at(0, 2, ct_unit.dot(sv))?;

                jac.set_at(1, 0, -su_unit.dot(ct))?;
                jac.set_at(1, 1, su_unit_u.dot(gap) + su_unit.dot(su))?;
                jac.set_at(1, 2, su_unit_v.dot(gap) + su_unit.dot(sv))?;

                jac.set_at(2, 0, -sv_unit.dot(ct))?;
                jac.set_at(2, 1, sv_unit_u.dot(gap) + sv_unit.dot(su))?;
                jac.set_at(2, 2, sv_unit_v.dot(gap) + sv_unit.dot(sv))?;
            }

            set_all(f, &[ct_unit.dot(gap), su_unit.dot(gap), sv_unit.dot(gap)]);
        }

        Ok(false)
    }
}

/// Newton equations of the surface/surface/plane intersection by tangent planes
/// (`SmSSPIntersectENFO`): `x = [u1 v1 u2 v2]`. With `steepest_descent` the plane of
/// each evaluation passes through the middle of the two surface points with the
/// cross product of the surface normals as its normal; otherwise it is `plane`.
#[derive(Debug)]
pub struct SspIntersect<'a> {
    surface: &'a Surface,
    other_surface: &'a Surface,
    steepest_descent: bool,
    plane: Option<&'a Plane>,
}

impl<'a> SspIntersect<'a> {
    /// The equations for `surface` against `other_surface` (the constructor; the
    /// squared tolerance is unused upstream).
    #[must_use]
    pub const fn new(
        surface: &'a Surface,
        other_surface: &'a Surface,
        steepest_descent: bool,
        plane: Option<&'a Plane>,
    ) -> Self {
        Self {
            surface,
            other_surface,
            steepest_descent,
            plane,
        }
    }
}

impl EvalNFunctions for SspIntersect<'_> {
    /// `Evaluate`: true when both surface points lie on the plane and coincide.
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        mut jacobian: Option<&mut Matrix>,
    ) -> Result<bool> {
        check_len(x, f, 4, OP_SSP_INTERSECT)?;
        let (Some(&u1), Some(&v1), Some(&u2), Some(&v2)) =
            (x.first(), x.get(1), x.get(2), x.get(3))
        else {
            return Err(SolidError::sm_err(OP_SSP_INTERSECT));
        };

        // get surface 3d points - with 1st derivatives
        let uv = Point2::new(u1, v1);
        let uv2 = Point2::new(u2, v2);
        let (s_f, duf, dvf) = evaluate_1st_derivatives(self.surface, uv, Side::Left, Side::Left)?;
        let (s_g, dug, dvg) =
            evaluate_1st_derivatives(self.other_surface, uv2, Side::Left, Side::Left)?;

        // pick the plane of this evaluation
        let temp_plane;
        let plane: &Plane = if self.steepest_descent {
            // plane = plane(surf1Point, surf2Point)
            let through = (s_f + s_g) / 2.0;
            let mut norm1 = duf.cross(dvf);
            let mut norm2 = dug.cross(dvg);
            if norm1.length_squared() < EFF_ZERO {
                norm1 = self.surface.evaluate_normal(uv, Side::Left, Side::Left)?;
            }
            if norm2.length_squared() < EFF_ZERO {
                norm2 = self
                    .other_surface
                    .evaluate_normal(uv2, Side::Left, Side::Left)?;
            }
            let plane_norm = norm1.cross(norm2);
            if plane_norm.length_squared() < EFF_ZERO_SQ {
                return Err(SolidError::sm_err(OP_SSP_INTERSECT));
            }
            temp_plane = Plane::from_normal(through, plane_norm.unitize_or_z());
            &temp_plane
        } else {
            self.plane
                .ok_or_else(|| SolidError::null_pointer(OP_SSP_INTERSECT))?
        };

        // evaluate Plane derivatives
        let pe = plane.evaluate(Point2::new(0.0, 0.0), 1, 1, true)?;
        let h = nth(
            pe.first().map_or(&[][..], Vec::as_slice),
            0,
            OP_SSP_INTERSECT,
        )?;
        let duh = nth(
            pe.get(1).map_or(&[][..], Vec::as_slice),
            0,
            OP_SSP_INTERSECT,
        )?;
        let dvh = nth(
            pe.first().map_or(&[][..], Vec::as_slice),
            1,
            OP_SSP_INTERSECT,
        )?;

        // Utilize plane equation aX + bY + cZ + d = 0 in solver
        let abc = duh.cross(dvh).unitize()?;
        let d = -abc.dot(h);

        // if we have a zero length tangent then we must fail
        if duf.length_squared() < EFF_ZERO_SQ
            || dvf.length_squared() < EFF_ZERO_SQ
            || dug.length_squared() < EFF_ZERO_SQ
            || dvg.length_squared() < EFF_ZERO_SQ
        {
            return Err(SolidError::sm_err(OP_SSP_INTERSECT));
        }

        // get gap vector
        let diff: Vec3 = s_f - s_g;

        // evaluate F[3]
        let f3 = abc.dot(s_f) + d + abc.dot(s_g) + d;
        if let Some(slot) = f.get_mut(3) {
            *slot = f3;
        }

        // exit criteria - Surface points on plane and close to each other
        if diff.length_squared() < EFF_ZERO_SQ && f3 < EFF_ZERO {
            return Ok(true);
        }

        if let Some(jac) = jacobian.as_deref_mut() {
            jac.set_at(0, 0, duf.x)?;
            jac.set_at(0, 1, dvf.x)?;
            jac.set_at(0, 2, -dug.x)?;
            jac.set_at(0, 3, -dvg.x)?;

            jac.set_at(1, 0, duf.y)?;
            jac.set_at(1, 1, dvf.y)?;
            jac.set_at(1, 2, -dug.y)?;
            jac.set_at(1, 3, -dvg.y)?;

            jac.set_at(2, 0, duf.z)?;
            jac.set_at(2, 1, dvf.z)?;
            jac.set_at(2, 2, -dug.z)?;
            jac.set_at(2, 3, -dvg.z)?;

            jac.set_at(3, 0, abc.dot(duf))?;
            jac.set_at(3, 1, abc.dot(dvf))?;
            jac.set_at(3, 2, abc.dot(dug))?;
            jac.set_at(3, 3, abc.dot(dvg))?;
        }

        // Compute function values (F[3] above)
        for (slot, value) in f.iter_mut().zip([diff.x, diff.y, diff.z]) {
            *slot = value;
        }

        Ok(false)
    }
}

/// Newton equations of the surface/surface intersection on a given plane
/// (`SmSSINormalizeENFO`): `x = [u1 v1 u2 v2]`, the plane `N . p + D = 0`.
#[derive(Debug)]
pub struct SsiNormalize<'a> {
    surface: &'a Surface,
    other_surface: &'a Surface,
    plane_normal: Vec3,
    plane_d: f64,
}

impl<'a> SsiNormalize<'a> {
    /// The equations for `surface` against `other_surface` on the plane with normal
    /// `plane_normal` and offset `plane_d` (the constructor; the squared tolerance
    /// is unused upstream).
    #[must_use]
    pub const fn new(
        surface: &'a Surface,
        other_surface: &'a Surface,
        plane_normal: Vec3,
        plane_d: f64,
    ) -> Self {
        Self {
            surface,
            other_surface,
            plane_normal,
            plane_d,
        }
    }
}

impl EvalNFunctions for SsiNormalize<'_> {
    /// `Evaluate`: never reports an answer found; the Jacobian is the finite
    /// difference one (upstream replaces its analytic one).
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
    ) -> Result<bool> {
        check_len(x, f, 4, OP_SSI_NORMALIZE)?;
        let (Some(&u1), Some(&v1), Some(&u2), Some(&v2)) =
            (x.first(), x.get(1), x.get(2), x.get(3))
        else {
            return Err(SolidError::sm_err(OP_SSI_NORMALIZE));
        };

        let uv = Point2::new(u1, v1);
        let uv2 = Point2::new(u2, v2);
        let plane_norm = self.plane_normal;

        // get surface 3d points - with 2nd derivatives
        let d1 = evaluate_2nd_derivatives(self.surface, uv, Side::Left, Side::Left)?;
        let d2 = evaluate_2nd_derivatives(self.other_surface, uv2, Side::Left, Side::Left)?;
        let mut srf_point = d1.point;
        let srf_point2 = d2.point;
        let (du, dv) = (d1.du, d1.dv);
        let (du2, dv2) = (d2.du, d2.dv);

        // if we have a zero length tangent then we must fail
        if du.length_squared() < EFF_ZERO_SQ
            || dv.length_squared() < EFF_ZERO_SQ
            || du2.length_squared() < EFF_ZERO_SQ
            || dv2.length_squared() < EFF_ZERO_SQ
        {
            return Err(SolidError::sm_err(OP_SSI_NORMALIZE));
        }

        // when the gap is small
        let mut gap: Vec3 = srf_point2 - srf_point;
        let tol = EFF_ZERO * 1000.0 * (1.0 + srf_point.max_dimension());
        if gap.length_squared() < tol {
            // translate one surface point by its normal and solve a little better
            let norm = self.surface.evaluate_normal(uv, Side::Left, Side::Left)?;
            // (upstream evaluates the second normal on the *first* surface)
            let norm2 = self.surface.evaluate_normal(uv2, Side::Left, Side::Left)?;
            let ave_norm = if norm.dot(norm2) < 0.0 {
                (norm - norm2) / 2.0
            } else {
                (norm + norm2) / 2.0
            };
            // project AveNorm to Plane
            let move_in_plane = plane_norm.cross(ave_norm).cross(plane_norm);
            if norm.dot(gap) > 0.0 {
                srf_point = srf_point - move_in_plane;
                gap = srf_point2 - srf_point;
            } else {
                srf_point = srf_point + move_in_plane;
                gap = srf_point2 - srf_point;
            }
        }

        // get surface normals
        let n = du.cross(dv);
        let n2 = du2.cross(dv2);

        // unitize surf1 Normal
        let mut nu = n;
        if nu.length_squared() < EFF_ZERO_SQ {
            nu = self.surface.evaluate_normal(uv, Side::Left, Side::Left)?;
            if nu.length_squared() < EFF_ZERO_SQ {
                return Err(SolidError::sm_err(OP_SSI_NORMALIZE));
            }
        }
        let nu = nu.unitize_or_z();

        // unitize surf2 Normal (upstream evaluates it on the first surface)
        let mut n2u = n2;
        if n2u.length_squared() < EFF_ZERO_SQ {
            n2u = self.surface.evaluate_normal(uv2, Side::Left, Side::Left)?;
            if n2u.length_squared() < EFF_ZERO_SQ {
                return Err(SolidError::sm_err(OP_SSI_NORMALIZE));
            }
        }
        let n2u = n2u.unitize_or_z();

        // Get two plane vectors equal to surf/plane intersection directions
        let v_1 = nu.cross(plane_norm);
        let v_2 = n2u.cross(plane_norm);

        // Compute function values
        set_all(
            f,
            &[
                v_1.dot(gap),
                v_2.dot(gap),
                srf_point.dot(plane_norm) + self.plane_d,
                srf_point2.dot(plane_norm) + self.plane_d,
            ],
        );

        if let Some(jac) = jacobian {
            // The analytic Jacobian upstream builds is overwritten by the finite
            // difference one over the natural domains of both surfaces.
            let dom = self.surface.natural_uv_domain();
            let dom2 = self.other_surface.natural_uv_domain();
            let domains = [
                dom.u_interval(),
                dom.v_interval(),
                dom2.u_interval(),
                dom2.v_interval(),
            ];
            let mut tmp_x = x.to_vec();
            let f_now = f.to_vec();
            compute_jacobian(self, &mut tmp_x, &f_now, jac, Some(&domains))?;
        }

        Ok(false)
    }
}

/// Newton equations of the surface/surface normalize, minimize and maximize solves,
/// and of the 3D signed directed minimize (`SmSSNormalizeENFO`): `x = [u1 v1 u2 v2]`
/// (and `w` for the signed directed minimize).
#[derive(Debug)]
pub struct SsNormalize<'a> {
    surface: &'a Surface,
    other_surface: &'a Surface,
    operation: SolverOperation,
    opt_vectors: &'a [Vec3],
}

impl<'a> SsNormalize<'a> {
    /// The equations for `surface` against `other_surface` (the constructor; the
    /// squared tolerance is unused upstream). `opt_vectors[0]` is the direction of
    /// a signed directed minimize.
    #[must_use]
    pub const fn new(
        surface: &'a Surface,
        other_surface: &'a Surface,
        operation: SolverOperation,
        opt_vectors: &'a [Vec3],
    ) -> Self {
        Self {
            surface,
            other_surface,
            operation,
            opt_vectors,
        }
    }
}

impl EvalNFunctions for SsNormalize<'_> {
    /// `Evaluate`: never reports an answer found.
    fn evaluate(
        &mut self,
        x: &[f64],
        f: &mut [f64],
        mut jacobian: Option<&mut Matrix>,
    ) -> Result<bool> {
        let directed = self.operation == SolverOperation::Signed3dDirectedMinimize;
        let count = if directed { 5 } else { 4 };
        check_len(x, f, count, OP_SS_NORMALIZE)?;
        let (Some(&u1), Some(&v1), Some(&u2), Some(&v2)) =
            (x.first(), x.get(1), x.get(2), x.get(3))
        else {
            return Err(SolidError::sm_err(OP_SS_NORMALIZE));
        };

        // Surf1Point and Surf2Point current solution values
        let uv1 = Point2::new(u1, v1);
        let uv2 = Point2::new(u2, v2);
        let d1 = evaluate_2nd_derivatives(self.surface, uv1, Side::Left, Side::Left)?;
        let d2 = evaluate_2nd_derivatives(self.other_surface, uv2, Side::Left, Side::Left)?;
        let (p1, s1u, s1v, s1uu, s1uv, s1vv) = (d1.point, d1.du, d1.dv, d1.duu, d1.duv, d1.dvv);
        let (p2, s2u, s2v, s2uu, s2uv, s2vv) = (d2.point, d2.du, d2.dv, d2.duu, d2.duv, d2.dvv);

        // if we have a zero length tangent then we must fail
        if s1u.length_squared() < EFF_ZERO_SQ
            || s1v.length_squared() < EFF_ZERO_SQ
            || s2u.length_squared() < EFF_ZERO_SQ
            || s2v.length_squared() < EFF_ZERO_SQ
        {
            return Err(SolidError::sm_err(OP_SS_NORMALIZE));
        }

        // when the gap is small
        let mut gap: Vec3 = p2 - p1;
        let tol = EFF_ZERO * 1000.0 * (1.0 + p1.max_dimension());
        if gap.length_squared() < tol {
            // translate one surface point by its normal and solve a little better
            let norm = if directed {
                nth(self.opt_vectors, 0, OP_SS_NORMALIZE)?
            } else {
                self.surface.evaluate_normal(uv1, Side::Left, Side::Left)?
            };
            gap = if norm.dot(gap) > 0.0 {
                p2 - (p1 - norm)
            } else {
                p2 - (p1 + norm)
            };
        }

        // Unitize S1u and S1v and get their derivatives.
        let s1u_unit = s1u.unitize_or_z();
        let s1u_unit_u = s1u.unitized_derivative(s1uu);
        let s1u_unit_v = s1u.unitized_derivative(s1uv);
        let s1v_unit = s1v.unitize_or_z();
        let s1v_unit_u = s1v.unitized_derivative(s1uv); // not the same as s1u_unit_v
        let s1v_unit_v = s1v.unitized_derivative(s1vv);

        if directed {
            if self.opt_vectors.is_empty() {
                return Err(SolidError::null_pointer(OP_SS_NORMALIZE));
            }
            let dir = nth(self.opt_vectors, 0, OP_SS_NORMALIZE)?;

            // Surf2 unit Normal and Normal derivatives
            let cross2 = s2u.cross(s2v);
            let cross2_u = s2u.derivative_of_cross_product(s2v, s2uu, s2uv);
            let cross2_v = s2u.derivative_of_cross_product(s2v, s2uv, s2vv);

            let n2 = cross2.unitize_or_z();
            let n2_u = cross2.unitized_derivative(cross2_u);
            let n2_v = cross2.unitized_derivative(cross2_v);

            if let Some(jac) = jacobian.as_deref_mut() {
                jac.set_at(0, 0, s1u_unit_u.dot(n2))?;
                jac.set_at(0, 1, s1u_unit_v.dot(n2))?;
                jac.set_at(0, 2, s1u_unit.dot(n2_u))?;
                jac.set_at(0, 3, s1u_unit.dot(n2_v))?;
                jac.set_at(0, 4, 0.0)?;

                jac.set_at(1, 0, s1v_unit_u.dot(n2))?;
                jac.set_at(1, 1, s1v_unit_v.dot(n2))?;
                jac.set_at(1, 2, s1v_unit.dot(n2_u))?;
                jac.set_at(1, 3, s1v_unit.dot(n2_v))?;
                jac.set_at(1, 4, 0.0)?;

                jac.set_at(2, 0, s1u.x)?;
                jac.set_at(2, 1, s1v.x)?;
                jac.set_at(2, 2, -s2u.x)?;
                jac.set_at(2, 3, -s2v.x)?;
                jac.set_at(2, 4, -dir.x)?;

                jac.set_at(3, 0, s1u.y)?;
                jac.set_at(3, 1, s1v.y)?;
                jac.set_at(3, 2, -s2u.y)?;
                jac.set_at(3, 3, -s2v.y)?;
                jac.set_at(3, 4, -dir.y)?;

                jac.set_at(4, 0, s1u.z)?;
                jac.set_at(4, 1, s1v.z)?;
                jac.set_at(4, 2, -s2u.z)?;
                jac.set_at(4, 3, -s2v.z)?;
                jac.set_at(4, 4, -dir.z)?;
            }

            // Compute function values
            let w = x.get(4).copied().unwrap_or(0.0);
            let vec: Vec3 = (p1 - p2) - w * dir;
            set_all(
                f,
                &[s1u_unit.dot(n2), s1v_unit.dot(n2), vec.x, vec.y, vec.z],
            );
        } else {
            // Unitize S2u and S2v and get their derivatives.
            let s2u_unit = s2u.unitize_or_z();
            let s2u_unit_u = s2u.unitized_derivative(s2uu);
            let s2u_unit_v = s2u.unitized_derivative(s2uv);
            let s2v_unit = s2v.unitize_or_z();
            let s2v_unit_u = s2v.unitized_derivative(s2uv); // not the same as s2u_unit_v
            let s2v_unit_v = s2v.unitized_derivative(s2vv);

            if let Some(jac) = jacobian.as_deref_mut() {
                jac.set_at(0, 0, s1u_unit_u.dot(gap) - s1u_unit.dot(s1u))?;
                jac.set_at(0, 1, s1u_unit_v.dot(gap) - s1u_unit.dot(s1v))?;
                jac.set_at(0, 2, s1u_unit.dot(s2u))?;
                jac.set_at(0, 3, s1u_unit.dot(s2v))?;

                jac.set_at(1, 0, s1v_unit_u.dot(gap) - s1v_unit.dot(s1u))?;
                jac.set_at(1, 1, s1v_unit_v.dot(gap) - s1v_unit.dot(s1v))?;
                jac.set_at(1, 2, s1v_unit.dot(s2u))?;
                jac.set_at(1, 3, s1v_unit.dot(s2v))?;

                jac.set_at(2, 0, -s2u_unit.dot(s1u))?;
                jac.set_at(2, 1, -s2u_unit.dot(s1v))?;
                jac.set_at(2, 2, s2u_unit_u.dot(gap) + s2u_unit.dot(s2u))?;
                jac.set_at(2, 3, s2u_unit_v.dot(gap) + s2u_unit.dot(s2v))?;

                jac.set_at(3, 0, -s2v_unit.dot(s1u))?;
                jac.set_at(3, 1, -s2v_unit.dot(s1v))?;
                jac.set_at(3, 2, s2v_unit_u.dot(gap) + s2v_unit.dot(s2v))?;
                jac.set_at(3, 3, s2v_unit_v.dot(gap) + s2v_unit.dot(s2v))?;
            }

            // Compute function values
            set_all(
                f,
                &[
                    s1u_unit.dot(gap),
                    s1v_unit.dot(gap),
                    s2u_unit.dot(gap),
                    s2v_unit.dot(gap),
                ],
            );
        }

        Ok(false)
    }
}
