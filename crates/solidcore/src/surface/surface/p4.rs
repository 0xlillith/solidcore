//! `SmSurface`, part 4 (work in progress).

use crate::curve::curve::Curve;
use crate::error::{Result, SolidError};
use crate::geom::global_solver::SolverOperation;
use crate::geom::local_solve_nd::EvalNFunctions;
use crate::geom::matrix::Matrix;
use crate::math::{Point2, Vec3, EFF_ZERO, EFF_ZERO_SQ};
use crate::nurbs::knots::Side;

use super::p7::{evaluate_1st_derivatives, evaluate_2nd_derivatives};
use super::Surface;

const OP_CS_INTERSECT: &str = "surface::surface::p4::CsIntersect::evaluate";
const OP_CS_NORMALIZE: &str = "surface::surface::p4::CsNormalize::evaluate";

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
