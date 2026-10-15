//! B-rep queries and checked surface evaluation entry points.
//!
//! Upstream: `source/SM_API/src/SmApiSurfaces.cpp`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`evaluate_surface_point`] | `SmApiEvaluateSurfacePoint` |
//! | [`evaluate_surface_normal`] | `SmApiEvaluateSurfaceNormal` |
//! | [`evaluate_surface_derivatives`] | `SmApiEvaluateSurfaceDerivatives` |
//!
//! Deviations: a nullable surface is `Option<&Surface>`; output parameters become
//! return values, so failure cannot modify caller-owned outputs. The remaining
//! functions assigned to query.brep are not implemented yet.

use crate::geom::vector2d::Point2;
use crate::geom::vector3d::{Point3, Vec3};
use crate::nurbs::knots::Side;
use crate::surface::surface::Surface;
use crate::{Result, SolidError};

fn checked_surface<'a>(
    surface: Option<&'a Surface>,
    uv: Point2,
    op: &'static str,
) -> Result<&'a Surface> {
    let surface = surface.ok_or_else(|| SolidError::invalid_input(op, "null surface"))?;
    let domain = surface.natural_uv_domain();
    let lo = domain.min();
    let hi = domain.max();
    // Keep the negated conjunction: a NaN parameter must fail the API range test.
    if !(uv.x >= lo.x && uv.x <= hi.x && uv.y >= lo.y && uv.y <= hi.y) {
        return Err(SolidError::invalid_input(
            op,
            "parameter outside natural domain",
        ));
    }
    Ok(surface)
}

/// Evaluate a surface point after the API's inclusive natural-domain check.
///
/// # Errors
/// Returns invalid input for a null surface or an out-of-domain (including NaN)
/// parameter, otherwise propagates the surface evaluation error unchanged.
pub fn evaluate_surface_point(surface: Option<&Surface>, uv: Point2) -> Result<Point3> {
    checked_surface(surface, uv, "query::brep::evaluate_surface_point")?.evaluate_point(uv)
}

/// Evaluate the unit normal, with both derivative sides taken from the left.
///
/// # Errors
/// Those of [`evaluate_surface_point`] and the surface's normal evaluation.
pub fn evaluate_surface_normal(surface: Option<&Surface>, uv: Point2) -> Result<Vec3> {
    checked_surface(surface, uv, "query::brep::evaluate_surface_normal")?.evaluate_normal(
        uv,
        Side::Left,
        Side::Left,
    )
}

/// Evaluate `(point, du, dv)`, with both derivative sides taken from the left.
///
/// # Errors
/// Those of [`evaluate_surface_point`] and the surface's first derivatives.
pub fn evaluate_surface_derivatives(
    surface: Option<&Surface>,
    uv: Point2,
) -> Result<(Point3, Vec3, Vec3)> {
    checked_surface(surface, uv, "query::brep::evaluate_surface_derivatives")?
        .evaluate_1st_derivatives(uv, true, true)
}
