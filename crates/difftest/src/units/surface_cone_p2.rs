//! Differential coverage of real cone ellipse trims; setup failures fail cases.
use crate::{Case, HarnessError, Output, UnitDiff};
use serde_json::Value;
use solidcore::curve::ellipse::{Ellipse, EllipseOptions};
use solidcore::geom::{extent1d::Extent1d, extent2d::Extent2d};
use solidcore::math::{Point3, Vec3};
use solidcore::surface::cone::p1::Cone;
use solidcore::surface::cone::p2::ellipse_trim_intervals;
use solidcore_ffi::units::surface_cone_p2 as ffi;

/// Cone part-2 differential runner.
pub struct Runner;

fn input(case: &Case) -> Result<Vec<f64>, HarnessError> {
    let values: Vec<f64> =
        serde_json::from_value(case.inputs.get("trim").cloned().unwrap_or(Value::Null))
            .map_err(|e| HarnessError::new(e.to_string()))?;
    ffi::validate(&values).map_err(|e| HarnessError::new(e.to_string()))?;
    Ok(values)
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "surface.cone_p2"
    }
    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        ffi::eval_trim(&input(case)?)
            .map(|v| Output::new().values("intervals", v))
            .map_err(|e| HarnessError::new(e.to_string()))
    }
    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let values = input(case)?;
        let [bottom, top, height, start, end, swap, inside, tilt, z, rx, ry, a, b, u0, v0, u1, v1, tol] =
            values.as_slice()
        else {
            return Err(HarnessError::new("invalid cone trim record"));
        };
        let kernel = solidcore::Kernel::default();
        let mut cone = Cone::new_bounded(
            Point3::ZERO,
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            *bottom,
            *top,
            *start,
            *end,
            *height,
            *swap != 0.0,
            *inside != 0.0,
            false,
        );
        cone.make_nurb(&kernel)
            .map_err(|e| HarnessError::new(e.to_string()))?;
        let ellipse = Ellipse::new(
            &kernel,
            Point3::new(0.0, 0.0, *z),
            Vec3::new(tilt.cos(), 0.0, tilt.sin()),
            Vec3::new(0.0, 1.0, 0.0),
            Extent1d::from_min_max(0.0, 360.0),
            *rx,
            *ry,
            &EllipseOptions {
                dim: 3,
                ..EllipseOptions::default()
            },
        )
        .map_err(|e| HarnessError::new(e.to_string()))?;
        let intervals = ellipse_trim_intervals(
            &kernel,
            &ellipse,
            &Extent1d::from_min_max(*a, *b),
            &cone,
            &Extent2d::from_bounds(*u0, *v0, *u1, *v1),
            *tol,
        )
        .map_err(|e| HarnessError::new(e.to_string()))?;
        let mut out = vec![intervals.len() as f64];
        for interval in intervals {
            out.push(interval.min());
            out.push(interval.max());
        }
        Ok(Output::new().values("intervals", out))
    }
}
