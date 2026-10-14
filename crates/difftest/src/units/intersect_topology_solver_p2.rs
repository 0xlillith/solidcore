//! Differential coverage of topology p2's local tree traversal.
//!
//! The shared p1 fixture adapter builds the same objects and trees on both sides;
//! the p2 C++ entry reaches the pinned `SmTopologySolver::LocalSolve`.
use super::intersect_topology_solver_p1 as shared;
use crate::{Case, CaseContext, HarnessError, Output, Rng, UnitDiff};
use serde_json::Value;
use solidcore_ffi::units::intersect_topology_solver_p2 as ffi;

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
        shared::Runner.run_rust(case)
    }
    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        shared::p2_tree_generate(rng)
    }
}
