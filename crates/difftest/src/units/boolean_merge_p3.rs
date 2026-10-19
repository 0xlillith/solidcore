//! Differential tests for merge orchestration. Leaf-only trees move stable
//! handles into an immutable input pool; each handle identifies a complete real
//! native fixture imported before running the port. No merge result is imported.
//! Geometric merges remain deferred at the exact missing dependency call.
use crate::{Case, CaseContext, HarnessError, Output, Rng, UnitDiff};
use port::{BooleanDriver, PiecewiseAccess};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use solidcore::boolean::merge::{p1::MergeOptions, p3 as port};
use solidcore::geom::a_object::BooleanOperation;
use solidcore::geom::context::MarkLock;
use solidcore::kernel::Kernel;
use solidcore::topology::{brep::p1 as bp, Brep};
use solidcore::SolidError;
use solidcore_ffi::units::boolean_merge_p3 as ffi;

/// The five merge functions of `boolean.merge_p3`.
pub struct Runner;
#[derive(Debug, Clone, Serialize, Deserialize)]
struct Operand {
    point: [f64; 3],
    #[serde(default)]
    size: Option<[f64; 3]>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
struct Inputs {
    family: u8,
    operation: i32,
    #[serde(default)]
    options: Option<u8>,
    #[serde(default)]
    a: Vec<Operand>,
    #[serde(default)]
    b: Vec<Operand>,
    #[serde(default)]
    tree: Vec<i64>,
}
fn harness(error: &(impl ToString + ?Sized)) -> HarnessError {
    HarnessError::new(error.to_string())
}
fn encode(i: &Inputs) -> Vec<f64> {
    let mut data = vec![
        f64::from(i.family),
        f64::from(i.operation),
        f64::from(i.options.unwrap_or(0)),
        f64::from(u8::from(i.options.is_some())),
        i.a.len() as f64,
        i.b.len() as f64,
    ];
    for operand in i.a.iter().chain(&i.b) {
        data.push(f64::from(u8::from(operand.size.is_some())));
        data.extend(operand.point);
        if let Some(size) = operand.size {
            data.extend(size);
        }
    }
    data.push(i.tree.len() as f64);
    data.extend(i.tree.iter().map(|&v| v as f64));
    data
}
fn options(mask: u8) -> MergeOptions {
    MergeOptions::new(
        mask & 1 != 0,
        mask & 2 != 0,
        mask & 4 != 0,
        mask & 8 != 0,
        mask & 16 != 0,
    )
}
fn status(error: &SolidError) -> i64 {
    match error {
        SolidError::InvalidInput { .. } => 1007,
        SolidError::Upstream { code, .. } => i64::from(*code),
        _ => -1,
    }
}
fn output(raw: &[f64], family: u8) -> Result<Output, HarnessError> {
    let (&status, rest) = raw.split_first().ok_or_else(|| harness("empty result"))?;
    if status != 0.0 {
        return Ok(Output::err(format!("status={status}")));
    }
    let mut identities = Vec::new();
    if family == 0 {
        let (&count, mut remaining) = rest
            .split_first()
            .ok_or_else(|| harness("missing tree count"))?;
        for _ in 0..count as usize {
            let (&id, next) = remaining
                .split_first()
                .ok_or_else(|| harness("missing identity"))?;
            let (&length, next) = next
                .split_first()
                .ok_or_else(|| harness("missing graph size"))?;
            identities.push(id as i64);
            remaining = next
                .get(length as usize..)
                .ok_or_else(|| harness("truncated graph"))?;
        }
        if !remaining.is_empty() {
            return Err(harness("trailing tree output"));
        }
    }
    Ok(Output::new()
        .int("status", status as i64)
        .ints("operand_identity", identities)
        .values("record", rest.to_vec()))
}
struct Pool {
    breps: Vec<Brep>,
    raw: Vec<Vec<f64>>,
    primary: Option<usize>,
    other: Option<usize>,
    manifold: bool,
}
impl Pool {
    fn load(i: &Inputs) -> Result<Self, HarnessError> {
        let mut breps = Vec::new();
        let mut raw = Vec::new();
        for operand in i.a.iter().chain(&i.b) {
            let exported = ffi::eval(&encode(&Inputs {
                family: 5,
                operation: 0,
                options: None,
                a: vec![operand.clone()],
                b: Vec::new(),
                tree: Vec::new(),
            }))
            .map_err(|error| harness(&error))?;
            // export_graph starts at validate_only; vertexuse import expects a
            // leading success tag. Geometry comes only from fixture construction.
            let mut tagged = vec![0.0];
            tagged.extend(&exported);
            let (brep, _) = super::topology_vertexuse::import(&tagged)?;
            breps.push(brep);
            raw.push(exported);
        }
        Ok(Self {
            breps,
            raw,
            primary: Some(0),
            other: Some(i.a.len()),
            manifold: true,
        })
    }
    fn record(&self, id: usize, values: &mut Vec<f64>) -> solidcore::Result<()> {
        let raw = self
            .raw
            .get(id)
            .ok_or_else(|| SolidError::invalid_input("harness::merge_p3", "bad pool handle"))?;
        values.extend([id as f64, raw.len() as f64]);
        values.extend(raw);
        Ok(())
    }
}
impl BooleanDriver for Pool {
    type Operand = usize;
    fn same_operand(&self, left: &usize, right: &usize) -> bool {
        left == right
    }
    fn construct(
        &mut self,
        _k: &mut Kernel,
        _primary: &mut Option<usize>,
        _other: &mut Option<usize>,
        _approx_tol: f64,
        _angle_tol: f64,
        _options: Option<&MergeOptions>,
    ) -> solidcore::Result<()> {
        // UNPORTED(boolean.merge_p2): SmMerge::SmMerge
        Err(SolidError::unported("SmMerge::SmMerge", "boolean.merge_p2"))
    }
    fn manifold_boolean(
        &mut self,
        _k: &mut Kernel,
        _primary: &mut Option<usize>,
        _other: &mut Option<usize>,
        _manifold: &mut bool,
        _lamina: &mut u64,
        _operation: BooleanOperation,
        _lock: Option<&mut MarkLock>,
    ) -> solidcore::Result<Option<usize>> {
        // UNPORTED(boolean.merge_p4): SmMerge::ManifoldBoolean
        Err(SolidError::unported(
            "SmMerge::ManifoldBoolean",
            "boolean.merge_p4",
        ))
    }
}
fn absent<T>() -> solidcore::Result<T> {
    Err(SolidError::invalid_input(
        "harness::merge_p3",
        "surface path outside this fixture adapter",
    ))
}
impl PiecewiseAccess for Pool {
    type Surface = (usize, usize);
    fn operands(&self) -> (Option<usize>, Option<usize>) {
        (self.primary, self.other)
    }
    fn set_operands(&mut self, a: Option<usize>, b: Option<usize>) {
        self.primary = a;
        self.other = b;
    }
    fn surfaces(&self, id: usize) -> solidcore::Result<Vec<Self::Surface>> {
        let brep = self
            .breps
            .get(id)
            .ok_or_else(|| SolidError::null_pointer("harness::merge_p3"))?;
        Ok(bp::surfaces(brep, None)?
            .iter()
            .enumerate()
            .map(|(s, _)| (id, s))
            .collect())
    }
    fn intersect_insert_relate(&mut self, _k: &mut Kernel) -> solidcore::Result<()> {
        // UNPORTED(boolean.topology_intersector_p2): SmTopologyIntersector::IntersectInsertRelate
        Err(SolidError::unported(
            "SmTopologyIntersector::IntersectInsertRelate",
            "boolean.topology_intersector_p2",
        ))
    }
    // No fixture reaches these methods: all piecewise cases defer at the real
    // IntersectInsertRelate dependency above. Unsupported adapter paths fail loudly.
    fn clear_relationships(&mut self) -> solidcore::Result<()> {
        absent()
    }
    fn set_subset(&mut self, _surfaces: &[Self::Surface]) {}
    fn surface_map(&mut self, _s: Self::Surface, _missing: bool) -> solidcore::Result<()> {
        absent()
    }
    fn surface_map_users(
        &self,
        _s: Self::Surface,
    ) -> solidcore::Result<Option<Vec<Self::Surface>>> {
        absent()
    }
    fn new_face_brep(&mut self, _k: &mut Kernel, _s: usize) -> solidcore::Result<usize> {
        absent()
    }
    fn manifold_flag(&self) -> bool {
        self.manifold
    }
    fn set_controls(&mut self, manifold: bool, _wires: bool, _lamina: u64) {
        self.manifold = manifold;
    }
    fn process_wires(&mut self, _k: &mut Kernel, _lock: &MarkLock) -> solidcore::Result<()> {
        absent()
    }
    fn release_marks(&mut self, _k: &mut Kernel, _lock: &mut MarkLock) -> solidcore::Result<()> {
        absent()
    }
    fn register_empty_ssi(
        &mut self,
        _a: Self::Surface,
        _b: Self::Surface,
    ) -> solidcore::Result<()> {
        absent()
    }
    fn clear_existing_ssi(&mut self) -> solidcore::Result<()> {
        absent()
    }
    fn notify_pre_edit(&mut self, _k: &mut Kernel, _a: usize) -> solidcore::Result<()> {
        absent()
    }
    fn delete_other(&mut self) -> solidcore::Result<()> {
        absent()
    }
}
fn rust(i: &Inputs) -> Result<Output, HarnessError> {
    if i.family > 4 {
        return Err(harness("fixture-only family is not a unit operation"));
    }
    ffi::validate(&encode(i)).map_err(|error| harness(&error))?;
    let mut pool = Pool::load(i)?;
    let mut k = Kernel::default();
    let opts = i.options.map(options);
    let mut values = vec![0.0];
    let result: solidcore::Result<()> = match i.family {
        0 => {
            let mut ids: Vec<_> = (0..i.a.len()).map(Some).collect();
            port::boolean_trees(&mut k, &mut ids, &i.tree, opts.as_ref(), &mut pool).and_then(
                |()| {
                    values.push(ids.len() as f64);
                    for id in ids {
                        if let Some(id) = id {
                            pool.record(id, &mut values)?;
                        } else {
                            values.extend([-1.0, 0.0]);
                        }
                    }
                    Ok(())
                },
            )
        }
        1 => {
            let mut ids: Vec<_> = (0..i.a.len()).collect();
            let op = u32::try_from(i.operation).unwrap_or(u32::MAX);
            port::merge_breps(&mut k, &mut ids, op, opts.as_ref(), &mut pool).and_then(|id| {
                values.push(1.0);
                pool.record(id, &mut values)
            })
        }
        2 => port::boolean_lists(
            &mut k,
            pool.breps
                .get(..i.a.len())
                .ok_or_else(|| harness("bad first list"))?,
            &[],
            pool.breps
                .get(i.a.len()..)
                .ok_or_else(|| harness("bad second list"))?,
            &[],
            i.operation,
            opts.as_ref(),
        )
        .map(|_| ()),
        3 => {
            let mut primary = Some(0);
            let mut other = Some(1);
            let mut manifold = true;
            let mut lamina = 0;
            let operation = match i.operation {
                0 => BooleanOperation::Unknown,
                1 => BooleanOperation::Union,
                2 => BooleanOperation::Intersection,
                3 => BooleanOperation::Difference,
                4 => BooleanOperation::ExclusiveOr,
                5 => BooleanOperation::Merge,
                6 => BooleanOperation::PartialMerge,
                7 => BooleanOperation::Imprint,
                8 => BooleanOperation::ImprintClassify,
                9 => BooleanOperation::ExtractSeparate,
                10 => BooleanOperation::Slice,
                _ => return Err(harness("invalid operation")),
            };
            port::non_manifold_boolean(
                &mut k,
                &mut pool,
                &mut primary,
                &mut other,
                &mut manifold,
                &mut lamina,
                operation,
            )
            .map(|_| ())
        }
        4 => port::piecewise_merge(
            &mut k,
            &mut pool,
            None,
            i.options.unwrap_or(0) & 1 != 0,
            i.options.unwrap_or(0) & 2 != 0,
        )
        .map(|_| ()),
        _ => return Err(harness("fixture-only family")),
    };
    match result {
        Ok(()) => output(&values, i.family),
        Err(error) => Ok(Output::unported(&error)
            .unwrap_or_else(|| Output::err(format!("status={}", status(&error))))),
    }
}
impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "boolean.merge_p3"
    }
    fn cases(&self, ctx: &CaseContext) -> Result<Vec<Case>, HarnessError> {
        crate::runner::default_cases(
            self,
            &CaseContext {
                seeds: ctx.seeds.max(256),
                ..ctx.clone()
            },
        )
    }
    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        rust(&serde_json::from_value(case.inputs.clone()).map_err(|error| harness(&error))?)
    }
    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs =
            serde_json::from_value(case.inputs.clone()).map_err(|error| harness(&error))?;
        output(
            &ffi::eval(&encode(&i)).map_err(|error| harness(&error))?,
            i.family,
        )
    }
    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        let family = match rng.below(10) {
            0..=5 => 0,
            6 => 1,
            7 => 2,
            8 => 3,
            _ => 4,
        };
        let n = if family >= 3 {
            1
        } else {
            rng.below(6) as usize
        };
        let m = if family == 0 || family == 1 {
            0
        } else if family >= 3 {
            1
        } else {
            rng.below(3) as usize
        };
        let scale = *rng.choose(&[1.0e-9, 1.0e-3, 1.0, 1.0e6, 1.0e12])?;
        let mut operand = || Operand {
            point: rng.point3(-scale, scale),
            size: None,
        };
        let a: Vec<_> = (0..n).map(|_| operand()).collect();
        let b: Vec<_> = (0..m).map(|_| operand()).collect();
        let mut tree: Vec<i64> = Vec::new();
        if family == 0 {
            let mut ids: Vec<_> = (0..n as i64).collect();
            for j in (1..ids.len()).rev() {
                ids.swap(j, rng.below((j + 1) as u64) as usize);
            }
            tree.extend(ids.into_iter().take(rng.below((n + 1) as u64) as usize));
            if rng.below(3) == 0 {
                tree.push(-rng.range_i64(1, 7));
            }
        }
        Some(json!(Inputs {
            family,
            operation: if family >= 3 {
                rng.below(11) as i32
            } else {
                rng.range_i64(-1, 5) as i32
            },
            options: if rng.below(2) == 0 {
                None
            } else {
                Some(rng.below(32) as u8)
            },
            a,
            b,
            tree
        }))
    }
}
