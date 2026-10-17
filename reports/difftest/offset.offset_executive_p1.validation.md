# offset.offset_executive_p1 — blocked validation, PR #175

This is an **isolated-case investigation**, not a normal aggregate difftest
report, and does not certify this unit as green or partial.

## Preserved implementation

- Repository: [altvnk/solidcore](https://github.com/altvnk/solidcore).
- Branch: `wip/offset.offset_executive_p1-pr175`.
- Tested commit: [`e173a72239c661892c167523528eb5301bfaee29`](https://github.com/altvnk/solidcore/commit/e173a72239c661892c167523528eb5301bfaee29).
- Kernel: `crates/solidcore/src/offset/offset_executive/p1.rs`.
- Adapter: `harness/cpp/shim/offset_offset_executive_p1.cpp` and
  `crates/solidcore-ffi/src/units/offset_offset_executive_p1.rs`.
- Runner: `crates/difftest/src/units/offset_offset_executive_p1.rs`.
- Inputs: `fixtures/offset.offset_executive_p1/*.json`, plus 600 fixed-seed
  generated cases (`gen0000` through `gen0599`).

The findings-only branch deliberately does not register this failing runner or
ship the uncertified implementation. Recover the WIP branch to reproduce it.
The WIP includes settings, custom geometry-creator scalar snapshots,
correspondence, synchronization, surface-extension caching, and reachable
construction prefixes. Unavailable callees return owner-tagged `Unported`;
continuations behind those calls are documented, not claimed as fully ported.

## Environment and method

Validated on 2026-10-17 using Linux x86_64, GNU g++ 16.2.1 and Rust 1.98.0.
The reference source is pinned to `c9c5979f04cea8c38213004d4f04df0dcfb40d05`
and built by the repository harness (`-ffp-contract=off`, no fast math).

Every original case ran in its own native process with an eight-second limit.
A temporary selector chose one exact case index and delegated both engines,
absolute floors and comparison to the existing runner; no input, output order,
generator or tolerance was changed. Progress markers distinguished a timeout
inside C++ from a Rust or comparison timeout. Only owned subprocesses were
terminated and waited for. No absolute floors are defined by this runner.

| Inputs | Passed | Failed (order) | Deferred | Rejected | Native timeout | Total |
|---|---:|---:|---:|---:|---:|---:|
| Fixtures | 73 | 2 | 8 | 4 | 4 | 91 |
| Generated | 450 | 0 | 9 | 37 | 104 | 600 |
| **Total** | **523** | **2** | **17** | **41** | **108** | **691** |

No completed comparison had a floating-point mismatch. Deferred, rejected and
timed-out cases do **not** count as parity passes. No crashes or harness errors
were observed in this scan. An ordinary whole-unit invocation cannot complete
within its budget because it reaches the negative-distance native cases.

### Observable order failures (ND-2)

Both `equal_minima_public_order_tie` and `equal_minima_matched_order_tie` in
`synchronize.json` fail with:

```text
ints.matched_first[0]: cpp=2 rust=0
ints.matched_first[1]: cpp=0 rust=2
```

The first fixture maps two analytic lines, from the origin to `(1,0,0)` and
`(2,0,0)`, each with natural and trimmed interval `[0,1]`, to one original
edgeuse. There are no other-side mappings. The second uses real quadratic
B-splines and also exercises matching against the radial-side range.

`SmRelation::GetSeconds` enumerates pointer-hashed unordered membership.
Synchronization's stable bubble sort swaps only strictly greater minima, so
its public tie order retains that unspecified order. We did not canonicalize
public results or assert that address-dependent variation was experimentally
measured within this binary. A portable reference-order policy is needed.

### Negative extension timeouts (HANG-16)

These fixtures each exceeded eight seconds inside C++:

- `extend_plane_negative`
- `extend_bilinear_negative`
- `extend_quadratic_negative`
- `extend_rational_quadratic_negative`

Their offset distance is `-0.125`, so the executive requests extension by
`-0.625`. All 104 generated timeouts also use negative distances, including
`-1e-12` and `-1e7`. Example: `gen0004`, seed `13762035380551339612`,
distance `-1e7`. The extension loop at `SmBSplineSurfaceNLib.cpp:4583-4591`
doubles a negative value while testing it against half the negative request,
which moves away from termination. The experiment establishes native timeouts,
not infinite execution or an experimentally identified stuck instruction.

## Reproduce from the WIP branch

After obtaining the preserved WIP and building the pinned C++ reference:

```sh
CARGO_BUILD_JOBS=2 nice -n 10 cargo build --locked --profile difftest -p difftest --features cpp

# Both observable-order regressions: exits nonzero, not a tolerance issue.
timeout --kill-after=2s 8 nice -n 10 bin/difftest/difftest \
  --unit offset.offset_executive_p1 --filter equal_minima --no-report

# Native negative-distance timeout. Exit 124 is an external timeout, NOT parity.
timeout --kill-after=2s 8 nice -n 10 bin/difftest/difftest \
  --unit offset.offset_executive_p1 --filter extend_plane_negative --no-report

# One fixed-seed negative-distance case, without an aggregate run.
timeout --kill-after=2s 8 nice -n 10 bin/difftest/difftest \
  --unit offset.offset_executive_p1 --filter gen0004 --no-report
```

The CLI uses substring filtering; the plane-negative filter also includes its
negative-zero neighbor, but the negative-distance fixture reaches the timeout
first. Compile before starting the eight-second execution limit.

## Other checks on the preserved WIP

| Check | Result |
|---|---|
| `cargo fmt --all --check` | PASS |
| `cargo clippy --locked --workspace --all-targets --features cpp -- -D warnings` | PASS |
| `cargo test --locked -p solidcore` | PASS: 1288 tests, 0 failed; doc-tests pass |
| `cargo check --locked -p solidcore --target wasm32-unknown-unknown` | PASS |
| Native `harness.selftest` | PASS: 28/28 |
| Optional `--green --no-report` | Incomplete at 300 s; 101 PASS / 32 PARTIAL summaries, no reported failures before timeout |

No full green-baseline pass is claimed. All builds ran under `nice -n 10`, with
two Cargo jobs and six reference-build jobs. Colima was not restarted or resized,
and no additional VM was created.

Raw logs, exact inputs, source hashes, all 691 outcomes and the 108-case timeout
list remain in the contributor checkout under
`bin/remote-pr175/latest-validation/`, and on the isolated Linux builder under
`/home/altvnk/solidcore-pr175-zed/bin/zed-case-validation/`. These ignored local
artifacts are not a machine-readable certification report. The immutable WIP
commit and named/seeded inputs above are the reproducible handoff.
