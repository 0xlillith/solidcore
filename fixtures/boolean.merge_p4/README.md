The oracle calls the real `SmMerge::ManifoldBoolean`, with independently built
native B-reps. The four procedural operands are an empty B-rep, an isolated shell
vertex, a B-spline line wire, and an analytic plane sheet. Directed cases include reversed/crossing/touching/coincident wires, tolerance gaps,
coincident/overlapping sheets, tiny geometry, far-origin translation, and each
option flag with non-manifold mode. There are no file-reader
inputs or analytic-to-B-spline substitutions for plane sheets. The constructor's
automatically calculated tolerances are outside this unit: the shim supplies
1e-9 distance and 0.3 radians, initializes the manifold/loop flags, and configures
all five options explicitly. Each call owns its contexts and native topology.

Coverage limitations: the shared merge_p1 temporary-transform helper currently
returns `Unported(query.brep)` for every non-null operand pair. Those cases run the
whole C++ operation but are deferred on the Rust side, including C++ success paths.
Only null-operand errors can establish parity at present; they are both-errored,
zero-value cases and are not success-path evidence. The classification/copy/delete
and region-bookkeeping tail is present but has no end-to-end parity certification.
The runner's later adapter intentionally fails if called instead of supplying
invented classifications or relabeling missing adapter integration as a dependency.
Do not promote this unit to green from these cases.

The null-input cases exercise the same input-validation concern as
`TestKernelMergeBrepsInvalidInput`, directly at this unit's NER checks; they do not
claim coverage of the higher-level list API's singleton/duplicate-list contracts.
Material census, BooleanLists, PiecewiseMerge and real solid volume results remain
outside the tested surface until the dependencies and canonical adapter are wired.

`ManifoldIntersector` is an adapter interface, not a second implementation of its
callees. Required methods delegate existing merge/intersector ports; defaults
replace only calls into unported dependencies. `CloseIntersectionLoop` is listed as unused and absent from all unit assignments,
although `SmMerge.cpp:2958` calls `SmTopologyIntersector.cpp:6108-6354`. Per overseer
direction it is stubbed under its proposed owner `boolean.topology_intersector_p2`;
the overseer will co-port it in that unit in a follow-up. No implementation is
duplicated here, and the caller's retry/cleanup tail is retained.
