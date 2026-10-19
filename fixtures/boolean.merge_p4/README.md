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
Only null-operand errors can establish parity at present. Their status is data,
so both sides compare the post-operation doing-Boolean flags, editing flags,
surviving operand graphs and intersection/component counts even on errors.
The native shim exports each full graph; the runner imports only the native
post-state on the C++ side and serializes it through the shared brep_p10 traversal
serializer (exact integer topology, geometry samples under the normal 1e-9 rule,
plus lazy gap/UV state). The Rust side serializes its independently mutated arenas.
Error-state matches remain **zero Boolean success paths**, even though they now
compare values. All eleven operations lack success-path parity. The complete
classification/copy/delete/region-bookkeeping tail has source-review evidence only.
Do not promote this unit to green from these cases.

The null-input cases exercise the same input-validation concern as
`TestKernelMergeBrepsInvalidInput`, directly at this unit's NER checks; they do not
claim coverage of the higher-level list API's singleton/duplicate-list contracts.
Material census, BooleanLists, PiecewiseMerge and real solid volume results remain
outside the tested surface until the dependencies and canonical adapter are wired.

No concrete kernel implementation of `ManifoldIntersector` exists. Its required
methods do not yet delegate to the already-ported callees:

- merge_p1: fix_found_face_coincidence, find_brep_region_of_other_face,
  propagate_other_src_region_attribs, boolean_post_process;
- merge_p5: process_wires_and_shell_vertices;
- topology_intersector_p1: Relations mates/relate and TopologyIntersector common
  and mapped object lists, surface subset and topology-deleted state.

Every required method has an `INTEGRATION-GAP(boolean.merge_p4)` marker. Owner:
the boolean.merge_p4 author follow-up, including reuse of p1::MergeState,
canonical relations, notification records and the caller AttributeStore.
The harness DeferredAdapter errors/panics if reached; wiring query.brep must not
turn an adapter gap into a pass or an Unported label for an already-ported unit.
Only genuinely missing dependency calls use the trait's Unported defaults.

`CloseIntersectionLoop` is listed as unused and absent from all unit assignments,
although `SmMerge.cpp:2958` calls `SmTopologyIntersector.cpp:6108-6354`. Per overseer
direction it is stubbed under its proposed owner `boolean.topology_intersector_p2`;
the overseer will co-port it in that unit in a follow-up. No implementation is
duplicated here, and the caller's retry/cleanup tail is retained.

Scope and ownership: the temporary transform now returns an inverse token whose
scope covers the complete operation and restores only surviving operands on
error and imprint exits. The wire/face-copy adapter calls borrow the caller's
same AttributeStore; canonical integration remains pending. Native snapshots
include every extracted component graph (record length, editing flag, graph),
in primary-then-other component order. A plane-domain constructor rejection
(width 2000000) exercises cleanup of geometry not yet owned by a face.

The #203 coordination explicitly requires rechecking CloseIntersectionLoop when
that unit lands/goes green; the proposed-owner tag must not outlive its co-port.
The FFI supports coordinates/dimensions up to 1e100 via shared reader bounds.
Its operand families avoid the closed-edge UB-175 path; adding closed curves
requires a new guard review.
