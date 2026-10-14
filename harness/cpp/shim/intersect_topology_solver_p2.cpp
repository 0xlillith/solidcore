// SPDX-License-Identifier: Apache-2.0
// Unit intersect.topology_solver_p2: real SmTopologySolver::LocalSolve, reached
// through SolveTrees / SolveTwoTrees. Reuse the p1 object-construction protocol
// and reference driver; its tree operation invokes the pinned C++ LocalSolve.
#include "sc_shim.h"

extern "C" int32_t sc_intersect_topology_solver_p1_eval(int32_t op,
    const double* input, int64_t n_in, double* out, int64_t cap, int64_t* out_len);

SC_EXPORT int32_t sc_intersect_topology_solver_p2_trees(const double* input,
    int64_t n_in, double* out, int64_t cap, int64_t* out_len)
{
    return sc_intersect_topology_solver_p1_eval(3, input, n_in, out, cap, out_len);
}

// Fixture preflight: unlike a managed node, an external point node is invisible
// to SmTree::GetObjectList, which walks the node manager rather than the root.
#include "StdAfx.h"
#include <SmTree.h>
#include <SmContext.h>
SC_EXPORT int32_t sc_intersect_topology_solver_p2_point_storage(int64_t* counts)
{
    SC_TRY
    SmContext context;
    SmTreeNode external_node;
    SmObjectList external_list, managed_list;
    SmTree external(SmPoint3d(1,2,3), 0.0, &external_node, &external_list);
    SmTree managed(SmPoint3d(1,2,3), 0.0, NULL, &managed_list);
    SmTArray<SmObjectList*> objects;
    SmStatus status = external.GetObjectList(objects);
    if (status != SM_SUCCESS) return SC_ERR_ARG;
    counts[0] = (int64_t)objects.GetSize();
    status = managed.GetObjectList(objects);
    if (status != SM_SUCCESS) return SC_ERR_ARG;
    counts[1] = (int64_t)objects.GetSize();
    return 0;
    SC_CATCH
}

SC_EXPORT int32_t sc_intersect_topology_solver_p2_local(const double* input,
    int64_t n_in, double* out, int64_t cap, int64_t* out_len)
{
    return sc_intersect_topology_solver_p1_eval(2, input, n_in, out, cap, out_len);
}
