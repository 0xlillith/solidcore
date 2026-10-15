// SPDX-License-Identifier: Apache-2.0
// topology.shell: real native shell/region queries on shared B-rep fixtures.
// No native pointer crosses the C ABI; identities are traversal-order indices.
#include "sc_shim.h"
#include "support/topology_fixture.h"
#include "StdAfx.h"
#include <SmBrep.h>
#include <SmRegion.h>
#include <SmShell.h>
#include <SmFaceuse.h>
#include <SmTArray.h>
#include <algorithm>

namespace {
template<class T> double identity(const std::vector<T*>& values, T* p) {
    if (!p) return -1.0;
    auto found = std::find(values.begin(), values.end(), p);
    if (found == values.end()) throw 1;
    return (double)(found - values.begin());
}
template<class T> void list(std::vector<double>& out, const SmTArray<T*>& values,
                            const std::vector<T*>& identities) {
    out.push_back((double)values.GetSize());
    for (ULONG i = 0; i < values.GetSize(); ++i)
        out.push_back(identity(identities, values[i]));
}
}

SC_EXPORT int32_t sc_topology_shell_query(const char* path, const double* input,
    int64_t count, int32_t set_void, int32_t void_state,
    double* output, int64_t cap, int64_t* output_len) {
    if (!output_len || cap < 0 || (cap && !output) || (set_void != 0 && set_void != 1))
        return SC_ERR_ARG;
    SC_TRY
    return sc_topology_fixture::apply(path, input, count,
        [&](const sc_topology_fixture::View& fixture) -> int32_t {
            std::vector<double> out;
            out.push_back((double)fixture.shells.size());
            out.push_back((double)fixture.regions.size());
            for (SmShell* shell : fixture.shells) {
                SmTArray<SmFaceuse*> faceuses;
                shell->GetFaceuses(faceuses);
                // GetMatedShells dereferences both pointers unconditionally.
                for (ULONG i = 0; i < faceuses.GetSize(); ++i)
                    if (!faceuses[i]->GetMate() || !faceuses[i]->GetMate()->GetShell())
                        return SC_ERR_ARG;
                SmTArray<SmShell*> mates;
                if (shell->GetMatedShells(mates) != SM_SUCCESS) return SC_ERR_ARG;
                out.push_back((double)shell->GetType());
                out.push_back((double)shell->IsVertexShell());
                out.push_back((double)shell->IsWireShell());
                out.push_back((double)shell->IsFaceuseShell());
                out.push_back((double)shell->IsUnknownShell());
                out.push_back((double)shell->IsClosed());
                out.push_back(identity(fixture.regions, shell->GetRegion()));
                out.push_back(shell->GetBrep() == &fixture.brep ? 1.0 : 0.0);
                out.push_back(identity(fixture.faceuses, shell->GetFirstFaceuse()));
                out.push_back(identity(fixture.edges, shell->GetWireEdge()));
                out.push_back(identity(fixture.vertices, shell->GetVertex()));
                for (SM_TYPE type : {16013L, 16012L, 16000L, 16999L, -1L})
                    out.push_back((double)shell->IsKindOf(type));
                list(out, faceuses, fixture.faceuses);
                list(out, mates, fixture.shells);
            }
            for (SmRegion* region : fixture.regions) {
                if (set_void) region->SetIsVoid(void_state);
                out.push_back((double)region->GetType());
                out.push_back((double)region->IsVoid());
                out.push_back((double)region->IsInfiniteRegion());
                out.push_back(region->GetBrep() == &fixture.brep ? 1.0 : 0.0);
                for (SM_TYPE type : {16013L, 16012L, 16000L, 16999L, -1L})
                    out.push_back((double)region->IsKindOf(type));
                SmTArray<SmShell*> shells;
                region->GetShells(shells);
                list(out, shells, fixture.shells);
            }
            return sc_copy_out(out.data(), (int64_t)out.size(), output, cap, output_len);
        });
    SC_CATCH
}

SC_EXPORT int32_t sc_topology_shell_wire_edges(const char* path, const double* input,
    int64_t count, double* output, int64_t cap, int64_t* output_len) {
    if (!output_len || cap < 0 || (cap && !output)) return SC_ERR_ARG;
    SC_TRY
    return sc_topology_fixture::apply(path, input, count,
        [&](const sc_topology_fixture::View& fixture) -> int32_t {
            std::vector<double> out;
            out.push_back((double)fixture.shells.size());
            for (SmShell* shell : fixture.shells) {
                if (shell->IsWireShell() && !shell->GetWireEdge()) return SC_ERR_ARG;
                SmTArray<SmEdge*> wires;
                const SmStatus status = shell->GetWireEdges(wires);
                out.push_back(status == SM_SUCCESS ? 0.0 : (double)status);
                list(out, wires, fixture.edges);
            }
            return sc_copy_out(out.data(), (int64_t)out.size(), output, cap, output_len);
        });
    SC_CATCH
}
