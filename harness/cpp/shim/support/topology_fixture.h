// SPDX-License-Identifier: Apache-2.0
// Internal C++ fixture sharing. No native object crosses the Rust FFI boundary.
#pragma once
#include <cstdint>
#include <functional>
#include <vector>

class SmBrep;
class SmRegion;
class SmShell;
class SmFaceuse;
class SmFace;
class SmLoopuse;
class SmLoop;
class SmEdgeuse;
class SmEdge;
class SmVertexuse;
class SmVertex;

namespace sc_topology_fixture {
// Traversal-order identities of an independently constructed native B-rep.
// All references are valid only during apply's callback.
struct View {
    SmBrep& brep;
    const std::vector<SmRegion*>& regions;
    const std::vector<SmShell*>& shells;
    const std::vector<SmFaceuse*>& faceuses;
    const std::vector<SmFace*>& faces;
    const std::vector<SmLoopuse*>& loopuses;
    const std::vector<SmLoop*>& loops;
    const std::vector<SmEdgeuse*>& edgeuses;
    const std::vector<SmEdge*>& edges;
    const std::vector<SmVertexuse*>& vertexuses;
    const std::vector<SmVertex*>& vertices;
};

// Uses the existing edge-use fixture protocol and native constructors. Setup
// errors fail the harness; the B-rep is destroyed before this function returns.
int32_t apply(const char* path, const double* input, int64_t count,
              const std::function<int32_t(const View&)>& operation);
}
