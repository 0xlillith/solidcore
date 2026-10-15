// SPDX-License-Identifier: Apache-2.0
// Disposable oracle process for surface.offset_surface. Binary native-endian protocol
// is private to the same-host Rust wrapper; never load an input-provided file/path.
#include "shim/sc_shim.h"
#include "shim/sc_offset_guard.h"
#include "StdAfx.h"
#include <nurbs.h>
#include <NL_Globals.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <limits>

extern "C" int32_t sc_surface_offset_surface_eval(int32_t, const double*, int64_t,
                                                   double*, int64_t, int64_t*);

// GNU link wrapping is confined to this disposable worker. All finite results are
// returned verbatim by the original SMLib functions. Reject only an extension
// search whose initial length is non-finite or whose first parameter step cannot
// advance, rather than changing tolerances or banning a geometry/continuity kind.
extern "C" NL_FLAG __real_N_CrvExtendByDist(NL_CURVE*, NL_REAL, NL_FLAG, NL_FLAG,
                                            NL_FLAG, NL_CURVE*, NL_STACKS*, NL_STACKS*);
extern "C" NL_FLAG __real_N_CrvArcLength(NL_CURVE*, NL_PARAMETER, NL_PARAMETER,
                                        NL_REAL, NL_FLAG, NL_REAL*);
namespace {
struct Extension {
    NL_CURVE* curve;
    NL_REAL distance;
    NL_FLAG relative, end, continuity;
};
// The worker is single-threaded. The scoped pointer also handles nested extensions.
const Extension* active_extension = nullptr;
struct ExtensionScope {
    const Extension* previous;
    explicit ExtensionScope(const Extension& current) : previous(active_extension) {
        active_extension = &current;
    }
    ~ExtensionScope() { active_extension = previous; }
};
}
extern "C" NL_FLAG __wrap_N_CrvExtendByDist(NL_CURVE* curve, NL_REAL distance,
    NL_FLAG relative, NL_FLAG end, NL_FLAG continuity, NL_CURVE* result,
    NL_STACKS* input_stack, NL_STACKS* output_stack)
{
    const Extension extension{curve, distance, relative, end, continuity};
    ExtensionScope scope(extension);
    return __real_N_CrvExtendByDist(curve, distance, relative, end, continuity,
                                   result, input_stack, output_stack);
}
extern "C" NL_FLAG __wrap_N_CrvArcLength(NL_CURVE* curve, NL_PARAMETER lo,
    NL_PARAMETER hi, NL_REAL tolerance, NL_FLAG flag, NL_REAL* length)
{
    const auto status = __real_N_CrvArcLength(curve, lo, hi, tolerance, flag, length);
    const Extension* ext = active_extension;
    if (!ext || status != NL_NO) return status;
    if (!std::isfinite(*length)) throw ScOffsetUnsafeInput{};
    // ST_curxcdr's initial arc-length call and precisely its subsequent arithmetic
    // (NL_CrvAdv.cpp:683-742). The small-curve and recursive-large-distance branches
    // return/recurse before this search, so do not reject those legitimate paths.
    if (ext->curve == curve && ext->continuity == NL_G2 && tolerance == 0.01 &&
        flag == NL_RELATIVE && !(*length <= NL_MTOL && *length <= NL_PTOL))
    {
        NL_REAL distance = ext->distance;
        if (ext->relative == NL_RELATIVE) distance *= *length;
        if (!(distance > 0.9313 * *length))
        {
            NL_REAL step = (hi - lo) * (distance / *length);
            NL_REAL endpoint = ext->end == NL_START ? lo : hi;
            NL_REAL next = ext->end == NL_START ? endpoint + step : endpoint - step;
            if (!std::isfinite(step) || next == endpoint) throw ScOffsetUnsafeInput{};
        }
    }
    return status;
}

template <class T> bool read(T& value)
{
    return bool(std::cin.read(reinterpret_cast<char*>(&value), sizeof(value)));
}
template <class T> void write(const T& value)
{
    std::cout.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

int main()
{
    int32_t op;
    uint64_t count;
    if (!read(op) || !read(count) || count == 0 || count > 1048576) return 2;
    std::vector<double> input(count);
    if (!std::cin.read(reinterpret_cast<char*>(input.data()), count * sizeof(double))) return 2;
    std::vector<double> output(4096);
    int64_t length = 0;
    int32_t status = sc_surface_offset_surface_eval(op, input.data(), count,
                                                   output.data(), output.size(), &length);
    if (status == SC_ERR_BUFFER && length > 0 && length <= 8388608)
    {
        output.resize(length);
        status = sc_surface_offset_surface_eval(op, input.data(), count,
                                               output.data(), output.size(), &length);
    }
    if (status != SC_OK) length = 0;
    if (length < 0 || uint64_t(length) > output.size()) return 2;
    write(status);
    write(uint64_t(length));
    std::cout.write(reinterpret_cast<const char*>(output.data()), length * sizeof(double));
    return std::cout ? 0 : 2;
}
