// SPDX-License-Identifier: Apache-2.0
#ifndef SC_OFFSET_GUARD_H
#define SC_OFFSET_GUARD_H
// Rejected oracle input/state, never an SMLib numerical result. Shared between the
// offset shim and its worker-only link wrappers; unwinds into the shim boundary.
struct ScOffsetUnsafeInput {};
#endif
