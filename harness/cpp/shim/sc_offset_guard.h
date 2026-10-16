// SPDX-License-Identifier: Apache-2.0
#ifndef SC_OFFSET_GUARD_H
#define SC_OFFSET_GUARD_H
// Rejected oracle input/state, never an SMLib numerical result. Shared between the
// offset shim and its worker-only link wrappers; unwinds into the shim boundary.
struct ScOffsetUnsafeInput {};

// Nonzero when the offset worker is linked with its NLib wrappers (GNU ld --wrap): the
// guards for UB-107 and the HANG-6/HANG-15 extension searches. The shim's weak default
// returns 0; the worker's strong definition overrides it where the wrappers exist.
extern "C" int sc_offset_extension_guards();
inline bool sc_offset_extensions_guarded() { return sc_offset_extension_guards() != 0; }
#endif
