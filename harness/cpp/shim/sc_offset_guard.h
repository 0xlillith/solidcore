// SPDX-License-Identifier: Apache-2.0
#ifndef SC_OFFSET_GUARD_H
#define SC_OFFSET_GUARD_H
// Rejected oracle input/state, never an SMLib numerical result. Shared between the
// offset shim and its worker-only link wrappers; unwinds into the shim boundary.
struct ScOffsetUnsafeInput {};

// Set to 1 by the offset worker's main() when it is linked with its NLib wrappers (GNU
// ld --wrap): the guards for UB-107 and the HANG-6/HANG-15 extension searches. Defined
// as 0 in the shim (a plain variable, no weak symbols: the shim builds everywhere).
extern "C" int sc_offset_extension_guards_enabled;
inline bool sc_offset_extensions_guarded() { return sc_offset_extension_guards_enabled != 0; }
#endif
