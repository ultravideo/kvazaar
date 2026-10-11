#ifndef STRATEGIES_CPP_H_
#define STRATEGIES_CPP_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 *
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

/**
 * \ingroup Optimization
 * \file
 * Prologue for strategy implementations written in C++17.
 *
 * The encoder is C. Some SIMD strategies are C++ so that block sizes, filter
 * taps and output modes can be template parameters, which makes their
 * specialization a compile-time guarantee instead of depending on inlining
 * heuristics, and lets intrinsics that need immediate operands take values
 * derived from them. Rules for these files:
 *
 * - Include this header first. It includes the C library and intrinsic
 *   headers, which have C++-specific declarations, before any encoder header
 *   can include them inside an extern "C" block.
 * - Include encoder headers inside extern "C" { ... } so their declarations,
 *   including the strategy register function, keep C linkage.
 * - Put everything except the register function in an anonymous namespace,
 *   so that no C++ symbols (template instantiations) are exported from the
 *   library.
 * - Do not use exceptions, RTTI, new/delete or anything else that needs the
 *   C++ runtime library; the library is linked as C. The build compiles these
 *   files with -std=c++17 -fno-exceptions -fno-rtti.
 */

#ifndef __cplusplus
#error "strategies-cpp.h is only for strategy implementations compiled as C++"
#endif

#include <cassert>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__aarch64__) || defined(_M_ARM64)
#include <arm_neon.h>
#endif
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#include <immintrin.h>
#endif

// Forces inlining where the specialization relies on it (for example so that
// arrays of vector accumulators stay in registers).
#if defined(_MSC_VER)
#define KVZ_FORCE_INLINE __forceinline
#else
#define KVZ_FORCE_INLINE inline __attribute__((always_inline))
#endif

// Fully unrolls the following loop over a compile-time trip count, for
// loops whose vector state must stay in registers.
#if defined(__clang__)
#define KVZ_UNROLL _Pragma("unroll")
#elif defined(__GNUC__) && __GNUC__ >= 8
#define KVZ_UNROLL _Pragma("GCC unroll 64")
#else
#define KVZ_UNROLL
#endif

// The strategy selector takes function pointers as void *; C++ requires the
// conversion to be explicit.
#define KVZ_STRATEGY_FPTR(f) (reinterpret_cast<void *>(&(f)))

#endif // STRATEGIES_CPP_H_
