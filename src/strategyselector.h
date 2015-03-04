#ifndef STRATEGYSELECTOR_H_
#define STRATEGYSELECTOR_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include "global.h"

#if defined(_DEBUG) && !defined(DEBUG_STRATEGYSELECTOR)
# define DEBUG_STRATEGYSELECTOR
#endif

//Hardware data (abstraction of defines). Extend for other compilers

#if defined(_M_IX86) || defined(__i586__) || defined(__i686__) || defined(_M_X64) || defined(_M_AMD64) || defined(__amd64__) || defined(__x86_64__)
#  define COMPILE_INTEL 1
#else
#  define COMPILE_INTEL 0
#endif

// Visual Studio note:
// Because these macros are only used to guard code that is guarded by CPUID
// at runtime, use /arch parameter to disable them, but enable all intrinsics
// supported by VisualStudio if SSE2 (highest) is enabled.
// AVX and AVX2 are handled by /arch directly and sse intrinsics will use VEX
// versions if they are defined.
#define MSC_X86_SIMD(level) (_M_X64 || (_M_IX86_FP >= (level)))

#if COMPILE_INTEL
#  if defined(__MMX__) || MSC_X86_SIMD(1)
#    define COMPILE_INTEL_MMX 1
#  endif
#  if defined(__SSE__) || MSC_X86_SIMD(1)
#    define COMPILE_INTEL_SSE 1
#  endif
#  if defined(__SSE2__) || MSC_X86_SIMD(2)
#    define COMPILE_INTEL_SSE2 1
#  endif
#  if defined(__SSE3__)
#    define COMPILE_INTEL_SSE3 1
#  endif
#  if defined(__SSSE3__) || MSC_X86_SIMD(2)
#    define COMPILE_INTEL_SSSE3 1
#  endif
#  if defined(__SSE4_1__) || MSC_X86_SIMD(2)
#    define COMPILE_INTEL_SSE41 1
#  endif
#  if defined(__SSE4_2__) || MSC_X86_SIMD(2)
#    define COMPILE_INTEL_SSE42 1
#  endif
#  if defined(__AVX__)
#    define COMPILE_INTEL_AVX 1
#   endif
#  if defined(__AVX2__)
#    define COMPILE_INTEL_AVX2 1
#   endif
#endif

#if defined (_M_PPC) || defined(__powerpc64__) || defined(__powerpc__)
#  define COMPILE_POWERPC 1
#  ifdef __ALTIVEC__
#    define COMPILE_POWERPC_ALTIVEC 1
#  else
#    define COMPILE_POWERPC_ALTIVEC 0
#  endif
#else
#  define COMPILE_POWERPC 0
#endif

#if defined (_M_ARM) || defined(__arm__) || defined(__thumb__)
#  define COMPILE_ARM 1
#else
#  define COMPILE_ARM 0
#endif



typedef struct {
  const char *type; //Type of the function, usually its name
  const char *strategy_name; //Name of the strategy (e.g. sse2)
  unsigned int priority; //Priority. 0 = lowest (default strategy)
  void *fptr; //Pointer to the function
} strategy_t;

typedef struct {
  unsigned int count;
  unsigned int allocated;
  strategy_t* strategies;
} strategy_list_t;

#define STRATEGY_LIST_ALLOC_SIZE 16

typedef struct {
  const char *strategy_type;
  void **fptr;
} strategy_to_select;

typedef struct {
  int intel;
  struct {
    int mmx;
    int sse;
    int sse2;
    int sse3;
    int ssse3;
    int sse41;
    int sse42;
    int avx;
    int avx2;
  } intel_flags;
  
  int powerpc;
  struct {
    int altivec;
  } powerpc_flags;
  
  int arm;
  struct {
    int neon;
  } arm_flags;
} hardware_flags;

extern hardware_flags g_hardware_flags;


int strategyselector_init(int32_t cpuid);
void strategyselector_free();
int strategyselector_register(void *opaque, const char *type, const char *strategy_name, int priority, void *fptr);


//Strategy to include
#include "strategies/strategies-nal.h"
#include "strategies/strategies-picture.h"
#include "strategies/strategies-dct.h"
#include "strategies/strategies-ipol.h"

static const strategy_to_select strategies_to_select[] = {
  STRATEGIES_NAL_EXPORTS
  STRATEGIES_PICTURE_EXPORTS
  STRATEGIES_DCT_EXPORTS
  STRATEGIES_IPOL_EXPORTS
  { NULL, NULL },
};

unsigned satd_8bit_8x8_generic(const pixel * const block1, const pixel * const block2);


#endif //STRATEGYSELECTOR_H_
