#ifndef STRATEGYSELECTOR_H_
#define STRATEGYSELECTOR_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2014 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

//Hardware data (abstraction of defines). Extend for other compilers

#if defined(_M_IX86) || defined(__i586__) || defined(__i686__) || defined(_M_X64) || defined(_M_AMD64) || defined(__amd64__) || defined(__x86_64__)
#define COMPILE_INTEL 1

#if defined(__MMX__)
#define COMPILE_INTEL_MMX 1
#endif

#if defined(__SSE__)
#define COMPILE_INTEL_SSE 1
#endif
    
#if defined(__SSE2__)
#define COMPILE_INTEL_SSE2 1
#endif

#if defined(__SSE3__)
#define COMPILE_INTEL_SSE3 1
#endif

#if defined(__SSSE3__)
#define COMPILE_INTEL_SSSE3 1
#endif

#if defined(__SSE4_1__)
#define COMPILE_INTEL_SSE41 1
#endif

#if defined(__SSE4_2__)
#define COMPILE_INTEL_SSE42 1
#endif

#if defined(__AVX__)
#define COMPILE_INTEL_AVX 1
#endif

#else
#define COMPILE_INTEL 0
#endif

#if defined (_M_PPC) || defined(__powerpc64__) || defined(__powerpc__)
#define COMPILE_POWERPC 1
#else
#define COMPILE_POWERPC 0
#endif

#if defined (_M_ARM) || defined(__arm__) || defined(__thumb__)
#define COMPILE_ARM 1
#else
#define COMPILE_ARM 0
#endif



typedef struct {
  const char *type; //Type of the function, usually its name
  const char *strategy_name; //Name of the strategy (e.g. sse2)
  unsigned int priority; //Priority. 0 = lowest (default strategy)
  void *fptr; //Pointer to the function
} strategy;

typedef struct {
  unsigned int count;
  unsigned int allocated;
  strategy* strategies;
} strategy_list;

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


int strategyselector_init();
void strategyselector_free();
int strategyselector_register(void *opaque, const char *type, const char *strategy_name, int priority, void *fptr);


//Strategy to include
#include "strategies/picture.h"

static const strategy_to_select strategies_to_select[] = {
  STRATEGIES_PICTURE_EXPORTS,
  {NULL, NULL},
};





#endif //STRATEGYSELECTOR_H_
