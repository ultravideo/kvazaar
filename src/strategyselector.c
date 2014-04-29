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

/*
 * \file
 */

#include <string.h>

#include "strategyselector.h"

hardware_flags g_hardware_flags;

static void set_hardware_flags();
static void* strategyselector_choose_for(const strategy_list * const strategies, const char * const strategy_type);

//Strategies to include (add new file here)
#include "strategies/picture.c"

//Returns 1 if successful
int strategyselector_init() {
  const strategy_to_select *cur_strategy_to_select = strategies_to_select;
  strategy_list strategies;
  
  strategies.allocated = 0;
  strategies.count = 0;
  strategies.strategies = NULL;
  
  set_hardware_flags();
  
  //Add new register function here
  if (!strategy_register_picture(&strategies)) {
    fprintf(stderr, "strategy_register_picture failed!\n");
    return 0;
  }
  
  while(cur_strategy_to_select->fptr) {
    *(cur_strategy_to_select->fptr) = strategyselector_choose_for(&strategies, cur_strategy_to_select->strategy_type);
    
    if (!(*(cur_strategy_to_select->fptr))) {
      fprintf(stderr, "Could not find a strategy for %s!\n", cur_strategy_to_select->strategy_type);
      return 0;
    }
    ++cur_strategy_to_select;
  }
  
  //We can free the structure now, as all strategies are statically set to pointers
  if (strategies.allocated) {
    free(strategies.strategies);
  }

  return 1;
}

void strategyselector_free() {
  //Do nothing (yet)
}

//Returns 1 if successful, 0 otherwise
int strategyselector_register(void * const opaque, const char * const type, const char * const strategy_name, int priority, void * const fptr) {
  strategy_list * const strategies = opaque;
  
  if (strategies->allocated == strategies->count) {
    strategy* new_strategies = realloc(strategies->strategies, sizeof(strategy) * (strategies->allocated + STRATEGY_LIST_ALLOC_SIZE));
    if (!new_strategies) {
      fprintf(stderr, "Could not increase strategies list size!\n");
      return 0;
    }
    strategies->strategies = new_strategies;
    strategies->allocated += STRATEGY_LIST_ALLOC_SIZE;
  }
  
  {
    strategy *new_strategy = &strategies->strategies[strategies->count++];
    new_strategy->type = type;
    new_strategy->strategy_name = strategy_name;
    new_strategy->priority = priority;
    new_strategy->fptr = fptr;
  }
#ifdef _DEBUG
  fprintf(stderr, "Registered strategy %s:%s with priority %d (%p)\n", type, strategy_name, priority, fptr);
#endif //_DEBUG
  
  return 1;
}

static void* strategyselector_choose_for(const strategy_list * const strategies, const char * const strategy_type) {
  unsigned int max_priority = 0;
  int max_priority_i = -1;
  char buffer[256];
  char *override = NULL;
  int i = 0;
  
  snprintf(buffer, 255, "KVAZAAR_OVERRIDE_%s", strategy_type);
  override = getenv(buffer);
  
  for (i=0; i < strategies->count; ++i) {
    if (strcmp(strategies->strategies[i].type, strategy_type) == 0) {
      if (override && strcmp(strategies->strategies[i].strategy_name, override) == 0) {
        fprintf(stderr, "%s environment variable present, choosing %s:%s\n", buffer, strategy_type, strategies->strategies[i].strategy_name);
        return strategies->strategies[i].fptr;
      }
      if (strategies->strategies[i].priority >= max_priority) {
        max_priority_i = i;
        max_priority = strategies->strategies[i].priority;
      }
    }
  }
  
  if (override) {
    fprintf(stderr, "%s environment variable present, but no strategy %s was found!\n", buffer, override);
    return NULL;
  }

#ifdef _DEBUG
  fprintf(stderr, "Choosing strategy for %s:\n", strategy_type);
  for (i=0; i < strategies->count; ++i) {
    if (strcmp(strategies->strategies[i].type, strategy_type) == 0) {
      if (i != max_priority_i) {
        fprintf(stderr, "- %s (%d, %p)\n", strategies->strategies[i].strategy_name, strategies->strategies[i].priority, strategies->strategies[i].fptr);
      } else {
        fprintf(stderr, "> %s (%d, %p)\n", strategies->strategies[i].strategy_name, strategies->strategies[i].priority, strategies->strategies[i].fptr);
      }
    }
  }
#endif //_DEBUG
  
  
  if (max_priority_i == -1) {
    return NULL;
  }
  
  return strategies->strategies[max_priority_i].fptr;
}

#if COMPILE_INTEL
#include "x86/cpu.h"
#endif

static void set_hardware_flags() {
  memset(&g_hardware_flags, 0, sizeof(g_hardware_flags));
  
  g_hardware_flags.arm = COMPILE_ARM;
  g_hardware_flags.intel = COMPILE_INTEL;
  g_hardware_flags.powerpc = COMPILE_POWERPC;
  
#if COMPILE_INTEL
  {
    int ecx = 0,edx =0;
    /* CPU feature bits */
    enum { BIT_SSE3 = 0,BIT_SSSE3 = 9, BIT_SSE41 = 19, BIT_SSE42 = 20, BIT_MMX = 24, BIT_SSE = 25, BIT_SSE2 = 26, BIT_AVX = 28};

    // Dig CPU features with cpuid
    kvz_cpu_cpuid(&ecx,&edx);
    
    // EDX
    if (edx & (1<<BIT_MMX))   g_hardware_flags.intel_flags.mmx = 1;
    if (edx & (1<<BIT_SSE))   g_hardware_flags.intel_flags.sse = 1;
    if (edx & (1<<BIT_SSE2))  g_hardware_flags.intel_flags.sse2 = 1;
    // ECX
    if (ecx & (1<<BIT_SSE3))  g_hardware_flags.intel_flags.sse3 = 1;;
    if (ecx & (1<<BIT_SSSE3)) g_hardware_flags.intel_flags.ssse3 = 1;
    if (ecx & (1<<BIT_SSE41)) g_hardware_flags.intel_flags.sse41 = 1;
    if (ecx & (1<<BIT_SSE42)) g_hardware_flags.intel_flags.sse42 = 1;
    if (ecx & (1<<BIT_AVX))   g_hardware_flags.intel_flags.avx = 1;
    
    fprintf(stderr, "Compiled: INTEL, flags:");
#if COMPILE_INTEL_MMX
    fprintf(stderr, " MMX");
#endif
#if COMPILE_INTEL_SSE
    fprintf(stderr, " SSE");
#endif
#if COMPILE_INTEL_SSE2
    fprintf(stderr, " SSE2");
#endif
#if COMPILE_INTEL_SSE3
    fprintf(stderr, " SSE3");
#endif
#if COMPILE_INTEL_SSSE3
    fprintf(stderr, " SSSE3");
#endif
#if COMPILE_INTEL_SSE41
    fprintf(stderr, " SSE41");
#endif
#if COMPILE_INTEL_SSE42
    fprintf(stderr, " SSE42");
#endif
#if COMPILE_INTEL_AVX
    fprintf(stderr, " AVX");
#endif
    fprintf(stderr, "\nRun on  : INTEL, flags:");
    if (g_hardware_flags.intel_flags.mmx) fprintf(stderr, " MMX");
    if (g_hardware_flags.intel_flags.sse) fprintf(stderr, " SSE");
    if (g_hardware_flags.intel_flags.sse2) fprintf(stderr, " SSE2");
    if (g_hardware_flags.intel_flags.sse3) fprintf(stderr, " SSE3");
    if (g_hardware_flags.intel_flags.ssse3) fprintf(stderr, " SSSE3");
    if (g_hardware_flags.intel_flags.sse41) fprintf(stderr, " SSE41");
    if (g_hardware_flags.intel_flags.sse42) fprintf(stderr, " SSE42");
    if (g_hardware_flags.intel_flags.avx) fprintf(stderr, " AVX");
    fprintf(stderr, "\n");
  }
#endif //COMPILE_INTEL
}
