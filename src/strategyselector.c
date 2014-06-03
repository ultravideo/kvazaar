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

#include <assert.h>
#include <string.h>

#include "strategyselector.h"

hardware_flags g_hardware_flags;

static void set_hardware_flags();
static void* strategyselector_choose_for(const strategy_list * const strategies, const char * const strategy_type);

//Strategies to include (add new file here)
#include "strategies/picture.c"
#include "strategies/nal.c"

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
  
  if (!strategy_register_nal(&strategies)) {
    fprintf(stderr, "strategy_register_nal failed!\n");
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
  
  // Because VS doesn't support snprintf, let's just assert that there is more
  // than enough room.
  assert(strnlen(strategy_type, 200));
  sprintf(buffer, "KVAZAAR_OVERRIDE_%s", strategy_type);

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

#if defined(__GNUC__)
#include <cpuid.h>
INLINE int get_cpuid(unsigned int level, unsigned int *eax, unsigned int *ebx, unsigned int *ecx, unsigned int *edx) {
  return __get_cpuid(level, eax, ebx, ecx, edx);
}
#else
#include <intrin.h>
//Adapter from __cpuid (VS) to __get_cpuid (GNU C).
INLINE int get_cpuid(unsigned int level, unsigned int *eax, unsigned int *ebx, unsigned int *ecx, unsigned int *edx) {
  int CPUInfo[4] = {*eax, *ebx, *ecx, *edx};
  __cpuid(CPUInfo, 0);
  // check if the CPU supports the cpuid instruction.
  if (CPUInfo[0] != 0) {
    __cpuid(CPUInfo, level);
    *eax = CPUInfo[0];
    *ebx = CPUInfo[1];
    *ecx = CPUInfo[2];
    *edx = CPUInfo[3];
    return 1;
  }
  return 0;
}
#endif //defined(__GNUC__)

#endif

#if COMPILE_POWERPC
#include <unistd.h>
#include <fcntl.h>
#include <linux/auxvec.h>
#include <asm/cputable.h>

//Source: http://freevec.org/function/altivec_runtime_detection_linux
static int altivec_available(void)
{
    int result = 0;
    unsigned long buf[64];
    ssize_t count;
    int fd, i;
 
    fd = open("/proc/self/auxv", O_RDONLY);
    if (fd < 0) {
        return 0;
    }
    // loop on reading
    do {
        count = read(fd, buf, sizeof(buf));
        if (count < 0)
            break;
        for (i=0; i < (count / sizeof(unsigned long)); i += 2) {
            if (buf[i] == AT_HWCAP) {
                result = !!(buf[i+1] & PPC_FEATURE_HAS_ALTIVEC);
                goto out_close;
            } else if (buf[i] == AT_NULL)
                goto out_close;
        }
    } while (count == sizeof(buf));
out_close:
    close(fd);
    return result;
}
#endif //COMPILE_POWERPC

static void set_hardware_flags() {
  memset(&g_hardware_flags, 0, sizeof(g_hardware_flags));
  
  g_hardware_flags.arm = COMPILE_ARM;
  g_hardware_flags.intel = COMPILE_INTEL;
  g_hardware_flags.powerpc = COMPILE_POWERPC;
  
#if COMPILE_INTEL
  {
    unsigned int eax = 0, ebx = 0, ecx = 0, edx =0;
    /* CPU feature bits */
    enum { BIT_SSE3 = 0,BIT_SSSE3 = 9, BIT_SSE41 = 19, BIT_SSE42 = 20, BIT_MMX = 24, BIT_SSE = 25, BIT_SSE2 = 26, BIT_AVX = 28};

    // Dig CPU features with cpuid
    get_cpuid(1, &eax, &ebx, &ecx, &edx);
    
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

#if COMPILE_POWERPC
  g_hardware_flags.powerpc_flags.altivec = altivec_available();
  
  fprintf(stderr, "Compiled: PowerPC, flags:");
#if COMPILE_POWERPC_ALTIVEC
  fprintf(stderr, " AltiVec");
#endif
  fprintf(stderr, "\nRun on  : PowerPC, flags:");
  if (g_hardware_flags.powerpc_flags.altivec) fprintf(stderr, " AltiVec");
  fprintf(stderr, "\n");
#endif
  
}
