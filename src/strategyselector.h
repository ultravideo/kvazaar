#ifndef STRATEGYSELECTOR_H_
#define STRATEGYSELECTOR_H_
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
 * Dynamic dispatch based on cpuid.
 */

#include "global.h" // IWYU pragma: keep

#if defined(KVZ_DEBUG) && !defined(DEBUG_STRATEGYSELECTOR)
# define DEBUG_STRATEGYSELECTOR
#endif

typedef struct {
  const char *type; //Type of the function, usually its name
  const char *strategy_name; //Name of the strategy (e.g. sse2)
  unsigned int priority; //Priority. 0 = lowest (default strategy)
  void *fptr; //Pointer to the function
} strategy_t;

typedef struct {
  unsigned int count;
  unsigned int allocated;//How much memory is allocated
  strategy_t* strategies;
} strategy_list_t;

#define STRATEGY_LIST_ALLOC_SIZE 16

typedef struct {
  const char *strategy_type;
  void **fptr;
} strategy_to_select_t;

typedef struct {
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

    bool hyper_threading;
  } intel_flags;
  
  struct {
    int altivec;
  } powerpc_flags;
  
  struct {
    int neon;
  } arm_flags;

  int logical_cpu_count;
  int physical_cpu_count;
} hardware_flags_t;

extern hardware_flags_t kvz_g_hardware_flags;
extern hardware_flags_t kvz_g_strategies_in_use;
extern hardware_flags_t kvz_g_strategies_available;

int kvz_strategyselector_init(int32_t cpuid, uint8_t bitdepth, uint8_t enable_logging_output);
int kvz_strategyselector_register(void *opaque, const char *type, const char *strategy_name, int priority, void *fptr);


//Strategy to include
#include "strategies/strategies-nal.h"
#include "strategies/strategies-picture.h"
#include "strategies/strategies-dct.h"
#include "strategies/strategies-ipol.h"
#include "strategies/strategies-quant.h"
#include "strategies/strategies-intra.h"
#include "strategies/strategies-sao.h"
#include "strategies/strategies-encode.h"

static const strategy_to_select_t strategies_to_select[] = {
  STRATEGIES_NAL_EXPORTS
  STRATEGIES_PICTURE_EXPORTS
  STRATEGIES_DCT_EXPORTS
  STRATEGIES_IPOL_EXPORTS
  STRATEGIES_QUANT_EXPORTS
  STRATEGIES_INTRA_EXPORTS
  STRATEGIES_SAO_EXPORTS
  STRATEGIES_ENCODE_EXPORTS
  { NULL, NULL },
};

#endif //STRATEGYSELECTOR_H_
