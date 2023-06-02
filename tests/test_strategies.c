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

#include "test_strategies.h"

#include "src/strategyselector.h"


strategy_list_t strategies;


void init_test_strategies()
{
  strategies.allocated = 0;
  strategies.count = 0;
  strategies.strategies = NULL;

  // Init strategyselector because it sets hardware flags.
  kvz_strategyselector_init(1, KVZ_BIT_DEPTH, 1);

  // Collect all strategies to be tested.
  if (!kvz_strategy_register_picture(&strategies, KVZ_BIT_DEPTH)) {
    fprintf(stderr, "strategy_register_picture failed!\n");
    return;
  }

  if (!kvz_strategy_register_dct(&strategies, KVZ_BIT_DEPTH)) {
    fprintf(stderr, "strategy_register_dct failed!\n");
    return;
  }

  if (!kvz_strategy_register_quant(&strategies, KVZ_BIT_DEPTH)) {
    fprintf(stderr, "strategy_register_quant failed!\n");
    return;
  }
}
