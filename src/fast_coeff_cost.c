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
 
#include "fast_coeff_cost.h"
#include "kvazaar.h"
#include "encoderstate.h"

// Note: Assumes that costs are non-negative, for pretty obvious reasons
static uint16_t to_q88(float f)
{
  return (uint16_t)(f * 256.0f + 0.5f);
}

static uint64_t to_4xq88(const double f[4])
{
  int i;
  uint64_t result = 0;

  for (i = 3; i >= 0; i--) {
    result <<= 16;
    result |= to_q88(f[i]);
  }
  return result;
}

int kvz_fast_coeff_table_parse(fast_coeff_table_t *fast_coeff_table, FILE *fast_coeff_table_f)
{
  int i;
  uint64_t *wts_by_qp = fast_coeff_table->wts_by_qp;

  for (i = 0; i < MAX_FAST_COEFF_COST_QP; i++) {
    double curr_wts[4];

    if (fscanf(fast_coeff_table_f, "%lf %lf %lf %lf\n", curr_wts + 0,
                                                    curr_wts + 1,
                                                    curr_wts + 2,
                                                    curr_wts + 3) != 4) {
      return 1;
    }
    wts_by_qp[i] = to_4xq88(curr_wts);
  }
  return 0;
}

void kvz_fast_coeff_use_default_table(fast_coeff_table_t *fast_coeff_table)
{
  int i;
  uint64_t *wts_by_qp = fast_coeff_table->wts_by_qp;

  for (i = 0; i < MAX_FAST_COEFF_COST_QP; i++) {
    wts_by_qp[i] = to_4xq88(default_fast_coeff_cost_wts[i]);
  }
}

uint64_t kvz_fast_coeff_get_weights(const encoder_state_t *state)
{
  const fast_coeff_table_t *table = &(state->encoder_control->fast_coeff_table);
  return table->wts_by_qp[state->qp];
}
