#ifndef RDO_H_
#define RDO_H_
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
 * \ingroup Compression
 * \file
 * Rate-Distortion Optimization related functionality.
 */

#include "cabac.h"
#include "cu.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "search_inter.h"


// This struct is for passing data to kvz_rdoq_sign_hiding
struct kvz_sh_rates_t {
  // Bit cost of increasing rate by one.
  int32_t inc[32 * 32];
  // Bit cost of decreasing rate by one.
  int32_t dec[32 * 32];
  // Bit cost of going from zero to one.
  int32_t sig_coeff_inc[32 * 32];
  // Coeff minus quantized coeff.
  int32_t quant_delta[32 * 32];
};


extern const uint32_t kvz_g_go_rice_range[5];
extern const uint32_t kvz_g_go_rice_prefix_len[5];

int kvz_init_rdcost_outfiles(const char *fn_template);
void kvz_close_rdcost_outfiles(void);

void  kvz_rdoq(encoder_state_t *state, coeff_t *coef, coeff_t *dest_coeff, int32_t width,
           int32_t height, int8_t type, int8_t scan_mode, int8_t block_type, int8_t tr_depth);

double kvz_get_coeff_cost(const encoder_state_t * const state,
                            const coeff_t *coeff,
                            int32_t width,
                            int32_t type,
                            int8_t scan_mode);

int32_t kvz_get_ic_rate(encoder_state_t *state, uint32_t abs_level, uint16_t ctx_num_one, uint16_t ctx_num_abs,
                    uint16_t abs_go_rice, uint32_t c1_idx, uint32_t c2_idx, int8_t type);
uint32_t kvz_get_coded_level(encoder_state_t * state, double* coded_cost, double* coded_cost0, double* coded_cost_sig,
                         int32_t level_double, uint32_t max_abs_level,
                         uint16_t ctx_num_sig, uint16_t ctx_num_one, uint16_t ctx_num_abs,
                         uint16_t abs_go_rice,
                         uint32_t c1_idx, uint32_t c2_idx,
                         int32_t q_bits,double temp, int8_t last, int8_t type);

kvz_mvd_cost_func kvz_calc_mvd_cost_cabac;

double kvz_get_mvd_coding_cost_cabac(const encoder_state_t* state,
                                     const cabac_data_t* cabac,
                                     int32_t mvd_hor,
                                     int32_t mvd_ver);

// Number of fixed point fractional bits used in the fractional bit table.
#define CTX_FRAC_BITS 15
#define CTX_FRAC_ONE_BIT (1 << CTX_FRAC_BITS)
#define CTX_FRAC_HALF_BIT (1 << (CTX_FRAC_BITS - 1))

extern const uint32_t kvz_entropy_bits[128];
#define CTX_ENTROPY_BITS(ctx, val) kvz_entropy_bits[(ctx)->uc_state ^ (val)]

// Floating point fractional bits, derived from kvz_entropy_bits

#endif
