#ifndef CABAC_H_
#define CABAC_H_
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
 * \ingroup CABAC
 * \file
 * Coding bins using CABAC.
 */

#include "global.h" // IWYU pragma: keep

#include "bitstream.h"

struct encoder_state_t;

// Types
typedef struct
{
  uint8_t  uc_state;
} cabac_ctx_t;

typedef struct
{
  cabac_ctx_t *cur_ctx;
  uint32_t   low;
  uint32_t   range;
  uint32_t   buffered_byte;
  int32_t    num_buffered_bytes;
  int32_t    bits_left;
  int8_t     only_count : 4;
  int8_t     update : 4;
  bitstream_t *stream;

  // CONTEXTS
  struct {
    cabac_ctx_t sao_merge_flag_model;
    cabac_ctx_t sao_type_idx_model;
    cabac_ctx_t split_flag_model[3]; //!< \brief split flag context models
    cabac_ctx_t intra_mode_model;    //!< \brief intra mode context models
    cabac_ctx_t chroma_pred_model[2];
    cabac_ctx_t inter_dir[5];
    cabac_ctx_t trans_subdiv_model[3]; //!< \brief intra mode context models
    cabac_ctx_t qt_cbf_model_luma[4];
    cabac_ctx_t qt_cbf_model_chroma[4];
    cabac_ctx_t cu_qp_delta_abs[4];
    cabac_ctx_t part_size_model[4];
    cabac_ctx_t cu_sig_coeff_group_model[4];
    cabac_ctx_t cu_sig_model_luma[27];
    cabac_ctx_t cu_sig_model_chroma[15];
    cabac_ctx_t cu_ctx_last_y_luma[15];
    cabac_ctx_t cu_ctx_last_y_chroma[15];
    cabac_ctx_t cu_ctx_last_x_luma[15];
    cabac_ctx_t cu_ctx_last_x_chroma[15];
    cabac_ctx_t cu_one_model_luma[16];
    cabac_ctx_t cu_one_model_chroma[8];
    cabac_ctx_t cu_abs_model_luma[4];
    cabac_ctx_t cu_abs_model_chroma[2];
    cabac_ctx_t cu_pred_mode_model;
    cabac_ctx_t cu_skip_flag_model[3];
    cabac_ctx_t cu_merge_idx_ext_model;
    cabac_ctx_t cu_merge_flag_ext_model;
    cabac_ctx_t cu_transquant_bypass;
    cabac_ctx_t cu_mvd_model[2];
    cabac_ctx_t cu_ref_pic_model[2];
    cabac_ctx_t mvp_idx_model[2];
    cabac_ctx_t cu_qt_root_cbf_model;
    cabac_ctx_t transform_skip_model_luma;
    cabac_ctx_t transform_skip_model_chroma;
  } ctx;
} cabac_data_t;


// Globals
extern const uint8_t kvz_g_auc_next_state_mps[128];
extern const uint8_t kvz_g_auc_next_state_lps[128];
extern const uint8_t kvz_g_auc_lpst_table[64][4];
extern const uint8_t kvz_g_auc_renorm_table[32];


// Functions
void kvz_cabac_start(cabac_data_t *data);
void kvz_cabac_encode_bin(cabac_data_t *data, uint32_t bin_value);
void kvz_cabac_encode_bin_ep(cabac_data_t *data, uint32_t bin_value);
void kvz_cabac_encode_bins_ep(cabac_data_t *data, uint32_t bin_values, int num_bins);
void kvz_cabac_encode_bin_trm(cabac_data_t *data, uint8_t bin_value);
void kvz_cabac_write(cabac_data_t *data);
void kvz_cabac_finish(cabac_data_t *data);
int kvz_cabac_write_coeff_remain(cabac_data_t* cabac, uint32_t symbol,
                                 uint32_t r_param);
void kvz_cabac_write_coeff_remain_encry(struct encoder_state_t * const state, cabac_data_t * const cabac, const uint32_t symbol,
                                        const uint32_t r_param, int32_t base_level);
uint32_t kvz_cabac_write_ep_ex_golomb(struct encoder_state_t * const state, cabac_data_t *data,
                                  uint32_t symbol, uint32_t count);
void kvz_cabac_write_unary_max_symbol(cabac_data_t *data, cabac_ctx_t *ctx,
                                      uint32_t symbol, int32_t offset,
                                      uint32_t max_symbol, double* bits_out);
void kvz_cabac_write_unary_max_symbol_ep(cabac_data_t *data, unsigned int symbol, unsigned int max_symbol);

extern const float kvz_f_entropy_bits[128];
#define CTX_ENTROPY_FBITS(ctx, val) kvz_f_entropy_bits[(ctx)->uc_state ^ (val)]

#define CABAC_FBITS_UPDATE(cabac, ctx, val, bits, name) do { \
  if((cabac)->only_count) (bits) += kvz_f_entropy_bits[(ctx)->uc_state ^ (val)]; \
  if((cabac)->update) {\
    (cabac)->cur_ctx = ctx;\
    CABAC_BIN((cabac), (val), (name));\
  } \
} while(0)

// Macros
#define CTX_STATE(ctx) ((ctx)->uc_state >> 1)
#define CTX_MPS(ctx) ((ctx)->uc_state & 1)
#define CTX_UPDATE_LPS(ctx) { (ctx)->uc_state = kvz_g_auc_next_state_lps[ (ctx)->uc_state ]; }
#define CTX_UPDATE_MPS(ctx) { (ctx)->uc_state = kvz_g_auc_next_state_mps[ (ctx)->uc_state ]; }


#ifdef VERBOSE
  #define CABAC_BIN(data, value, name) { \
    uint32_t prev_state = (data)->cur_ctx->uc_state; \
    kvz_cabac_encode_bin((data), (value)); \
    if(!(data)->only_count)  printf("%s = %u, state = %u -> %u MPS = %u\n", \
           (name), (uint32_t)(value), prev_state, (data)->cur_ctx->uc_state, CTX_MPS((data)->cur_ctx)); }

  #define CABAC_BINS_EP(data, value, bins, name) { \
    uint32_t prev_state = (data)->cur_ctx->uc_state; \
    kvz_cabac_encode_bins_ep((data), (value), (bins)); \
    if(!(data)->only_count) printf("%s = %u(%u bins), state = %u -> %u\n", \
           (name), (uint32_t)(value), (bins), prev_state, (data)->cur_ctx->uc_state); }

  #define CABAC_BIN_EP(data, value, name) { \
    uint32_t prev_state = (data)->cur_ctx->uc_state; \
    kvz_cabac_encode_bin_ep((data), (value)); \
    if(!(data)->only_count) printf("%s = %u, state = %u -> %u\n", \
           (name), (uint32_t)(value), prev_state, (data)->cur_ctx->uc_state); }
#else
  #define CABAC_BIN(data, value, name) \
    kvz_cabac_encode_bin((data), (value));
  #define CABAC_BINS_EP(data, value, bins, name) \
    kvz_cabac_encode_bins_ep((data), (value), (bins));
  #define CABAC_BIN_EP(data, value, name) \
    kvz_cabac_encode_bin_ep((data), (value));
#endif

#endif
