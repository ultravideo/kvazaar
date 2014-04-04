#ifndef CONTEXT_H_
#define CONTEXT_H_
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
 * \brief Context derivation for CABAC.
 */

#include "global.h"

#include "encoder.h"


// Types
typedef struct
{
  uint8_t  uc_state;
} cabac_ctx;


// Functions
void ctx_init(cabac_ctx* ctx, uint32_t qp, uint32_t init_value);
void init_contexts(encoder_control *encoder, int8_t slice);
int32_t context_calc_pattern_sig_ctx( const uint32_t *sig_coeff_group_flag, uint32_t pos_x, uint32_t pos_y, int32_t width);

uint32_t context_get_sig_coeff_group( uint32_t *sig_coeff_group_flag,uint32_t pos_x, uint32_t pos_y,int32_t width);


int32_t context_get_sig_ctx_inc(int32_t pattern_sig_ctx,uint32_t scan_idx,int32_t pos_x,
                                int32_t pos_y,int32_t block_type, int8_t texture_type);

// CONTEXTS
extern cabac_ctx g_sao_merge_flag_model;
extern cabac_ctx g_sao_type_idx_model;
extern cabac_ctx g_split_flag_model[3];
extern cabac_ctx g_intra_mode_model;
extern cabac_ctx g_chroma_pred_model[2];
extern cabac_ctx g_trans_subdiv_model[3];
extern cabac_ctx g_qt_cbf_model_luma[3];
extern cabac_ctx g_qt_cbf_model_chroma[3];
extern cabac_ctx g_part_size_model[4];
extern cabac_ctx g_cu_sig_coeff_group_model[4];
extern cabac_ctx g_cu_sig_model_luma[27];
extern cabac_ctx g_cu_sig_model_chroma[15];
extern cabac_ctx g_cu_ctx_last_y_luma[15];
extern cabac_ctx g_cu_ctx_last_y_chroma[15];
extern cabac_ctx g_cu_ctx_last_x_luma[15];
extern cabac_ctx g_cu_ctx_last_x_chroma[15];
extern cabac_ctx g_cu_one_model_luma[16];
extern cabac_ctx g_cu_one_model_chroma[8];
extern cabac_ctx g_cu_abs_model_luma[4];
extern cabac_ctx g_cu_abs_model_chroma[2];
extern cabac_ctx g_cu_pred_mode_model;
extern cabac_ctx g_cu_skip_flag_model[3];
extern cabac_ctx g_cu_merge_idx_ext_model;
extern cabac_ctx g_cu_merge_flag_ext_model;
extern cabac_ctx g_cu_mvd_model[2];
extern cabac_ctx g_cu_ref_pic_model[2];
extern cabac_ctx g_mvp_idx_model[2];
extern cabac_ctx g_cu_qt_root_cbf_model;
extern cabac_ctx g_transform_skip_model_luma;
extern cabac_ctx g_transform_skip_model_chroma;
#define CNU 154

#endif
