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

#include "context.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"


// CONTEXTS
cabac_ctx g_sao_merge_flag_model;
cabac_ctx g_sao_type_idx_model;
cabac_ctx g_split_flag_model[3]; //!< \brief split flag context models
cabac_ctx g_intra_mode_model;    //!< \brief intra mode context models
cabac_ctx g_chroma_pred_model[2];
cabac_ctx g_trans_subdiv_model[3]; //!< \brief intra mode context models
cabac_ctx g_qt_cbf_model_luma[3];
cabac_ctx g_qt_cbf_model_chroma[3];
cabac_ctx g_part_size_model[4];
cabac_ctx g_cu_sig_coeff_group_model[4];
cabac_ctx g_cu_sig_model_luma[27];
cabac_ctx g_cu_sig_model_chroma[15];
cabac_ctx g_cu_ctx_last_y_luma[15];
cabac_ctx g_cu_ctx_last_y_chroma[15];
cabac_ctx g_cu_ctx_last_x_luma[15];
cabac_ctx g_cu_ctx_last_x_chroma[15];
cabac_ctx g_cu_one_model_luma[16];
cabac_ctx g_cu_one_model_chroma[8];
cabac_ctx g_cu_abs_model_luma[4];
cabac_ctx g_cu_abs_model_chroma[2];
cabac_ctx g_cu_pred_mode_model;
cabac_ctx g_cu_skip_flag_model[3];
cabac_ctx g_cu_merge_idx_ext_model;
cabac_ctx g_cu_merge_flag_ext_model;
cabac_ctx g_cu_mvd_model[2];
cabac_ctx g_cu_ref_pic_model[2];
cabac_ctx g_mvp_idx_model[2];
cabac_ctx g_cu_qt_root_cbf_model;


/**
 * \brief Initialize struct cabac_ctx.
 */
void ctx_init(cabac_ctx *ctx, uint32_t qp, uint32_t init_value)
{
  int slope = (init_value >> 4) * 5 - 45;
  int offset = ((init_value & 15) << 3) - 16;
  int init_state = MIN(MAX(1, ((slope * (int)qp) >> 4) + offset), 126);

  if (init_state >= 64) {
    ctx->uc_state = ((init_state - 64) << 1) + 1;
  } else {
    ctx->uc_state = (63 - init_state) << 1;
  }
  ctx->bins_coded = 0;
}

/**
 * \brief Initialize cabac context to be used for coding
 * \param encoder encoder control struct
 * \param slice type of slice we are coding (P/B/I)
 */
void init_contexts(encoder_control *encoder, int8_t slice)
{
  uint16_t i;

  // Initialize contexts
  ctx_init(&g_sao_merge_flag_model, encoder->QP, INIT_SAO_MERGE_FLAG[slice]);
  ctx_init(&g_sao_type_idx_model, encoder->QP, INIT_SAO_TYPE_IDX[slice]);

  ctx_init(&g_cu_merge_flag_ext_model, encoder->QP, INIT_MERGE_FLAG_EXT[slice][0]);
  ctx_init(&g_cu_merge_idx_ext_model, encoder->QP, INIT_MERGE_IDX_EXT[slice][0]);
  ctx_init(&g_cu_pred_mode_model, encoder->QP, INIT_PRED_MODE[slice][0]);

  ctx_init(&g_cu_skip_flag_model[0], encoder->QP, INIT_SKIP_FLAG[slice][0]);
  ctx_init(&g_cu_skip_flag_model[1], encoder->QP, INIT_SKIP_FLAG[slice][1]);
  ctx_init(&g_cu_skip_flag_model[2], encoder->QP, INIT_SKIP_FLAG[slice][2]);

  ctx_init(&g_split_flag_model[0], encoder->QP, INIT_SPLIT_FLAG[slice][0]);
  ctx_init(&g_split_flag_model[1], encoder->QP, INIT_SPLIT_FLAG[slice][1]);
  ctx_init(&g_split_flag_model[2], encoder->QP, INIT_SPLIT_FLAG[slice][2]);

  ctx_init(&g_intra_mode_model, encoder->QP, INIT_INTRA_PRED_MODE[slice]);  

  ctx_init(&g_chroma_pred_model[0], encoder->QP, INIT_CHROMA_PRED_MODE[slice][0]);
  ctx_init(&g_chroma_pred_model[1], encoder->QP, INIT_CHROMA_PRED_MODE[slice][1]);
  
  ctx_init(&g_cu_abs_model_chroma[0], encoder->QP, INIT_ABS_FLAG[slice][4]);
  ctx_init(&g_cu_abs_model_chroma[1], encoder->QP, INIT_ABS_FLAG[slice][5]);

  //TODO: ignore P/B contexts on intra frame
  ctx_init(&g_cu_qt_root_cbf_model, encoder->QP, INIT_QT_ROOT_CBF[slice][0]);

  ctx_init(&g_cu_mvd_model[0], encoder->QP, INIT_MVD[slice][0]);
  ctx_init(&g_cu_mvd_model[1], encoder->QP, INIT_MVD[slice][1]);
  ctx_init(&g_cu_ref_pic_model[0], encoder->QP, INIT_REF_PIC[slice][0]);
  ctx_init(&g_cu_ref_pic_model[1], encoder->QP, INIT_REF_PIC[slice][1]);
  ctx_init(&g_mvp_idx_model[0], encoder->QP, INIT_MVP_IDX[slice][0]);
  ctx_init(&g_mvp_idx_model[1], encoder->QP, INIT_MVP_IDX[slice][1]);
  
  for (i = 0; i < 4; i++) {    
    ctx_init(&g_cu_sig_coeff_group_model[i], encoder->QP, INIT_SIG_CG_FLAG[slice][i]);
    ctx_init(&g_cu_abs_model_luma[i], encoder->QP, INIT_ABS_FLAG[slice][i]);
    ctx_init(&g_part_size_model[i], encoder->QP, INIT_PART_SIZE[slice][i]);
  }
  for (i = 0; i < 3; i++) {
    ctx_init(&g_trans_subdiv_model[i], encoder->QP, INIT_TRANS_SUBDIV_FLAG[slice][i]);
    ctx_init(&g_qt_cbf_model_luma[i], encoder->QP, INIT_QT_CBF[slice][i]);
    ctx_init(&g_qt_cbf_model_chroma[i], encoder->QP, INIT_QT_CBF[slice][i+3]);
  }

  for (i = 0; i < 8; i++) {
    ctx_init(&g_cu_one_model_chroma[i], encoder->QP, INIT_ONE_FLAG[slice][i+16]);
  }

  for (i = 0; i < 15; i++) {
    ctx_init(&g_cu_ctx_last_y_luma[i], encoder->QP, INIT_LAST[slice][i] );
    ctx_init(&g_cu_ctx_last_x_luma[i], encoder->QP, INIT_LAST[slice][i] );

    ctx_init(&g_cu_ctx_last_y_chroma[i], encoder->QP, INIT_LAST[slice][i+15] );
    ctx_init(&g_cu_ctx_last_x_chroma[i], encoder->QP, INIT_LAST[slice][i+15] );

    ctx_init(&g_cu_one_model_luma[i], encoder->QP, INIT_ONE_FLAG[slice][i]);
  }
  ctx_init(&g_cu_one_model_luma[15], encoder->QP, INIT_ONE_FLAG[slice][15]);
  
  for (i = 0; i < 27; i++) {
    ctx_init(&g_cu_sig_model_luma[i], encoder->QP, INIT_SIG_FLAG[slice][i]);
    if(i < 15) ctx_init(&g_cu_sig_model_chroma[i], encoder->QP, INIT_SIG_FLAG[slice][i+27]);
  }

}


uint32_t context_get_sig_coeff_group( uint32_t *sig_coeff_group_flag,
                                      uint32_t pos_x,
                                      uint32_t pos_y,
                                      int32_t width)
{
  uint32_t uiRight = 0;
  uint32_t uiLower = 0;
  width >>= 2;
  if (pos_x < (uint32_t)width - 1) uiRight = (sig_coeff_group_flag[pos_y * width + pos_x + 1] != 0);
  if (pos_y < (uint32_t)width - 1) uiLower = (sig_coeff_group_flag[(pos_y  + 1 ) * width + pos_x] != 0);

  return uiRight || uiLower;
}


/** 
 * \brief Pattern decision for context derivation process of significant_coeff_flag
 * \param sig_coeff_group_flag pointer to prior coded significant coeff group
 * \param pos_x column of current coefficient group
 * \param pos_y row of current coefficient group
 * \param width width of the block
 * \returns pattern for current coefficient group
 */
 
int32_t context_calc_pattern_sig_ctx(const uint32_t *sig_coeff_group_flag, uint32_t pos_x, uint32_t pos_y, int32_t width)
{
  uint32_t sigRight = 0;
  uint32_t sigLower = 0;

  if (width == 4) return -1;

  width >>= 2;
  if (pos_x < (uint32_t)width - 1) sigRight = (sig_coeff_group_flag[pos_y * width + pos_x + 1] != 0);
  if (pos_y < (uint32_t)width - 1) sigLower = (sig_coeff_group_flag[(pos_y  + 1 ) * width + pos_x] != 0);
  
  return sigRight + (sigLower<<1);
}


/**
 * \brief Context derivation process of coeff_abs_significant_flag
 * \param pattern_sig_ctx pattern for current coefficient group
 * \param scan_idx pixel scan type in use
 * \param pos_x column of current scan position
 * \param pos_y row of current scan position
 * \param block_type log2 value of block size if square block, or 4 otherwise
 * \param width width of the block
 * \param texture_type texture type (TEXT_LUMA...)
 * \returns ctx_inc for current scan position
 */

int32_t context_get_sig_ctx_inc(int32_t pattern_sig_ctx, uint32_t scan_idx, int32_t pos_x,
                                int32_t pos_y, int32_t block_type, int32_t width, int8_t texture_type)
{
  const int32_t ctx_ind_map[16] =
  {
    0, 1, 4, 5,
    2, 3, 4, 5,
    6, 6, 8, 8,
    7, 7, 8, 8
  };

  int32_t cnt,offset,pos_x_in_subset,pos_y_in_subset;

  if (pos_x + pos_y == 0) return 0;

  if (block_type == 2) return ctx_ind_map[4 * pos_y + pos_x];

  cnt = 0;
  offset = (block_type == 3) ? ((scan_idx == SCAN_DIAG) ? 9 : 15) : ((texture_type == 0) ? 21 : 12);
  pos_x_in_subset = pos_x - ((pos_x>>2)<<2);
  pos_y_in_subset = pos_y - ((pos_y>>2)<<2);
  
  if (pattern_sig_ctx == 0) {
    cnt = (pos_x_in_subset + pos_y_in_subset <= 2) ? ((pos_x_in_subset + pos_y_in_subset==0) ? 2 : 1) : 0;  
  } else if (pattern_sig_ctx==1) {
    cnt = (pos_y_in_subset <= 1) ? ((pos_y_in_subset == 0) ? 2 : 1) : 0;
  } else if (pattern_sig_ctx==2) {
    cnt = (pos_x_in_subset <= 1) ? ((pos_x_in_subset == 0) ? 2 : 1) : 0;
  } else {
    cnt = 2;
  }
  return (( texture_type == 0 && ((pos_x>>2) + (pos_y>>2)) > 0 ) ? 3 : 0) + offset + cnt;
}

/*
 * Entropy bits to estimate coded bits in RDO / RDOQ (From HM 12.0)
 */
const uint32_t entropy_bits[128] =
{
  0x08000, 0x08000, 0x076da, 0x089a0, 0x06e92, 0x09340, 0x0670a, 0x09cdf, 0x06029, 0x0a67f, 0x059dd, 0x0b01f, 0x05413, 0x0b9bf, 0x04ebf, 0x0c35f,
  0x049d3, 0x0ccff, 0x04546, 0x0d69e, 0x0410d, 0x0e03e, 0x03d22, 0x0e9de, 0x0397d, 0x0f37e, 0x03619, 0x0fd1e, 0x032ee, 0x106be, 0x02ffa, 0x1105d,
  0x02d37, 0x119fd, 0x02aa2, 0x1239d, 0x02836, 0x12d3d, 0x025f2, 0x136dd, 0x023d1, 0x1407c, 0x021d2, 0x14a1c, 0x01ff2, 0x153bc, 0x01e2f, 0x15d5c,
  0x01c87, 0x166fc, 0x01af7, 0x1709b, 0x0197f, 0x17a3b, 0x0181d, 0x183db, 0x016d0, 0x18d7b, 0x01595, 0x1971b, 0x0146c, 0x1a0bb, 0x01354, 0x1aa5a,
  0x0124c, 0x1b3fa, 0x01153, 0x1bd9a, 0x01067, 0x1c73a, 0x00f89, 0x1d0da, 0x00eb7, 0x1da79, 0x00df0, 0x1e419, 0x00d34, 0x1edb9, 0x00c82, 0x1f759,
  0x00bda, 0x200f9, 0x00b3c, 0x20a99, 0x00aa5, 0x21438, 0x00a17, 0x21dd8, 0x00990, 0x22778, 0x00911, 0x23118, 0x00898, 0x23ab8, 0x00826, 0x24458,
  0x007ba, 0x24df7, 0x00753, 0x25797, 0x006f2, 0x26137, 0x00696, 0x26ad7, 0x0063f, 0x27477, 0x005ed, 0x27e17, 0x0059f, 0x287b6, 0x00554, 0x29156,
  0x0050e, 0x29af6, 0x004cc, 0x2a497, 0x0048d, 0x2ae35, 0x00451, 0x2b7d6, 0x00418, 0x2c176, 0x003e2, 0x2cb15, 0x003af, 0x2d4b5, 0x0037f, 0x2de55
};
