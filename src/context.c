/**
 * \file
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
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

