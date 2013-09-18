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


/* CONTEXTS */
cabac_ctx g_split_flag_model[3]; /*<! \brief split flag context models */
cabac_ctx g_intra_mode_model;    /*<! \brief intra mode context models */
cabac_ctx g_chroma_pred_model[2];
cabac_ctx g_trans_subdiv_model[3];    /*<! \brief intra mode context models */
cabac_ctx g_qt_cbf_model_luma[3];
cabac_ctx g_qt_cbf_model_chroma[3];
//cabac_ctx g_QtCbfSCModelV[3];
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

void init_contexts(encoder_control *encoder, int8_t slice)
{
  uint16_t i;

  /* Initialize contexts */
  /* TODO: add P/B slice */  
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
  


  for(i = 0; i < 4; i++)
  {    
    ctx_init(&g_cu_sig_coeff_group_model[i], encoder->QP, INIT_SIG_CG_FLAG[slice][i]);
    ctx_init(&g_cu_abs_model_luma[i], encoder->QP, INIT_ABS_FLAG[slice][i]);
    ctx_init(&g_part_size_model[i], encoder->QP, INIT_PART_SIZE[slice][i]);
  }
  for(i = 0; i < 3; i++)
  {
    ctx_init(&g_trans_subdiv_model[i], encoder->QP, INIT_TRANS_SUBDIV_FLAG[slice][i]);
    ctx_init(&g_qt_cbf_model_luma[i], encoder->QP, INIT_QT_CBF[slice][i]);
    ctx_init(&g_qt_cbf_model_chroma[i], encoder->QP, INIT_QT_CBF[slice][i+3]);
    //cxt_init(&g_QtCbfSCModelV[i], encoder->QP, INIT_QT_CBF[SLICE][i]);
  }

  for(i = 0; i < 8; i++)
  {
    ctx_init(&g_cu_one_model_chroma[i], encoder->QP, INIT_ONE_FLAG[slice][i+16]);
  }

  for(i = 0; i < 15; i++)
  {
    ctx_init(&g_cu_ctx_last_y_luma[i], encoder->QP, INIT_LAST[slice][i] );
    ctx_init(&g_cu_ctx_last_x_luma[i], encoder->QP, INIT_LAST[slice][i] );

    ctx_init(&g_cu_ctx_last_y_chroma[i], encoder->QP, INIT_LAST[slice][i+15] );
    ctx_init(&g_cu_ctx_last_x_chroma[i], encoder->QP, INIT_LAST[slice][i+15] );

    ctx_init(&g_cu_one_model_luma[i], encoder->QP, INIT_ONE_FLAG[slice][i]);
  }
  ctx_init(&g_cu_one_model_luma[15], encoder->QP, INIT_ONE_FLAG[slice][15]);
  
  for(i = 0; i < 27; i++)
  { 

    ctx_init(&g_cu_sig_model_luma[i], encoder->QP, INIT_SIG_FLAG[slice][i]);
    if(i < 15)
    {
      ctx_init(&g_cu_sig_model_chroma[i], encoder->QP, INIT_SIG_FLAG[slice][i+27]);
    }
  }

}

uint32_t context_get_sig_coeff_group( uint32_t* sig_coeff_group_flag,
                                    uint32_t pos_x,
                                    uint32_t pos_y,
                                    int32_t width)
{
  uint32_t uiRight = 0;
  uint32_t uiLower = 0;
  width >>= 2;
  if( pos_x < (uint32_t)width - 1 )
    uiRight = (sig_coeff_group_flag[ pos_y * width + pos_x + 1 ] != 0);

  if (pos_y < (uint32_t)width - 1 )
    uiLower = (sig_coeff_group_flag[ (pos_y  + 1 ) * width + pos_x ] != 0);

  return (uiRight || uiLower);
}


/*! 
 \brief Pattern decision for context derivation process of significant_coeff_flag
 \param sigCoeffGroupFlag pointer to prior coded significant coeff group
 \param posXCG column of current coefficient group
 \param posYCG row of current coefficient group
 \param width width of the block
 \param height height of the block
 \returns pattern for current coefficient group
*/
 
int32_t context_calc_pattern_sig_ctx( const uint32_t* sig_coeff_group_flag, uint32_t pos_x, uint32_t pos_y, int32_t width)
{
  if( width == 4) return -1;

  {
  uint32_t sigRight = 0;
  uint32_t sigLower = 0;
  width >>= 2;
  if( pos_x < (uint32_t)width - 1 )
  {
    sigRight = (sig_coeff_group_flag[ pos_y * width + pos_x + 1 ] != 0);
  }
  if (pos_y < (uint32_t)width - 1 )
  {
    sigLower = (sig_coeff_group_flag[ (pos_y  + 1 ) * width + pos_x ] != 0);
  }
  
  return sigRight + (sigLower<<1);
  }
}


/*! 
 \brief Context derivation process of coeff_abs_significant_flag
 \param patternSigCtx pattern for current coefficient group
 \param posX column of current scan position
 \param posY row of current scan position
 \param blockType log2 value of block size if square block, or 4 otherwise
 \param width width of the block
 \param textureType texture type (TEXT_LUMA...)
 \returns ctxInc for current scan position
*/

int32_t context_get_sig_ctx_inc(int32_t pattern_sig_ctx,uint32_t scan_idx,int32_t pos_x,
                             int32_t pos_y,int32_t block_type,int32_t width,
                             int8_t texture_type)
{
  const int32_t ctx_ind_map[16] =
  {
    0, 1, 4, 5,
    2, 3, 4, 5,
    6, 6, 8, 8,
    7, 7, 8, 8
  };

  if( pos_x + pos_y == 0 )
    return 0;

  if ( block_type == 2 )
    return ctx_ind_map[ 4 * pos_y + pos_x ];

  {
  int32_t cnt = 0;
  int32_t offset = block_type == 3 ? (scan_idx==SCAN_DIAG ? 9 : 15) : (texture_type == 0 ? 21 : 12);
  int32_t posXinSubset = pos_x-((pos_x>>2)<<2);
  int32_t posYinSubset = pos_y-((pos_y>>2)<<2);
  
  if(pattern_sig_ctx==0)
  {
    cnt = posXinSubset+posYinSubset<=2 ? (posXinSubset+posYinSubset==0 ? 2 : 1) : 0;
  }
  else if(pattern_sig_ctx==1)
  {
    cnt = posYinSubset<=1 ? (posYinSubset==0 ? 2 : 1) : 0;
  }
  else if(pattern_sig_ctx==2)
  {
    cnt = posXinSubset<=1 ? (posXinSubset==0 ? 2 : 1) : 0;
  }
  else
  {
    cnt = 2;
  }
  return (( texture_type == 0 && ((pos_x>>2) + (pos_y>>2)) > 0 ) ? 3 : 0) + offset + cnt;
  }
}

