/**
 *  Part of HEVC Encoder
 *  By Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file context.c
    \brief Functions for context derication
    \author Marko Viitanen
    \date 2012-08    
  This file contains context derivation functions
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "config.h"
#include "encoder.h"
#include "cabac.h"
#include "context.h"

/* CONTEXTS */
/* ToDo: move somewhere else */
cabac_ctx *SplitFlagSCModel;
cabac_ctx g_SplitFlagSCModel[3]; /*<! \brief split flag context models */
cabac_ctx g_IntraModeSCModel;    /*<! \brief intra mode context models */
cabac_ctx g_ChromaPredSCModel[2];
cabac_ctx g_TransSubdivSCModel[4];    /*<! \brief intra mode context models */
cabac_ctx g_QtCbfSCModelY[3];
cabac_ctx g_QtCbfSCModelU[3];
//cabac_ctx g_QtCbfSCModelV[3];
cabac_ctx g_PartSizeSCModel;
cabac_ctx g_CUSigCoeffGroupSCModel[4];
cabac_ctx g_CUSigSCModel_luma[24];
cabac_ctx g_CUSigSCModel_chroma[24];
cabac_ctx g_CuCtxLastY_luma[15];
cabac_ctx g_CuCtxLastY_chroma[15];
cabac_ctx g_CuCtxLastX_luma[15];
cabac_ctx g_CuCtxLastX_chroma[15];
cabac_ctx g_CUOneSCModel_luma[24];


void init_contexts(encoder_control *encoder)
{
  uint16_t i;
  /* Initialize contexts */
  /* ToDo: add P/B slice */
  ctx_init(&g_SplitFlagSCModel[0], encoder->QP, INIT_SPLIT_FLAG[SLICE_I][0]);
  ctx_init(&g_SplitFlagSCModel[1], encoder->QP, INIT_SPLIT_FLAG[SLICE_I][1]);
  ctx_init(&g_SplitFlagSCModel[2], encoder->QP, INIT_SPLIT_FLAG[SLICE_I][2]);

  ctx_init(&g_IntraModeSCModel, encoder->QP, INIT_INTRA_PRED_MODE[SLICE_I]);

  ctx_init(&g_ChromaPredSCModel[0], encoder->QP, INIT_CHROMA_PRED_MODE[SLICE_I][0]);
  ctx_init(&g_ChromaPredSCModel[1], encoder->QP, INIT_CHROMA_PRED_MODE[SLICE_I][1]);
  

  for(i = 0; i < 4; i++)
  {
    ctx_init(&g_TransSubdivSCModel[i], encoder->QP, INIT_TRANS_SUBDIV_FLAG[SLICE_I][i]);
    ctx_init(&g_CUSigCoeffGroupSCModel[i], encoder->QP, INIT_SIG_CG_FLAG[SLICE_I][i]);
  }
  for(i = 0; i < 3; i++)
  {
    ctx_init(&g_QtCbfSCModelY[i], encoder->QP, INIT_QT_CBF[SLICE_I][i]);
    ctx_init(&g_QtCbfSCModelU[i], encoder->QP, INIT_QT_CBF[SLICE_I][i+3]);
    //cxt_init(&g_QtCbfSCModelV[i], encoder->QP, INIT_QT_CBF[SLICE_I][i]);
  }
  for(i = 0; i < 15; i++)
  {
    ctx_init(&g_CuCtxLastY_luma[i], encoder->QP, INIT_LAST[SLICE_I][i] );
    ctx_init(&g_CuCtxLastX_luma[i], encoder->QP, INIT_LAST[SLICE_I][i] );

    ctx_init(&g_CuCtxLastY_chroma[i], encoder->QP, INIT_LAST[SLICE_I][i+15] );
    ctx_init(&g_CuCtxLastX_chroma[i], encoder->QP, INIT_LAST[SLICE_I][i+15] );
  }
  
  for(i = 0; i < 24; i++)
  {
    ctx_init(&g_CUOneSCModel_luma[i], encoder->QP, INIT_ONE_FLAG[SLICE_I][i]);

    ctx_init(&g_CUSigSCModel_luma[i], encoder->QP, INIT_SIG_FLAG[SLICE_I][i]);
    if(i < 21)
    {
      ctx_init(&g_CUSigSCModel_chroma[i], encoder->QP, INIT_SIG_FLAG[SLICE_I][i+24]);
    }
  }

}



//uint8_t get_context_coeff_abs_significant_flag(uint8_t 


  /** Pattern decision for context derivation process of significant_coeff_flag
 * \param sigCoeffGroupFlag pointer to prior coded significant coeff group
 * \param posXCG column of current coefficient group
 * \param posYCG row of current coefficient group
 * \param width width of the block
 * \param height height of the block
 * \returns pattern for current coefficient group
 */
 /*
Int  TComTrQuant::calcPatternSigCtx( const UInt* sigCoeffGroupFlag, UInt posXCG, UInt posYCG, Int width, Int height )
{
  if( width == 4 && height == 4 ) return -1;

  UInt sigRight = 0;
  UInt sigLower = 0;

  width >>= 2;
  height >>= 2;
  if( posXCG < width - 1 )
  {
    sigRight = (sigCoeffGroupFlag[ posYCG * width + posXCG + 1 ] != 0);
  }
  if (posYCG < height - 1 )
  {
    sigLower = (sigCoeffGroupFlag[ (posYCG  + 1 ) * width + posXCG ] != 0);
  }
  return sigRight + (sigLower<<1);
}
*/

/** Context derivation process of coeff_abs_significant_flag
 * \param patternSigCtx pattern for current coefficient group
 * \param posX column of current scan position
 * \param posY row of current scan position
 * \param blockType log2 value of block size if square block, or 4 otherwise
 * \param width width of the block
 * \param height height of the block
 * \param textureType texture type (TEXT_LUMA...)
 * \returns ctxInc for current scan position
 */
/*
Int TComTrQuant::getSigCtxInc    (
                                   Int                             patternSigCtx,
                                   UInt                            scanIdx,
                                   Int                             posX,
                                   Int                             posY,
                                   Int                             blockType,
                                   Int                             width
                                  ,Int                             height
                                  ,TextType                        textureType
                                  )
{
  const Int ctxIndMap[16] =
  {
    0, 1, 4, 5,
    2, 3, 4, 5,
    6, 6, 8, 8,
    7, 7, 8, 8
  };

  if( posX + posY == 0 )
  {
    return 0;
  }

  if ( blockType == 2 )
  {
    return ctxIndMap[ 4 * posY + posX ];
  }

  Int offset = blockType == 3 ? (scanIdx==SCAN_DIAG ? 9 : 15) : (textureType == TEXT_LUMA ? 21 : 12);


  Int posXinSubset = posX-((posX>>2)<<2);
  Int posYinSubset = posY-((posY>>2)<<2);
  Int cnt = 0;
  if(patternSigCtx==0)
  {
    cnt = posXinSubset+posYinSubset<=2 ? (posXinSubset+posYinSubset==0 ? 2 : 1) : 0;
  }
  else if(patternSigCtx==1)
  {
    cnt = posYinSubset<=1 ? (posYinSubset==0 ? 2 : 1) : 0;
  }
  else if(patternSigCtx==2)
  {
    cnt = posXinSubset<=1 ? (posXinSubset==0 ? 2 : 1) : 0;
  }
  else
  {
    cnt = 2;
  }

  return (( textureType == TEXT_LUMA && ((posX>>2) + (posY>>2)) > 0 ) ? 3 : 0) + offset + cnt;
}
*/
