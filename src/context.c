/**
 *  Part of HEVC Encoder
 *  By Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file context.c
    \brief Functions for context derication
    \author Marko Viitanen
    \date 2013-04    
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
cabac_ctx g_SplitFlagSCModel[3]; /*<! \brief split flag context models */
cabac_ctx g_IntraModeSCModel;    /*<! \brief intra mode context models */
cabac_ctx g_ChromaPredSCModel[2];
cabac_ctx g_TransSubdivSCModel[3];    /*<! \brief intra mode context models */
cabac_ctx g_QtCbfSCModelY[3];
cabac_ctx g_QtCbfSCModelU[3];
//cabac_ctx g_QtCbfSCModelV[3];
cabac_ctx g_PartSizeSCModel[4];
cabac_ctx g_CUSigCoeffGroupSCModel[4];
cabac_ctx g_CUSigSCModel_luma[27];
cabac_ctx g_CUSigSCModel_chroma[15];
cabac_ctx g_CuCtxLastY_luma[15];
cabac_ctx g_CuCtxLastY_chroma[15];
cabac_ctx g_CuCtxLastX_luma[15];
cabac_ctx g_CuCtxLastX_chroma[15];
cabac_ctx g_CUOneSCModel_luma[16];
cabac_ctx g_CUOneSCModel_chroma[8];
cabac_ctx g_cCUAbsSCModel_luma[4];
cabac_ctx g_cCUAbsSCModel_chroma[2];
cabac_ctx g_cCUPredModeSCModel;
cabac_ctx g_cCUSkipFlagSCModel[3];
cabac_ctx g_cCUMergeIdxExtSCModel;
cabac_ctx g_cCUMergeFlagExtSCModel;
cabac_ctx g_cCUMvdSCModel[2];
cabac_ctx g_cCURefPicSCModel[2];
cabac_ctx g_cMVPIdxSCModel[2];


void init_contexts(encoder_control *encoder, int8_t SLICE)
{
  uint16_t i;

  /* Initialize contexts */
  /* ToDo: add P/B slice */  
  ctx_init(&g_cCUMergeFlagExtSCModel, encoder->QP, INIT_MERGE_FLAG_EXT[SLICE][0]);
  ctx_init(&g_cCUMergeIdxExtSCModel, encoder->QP, INIT_MERGE_IDX_EXT[SLICE][0]);
  ctx_init(&g_cCUPredModeSCModel, encoder->QP, INIT_PRED_MODE[SLICE][0]);

  ctx_init(&g_cCUSkipFlagSCModel[0], encoder->QP, INIT_SKIP_FLAG[SLICE][0]);
  ctx_init(&g_cCUSkipFlagSCModel[1], encoder->QP, INIT_SKIP_FLAG[SLICE][1]);
  ctx_init(&g_cCUSkipFlagSCModel[2], encoder->QP, INIT_SKIP_FLAG[SLICE][2]);


  ctx_init(&g_SplitFlagSCModel[0], encoder->QP, INIT_SPLIT_FLAG[SLICE][0]);
  ctx_init(&g_SplitFlagSCModel[1], encoder->QP, INIT_SPLIT_FLAG[SLICE][1]);
  ctx_init(&g_SplitFlagSCModel[2], encoder->QP, INIT_SPLIT_FLAG[SLICE][2]);

  ctx_init(&g_IntraModeSCModel, encoder->QP, INIT_INTRA_PRED_MODE[SLICE]);  

  ctx_init(&g_ChromaPredSCModel[0], encoder->QP, INIT_CHROMA_PRED_MODE[SLICE][0]);
  ctx_init(&g_ChromaPredSCModel[1], encoder->QP, INIT_CHROMA_PRED_MODE[SLICE][1]);
  
  ctx_init(&g_cCUAbsSCModel_chroma[0], encoder->QP, INIT_ABS_FLAG[SLICE][4]);
  ctx_init(&g_cCUAbsSCModel_chroma[1], encoder->QP, INIT_ABS_FLAG[SLICE][5]);

  //ToDo: ignore P/B contexts on intra frame
  ctx_init(&g_cCUMvdSCModel[0], encoder->QP, INIT_MVD[SLICE][0]);
  ctx_init(&g_cCUMvdSCModel[1], encoder->QP, INIT_MVD[SLICE][1]);
  ctx_init(&g_cCURefPicSCModel[0], encoder->QP, INIT_REF_PIC[SLICE][0]);
  ctx_init(&g_cCURefPicSCModel[1], encoder->QP, INIT_REF_PIC[SLICE][1]);
  ctx_init(&g_cMVPIdxSCModel[0], encoder->QP, INIT_MVP_IDX[SLICE][0]);
  ctx_init(&g_cMVPIdxSCModel[1], encoder->QP, INIT_MVP_IDX[SLICE][1]);
  


  for(i = 0; i < 4; i++)
  {    
    ctx_init(&g_CUSigCoeffGroupSCModel[i], encoder->QP, INIT_SIG_CG_FLAG[SLICE][i]);
    ctx_init(&g_cCUAbsSCModel_luma[i], encoder->QP, INIT_ABS_FLAG[SLICE][i]);
    ctx_init(&g_PartSizeSCModel[i], encoder->QP, INIT_PART_SIZE[SLICE][i]);
  }
  for(i = 0; i < 3; i++)
  {
    ctx_init(&g_TransSubdivSCModel[i], encoder->QP, INIT_TRANS_SUBDIV_FLAG[SLICE][i]);
    ctx_init(&g_QtCbfSCModelY[i], encoder->QP, INIT_QT_CBF[SLICE][i]);
    ctx_init(&g_QtCbfSCModelU[i], encoder->QP, INIT_QT_CBF[SLICE][i+3]);
    //cxt_init(&g_QtCbfSCModelV[i], encoder->QP, INIT_QT_CBF[SLICE][i]);
  }

  for(i = 0; i < 8; i++)
  {
    ctx_init(&g_CUOneSCModel_chroma[i], encoder->QP, INIT_ONE_FLAG[SLICE][i+16]);
  }

  for(i = 0; i < 15; i++)
  {
    ctx_init(&g_CuCtxLastY_luma[i], encoder->QP, INIT_LAST[SLICE][i] );
    ctx_init(&g_CuCtxLastX_luma[i], encoder->QP, INIT_LAST[SLICE][i] );

    ctx_init(&g_CuCtxLastY_chroma[i], encoder->QP, INIT_LAST[SLICE][i+15] );
    ctx_init(&g_CuCtxLastX_chroma[i], encoder->QP, INIT_LAST[SLICE][i+15] );

    ctx_init(&g_CUOneSCModel_luma[i], encoder->QP, INIT_ONE_FLAG[SLICE][i]);
  }
  ctx_init(&g_CUOneSCModel_luma[15], encoder->QP, INIT_ONE_FLAG[SLICE][15]);
  
  for(i = 0; i < 27; i++)
  { 

    ctx_init(&g_CUSigSCModel_luma[i], encoder->QP, INIT_SIG_FLAG[SLICE][i]);
    if(i < 15)
    {
      ctx_init(&g_CUSigSCModel_chroma[i], encoder->QP, INIT_SIG_FLAG[SLICE][i+27]);
    }
  }

}

uint32_t context_get_sigCoeffGroup( uint32_t* uiSigCoeffGroupFlag,
                                    uint32_t uiCGPosX,
                                    uint32_t uiCGPosY,                                    
                                    int32_t width)
{
  uint32_t uiRight = 0;
  uint32_t uiLower = 0;
  width >>= 2;
  if( uiCGPosX < (uint32_t)width - 1 )
    uiRight = (uiSigCoeffGroupFlag[ uiCGPosY * width + uiCGPosX + 1 ] != 0);

  if (uiCGPosY < (uint32_t)width - 1 )
    uiLower = (uiSigCoeffGroupFlag[ (uiCGPosY  + 1 ) * width + uiCGPosX ] != 0);

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
 
int32_t  context_calcPatternSigCtx( const uint32_t* sigCoeffGroupFlag, uint32_t posXCG, uint32_t posYCG, int32_t width)
{
  if( width == 4) return -1;

  {
  uint32_t sigRight = 0;
  uint32_t sigLower = 0;
  width >>= 2;
  if( posXCG < (uint32_t)width - 1 )
  {
    sigRight = (sigCoeffGroupFlag[ posYCG * width + posXCG + 1 ] != 0);
  }
  if (posYCG < (uint32_t)width - 1 )
  {
    sigLower = (sigCoeffGroupFlag[ (posYCG  + 1 ) * width + posXCG ] != 0);
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

int32_t context_getSigCtxInc(int32_t patternSigCtx,uint32_t scanIdx,int32_t posX,
                             int32_t posY,int32_t blockType,int32_t width,
                             int8_t textureType)
{
  const int32_t ctxIndMap[16] =
  {
    0, 1, 4, 5,
    2, 3, 4, 5,
    6, 6, 8, 8,
    7, 7, 8, 8
  };

  if( posX + posY == 0 )
    return 0;

  if ( blockType == 2 )
    return ctxIndMap[ 4 * posY + posX ];

  {
  int32_t cnt = 0;
  int32_t offset = blockType == 3 ? (scanIdx==SCAN_DIAG ? 9 : 15) : (textureType == 0 ? 21 : 12);
  int32_t posXinSubset = posX-((posX>>2)<<2);
  int32_t posYinSubset = posY-((posY>>2)<<2);
  
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
  return (( textureType == 0 && ((posX>>2) + (posY>>2)) > 0 ) ? 3 : 0) + offset + cnt;
  }
}

