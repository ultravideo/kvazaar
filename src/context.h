/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file context.h
    \brief Context
    \author Marko Viitanen
    \date 2012-08
    
    Context derivation function headers
*/
#ifndef __CONTEXT_H
#define __CONTEXT_H

#include "global.h"

#include "encoder.h"
#include "cabac.h"


void init_contexts(encoder_control *encoder, int8_t SLICE);

int32_t  context_calcPatternSigCtx( const uint32_t* sigCoeffGroupFlag, uint32_t posXCG, uint32_t posYCG, int32_t width);

uint32_t context_get_sigCoeffGroup( uint32_t* uiSigCoeffGroupFlag,uint32_t uiCGPosX,
                                    uint32_t uiCGPosY,int32_t width);

int32_t context_getSigCtxInc(int32_t patternSigCtx,uint32_t scanIdx,int32_t posX,
                             int32_t posY,int32_t blockType,int32_t width,
                             int8_t textureType);


/* CONTEXTS */
extern cabac_ctx g_SplitFlagSCModel[3];
extern cabac_ctx g_IntraModeSCModel;
extern cabac_ctx g_ChromaPredSCModel[2];
extern cabac_ctx g_TransSubdivSCModel[3];
extern cabac_ctx g_QtCbfSCModelY[3];
extern cabac_ctx g_QtCbfSCModelU[3];
extern cabac_ctx g_PartSizeSCModel[4];
extern cabac_ctx g_CUSigCoeffGroupSCModel[4];
extern cabac_ctx g_CUSigSCModel_luma[27];
extern cabac_ctx g_CUSigSCModel_chroma[15];
extern cabac_ctx g_CuCtxLastY_luma[15];
extern cabac_ctx g_CuCtxLastY_chroma[15];
extern cabac_ctx g_CuCtxLastX_luma[15];
extern cabac_ctx g_CuCtxLastX_chroma[15];
extern cabac_ctx g_CUOneSCModel_luma[16];
extern cabac_ctx g_CUOneSCModel_chroma[8];
extern cabac_ctx g_cCUAbsSCModel_luma[4];
extern cabac_ctx g_cCUAbsSCModel_chroma[2];
extern cabac_ctx g_cCUPredModeSCModel;
extern cabac_ctx g_cCUSkipFlagSCModel[3];
extern cabac_ctx g_cCUMergeIdxExtSCModel;
extern cabac_ctx g_cCUMergeFlagExtSCModel;
extern cabac_ctx g_cCUMvdSCModel[2];
extern cabac_ctx g_cCURefPicSCModel[2];
extern cabac_ctx g_cMVPIdxSCModel[2];
extern cabac_ctx g_cCUQtRootCbfSCModel;
#define CNU 154

static const uint8_t 
INIT_QT_ROOT_CBF[3][1] = 
{
  {  79, }, 
  {  79, }, 
  { CNU, }, 
};

static const uint8_t 
INIT_MVP_IDX[3][2] =  
{
  { 168,  CNU, }, 
  { 168,  CNU, }, 
  { CNU,  CNU, }, 
};

static const uint8_t 
INIT_REF_PIC[3][2] =  
{
  { 153,  153 }, 
  { 153,  153 }, 
  { CNU,  CNU }, 
};

static const uint8_t 
INIT_MVD[3][2] =  
{
  { 169,  198, }, 
  { 140,  198, }, 
  { CNU,  CNU, }, 
};

static const uint8_t
INIT_MERGE_FLAG_EXT[3][1] = 
{
  { 154, }, 
  { 110, }, 
  { CNU, }
};

static const uint8_t 
INIT_MERGE_IDX_EXT[3][1] =  
{
  { 137, }, 
  { 122, }, 
  { CNU, }
};

static const uint8_t 
INIT_SKIP_FLAG[3][3] =  
{
  { 197,  185,  201, }, 
  { 197,  185,  201, }, 
  { CNU,  CNU,  CNU, }
};

static const uint8_t 
INIT_PRED_MODE[3][1] = 
{
  { 134, }, 
  { 149, }, 
  { CNU, }
};


static const uint8_t 
INIT_PART_SIZE[3][4] =  
{
  { 154,  139,  CNU,  CNU, }, 
  { 154,  139,  CNU,  CNU, }, 
  { 184,  CNU,  CNU,  CNU, }, 
};

static const uint8_t  INIT_SPLIT_FLAG[3][3] =  
                       { { 107,  139,  126 },
                         { 107,  139,  126 },
                         { 139,  141,  157 } };

static const uint8_t INIT_INTRA_PRED_MODE[3] = { 183,154,184 };

static const uint8_t INIT_CHROMA_PRED_MODE[3][2] = { { 152,  139 }, { 152,  139 }, {  63,  139 } };


static const uint8_t INIT_TRANS_SUBDIV_FLAG[3][3] = 
{
  { 224,  167,  122 }, 
  { 124,  138,   94 }, 
  { 153,  138,  138 }
};

static const uint8_t INIT_QT_CBF[3][6] =  
{
  { 153,  111,  CNU,  149,   92,  167 },
  { 153,  111,  CNU,  149,  107,  167 },
  { 111,  141,  CNU,   94,  138,  182 }
};

static const uint8_t INIT_SIG_CG_FLAG[3][4] =  
 {  { 121,  140,  61,  154  },  { 121,  140,  61,  154 }, {  91,  171,  134,  141  } };

static const uint8_t INIT_SIG_FLAG[3][42] = 
  {
   {170,154,139,153,139,123,123, 63,124,166,183,140,136,153,154,166,183,140,136,153,154,166,183,140,136,153,154,170,153,138,138,122,121,122,121,167,151,183,140,151,183,140,},
   {155,154,139,153,139,123,123,63,153,166,183,140,136,153,154,166,183,140,136,153,154,166,183,140,136,153,154,170,153,123,123,107,121,107,121,167,151,183,140,151,183,140,},
   {111,111,125,110,110,94,124,108,124,107,125,141,179,153,125,107,125,141,179,153,125,107,125,141,179,153,125,140,139,182,182,152,136,152,136,153,136,139,111,136,139,111,}
  };

static const uint8_t INIT_LAST[3][30] =  
{
  { 125,  110,  124,  110,   95,   94,  125,  111,  111,   79,  125,  126,  111,  111,   79,
    108,  123,   93,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU  }, 
  { 125,  110,   94,  110,   95,   79,  125,  111,  110,   78,  110,  111,  111,   95,   94,
    108,  123,  108,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU  }, 
  { 110,  110,  124,  125,  140,  153,  125,  127,  140,  109,  111,  143,  127,  111,   79,
    108,  123,   63,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU  }
};

static const uint8_t INIT_ONE_FLAG[3][24] = 
{
  {154,196,167,167,154,152,167,182,182,134,149,136,153,121,136,122,169,208,166,167,154,152,167,182},
  {154,196,196,167,154,152,167,182,182,134,149,136,153,121,136,137,169,194,166,167,154,167,137,182},
  {140, 92,137,138,140,152,138,139,153, 74,149, 92,139,107,122,152,140,179,166,182,140,227,122,197}
};

static const uint8_t INIT_ABS_FLAG[3][6] =  
{
  { 107,167, 91,107,107,167}, 
  { 107,167, 91,122,107,167}, 
  { 138,153,136,167,152,152}, 
};


#endif
