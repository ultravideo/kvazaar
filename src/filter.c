/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file filter.c
    \brief filtering
    \author Marko Viitanen
    \date 2013-04
    
    Filtering functions
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "config.h"
#include "bitstream.h"
#include "picture.h"
#include "cabac.h"
#include "encoder.h"
#include "filter.h"

const uint8_t tctable_8x8[54] =
{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,5,5,6,6,7,8,9,10,11,13,14,16,18,20,22,24
};

const uint8_t betatable_8x8[52] =
{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,7,8,9,10,11,12,13,14,15,16,17,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64
};


void filter_luma( uint8_t* piSrc, int32_t iOffset, int32_t tc , int8_t sw, int8_t bPartPNoFilter, int8_t bPartQNoFilter, int32_t iThrCut, int8_t bFilterSecondP, int8_t bFilterSecondQ)
{
  int32_t delta;
  
  uint8_t m4  = piSrc[0];
  uint8_t m3  = piSrc[-iOffset];
  uint8_t m5  = piSrc[ iOffset];
  uint8_t m2  = piSrc[-iOffset*2];
  uint8_t m6  = piSrc[ iOffset*2];
  uint8_t m1  = piSrc[-iOffset*3];
  uint8_t m7  = piSrc[ iOffset*3];
  uint8_t m0  = piSrc[-iOffset*4];

  if (sw)
  {
    piSrc[-iOffset]   = CLIP(m3-2*tc, m3+2*tc, ((m1 + 2*m2 + 2*m3 + 2*m4 + m5 + 4) >> 3));
    piSrc[0]          = CLIP(m4-2*tc, m4+2*tc, ((m2 + 2*m3 + 2*m4 + 2*m5 + m6 + 4) >> 3));
    piSrc[-iOffset*2] = CLIP(m2-2*tc, m2+2*tc, ((m1 + m2 + m3 + m4 + 2)>>2));
    piSrc[ iOffset]   = CLIP(m5-2*tc, m5+2*tc, ((m3 + m4 + m5 + m6 + 2)>>2));
    piSrc[-iOffset*3] = CLIP(m1-2*tc, m1+2*tc, ((2*m0 + 3*m1 + m2 + m3 + m4 + 4 )>>3));
    piSrc[ iOffset*2] = CLIP(m6-2*tc, m6+2*tc, ((m3 + m4 + m5 + 3*m6 + 2*m7 +4 )>>3));
  }
  else
  {
    // Weak filter
    delta = (9*(m4-m3) -3*(m5-m2) + 8)>>4 ;

    if ( abs(delta) < iThrCut )
    {
      int32_t tc2 = tc>>1;
      delta = CLIP(-tc, tc, delta);        
      piSrc[-iOffset] = CLIP(0,(1 << g_bitDepth)-1,(m3+delta));
      piSrc[0] = CLIP(0,(1 << g_bitDepth)-1,(m4-delta));

      
      if(bFilterSecondP)
      {
        int32_t delta1 = CLIP(-tc2, tc2, (( ((m1+m3+1)>>1)- m2+delta)>>1));
        piSrc[-iOffset*2] = CLIP(0,(1 << g_bitDepth)-1,(m2+delta1));
      }
      if(bFilterSecondQ)
      {
        int32_t delta2 = CLIP(-tc2, tc2, (( ((m6+m4+1)>>1)- m5-delta)>>1));
        piSrc[ iOffset] = CLIP(0,(1 << g_bitDepth)-1,(m5+delta2));
      }
    }
  }

  if(bPartPNoFilter)
  {
    piSrc[-iOffset] = m3;
    piSrc[-iOffset*2] = m2;
    piSrc[-iOffset*3] = m1;
  }
  if(bPartQNoFilter)
  {
    piSrc[0] = m4;
    piSrc[ iOffset] = m5;
    piSrc[ iOffset*2] = m6;
  }
}

void filter_deblock_edge_luma(encoder_control* encoder, int32_t xpos, int32_t ypos, int8_t depth, int32_t edge, int8_t dir)
{
  int i,iIdx;
  int32_t iStride = encoder->in.cur_pic.width;
  int32_t iOffset = 0;
  int32_t betaOffsetDiv2 = encoder->betaOffsetdiv2;
  int32_t tcOffsetDiv2   = encoder->tcOffsetdiv2;
  const int8_t scu_width       = (LCU_WIDTH>>MAX_DEPTH);  
  const int8_t scu_width_log2  = TOBITS(scu_width);
  int8_t uiNumParts = 1;
  int8_t uiBs       = 1; /* Filter strength */
  /* ToDo: support 10+bits */
  uint8_t* src      = &encoder->in.cur_pic.yRecData[xpos+ypos*iStride];
  uint8_t* piTmpSrc = src;
  int32_t iSrcStep;
  CU_info* cu = &encoder->in.cur_pic.CU[0][(xpos>>scu_width_log2) + (ypos>>scu_width_log2)*(encoder->in.width>>scu_width_log2)];
  
  if(dir == EDGE_VER)
  {
    iOffset = 1;
    iSrcStep = iStride;
    piTmpSrc += edge*scu_width;
  }
  else
  {
    iOffset = iStride;
    iSrcStep = 1;
    piTmpSrc += edge*scu_width*iStride;
  }
  
  /* For each subpart */
  for(iIdx = 0; iIdx < uiNumParts; iIdx++)
  {
    int32_t iQP            = encoder->QP;
    int32_t iBitdepthScale = 1 << (g_bitDepth-8);
    int32_t iIndexTC       = CLIP(0, 51+2, (int32_t)(iQP + 2*(uiBs-1) + (tcOffsetDiv2 << 1)));
    int32_t iIndexB        = CLIP(0, 51, iQP + (betaOffsetDiv2 << 1));
    int32_t iTc            = tctable_8x8[iIndexTC]*iBitdepthScale;
    int32_t iBeta          = betatable_8x8[iIndexB]*iBitdepthScale;
    int32_t iSideThreshold = (iBeta+(iBeta>>1))>>3;
    int32_t iThrCut        = iTc*10;
    uint32_t uiBlocksInPart= scu_width / 4 ? scu_width / 4 : 1;
    uint32_t iBlkIdx;

    for (iBlkIdx = 0; iBlkIdx < uiBlocksInPart; iBlkIdx++)
    {
      uint8_t* piTmpSrcShift;
      int32_t dp0,dq0,dp3,dq3,d0,d3,dp,dq,d;

      /* Check conditions for filtering */
      piTmpSrcShift = piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+0);
      dp0 = abs( piTmpSrcShift[-iOffset*3] - 2*piTmpSrcShift[-iOffset*2] + piTmpSrcShift[-iOffset] );
      piTmpSrcShift = piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+0);
      dq0 = abs( piTmpSrcShift[0] - 2*piTmpSrcShift[iOffset] + piTmpSrcShift[iOffset*2] );
      piTmpSrcShift = piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+3);
      dp3 = abs( piTmpSrcShift[-iOffset*3] - 2*piTmpSrcShift[-iOffset*2] + piTmpSrcShift[-iOffset] );
      piTmpSrcShift = piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+3);
      dq3 = abs( piTmpSrcShift[0] - 2*piTmpSrcShift[iOffset] + piTmpSrcShift[iOffset*2] );
      d0 = dp0 + dq0;
      d3 = dp3 + dq3;        
      dp = dp0 + dp3;
      dq = dq0 + dq3;
      d  =  d0 + d3;

      #if ENABLE_PCM == 1
      //ToDo: add PCM deblocking
      #endif
      if (d < iBeta)
      { 
        #define useStrongFiltering(offset,d,beta,tc,src) ( ((abs(src[-offset*4]-src[-offset]) + abs(src[-offset*3]-src[0])) < (beta>>3)) && (d<(beta>>2)) && ( abs(src[-offset]-src[0]) < ((tc*5+1)>>1)) )
        int8_t bFilterP = (dp < iSideThreshold)?1:0;
        int8_t bFilterQ = (dq < iSideThreshold)?1:0;          
        int8_t sw = useStrongFiltering( iOffset, 2*d0, iBeta, iTc, (piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+0)))
                            && useStrongFiltering( iOffset, 2*d3, iBeta, iTc, (piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+3)));
        for (i = 0; i < 8/2; i++)
        {
          filter_luma( piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+i), iOffset, iTc, sw, 0, 0, iThrCut, bFilterP, bFilterQ);
        }
      }
    }
  }
}

void filter_deblock_edge_chroma(encoder_control* encoder,int32_t idx, int32_t xpos, int32_t ypos, int8_t depth, int32_t edge, int8_t dir)
{
  //int i,iIdx;
  int32_t iStride = encoder->in.cur_pic.width;
  int32_t iOffset = 0;  
  int32_t tcOffsetDiv2   = encoder->betaOffsetdiv2;
  const int8_t scu_width       = (LCU_WIDTH>>(MAX_DEPTH+1));
  const int8_t scu_width_log2  = TOBITS(scu_width);
  int8_t uiNumParts = 1;
  int8_t uiBs       = 1; /* Filter strength */
  /* ToDo: support 10+bits */
  uint8_t* src      = &encoder->in.cur_pic.yRecData[xpos+ypos*iStride];
  uint8_t* piTmpSrc = src;
  //int32_t iSrcStep;
  CU_info* cu = &encoder->in.cur_pic.CU[0][(xpos>>scu_width_log2) + (ypos>>scu_width_log2)*(encoder->in.width>>scu_width_log2)];

  uint32_t uiEdgeNumInLCUVert = g_auiZscanToRaster[idx]%(1<<MAX_DEPTH) + edge;
  uint32_t uiEdgeNumInLCUHor = g_auiZscanToRaster[idx]/(1<<MAX_DEPTH) + edge;
  
  if ( (scu_width < 8) && (( (uiEdgeNumInLCUVert%(8/scu_width))&&(dir==0) ) || ( (uiEdgeNumInLCUHor%(8/scu_width))&& dir ) ))
  {
    return;
  }
  /*
  
  if(dir == EDGE_VER)
  {
    iOffset = 1;
    iSrcStep = iStride;
    piTmpSrc += edge*scu_width;
  }
  else
  {
    iOffset = iStride;
    iSrcStep = 1;
    piTmpSrc += edge*scu_width*iStride;
  }
  
  // For each subpart
  for(iIdx = 0; iIdx < uiNumParts; iIdx++)
  {
    int32_t iQP            = encoder->QP;
    int32_t iBitdepthScale = 1 << (g_bitDepth-8);
    int32_t iIndexTC       = CLIP(0, 51+2, (int32_t)(iQP + 2*(uiBs-1) + (tcOffsetDiv2 << 1)));
    int32_t iIndexB        = CLIP(0, 51, iQP + (betaOffsetDiv2 << 1));
    int32_t iTc            = tctable_8x8[iIndexTC]*iBitdepthScale;
    int32_t iBeta          = betatable_8x8[iIndexB]*iBitdepthScale;
    int32_t iSideThreshold = (iBeta+(iBeta>>1))>>3;
    int32_t iThrCut        = iTc*10;
    uint32_t uiBlocksInPart= scu_width / 4 ? scu_width / 4 : 1;
    uint32_t iBlkIdx;

    for (iBlkIdx = 0; iBlkIdx < uiBlocksInPart; iBlkIdx++)
    {
      uint8_t* piTmpSrcShift;
      int32_t dp0,dq0,dp3,dq3,d0,d3,dp,dq,d;


      piTmpSrcShift = piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+0);
      dp0 = abs( piTmpSrcShift[-iOffset*3] - 2*piTmpSrcShift[-iOffset*2] + piTmpSrcShift[-iOffset] );
      piTmpSrcShift = piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+0);
      dq0 = abs( piTmpSrcShift[0] - 2*piTmpSrcShift[iOffset] + piTmpSrcShift[iOffset*2] );
      piTmpSrcShift = piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+3);
      dp3 = abs( piTmpSrcShift[-iOffset*3] - 2*piTmpSrcShift[-iOffset*2] + piTmpSrcShift[-iOffset] );
      piTmpSrcShift = piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+3);
      dq3 = abs( piTmpSrcShift[0] - 2*piTmpSrcShift[iOffset] + piTmpSrcShift[iOffset*2] );
      d0 = dp0 + dq0;
      d3 = dp3 + dq3;        
      dp = dp0 + dp3;
      dq = dq0 + dq3;
      d  =  d0 + d3;

      #if ENABLE_PCM == 1
      //ToDo: add PCM deblocking
      #endif
      if (d < iBeta)
      { 
        int8_t bFilterP = (dp < iSideThreshold)?1:0;
        int8_t bFilterQ = (dq < iSideThreshold)?1:0;          
        int8_t sw = 0;// xUseStrongFiltering( iOffset, 2*d0, iBeta, iTc, piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+0))
                            //&& xUseStrongFiltering( iOffset, 2*d3, iBeta, iTc, piTmpSrc+iSrcStep*(iIdx*uiPelsInPart+iBlkIdx*4+3));          
        for (i = 0; i < 8/2; i++)
        {
          filter_luma( piTmpSrc+iSrcStep*(iIdx*scu_width+iBlkIdx*4+i), iOffset, iTc, sw, 0, 0, iThrCut, bFilterP, bFilterQ);
        }
      }
    }
  }
  */
}

void filter_deblock_CU(encoder_control* encoder, int32_t xpos, int32_t ypos, int8_t depth, int32_t edge, int8_t dir)
{
  /*
  if(pcCU->getPic()==0||pcCU->getPartitionSize(uiAbsZorderIdx)==SIZE_NONE)
  {
    return;
  }
  TComPic* pcPic     = pcCU->getPic();
  UInt uiCurNumParts = pcPic->getNumPartInCU() >> (uiDepth<<1);
  UInt uiQNumParts   = uiCurNumParts>>2;
  
  if( pcCU->getDepth(uiAbsZorderIdx) > uiDepth )
  {
    for ( UInt uiPartIdx = 0; uiPartIdx < 4; uiPartIdx++, uiAbsZorderIdx+=uiQNumParts )
    {
      UInt uiLPelX   = pcCU->getCUPelX() + g_auiRasterToPelX[ g_auiZscanToRaster[uiAbsZorderIdx] ];
      UInt uiTPelY   = pcCU->getCUPelY() + g_auiRasterToPelY[ g_auiZscanToRaster[uiAbsZorderIdx] ];
      if( ( uiLPelX < pcCU->getSlice()->getSPS()->getPicWidthInLumaSamples() ) && ( uiTPelY < pcCU->getSlice()->getSPS()->getPicHeightInLumaSamples() ) )
      {
        filter_deblock_CU( pcCU, uiAbsZorderIdx, uiDepth+1, Edge );
      }
    }
    return;
  }
  
  xSetLoopfilterParam( pcCU, uiAbsZorderIdx );
  
  xSetEdgefilterTU   ( pcCU, uiAbsZorderIdx , uiAbsZorderIdx, uiDepth );
  xSetEdgefilterPU   ( pcCU, uiAbsZorderIdx );
  
  Int iDir = Edge;
  for( UInt uiPartIdx = uiAbsZorderIdx; uiPartIdx < uiAbsZorderIdx + uiCurNumParts; uiPartIdx++ )
  {
    UInt uiBSCheck;
    if( (g_uiMaxCUWidth >> g_uiMaxCUDepth) == 4 ) 
    {
      uiBSCheck = (iDir == EDGE_VER && uiPartIdx%2 == 0) || (iDir == EDGE_HOR && (uiPartIdx-((uiPartIdx>>2)<<2))/2 == 0);
    }
    else
    {
      uiBSCheck = 1;
    }
    
    if ( m_aapbEdgeFilter[iDir][uiPartIdx] && uiBSCheck )
    {
      xGetBoundaryStrengthSingle ( pcCU, iDir, uiPartIdx );
    }
  }
  
  UInt uiPelsInPart = g_uiMaxCUWidth >> g_uiMaxCUDepth;
  UInt PartIdxIncr = 8 / uiPelsInPart ? 8 / uiPelsInPart : 1 ;
  
  UInt uiSizeInPU = pcPic->getNumPartInWidth()>>(uiDepth);
  
  for ( UInt iEdge = 0; iEdge < uiSizeInPU ; iEdge+=PartIdxIncr)
  {
    xEdgeFilterLuma     ( pcCU, uiAbsZorderIdx, uiDepth, iDir, iEdge );
    if ( (uiPelsInPart>8) || (iEdge % ( (8<<1)/uiPelsInPart ) ) == 0 )
    {
      xEdgeFilterChroma   ( pcCU, uiAbsZorderIdx, uiDepth, iDir, iEdge );
    }
  }
  */
}
void filter_deblock(encoder_control* encoder)
{
  int x,y;
  const int8_t scu_width = (LCU_WIDTH>>(MAX_DEPTH));
  int16_t width = 16;
  /*
  // Horizontal filtering
  for ( UInt uiCUAddr = 0; uiCUAddr < pcPic->getNumCUsInFrame(); uiCUAddr++ )
  {
    TComDataCU* pcCU = pcPic->getCU( uiCUAddr );

    ::memset( m_aapucBS       [EDGE_VER], 0, sizeof( UChar ) * m_uiNumPartitions );
    ::memset( m_aapbEdgeFilter[EDGE_VER], 0, sizeof( Bool  ) * m_uiNumPartitions );

    // CU-based deblocking
    filter_deblock_CU( pcCU, 0, 0, EDGE_VER );
  }

  // Vertical filtering
  for ( UInt uiCUAddr = 0; uiCUAddr < pcPic->getNumCUsInFrame(); uiCUAddr++ )
  {
    TComDataCU* pcCU = pcPic->getCU( uiCUAddr );

    ::memset( m_aapucBS       [EDGE_HOR], 0, sizeof( UChar ) * m_uiNumPartitions );
    ::memset( m_aapbEdgeFilter[EDGE_HOR], 0, sizeof( Bool  ) * m_uiNumPartitions );

    // CU-based deblocking
    filter_deblock_CU( pcCU, 0, 0, EDGE_HOR );
  }
  */
  for(y = width-1; y < encoder->in.height-1; y+=width)
  {
    for(x = width-1; x < encoder->in.width-1; x+=width)
    {
      filter_deblock_edge_luma(encoder, x, y, 2, EDGE_VER, 0);
    }
  }

  for(y = width-1; y < encoder->in.height-1; y+=width)
  {
    for(x = width-1; x < encoder->in.width-1; x+=width)
    {
      filter_deblock_edge_luma(encoder, x, y, 2, EDGE_HOR, 0);
    }
  }
  
}