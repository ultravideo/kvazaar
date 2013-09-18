/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file filter.c
    \brief filtering
    \author Marko Viitanen
    \date 2013-05
    
    Filtering functions
*/

#include "filter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "bitstream.h"
#include "picture.h"
#include "cabac.h"


extern const uint8_t g_aucChromaScale[58];
const uint8_t tctable_8x8[54] =
{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,5,5,6,6,7,8,9,10,11,13,14,16,18,20,22,24
};

const uint8_t betatable_8x8[52] =
{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,7,8,9,10,11,12,13,14,15,16,17,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64
};


INLINE void filter_deblock_luma( uint8_t* src, int32_t offset, int32_t tc , int8_t sw, int8_t part_P_nofilter, int8_t part_Q_nofilter, int32_t thr_cut, int8_t filter_second_P, int8_t filter_second_Q)
{
  int32_t delta;
  
  int16_t m0  = src[-offset*4];
  int16_t m1  = src[-offset*3];
  int16_t m2  = src[-offset*2];
  int16_t m3  = src[-offset];
  int16_t m4  = src[0];  
  int16_t m5  = src[ offset];  
  int16_t m6  = src[ offset*2];  
  int16_t m7  = src[ offset*3];  

  if (sw)
  {
    src[-offset*3] = CLIP(m1-2*tc, m1+2*tc, ((2*m0 + 3*m1 + m2 + m3 + m4 + 4 )>>3));
    src[-offset*2] = CLIP(m2-2*tc, m2+2*tc, ((m1 + m2 + m3 + m4 + 2)>>2));
    src[-offset]   = CLIP(m3-2*tc, m3+2*tc, ((m1 + 2*m2 + 2*m3 + 2*m4 + m5 + 4) >> 3));
    src[0]         = CLIP(m4-2*tc, m4+2*tc, ((m2 + 2*m3 + 2*m4 + 2*m5 + m6 + 4) >> 3));
    src[ offset]   = CLIP(m5-2*tc, m5+2*tc, ((m3 + m4 + m5 + m6 + 2)>>2));
    src[ offset*2] = CLIP(m6-2*tc, m6+2*tc, ((m3 + m4 + m5 + 3*m6 + 2*m7 + 4 )>>3));
  }
  else
  {
    // Weak filter
    delta = (9*(m4-m3) -3*(m5-m2) + 8)>>4;

    if ( abs(delta) < thr_cut )
    {
      int32_t tc2 = tc>>1;
      delta = CLIP(-tc, tc, delta);        
      src[-offset] = CLIP(0,(1 << g_bitdepth)-1,(m3+delta));
      src[0] = CLIP(0,(1 << g_bitdepth)-1,(m4-delta));

      
      if(filter_second_P)
      {
        int32_t delta1 = CLIP(-tc2, tc2, (( ((m1+m3+1)>>1)- m2+delta)>>1));
        src[-offset*2] = CLIP(0,(1 << g_bitdepth)-1,(m2+delta1));
      }
      if(filter_second_Q)
      {
        int32_t delta2 = CLIP(-tc2, tc2, (( ((m6+m4+1)>>1)- m5-delta)>>1));
        src[ offset] = CLIP(0,(1 << g_bitdepth)-1,(m5+delta2));
      }
    }
  }

  if(part_P_nofilter)
  {
    src[-offset]   = (uint8_t)m3;
    src[-offset*2] = (uint8_t)m2;
    src[-offset*3] = (uint8_t)m1;
  }
  if(part_Q_nofilter)
  {
    src[0]          = (uint8_t)m4;
    src[ offset]   = (uint8_t)m5;
    src[ offset*2] = (uint8_t)m6;
  }
}

INLINE void filter_deblock_chroma( uint8_t* src, int32_t offset, int32_t tc ,int8_t part_P_nofilter, int8_t part_Q_nofilter)
{
  int32_t delta;
  
  int16_t m2  = src[-offset*2];
  int16_t m3  = src[-offset];
  int16_t m4  = src[0];
  int16_t m5  = src[ offset];  
  
  delta = CLIP(-tc,tc, (((( m4 - m3 ) << 2 ) + m2 - m5 + 4 ) >> 3) );
  if(!part_P_nofilter)
  {
    src[-offset] = CLIP(0,(1 << g_bitdepth)-1,m3+delta);
  }
  if(!part_Q_nofilter)
  {
    src[0] = CLIP(0,(1 << g_bitdepth)-1,m4-delta);
  }
}

void filter_deblock_edge_luma(encoder_control* encoder, int32_t xpos, int32_t ypos, int8_t depth, int8_t dir)
{
  int32_t stride = encoder->in.cur_pic->width;
  int32_t offset = stride;
  int32_t betaOffsetDiv2 = encoder->beta_offset_div2;
  int32_t tcOffsetDiv2   = encoder->tc_offset_div2;
  int8_t uiBs       = 2; /* Filter strength */
  /* TODO: support 10+bits */
  uint8_t* origsrc      = &encoder->in.cur_pic->yRecData[xpos+ypos*stride];
  uint8_t* src = origsrc;
  int32_t step = 1;
  //CU_info* cu = &encoder->in.cur_pic->CU[depth][(xpos>>scu_width_log2) + (ypos>>scu_width_log2)*(encoder->in.width>>scu_width_log2)];
  
  if(dir == EDGE_VER)
  {
    offset = 1;
    step = stride;
  }
  
  {
    int32_t QP             = encoder->QP;
    int32_t bitdepth_scale = 1 << (g_bitdepth-8);
    int32_t TC_index       = CLIP(0, 51+2, (int32_t)(QP + 2*(uiBs-1) + (tcOffsetDiv2 << 1)));
    int32_t B_index        = CLIP(0, 51, QP + (betaOffsetDiv2 << 1));
    int32_t Tc             = tctable_8x8[TC_index]*bitdepth_scale;
    int32_t Beta           = betatable_8x8[B_index]*bitdepth_scale;
    int32_t side_threshold = (Beta+(Beta>>1))>>3;
    int32_t iThrCut        = Tc*10;
    uint32_t blocks_in_part= (LCU_WIDTH>>depth) / 4;
    uint32_t block_idx;

    /* TODO: add CU based QP calculation */

    /* For each 4-pixel part in the edge */
    for (block_idx = 0; block_idx < blocks_in_part; block_idx++)
    {
      int32_t dp0,dq0,dp3,dq3,d0,d3,dp,dq,d;

      /* Check conditions for filtering */
      #define calc_DP(s,o) abs( (int16_t)s[-o*3] - (int16_t)2*s[-o*2] + (int16_t)s[-o] )
      #define calc_DQ(s,o) abs( (int16_t)s[0]    - (int16_t)2*s[o]    + (int16_t)s[o*2] )

      dp0 = calc_DP( (src+step*(block_idx*4+0)), offset);      
      dq0 = calc_DQ( (src+step*(block_idx*4+0)), offset);
      dp3 = calc_DP( (src+step*(block_idx*4+3)), offset);      
      dq3 = calc_DQ( (src+step*(block_idx*4+3)), offset);
      d0 = dp0 + dq0;
      d3 = dp3 + dq3;
      dp = dp0 + dp3;
      dq = dq0 + dq3;
      d  =  d0 + d3;

      #if ENABLE_PCM == 1
      //TODO: add PCM deblocking
      #endif
      if (d < Beta)
      { 
        int8_t filter_P = (dp < side_threshold)?1:0;
        int8_t filter_Q = (dq < side_threshold)?1:0;

        /* Strong filtering flag checking */
        #define useStrongFiltering(o,d,s) ( ((abs(s[-o*4]-s[-o]) + abs(s[o*3]-s[0])) < (Beta>>3)) && (d<(Beta>>2)) && ( abs(s[-o]-s[0]) < ((Tc*5+1)>>1)) )      
        int8_t sw = useStrongFiltering( offset, 2*d0, (src+step*(block_idx*4+0))) &&
                    useStrongFiltering( offset, 2*d3, (src+step*(block_idx*4+3)));

        /* Filter four rows/columns */
        filter_deblock_luma( src+step*(block_idx*4+0), offset, Tc, sw, 0, 0, iThrCut, filter_P, filter_Q);
        filter_deblock_luma( src+step*(block_idx*4+1), offset, Tc, sw, 0, 0, iThrCut, filter_P, filter_Q);
        filter_deblock_luma( src+step*(block_idx*4+2), offset, Tc, sw, 0, 0, iThrCut, filter_P, filter_Q);
        filter_deblock_luma( src+step*(block_idx*4+3), offset, Tc, sw, 0, 0, iThrCut, filter_P, filter_Q);
      }
    }
  }
}

void filter_deblock_edge_chroma(encoder_control* encoder,int32_t xpos, int32_t ypos, int8_t depth, int8_t dir)
{
  int32_t stride = encoder->in.cur_pic->width>>1;
  int32_t tcOffsetDiv2   = encoder->tc_offset_div2;
  int8_t uiNumParts = 1;
  /* TODO: support 10+bits */
  uint8_t* srcU      = &encoder->in.cur_pic->uRecData[xpos+ypos*stride];
  uint8_t* srcV      = &encoder->in.cur_pic->vRecData[xpos+ypos*stride];
  /* Init offset and step to EDGE_HOR */
  int32_t offset = stride;
  int32_t step = 1;

  /* We cannot filter edges not on 8x8 grid */
  if( depth == MAX_DEPTH && (( (ypos & 0x7) && dir == EDGE_HOR ) || ( (xpos & 0x7) && dir == EDGE_VER ) ) )
  {
    return;
  }

  if(dir == EDGE_VER)
  {
    offset = 1;
    step = stride;
  }

  // For each subpart
  {
    int32_t QP             = g_aucChromaScale[encoder->QP];
    int32_t bitdepth_scale = 1 << (g_bitdepth-8);
    int32_t TC_index       = CLIP(0, 51+2, (int32_t)(QP + 2 + (tcOffsetDiv2 << 1)));    
    int32_t Tc             = tctable_8x8[TC_index]*bitdepth_scale;
    uint32_t blocks_in_part= (LCU_WIDTH>>(depth+1)) / 4;
    uint32_t blk_idx;

    for (blk_idx = 0; blk_idx < blocks_in_part; blk_idx++)
    {
      /* Chroma U */
      filter_deblock_chroma( srcU+step*(blk_idx*4+0), offset, Tc,0, 0);
      filter_deblock_chroma( srcU+step*(blk_idx*4+1), offset, Tc,0, 0);
      filter_deblock_chroma( srcU+step*(blk_idx*4+2), offset, Tc,0, 0);
      filter_deblock_chroma( srcU+step*(blk_idx*4+3), offset, Tc,0, 0);
      /* Chroma V */
      filter_deblock_chroma( srcV+step*(blk_idx*4+0), offset, Tc,0, 0);
      filter_deblock_chroma( srcV+step*(blk_idx*4+1), offset, Tc,0, 0);
      filter_deblock_chroma( srcV+step*(blk_idx*4+2), offset, Tc,0, 0);
      filter_deblock_chroma( srcV+step*(blk_idx*4+3), offset, Tc,0, 0);
    }
  }
}

/*! \brief function to split LCU into smaller CU blocks
  \param encoder the encoder info structure
  \param xCtb block x-position (as SCU)
  \param yCtb block y-position (as SCU)
  \param depth block depth
  \param edge which edge we are filtering

  This function takes (SCU) block position as input and splits the block
  until the coded block size has been achived. Calls luma and chroma filtering
  functions for each coded CU size
*/
void filter_deblock_CU(encoder_control* encoder, int32_t xCtb, int32_t yCtb, int8_t depth, int32_t edge)
{
  CU_info *cur_CU = &encoder->in.cur_pic->CU[depth][xCtb+yCtb*(encoder->in.width_in_lcu<<MAX_DEPTH)];
  uint8_t split_flag = (cur_CU->depth > depth)?1:0;
  uint8_t border_x = ((encoder->in.width)<( xCtb*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>depth) ))?1:0;
  uint8_t border_y = ((encoder->in.height)<( yCtb*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>depth) ))?1:0;
  uint8_t border_split_x = ((encoder->in.width)  < ( (xCtb+1)*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>(depth+1)) ))?0:1;
  uint8_t border_split_y = ((encoder->in.height) < ( (yCtb+1)*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>(depth+1)) ))?0:1;

  uint8_t border = border_x | border_y; /*!< are we in any border CU */

  /* split 64x64, on split flag and on border */
  if(!depth || split_flag || border)
  {
    /* Split blocks and remember to change x and y block positions */
    uint8_t change = 1<<(MAX_DEPTH-1-depth);
    filter_deblock_CU(encoder,xCtb,yCtb,depth+1,edge); /* x,y */

    if(!border_x || border_split_x)
    {
      filter_deblock_CU(encoder,xCtb+change,yCtb,depth+1,edge); /* x+1,y */
    }
    if(!border_y || border_split_y)
    {
      filter_deblock_CU(encoder,xCtb,yCtb+change,depth+1,edge); /* x,y+1 */
    }
    if((!border_x && !border_y) || (border_split_x && border_split_y) )
    {
      filter_deblock_CU(encoder,xCtb+change,yCtb+change,depth+1,edge); /* x+1,y+1 */
    }
    return;
  }
  /* no filtering on borders (where filter would use pixels outside the picture) */
  if((!xCtb && edge == EDGE_VER) || (!yCtb && edge == EDGE_HOR)) return;

  /* do the filtering for block edge */
  filter_deblock_edge_luma(encoder, xCtb*(LCU_WIDTH>>MAX_DEPTH), yCtb*(LCU_WIDTH>>MAX_DEPTH), depth, edge);
  filter_deblock_edge_chroma(encoder, xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)), yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)), depth, edge);
}



/*! \brief Main function for Deblocking filtering
  \param encoder the encoder info structure

  This is the main function for deblocking filter, it will loop through all the
  Largest Coding Units and call filter_deblock_CU with absolute X and Y coordinate
  of the LCU.
*/
void filter_deblock(encoder_control* encoder)
{
  int16_t xCtb,yCtb;

  /* TODO: Optimization: add thread for each LCU */
  /* Loop through every LCU in the slice */
  for(yCtb = 0; yCtb < encoder->in.height_in_lcu; yCtb++)
  {
    for(xCtb = 0; xCtb < encoder->in.width_in_lcu; xCtb++)
    {
      filter_deblock_CU(encoder, xCtb<<MAX_DEPTH, yCtb<<MAX_DEPTH, 0, EDGE_VER);
    }
  }

  /* Loop through every LCU in the slice */
  for(yCtb = 0; yCtb < encoder->in.height_in_lcu; yCtb++)
  {
    for(xCtb = 0; xCtb < encoder->in.width_in_lcu; xCtb++)
    {
      filter_deblock_CU(encoder, xCtb<<MAX_DEPTH, yCtb<<MAX_DEPTH, 0, EDGE_HOR);
    }
  }
  
}