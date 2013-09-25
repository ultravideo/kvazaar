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

#include "filter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "bitstream.h"
#include "picture.h"
#include "cabac.h"
#include "transform.h"

//////////////////////////////////////////////////////////////////////////
// INITIALIZATIONS
const uint8_t g_tc_table_8x8[54] =
{
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  1,  1,
   1,  1,  1,  1,  1,  1,  1,  2,  2,  2,
   2,  3,  3,  3,  3,  4,  4,  4,  5,  5,
   6,  6,  7,  8,  9, 10, 11, 13, 14, 16,
  18, 20, 22, 24
};

const uint8_t g_beta_table_8x8[52] =
{
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  6,  7,  8,  9,
  10, 11, 12, 13, 14, 15, 16, 17, 18, 20,
  22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
  42, 44, 46, 48, 50, 52, 54, 56, 58, 60,
  62, 64
};

const int16_t g_luma_filter[4][8] =
{
  {  0, 0,   0, 64,  0,   0, 0,  0 },
  { -1, 4, -10, 58, 17,  -5, 1,  0 },
  { -1, 4, -11, 40, 40, -11, 4, -1 },
  {  0, 1,  -5, 17, 58, -10, 4, -1 }
};

const int16_t g_chroma_filter[8][4] =
{
  {  0, 64,  0,  0 },
  { -2, 58, 10, -2 },
  { -4, 54, 16, -2 },
  { -6, 46, 28, -4 },
  { -4, 36, 36, -4 },
  { -4, 28, 46, -6 },
  { -2, 16, 54, -4 },
  { -2, 10, 58, -2 }
};

//////////////////////////////////////////////////////////////////////////
// FUNCTIONS

/**
 * \brief
 */
INLINE void filter_deblock_luma(uint8_t *src, int32_t offset,
                                int32_t tc, int8_t sw,
                                int8_t part_P_nofilter, int8_t part_Q_nofilter,
                                int32_t thr_cut,
                                int8_t filter_second_P, int8_t filter_second_Q)
{
  int32_t delta;
  
  int16_t m0 = src[-offset * 4];
  int16_t m1 = src[-offset * 3];
  int16_t m2 = src[-offset * 2];
  int16_t m3 = src[-offset];
  int16_t m4 = src[0];
  int16_t m5 = src[offset];
  int16_t m6 = src[offset * 2];
  int16_t m7 = src[offset * 3];

  if (sw) {
    src[-offset * 3] = CLIP(m1 - 2*tc, m1 + 2*tc, (2*m0 + 3*m1 +   m2 +   m3 +   m4 + 4) >> 3);
    src[-offset * 2] = CLIP(m2 - 2*tc, m2 + 2*tc, (  m1 +   m2 +   m3 +   m4        + 2) >> 2);
    src[-offset]     = CLIP(m3 - 2*tc, m3 + 2*tc, (  m1 + 2*m2 + 2*m3 + 2*m4 +   m5 + 4) >> 3);
    src[0]           = CLIP(m4 - 2*tc, m4 + 2*tc, (  m2 + 2*m3 + 2*m4 + 2*m5 +   m6 + 4) >> 3);
    src[offset]      = CLIP(m5 - 2*tc, m5 + 2*tc, (  m3 +   m4 +   m5 +   m6        + 2) >> 2);
    src[offset * 2]  = CLIP(m6 - 2*tc, m6 + 2*tc, (  m3 +   m4 +   m5 + 3*m6 + 2*m7 + 4) >> 3);
  } else {
    // Weak filter
    delta = (9*(m4 - m3) - 3*(m5 - m2) + 8) >> 4;

    if (abs(delta) < thr_cut) {
      int32_t tc2  = tc >> 1;
      delta        = CLIP(-tc, tc, delta);
      src[-offset] = CLIP(0, (1 << g_bitdepth) - 1, (m3 + delta));
      src[0]       = CLIP(0, (1 << g_bitdepth) - 1, (m4 - delta));

      if(filter_second_P) {
        int32_t delta1   = CLIP(-tc2, tc2, (((m1 + m3 + 1) >> 1) - m2 + delta) >> 1);
        src[-offset * 2] = CLIP(0, (1 << g_bitdepth) - 1, m2 + delta1);
      }
      if(filter_second_Q) {
        int32_t delta2 = CLIP(-tc2, tc2, (((m6 + m4 + 1) >> 1) - m5 - delta) >> 1);
        src[offset]    = CLIP(0, (1 << g_bitdepth) - 1, m5 + delta2);
      }
    }
  }

  if(part_P_nofilter) {
    src[-offset]   = (uint8_t)m3;
    src[-offset*2] = (uint8_t)m2;
    src[-offset*3] = (uint8_t)m1;
  }
  if(part_Q_nofilter) {
    src[0]        = (uint8_t)m4;
    src[offset]   = (uint8_t)m5;
    src[offset*2] = (uint8_t)m6;
  }
}

/**
 * \brief
 */
INLINE void filter_deblock_chroma(uint8_t *src, int32_t offset, int32_t tc,
                                  int8_t part_P_nofilter, int8_t part_Q_nofilter)
{
  int32_t delta;
  int16_t m2 = src[-offset * 2];
  int16_t m3 = src[-offset];
  int16_t m4 = src[0];
  int16_t m5 = src[offset];  
  
  delta = CLIP(-tc,tc, (((m4 - m3) << 2) + m2 - m5 + 4 ) >> 3);
  if(!part_P_nofilter) {
    src[-offset] = CLIP(0, (1 << g_bitdepth) - 1, m3 + delta);
  }
  if(!part_Q_nofilter) {
    src[0] = CLIP(0, (1 << g_bitdepth) - 1, m4 - delta);
  }
}

/**
 * \brief
 */
void filter_deblock_edge_luma(encoder_control *encoder, 
                              int32_t xpos, int32_t ypos, 
                              int8_t depth, int8_t dir)
{
  int32_t stride = encoder->in.cur_pic->width;
  int32_t offset = stride;
  int32_t beta_offset_div2 = encoder->beta_offset_div2;
  int32_t tc_offset_div2   = encoder->tc_offset_div2;
  int8_t strength = 2; // Filter strength
  // TODO: support 10+bits
  uint8_t *orig_src = &encoder->in.cur_pic->y_recdata[xpos + ypos*stride];
  uint8_t *src = orig_src;
  int32_t step = 1;
  
  if(dir == EDGE_VER) {
    offset = 1;
    step = stride;
  }
  
  {
    int32_t qp              = encoder->QP;
    int32_t bitdepth_scale  = 1 << (g_bitdepth - 8);
    int32_t tc_index        = CLIP(0, 51 + 2, (int32_t)(qp + 2*(strength - 1) + (tc_offset_div2 << 1)));
    int32_t b_index         = CLIP(0, 51, qp + (beta_offset_div2 << 1));
    int32_t tc              = g_tc_table_8x8[tc_index] * bitdepth_scale;
    int32_t beta            = g_beta_table_8x8[b_index] * bitdepth_scale;
    int32_t side_threshold  = (beta + (beta >>1 )) >> 3;
    int32_t thr_cut         = tc * 10;
    uint32_t blocks_in_part = (LCU_WIDTH >> depth) / 4;
    uint32_t block_idx;

    // TODO: add CU based QP calculation 

    // For each 4-pixel part in the edge
    for (block_idx = 0; block_idx < blocks_in_part; ++block_idx)
    {
      int32_t dp0, dq0, dp3, dq3, d0, d3, dp, dq, d;

      // Check conditions for filtering
      // TODO: Get rid of these inline defines.
      #define calc_DP(s,o) abs( (int16_t)s[-o*3] - (int16_t)2*s[-o*2] + (int16_t)s[-o] )
      #define calc_DQ(s,o) abs( (int16_t)s[0]    - (int16_t)2*s[o]    + (int16_t)s[o*2] )

      dp0 = calc_DP((src+step*(block_idx*4+0)), offset);      
      dq0 = calc_DQ((src+step*(block_idx*4+0)), offset);
      dp3 = calc_DP((src+step*(block_idx*4+3)), offset);      
      dq3 = calc_DQ((src+step*(block_idx*4+3)), offset);
      d0 = dp0 + dq0;
      d3 = dp3 + dq3;
      dp = dp0 + dp3;
      dq = dq0 + dq3;
      d  =  d0 + d3;

      #if ENABLE_PCM
      // TODO: add PCM deblocking
      #endif

      if (d < beta) { 
        int8_t filter_P = (dp < side_threshold) ? 1 : 0;
        int8_t filter_Q = (dq < side_threshold) ? 1 : 0;

        // Strong filtering flag checking
        #define useStrongFiltering(o,d,s) ( ((abs(s[-o*4]-s[-o]) + abs(s[o*3]-s[0])) < (beta>>3)) && (d<(beta>>2)) && ( abs(s[-o]-s[0]) < ((tc*5+1)>>1)) )      
        int8_t sw = useStrongFiltering(offset, 2*d0, (src+step*(block_idx*4+0))) &&
                    useStrongFiltering(offset, 2*d3, (src+step*(block_idx*4+3)));

        // Filter four rows/columns
        filter_deblock_luma(src + step * (4*block_idx + 0), offset, tc, sw, 0, 0, thr_cut, filter_P, filter_Q);
        filter_deblock_luma(src + step * (4*block_idx + 1), offset, tc, sw, 0, 0, thr_cut, filter_P, filter_Q);
        filter_deblock_luma(src + step * (4*block_idx + 2), offset, tc, sw, 0, 0, thr_cut, filter_P, filter_Q);
        filter_deblock_luma(src + step * (4*block_idx + 3), offset, tc, sw, 0, 0, thr_cut, filter_P, filter_Q);
      }
    }
  }
}

/**
 * \brief
 */
void filter_deblock_edge_chroma(encoder_control *encoder, 
                                int32_t x, int32_t y, 
                                int8_t depth, int8_t dir)
{
  int32_t stride = encoder->in.cur_pic->width >> 1;
  int32_t tc_offset_div2 = encoder->tc_offset_div2;
  int8_t num_parts = 1;
  // TODO: support 10+bits
  uint8_t* src_u = &encoder->in.cur_pic->u_recdata[x + y*stride];
  uint8_t* src_v = &encoder->in.cur_pic->v_recdata[x + y*stride];
  // Init offset and step to EDGE_HOR
  int32_t offset = stride;
  int32_t step = 1;

  // We cannot filter edges not on 8x8 grid
  if(depth == MAX_DEPTH && (( (y & 0x7) && dir == EDGE_HOR ) || ( (x & 0x7) && dir == EDGE_VER ) ) )
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
    int32_t QP             = g_chroma_scale[encoder->QP];
    int32_t bitdepth_scale = 1 << (g_bitdepth-8);
    int32_t TC_index       = CLIP(0, 51+2, (int32_t)(QP + 2 + (tc_offset_div2 << 1)));    
    int32_t Tc             = g_tc_table_8x8[TC_index]*bitdepth_scale;
    uint32_t blocks_in_part= (LCU_WIDTH>>(depth+1)) / 4;
    uint32_t blk_idx;

    for (blk_idx = 0; blk_idx < blocks_in_part; ++blk_idx)
    {
      // Chroma U
      filter_deblock_chroma(src_u + step * (4*blk_idx + 0), offset, Tc, 0, 0);
      filter_deblock_chroma(src_u + step * (4*blk_idx + 1), offset, Tc, 0, 0);
      filter_deblock_chroma(src_u + step * (4*blk_idx + 2), offset, Tc, 0, 0);
      filter_deblock_chroma(src_u + step * (4*blk_idx + 3), offset, Tc, 0, 0);
      // Chroma V
      filter_deblock_chroma(src_v + step * (4*blk_idx + 0), offset, Tc, 0, 0);
      filter_deblock_chroma(src_v + step * (4*blk_idx + 1), offset, Tc, 0, 0);
      filter_deblock_chroma(src_v + step * (4*blk_idx + 2), offset, Tc, 0, 0);
      filter_deblock_chroma(src_v + step * (4*blk_idx + 3), offset, Tc, 0, 0);
    }
  }
}

/** 
 * \brief function to split LCU into smaller CU blocks
 * \param encoder the encoder info structure
 * \param xCtb block x-position (as SCU)
 * \param yCtb block y-position (as SCU)
 * \param depth block depth
 * \param edge which edge we are filtering
 *
 * This function takes (SCU) block position as input and splits the block
 * until the coded block size has been achived. Calls luma and chroma filtering
 * functions for each coded CU size.
 */
void filter_deblock_cu(encoder_control *encoder, int32_t x, int32_t y, int8_t depth, int32_t edge)
{
  cu_info *cur_cu = &encoder->in.cur_pic->cu_array[depth][x + y*(encoder->in.width_in_lcu << MAX_DEPTH)];
  uint8_t split_flag = (cur_cu->depth > depth) ? 1 : 0;
  uint8_t border_x = (encoder->in.width  < x*(LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth)) ? 1 : 0;
  uint8_t border_y = (encoder->in.height < y*(LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth)) ? 1 : 0;
  uint8_t border_split_x = (encoder->in.width  < ((x + 1) * (LCU_WIDTH >> MAX_DEPTH)) + (LCU_WIDTH >> (depth + 1))) ? 0 : 1;
  uint8_t border_split_y = (encoder->in.height < ((y + 1) * (LCU_WIDTH >> MAX_DEPTH)) + (LCU_WIDTH >> (depth + 1))) ? 0 : 1;

  uint8_t border = border_x | border_y; // are we in any border CU?

  // split 64x64, on split flag and on border
  if (depth == 0 || split_flag || border) {
    // Split the four sub-blocks of this block recursively.
    uint8_t change = 1 << (MAX_DEPTH - 1 - depth);
    filter_deblock_cu(encoder, x, y, depth + 1, edge);
    if(!border_x || border_split_x) {
      filter_deblock_cu(encoder, x + change, y, depth + 1, edge);
    }
    if(!border_y || border_split_y) {
      filter_deblock_cu(encoder, x , y + change, depth + 1, edge);
    }
    if((!border_x && !border_y) || (border_split_x && border_split_y)) {
      filter_deblock_cu(encoder, x + change, y + change, depth + 1, edge);
    }
    return;
  }

  // no filtering on borders (where filter would use pixels outside the picture)
  if ((x == 0 && edge == EDGE_VER) || (y == 0 && edge == EDGE_HOR)) return;

  // do the filtering for block edge
  filter_deblock_edge_luma(encoder,   x*(LCU_WIDTH >> MAX_DEPTH),       y*(LCU_WIDTH >> MAX_DEPTH),       depth, edge);
  filter_deblock_edge_chroma(encoder, x*(LCU_WIDTH >> (MAX_DEPTH + 1)), y*(LCU_WIDTH >> (MAX_DEPTH + 1)), depth, edge);
}

/** 
 * \brief Main function for Deblocking filtering
 * \param encoder the encoder info structure
 *
 * This is the main function for deblocking filter, it will loop through all
 * the Largest Coding Units (LCU) and call filter_deblock_cu with absolute
 * X and Y coordinates of the LCU.
 */
void filter_deblock(encoder_control* encoder)
{
  int16_t x, y;

  // TODO: Optimization: add thread for each LCU
  // Filter vertically.
  for (y = 0; y < encoder->in.height_in_lcu; y++)
  {
    for (x = 0; x < encoder->in.width_in_lcu; x++)
    {
      filter_deblock_cu(encoder, x << MAX_DEPTH, y << MAX_DEPTH, 0, EDGE_VER);
    }
  }

  // Filter horizontally.
  for (y = 0; y < encoder->in.height_in_lcu; y++)
  {
    for (x = 0; x < encoder->in.width_in_lcu; x++)
    {
      filter_deblock_cu(encoder, x << MAX_DEPTH, y << MAX_DEPTH, 0, EDGE_HOR);
    }
  }
}


/** 
 * \brief Interpolation for chroma half-pixel
 * \param src source image in integer pels (-2..width+3, -2..height+3)
 * \param src_stride stride of source image
 * \param width width of source image block
 * \param height height of source image block
 * \param dst destination image in half-pixel resolution
 * \param dst_stride stride of destination image
 *
 */
void filter_inter_halfpel_chroma(int16_t *src, int16_t src_stride, int width, int height, int16_t *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag)
{
  /* ____________
   * | B0,0|ae0,0|
   * |ea0,0|ee0,0|
   * 
   * ae0,0 = (-4*B-1,0  + 36*B0,0  + 36*B1,0  - 4*B2,0)  >> shift1
   * ea0,0 = (-4*B0,-1  + 36*B0,0  + 36*B0,1  - 4*B0,2)  >> shift1
   * ee0,0 = (-4*ae0,-1 + 36*ae0,0 + 36*ae0,1 - 4*ae0,2) >> shift2
   */

  int32_t x, y;
  int32_t shift1 = g_bitdepth-8;
  int32_t shift2 = 6;
  int32_t shift3 = 14-g_bitdepth;
  int32_t offset = 1<<(shift2-1); //!< offset for rounding purposes

  // Loop source pixels and generate four filtered half-pel pixels on each round
  for (y = 0; y < height; y++) {
    int dst_pos_y = (y<<1)*dst_stride;
    int src_pos_y = y*src_stride;
    for (x = 0; x < width; x++) {

      // Calculate current dst and src pixel positions
      int dst_pos = dst_pos_y+(x<<1);
      int src_pos = src_pos_y+x;

      // Temporary variables..
      int32_t ae_temp, ae_temp1, ae_temp2, ae_temp3;

      // Original pixel (not really needed)
      dst[dst_pos] = src[src_pos]; //B0,0

      // ae0,0 - We need this only when hor_flag and for ee0,0
      if (hor_flag) {
        ae_temp = ((-4*src[src_pos - 1] + 36*src[src_pos] + 36*src[src_pos + 1] - 4*src[src_pos + 2] ) >> shift1) + offset; // ae0,0
      }
      // ea0,0 - needed only when ver_flag
      if(ver_flag) {
        dst[dst_pos + 1*dst_stride] = (((-4*src[src_pos - src_stride] + 36*src[src_pos] + 36*src[src_pos + src_stride] 
                                        - 4*src[src_pos + 2*src_stride]  ) >> shift1) + offset) >> shift3; // ea0,0
      }

      // When both flags, we use _only_ this pixel (but still need ae0,0 for it)
      if (hor_flag && ver_flag) {
        // Calculate temporary values..
        //TODO: optimization, store these values
        src_pos -= src_stride;  //0,-1
        ae_temp1 = ((-4*src[src_pos - 1] + 36*src[src_pos] + 36*src[src_pos + 1] - 4*src[src_pos + 2] ) >> shift1) + offset; // ae0,-1
        src_pos += src_stride;  //0,1
        ae_temp2 = ((-4*src[src_pos - 1] + 36*src[src_pos] + 36*src[src_pos + 1] - 4*src[src_pos + 2] ) >> shift1) + offset; // ae0,1
        src_pos += src_stride;  //0,2
        ae_temp3 = ((-4*src[src_pos - 1] + 36*src[src_pos] + 36*src[src_pos + 1] - 4*src[src_pos + 2] ) >> shift1) + offset; // ae0,2

        dst[dst_pos + 1*dst_stride + 1] = ((( -4*ae_temp1 + 36*dst[dst_pos + 1] + 36*ae_temp2 - 4*ae_temp3 ) >> shift2) + offset) >> shift3; // ee0,0
      }

      if(hor_flag) {
        dst[dst_pos + 1] = ae_temp >> shift3;
      }
    }
  }
}
