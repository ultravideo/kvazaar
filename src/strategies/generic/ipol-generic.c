/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

/*
 * \file
 */

#include <stdlib.h>

#include "ipol-generic.h"
#include "strategyselector.h"
#include "encoder.h"
#include "picture-generic.h"

extern int8_t g_luma_filter[4][8];
extern int8_t g_chroma_filter[8][4];

int16_t eight_tap_filter_hor_generic(int8_t *filter, kvz_pixel *data)
{
  int16_t temp = 0;
  for (int i = 0; i < 8; ++i)
  {
    temp += filter[i] * data[i];
  }

  return temp;
}

int32_t eight_tap_filter_hor_16bit_generic(int8_t *filter, int16_t *data)
{
  int32_t temp = 0;
  for (int i = 0; i < 8; ++i)
  {
    temp += filter[i] * data[i];
  }

  return temp;
}

int16_t eight_tap_filter_ver_generic(int8_t *filter, kvz_pixel *data, int16_t stride)
{
  int16_t temp = 0;
  for (int i = 0; i < 8; ++i)
  {
    temp += filter[i] * data[stride * i];
  }

  return temp;
}

int32_t eight_tap_filter_ver_16bit_generic(int8_t *filter, int16_t *data, int16_t stride)
{
  int32_t temp = 0;
  for (int i = 0; i < 8; ++i)
  {
    temp += filter[i] * data[stride * i];
  }

  return temp;
}

int16_t four_tap_filter_hor_generic(int8_t *filter, kvz_pixel *data)
{
  int16_t temp = 0;
  for (int i = 0; i < 4; ++i)
  {
    temp += filter[i] * data[i];
  }

  return temp;
}

int32_t four_tap_filter_hor_16bit_generic(int8_t *filter, int16_t *data)
{
  int32_t temp = 0;
  for (int i = 0; i < 4; ++i)
  {
    temp += filter[i] * data[i];
  }

  return temp;
}

int16_t four_tap_filter_ver_generic(int8_t *filter, kvz_pixel *data, int16_t stride)
{
  int16_t temp = 0;
  for (int i = 0; i < 4; ++i)
  {
    temp += filter[i] * data[stride * i];
  }

  return temp;
}

int32_t four_tap_filter_ver_16bit_generic(int8_t *filter, int16_t *data, int16_t stride)
{
  int32_t temp = 0;
  for (int i = 0; i < 4; ++i)
  {
    temp += filter[i] * data[stride * i];
  }

  return temp;
}

void filter_inter_quarterpel_luma_generic(const encoder_control_t * const encoder, kvz_pixel *src, int16_t src_stride, int width, int height, kvz_pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag)
{
  //TODO: horizontal and vertical only filtering
  int32_t x, y;
  int16_t shift1 = KVZ_BIT_DEPTH - 8;
  int32_t shift2 = 6;
  int32_t shift3 = 14 - KVZ_BIT_DEPTH;
  int32_t offset23 = 1 << (shift2 + shift3 - 1);

  //coefficients for 1/4, 2/4 and 3/4 positions
  int8_t *c0, *c1, *c2, *c3;

  c0 = g_luma_filter[0];
  c1 = g_luma_filter[1];
  c2 = g_luma_filter[2];
  c3 = g_luma_filter[3];

  #define FILTER_OFFSET 3
  #define FILTER_SIZE 8

  int16_t flipped_hor_filtered[4 * (LCU_WIDTH + 1) + FILTER_SIZE][(LCU_WIDTH + 1) + FILTER_SIZE];

  // Filter horizontally and flip x and y
  for (x = 0; x < width; ++x) {
    for (y = 0; y < height + FILTER_SIZE; ++y) {
      int ypos = y - FILTER_OFFSET;
      int xpos = x - FILTER_OFFSET;
      // Original pixel
      flipped_hor_filtered[4 * x + 0][y] = (c0[FILTER_OFFSET] * src[src_stride*ypos + xpos + FILTER_OFFSET]) << shift1;
      flipped_hor_filtered[4 * x + 1][y] = eight_tap_filter_hor_generic(c1, &src[src_stride*ypos + xpos]) << shift1;
      flipped_hor_filtered[4 * x + 2][y] = eight_tap_filter_hor_generic(c2, &src[src_stride*ypos + xpos]) << shift1;
      flipped_hor_filtered[4 * x + 3][y] = eight_tap_filter_hor_generic(c3, &src[src_stride*ypos + xpos]) << shift1;

    }
  }

  // Filter vertically and flip x and y
  for (x = 0; x < 4 * width; ++x) {
    for (y = 0; y < height; ++y) {
      int ypos = y;
      int xpos = x;
      dst[(4 * y + 0)*dst_stride + x] = fast_clip_32bit_to_pixel(((c0[FILTER_OFFSET] * flipped_hor_filtered[xpos][ypos + FILTER_OFFSET] + offset23) >> shift2) >> shift3); 
      dst[(4 * y + 1)*dst_stride + x] = fast_clip_32bit_to_pixel(((eight_tap_filter_hor_16bit_generic(c1, &flipped_hor_filtered[xpos][ypos]) + offset23) >> shift2) >> shift3);
      dst[(4 * y + 2)*dst_stride + x] = fast_clip_32bit_to_pixel(((eight_tap_filter_hor_16bit_generic(c2, &flipped_hor_filtered[xpos][ypos]) + offset23) >> shift2) >> shift3);
      dst[(4 * y + 3)*dst_stride + x] = fast_clip_32bit_to_pixel(((eight_tap_filter_hor_16bit_generic(c3, &flipped_hor_filtered[xpos][ypos]) + offset23) >> shift2) >> shift3);

    }
  }
}

void sample_quarterpel_luma_generic(const encoder_control_t * const encoder, kvz_pixel *src, int16_t src_stride, int width, int height, kvz_pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag, const int16_t mv[2])
{
  //TODO: horizontal and vertical only filtering
  int32_t x, y;
  int16_t shift1 = KVZ_BIT_DEPTH - 8;
  int32_t shift2 = 6;
  int32_t shift3 = 14 - KVZ_BIT_DEPTH;
  int32_t offset23 = 1 << (shift2 + shift3 - 1);

  //coefficients for 1/4, 2/4 and 3/4 positions
  int8_t *hor_filter = g_luma_filter[mv[0]&3];
  int8_t *ver_filter = g_luma_filter[mv[1]&3];

  int16_t flipped_hor_filtered[(LCU_WIDTH + 1) + FILTER_SIZE][(LCU_WIDTH + 1) + FILTER_SIZE];

  // Filter horizontally and flip x and y
  for (x = 0; x < width; ++x) {
    for (y = 0; y < height + FILTER_SIZE; ++y) {
      int ypos = y - FILTER_OFFSET;
      int xpos = x - FILTER_OFFSET;
      flipped_hor_filtered[x][y] = eight_tap_filter_hor_generic(hor_filter, &src[src_stride*ypos + xpos]) >> shift1;
    }
  }

  // Filter vertically and flip x and y
  for (x = 0; x < width; ++x) {
    for (y = 0; y < height; ++y) {
      int ypos = y;
      int xpos = x;
      dst[y*dst_stride + x] = fast_clip_32bit_to_pixel(((eight_tap_filter_hor_16bit_generic(ver_filter, &flipped_hor_filtered[xpos][ypos]) + offset23) >> shift2) >> shift3);
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
void filter_inter_halfpel_chroma_generic(const encoder_control_t * const encoder, kvz_pixel *src, int16_t src_stride, int width, int height, kvz_pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag)
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
  int32_t shift1 = KVZ_BIT_DEPTH - 8;
  int32_t shift2 = 6;
  int32_t shift3 = 14 - KVZ_BIT_DEPTH;
  int32_t offset3 = 1 << (shift3 - 1);
  int32_t offset23 = 1 << (shift2 + shift3 - 1);

  int8_t* c = g_chroma_filter[4];
  int16_t temp[4] = {0,0,0,0};

  // Loop source pixels and generate four filtered half-pel pixels on each round
  for (y = 0; y < height; y++) {
    int dst_pos_y = (y << 1)*dst_stride;
    int src_pos_y = y*src_stride;
    for (x = 0; x < width; x++) {
      // Calculate current dst and src pixel positions
      int dst_pos = dst_pos_y + (x << 1);
      int src_pos = src_pos_y + x;

      // Original pixel (not really needed)
      dst[dst_pos] = src[src_pos]; //B0,0

      // ae0,0 - We need this only when hor_flag and for ee0,0
      if (hor_flag) {
        temp[1] = four_tap_filter_hor_generic(c, &src[src_pos - 1]) >> shift1; // ae0,0
      }
      // ea0,0 - needed only when ver_flag
      if (ver_flag) {
        dst[dst_pos + 1 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_generic(c, &src[src_pos - src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3); // ea0,0
      }

      // When both flags, we use _only_ this pixel (but still need ae0,0 for it)
      if (hor_flag && ver_flag) {
        // Calculate temporary values..
        src_pos -= src_stride;  //0,-1
        temp[0] = (four_tap_filter_hor_generic(c, &src[src_pos - 1]) >> shift1); // ae0,-1
        src_pos += 2 * src_stride;  //0,1
        temp[2] = (four_tap_filter_hor_generic(c, &src[src_pos - 1]) >> shift1); // ae0,1
        src_pos += src_stride;  //0,2
        temp[3] = (four_tap_filter_hor_generic(c, &src[src_pos - 1]) >> shift1); // ae0,2
        
        dst[dst_pos + 1 * dst_stride + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_generic(c, temp) + offset23) >> shift2) >> shift3); // ee0,0
      }

      if (hor_flag) {
        dst[dst_pos + 1] = fast_clip_32bit_to_pixel((temp[1] + offset3) >> shift3);
      }
    }
  }
}

void filter_inter_octpel_chroma_generic(const encoder_control_t * const encoder, kvz_pixel *src, int16_t src_stride, int width, int height, kvz_pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag)
{

  int32_t x, y;
  int32_t shift1 = KVZ_BIT_DEPTH - 8;
  int32_t shift2 = 6;
  int32_t shift3 = 14 - KVZ_BIT_DEPTH;
  int32_t offset3 = 1 << (shift3 - 1);
  int32_t offset23 = 1 << (shift2 + shift3 - 1);

  //coefficients for 1/8, 2/8, 3/8, 4/8, 5/8, 6/8 and 7/8 positions
  int8_t *c1, *c2, *c3, *c4, *c5, *c6, *c7;

  int i;
  c1 = g_chroma_filter[1];
  c2 = g_chroma_filter[2];
  c3 = g_chroma_filter[3];
  c4 = g_chroma_filter[4];
  c5 = g_chroma_filter[5];
  c6 = g_chroma_filter[6];
  c7 = g_chroma_filter[7];

  int16_t temp[7][4]; // Temporary horizontal values calculated from integer pixels


  // Loop source pixels and generate 64 filtered 1/8-pel pixels on each round
  for (y = 0; y < height; y++) {
    int dst_pos_y = (y << 3)*dst_stride;
    int src_pos_y = y*src_stride;
    for (x = 0; x < width; x++) {
      // Calculate current dst and src pixel positions
      int dst_pos = dst_pos_y + (x << 3);
      int src_pos = src_pos_y + x;
      
      // Original pixel
      dst[dst_pos] = src[src_pos];

      // Horizontal 1/8-values
      if (hor_flag && !ver_flag) {

        temp[0][1] = (four_tap_filter_hor_generic(c1, &src[src_pos - 1]) >> shift1); // ae0,0 h0
        temp[1][1] = (four_tap_filter_hor_generic(c2, &src[src_pos - 1]) >> shift1);
        temp[2][1] = (four_tap_filter_hor_generic(c3, &src[src_pos - 1]) >> shift1);
        temp[3][1] = (four_tap_filter_hor_generic(c4, &src[src_pos - 1]) >> shift1);
        temp[4][1] = (four_tap_filter_hor_generic(c5, &src[src_pos - 1]) >> shift1);
        temp[5][1] = (four_tap_filter_hor_generic(c6, &src[src_pos - 1]) >> shift1);
        temp[6][1] = (four_tap_filter_hor_generic(c7, &src[src_pos - 1]) >> shift1);

      }

      // Vertical 1/8-values
      if (ver_flag) {
        dst[dst_pos + 1 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_generic(c1, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3); //
        dst[dst_pos + 2 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_generic(c2, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
        dst[dst_pos + 3 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_generic(c3, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
        dst[dst_pos + 4 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_generic(c4, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
        dst[dst_pos + 5 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_generic(c5, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
        dst[dst_pos + 6 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_generic(c6, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
        dst[dst_pos + 7 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_generic(c7, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
      }

      // When both flags, interpolate values from temporary horizontal values
      if (hor_flag && ver_flag) {

        // Calculate temporary values
        src_pos -= 1 * src_stride;  //0,-3
        for (i = 0; i < 4; ++i) {

          temp[0][i] = (four_tap_filter_hor_generic(c1, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[1][i] = (four_tap_filter_hor_generic(c2, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[2][i] = (four_tap_filter_hor_generic(c3, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[3][i] = (four_tap_filter_hor_generic(c4, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[4][i] = (four_tap_filter_hor_generic(c5, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[5][i] = (four_tap_filter_hor_generic(c6, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[6][i] = (four_tap_filter_hor_generic(c7, &src[src_pos + i * src_stride - 1]) >> shift1);

        }


        //Calculate values from temporary horizontal 1/8-values
        for (i = 0; i<7; ++i){
          dst[dst_pos + 1 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_generic(c1, &temp[i][0]) + offset23) >> shift2) >> shift3); // ee0,0
          dst[dst_pos + 2 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_generic(c2, &temp[i][0]) + offset23) >> shift2) >> shift3);
          dst[dst_pos + 3 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_generic(c3, &temp[i][0]) + offset23) >> shift2) >> shift3);
          dst[dst_pos + 4 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_generic(c4, &temp[i][0]) + offset23) >> shift2) >> shift3);
          dst[dst_pos + 5 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_generic(c5, &temp[i][0]) + offset23) >> shift2) >> shift3);
          dst[dst_pos + 6 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_generic(c6, &temp[i][0]) + offset23) >> shift2) >> shift3);
          dst[dst_pos + 7 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_generic(c7, &temp[i][0]) + offset23) >> shift2) >> shift3);
          
        }

      }

      if (hor_flag) {
        dst[dst_pos + 1] = fast_clip_32bit_to_pixel((temp[0][1] + offset3) >> shift3);
        dst[dst_pos + 2] = fast_clip_32bit_to_pixel((temp[1][1] + offset3) >> shift3);
        dst[dst_pos + 3] = fast_clip_32bit_to_pixel((temp[2][1] + offset3) >> shift3);
        dst[dst_pos + 4] = fast_clip_32bit_to_pixel((temp[3][1] + offset3) >> shift3);
        dst[dst_pos + 5] = fast_clip_32bit_to_pixel((temp[4][1] + offset3) >> shift3);
        dst[dst_pos + 6] = fast_clip_32bit_to_pixel((temp[5][1] + offset3) >> shift3);
        dst[dst_pos + 7] = fast_clip_32bit_to_pixel((temp[6][1] + offset3) >> shift3);
      }


    }
  }
}

void extend_borders_generic(int xpos, int ypos, int mv_x, int mv_y, int off_x, int off_y, kvz_pixel *ref, int ref_width, int ref_height,
  int filterSize, int width, int height, kvz_pixel *dst) {

  int16_t mv[2] = { mv_x, mv_y };
  int halfFilterSize = filterSize >> 1;

  int dst_y; int y; int dst_x; int x; int coord_x; int coord_y;
  int8_t overflow_neg_y_temp, overflow_pos_y_temp, overflow_neg_x_temp, overflow_pos_x_temp;

  for (dst_y = 0, y = ypos - halfFilterSize; y < ((ypos + height)) + halfFilterSize; dst_y++, y++) {

    // calculate y-pixel offset
    coord_y = y + off_y + mv[1];

    // On y-overflow set coord_y accordingly
    overflow_neg_y_temp = (coord_y < 0) ? 1 : 0;
    overflow_pos_y_temp = (coord_y >= ref_height) ? 1 : 0;
    if (overflow_neg_y_temp)      coord_y = 0;
    else if (overflow_pos_y_temp) coord_y = (ref_height)-1;
    coord_y *= ref_width;

    for (dst_x = 0, x = (xpos)-halfFilterSize; x < ((xpos + width)) + halfFilterSize; dst_x++, x++) {
      coord_x = x + off_x + mv[0];

      // On x-overflow set coord_x accordingly
      overflow_neg_x_temp = (coord_x < 0) ? 1 : 0;
      overflow_pos_x_temp = (coord_x >= ref_width) ? 1 : 0;
      if (overflow_neg_x_temp)      coord_x = 0;
      else if (overflow_pos_x_temp) coord_x = ref_width - 1;

      // Store source block data (with extended borders)
      dst[dst_y*(width + filterSize) + dst_x] = ref[coord_y + coord_x];
    }
  }
}


int strategy_register_ipol_generic(void* opaque)
{
  bool success = true;

  success &= strategyselector_register(opaque, "filter_inter_quarterpel_luma", "generic", 0, &filter_inter_quarterpel_luma_generic);
  success &= strategyselector_register(opaque, "filter_inter_halfpel_chroma", "generic", 0, &filter_inter_halfpel_chroma_generic);
  success &= strategyselector_register(opaque, "filter_inter_octpel_chroma", "generic", 0, &filter_inter_octpel_chroma_generic);
  success &= strategyselector_register(opaque, "extend_borders", "generic", 0, &extend_borders_generic);

  return success;
}
