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

#include <stdlib.h>

#include "ipol-generic.h"
#include "strategyselector.h"
#include "encoder.h"

extern int8_t g_luma_filter[4][8];
extern int8_t g_chroma_filter[8][4];

int16_t eight_tap_filter_hor_generic(int8_t *filter, pixel *data)
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

int16_t eight_tap_filter_ver_generic(int8_t *filter, pixel *data, int16_t stride)
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

void filter_inter_quarterpel_luma_generic(const encoder_control * const encoder, pixel *src, int16_t src_stride, int width, int height, pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag)
{

  int32_t x, y;
  int32_t shift1 = encoder->bitdepth - 8;
  int32_t shift2 = 6;
  int32_t shift3 = 14 - encoder->bitdepth;
  int32_t offset3 = 1 << (shift3 - 1);
  int32_t offset23 = 1 << (shift2 + shift3 - 1);

  //coefficients for 1/4, 2/4 and 3/4 positions
  int8_t *c1, *c2, *c3;

  int i;
  c1 = g_luma_filter[1];
  c2 = g_luma_filter[2];
  c3 = g_luma_filter[3];

  int16_t temp[8][3];

  // Loop source pixels and generate sixteen filtered quarter-pel pixels on each round
  for (y = 0; y < height; y++) {
    int dst_pos_y = (y << 2)*dst_stride;
    int src_pos_y = y*src_stride;
    for (x = 0; x < width; x++) {
      // Calculate current dst and src pixel positions
      int dst_pos = dst_pos_y + (x << 2);
      int src_pos = src_pos_y + x;

      // Original pixel
      dst[dst_pos] = src[src_pos];

      //
      if (hor_flag && !ver_flag) {

        temp[3][0] = eight_tap_filter_hor_generic(c1, &src[src_pos - 3]) >> shift1;
        temp[3][1] = eight_tap_filter_hor_generic(c2, &src[src_pos - 3]) >> shift1;
        temp[3][2] = eight_tap_filter_hor_generic(c3, &src[src_pos - 3]) >> shift1;
      }
      // ea0,0 - needed only when ver_flag
      if (ver_flag) {
        dst[dst_pos + 1 * dst_stride] = ((eight_tap_filter_ver_generic(c1, &src[src_pos - 3 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3;
        dst[dst_pos + 2 * dst_stride] = ((eight_tap_filter_ver_generic(c2, &src[src_pos - 3 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3;
        dst[dst_pos + 3 * dst_stride] = ((eight_tap_filter_ver_generic(c3, &src[src_pos - 3 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3;
      }

      // When both flags, we use _only_ this pixel (but still need ae0,0 for it)
      if (hor_flag && ver_flag) {

        // Calculate temporary values..
        src_pos -= 3 * src_stride;  //0,-3
        for (i = 0; i < 8; ++i) {

          temp[i][0] = eight_tap_filter_hor_generic(c1, &src[src_pos - 3]) >> shift1; // h0(0,-3+i)
          temp[i][1] = eight_tap_filter_hor_generic(c2, &src[src_pos - 3]) >> shift1; // h1(0,-3+i)
          temp[i][2] = eight_tap_filter_hor_generic(c3, &src[src_pos - 3]) >> shift1; // h2(0,-3+i)
          src_pos += src_stride;
        }



        for (i = 0; i<3; ++i){
          dst[dst_pos + 1 * dst_stride + i + 1] = ((eight_tap_filter_ver_16bit_generic(c1, &temp[0][i], 3) + offset23) >> shift2) >> shift3;
          dst[dst_pos + 2 * dst_stride + i + 1] = ((eight_tap_filter_ver_16bit_generic(c2, &temp[0][i], 3) + offset23) >> shift2) >> shift3;
          dst[dst_pos + 3 * dst_stride + i + 1] = ((eight_tap_filter_ver_16bit_generic(c3, &temp[0][i], 3) + offset23) >> shift2) >> shift3;

        }

      }

      if (hor_flag) {
        dst[dst_pos + 1] = (temp[3][0] + offset3) >> shift3;
        dst[dst_pos + 2] = (temp[3][1] + offset3) >> shift3;
        dst[dst_pos + 3] = (temp[3][2] + offset3) >> shift3;
      }


    }
  }

  //Clamp values to bitdepth
  for (i = 0; i < width*height * 16; ++i) {
    if (dst[i] >((1 << encoder->bitdepth) - 1)) dst[i] = (pixel)((1 << encoder->bitdepth) - 1);
    if (dst[i] < 0) dst[i] = 0;
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
void filter_inter_halfpel_chroma_generic(const encoder_control * const encoder, pixel *src, int16_t src_stride, int width, int height, pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag)
{
  /* ____________
  * | B0,0|ae0,0|
  * |ea0,0|ee0,0|
  *
  * ae0,0 = (-4*B-1,0  + 36*B0,0  + 36*B1,0  - 4*B2,0)  >> shift1
  * ea0,0 = (-4*B0,-1  + 36*B0,0  + 36*B0,1  - 4*B0,2)  >> shift1
  * ee0,0 = (-4*ae0,-1 + 36*ae0,0 + 36*ae0,1 - 4*ae0,2) >> shift2
  */
  int i = 0;
  int32_t x, y;
  int32_t shift1 = encoder->bitdepth - 8;
  int32_t shift2 = 6;
  int32_t shift3 = 14 - encoder->bitdepth;
  int32_t offset3 = 1 << (shift3 - 1);
  int32_t offset23 = 1 << (shift2 + shift3 - 1);

  // Loop source pixels and generate four filtered half-pel pixels on each round
  for (y = 0; y < height; y++) {
    int dst_pos_y = (y << 1)*dst_stride;
    int src_pos_y = y*src_stride;
    for (x = 0; x < width; x++) {
      // Calculate current dst and src pixel positions
      int dst_pos = dst_pos_y + (x << 1);
      int src_pos = src_pos_y + x;

      // Temporary variables..
      int32_t ae_temp = 0;

      // Original pixel (not really needed)
      dst[dst_pos] = src[src_pos]; //B0,0

      // ae0,0 - We need this only when hor_flag and for ee0,0
      if (hor_flag) {
        ae_temp = ((-4 * src[src_pos - 1] + 36 * src[src_pos] + 36 * src[src_pos + 1] - 4 * src[src_pos + 2]) >> shift1); // ae0,0
      }
      // ea0,0 - needed only when ver_flag
      if (ver_flag) {
        dst[dst_pos + 1 * dst_stride] = (((-4 * src[src_pos - src_stride] + 36 * src[src_pos] + 36 * src[src_pos + src_stride]
          - 4 * src[src_pos + 2 * src_stride]) >> shift1) + (1 << (shift3 - 1))) >> shift3; // ea0,0
      }

      // When both flags, we use _only_ this pixel (but still need ae0,0 for it)
      if (hor_flag && ver_flag) {
        int32_t ae_temp1, ae_temp2, ae_temp3;
        // Calculate temporary values..
        //TODO: optimization, store these values
        src_pos -= src_stride;  //0,-1
        ae_temp1 = ((-4 * src[src_pos - 1] + 36 * src[src_pos] + 36 * src[src_pos + 1] - 4 * src[src_pos + 2]) >> shift1); // ae0,-1
        src_pos += 2 * src_stride;  //0,1
        ae_temp2 = ((-4 * src[src_pos - 1] + 36 * src[src_pos] + 36 * src[src_pos + 1] - 4 * src[src_pos + 2]) >> shift1); // ae0,1
        src_pos += src_stride;  //0,2
        ae_temp3 = ((-4 * src[src_pos - 1] + 36 * src[src_pos] + 36 * src[src_pos + 1] - 4 * src[src_pos + 2]) >> shift1); // ae0,2

        dst[dst_pos + 1 * dst_stride + 1] = (((-4 * ae_temp1 + 36 * ae_temp + 36 * ae_temp2 - 4 * ae_temp3) + offset23) >> shift2) >> shift3; // ee0,0
      }

      if (hor_flag) {
        dst[dst_pos + 1] = (ae_temp + offset3) >> shift3;
      }
    }
  }
  //Clamp values to bitdepth
  for (i = 0; i < width*height * 4; ++i) {
    if (dst[i] >((1 << encoder->bitdepth) - 1)) dst[i] = (pixel)((1 << encoder->bitdepth) - 1);
    if (dst[i] < 0) dst[i] = 0;
  }
}

void filter_inter_octpel_chroma_generic(const encoder_control * const encoder, pixel *src, int16_t src_stride, int width, int height, pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag)
{

  int32_t x, y;
  int32_t shift1 = encoder->bitdepth - 8;
  int32_t shift2 = 6;
  int32_t shift3 = 14 - encoder->bitdepth;
  int32_t offset3 = 1 << (shift3 - 1);
  int32_t offset23 = 1 << (shift2 + shift3 - 1);

  //coefficients for 1/8, 2/8, 3/8, 4/8, 5/8, 6/8 and 7/8 positions
  int8_t c1[4], c2[4], c3[4], c4[4], c5[4], c6[4], c7[4];

  int i;
  for (i = 0; i < 4; ++i) {
    c1[i] = g_chroma_filter[1][i];
    c2[i] = g_chroma_filter[2][i];
    c3[i] = g_chroma_filter[3][i];
    c4[i] = g_chroma_filter[4][i];
    c5[i] = g_chroma_filter[5][i];
    c6[i] = g_chroma_filter[6][i];
    c7[i] = g_chroma_filter[7][i];
  }

  // Loop source pixels and generate 64 filtered 1/8-pel pixels on each round
  for (y = 0; y < height; y++) {
    int dst_pos_y = (y << 3)*dst_stride;
    int src_pos_y = y*src_stride;
    for (x = 0; x < width; x++) {
      // Calculate current dst and src pixel positions
      int dst_pos = dst_pos_y + (x << 3);
      int src_pos = src_pos_y + x;

      // Temporary horizontally interpolated postions
      int32_t h_temp[7] = { 0, 0, 0, 0, 0, 0, 0 };

      // Original pixel
      dst[dst_pos] = src[src_pos];

      // Horizontal 1/8-values
      if (hor_flag) {

        h_temp[0] = ((c1[0] * src[src_pos - 1]
          + c1[1] * src[src_pos]
          + c1[2] * src[src_pos + 1]
          + c1[3] * src[src_pos + 2]) >> shift1); // ae0,0 h0

        h_temp[1] = ((c2[0] * src[src_pos - 1]
          + c2[1] * src[src_pos]
          + c2[2] * src[src_pos + 1]
          + c2[3] * src[src_pos + 2]) >> shift1); // ae0,0 h1

        h_temp[2] = ((c3[0] * src[src_pos - 1]
          + c3[1] * src[src_pos]
          + c3[2] * src[src_pos + 1]
          + c3[3] * src[src_pos + 2]) >> shift1); // ae0,0 h2

        h_temp[3] = ((c4[0] * src[src_pos - 1]
          + c4[1] * src[src_pos]
          + c4[2] * src[src_pos + 1]
          + c4[3] * src[src_pos + 2]) >> shift1); // ae0,0 h2

        h_temp[4] = ((c5[0] * src[src_pos - 1]
          + c5[1] * src[src_pos]
          + c5[2] * src[src_pos + 1]
          + c5[3] * src[src_pos + 2]) >> shift1); // ae0,0 h2

        h_temp[5] = ((c6[0] * src[src_pos - 1]
          + c6[1] * src[src_pos]
          + c6[2] * src[src_pos + 1]
          + c6[3] * src[src_pos + 2]) >> shift1); // ae0,0 h2

        h_temp[6] = ((c7[0] * src[src_pos - 1]
          + c7[1] * src[src_pos]
          + c7[2] * src[src_pos + 1]
          + c7[3] * src[src_pos + 2]) >> shift1); // ae0,0 h2
      }

      // Vertical 1/8-values
      if (ver_flag) {
        dst[dst_pos + 1 * dst_stride] = (((c1[0] * src[src_pos - 1 * src_stride]
          + c1[1] * src[src_pos]
          + c1[2] * src[src_pos + 1 * src_stride]
          + c1[3] * src[src_pos + 2 * src_stride]) >> shift1)
          + (1 << (shift3 - 1))) >> shift3; //

        dst[dst_pos + 2 * dst_stride] = (((c2[0] * src[src_pos - 1 * src_stride]
          + c2[1] * src[src_pos]
          + c2[2] * src[src_pos + 1 * src_stride]
          + c2[3] * src[src_pos + 2 * src_stride]) >> shift1)
          + (1 << (shift3 - 1))) >> shift3; //

        dst[dst_pos + 3 * dst_stride] = (((c3[0] * src[src_pos - 1 * src_stride]
          + c3[1] * src[src_pos]
          + c3[2] * src[src_pos + 1 * src_stride]
          + c3[3] * src[src_pos + 2 * src_stride]) >> shift1)
          + (1 << (shift3 - 1))) >> shift3; //

        dst[dst_pos + 4 * dst_stride] = (((c4[0] * src[src_pos - 1 * src_stride]
          + c4[1] * src[src_pos]
          + c4[2] * src[src_pos + 1 * src_stride]
          + c4[3] * src[src_pos + 2 * src_stride]) >> shift1)
          + (1 << (shift3 - 1))) >> shift3; //

        dst[dst_pos + 5 * dst_stride] = (((c5[0] * src[src_pos - 1 * src_stride]
          + c5[1] * src[src_pos]
          + c5[2] * src[src_pos + 1 * src_stride]
          + c5[3] * src[src_pos + 2 * src_stride]) >> shift1)
          + (1 << (shift3 - 1))) >> shift3; //

        dst[dst_pos + 6 * dst_stride] = (((c6[0] * src[src_pos - 1 * src_stride]
          + c6[1] * src[src_pos]
          + c6[2] * src[src_pos + 1 * src_stride]
          + c6[3] * src[src_pos + 2 * src_stride]) >> shift1)
          + (1 << (shift3 - 1))) >> shift3; //

        dst[dst_pos + 7 * dst_stride] = (((c7[0] * src[src_pos - 1 * src_stride]
          + c7[1] * src[src_pos]
          + c7[2] * src[src_pos + 1 * src_stride]
          + c7[3] * src[src_pos + 2 * src_stride]) >> shift1)
          + (1 << (shift3 - 1))) >> shift3; //
      }

      // When both flags, interpolate values from temporary horizontal values
      if (hor_flag && ver_flag) {

        int32_t temp[3][7]; // Temporary horizontal values calculated from integer pixels

        // Calculate temporary values
        src_pos -= 1 * src_stride;  //0,-3
        for (i = 0; i < 3; ++i) {

          temp[i][0] = ((c1[0] * src[src_pos - 1] + c1[1] * src[src_pos]
            + c1[2] * src[src_pos + 1] + c1[3] * src[src_pos + 2])
            >> shift1); // h0(0,-3+i)

          temp[i][1] = ((c2[0] * src[src_pos - 1] + c2[1] * src[src_pos]
            + c2[2] * src[src_pos + 1] + c2[3] * src[src_pos + 2])
            >> shift1); // h1(0,-3+i)

          temp[i][2] = ((c3[0] * src[src_pos - 1] + c3[1] * src[src_pos]
            + c3[2] * src[src_pos + 1] + c3[3] * src[src_pos + 2])
            >> shift1); // h2(0,-3+i)

          temp[i][3] = ((c4[0] * src[src_pos - 1] + c4[1] * src[src_pos]
            + c4[2] * src[src_pos + 1] + c4[3] * src[src_pos + 2])
            >> shift1); // h2(0,-3+i)

          temp[i][4] = ((c5[0] * src[src_pos - 1] + c5[1] * src[src_pos]
            + c5[2] * src[src_pos + 1] + c5[3] * src[src_pos + 2])
            >> shift1); // h2(0,-3+i)

          temp[i][5] = ((c6[0] * src[src_pos - 1] + c6[1] * src[src_pos]
            + c6[2] * src[src_pos + 1] + c6[3] * src[src_pos + 2])
            >> shift1); // h2(0,-3+i)

          temp[i][6] = ((c7[0] * src[src_pos - 1] + c7[1] * src[src_pos]
            + c7[2] * src[src_pos + 1] + c7[3] * src[src_pos + 2])
            >> shift1); // h2(0,-3+i)

          if (i == 0) {
            //Skip calculating h_temp again
            src_pos += 2 * src_stride;
          }
          else {
            src_pos += src_stride;
          }
        }


        //Calculate values from temporary horizontal 1/8-values
        for (i = 0; i<7; ++i){
          dst[dst_pos + 1 * dst_stride + i + 1] = (((c1[0] * temp[0][i] + c1[1] * h_temp[i]
            + c1[2] * temp[1][i] + c1[3] * temp[2][i])
            + offset23) >> shift2) >> shift3; // ee0,0

          dst[dst_pos + 2 * dst_stride + i + 1] = (((c2[0] * temp[0][i] + c2[1] * h_temp[i]
            + c2[2] * temp[1][i] + c2[3] * temp[2][i])
            + offset23) >> shift2) >> shift3; // ee0,0

          dst[dst_pos + 3 * dst_stride + i + 1] = (((c3[0] * temp[0][i] + c3[1] * h_temp[i]
            + c3[2] * temp[1][i] + c3[3] * temp[2][i])
            + offset23) >> shift2) >> shift3; // ee0,0

          dst[dst_pos + 4 * dst_stride + i + 1] = (((c4[0] * temp[0][i] + c4[1] * h_temp[i]
            + c4[2] * temp[1][i] + c4[3] * temp[2][i])
            + offset23) >> shift2) >> shift3; // ee0,0

          dst[dst_pos + 5 * dst_stride + i + 1] = (((c5[0] * temp[0][i] + c5[1] * h_temp[i]
            + c5[2] * temp[1][i] + c5[3] * temp[2][i])
            + offset23) >> shift2) >> shift3; // ee0,0

          dst[dst_pos + 6 * dst_stride + i + 1] = (((c6[0] * temp[0][i] + c6[1] * h_temp[i]
            + c6[2] * temp[1][i] + c6[3] * temp[2][i])
            + offset23) >> shift2) >> shift3; // ee0,0

          dst[dst_pos + 7 * dst_stride + i + 1] = (((c7[0] * temp[0][i] + c7[1] * h_temp[i]
            + c7[2] * temp[1][i] + c7[3] * temp[2][i])
            + offset23) >> shift2) >> shift3; // ee0,0

        }

      }

      if (hor_flag) {
        dst[dst_pos + 1] = (h_temp[0] + offset3) >> shift3;
        dst[dst_pos + 2] = (h_temp[1] + offset3) >> shift3;
        dst[dst_pos + 3] = (h_temp[2] + offset3) >> shift3;
        dst[dst_pos + 4] = (h_temp[3] + offset3) >> shift3;
        dst[dst_pos + 5] = (h_temp[4] + offset3) >> shift3;
        dst[dst_pos + 6] = (h_temp[5] + offset3) >> shift3;
        dst[dst_pos + 7] = (h_temp[6] + offset3) >> shift3;
      }


    }
  }

  //Clamp values to bitdepth
  for (i = 0; i < width*height * 64; ++i) {
    if (dst[i] >((1 << encoder->bitdepth) - 1)) dst[i] = (pixel)((1 << encoder->bitdepth) - 1);
    if (dst[i] < 0) dst[i] = 0;
  }
}

void extend_borders_generic(int xpos, int ypos, int mv_x, int mv_y, int off_x, int off_y, pixel *ref, int ref_width, int ref_height,
  int filterSize, int width, int height, pixel *dst) {

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
