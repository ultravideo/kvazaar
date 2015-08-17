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

#include "ipol-avx2.h"
#include "strategyselector.h"

#if COMPILE_INTEL_AVX2
#include <stdlib.h>

#include <immintrin.h>


#include "encoder.h"
#include "strategies/generic/picture-generic.h"


#define FILTER_OFFSET 3
#define FILTER_SIZE 8

#define MAX_HEIGHT (4 * (LCU_WIDTH + 1) + FILTER_SIZE)
#define MAX_WIDTH ((LCU_WIDTH + 1) + FILTER_SIZE)

extern int8_t g_luma_filter[4][8];
extern int8_t g_chroma_filter[8][4];

void eight_tap_filter_x8_and_flip(__m128i data01, __m128i data23, __m128i data45, __m128i data67, __m128i* filter, __m128i* dst)
{
  __m128i a, b, c, d;
  __m128i fir = _mm_broadcastq_epi64(_mm_loadl_epi64(filter));

  a = _mm_maddubs_epi16(data01, fir);
  b = _mm_maddubs_epi16(data23, fir);
  a = _mm_hadd_epi16(a, b);

  c = _mm_maddubs_epi16(data45, fir);
  d = _mm_maddubs_epi16(data67, fir);
  c = _mm_hadd_epi16(c, d);

  a = _mm_hadd_epi16(a, c);

  _mm_storeu_si128(dst, a);
}

__m128i eight_tap_filter_x4_and_flip_16bit(__m128i data0, __m128i data1, __m128i data2, __m128i data3, __m128i* filter)
{
  __m128i a, b, c, d;
  __m128i fir = _mm_cvtepi8_epi16(_mm_loadu_si128((__m128i*)(filter)));

  a = _mm_madd_epi16(data0, fir);
  b = _mm_madd_epi16(data1, fir);
  a = _mm_hadd_epi32(a, b);

  c = _mm_madd_epi16(data2, fir);
  d = _mm_madd_epi16(data3, fir);
  c = _mm_hadd_epi32(c, d);

  a = _mm_hadd_epi32(a, c);

  return a;
}

void eight_tap_filter_and_flip_avx2(int8_t filter[4][8], kvz_pixel *src, int16_t src_stride, int16_t* __restrict dst)
{

  //Load 2 rows per xmm register
  __m128i rows01 = _mm_loadl_epi64((__m128i*)(src + 0 * src_stride));
  rows01 = _mm_castpd_si128(_mm_loadh_pd(_mm_castsi128_pd(rows01), (double*)(src + 1 * src_stride)));

  __m128i rows23 = _mm_loadl_epi64((__m128i*)(src + 2 * src_stride));
  rows23 = _mm_castpd_si128(_mm_loadh_pd(_mm_castsi128_pd(rows23), (double*)(src + 3 * src_stride)));

  __m128i rows45 = _mm_loadl_epi64((__m128i*)(src + 4 * src_stride));
  rows45 = _mm_castpd_si128(_mm_loadh_pd(_mm_castsi128_pd(rows45), (double*)(src + 5 * src_stride)));

  __m128i rows67 = _mm_loadl_epi64((__m128i*)(src + 6 * src_stride));
  rows67 = _mm_castpd_si128(_mm_loadh_pd(_mm_castsi128_pd(rows67), (double*)(src + 7 * src_stride)));

  //Filter rows
  const int dst_stride = MAX_WIDTH;
  eight_tap_filter_x8_and_flip(rows01, rows23, rows45, rows67, (__m128i*)(&filter[0]), (__m128i*)(dst + 0));
  eight_tap_filter_x8_and_flip(rows01, rows23, rows45, rows67, (__m128i*)(&filter[1]), (__m128i*)(dst + 1 * dst_stride));
  eight_tap_filter_x8_and_flip(rows01, rows23, rows45, rows67, (__m128i*)(&filter[2]), (__m128i*)(dst + 2 * dst_stride));
  eight_tap_filter_x8_and_flip(rows01, rows23, rows45, rows67, (__m128i*)(&filter[3]), (__m128i*)(dst + 3 * dst_stride));
}

static INLINE void eight_tap_filter_and_flip_16bit_avx2(int8_t filter[4][8], int16_t *src, int16_t src_stride, int offset, int combined_shift, kvz_pixel* __restrict dst, int16_t dst_stride)
{

  //Load a row per xmm register
  __m128i row0 = _mm_loadu_si128((__m128i*)(src + 0 * src_stride));
  __m128i row1 = _mm_loadu_si128((__m128i*)(src + 1 * src_stride));
  __m128i row2 = _mm_loadu_si128((__m128i*)(src + 2 * src_stride));
  __m128i row3 = _mm_loadu_si128((__m128i*)(src + 3 * src_stride));

  //Filter rows
  union {
    __m128i vector;
    int32_t array[4];
  } temp[4];

  temp[0].vector = eight_tap_filter_x4_and_flip_16bit(row0, row1, row2, row3, (__m128i*)(&filter[0]));
  temp[1].vector = eight_tap_filter_x4_and_flip_16bit(row0, row1, row2, row3, (__m128i*)(&filter[1]));
  temp[2].vector = eight_tap_filter_x4_and_flip_16bit(row0, row1, row2, row3, (__m128i*)(&filter[2]));
  temp[3].vector = eight_tap_filter_x4_and_flip_16bit(row0, row1, row2, row3, (__m128i*)(&filter[3]));

  __m128i packed_offset = _mm_set1_epi32(offset);

  temp[0].vector = _mm_add_epi32(temp[0].vector, packed_offset);
  temp[0].vector = _mm_srai_epi32(temp[0].vector, combined_shift);
  temp[1].vector = _mm_add_epi32(temp[1].vector, packed_offset);
  temp[1].vector = _mm_srai_epi32(temp[1].vector, combined_shift);

  temp[0].vector = _mm_packus_epi32(temp[0].vector, temp[1].vector);

  temp[2].vector = _mm_add_epi32(temp[2].vector, packed_offset);
  temp[2].vector = _mm_srai_epi32(temp[2].vector, combined_shift);
  temp[3].vector = _mm_add_epi32(temp[3].vector, packed_offset);
  temp[3].vector = _mm_srai_epi32(temp[3].vector, combined_shift);

  temp[2].vector = _mm_packus_epi32(temp[2].vector, temp[3].vector);

  temp[0].vector = _mm_packus_epi16(temp[0].vector, temp[2].vector);

  int32_t* four_pixels = (int32_t*)&(dst[0 * dst_stride]);
  *four_pixels = temp[0].array[0];

  four_pixels = (int32_t*)&(dst[1 * dst_stride]);
  *four_pixels = _mm_extract_epi32(temp[0].vector, 1);

  four_pixels = (int32_t*)&(dst[2 * dst_stride]);
  *four_pixels = _mm_extract_epi32(temp[0].vector, 2);

  four_pixels = (int32_t*)&(dst[3 * dst_stride]);
  *four_pixels = _mm_extract_epi32(temp[0].vector, 3);


}

int16_t eight_tap_filter_hor_avx2(int8_t *filter, kvz_pixel *data)
{
  union {
    __m128i vector;
    int16_t array[8];
  } sample;

  __m128i packed_data = _mm_loadu_si128((__m128i*)data);
  __m128i packed_filter = _mm_loadu_si128((__m128i*)filter);

  sample.vector = _mm_maddubs_epi16(packed_data, packed_filter);
  sample.vector = _mm_hadd_epi16(sample.vector, sample.vector);
  sample.vector = _mm_hadd_epi16(sample.vector, sample.vector);

  return sample.array[0];
}

int32_t eight_tap_filter_hor_16bit_avx2(int8_t *filter, int16_t *data)
{
  int32_t temp = 0;
  for (int i = 0; i < 8; ++i)
  {
    temp += filter[i] * data[i];
  }

  return temp;
}

int16_t eight_tap_filter_ver_avx2(int8_t *filter, kvz_pixel *data, int16_t stride)
{
  int16_t temp = 0;
  for (int i = 0; i < 8; ++i)
  {
    temp += filter[i] * data[stride * i];
  }

  return temp;
}

int32_t eight_tap_filter_ver_16bit_avx2(int8_t *filter, int16_t *data, int16_t stride)
{
  int32_t temp = 0;
  for (int i = 0; i < 8; ++i)
  {
    temp += filter[i] * data[stride * i];
  }

  return temp;
}

int16_t four_tap_filter_hor_avx2(int8_t *filter, kvz_pixel *data)
{
  int16_t temp = 0;
  for (int i = 0; i < 4; ++i)
  {
    temp += filter[i] * data[i];
  }

  return temp;
}

int32_t four_tap_filter_hor_16bit_avx2(int8_t *filter, int16_t *data)
{
  int32_t temp = 0;
  for (int i = 0; i < 4; ++i)
  {
    temp += filter[i] * data[i];
  }

  return temp;
}

int16_t four_tap_filter_ver_avx2(int8_t *filter, kvz_pixel *data, int16_t stride)
{
  int16_t temp = 0;
  for (int i = 0; i < 4; ++i)
  {
    temp += filter[i] * data[stride * i];
  }

  return temp;
}

int32_t four_tap_filter_ver_16bit_avx2(int8_t *filter, int16_t *data, int16_t stride)
{
  int32_t temp = 0;
  for (int i = 0; i < 4; ++i)
  {
    temp += filter[i] * data[stride * i];
  }

  return temp;
}

void filter_inter_quarterpel_luma_avx2(const encoder_control_t * const encoder, kvz_pixel *src, int16_t src_stride, int width, int height, kvz_pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag)
{

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

  int16_t flipped_hor_filtered[MAX_HEIGHT][MAX_WIDTH];

  // Filter horizontally and flip x and y
  for (x = 0; x < width; ++x) {
    for (y = 0; y < height; y += 8) {
      int ypos = y - FILTER_OFFSET;
      int xpos = x - FILTER_OFFSET;

      eight_tap_filter_and_flip_avx2(g_luma_filter, &src[src_stride*ypos + xpos], src_stride, (int16_t*)&(flipped_hor_filtered[4 * x + 0][y]));
    
    }

    for (; y < height + FILTER_SIZE - 1; ++y) {
      int ypos = y - FILTER_OFFSET;
      int xpos = x - FILTER_OFFSET;
      flipped_hor_filtered[4 * x + 0][y] = eight_tap_filter_hor_avx2(c0, &src[src_stride*ypos + xpos]) << shift1;
      flipped_hor_filtered[4 * x + 1][y] = eight_tap_filter_hor_avx2(c1, &src[src_stride*ypos + xpos]) << shift1;
      flipped_hor_filtered[4 * x + 2][y] = eight_tap_filter_hor_avx2(c2, &src[src_stride*ypos + xpos]) << shift1;
      flipped_hor_filtered[4 * x + 3][y] = eight_tap_filter_hor_avx2(c3, &src[src_stride*ypos + xpos]) << shift1;
    }
  }

  // Filter vertically and flip x and y
  for (y = 0; y < height; ++y) {
    for (x = 0; x < 4 * width - 3; x += 4) {

     eight_tap_filter_and_flip_16bit_avx2(g_luma_filter, &flipped_hor_filtered[x][y], MAX_WIDTH, offset23, shift2 + shift3, &(dst[(4 * y + 0)*dst_stride + x]), dst_stride);

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
void filter_inter_halfpel_chroma_avx2(const encoder_control_t * const encoder, kvz_pixel *src, int16_t src_stride, int width, int height, kvz_pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag)
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
        temp[1] = four_tap_filter_hor_avx2(c, &src[src_pos - 1]) >> shift1; // ae0,0
      }
      // ea0,0 - needed only when ver_flag
      if (ver_flag) {
        dst[dst_pos + 1 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_avx2(c, &src[src_pos - src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3); // ea0,0
      }

      // When both flags, we use _only_ this pixel (but still need ae0,0 for it)
      if (hor_flag && ver_flag) {
        // Calculate temporary values..
        src_pos -= src_stride;  //0,-1
        temp[0] = (four_tap_filter_hor_avx2(c, &src[src_pos - 1]) >> shift1); // ae0,-1
        src_pos += 2 * src_stride;  //0,1
        temp[2] = (four_tap_filter_hor_avx2(c, &src[src_pos - 1]) >> shift1); // ae0,1
        src_pos += src_stride;  //0,2
        temp[3] = (four_tap_filter_hor_avx2(c, &src[src_pos - 1]) >> shift1); // ae0,2

        dst[dst_pos + 1 * dst_stride + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_avx2(c, temp) + offset23) >> shift2) >> shift3); // ee0,0
      }

      if (hor_flag) {
        dst[dst_pos + 1] = fast_clip_32bit_to_pixel((temp[1] + offset3) >> shift3);
      }
    }
  }
}

void filter_inter_octpel_chroma_avx2(const encoder_control_t * const encoder, kvz_pixel *src, int16_t src_stride, int width, int height, kvz_pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag)
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

        temp[0][1] = (four_tap_filter_hor_avx2(c1, &src[src_pos - 1]) >> shift1); // ae0,0 h0
        temp[1][1] = (four_tap_filter_hor_avx2(c2, &src[src_pos - 1]) >> shift1);
        temp[2][1] = (four_tap_filter_hor_avx2(c3, &src[src_pos - 1]) >> shift1);
        temp[3][1] = (four_tap_filter_hor_avx2(c4, &src[src_pos - 1]) >> shift1);
        temp[4][1] = (four_tap_filter_hor_avx2(c5, &src[src_pos - 1]) >> shift1);
        temp[5][1] = (four_tap_filter_hor_avx2(c6, &src[src_pos - 1]) >> shift1);
        temp[6][1] = (four_tap_filter_hor_avx2(c7, &src[src_pos - 1]) >> shift1);

      }

      // Vertical 1/8-values
      if (ver_flag) {
        dst[dst_pos + 1 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_avx2(c1, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3); //
        dst[dst_pos + 2 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_avx2(c2, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
        dst[dst_pos + 3 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_avx2(c3, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
        dst[dst_pos + 4 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_avx2(c4, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
        dst[dst_pos + 5 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_avx2(c5, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
        dst[dst_pos + 6 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_avx2(c6, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
        dst[dst_pos + 7 * dst_stride] = fast_clip_32bit_to_pixel(((four_tap_filter_ver_avx2(c7, &src[src_pos - 1 * src_stride], src_stride) >> shift1) + (1 << (shift3 - 1))) >> shift3);
      }

      // When both flags, interpolate values from temporary horizontal values
      if (hor_flag && ver_flag) {

        // Calculate temporary values
        src_pos -= 1 * src_stride;  //0,-3
        for (i = 0; i < 4; ++i) {

          temp[0][i] = (four_tap_filter_hor_avx2(c1, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[1][i] = (four_tap_filter_hor_avx2(c2, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[2][i] = (four_tap_filter_hor_avx2(c3, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[3][i] = (four_tap_filter_hor_avx2(c4, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[4][i] = (four_tap_filter_hor_avx2(c5, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[5][i] = (four_tap_filter_hor_avx2(c6, &src[src_pos + i * src_stride - 1]) >> shift1);
          temp[6][i] = (four_tap_filter_hor_avx2(c7, &src[src_pos + i * src_stride - 1]) >> shift1);

        }


        //Calculate values from temporary horizontal 1/8-values
        for (i = 0; i<7; ++i){
          dst[dst_pos + 1 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_avx2(c1, &temp[i][0]) + offset23) >> shift2) >> shift3); // ee0,0
          dst[dst_pos + 2 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_avx2(c2, &temp[i][0]) + offset23) >> shift2) >> shift3);
          dst[dst_pos + 3 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_avx2(c3, &temp[i][0]) + offset23) >> shift2) >> shift3);
          dst[dst_pos + 4 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_avx2(c4, &temp[i][0]) + offset23) >> shift2) >> shift3);
          dst[dst_pos + 5 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_avx2(c5, &temp[i][0]) + offset23) >> shift2) >> shift3);
          dst[dst_pos + 6 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_avx2(c6, &temp[i][0]) + offset23) >> shift2) >> shift3);
          dst[dst_pos + 7 * dst_stride + i + 1] = fast_clip_32bit_to_pixel(((four_tap_filter_hor_16bit_avx2(c7, &temp[i][0]) + offset23) >> shift2) >> shift3);

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

void get_extended_block_avx2(int xpos, int ypos, int mv_x, int mv_y, int off_x, int off_y, kvz_pixel *ref, int ref_width, int ref_height,
  int filter_size, int width, int height, kvz_extended_block *out) {

  int half_filter_size = filter_size >> 1;

  out->buffer = ref + (ypos - half_filter_size + off_y + mv_y) * ref_width + (xpos - half_filter_size + off_x + mv_x);
  out->stride = ref_width;
  out->orig_topleft = out->buffer + out->stride * half_filter_size + half_filter_size;
  out->malloc_used = 0;

  int min_y = ypos - half_filter_size + off_y + mv_y;
  int max_y = min_y + height + filter_size;
  int out_of_bounds_y = (min_y < 0) || (max_y >= ref_height);

  int min_x = xpos - half_filter_size + off_x + mv_x;
  int max_x = min_x + width + filter_size;
  int out_of_bounds_x = (min_x < 0) || (max_x >= ref_width);

  int sample_out_of_bounds = out_of_bounds_y || out_of_bounds_x;

  if (sample_out_of_bounds){
    out->buffer = MALLOC(kvz_pixel, (width + filter_size) * (width + filter_size));
    if (!out->buffer){
      fprintf(stderr, "Memory allocation failed!\n");
      assert(0);
    }
    out->stride = width + filter_size;
    out->orig_topleft = out->buffer + out->stride * half_filter_size + half_filter_size;
    out->malloc_used = 1;

    int dst_y; int y; int dst_x; int x; int coord_x; int coord_y;

    for (dst_y = 0, y = ypos - half_filter_size; y < ((ypos + height)) + half_filter_size; dst_y++, y++) {

      // calculate y-pixel offset
      coord_y = y + off_y + mv_y;
      coord_y = CLIP(0, (ref_height)-1, coord_y);
      coord_y *= ref_width;

      if (!out_of_bounds_x){
        memcpy(&out->buffer[dst_y*(width + filter_size) + 0], &ref[coord_y + min_x], (width + filter_size) * sizeof(kvz_pixel));
      } else {
        for (dst_x = 0, x = (xpos)-half_filter_size; x < ((xpos + width)) + half_filter_size; dst_x++, x++) {

          coord_x = x + off_x + mv_x;
          coord_x = CLIP(0, (ref_width)-1, coord_x);

          // Store source block data (with extended borders)
          out->buffer[dst_y*(width + filter_size) + dst_x] = ref[coord_y + coord_x];
        }
      }
    }
  }
}

#endif //COMPILE_INTEL_AVX2

int strategy_register_ipol_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_INTEL_AVX2
  if (bitdepth == 8){
    success &= strategyselector_register(opaque, "filter_inter_quarterpel_luma", "avx2", 40, &filter_inter_quarterpel_luma_avx2);
    success &= strategyselector_register(opaque, "filter_inter_halfpel_chroma", "avx2", 40, &filter_inter_halfpel_chroma_avx2);
    success &= strategyselector_register(opaque, "filter_inter_octpel_chroma", "avx2", 40, &filter_inter_octpel_chroma_avx2);
  }
  success &= strategyselector_register(opaque, "get_extended_block", "avx2", 40, &get_extended_block_avx2);
#endif //COMPILE_INTEL_AVX2
  return success;
}
