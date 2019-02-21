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

#include "strategies/avx2/sao-avx2.h"

#if COMPILE_INTEL_AVX2
#include <immintrin.h>

#include "cu.h"
#include "encoder.h"
#include "encoderstate.h"
#include "kvazaar.h"
#include "sao.h"
#include "strategyselector.h"


// These optimizations are based heavily on sao-generic.c.
// Might be useful to check that if (when) this file
// is difficult to understand.

// Mapping of edge_idx values to eo-classes.
static int sao_calc_eo_cat(kvz_pixel a, kvz_pixel b, kvz_pixel c)
{
 // Mapping relationships between a, b and c to eo_idx.
 static const int sao_eo_idx_to_eo_category[] = { 1, 2, 0, 3, 4 };

 int eo_idx = 2 + SIGN3((int)c - (int)a) + SIGN3((int)c - (int)b);

 //printf("%d ", SIGN3((int)c - (int)a));
 return sao_eo_idx_to_eo_category[eo_idx];
}

static int sao_calc_eo_cat_test(kvz_pixel a, kvz_pixel b, kvz_pixel c)
{
 // Mapping relationships between a, b and c to eo_idx.
 static const int sao_eo_idx_to_eo_category[] = { 1, 2, 0, 3, 4 };

 int eo_idx = 2 + SIGN3((int)c - (int)a) + SIGN3((int)c - (int)b);

 //printf("%d ", sao_eo_idx_to_eo_category[eo_idx]);
 return sao_eo_idx_to_eo_category[eo_idx];
}

// Mapping of edge_idx values to eo-classes.
static __m256i sao_calc_eo_cat_avx2(__m128i vector_a_epi8, __m128i vector_b_epi8, __m128i vector_c_epi8)
{
 // Mapping relationships between a, b and c to eo_idx.
 __m256i vector_sao_eo_idx_to_eo_category_epi32 = _mm256_setr_epi32(1, 2, 0, 3, 4, 0, 0, 0);

 __m256i eo_idx_epi32 = _mm256_set1_epi32(2);
 __m256i vector_a_epi32 = _mm256_cvtepu8_epi32(vector_a_epi8);
 __m256i vector_b_epi32 = _mm256_cvtepu8_epi32(vector_b_epi8);
 __m256i vector_c_epi32 = _mm256_cvtepu8_epi32(vector_c_epi8);

 __m256i temp1_epi32 = _mm256_sign_epi32(_mm256_set1_epi32(1), _mm256_sub_epi32(vector_c_epi32, vector_a_epi32));
 __m256i temp2_epi32 = _mm256_sign_epi32(_mm256_set1_epi32(1), _mm256_sub_epi32(vector_c_epi32, vector_b_epi32));


 eo_idx_epi32 = _mm256_add_epi32(eo_idx_epi32, temp1_epi32);
 eo_idx_epi32 = _mm256_add_epi32(eo_idx_epi32, temp2_epi32);

 __m256i v_cat_epi32 = _mm256_permutevar8x32_epi32(vector_sao_eo_idx_to_eo_category_epi32, eo_idx_epi32);
 return v_cat_epi32;
}


// Mapping of edge_idx values to eo-classes.
static __m256i sao_calc_eo_cat_6_pixels_avx2(__m128i vector_a_epi8, __m128i vector_b_epi8, __m128i vector_c_epi8)
{
 // Mapping relationships between a, b and c to eo_idx.
 __m256i vector_sao_eo_idx_to_eo_category_epi32 = _mm256_setr_epi32(1, 2, 0, 3, 4, 0, 0, 0);
 __m256i use_only_six = _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, 0, 0);

 __m256i eo_idx_epi32 = _mm256_set1_epi32(2);
 __m256i vector_a_epi32 = _mm256_cvtepu8_epi32(vector_a_epi8);
 __m256i vector_b_epi32 = _mm256_cvtepu8_epi32(vector_b_epi8);
 __m256i vector_c_epi32 = _mm256_cvtepu8_epi32(vector_c_epi8);

 __m256i temp1_epi32 = _mm256_sign_epi32(_mm256_set1_epi32(1), _mm256_sub_epi32(vector_c_epi32, vector_a_epi32));
 __m256i temp2_epi32 = _mm256_sign_epi32(_mm256_set1_epi32(1), _mm256_sub_epi32(vector_c_epi32, vector_b_epi32));

 eo_idx_epi32 = _mm256_add_epi32(eo_idx_epi32, temp1_epi32);
 eo_idx_epi32 = _mm256_add_epi32(eo_idx_epi32, temp2_epi32);

 __m256i v_cat_epi32 = _mm256_permutevar8x32_epi32(vector_sao_eo_idx_to_eo_category_epi32, eo_idx_epi32);
 v_cat_epi32 = _mm256_and_si256(v_cat_epi32, use_only_six);
 return v_cat_epi32;
}


static int sao_edge_ddistortion_avx2(const kvz_pixel *orig_data,
 const kvz_pixel *rec_data,
 int block_width,
 int block_height,
 int eo_class,
 int offsets[NUM_SAO_EDGE_CATEGORIES])
{
 int y, x;
 int sum = 0;
 vector2d_t a_ofs = g_sao_edge_offsets[eo_class][0];
 vector2d_t b_ofs = g_sao_edge_offsets[eo_class][1];

 __m256i offsets_epi32 = _mm256_setr_epi32(offsets[0], offsets[1], offsets[2], offsets[3], offsets[4], 0, 0, 0);
 __m256i tmp_diff_epi32;
 __m256i tmp_sum_epi32 = _mm256_setzero_si256();
 __m256i tmp_offset_epi32;
 __m256i tmp1_vec_epi32;
 __m256i tmp2_vec_epi32;
 __m256i zeros_epi32 = _mm256_setzero_si256();
 __m256i offset_zeros_epi32;

 for (y = 1; y < block_height - 1; ++y) {
  for (x = 1; x < block_width - 8; x+=8) {
   const kvz_pixel *c_data = &rec_data[y * block_width + x];

   kvz_pixel c = c_data[0];
   
   __m128i vector_a_epi8 = _mm_loadl_epi64((__m128i*)&c_data[a_ofs.y * block_width + a_ofs.x]);
   __m128i vector_c_epi8 = _mm_loadl_epi64((__m128i*)&c_data[0]);
   __m128i vector_b_epi8 = _mm_loadl_epi64((__m128i*)&c_data[b_ofs.y * block_width + b_ofs.x]);


   __m256i v_cat_epi32 = sao_calc_eo_cat_avx2(vector_a_epi8, vector_b_epi8, vector_c_epi8);

   tmp_diff_epi32 = _mm256_load_si256((__m256i*)&orig_data[y * block_width + x] - c);

   tmp_offset_epi32 = _mm256_permutevar8x32_epi32(offsets_epi32, v_cat_epi32);

   offset_zeros_epi32 = _mm256_cmpeq_epi32(zeros_epi32, tmp_offset_epi32);
   

   // (diff - offset) * (diff - offset)
   tmp1_vec_epi32 = _mm256_mul_epi32(_mm256_sub_epi32(tmp_diff_epi32, tmp_offset_epi32), _mm256_sub_epi32(tmp_diff_epi32, tmp_offset_epi32));

   // diff * diff
   tmp2_vec_epi32 = _mm256_mul_epi32(tmp_diff_epi32, tmp_diff_epi32);

   // Offset is applied to reconstruction, so it is subtracted from diff.
   // sum += (diff - offset) * (diff - offset) - diff * diff;

   tmp_sum_epi32 = _mm256_add_epi32(tmp_sum_epi32, _mm256_andnot_si256(offset_zeros_epi32, _mm256_sub_epi32(tmp1_vec_epi32, tmp2_vec_epi32)));
  }

  // Load the last 6 pixels to use

  const kvz_pixel *c_data = &rec_data[y * block_width + x];
  const kvz_pixel *c_data2 = &rec_data[y * block_width + x +2];
  const kvz_pixel *c_data4 = &rec_data[y * block_width + x + 4];

  kvz_pixel c = c_data[0];

  __m128i vector_a_epi8 = _mm_setr_epi16(c_data[a_ofs.y * block_width + a_ofs.x], c_data2[a_ofs.y * block_width + a_ofs.x], c_data4[a_ofs.y * block_width + a_ofs.x], 0, 0, 0, 0, 0);
  __m128i vector_c_epi8 = _mm_setr_epi16(c_data[0], c_data2[0], c_data4[0], 0, 0, 0, 0, 0);
  __m128i vector_b_epi8 = _mm_setr_epi16(c_data[b_ofs.y * block_width + b_ofs.x], c_data2[b_ofs.y * block_width + b_ofs.x], c_data4[b_ofs.y * block_width + b_ofs.x], 0, 0, 0, 0, 0);

  __m256i v_cat_epi32 = sao_calc_eo_cat_6_pixels_avx2(vector_a_epi8, vector_b_epi8, vector_c_epi8);
  
  tmp_diff_epi32 = _mm256_setr_epi32(orig_data[y * block_width + x] - c, orig_data[y * block_width + x + 1] - c, orig_data[y * block_width + x + 2] - c, orig_data[y * block_width + x + 3] - c, orig_data[y * block_width + x + 3] - c, orig_data[y * block_width + x + 4] - c, 0, 0);

  tmp_offset_epi32 = _mm256_permutevar8x32_epi32(offsets_epi32, v_cat_epi32);

  offset_zeros_epi32 = _mm256_cmpeq_epi32(zeros_epi32, tmp_offset_epi32);

  // (diff - offset) * (diff - offset)
  tmp1_vec_epi32 = _mm256_mul_epi32(_mm256_sub_epi32(tmp_diff_epi32, tmp_offset_epi32), _mm256_sub_epi32(tmp_diff_epi32, tmp_offset_epi32));

  // diff * diff
  tmp2_vec_epi32 = _mm256_mul_epi32(tmp_diff_epi32, tmp_diff_epi32);

  // Offset is applied to reconstruction, so it is subtracted from diff.
  // sum += (diff - offset) * (diff - offset) - diff * diff;

  tmp_sum_epi32 = _mm256_add_epi32(tmp_sum_epi32, _mm256_andnot_si256(offset_zeros_epi32, _mm256_sub_epi32(tmp1_vec_epi32, tmp2_vec_epi32)));

  tmp_sum_epi32 = _mm256_hadd_epi32(tmp_sum_epi32, tmp_sum_epi32);
  tmp_sum_epi32 = _mm256_hadd_epi32(tmp_sum_epi32, tmp_sum_epi32);
  int* pointer = (int*)&tmp_sum_epi32;

  sum += (pointer[0] + pointer[4]);
  
 }

 return sum;
}

/**
* \param orig_data  Original pixel data. 64x64 for luma, 32x32 for chroma.
* \param rec_data  Reconstructed pixel data. 64x64 for luma, 32x32 for chroma.
* \param dir_offsets
* \param is_chroma  0 for luma, 1 for chroma. Indicates
*/
static void calc_sao_edge_dir_avx2(const kvz_pixel *orig_data,
 const kvz_pixel *rec_data,
 int eo_class,
 int block_width,
 int block_height,
 int cat_sum_cnt[2][NUM_SAO_EDGE_CATEGORIES])
{
 int y, x;
 vector2d_t a_ofs = g_sao_edge_offsets[eo_class][0];
 vector2d_t b_ofs = g_sao_edge_offsets[eo_class][1];
 // Arrays orig_data and rec_data are quarter size for chroma.

 // Don't sample the edge pixels because this function doesn't have access to
 // their neighbours.

 __m256i zeros_epi32 = _mm256_setzero_si256();
 __m256i ones_epi32 = _mm256_set1_epi32(1);
 __m256i twos_epi32 = _mm256_set1_epi32(2);
 __m256i threes_epi32 = _mm256_set1_epi32(3);
 __m256i fours_epi32 = _mm256_set1_epi32(4);

 __m256i tmp_zero_values_epi32 = _mm256_setzero_si256();
 __m256i tmp_one_values_epi32 = _mm256_setzero_si256();
 __m256i tmp_two_values_epi32 = _mm256_setzero_si256();
 __m256i tmp_three_values_epi32 = _mm256_setzero_si256();
 __m256i tmp_four_values_epi32 = _mm256_setzero_si256();



 __m256i temp_epi32 = _mm256_setzero_si256();
 __m256i temp_mem_epi32 = _mm256_setzero_si256();

 for (y = 1; y < block_height - 1; ++y) {
  for (x = 1; x < block_width - 8; x += 8) {
   const kvz_pixel *c_data = &rec_data[y * block_width + x];

   kvz_pixel c = c_data[0];

   __m128i vector_a_epi8 = _mm_loadl_epi64((__m128i*)&c_data[a_ofs.y * block_width + a_ofs.x]);
   __m128i vector_c_epi8 = _mm_loadl_epi64((__m128i*)&c);
   __m128i vector_b_epi8 = _mm_loadl_epi64((__m128i*)&c_data[b_ofs.y * block_width + b_ofs.x]);


   __m256i v_cat_epi32 = sao_calc_eo_cat_avx2(vector_a_epi8, vector_b_epi8, vector_c_epi8);


   // Check wich values are right for specific cat amount.
   // It's done for every single value that cat could get {1, 2, 0, 3, 4}

   //--------------------------------------------------------------------------
   __m256i mask_epi32 = _mm256_cmpeq_epi32(zeros_epi32, v_cat_epi32);
   int temp_cnt = __popcnt(_mm256_movemask_epi8(mask_epi32))/4;
   cat_sum_cnt[1][0] += temp_cnt;
   temp_mem_epi32 = _mm256_load_si256((__m256i*)&orig_data[y * block_width + x] - c);
   temp_epi32 = _mm256_and_si256(mask_epi32, temp_mem_epi32);
   tmp_zero_values_epi32 = _mm256_add_epi32(tmp_zero_values_epi32, temp_epi32);
   //--------------------------------------------------------------------------

   mask_epi32 = _mm256_cmpeq_epi32(ones_epi32, v_cat_epi32);
   temp_cnt = __popcnt(_mm256_movemask_epi8(mask_epi32)) / 4;
   cat_sum_cnt[1][1] += temp_cnt;
   temp_epi32 = _mm256_and_si256(mask_epi32, temp_mem_epi32);
   tmp_one_values_epi32 = _mm256_add_epi32(tmp_zero_values_epi32, temp_epi32);
   //--------------------------------------------------------------------------

   mask_epi32 = _mm256_cmpeq_epi32(twos_epi32, v_cat_epi32);
   temp_cnt = __popcnt(_mm256_movemask_epi8(mask_epi32)) / 4;
   cat_sum_cnt[1][2] += temp_cnt;
   temp_epi32 = _mm256_and_si256(mask_epi32, temp_mem_epi32);
   tmp_two_values_epi32 = _mm256_add_epi32(tmp_zero_values_epi32, temp_epi32);
   //--------------------------------------------------------------------------

   mask_epi32 = _mm256_cmpeq_epi32(threes_epi32, v_cat_epi32);
   temp_cnt = __popcnt(_mm256_movemask_epi8(mask_epi32)) / 4;
   cat_sum_cnt[1][3] += temp_cnt;
   temp_epi32 = _mm256_and_si256(mask_epi32, temp_mem_epi32);
   tmp_three_values_epi32 = _mm256_add_epi32(tmp_zero_values_epi32, temp_epi32);
   //--------------------------------------------------------------------------

   mask_epi32 = _mm256_cmpeq_epi32(fours_epi32, v_cat_epi32);
   temp_cnt = __popcnt(_mm256_movemask_epi8(mask_epi32)) / 4;
   cat_sum_cnt[1][4] += temp_cnt;
   temp_epi32 = _mm256_and_si256(mask_epi32, temp_mem_epi32);
   tmp_four_values_epi32 = _mm256_add_epi32(tmp_zero_values_epi32, temp_epi32);
   //--------------------------------------------------------------------------


  }
  temp_epi32 = _mm256_hadd_epi32(tmp_zero_values_epi32, tmp_one_values_epi32);
  temp_mem_epi32 = _mm256_hadd_epi32(tmp_two_values_epi32, tmp_three_values_epi32);

  temp_epi32 = _mm256_hadd_epi32(temp_epi32, temp_mem_epi32);

  int*temp = (int*)&temp_epi32;

  cat_sum_cnt[0][0] += (temp[0] + temp[4]);
  cat_sum_cnt[0][1] += (temp[1] + temp[5]);
  cat_sum_cnt[0][2] += (temp[2] + temp[6]);
  cat_sum_cnt[0][3] += (temp[3] + temp[7]);

  tmp_four_values_epi32 = _mm256_hadd_epi32(tmp_four_values_epi32, tmp_four_values_epi32);
  tmp_four_values_epi32 = _mm256_hadd_epi32(tmp_four_values_epi32, tmp_four_values_epi32);
  tmp_four_values_epi32 = _mm256_hadd_epi32(tmp_four_values_epi32, tmp_four_values_epi32);

  temp = (int*)&tmp_four_values_epi32;
  cat_sum_cnt[0][4] += (temp[0] + temp[4]);


  // Load the last 6 pixels to use

  const kvz_pixel *c_data = &rec_data[y * block_width + x];
  const kvz_pixel *c_data2 = &rec_data[y * block_width + x + 2];
  const kvz_pixel *c_data4 = &rec_data[y * block_width + x + 4];

  kvz_pixel c = c_data[0];

  __m128i vector_a_epi8 = _mm_setr_epi16(c_data[a_ofs.y * block_width + a_ofs.x], c_data2[a_ofs.y * block_width + a_ofs.x], c_data4[a_ofs.y * block_width + a_ofs.x], 0, 0, 0, 0, 0);
  __m128i vector_c_epi8 = _mm_setr_epi16(c_data[0], c_data2[0], c_data4[0], 0, 0, 0, 0, 0);
  __m128i vector_b_epi8 = _mm_setr_epi16(c_data[b_ofs.y * block_width + b_ofs.x], c_data2[b_ofs.y * block_width + b_ofs.x], c_data4[b_ofs.y * block_width + b_ofs.x], 0, 0, 0, 0, 0);

  __m256i v_cat_epi32 = sao_calc_eo_cat_6_pixels_avx2(vector_a_epi8, vector_b_epi8, vector_c_epi8);

  __m256i temp_mem_epi32 = _mm256_setr_epi32(orig_data[y * block_width + x] - c, orig_data[y * block_width + x + 1] - c, orig_data[y * block_width + x + 2] - c, orig_data[y * block_width + x + 3] - c, orig_data[y * block_width + x + 3] - c, orig_data[y * block_width + x + 4] - c, 0, 0);

  // Check wich values are right for specific cat amount.
  // It's done for every single value that cat could get {1, 2, 0, 3, 4}
  //--------------------------------------------------------------------------
  __m256i mask_epi32 = _mm256_cmpeq_epi32(zeros_epi32, v_cat_epi32);
  int temp_cnt = __popcnt(_mm256_movemask_epi8(mask_epi32)) / 4 - 2;
  cat_sum_cnt[1][0] += temp_cnt;
  temp_epi32 = _mm256_and_si256(mask_epi32, temp_mem_epi32);
  tmp_zero_values_epi32 = _mm256_add_epi32(tmp_zero_values_epi32, temp_epi32);
  //--------------------------------------------------------------------------

  mask_epi32 = _mm256_cmpeq_epi32(ones_epi32, v_cat_epi32);
  temp_cnt = __popcnt(_mm256_movemask_epi8(mask_epi32)) / 4;
  cat_sum_cnt[1][1] += temp_cnt;
  temp_epi32 = _mm256_and_si256(mask_epi32, temp_mem_epi32);
  tmp_one_values_epi32 = _mm256_add_epi32(tmp_zero_values_epi32, temp_epi32);
  //--------------------------------------------------------------------------

  mask_epi32 = _mm256_cmpeq_epi32(twos_epi32, v_cat_epi32);
  temp_cnt = __popcnt(_mm256_movemask_epi8(mask_epi32)) / 4;
  cat_sum_cnt[1][2] += temp_cnt;
  temp_epi32 = _mm256_and_si256(mask_epi32, temp_mem_epi32);
  tmp_two_values_epi32 = _mm256_add_epi32(tmp_zero_values_epi32, temp_epi32);
  //--------------------------------------------------------------------------

  mask_epi32 = _mm256_cmpeq_epi32(threes_epi32, v_cat_epi32);
  temp_cnt = __popcnt(_mm256_movemask_epi8(mask_epi32)) / 4;
  cat_sum_cnt[1][3] += temp_cnt;
  temp_epi32 = _mm256_and_si256(mask_epi32, temp_mem_epi32);
  tmp_three_values_epi32 = _mm256_add_epi32(tmp_zero_values_epi32, temp_epi32);
  //--------------------------------------------------------------------------

  mask_epi32 = _mm256_cmpeq_epi32(fours_epi32, v_cat_epi32);
  temp_cnt = __popcnt(_mm256_movemask_epi8(mask_epi32)) / 4;
  cat_sum_cnt[1][4] += temp_cnt;
  temp_epi32 = _mm256_and_si256(mask_epi32, temp_mem_epi32);
  tmp_four_values_epi32 = _mm256_add_epi32(tmp_zero_values_epi32, temp_epi32);
  //--------------------------------------------------------------------------

  temp_epi32 = _mm256_hadd_epi32(tmp_zero_values_epi32, tmp_one_values_epi32);

  int*remove = (int*)&temp_epi32;

  temp_mem_epi32 = _mm256_hadd_epi32(tmp_two_values_epi32, tmp_three_values_epi32);

  temp_mem_epi32 = _mm256_hadd_epi32(temp_epi32, temp_mem_epi32);

  temp = (int*)&temp_mem_epi32;

  cat_sum_cnt[0][0] += (temp[0] + temp[4] - remove[5]);
  cat_sum_cnt[0][1] += (temp[1] + temp[5]);
  cat_sum_cnt[0][2] += (temp[2] + temp[6]);
  cat_sum_cnt[0][3] += (temp[3] + temp[7]);

  tmp_four_values_epi32 = _mm256_hadd_epi32(tmp_four_values_epi32, tmp_four_values_epi32);
  tmp_four_values_epi32 = _mm256_hadd_epi32(tmp_four_values_epi32, tmp_four_values_epi32);
  tmp_four_values_epi32 = _mm256_hadd_epi32(tmp_four_values_epi32, tmp_four_values_epi32);

  temp = (int*)&tmp_four_values_epi32;
  cat_sum_cnt[0][4] += (temp[0] + temp[4]);

 }
}


static void sao_reconstruct_color_avx(const encoder_control_t * const encoder,
 const kvz_pixel *rec_data,
 kvz_pixel *new_rec_data,
 const sao_info_t *sao,
 int stride,
 int new_stride,
 int block_width,
 int block_height,
 color_t color_i)
{
 // Arrays orig_data and rec_data are quarter size for chroma.
 int offset_v = color_i == COLOR_V ? 5 : 0;

 if (sao->type == SAO_TYPE_BAND) {
  int offsets[1 << KVZ_BIT_DEPTH];
  kvz_calc_sao_offset_array(encoder, sao, offsets, color_i);
  for (int y = 0; y < block_height; ++y) {
   for (int x = 0; x < block_width; ++x) {
    new_rec_data[y * new_stride + x] = offsets[rec_data[y * stride + x]];
   }
  }
 }
 else {
  // Don't sample the edge pixels because this function doesn't have access to
  // their neighbours.
  for (int y = 0; y < block_height; ++y) {
   for (int x = 0; x < block_width; x += 8) {

    for (int i = 0; i < 8; ++i) {

     int test = x + i;
     vector2d_t a_ofs = g_sao_edge_offsets[sao->eo_class][0];
     vector2d_t b_ofs = g_sao_edge_offsets[sao->eo_class][1];
     const kvz_pixel *c_data = &rec_data[y * stride + test];
     kvz_pixel *new_data = &new_rec_data[y * new_stride + test];
     kvz_pixel a = c_data[a_ofs.y * stride + a_ofs.x];
     kvz_pixel c = c_data[0];
     kvz_pixel b = c_data[b_ofs.y * stride + b_ofs.x];

     int eo_cat = sao_calc_eo_cat(a, b, c);

     new_data[0] = (kvz_pixel)CLIP(0, (1 << KVZ_BIT_DEPTH) - 1, c_data[0] + sao->offsets[eo_cat + offset_v]);

    }

   }
  }
 }
}

static void sao_reconstruct_color_avx2(const encoder_control_t * const encoder,
 const kvz_pixel *rec_data,
 kvz_pixel *new_rec_data,
 const sao_info_t *sao,
 int stride,
 int new_stride,
 int block_width,
 int block_height,
 color_t color_i)
{
 // Arrays orig_data and rec_data are quarter size for chroma.
 int offset_v = color_i == COLOR_V ? 5 : 0;
 

 if (sao->type == SAO_TYPE_BAND) {
  int offsets[1 << KVZ_BIT_DEPTH];
  kvz_calc_sao_offset_array(encoder, sao, offsets, color_i);
  for (int y = 0; y < block_height; ++y) {
   for (int x = 0; x < block_width; ++x) {
    new_rec_data[y * new_stride + x] = offsets[rec_data[y * stride + x]];
   }
  }
 }
 else {

  // Don't sample the edge pixels because this function doesn't have access to
  // their neighbours.

  __m256i offset_v_epi32 = _mm256_set1_epi32(offset_v);

  vector2d_t a_ofs = g_sao_edge_offsets[sao->eo_class][0];
  vector2d_t b_ofs = g_sao_edge_offsets[sao->eo_class][1];

  for (int y = 0; y < block_height; ++y) {
   int test = 0;

   for (int x = 0; x < block_width - 8; x+=8) {


    const kvz_pixel *c_data = &rec_data[y * stride + x];

    __m128i vector_a_epi8 = _mm_loadl_epi64((__m128i*)&c_data[a_ofs.y * stride + a_ofs.x]);
    __m128i vector_c_epi8 = _mm_loadl_epi64((__m128i*)&c_data[0]);
    __m128i vector_b_epi8 = _mm_loadl_epi64((__m128i*)&c_data[b_ofs.y * stride + b_ofs.x]);


    __m256i v_cat_epi32 = sao_calc_eo_cat_avx2(vector_a_epi8, vector_b_epi8, vector_c_epi8);


    v_cat_epi32 = _mm256_add_epi32(v_cat_epi32, offset_v_epi32);

    __m256i vector_c_data0_epi32 = _mm256_cvtepu8_epi32(vector_c_epi8);


    int*temp = (int*)&v_cat_epi32;
    __m256i vector_sao_offsets_epi32 = _mm256_set_epi32(sao->offsets[temp[7]], sao->offsets[temp[6]], sao->offsets[temp[5]], sao->offsets[temp[4]], sao->offsets[temp[3]], sao->offsets[temp[2]], sao->offsets[temp[1]], sao->offsets[temp[0]]);
    vector_sao_offsets_epi32 = _mm256_add_epi32(vector_sao_offsets_epi32, vector_c_data0_epi32);

    __m256i temp_epi16 = _mm256_packus_epi32(vector_sao_offsets_epi32, vector_sao_offsets_epi32);
    __m256i temp_epi8 = _mm256_packus_epi16(temp_epi16, temp_epi16);



    int*temp2 = (int*)&vector_sao_offsets_epi32;
    
    for (int i = 0; i < 8; ++i) {

     const kvz_pixel *c_data = &rec_data[y * stride + x + i];

     kvz_pixel *new_data = &new_rec_data[y * new_stride + x + i];

     //printf("%d ", c_data[0] + sao->offsets[temp[i]]);
     //printf("%d \n", temp2[i]);


     new_data[0] = (kvz_pixel)CLIP(0, (1 << KVZ_BIT_DEPTH) - 1, temp2[i]);//c_data[0] + sao->offsets[temp[i]]);
     test = x;
    }
    //Low = 0
    //High = (1 << KVZ_BIT_DEPTH)
    //Value = c_data[0] + sao->offsets[eo_cat + offset_v]
    //new_data[0] = (kvz_pixel)CLIP(0, (1 << KVZ_BIT_DEPTH) - 1, c_data[0] + sao->offsets[eo_cat + offset_v]);
   }


   for (int i = 0; i < (block_width - test); ++i) {

    const kvz_pixel *c_data = &rec_data[y * stride + test + i];

    kvz_pixel *new_data = &new_rec_data[y * new_stride + test + i];
    kvz_pixel a = c_data[a_ofs.y * stride + a_ofs.x];
    kvz_pixel c = c_data[0];
    kvz_pixel b = c_data[b_ofs.y * stride + b_ofs.x];


    int eo_cat = sao_calc_eo_cat(a, b, c);

    new_data[0] = (kvz_pixel)CLIP(0, (1 << KVZ_BIT_DEPTH) - 1, c_data[0] + sao->offsets[eo_cat + offset_v]);

   }


  }
 }
}


static int sao_band_ddistortion_avx2(const encoder_state_t * const state,
 const kvz_pixel *orig_data,
 const kvz_pixel *rec_data,
 int block_width,
 int block_height,
 int band_pos,
 int sao_bands[4])
{
 int y, x;
 int shift = state->encoder_control->bitdepth - 5;
 int sum = 0;

 for (y = 0; y < block_height; ++y) {
  for (x = 0; x < block_width; ++x) {
   int band = (rec_data[y * block_width + x] >> shift) - band_pos;
   int offset = 0;
   if (band >= 0 && band < 4) {
    offset = sao_bands[band];
   }
   if (offset != 0) {
    int diff = orig_data[y * block_width + x] - rec_data[y * block_width + x];
    // Offset is applied to reconstruction, so it is subtracted from diff.
    sum += (diff - offset) * (diff - offset) - diff * diff;
   }
  }
 }

 return sum;
}

#endif //COMPILE_INTEL_AVX2

int kvz_strategy_register_sao_avx2(void* opaque, uint8_t bitdepth)
{
  bool success = true;
#if COMPILE_INTEL_AVX2
  if (bitdepth == 8) {
    success &= kvz_strategyselector_register(opaque, "sao_edge_ddistortion", "avx2", 40, &sao_edge_ddistortion_avx2);
    success &= kvz_strategyselector_register(opaque, "calc_sao_edge_dir", "avx2", 40, &calc_sao_edge_dir_avx2);
    success &= kvz_strategyselector_register(opaque, "sao_reconstruct_color", "avx2", 40, &sao_reconstruct_color_avx2);
    success &= kvz_strategyselector_register(opaque, "sao_band_ddistortion", "avx2", 40, &sao_band_ddistortion_avx2);
  }
#endif //COMPILE_INTEL_AVX2
  return success;
}
