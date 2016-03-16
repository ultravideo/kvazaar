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

#include <stdlib.h>

#include "sao-avx2.h"
#include "sao.h"
#include "strategyselector.h"
#include "rdo.h"
#include "strategies/strategies-common.h"

#if COMPILE_INTEL_AVX2
#include <immintrin.h>
// These optimizations are based heavily on sao-generic.c.
// Might be useful to check that if (when) this file
// is difficult to understand.


static INLINE __m256i load_6_offsets(const int* offsets){

  return _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((__m128i*) offsets)), _mm_loadl_epi64((__m128i*)&(offsets[4])), 1);
}

static INLINE __m128i load_6_pixels(const kvz_pixel* data){

  return _mm_insert_epi16(_mm_cvtsi32_si128(*(int32_t*)&(data[0])), *(int16_t*)&(data[4]), 2);
}

static INLINE __m256i load_5_offsets(const int* offsets){

  return _mm256_inserti128_si256(_mm256_castsi128_si256(_mm_loadu_si128((__m128i*) offsets)), _mm_insert_epi32(_mm_setzero_si128(), offsets[4], 0), 1);
}


// Mapping of edge_idx values to eo-classes.
static int sao_calc_eo_cat(kvz_pixel a, kvz_pixel b, kvz_pixel c)
{
  // Mapping relationships between a, b and c to eo_idx.
  static const int sao_eo_idx_to_eo_category[] = { 1, 2, 0, 3, 4 };

  int eo_idx = 2 + SIGN3((int)c - (int)a) + SIGN3((int)c - (int)b);

  return sao_eo_idx_to_eo_category[eo_idx];
}


static __m128i sao_calc_eo_cat_avx2(__m128i* a, __m128i* b, __m128i* c)
{
  __m128i v_eo_idx = _mm_set1_epi16(2);
  __m128i v_a = _mm_cvtepu8_epi16(*a);
  __m128i v_c = _mm_cvtepu8_epi16(*c);
  __m128i v_b = _mm_cvtepu8_epi16(*b);
  
  __m128i temp_a = _mm_sign_epi16(_mm_set1_epi16(1), _mm_sub_epi16(v_c, v_a));
  __m128i temp_b = _mm_sign_epi16(_mm_set1_epi16(1), _mm_sub_epi16(v_c, v_b));
  v_eo_idx = _mm_add_epi16(v_eo_idx, temp_a);
  v_eo_idx = _mm_add_epi16(v_eo_idx, temp_b);
  
  v_eo_idx = _mm_packus_epi16(v_eo_idx, v_eo_idx);
  __m128i v_cat_lookup = _mm_setr_epi8(1,2,0,3,4,0,0,0,0,0,0,0,0,0,0,0);
  __m128i v_cat = _mm_shuffle_epi8(v_cat_lookup, v_eo_idx);


  return v_cat;
}


int kvz_sao_edge_ddistortion_avx2(const kvz_pixel *orig_data, const kvz_pixel *rec_data,
                         int block_width, int block_height,
                         int eo_class, int offsets[NUM_SAO_EDGE_CATEGORIES])
{
  int y, x;
  int sum = 0;
  vector2d_t a_ofs = g_sao_edge_offsets[eo_class][0];
  vector2d_t b_ofs = g_sao_edge_offsets[eo_class][1];

  __m256i v_accum = { 0 };

  for (y = 1; y < block_height - 1; ++y) {

    for (x = 1; x < block_width - 8; x+=8) {
      const kvz_pixel *c_data = &rec_data[y * block_width + x];

      __m128i v_c_data = _mm_loadl_epi64((__m128i*)c_data);
      __m128i v_a = _mm_loadl_epi64((__m128i*)(&c_data[a_ofs.y * block_width + a_ofs.x]));
      __m128i v_c = v_c_data;
      __m128i v_b = _mm_loadl_epi64((__m128i*)(&c_data[b_ofs.y * block_width + b_ofs.x]));

      __m256i v_cat = _mm256_cvtepu8_epi32(sao_calc_eo_cat_avx2(&v_a, &v_b, &v_c));

      __m256i v_offset = _mm256_loadu_si256((__m256i*) offsets);
      v_offset = _mm256_permutevar8x32_epi32(v_offset, v_cat);
   
      __m256i v_diff = _mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i*)&(orig_data[y * block_width + x])));
      v_diff = _mm256_sub_epi32(v_diff, _mm256_cvtepu8_epi32(v_c));
      __m256i v_diff_minus_offset = _mm256_sub_epi32(v_diff, v_offset);
      __m256i v_temp_sum = _mm256_sub_epi32(_mm256_mullo_epi32(v_diff_minus_offset, v_diff_minus_offset), _mm256_mullo_epi32(v_diff, v_diff));
      v_accum = _mm256_add_epi32(v_accum, v_temp_sum);
    }

    //Handle last 6 pixels separately to prevent reading over boundary
    const kvz_pixel *c_data = &rec_data[y * block_width + x];
    __m128i v_c_data = load_6_pixels(c_data);
    const kvz_pixel* a_ptr = &c_data[a_ofs.y * block_width + a_ofs.x];
    const kvz_pixel* b_ptr = &c_data[b_ofs.y * block_width + b_ofs.x];
    __m128i v_a = load_6_pixels(a_ptr);
    __m128i v_c = v_c_data;
    __m128i v_b = load_6_pixels(b_ptr);

    __m256i v_cat = _mm256_cvtepu8_epi32(sao_calc_eo_cat_avx2(&v_a, &v_b, &v_c));

    __m256i v_offset = load_6_offsets(offsets);
    v_offset = _mm256_permutevar8x32_epi32(v_offset, v_cat);
   
    const kvz_pixel* orig_ptr = &(orig_data[y * block_width + x]);
    __m256i v_diff = _mm256_cvtepu8_epi32(load_6_pixels(orig_ptr));
    v_diff = _mm256_sub_epi32(v_diff, _mm256_cvtepu8_epi32(v_c));

    __m256i v_diff_minus_offset = _mm256_sub_epi32(v_diff, v_offset);
    __m256i v_temp_sum = _mm256_sub_epi32(_mm256_mullo_epi32(v_diff_minus_offset, v_diff_minus_offset), _mm256_mullo_epi32(v_diff, v_diff));
    v_accum = _mm256_add_epi32(v_accum, v_temp_sum);
  }

  //Full horizontal sum
  v_accum = _mm256_add_epi32(v_accum, _mm256_castsi128_si256(_mm256_extracti128_si256(v_accum, 1)));
  v_accum = _mm256_add_epi32(v_accum, _mm256_shuffle_epi32(v_accum, KVZ_PERMUTE(2, 3, 0, 1)));
  v_accum = _mm256_add_epi32(v_accum, _mm256_shuffle_epi32(v_accum, KVZ_PERMUTE(1, 0, 1, 0)));
  sum += _mm_cvtsi128_si32(_mm256_castsi256_si128(v_accum));

  return sum;
}


/**
 * \param orig_data  Original pixel data. 64x64 for luma, 32x32 for chroma.
 * \param rec_data  Reconstructed pixel data. 64x64 for luma, 32x32 for chroma.
 * \param dir_offsets
 * \param is_chroma  0 for luma, 1 for chroma. Indicates
 */
void kvz_calc_sao_edge_dir_avx2(const kvz_pixel *orig_data, const kvz_pixel *rec_data,
                              int eo_class, int block_width, int block_height,
                              int cat_sum_cnt[2][NUM_SAO_EDGE_CATEGORIES])
{
  int y, x;
  vector2d_t a_ofs = g_sao_edge_offsets[eo_class][0];
  vector2d_t b_ofs = g_sao_edge_offsets[eo_class][1];
  // Arrays orig_data and rec_data are quarter size for chroma.

  // Don't sample the edge pixels because this function doesn't have access to
  // their neighbours.
  for (y = 1; y < block_height - 1; ++y) {
    for (x = 1; x < block_width - 1; ++x) {
      const kvz_pixel *c_data = &rec_data[y * block_width + x];
      kvz_pixel a = c_data[a_ofs.y * block_width + a_ofs.x];
      kvz_pixel c = c_data[0];
      kvz_pixel b = c_data[b_ofs.y * block_width + b_ofs.x];

      int eo_cat = sao_calc_eo_cat(a, b, c);

      cat_sum_cnt[0][eo_cat] += orig_data[y * block_width + x] - c;
      cat_sum_cnt[1][eo_cat] += 1;
    }
  }
}


void kvz_sao_reconstruct_color_avx2(const encoder_control_t * const encoder, 
                                  const kvz_pixel *rec_data, kvz_pixel *new_rec_data,
                                  const sao_info_t *sao,
                                  int stride, int new_stride,
                                  int block_width, int block_height,
                                  color_t color_i)
{
  int y, x;
  // Arrays orig_data and rec_data are quarter size for chroma.
  int offset_v = color_i == COLOR_V ? 5 : 0;

  if(sao->type == SAO_TYPE_BAND) {
    int offsets[1<<KVZ_BIT_DEPTH];
    kvz_calc_sao_offset_array(encoder, sao, offsets, color_i);
    for (y = 0; y < block_height; ++y) {
      for (x = 0; x < block_width; ++x) {
        new_rec_data[y * new_stride + x] = offsets[rec_data[y * stride + x]];
      }
    }
  } else {
    // Don't sample the edge pixels because this function doesn't have access to
    // their neighbours.
    for (y = 0; y < block_height; ++y) {
      for (x = 0; x < block_width; ++x) {
        vector2d_t a_ofs = g_sao_edge_offsets[sao->eo_class][0];
        vector2d_t b_ofs = g_sao_edge_offsets[sao->eo_class][1];
        const kvz_pixel *c_data = &rec_data[y * stride + x];
        kvz_pixel *new_data = &new_rec_data[y * new_stride + x];
        kvz_pixel a = c_data[a_ofs.y * stride + a_ofs.x];
        kvz_pixel c = c_data[0];
        kvz_pixel b = c_data[b_ofs.y * stride + b_ofs.x];

        int eo_cat = sao_calc_eo_cat(a, b, c);

        new_data[0] = (kvz_pixel)CLIP(0, (1 << KVZ_BIT_DEPTH) - 1, c_data[0] + sao->offsets[eo_cat + offset_v]);
      }
    }
  }
}


int kvz_sao_band_ddistortion_avx2(const encoder_state_t * const state, const kvz_pixel *orig_data, const kvz_pixel *rec_data,
                         int block_width, int block_height,
                         int band_pos, int sao_bands[4])
{
  int y, x;
  int shift = state->encoder_control->bitdepth-5;
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
    success &= kvz_strategyselector_register(opaque, "sao_edge_ddistortion", "avx2", 40, &kvz_sao_edge_ddistortion_avx2);
    success &= kvz_strategyselector_register(opaque, "calc_sao_edge_dir", "avx2", 40, &kvz_calc_sao_edge_dir_avx2);
    success &= kvz_strategyselector_register(opaque, "sao_reconstruct_color", "avx2", 40, &kvz_sao_reconstruct_color_avx2);
    success &= kvz_strategyselector_register(opaque, "sao_band_ddistortion", "avx2", 40, &kvz_sao_band_ddistortion_avx2);
  }
#endif //COMPILE_INTEL_AVX2
  return success;
}
