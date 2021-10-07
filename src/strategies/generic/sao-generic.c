/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "strategies/generic/sao-generic.h"
#include "strategies/generic/sao_shared_generics.h"

#include "cu.h"
#include "encoder.h"
#include "encoderstate.h"
#include "kvazaar.h"
#include "sao.h"
#include "strategyselector.h"


/**
 * \param orig_data  Original pixel data. 64x64 for luma, 32x32 for chroma.
 * \param rec_data  Reconstructed pixel data. 64x64 for luma, 32x32 for chroma.
 * \param dir_offsets
 * \param is_chroma  0 for luma, 1 for chroma. Indicates
 */
static void calc_sao_edge_dir_generic(const kvz_pixel *orig_data,
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


static void sao_reconstruct_color_generic(const encoder_control_t * const encoder,
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
    int offsets[1<<KVZ_BIT_DEPTH];
    kvz_calc_sao_offset_array(encoder, sao, offsets, color_i);
    for (int y = 0; y < block_height; ++y) {
      for (int x = 0; x < block_width; ++x) {
        new_rec_data[y * new_stride + x] = offsets[rec_data[y * stride + x]];
      }
    }
  } else {
    // Don't sample the edge pixels because this function doesn't have access to
    // their neighbours.
    for (int y = 0; y < block_height; ++y) {
      for (int x = 0; x < block_width; ++x) {
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



int kvz_strategy_register_sao_generic(void* opaque, uint8_t bitdepth)
{
  bool success = true;

  success &= kvz_strategyselector_register(opaque, "sao_edge_ddistortion", "generic", 0, &sao_edge_ddistortion_generic);
  success &= kvz_strategyselector_register(opaque, "calc_sao_edge_dir", "generic", 0, &calc_sao_edge_dir_generic);
  success &= kvz_strategyselector_register(opaque, "sao_reconstruct_color", "generic", 0, &sao_reconstruct_color_generic);
  success &= kvz_strategyselector_register(opaque, "sao_band_ddistortion", "generic", 0, &sao_band_ddistortion_generic);

  return success;
}
