#ifndef STRATEGIES_SAO_H_
#define STRATEGIES_SAO_H_
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

/**
 * \ingroup Optimization
 * \file
 * Interface for sao functions.
 */

#include "encoder.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "sao.h"


// Declare function pointers.
typedef int (sao_edge_ddistortion_func)(const kvz_pixel *orig_data, const kvz_pixel *rec_data,
  int block_width, int block_height,
  int eo_class, int offsets[NUM_SAO_EDGE_CATEGORIES]);

typedef void (calc_sao_edge_dir_func)(const kvz_pixel *orig_data, const kvz_pixel *rec_data,
  int eo_class, int block_width, int block_height,
  int cat_sum_cnt[2][NUM_SAO_EDGE_CATEGORIES]);

typedef void (sao_reconstruct_color_func)(const encoder_control_t * const encoder,
  const kvz_pixel *rec_data, kvz_pixel *new_rec_data,
  const sao_info_t *sao,
  int stride, int new_stride,
  int block_width, int block_height,
  color_t color_i);

typedef int (sao_band_ddistortion_func)(const encoder_state_t * const state, const kvz_pixel *orig_data, const kvz_pixel *rec_data,
  int block_width, int block_height,
  int band_pos, const int sao_bands[4]);

// Declare function pointers.
extern sao_edge_ddistortion_func * kvz_sao_edge_ddistortion;
extern calc_sao_edge_dir_func * kvz_calc_sao_edge_dir;
extern sao_reconstruct_color_func * kvz_sao_reconstruct_color;
extern sao_band_ddistortion_func * kvz_sao_band_ddistortion;

int kvz_strategy_register_sao(void* opaque, uint8_t bitdepth);


#define STRATEGIES_SAO_EXPORTS \
  {"sao_edge_ddistortion", (void**) &kvz_sao_edge_ddistortion}, \
  {"calc_sao_edge_dir", (void**) &kvz_calc_sao_edge_dir}, \
  {"sao_reconstruct_color", (void**) &kvz_sao_reconstruct_color}, \
  {"sao_band_ddistortion", (void**) &kvz_sao_band_ddistortion}, \



#endif //STRATEGIES_SAO_H_
