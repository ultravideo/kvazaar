#ifndef SAO_H_
#define SAO_H_
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
 * \ingroup Reconstruction
 * \file
 * Sample Adaptive Offset filter.
 */

#include "checkpoint.h"
#include "cu.h"
#include "encoder.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "videoframe.h"


typedef enum { SAO_TYPE_NONE = 0, SAO_TYPE_BAND, SAO_TYPE_EDGE } sao_type;
typedef enum { SAO_EO0 = 0, SAO_EO1, SAO_EO2, SAO_EO3, SAO_NUM_EO } sao_eo_class;
typedef enum { SAO_EO_CAT0 = 0, SAO_EO_CAT1, SAO_EO_CAT2, SAO_EO_CAT3, SAO_EO_CAT4, NUM_SAO_EDGE_CATEGORIES } sao_eo_cat;


typedef struct sao_info_t {
  sao_type type;
  sao_eo_class eo_class;
  int ddistortion;
  int merge_left_flag;
  int merge_up_flag;
  int band_position[2];
  int offsets[NUM_SAO_EDGE_CATEGORIES * 2];
} sao_info_t;


// Offsets of a and b in relation to c.
// dir_offset[dir][a or b]
// |       |   a   | a     |     a |
// | a c b |   c   |   c   |   c   |
// |       |   b   |     b | b     |
static const vector2d_t g_sao_edge_offsets[SAO_NUM_EO][2] = {
  { { -1, 0 }, { 1, 0 } },
  { { 0, -1 }, { 0, 1 } },
  { { -1, -1 }, { 1, 1 } },
  { { 1, -1 }, { -1, 1 } }
};


#define CHECKPOINT_SAO_INFO(prefix_str, sao) CHECKPOINT(prefix_str " type=%d eo_class=%d ddistortion=%d " \
  "merge_left_flag=%d merge_up_flag=%d band_position=%d " \
  "offsets[0]=%d offsets[1]=%d offsets[2]=%d offsets[3]=%d offsets[4]=%d", \
  (sao).type, (sao).eo_class, (sao).ddistortion, \
  (sao).merge_left_flag, (sao).merge_up_flag, (sao).band_position[0], \
  (sao).offsets[0], (sao).offsets[1], (sao).offsets[2], (sao).offsets[3], (sao).offsets[4])


void kvz_sao_reconstruct(const encoder_state_t *state,
                         const kvz_pixel *buffer,
                         int stride,
                         int frame_x,
                         int frame_y,
                         int width,
                         int height,
                         const sao_info_t *sao,
                         color_t color);

void kvz_sao_search_lcu(const encoder_state_t* const state, int lcu_x, int lcu_y);
void kvz_calc_sao_offset_array(const encoder_control_t * const encoder, const sao_info_t *sao, int *offset, color_t color_i);

#endif
