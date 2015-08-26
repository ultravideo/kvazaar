#ifndef SAO_H_
#define SAO_H_
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
 * \brief Coding Unit (CU) and picture data related functions.
 */

#include "checkpoint.h"
#include "global.h"
#include "videoframe.h"
#include "encoder.h"
#include "encoderstate.h"
#include "math.h"

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

#define CHECKPOINT_SAO_INFO(prefix_str, sao) CHECKPOINT(prefix_str " type=%d eo_class=%d ddistortion=%d " \
  "merge_left_flag=%d merge_up_flag=%d band_position=%d " \
  "offsets[0]=%d offsets[1]=%d offsets[2]=%d offsets[3]=%d offsets[4]=%d", \
  (sao).type, (sao).eo_class, (sao).ddistortion, \
  (sao).merge_left_flag, (sao).merge_up_flag, (sao).band_position[0], \
  (sao).offsets[0], (sao).offsets[1], (sao).offsets[2], (sao).offsets[3], (sao).offsets[4])


void kvz_init_sao_info(sao_info_t *sao);
void kvz_sao_search_chroma(const encoder_state_t * state, const videoframe_t *frame, unsigned x_ctb, unsigned y_ctb, sao_info_t *sao, sao_info_t *sao_top, sao_info_t *sao_left, int32_t merge_cost[3]);
void kvz_sao_search_luma(const encoder_state_t * state, const videoframe_t *frame, unsigned x_ctb, unsigned y_ctb, sao_info_t *sao, sao_info_t *sao_top, sao_info_t *sao_left, int32_t merge_cost[3]);
void kvz_sao_reconstruct(const encoder_control_t * encoder, videoframe_t *frame, const kvz_pixel *old_rec,
                     unsigned x_ctb, unsigned y_ctb,
                     const sao_info_t *sao, color_t color_i);
void kvz_sao_reconstruct_frame(encoder_state_t *state);

#endif
