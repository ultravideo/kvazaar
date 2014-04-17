#ifndef SAO_H_
#define SAO_H_
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
 * \brief Coding Unit (CU) and picture data related functions.
 */

#include "global.h"
#include "picture.h"
#include "encoder.h"
#include "math.h"

typedef enum { SAO_TYPE_NONE = 0, SAO_TYPE_BAND, SAO_TYPE_EDGE } sao_type;
typedef enum { SAO_EO0 = 0, SAO_EO1, SAO_EO2, SAO_EO3, SAO_NUM_EO } sao_eo_class;
typedef enum { SAO_EO_CAT0 = 0, SAO_EO_CAT1, SAO_EO_CAT2, SAO_EO_CAT3, SAO_EO_CAT4, NUM_SAO_EDGE_CATEGORIES } sao_eo_cat;


typedef struct sao_info_struct {
  sao_type type;
  sao_eo_class eo_class;
  int ddistortion;
  int merge_left_flag;
  int merge_up_flag;
  int band_position;
  int offsets[NUM_SAO_EDGE_CATEGORIES];
} sao_info;


void init_sao_info(sao_info *sao);
void sao_search_chroma(const encoder_control * encoder, const picture *pic, unsigned x_ctb, unsigned y_ctb, sao_info *sao, sao_info *sao_top, sao_info *sao_left);
void sao_search_luma(const encoder_control * encoder, const picture *pic, unsigned x_ctb, unsigned y_ctb, sao_info *sao, sao_info *sao_top, sao_info *sao_left);
void sao_reconstruct(const encoder_control * encoder, picture *pic, const pixel *old_rec,
                     unsigned x_ctb, unsigned y_ctb,
                     const sao_info *sao, color_index color_i);
void sao_reconstruct_frame(const encoder_control * const encoder);

#endif
