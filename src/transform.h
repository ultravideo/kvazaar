#ifndef TRANSFORM_H_
#define TRANSFORM_H_
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
 * \brief Transformations, such as quantization and DST.
 */

#include "global.h"

#include "encoder.h"


extern int32_t* g_quant_coeff[4][6][6];
extern const int32_t g_quant_intra_default_8x8[64];
extern const uint8_t g_chroma_scale[58];


void quant(encoder_control *encoder, int16_t *coef, int16_t *q_coef, int32_t width,
           int32_t height, uint32_t *ac_sum, int8_t type, int8_t scan_idx, int8_t block_type);
void dequant(encoder_control *encoder, int16_t *q_coef, int16_t *coef, int32_t width, int32_t height,int8_t type, int8_t block_type);

void transform2d(int16_t *block,int16_t *coeff, int8_t block_size, int32_t mode);
void itransform2d(int16_t *block,int16_t *coeff, int8_t block_size, int32_t mode);

void scalinglist_init();
void scalinglist_process_enc( int32_t *coeff, int32_t *quant_coeff, int32_t quant_scales, 
                             uint32_t height,uint32_t width, uint32_t ratio, int32_t size_num, uint32_t dc, uint8_t flat);
void scalinglist_process();
void scalinglist_set(int32_t *coeff, uint32_t list_id, uint32_t size_id, uint32_t qp);
void scalinglist_destroy();

#endif
