#ifndef TRANSFORM_H_
#define TRANSFORM_H_
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
 * \brief Transformations, such as quantization and DST.
 */

#include "global.h"

#include "encoder.h"
#include "encoderstate.h"

extern const uint8_t g_chroma_scale[58];
extern const int16_t g_inv_quant_scales[6];



void quant(const encoder_state_t *encoder_state, int16_t *coef, int16_t *q_coef, int32_t width,
           int32_t height, int8_t type, int8_t scan_idx, int8_t block_type);
void dequant(const encoder_state_t *encoder_state, int16_t *q_coef, int16_t *coef, int32_t width, int32_t height,int8_t type, int8_t block_type);

void transformskip(const encoder_control_t *encoder, int16_t *block,int16_t *coeff, int8_t block_size);
void itransformskip(const encoder_control_t *encoder, int16_t *block,int16_t *coeff, int8_t block_size);

void transform2d(const encoder_control_t *encoder, int16_t *block,int16_t *coeff, int8_t block_size, int32_t mode);
void itransform2d(const encoder_control_t *encoder, int16_t *block,int16_t *coeff, int8_t block_size, int32_t mode);

int32_t get_scaled_qp(int8_t type, int8_t qp, int8_t qp_offset);

void quantize_lcu_luma_residual(encoder_state_t *encoder_state, int32_t x, int32_t y, uint8_t depth, cu_info_t *cur_cu, lcu_t* lcu);
void quantize_lcu_chroma_residual(encoder_state_t *encoder_state, int32_t x, int32_t y, uint8_t depth, cu_info_t *cur_cu, lcu_t* lcu);

#endif
