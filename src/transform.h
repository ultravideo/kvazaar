#ifndef TRANSFORM_H_
#define TRANSFORM_H_
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
 * Quantization and transform functions.
 */

#include "cu.h"
#include "encoder.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep

extern const uint8_t kvz_g_chroma_scale[58];
extern const int16_t kvz_g_inv_quant_scales[6];

void kvz_transformskip(const encoder_control_t *encoder, int16_t *block,int16_t *coeff, int8_t block_size);
void kvz_itransformskip(const encoder_control_t *encoder, int16_t *block,int16_t *coeff, int8_t block_size);

void kvz_transform2d(const encoder_control_t * const encoder,
                     int16_t *block,
                     int16_t *coeff,
                     int8_t block_size,
                     color_t color,
                     cu_type_t type);
void kvz_itransform2d(const encoder_control_t * const encoder,
                      int16_t *block,
                      int16_t *coeff,
                      int8_t block_size,
                      color_t color,
                      cu_type_t type);

int32_t kvz_get_scaled_qp(int8_t type, int8_t qp, int8_t qp_offset);

void kvz_quantize_lcu_residual(encoder_state_t *state,
                               bool luma,
                               bool chroma,
                               int32_t x,
                               int32_t y,
                               uint8_t depth,
                               cu_info_t *cur_cu,
                               lcu_t* lcu,
                               bool early_skip);

#endif
