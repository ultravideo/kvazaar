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

#include <math.h>

#define SCALING_LIST_4x4      0
#define SCALING_LIST_8x8      1
#define SCALING_LIST_16x16    2
#define SCALING_LIST_32x32    3
#define SCALING_LIST_SIZE_NUM 4
#define SCALING_LIST_NUM      6
#define MAX_MATRIX_COEF_NUM   64

extern uint8_t  g_scaling_list_enable;
extern int32_t  g_scaling_list_dc   [SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM];
extern int32_t* g_scaling_list_coeff[SCALING_LIST_SIZE_NUM][SCALING_LIST_NUM];
extern int32_t* g_quant_coeff[4][6][6];
extern double* g_error_scale[4][6][6];
extern const uint8_t g_scaling_list_num[4];
extern const uint16_t g_scaling_list_size[4];
extern const int32_t g_quant_intra_default_8x8[64];
extern const uint8_t g_chroma_scale[58];
extern const int16_t g_inv_quant_scales[6];


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
void scalinglist_set_err_scale(uint32_t list,uint32_t size, uint32_t qp);
void scalinglist_destroy();

int32_t *scalinglist_get_default(uint32_t size_id, uint32_t list_id);
int scalinglist_parse(FILE *fp);

int32_t get_scaled_qp(int8_t type, int8_t qp, int8_t qp_offset);

#endif
