#ifndef INTRA_H_
#define INTRA_H_
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
 * \brief Handling Coding Units (CU's) for intra frames.
 */

#include "global.h"

#include "picture.h"
#include "encoder.h"

void intra_set_block_mode(picture* pic,uint32_t x_ctb, uint32_t y_ctb, uint8_t depth, uint8_t mode, uint8_t part_mode);

int8_t intra_get_dir_luma_predictor(uint32_t x, uint32_t y, int8_t* preds,
                                    cu_info* cur_cu, cu_info* left_cu, cu_info* above_cu);
void intra_dc_pred_filtering(pixel* src, int32_t src_stride, pixel* dst, int32_t dst_stride, int32_t width, int32_t height );

void intra_build_reference_border(int32_t x_luma, int32_t y_luma, int16_t out_width, pixel *dst, int32_t dst_stride, int8_t chroma, int32_t pic_width, int32_t pic_height, lcu_t *lcu);
void intra_filter(pixel* ref, int32_t stride, int32_t width, int8_t mode);

/* Predictions */
int16_t intra_prediction(pixel *orig, int32_t origstride, pixel *rec, int16_t recstride,
                         uint8_t width, uint32_t *sad_out,
                         int8_t *intra_preds, uint32_t *bitcost_out);

pixel intra_get_dc_pred(pixel* pic, uint16_t pic_width, uint8_t width);
void intra_get_planar_pred(pixel* src,int32_t srcstride, uint32_t width, pixel* dst, int32_t dststride);
void intra_get_angular_pred(pixel* src, int32_t src_stride, pixel* p_dst, int32_t dst_stride, int32_t width, int32_t dir_mode, int8_t filter);

void intra_recon(pixel* rec, uint32_t rec_stride, uint32_t width, pixel* dst, int32_t dst_stride, int8_t mode, int8_t chroma);

void intra_recon_lcu(encoder_control* encoder, cabac_data *cabac, int x, int y, int depth, lcu_t *lcu, uint32_t pic_width, uint32_t pic_height);

#endif
