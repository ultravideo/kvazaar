#ifndef INTRA_H_
#define INTRA_H_
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
 * \brief Handling Coding Units (CU's) for intra frames.
 */

#include "global.h"

#include "image.h"
#include "encoder.h"
#include "encoderstate.h"

//void intra_set_block_mode(image* im,uint32_t x_ctb, uint32_t y_ctb, uint8_t depth, uint8_t mode, uint8_t part_mode);

int8_t intra_get_dir_luma_predictor(uint32_t x, uint32_t y, int8_t* preds,
                                    const cu_info* cur_cu, const cu_info* left_cu, const cu_info* above_cu);
void intra_dc_pred_filtering(const pixel* src, int32_t src_stride, pixel* dst, int32_t dst_stride, int32_t width, int32_t height );

void intra_build_reference_border(const encoder_control *encoder, int32_t x_luma, int32_t y_luma, int16_t out_width, pixel *dst, int32_t dst_stride, int8_t chroma, int32_t pic_width, int32_t pic_height, lcu_t *lcu);
void intra_filter(pixel* ref, int32_t stride, int32_t width, int8_t mode);

/* Predictions */
void intra_get_pred(const encoder_control * const encoder, const pixel *rec, const pixel *rec_filtered, int recstride, pixel *dst, int width, int mode, int is_chroma);

pixel intra_get_dc_pred(const pixel* pic, uint16_t pic_width, uint8_t width);
void intra_get_planar_pred(const pixel* src,int32_t srcstride, uint32_t width, pixel* dst, int32_t dststride);
void intra_get_angular_pred(const encoder_control *encoder, const pixel* src, int32_t src_stride, pixel* dst, int32_t dst_stride, int32_t width, int32_t dir_mode, int8_t filter);

void intra_recon(const encoder_control *encoder, pixel* rec, int32_t rec_stride, uint32_t width, pixel* dst, int32_t dst_stride, int8_t mode, int8_t chroma);

void intra_recon_lcu_luma(encoder_state *encoder_state, int x, int y, int depth, int8_t intra_mode, cu_info *cur_cu, lcu_t *lcu);
void intra_recon_lcu_chroma(encoder_state *encoder_state, int x, int y, int depth, int8_t intra_mode, cu_info *cur_cu, lcu_t *lcu);

#endif
