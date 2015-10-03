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

#include "encoder.h"
#include "encoderstate.h"

typedef struct {
  kvz_pixel left[2 * 32 + 1];
  kvz_pixel top[2 * 32 + 1];
} kvz_intra_ref;
typedef struct
{
  kvz_intra_ref ref;
  kvz_intra_ref filtered_ref;
  bool filtered_initialized;
} kvz_intra_references;

//void kvz_intra_set_block_mode(image* im,uint32_t x_ctb, uint32_t y_ctb, uint8_t depth, uint8_t mode, uint8_t part_mode);

int8_t kvz_intra_get_dir_luma_predictor(uint32_t x, uint32_t y, int8_t* preds,
                                    const cu_info_t* cur_cu, const cu_info_t* left_cu, const cu_info_t* above_cu);
void kvz_intra_dc_pred_filtering(const kvz_pixel* src, int32_t src_stride, kvz_pixel* dst, int32_t dst_stride, int32_t width, int32_t height );

void kvz_intra_build_reference_border(const encoder_control_t *encoder, int32_t x_luma, int32_t y_luma, int16_t out_width, kvz_pixel *dst, int32_t dst_stride, int8_t chroma, int32_t pic_width, int32_t pic_height, lcu_t *lcu);
void kvz_intra_filter(kvz_pixel* ref, int32_t stride, int32_t width, int8_t mode);

/* Predictions */
void kvz_intra_get_pred(const encoder_control_t * const encoder, const kvz_pixel *rec, const kvz_pixel *rec_filtered, int recstride, kvz_pixel *dst, int width, int mode, int is_chroma);

kvz_pixel kvz_intra_get_dc_pred(const kvz_pixel* pic, uint16_t pic_width, uint8_t width);
void kvz_intra_get_planar_pred(const kvz_pixel* src,int32_t srcstride, uint32_t width, kvz_pixel* dst, int32_t dststride);
void kvz_intra_get_angular_pred(const encoder_control_t *encoder, const kvz_pixel* src, int32_t src_stride, kvz_pixel* dst, int32_t dst_stride, int32_t width, int32_t dir_mode, int8_t filter);

/**
* \brief Generage angular predictions.
* \param width    Width in pixels, range 4..32.
* \param color    What color pixels to use.
* \param luma_px  Luma coordinates of the prediction block.
* \param pic_px   Picture dimensions in luma pixels.
* \param lcu      LCU struct.
* \param out_left_ref  Left reference pixels, index 0 is the top-left.
* \param out_top_ref   Top reference pixels, index 0 is the top-left.
*/
void kvz_intra_build_reference(
  const int_fast8_t log2_width,
  const color_t color,
  const vector2d_t *const luma_px,
  const vector2d_t *const pic_px,
  const lcu_t *const lcu,
  kvz_intra_references *const refs);

/**
 * \brief Generate intra predictions.
 * \param refs   Reference pixels used for the prediction.     
 * \param log2_width  Width of the predicted block.
 * \param mode   Intra mode used for the prediction.
 * \param color  Color of the prediction.
 * \param dst    Buffer for the predicted pixels.
 */
void kvz_intra_get_pred_new(
  kvz_intra_references *refs,
  int_fast8_t log2_width,
  int_fast8_t mode,
  color_t color,
  kvz_pixel *dst);


void kvz_intra_recon(const encoder_control_t *encoder, kvz_pixel* rec, int32_t rec_stride, uint32_t width, kvz_pixel* dst, int32_t dst_stride, int8_t mode, int8_t chroma);

void kvz_intra_recon_lcu_luma(encoder_state_t *state, int x, int y, int depth, int8_t intra_mode, cu_info_t *cur_cu, lcu_t *lcu);
void kvz_intra_recon_lcu_chroma(encoder_state_t *state, int x, int y, int depth, int8_t intra_mode, cu_info_t *cur_cu, lcu_t *lcu);

#endif
