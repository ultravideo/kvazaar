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

/**
* \ingroup Reconstruction
* \file
* Intra prediction.
*/

#include "global.h" // IWYU pragma: keep

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


/**
* \brief Function for deriving intra luma predictions
* \param pic picture to use
* \param x_cu x CU position (smallest CU)
* \param y_cu y CU position (smallest CU)
* \param preds output buffer for 3 predictions
* \returns (predictions are found)?1:0
*/
int8_t kvz_intra_get_dir_luma_predictor(
  const uint32_t x,
  const uint32_t y,
  int8_t *preds,
  const cu_info_t *const cur_cu,
  const cu_info_t *const left_cu,
  const cu_info_t *const above_cu);

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
void kvz_intra_predict(
  kvz_intra_references *refs,
  int_fast8_t log2_width,
  int_fast8_t mode,
  color_t color,
  kvz_pixel *dst);

/**
 * \brief Do a full intra prediction cycle on a CU in lcu for luma.
 */
void kvz_intra_recon_lcu_luma(
  encoder_state_t *const state,
  int x,
  int y,
  int depth,
  int8_t intra_mode,
  cu_info_t *cur_cu,
  lcu_t *lcu);

/**
* \brief Do a full intra prediction cycle on a CU in lcu for chroma.
*/
void kvz_intra_recon_lcu_chroma(
  encoder_state_t *const state,
  int x,
  int y,
  int depth,
  int8_t intra_mode,
  cu_info_t *cur_cu,
  lcu_t *lcu);

#endif
