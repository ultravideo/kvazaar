#ifndef INTRA_H_
#define INTRA_H_
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
* Intra prediction.
*/

#include "cu.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"


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
* \param x          x-coordinate of the PU in pixels
* \param y          y-coordinate of the PU in pixels
* \param preds      output buffer for 3 predictions
* \param cur_pu     PU to check
* \param left_pu    PU to the left of cur_pu
* \param above_pu   PU above cur_pu
* \returns          1 if predictions are found, otherwise 0
*/
int8_t kvz_intra_get_dir_luma_predictor(
  const uint32_t x,
  const uint32_t y,
  int8_t *preds,
  const cu_info_t *const cur_pu,
  const cu_info_t *const left_pu,
  const cu_info_t *const above_pu);

#if KVZ_SEL_ENCRYPTION
/**
* \brief Function for deriving intra luma predictions with encryption
* \param x          x-coordinate of the PU in pixels
* \param y          y-coordinate of the PU in pixels
* \param preds      output buffer for 3 predictions
* \param cur_pu     PU to check
* \param left_pu    PU to the left of cur_pu
* \param above_pu   PU above cur_pu
* \returns          1 if predictions are found, otherwise 0
*/
int8_t kvz_intra_get_dir_luma_predictor_encry(
const uint32_t x,
const uint32_t y,
int8_t *preds,
const cu_info_t *const cur_pu,
const cu_info_t *const left_pu,
const cu_info_t *const above_pu);
#endif

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
 * \param refs            Reference pixels used for the prediction.
 * \param log2_width      Width of the predicted block.
 * \param mode            Intra mode used for the prediction.
 * \param color           Color of the prediction.
 * \param dst             Buffer for the predicted pixels.
 * \param filter_boundary Whether to filter the boundary on modes 10 and 26.
 */
void kvz_intra_predict(
  kvz_intra_references *refs,
  int_fast8_t log2_width,
  int_fast8_t mode,
  color_t color,
  kvz_pixel *dst,
  bool filter_boundary);

void kvz_intra_recon_cu(
  encoder_state_t *const state,
  int x,
  int y,
  int depth,
  int8_t mode_luma,
  int8_t mode_chroma,
  cu_info_t *cur_cu,
  lcu_t *lcu);

#endif
