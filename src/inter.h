#ifndef INTER_H_
#define INTER_H_
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
 * Inter prediction.
 */

#include "cu.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "image.h"
#include "kvazaar.h"


typedef struct {
  uint8_t dir;
  uint8_t ref[2]; // index to L0/L1
  int16_t mv[2][2];

} inter_merge_cand_t;

void kvz_inter_recon_cu(const encoder_state_t * const state,
                        lcu_t *lcu,
                        int32_t x,
                        int32_t y,
                        int32_t width,
                        bool predict_luma,
                        bool predict_chroma);

void kvz_inter_pred_pu(const encoder_state_t * const state,
  lcu_t *lcu,
  int32_t x,
  int32_t y,
  int32_t width,
  bool predict_luma,
  bool predict_chroma,
  int i_pu);

void kvz_inter_recon_bipred(const encoder_state_t * const state,
                            const kvz_picture * ref1,
                            const kvz_picture * ref2,
                            int32_t xpos,
                            int32_t ypos,
                            int32_t width,
                            int32_t height,
                            int16_t mv_param[2][2],
                            lcu_t* lcu,
                            bool predict_luma,
                            bool predict_chroma);


void kvz_inter_get_mv_cand(const encoder_state_t * const state,
                           int32_t x,
                           int32_t y,
                           int32_t width,
                           int32_t height,
                           int16_t mv_cand[2][2],
                           const cu_info_t* cur_cu,
                           lcu_t *lcu,
                           int8_t reflist);

void kvz_inter_get_mv_cand_cua(const encoder_state_t * const state,
                               int32_t x,
                               int32_t y,
                               int32_t width,
                               int32_t height,
                               int16_t mv_cand[2][2],
                               const cu_info_t* cur_cu,
                               int8_t reflist);

uint8_t kvz_inter_get_merge_cand(const encoder_state_t * const state,
                                 int32_t x, int32_t y,
                                 int32_t width, int32_t height,
                                 bool use_a1, bool use_b1,
                                 inter_merge_cand_t mv_cand[MRG_MAX_NUM_CANDS],
                                 lcu_t *lcu);
#endif
