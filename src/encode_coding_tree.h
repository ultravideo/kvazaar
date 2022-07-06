#ifndef ENCODE_CODING_TREE_H_
#define ENCODE_CODING_TREE_H_

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
 * \file
 * Functions for writing the coding quadtree and related syntax.
 */

#include "encoderstate.h"
#include "global.h"

void kvz_encode_coding_tree(encoder_state_t * const state,
                            uint16_t x_ctb,
                            uint16_t y_ctb,
                            uint8_t depth);

void kvz_encode_mvd(encoder_state_t * const state,
                    cabac_data_t *cabac,
                    int32_t mvd_hor,
                    int32_t mvd_ver,
                    double* bits_out);

double kvz_mock_encode_coding_unit(
  encoder_state_t* const state,
  cabac_data_t* cabac,
  int x, int y, int depth,
  lcu_t* lcu, cu_info_t* cur_cu);

double kvz_encode_part_mode(encoder_state_t* const state,
  cabac_data_t* const cabac,
  const cu_info_t* const cur_cu,
  int depth);

void kvz_encode_inter_prediction_unit(encoder_state_t* const state,
                                      cabac_data_t* const cabac,
                                      const cu_info_t* const cur_cu,
                                      int x, int y, int width, int height,
                                      int depth, 
                                      lcu_t* lcu,
                                      double* bits_out);

void kvz_encode_last_significant_xy(cabac_data_t * const cabac,
                                    uint8_t lastpos_x, uint8_t lastpos_y,
                                    uint8_t width, uint8_t height,
                                    uint8_t type, uint8_t scan, double* bits_out);

#endif // ENCODE_CODING_TREE_H_
