#ifndef STRATEGIES_QUANT_H_
#define STRATEGIES_QUANT_H_
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
 * \ingroup Optimization
 * \file
 * Interface for quantization functions.
 */

#include "cu.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "tables.h"

#define MIN_FAST_COEFF_COST_QP 12
#define MAX_FAST_COEFF_COST_QP 42
#define FAST_COEFF_QP_COUNT ((MAX_FAST_COEFF_COST_QP) - (MIN_FAST_COEFF_COST_QP) + (1))

// Note: Assumes that costs are non-negative, for pretty obvious reasons
#define TO_Q88(f) ((uint16_t)((f) * (256.0f) + (0.5f)))

#define TO_4XQ88(f0,f1,f2,f3) ( \
                                ((uint64_t)TO_Q88((f0)) <<  0) | \
                                ((uint64_t)TO_Q88((f1)) << 16) | \
                                ((uint64_t)TO_Q88((f2)) << 32) | \
                                ((uint64_t)TO_Q88((f3)) << 48)   \
                              )

// Weights for 4 buckets (coeff 0, coeff 1, coeff 2, coeff >= 3), for QPs from
// MIN_FAST_COEFF_COST_QP to MAX_FAST_COEFF_COST_QP with ultrafast encoding
static const uint64_t fast_coeff_cost_wts[FAST_COEFF_QP_COUNT] = {
  TO_4XQ88(0.134012, 4.207784, 3.431633, 6.556149),
  TO_4XQ88(0.122972, 4.278606, 3.398710, 6.527168),
  TO_4XQ88(0.109453, 4.373356, 3.335091, 6.499565),
  TO_4XQ88(0.095884, 4.485212, 3.217406, 6.488617),
  TO_4XQ88(0.083252, 4.591962, 3.088931, 6.481644),
  TO_4XQ88(0.071878, 4.674704, 2.989839, 6.464058),
  TO_4XQ88(0.060957, 4.743912, 2.922084, 6.438406),
  TO_4XQ88(0.050599, 4.799877, 2.901414, 6.405108),
  TO_4XQ88(0.040942, 4.848476, 2.899774, 6.368717),
  TO_4XQ88(0.033205, 4.882974, 2.908347, 6.329815),
  TO_4XQ88(0.026834, 4.909299, 2.913517, 6.292657),
  TO_4XQ88(0.022367, 4.924819, 2.919600, 6.257252),
  TO_4XQ88(0.018591, 4.935092, 2.946776, 6.214822),
  TO_4XQ88(0.015312, 4.939410, 3.028364, 6.154816),
  TO_4XQ88(0.012358, 4.942173, 3.127025, 6.094022),
  TO_4XQ88(0.010188, 4.945157, 3.215646, 6.051293),
  TO_4XQ88(0.008442, 4.948889, 3.293383, 6.020484),
  TO_4XQ88(0.007136, 4.954426, 3.360035, 6.000430),
  TO_4XQ88(0.006015, 4.962144, 3.418237, 5.983726),
  TO_4XQ88(0.005135, 4.974654, 3.461126, 5.970288),
  TO_4XQ88(0.004360, 4.989681, 3.500489, 5.958140),
  TO_4XQ88(0.003711, 5.006930, 3.545568, 5.948645),
  TO_4XQ88(0.003128, 5.024501, 3.596080, 5.941973),
  TO_4XQ88(0.002656, 5.043468, 3.649575, 5.940536),
  TO_4XQ88(0.002246, 5.065988, 3.698400, 5.940656),
  TO_4XQ88(0.001924, 5.097480, 3.733661, 5.941015),
  TO_4XQ88(0.001638, 5.133176, 3.763985, 5.942503),
  TO_4XQ88(0.001392, 5.170478, 3.796842, 5.954164),
  TO_4XQ88(0.001166, 5.206007, 3.835210, 5.980734),
  TO_4XQ88(0.000987, 5.234321, 3.878463, 6.031444),
  TO_4XQ88(0.000853, 5.255265, 3.915359, 6.080584),
};

#undef TO_4XQ88
#undef TO_Q88

// Declare function pointers.
typedef unsigned (quant_func)(const encoder_state_t * const state, coeff_t *coef, coeff_t *q_coef, int32_t width,
  int32_t height, int8_t type, int8_t scan_idx, int8_t block_type);
typedef unsigned (quant_residual_func)(encoder_state_t *const state,
  const cu_info_t *const cur_cu, const int width, const color_t color,
  const coeff_scan_order_t scan_order, const int use_trskip,
  const int in_stride, const int out_stride,
  const kvz_pixel *const ref_in, const kvz_pixel *const pred_in,
  kvz_pixel *rec_out, coeff_t *coeff_out,
  bool early_skip);
typedef unsigned (dequant_func)(const encoder_state_t * const state, coeff_t *q_coef, coeff_t *coef, int32_t width,
  int32_t height, int8_t type, int8_t block_type);
typedef uint32_t (fast_coeff_cost_func)(const coeff_t *coeff, int32_t width, int32_t qp);

typedef uint32_t (coeff_abs_sum_func)(const coeff_t *coeffs, size_t length);

// Declare function pointers.
extern quant_func * kvz_quant;
extern quant_residual_func * kvz_quantize_residual;
extern dequant_func *kvz_dequant;
extern coeff_abs_sum_func *kvz_coeff_abs_sum;
extern fast_coeff_cost_func *kvz_fast_coeff_cost;

int kvz_strategy_register_quant(void* opaque, uint8_t bitdepth);


#define STRATEGIES_QUANT_EXPORTS \
  {"quant", (void**) &kvz_quant}, \
  {"quantize_residual", (void**) &kvz_quantize_residual}, \
  {"dequant", (void**) &kvz_dequant}, \
  {"coeff_abs_sum", (void**) &kvz_coeff_abs_sum}, \
  {"fast_coeff_cost", (void**) &kvz_fast_coeff_cost}, \



#endif //STRATEGIES_QUANT_H_
