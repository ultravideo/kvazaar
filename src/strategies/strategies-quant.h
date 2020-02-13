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

#define MIN_FAST_COEFF_COST_QP  0
#define MAX_FAST_COEFF_COST_QP 50
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
// 0 to 50 with ultrafast encoding
static const uint64_t fast_coeff_cost_wts[FAST_COEFF_QP_COUNT] = {
  TO_4XQ88(0.917891, 3.008218, 2.310740, 6.368842),
  TO_4XQ88(0.843318, 3.019534, 2.491460, 6.317865),
  TO_4XQ88(0.768744, 3.030851, 2.672180, 6.266887),
  TO_4XQ88(0.694171, 3.042167, 2.852900, 6.215910),
  TO_4XQ88(0.619597, 3.053484, 3.033620, 6.164932),
  TO_4XQ88(0.545024, 3.064800, 3.214340, 6.113955),
  TO_4XQ88(0.470522, 3.117827, 3.321078, 6.157121),
  TO_4XQ88(0.396021, 3.170854, 3.427816, 6.200287),
  TO_4XQ88(0.321520, 3.223881, 3.534554, 6.243453),
  TO_4XQ88(0.247018, 3.276909, 3.641292, 6.286618),
  TO_4XQ88(0.172517, 3.329936, 3.748030, 6.329784),
  TO_4XQ88(0.151807, 3.398547, 3.752390, 6.320054),
  TO_4XQ88(0.131096, 3.467158, 3.756750, 6.310324),
  TO_4XQ88(0.110386, 3.535770, 3.761109, 6.300594),
  TO_4XQ88(0.089676, 3.604381, 3.765469, 6.290864),
  TO_4XQ88(0.068966, 3.672992, 3.769829, 6.281134),
  TO_4XQ88(0.058852, 3.716881, 3.758086, 6.263672),
  TO_4XQ88(0.048739, 3.760770, 3.746343, 6.246210),
  TO_4XQ88(0.038625, 3.804658, 3.734600, 6.228749),
  TO_4XQ88(0.028512, 3.848547, 3.722858, 6.211287),
  TO_4XQ88(0.018398, 3.892435, 3.711115, 6.193825),
  TO_4XQ88(0.015652, 3.891003, 3.797404, 6.158838),
  TO_4XQ88(0.012906, 3.889571, 3.883692, 6.123852),
  TO_4XQ88(0.010159, 3.888139, 3.969981, 6.088865),
  TO_4XQ88(0.007413, 3.886707, 4.056270, 6.053878),
  TO_4XQ88(0.004667, 3.885275, 4.142559, 6.018891),
  TO_4XQ88(0.004128, 3.870089, 4.226282, 6.035303),
  TO_4XQ88(0.003588, 3.854903, 4.310005, 6.051715),
  TO_4XQ88(0.003049, 3.839716, 4.393728, 6.068127),
  TO_4XQ88(0.002510, 3.824530, 4.477452, 6.084539),
  TO_4XQ88(0.001971, 3.809344, 4.561175, 6.100951),
  TO_4XQ88(0.001763, 3.797488, 4.633904, 6.175179),
  TO_4XQ88(0.001555, 3.785632, 4.706633, 6.249407),
  TO_4XQ88(0.001347, 3.773777, 4.779361, 6.323634),
  TO_4XQ88(0.001139, 3.761921, 4.852090, 6.397862),
  TO_4XQ88(0.000931, 3.750065, 4.924819, 6.472090),
  TO_4XQ88(0.000816, 3.748832, 4.993478, 6.557443),
  TO_4XQ88(0.000701, 3.747598, 5.062136, 6.642795),
  TO_4XQ88(0.000586, 3.746364, 5.130795, 6.728148),
  TO_4XQ88(0.000470, 3.745131, 5.199454, 6.813500),
  TO_4XQ88(0.000355, 3.743898, 5.268112, 6.898853),
  TO_4XQ88(0.000299, 3.767407, 5.286614, 6.968463),
  TO_4XQ88(0.000244, 3.790916, 5.305115, 7.038073),
  TO_4XQ88(0.000188, 3.814425, 5.323616, 7.107684),
  TO_4XQ88(0.000133, 3.837935, 5.342118, 7.177295),
  TO_4XQ88(0.000077, 3.861444, 5.360619, 7.246905),
  TO_4XQ88(0.000062, 3.876310, 5.356118, 7.323495),
  TO_4XQ88(0.000046, 3.891176, 5.351617, 7.400084),
  TO_4XQ88(0.000031, 3.906042, 5.347116, 7.476674),
  TO_4XQ88(0.000015, 3.920908, 5.342615, 7.553263),
  TO_4XQ88(0.000000, 3.935774, 5.338114, 7.629853),
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
