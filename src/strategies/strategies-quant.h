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
#define TO_Q88(f) ((uint16_t)((f) * (256.0f)))

#define TO_4XQ88(f0,f1,f2,f3) ( \
                                ((uint64_t)TO_Q88((f0)) <<  0) | \
                                ((uint64_t)TO_Q88((f1)) << 16) | \
                                ((uint64_t)TO_Q88((f2)) << 32) | \
                                ((uint64_t)TO_Q88((f3)) << 48)   \
                              )

// Weights for 4 buckets (coeff 0, coeff 1, coeff 2, coeff >= 3), for QPs from
// 0 to 50 with ultrafast encoding
static const uint64_t fast_coeff_cost_wts[FAST_COEFF_QP_COUNT] = {
  TO_4XQ88(0.783627, 2.844662, 2.204238, 6.229038),
  TO_4XQ88(0.723039, 2.875346, 2.399723, 6.191952),
  TO_4XQ88(0.662451, 2.906030, 2.595208, 6.154866),
  TO_4XQ88(0.601863, 2.936714, 2.790693, 6.117780),
  TO_4XQ88(0.541275, 2.967398, 2.986178, 6.080694),
  TO_4XQ88(0.480687, 2.998082, 3.181662, 6.043608),
  TO_4XQ88(0.411525, 3.057255, 3.277193, 6.089339),
  TO_4XQ88(0.342363, 3.116427, 3.372724, 6.135070),
  TO_4XQ88(0.273201, 3.175600, 3.468255, 6.180802),
  TO_4XQ88(0.204039, 3.234772, 3.563786, 6.226533),
  TO_4XQ88(0.134877, 3.293944, 3.659317, 6.272264),
  TO_4XQ88(0.120459, 3.366286, 3.676365, 6.274461),
  TO_4XQ88(0.106042, 3.438629, 3.693413, 6.276658),
  TO_4XQ88(0.091625, 3.510971, 3.710461, 6.278855),
  TO_4XQ88(0.077208, 3.583313, 3.727508, 6.281052),
  TO_4XQ88(0.062791, 3.655655, 3.744556, 6.283249),
  TO_4XQ88(0.053615, 3.702390, 3.737735, 6.265943),
  TO_4XQ88(0.044440, 3.749125, 3.730914, 6.248637),
  TO_4XQ88(0.035264, 3.795861, 3.724093, 6.231331),
  TO_4XQ88(0.026089, 3.842596, 3.717272, 6.214025),
  TO_4XQ88(0.016913, 3.889331, 3.710451, 6.196719),
  TO_4XQ88(0.014510, 3.881348, 3.789645, 6.153557),
  TO_4XQ88(0.012106, 3.873365, 3.868839, 6.110395),
  TO_4XQ88(0.009703, 3.865382, 3.948033, 6.067233),
  TO_4XQ88(0.007299, 3.857399, 4.027226, 6.024072),
  TO_4XQ88(0.004896, 3.849416, 4.106420, 5.980910),
  TO_4XQ88(0.004573, 3.838888, 4.203887, 6.019835),
  TO_4XQ88(0.004251, 3.828359, 4.301355, 6.058760),
  TO_4XQ88(0.003928, 3.817830, 4.398822, 6.097685),
  TO_4XQ88(0.003606, 3.807302, 4.496290, 6.136611),
  TO_4XQ88(0.003283, 3.796773, 4.593757, 6.175536),
  TO_4XQ88(0.003089, 3.789503, 4.657216, 6.230311),
  TO_4XQ88(0.002894, 3.782233, 4.720675, 6.285086),
  TO_4XQ88(0.002700, 3.774963, 4.784134, 6.339861),
  TO_4XQ88(0.002505, 3.767693, 4.847592, 6.394636),
  TO_4XQ88(0.002311, 3.760423, 4.911051, 6.449411),
  TO_4XQ88(0.002032, 3.762067, 4.987400, 6.588120),
  TO_4XQ88(0.001753, 3.763711, 5.063748, 6.726829),
  TO_4XQ88(0.001473, 3.765354, 5.140096, 6.865538),
  TO_4XQ88(0.001194, 3.766999, 5.216444, 7.004246),
  TO_4XQ88(0.000915, 3.768642, 5.292793, 7.142955),
  TO_4XQ88(0.000798, 3.766159, 5.325907, 7.237017),
  TO_4XQ88(0.000681, 3.763676, 5.359022, 7.331078),
  TO_4XQ88(0.000563, 3.761193, 5.392137, 7.425140),
  TO_4XQ88(0.000446, 3.758710, 5.425252, 7.519202),
  TO_4XQ88(0.000329, 3.756227, 5.458367, 7.613263),
  TO_4XQ88(0.000273, 3.764147, 5.474367, 7.692495),
  TO_4XQ88(0.000217, 3.772068, 5.490368, 7.771726),
  TO_4XQ88(0.000160, 3.779989, 5.506369, 7.850957),
  TO_4XQ88(0.000104, 3.787909, 5.522369, 7.930189),
  TO_4XQ88(0.000048, 3.795830, 5.538370, 8.009420),
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
