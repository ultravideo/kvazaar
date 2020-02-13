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
// 0 to 50 with ultrafast encoding
static const uint64_t fast_coeff_cost_wts[FAST_COEFF_QP_COUNT] = {
  TO_4XQ88(0.139807, 4.172847, 3.442430, 6.573553),
  TO_4XQ88(0.124146, 4.267273, 3.413249, 6.526515),
  TO_4XQ88(0.108484, 4.361699, 3.384069, 6.479476),
  TO_4XQ88(0.095441, 4.487631, 3.222517, 6.486733),
  TO_4XQ88(0.082397, 4.613562, 3.060965, 6.493991),
  TO_4XQ88(0.071664, 4.681115, 2.974996, 6.467879),
  TO_4XQ88(0.060932, 4.748667, 2.889028, 6.441766),
  TO_4XQ88(0.050267, 4.803366, 2.892400, 6.406573),
  TO_4XQ88(0.039603, 4.858065, 2.895773, 6.371380),
  TO_4XQ88(0.032558, 4.886846, 2.909959, 6.328989),
  TO_4XQ88(0.025512, 4.915628, 2.924145, 6.286599),
  TO_4XQ88(0.022040, 4.927304, 2.909608, 6.260581),
  TO_4XQ88(0.018568, 4.938981, 2.895071, 6.234563),
  TO_4XQ88(0.015158, 4.940147, 3.020277, 6.155190),
  TO_4XQ88(0.011747, 4.941313, 3.145482, 6.075818),
  TO_4XQ88(0.009987, 4.944803, 3.220803, 6.045646),
  TO_4XQ88(0.008226, 4.948294, 3.296123, 6.015475),
  TO_4XQ88(0.007048, 4.953393, 3.364038, 5.998843),
  TO_4XQ88(0.005870, 4.958493, 3.431954, 5.982210),
  TO_4XQ88(0.005085, 4.973462, 3.462797, 5.969676),
  TO_4XQ88(0.004299, 4.988431, 3.493640, 5.957141),
  TO_4XQ88(0.003679, 5.006778, 3.542994, 5.947307),
  TO_4XQ88(0.003060, 5.025124, 3.592349, 5.937474),
  TO_4XQ88(0.002626, 5.041786, 3.651788, 5.939799),
  TO_4XQ88(0.002191, 5.058447, 3.711226, 5.942125),
  TO_4XQ88(0.001906, 5.095488, 3.735999, 5.940480),
  TO_4XQ88(0.001621, 5.132529, 3.760771, 5.938835),
  TO_4XQ88(0.001383, 5.171319, 3.794231, 5.947102),
  TO_4XQ88(0.001144, 5.210109, 3.827690, 5.955370),
  TO_4XQ88(0.000965, 5.237812, 3.881474, 6.032188),
  TO_4XQ88(0.000787, 5.265515, 3.935259, 6.109006),
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
