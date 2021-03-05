#ifndef STRATEGIES_IPOL_H_
#define STRATEGIES_IPOL_H_
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
 * Interface for subpixel interpolation functions.
 */

#include "encoder.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "search_inter.h"


typedef struct { kvz_pixel *buffer; kvz_pixel *orig_topleft; unsigned stride; unsigned malloc_used; } kvz_extended_block;

typedef void(ipol_blocks_func)(const encoder_control_t * encoder, kvz_pixel *src, int16_t src_stride, int width, int height,
  kvz_pixel filtered[4][LCU_WIDTH * LCU_WIDTH], int16_t hor_intermediate[5][(KVZ_EXT_BLOCK_W_LUMA + 1) * LCU_WIDTH], int8_t fme_level, int16_t hor_first_cols[5][KVZ_EXT_BLOCK_W_LUMA + 1], 
  int8_t sample_off_x, int8_t sample_off_y);

typedef struct {
  // Source samples
  kvz_pixel *src; // Top-left sample
  int src_w; // Width
  int src_h; // Height
  int src_s; // Stride

  // Requested sampling position, base dimensions, and padding
  int blk_x;
  int blk_y;
  int blk_w; // Width
  int blk_h; // Height
  int pad_l; // Left
  int pad_r; // Right
  int pad_t; // Top
  int pad_b; // Bottom

  // Buffer for possible extrapolation. Free memory provided by the caller.
  kvz_pixel *buf;

  // Extended block data. These are set by the function.
  kvz_pixel **ext; // Top-left sample with padding
  kvz_pixel **ext_origin; // Top-left sample without padding
  int *ext_s; // Stride
} kvz_epol_args;

typedef void(epol_func)(kvz_epol_args *args);


typedef void(kvz_sample_quarterpel_luma_func)(const encoder_control_t * const encoder, kvz_pixel *src, int16_t src_stride, int width, int height, kvz_pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag, const int16_t mv[2]);
typedef void(kvz_sample_octpel_chroma_func)(const encoder_control_t * const encoder, kvz_pixel *src, int16_t src_stride, int width, int height, kvz_pixel *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag, const int16_t mv[2]);

typedef void(kvz_sample_14bit_quarterpel_luma_func)(const encoder_control_t * const encoder, kvz_pixel *src, int16_t src_stride, int width, int height, int16_t *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag, const int16_t mv[2]);
typedef void(kvz_sample_14bit_octpel_chroma_func)(const encoder_control_t * const encoder, kvz_pixel *src, int16_t src_stride, int width, int height, int16_t *dst, int16_t dst_stride, int8_t hor_flag, int8_t ver_flag, const int16_t mv[2]);

// Declare function pointers.
extern ipol_blocks_func * kvz_filter_hpel_blocks_hor_ver_luma;
extern ipol_blocks_func * kvz_filter_hpel_blocks_diag_luma;
extern ipol_blocks_func * kvz_filter_qpel_blocks_hor_ver_luma;
extern ipol_blocks_func * kvz_filter_qpel_blocks_diag_luma;
extern epol_func * kvz_get_extended_block;
extern kvz_sample_quarterpel_luma_func * kvz_sample_quarterpel_luma;
extern kvz_sample_octpel_chroma_func * kvz_sample_octpel_chroma;
extern kvz_sample_14bit_quarterpel_luma_func * kvz_sample_14bit_quarterpel_luma;
extern kvz_sample_14bit_octpel_chroma_func * kvz_sample_14bit_octpel_chroma;


int kvz_strategy_register_ipol(void* opaque, uint8_t bitdepth);


#define STRATEGIES_IPOL_EXPORTS \
  {"filter_hpel_blocks_hor_ver_luma", (void**) &kvz_filter_hpel_blocks_hor_ver_luma}, \
  {"filter_hpel_blocks_diag_luma",    (void**) &kvz_filter_hpel_blocks_diag_luma}, \
  {"filter_qpel_blocks_hor_ver_luma", (void**) &kvz_filter_qpel_blocks_hor_ver_luma}, \
  {"filter_qpel_blocks_diag_luma",    (void**) &kvz_filter_qpel_blocks_diag_luma}, \
  {"sample_quarterpel_luma", (void**) &kvz_sample_quarterpel_luma}, \
  {"sample_octpel_chroma", (void**) &kvz_sample_octpel_chroma}, \
  {"sample_14bit_quarterpel_luma", (void**) &kvz_sample_14bit_quarterpel_luma}, \
  {"sample_14bit_octpel_chroma", (void**) &kvz_sample_14bit_octpel_chroma}, \
  {"get_extended_block", (void**) &kvz_get_extended_block}, \



#endif //STRATEGIES_IPOL_H_
