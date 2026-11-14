#ifndef STRATEGIES_PICTURE_GENERIC_H_
#define STRATEGIES_PICTURE_GENERIC_H_
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
 * \ingroup Optimization
 * \file
 * Generic C implementations of optimized functions.
 */

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"

int kvz_strategy_register_picture_generic(void* opaque, uint8_t bitdepth);

// Function to clip int16_t to pixel. (0-255 or 0-1023)
// Assumes PIXEL_MAX to be 2^n-1
kvz_pixel kvz_fast_clip_16bit_to_pixel(int16_t value);

// Function to clip int32_t to pixel. (0-255 or 0-1023)
// Assumes PIXEL_MAX to be 2^n-1
kvz_pixel kvz_fast_clip_32bit_to_pixel(int32_t value);

unsigned kvz_satd_4x4_subblock_generic(const kvz_pixel * buf1,
                                       const int32_t     stride1,
                                       const kvz_pixel * buf2,
                                       const int32_t     stride2);

void kvz_satd_4x4_subblock_quad_generic(const kvz_pixel *preds[4],
                                        const int stride,
                                        const kvz_pixel *orig,
                                        const int orig_stride,
                                        unsigned costs[4]);

// Temporal filter functions, also used for the cases the SIMD versions do
// not handle.
void kvz_bilinear_ssd_multi_generic(const kvz_pixel *org, int org_stride,
                                    const kvz_pixel *ref, int ref_stride,
                                    int size, const int *mv_x, const int *mv_y,
                                    int num, unsigned *out);
void kvz_tf_bilinear_pred_generic(kvz_pixel *dst, int dst_stride,
                                  const kvz_pixel *ref, int ref_stride,
                                  int width, int height,
                                  int w00, int w10, int w01, int w11, int shift);
void kvz_tf_block_stats_generic(const kvz_pixel *org, int org_stride,
                                const kvz_pixel *pred, int pred_stride,
                                int width, int height, int64_t out[5]);
void kvz_tf_filter_block_generic(kvz_pixel *dst, int dst_stride,
                                 const kvz_pixel *src, int src_stride,
                                 const kvz_pixel *const *preds, int pred_stride,
                                 int num_refs, const float *weights,
                                 const float *const *luts,
                                 int width, int height, int max_val);


#endif //STRATEGIES_PICTURE_GENERIC_H_
