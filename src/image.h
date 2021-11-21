#ifndef IMAGE_H_
#define IMAGE_H_
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
 * \ingroup DataStructures
 * \file
 * A reference counted YUV pixel buffer.
 */

#include "global.h" // IWYU pragma: keep

#include "kvazaar.h"
#include "strategies/optimized_sad_func_ptr_t.h"


typedef struct {
  kvz_pixel y[LCU_LUMA_SIZE];
  kvz_pixel u[LCU_CHROMA_SIZE];
  kvz_pixel v[LCU_CHROMA_SIZE];
  enum kvz_chroma_format chroma_format;
} lcu_yuv_t;

typedef struct {
  int size;
  kvz_pixel *y;
  kvz_pixel *u;
  kvz_pixel *v;
} yuv_t;

typedef struct {
  int size;
  kvz_pixel_im *y;
  kvz_pixel_im *u;
  kvz_pixel_im *v;
} yuv_im_t;

kvz_picture *kvz_image_alloc_420(const int32_t width, const int32_t height);
kvz_picture *kvz_image_alloc(enum kvz_chroma_format chroma_format, const int32_t width, const int32_t height);

void kvz_image_free(kvz_picture *im);

kvz_picture *kvz_image_copy_ref(kvz_picture *im);

kvz_picture *kvz_image_make_subimage(kvz_picture *const orig_image,
                             const unsigned x_offset,
                             const unsigned y_offset,
                             const unsigned width,
                             const unsigned height);

yuv_t * kvz_yuv_t_alloc(int luma_size, int chroma_size);
void kvz_yuv_t_free(yuv_t * yuv);


//Algorithms
unsigned kvz_image_calc_sad(const kvz_picture *pic,
                            const kvz_picture *ref,
                            int pic_x,
                            int pic_y,
                            int ref_x,
                            int ref_y,
                            int block_width,
                            int block_height,
                            optimized_sad_func_ptr_t optimized_sad);


unsigned kvz_image_calc_satd(const kvz_picture *pic,
                             const kvz_picture *ref,
                             int pic_x,
                             int pic_y,
                             int ref_x,
                             int ref_y,
                             int block_width,
                             int block_height);


void kvz_pixels_blit(const kvz_pixel* orig, kvz_pixel *dst,
                         unsigned width, unsigned height,
                         unsigned orig_stride, unsigned dst_stride);


#endif
