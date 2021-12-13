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

#include "inter.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "encoder.h"
#include "imagelist.h"
#include "strategies/generic/picture-generic.h"
#include "strategies/strategies-ipol.h"
#include "videoframe.h"
#include "strategies/strategies-picture.h"


typedef struct {
  const cu_info_t *a[2];
  const cu_info_t *b[3];
  const cu_info_t *c3;
  const cu_info_t *h;
} merge_candidates_t;


static void inter_recon_frac_luma(const encoder_state_t * const state,
                                  const kvz_picture * const ref,
                                  int32_t xpos,
                                  int32_t ypos,
                                  int32_t block_width,
                                  int32_t block_height,
                                  const int16_t mv_param[2],
                                  yuv_t *out,
                                  unsigned out_stride)
{
  int mv_frac_x = (mv_param[0] & 3);
  int mv_frac_y = (mv_param[1] & 3);

  // Space for extrapolated pixels and the part from the picture.
  // Some extra for AVX2.
  // The extrapolation function will set the pointers and stride.
  kvz_pixel ext_buffer[KVZ_IPOL_MAX_INPUT_SIZE_LUMA_SIMD];
  kvz_pixel *ext = NULL;
  kvz_pixel *ext_origin = NULL;
  int ext_s = 0;
  kvz_epol_args epol_args = {
    .src = ref->y,
    .src_w = ref->width,
    .src_h = ref->height,
    .src_s = ref->stride,
    .blk_x = state->tile->offset_x + xpos + (mv_param[0] >> 2),
    .blk_y = state->tile->offset_y + ypos + (mv_param[1] >> 2),
    .blk_w = block_width,
    .blk_h = block_height,
    .pad_l = KVZ_LUMA_FILTER_OFFSET,
    .pad_r = KVZ_EXT_PADDING_LUMA - KVZ_LUMA_FILTER_OFFSET,
    .pad_t = KVZ_LUMA_FILTER_OFFSET,
    .pad_b = KVZ_EXT_PADDING_LUMA - KVZ_LUMA_FILTER_OFFSET,
    .pad_b_simd = 1 // One row for AVX2
  };

  // Initialize separately. Gets rid of warning
  // about using nonstandard extension.
  epol_args.buf = ext_buffer;
  epol_args.ext = &ext;
  epol_args.ext_origin = &ext_origin;
  epol_args.ext_s = &ext_s;

  kvz_get_extended_block(&epol_args);
  kvz_sample_quarterpel_luma(state->encoder_control,
    ext_origin,
    ext_s,
    block_width,
    block_height,
    out->y,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);
}

static void inter_recon_frac_luma_hi(const encoder_state_t *const state,
  const kvz_picture *const ref,
  int32_t xpos,
  int32_t ypos,
  int32_t block_width,
  int32_t block_height,
  const int16_t mv_param[2],
  yuv_im_t *out,
  const unsigned out_stride)
{
  int mv_frac_x = (mv_param[0] & 3);
  int mv_frac_y = (mv_param[1] & 3);

  // Space for extrapolated pixels and the part from the picture.
  // Some extra for AVX2.
  // The extrapolation function will set the pointers and stride.
  kvz_pixel ext_buffer[KVZ_IPOL_MAX_INPUT_SIZE_LUMA_SIMD];
  kvz_pixel *ext = NULL;
  kvz_pixel *ext_origin = NULL;
  int ext_s = 0;
  kvz_epol_args epol_args = {
    .src = ref->y,
    .src_w = ref->width,
    .src_h = ref->height,
    .src_s = ref->stride,
    .blk_x = state->tile->offset_x + xpos + (mv_param[0] >> 2),
    .blk_y = state->tile->offset_y + ypos + (mv_param[1] >> 2),
    .blk_w = block_width,
    .blk_h = block_height,
    .pad_l = KVZ_LUMA_FILTER_OFFSET,
    .pad_r = KVZ_EXT_PADDING_LUMA - KVZ_LUMA_FILTER_OFFSET,
    .pad_t = KVZ_LUMA_FILTER_OFFSET,
    .pad_b = KVZ_EXT_PADDING_LUMA - KVZ_LUMA_FILTER_OFFSET,
    .pad_b_simd = 1 // One row for AVX2
  };

  // Initialize separately. Gets rid of warning
  // about using nonstandard extension.
  epol_args.buf = ext_buffer;
  epol_args.ext = &ext;
  epol_args.ext_origin = &ext_origin;
  epol_args.ext_s = &ext_s;

  kvz_get_extended_block(&epol_args);
  kvz_sample_quarterpel_luma_hi(state->encoder_control,
    ext_origin,
    ext_s,
    block_width,
    block_height,
    out->y,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);
}

static void inter_recon_frac_chroma(const encoder_state_t *const state,
  const kvz_picture *const ref,
  int32_t pu_x,
  int32_t pu_y,
  int32_t pu_w,
  int32_t pu_h,
  const int16_t mv_param[2],
  yuv_t *out,
  const unsigned out_stride)
{
  int mv_frac_x = (mv_param[0] & 7);
  int mv_frac_y = (mv_param[1] & 7);

  // Take into account chroma subsampling
  unsigned pb_w = pu_w / 2;
  unsigned pb_h = pu_h / 2;

  // Space for extrapolated pixels and the part from the picture.
  // Some extra for AVX2.
  // The extrapolation function will set the pointers and stride.
  kvz_pixel ext_buffer[KVZ_IPOL_MAX_INPUT_SIZE_CHROMA_SIMD];
  kvz_pixel *ext = NULL;
  kvz_pixel *ext_origin = NULL;
  int ext_s = 0;

  // Chroma U
  // Divisions by 2 due to 4:2:0 chroma subsampling
  kvz_epol_args epol_args = {
    .src = ref->u,
    .src_w = ref->width / 2,
    .src_h = ref->height / 2,
    .src_s = ref->stride / 2,
    .blk_x = (state->tile->offset_x + pu_x) / 2 + (mv_param[0] >> 3),
    .blk_y = (state->tile->offset_y + pu_y) / 2 + (mv_param[1] >> 3),
    .blk_w = pb_w,
    .blk_h = pb_h,
    .pad_l = KVZ_CHROMA_FILTER_OFFSET,
    .pad_r = KVZ_EXT_PADDING_CHROMA - KVZ_CHROMA_FILTER_OFFSET,
    .pad_t = KVZ_CHROMA_FILTER_OFFSET,
    .pad_b = KVZ_EXT_PADDING_CHROMA - KVZ_CHROMA_FILTER_OFFSET,
    .pad_b_simd = 3 // Three rows for AVX2
  };

  // Initialize separately. Gets rid of warning
  // about using nonstandard extension.
  epol_args.buf = ext_buffer;
  epol_args.ext = &ext;
  epol_args.ext_origin = &ext_origin;
  epol_args.ext_s = &ext_s;

  kvz_get_extended_block(&epol_args);
  kvz_sample_octpel_chroma(state->encoder_control,
    ext_origin,
    ext_s,
    pb_w,
    pb_h,
    out->u,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);

  // Chroma V
  epol_args.src = ref->v;
  kvz_get_extended_block(&epol_args);
  kvz_sample_octpel_chroma(state->encoder_control,
    ext_origin,
    ext_s,
    pb_w,
    pb_h,
    out->v,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);
}

static void inter_recon_frac_chroma_hi(const encoder_state_t *const state,
  const kvz_picture *const ref,
  int32_t pu_x,
  int32_t pu_y,
  int32_t pu_w,
  int32_t pu_h,
  const int16_t mv_param[2],
  yuv_im_t *out,
  const unsigned out_stride)
{
  int mv_frac_x = (mv_param[0] & 7);
  int mv_frac_y = (mv_param[1] & 7);

  // Take into account chroma subsampling
  unsigned pb_w = pu_w / 2;
  unsigned pb_h = pu_h / 2;

  // Space for extrapolated pixels and the part from the picture.
  // Some extra for AVX2.
  // The extrapolation function will set the pointers and stride.
  kvz_pixel ext_buffer[KVZ_IPOL_MAX_INPUT_SIZE_CHROMA_SIMD];
  kvz_pixel *ext = NULL;
  kvz_pixel *ext_origin = NULL;
  int ext_s = 0;

  // Chroma U
  // Divisions by 2 due to 4:2:0 chroma subsampling
  kvz_epol_args epol_args = {
    .src = ref->u,
    .src_w = ref->width / 2,
    .src_h = ref->height / 2,
    .src_s = ref->stride / 2,
    .blk_x = (state->tile->offset_x + pu_x) / 2 + (mv_param[0] >> 3),
    .blk_y = (state->tile->offset_y + pu_y) / 2 + (mv_param[1] >> 3),
    .blk_w = pb_w,
    .blk_h = pb_h,
    .pad_l = KVZ_CHROMA_FILTER_OFFSET,
    .pad_r = KVZ_EXT_PADDING_CHROMA - KVZ_CHROMA_FILTER_OFFSET,
    .pad_t = KVZ_CHROMA_FILTER_OFFSET,
    .pad_b = KVZ_EXT_PADDING_CHROMA - KVZ_CHROMA_FILTER_OFFSET,
    .pad_b_simd = 3 // Three rows for AVX2
  };

  // Initialize separately. Gets rid of warning
  // about using nonstandard extension.
  epol_args.buf = ext_buffer;
  epol_args.ext = &ext;
  epol_args.ext_origin = &ext_origin;
  epol_args.ext_s = &ext_s;

  kvz_get_extended_block(&epol_args);
  kvz_sample_octpel_chroma_hi(state->encoder_control,
    ext_origin,
    ext_s,
    pb_w,
    pb_h,
    out->u,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);

  // Chroma V
  epol_args.src = ref->v;
  kvz_get_extended_block(&epol_args);
  kvz_sample_octpel_chroma_hi(state->encoder_control,
    ext_origin,
    ext_s,
    pb_w,
    pb_h,
    out->v,
    out_stride,
    mv_frac_x,
    mv_frac_y,
    mv_param);
}


/**
* \brief Copy from frame with extended border.
*
* \param ref_buf      pointer to the start of ref buffer
* \param ref_stride   stride of ref buffer
* \param ref_width    width of frame
* \param ref_height   height of frame
* \param rec_buf      pointer to the start of pu in rec buffer
* \param rec_stride   stride of rec buffer
* \param width        width of copied block
* \param height       height of copied block
* \param mv_in_frame  coordinates of copied block in frame coordinates
*/
static void inter_cp_with_ext_border(const kvz_pixel *ref_buf, int ref_stride,
                                     int ref_width, int ref_height,
                                     kvz_pixel *rec_buf, int rec_stride,
                                     int width, int height,
                                     const vector2d_t *mv_in_frame)
{
  for (int y = mv_in_frame->y; y < mv_in_frame->y + height; ++y) {
    for (int x = mv_in_frame->x; x < mv_in_frame->x + width; ++x) {
      vector2d_t in_frame = {
        CLIP(0, ref_width - 1, x),
        CLIP(0, ref_height - 1, y),
      };
      vector2d_t in_pu = {
        x - mv_in_frame->x,
        y - mv_in_frame->y,
      };
      int pu_index = in_pu.y * rec_stride + in_pu.x;
      int frame_index = in_frame.y * ref_stride + in_frame.x;
      rec_buf[pu_index] = ref_buf[frame_index];
    }
  }
}


/**
 * \brief Reconstruct an inter PU using uniprediction.
 *
 * \param state          encoder state
 * \param ref            picture to copy the data from
 * \param pu_x           PU x position
 * \param pu_y           PU y position
 * \param width          PU width
 * \param height         PU height
 * \param mv_param       motion vector
 * \param yuv_px         destination buffer for pixel precision
 * \param yuv_im         destination buffer for high-precision, or NULL if not needed
 * \param predict_luma   Enable or disable luma prediction for this call.
 * \param predict_chroma Enable or disable chroma prediction for this call.
*/
static unsigned inter_recon_unipred(const encoder_state_t * const state,
                                    const kvz_picture * const ref,
                                    int32_t pu_x,
                                    int32_t pu_y,
                                    int32_t pu_w,
                                    int32_t pu_h,
                                    int32_t out_stride_luma,
                                    const int16_t mv_param[2],
                                    yuv_t *yuv_px,
                                    yuv_im_t *yuv_im,
                                    bool predict_luma,
                                    bool predict_chroma)
{
  const vector2d_t int_mv = { mv_param[0] >> 2, mv_param[1] >> 2 };
  const vector2d_t int_mv_in_frame = {
    int_mv.x + pu_x + state->tile->offset_x,
    int_mv.y + pu_y + state->tile->offset_y
  };

  const bool int_mv_outside_frame = int_mv_in_frame.x < 0 ||
    int_mv_in_frame.y < 0 ||
    int_mv_in_frame.x + pu_w > ref->width ||
    int_mv_in_frame.y + pu_h > ref->height;

  // With 420, odd coordinates need interpolation.
  const bool fractional_chroma = (int_mv.x & 1) || (int_mv.y & 1);
  const bool fractional_luma = (mv_param[0] & 3) || (mv_param[1] & 3);

  // Generate prediction for luma.
  if (predict_luma) {
    if (fractional_luma) {
      // With a fractional MV, do interpolation.
      if (state->encoder_control->cfg.bipred && yuv_im) {
        inter_recon_frac_luma_hi(state, ref,
          pu_x, pu_y,
          pu_w, pu_h,
          mv_param, yuv_im, out_stride_luma);
      }
      else {
        inter_recon_frac_luma(state, ref,
          pu_x, pu_y,
          pu_w, pu_h,
          mv_param, yuv_px, out_stride_luma);
      }
    }
    else {
      // With an integer MV, copy pixels directly from the reference.
      if (int_mv_outside_frame) {
        inter_cp_with_ext_border(ref->y, ref->width,
          ref->width, ref->height,
          yuv_px->y, out_stride_luma,
          pu_w, pu_h,
          &int_mv_in_frame);
      }
      else {
        const int frame_mv_index = int_mv_in_frame.y * ref->width + int_mv_in_frame.x;
        kvz_pixels_blit(&ref->y[frame_mv_index],
          yuv_px->y,
          pu_w, pu_h,
          ref->width, out_stride_luma);
      }
    }
  }

  if (!predict_chroma) {
    return fractional_luma;
  }

  const unsigned out_stride_c = out_stride_luma / 2;

  // Generate prediction for chroma.
  if (fractional_luma || fractional_chroma) {
    // With a fractional MV, do interpolation.
    if (state->encoder_control->cfg.bipred && yuv_im) {
      inter_recon_frac_chroma_hi(state, ref,
                                    pu_x, pu_y,
                                    pu_w, pu_h, 
                                    mv_param, yuv_im, out_stride_c);
    } else {
      inter_recon_frac_chroma(state, ref,
                              pu_x, pu_y,
                              pu_w, pu_h,
                              mv_param, yuv_px, out_stride_c);
    }
  } else {
    // With an integer MV, copy pixels directly from the reference.
    const vector2d_t int_mv_in_frame_c = { int_mv_in_frame.x / 2, int_mv_in_frame.y / 2 };

    if (int_mv_outside_frame) {
      inter_cp_with_ext_border(ref->u, ref->width / 2,
                               ref->width / 2, ref->height / 2,
                               yuv_px->u, out_stride_c,
                               pu_w / 2, pu_h / 2,
                               &int_mv_in_frame_c);
      inter_cp_with_ext_border(ref->v, ref->width / 2,
                               ref->width / 2, ref->height / 2,
                               yuv_px->v, out_stride_c,
                               pu_w / 2, pu_h / 2,
                               &int_mv_in_frame_c);
    } else {
      const int frame_mv_index = int_mv_in_frame_c.y * ref->width / 2 + int_mv_in_frame_c.x;

      kvz_pixels_blit(&ref->u[frame_mv_index],
                      yuv_px->u,
                      pu_w / 2, pu_h / 2,
                      ref->width / 2, out_stride_c);
      kvz_pixels_blit(&ref->v[frame_mv_index],
                      yuv_px->v,
                      pu_w / 2, pu_h / 2,
                      ref->width / 2, out_stride_c);
    }
  }

  return fractional_luma | ((fractional_luma || fractional_chroma) << 1);
}
/**
 * \brief Reconstruct bi-pred inter PU
 *
 * \param state          encoder state
 * \param ref1           reference picture to copy the data from
 * \param ref2           other reference picture to copy the data from
 * \param pu_x           PU x position
 * \param pu_y           PU y position
 * \param width          PU width
 * \param height         PU height
 * \param mv_param       motion vectors
 * \param lcu            destination lcu
 * \param predict_luma   Enable or disable luma prediction for this call.
 * \param predict_chroma Enable or disable chroma prediction for this call.
 */
void kvz_inter_recon_bipred(const encoder_state_t *const state,
  const kvz_picture *ref1,
  const kvz_picture *ref2,
  int32_t pu_x,
  int32_t pu_y,
  int32_t pu_w,
  int32_t pu_h,
  int16_t mv_param[2][2],
  lcu_t *lcu,
  bool predict_luma,
  bool predict_chroma)
{
  // Allocate maximum size arrays for interpolated and copied samples
  ALIGNED(64) kvz_pixel px_buf_L0[LCU_LUMA_SIZE + 2 * LCU_CHROMA_SIZE];
  ALIGNED(64) kvz_pixel px_buf_L1[LCU_LUMA_SIZE + 2 * LCU_CHROMA_SIZE];
  ALIGNED(64) kvz_pixel_im im_buf_L0[LCU_LUMA_SIZE + 2 * LCU_CHROMA_SIZE];
  ALIGNED(64) kvz_pixel_im im_buf_L1[LCU_LUMA_SIZE + 2 * LCU_CHROMA_SIZE];

  yuv_t px_L0;
  px_L0.size = pu_w * pu_h;
  px_L0.y = &px_buf_L0[0];
  px_L0.u = &px_buf_L0[LCU_LUMA_SIZE];
  px_L0.v = &px_buf_L0[LCU_LUMA_SIZE + LCU_CHROMA_SIZE];

  yuv_t px_L1;
  px_L1.size = pu_w * pu_h;
  px_L1.y = &px_buf_L1[0];
  px_L1.u = &px_buf_L1[LCU_LUMA_SIZE];
  px_L1.v = &px_buf_L1[LCU_LUMA_SIZE + LCU_CHROMA_SIZE];

  yuv_im_t im_L0;
  im_L0.size = pu_w * pu_h;
  im_L0.y = &im_buf_L0[0];
  im_L0.u = &im_buf_L0[LCU_LUMA_SIZE];
  im_L0.v = &im_buf_L0[LCU_LUMA_SIZE + LCU_CHROMA_SIZE];

  yuv_im_t im_L1;
  im_L1.size = pu_w * pu_h;
  im_L1.y = &im_buf_L1[0];
  im_L1.u = &im_buf_L1[LCU_LUMA_SIZE];
  im_L1.v = &im_buf_L1[LCU_LUMA_SIZE + LCU_CHROMA_SIZE];

  // Sample blocks from both reference picture lists.
  // Flags state if the outputs were written to high-precision / interpolated sample buffers.
  unsigned im_flags_L0 = inter_recon_unipred(state, ref1, pu_x, pu_y, pu_w, pu_h, pu_w, mv_param[0],
                                             &px_L0, &im_L0, predict_luma, predict_chroma);
  unsigned im_flags_L1 = inter_recon_unipred(state, ref2, pu_x, pu_y, pu_w, pu_h, pu_w, mv_param[1],
                                             &px_L1, &im_L1, predict_luma, predict_chroma);

  // After reconstruction, merge the predictors by taking an average of each pixel
  kvz_bipred_average(lcu, &px_L0, &px_L1, &im_L0, &im_L1,
                     pu_x, pu_y, pu_w, pu_h,
                     im_flags_L0, im_flags_L1,
                     predict_luma, predict_chroma);
}


/**
 * Reconstruct a single CU.
 *
 * The CU may consist of multiple PUs, each of which can use either
 * uniprediction or biprediction.
 *
 * \param state   encoder state
 * \param lcu     containing LCU
 * \param x       x-coordinate of the CU in pixels
 * \param y       y-coordinate of the CU in pixels
 * \param width   CU width
 * \param predict_luma   Enable or disable luma prediction for this call.
 * \param predict_chroma Enable or disable chroma prediction for this call.
 */
void kvz_inter_recon_cu(const encoder_state_t * const state,
                        lcu_t *lcu,
                        int32_t x,
                        int32_t y,
                        int32_t width,
                        bool predict_luma,
                        bool predict_chroma)
{
  cu_info_t *cu = LCU_GET_CU_AT_PX(lcu, SUB_SCU(x), SUB_SCU(y));
  const int num_pu = kvz_part_mode_num_parts[cu->part_size];
  for (int i = 0; i < num_pu; ++i) {
    kvz_inter_pred_pu(state, lcu, x, y, width, predict_luma, predict_chroma, i);
  }
}

/**
 * Predict a single PU.
 *
 * The PU may use either uniprediction or biprediction.
 *
 * \param state          encoder state
 * \param lcu            containing LCU
 * \param x              x-coordinate of the CU in pixels
 * \param y              y-coordinate of the CU in pixels
 * \param width          CU width
 * \param predict_luma   Enable or disable luma prediction for this call.
 * \param predict_chroma Enable or disable chroma prediction for this call.
 * \param i_pu           Index of the PU. Always zero for 2Nx2N. Used for SMP+AMP.
 */
void kvz_inter_pred_pu(const encoder_state_t * const state,
                       lcu_t *lcu,
                       int32_t x,
                       int32_t y,
                       int32_t width,
                       bool predict_luma,
                       bool predict_chroma,
                       int i_pu)

{
  cu_info_t *cu = LCU_GET_CU_AT_PX(lcu, SUB_SCU(x), SUB_SCU(y));
  const int pu_x = PU_GET_X(cu->part_size, width, x, i_pu);
  const int pu_y = PU_GET_Y(cu->part_size, width, y, i_pu);
  const int pu_w = PU_GET_W(cu->part_size, width, i_pu);
  const int pu_h = PU_GET_H(cu->part_size, width, i_pu);
  cu_info_t *pu = LCU_GET_CU_AT_PX(lcu, SUB_SCU(pu_x), SUB_SCU(pu_y));

  if (pu->inter.mv_dir == 3) {
    const kvz_picture *const refs[2] = {
      state->frame->ref->images[
        state->frame->ref_LX[0][
          pu->inter.mv_ref[0]]],
      state->frame->ref->images[
        state->frame->ref_LX[1][
          pu->inter.mv_ref[1]]],
    };
    kvz_inter_recon_bipred(state,
      refs[0], refs[1],
      pu_x, pu_y,
      pu_w, pu_h,
      pu->inter.mv,
      lcu,
      predict_luma, predict_chroma);
  }
  else {
    const int mv_idx = pu->inter.mv_dir - 1;
    const kvz_picture *const ref =
      state->frame->ref->images[
        state->frame->ref_LX[mv_idx][
          pu->inter.mv_ref[mv_idx]]];

    const unsigned offset_luma = SUB_SCU(pu_y) * LCU_WIDTH + SUB_SCU(pu_x);
    const unsigned offset_chroma = SUB_SCU(pu_y) / 2 * LCU_WIDTH_C + SUB_SCU(pu_x) / 2;
    yuv_t lcu_adapter;
    lcu_adapter.size = pu_w * pu_h;
    lcu_adapter.y = lcu->rec.y + offset_luma,
    lcu_adapter.u = lcu->rec.u + offset_chroma,
    lcu_adapter.v = lcu->rec.v + offset_chroma,

    inter_recon_unipred(state,
      ref,
      pu_x, pu_y,
      pu_w, pu_h,
      LCU_WIDTH,
      pu->inter.mv[mv_idx],
      &lcu_adapter,
      NULL,
      predict_luma, predict_chroma);
  }
}

/**
 * \brief Clear unused L0/L1 motion vectors and reference
 * \param cu coding unit to clear
 */
static void inter_clear_cu_unused(cu_info_t* cu)
{
  for (unsigned i = 0; i < 2; ++i) {
    if (cu->inter.mv_dir & (1 << i)) continue;

    cu->inter.mv[i][0] = 0;
    cu->inter.mv[i][1] = 0;
    cu->inter.mv_ref[i] = 255;
  }
}

/**
 * \brief Check whether a0 mv cand block is coded before the current block.
 * \param x       x-coordinate of the current block (in pixels)
 * \param y       y-coordinate of the current block (in pixels)
 * \param width   width of the current block (in pixels)
 * \param height  height of the current block (in pixels)
 * \return        True, if the a0 mv candidate block is coded before the
 *                current block. Otherwise false.
 */
static bool is_a0_cand_coded(int x, int y, int width, int height)
{
  int size = MIN(width & ~(width - 1), height & ~(height - 1));

  if (height != size) {
    // For SMP and AMP blocks the situation is equivalent to a square block
    // at the lower left corner of the PU.
    y = y + height - size;
  }

  while (size < LCU_WIDTH) {
    const int parent_size = 2 * size;
    const int cu_index    = (x % parent_size != 0) + 2 * (y % parent_size != 0);
    switch (cu_index) {
      case 0:
        // A0 is in the CU directly left of the parent CU so it has been
        // coded already.
        //    +---+---+
        //    | X |   |
        //    |---+---+
        // A0 |   |   |
        //    +---+---+
        return true;

      case 1:
        // A0 is in the CU that will be coded after the current CU.
        //    +---+---+
        //    |   | X |
        //    |---+---+
        //    |A0 |   |
        //    +---+---+
        return false;

      case 2:
        //    +---+---+
        //    |   |   |
        //    |---+---+
        //    | X |   |
        //    +---+---+
        // A0

        // Move to the parent block.
        y -= size;
        size = parent_size;
        break;

      case 3:
        // A0 is in the CU directly down of the parent CU so is has not
        // been coded yet.
        //    +---+---+
        //    |   |   |
        //    |---+---+
        //    |   | X |
        //    +---+---+
        //     A0
        return false;
    }
  }

  // For 64x64 blocks A0 candidate is located outside the LCU.
  return false;
}

/**
 * \brief Check whether b0 mv cand block is coded before the current block.
 * \param x       x-coordinate of the current block (in pixels)
 * \param y       y-coordinate of the current block (in pixels)
 * \param width   width of the current block (in pixels)
 * \param height  height of the current block (in pixels)
 * \return        True, if the b0 mv candidate block is coded before the
 *                current block. Otherwise false.
 */
static bool is_b0_cand_coded(int x, int y, int width, int height)
{
  int size = MIN(width & ~(width - 1), height & ~(height - 1));

  if (width != size) {
    // For SMP and AMP blocks the situation is equivalent to a square block
    // at the upper right corner of the PU.
    x = x + width - size;
  }

  while (size < LCU_WIDTH) {
    const int parent_size = 2 * size;
    const int cu_index    = (x % parent_size != 0) + 2 * (y % parent_size != 0);
    switch (cu_index) {
      case 0:
        // B0 is in the CU directly above the parent CU so it has been
        // coded already.
        //         B0
        //    +---+---+
        //    | X |   |
        //    |---+---+
        //    |   |   |
        //    +---+---+
        return true;

      case 1:
        //             B0
        //    +---+---+
        //    |   | X |
        //    |---+---+
        //    |   |   |
        //    +---+---+

        // Move to the parent block.
        x -= size;
        size = parent_size;
        break;

      case 2:
        //    +---+---+
        //    |   |B0 |
        //    |---+---+
        //    | X |   |
        //    +---+---+
        return true;

      case 3:
        // B0 is in the CU directly right of the parent CU so is has not
        // been coded yet.
        //    +---+---+
        //    |   |   | B0
        //    |---+---+
        //    |   | X |
        //    +---+---+
        return false;
    }
  }

  // The LCU to the right and up of the current LCU has been coded already.
  return true;
}


/**
 * \brief Get merge candidates for current block
 *
 * \param state     encoder control state to use
 * \param x         block x position in SCU
 * \param y         block y position in SCU
 * \param width     current block width
 * \param height    current block height
 * \param ref_list  which reference list, L0 is 1 and L1 is 2
 * \param ref_idx   index in the reference list
 * \param cand_out  will be filled with C3 and H candidates
 */
static void get_temporal_merge_candidates(const encoder_state_t * const state,
                                          int32_t x,
                                          int32_t y,
                                          int32_t width,
                                          int32_t height,
                                          uint8_t ref_list,
                                          uint8_t ref_idx,
                                          merge_candidates_t *cand_out)
{
  /*
  Predictor block locations
  _________
  |CurrentPU|
  | |C0|__  |
  |    |C3| |
  |_________|_
            |H|
  */

  cand_out->c3 = cand_out->h = NULL;

  // Find temporal reference
  if (state->frame->ref->used_size) {
    uint32_t colocated_ref;

    // Select L0/L1 ref_idx reference
    if (state->frame->ref_LX_size[ref_list-1] > ref_idx) {
      colocated_ref = state->frame->ref_LX[ref_list - 1][ref_idx];
    } else {
      // not found
      return;
    }

    cu_array_t *ref_cu_array = state->frame->ref->cu_arrays[colocated_ref];
    int cu_per_width = ref_cu_array->width / SCU_WIDTH;

    uint32_t xColBr = x + width;
    uint32_t yColBr = y + height;

    // H must be available
    if (xColBr < state->encoder_control->in.width &&
        yColBr < state->encoder_control->in.height) {
      int32_t H_offset = -1;

      // Y inside the current CTU / LCU
      if (yColBr % LCU_WIDTH != 0) {
        H_offset = ((xColBr >> 4) << 4) / SCU_WIDTH +
                  (((yColBr >> 4) << 4) / SCU_WIDTH) * cu_per_width;
      }

      if (H_offset >= 0) {
        // Only use when it's inter block
        if (ref_cu_array->data[H_offset].type == CU_INTER) {
          cand_out->h = &ref_cu_array->data[H_offset];
        }
      }
    }
    uint32_t xColCtr = x + (width / 2);
    uint32_t yColCtr = y + (height / 2);

    // C3 must be inside the LCU, in the center position of current CU
    if (xColCtr < state->encoder_control->in.width && yColCtr < state->encoder_control->in.height) {
      uint32_t C3_offset = ((xColCtr >> 4) << 4) / SCU_WIDTH + ((((yColCtr >> 4) << 4) / SCU_WIDTH) * cu_per_width);
      if (ref_cu_array->data[C3_offset].type == CU_INTER) {
        cand_out->c3 = &ref_cu_array->data[C3_offset];
      }
    }
  }
}

/**
 * \brief Get merge candidates for current block.
 *
 * The output parameters b0, b1, b2, a0, a1 are pointed to the
 * corresponding cu_info_t struct in lcu->cu, or set to NULL, if the
 * candidate is not available.
 *
 * \param x               block x position in pixels
 * \param y               block y position in pixels
 * \param width           block width in pixels
 * \param height          block height in pixels
 * \param picture_width   tile width in pixels
 * \param picture_height  tile height in pixels
 * \param lcu             current LCU
 * \param cand_out        will be filled with A and B candidates
 */
static void get_spatial_merge_candidates(int32_t x,
                                         int32_t y,
                                         int32_t width,
                                         int32_t height,
                                         int32_t picture_width,
                                         int32_t picture_height,
                                         lcu_t *lcu,
                                         merge_candidates_t *cand_out)
{
  /*
  Predictor block locations
  ____      _______
  |B2|______|B1|B0|
     |         |
     |  Cur CU |
   __|         |
  |A1|_________|
  |A0|
  */
  int32_t x_local = SUB_SCU(x); //!< coordinates from top-left of this LCU
  int32_t y_local = SUB_SCU(y);
  // A0 and A1 availability testing
  if (x != 0) {
    cu_info_t *a1 = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local + height - 1);
    // Do not check a1->coded because the block above is always coded before
    // the current one and the flag is not set when searching an SMP block.
    if (a1->type == CU_INTER) {
      inter_clear_cu_unused(a1);
      cand_out->a[1] = a1;
    }

    if (y_local + height < LCU_WIDTH && y + height < picture_height) {
      cu_info_t *a0 = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local + height);
      if (a0->type == CU_INTER && is_a0_cand_coded(x, y, width, height)) {
        inter_clear_cu_unused(a0);
        cand_out->a[0] = a0;
      }
    }
  }

  // B0, B1 and B2 availability testing
  if (y != 0) {
    cu_info_t *b0 = NULL;
    if (x + width < picture_width) {
      if (x_local + width < LCU_WIDTH) {
        b0 = LCU_GET_CU_AT_PX(lcu, x_local + width, y_local - 1);
      } else if (y_local == 0) {
        // Special case, top-right CU
        b0 = LCU_GET_TOP_RIGHT_CU(lcu);
      }
    }
    if (b0 && b0->type == CU_INTER && is_b0_cand_coded(x, y, width, height)) {
      inter_clear_cu_unused(b0);
      cand_out->b[0] = b0;
    }

    cu_info_t *b1 = LCU_GET_CU_AT_PX(lcu, x_local + width - 1, y_local - 1);
    // Do not check b1->coded because the block to the left is always coded
    // before the current one and the flag is not set when searching an SMP
    // block.
    if (b1->type == CU_INTER) {
      inter_clear_cu_unused(b1);
      cand_out->b[1] = b1;
    }

    if (x != 0) {
      cu_info_t *b2 = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local - 1);
      // Do not check b2->coded because the block above and to the left is
      // always coded before the current one.
      if (b2->type == CU_INTER) {
        inter_clear_cu_unused(b2);
        cand_out->b[2] = b2;
      }
    }
  }
}

/**
 * \brief Get merge candidates for current block.
 *
 * The output parameters b0, b1, b2, a0, a1 are pointed to the
 * corresponding cu_info_t struct in lcu->cu, or set to NULL, if the
 * candidate is not available.
 *
 * \param cua             cu information
 * \param x               block x position in pixels
 * \param y               block y position in pixels
 * \param width           block width in pixels
 * \param height          block height in pixels
 * \param picture_width   tile width in pixels
 * \param picture_height  tile height in pixels
 * \param cand_out        will be filled with A and B candidates
 */
static void get_spatial_merge_candidates_cua(const cu_array_t *cua,
                                             int32_t x,
                                             int32_t y,
                                             int32_t width,
                                             int32_t height,
                                             int32_t picture_width,
                                             int32_t picture_height,
                                             merge_candidates_t *cand_out)
{
  /*
  Predictor block locations
  ____      _______
  |B2|______|B1|B0|
     |         |
     |  Cur CU |
   __|         |
  |A1|_________|
  |A0|
  */
  int32_t x_local = SUB_SCU(x); //!< coordinates from top-left of this LCU
  int32_t y_local = SUB_SCU(y);
  // A0 and A1 availability testing
  if (x != 0) {
    const cu_info_t *a1 = kvz_cu_array_at_const(cua, x - 1, y + height - 1);
    // The block above is always coded before the current one.
    if (a1->type == CU_INTER) {
      cand_out->a[1] = a1;
    }

    if (y_local + height < LCU_WIDTH && y + height < picture_height) {
      const cu_info_t *a0 = kvz_cu_array_at_const(cua, x - 1, y + height);
      if (a0->type == CU_INTER && is_a0_cand_coded(x, y, width, height)) {
        cand_out->a[0] = a0;
      }
    }
  }

  // B0, B1 and B2 availability testing
  if (y != 0) {
    if (x + width < picture_width && (x_local + width < LCU_WIDTH || y_local == 0)) {
      const cu_info_t *b0 = kvz_cu_array_at_const(cua, x + width, y - 1);
      if (b0->type == CU_INTER && is_b0_cand_coded(x, y, width, height)) {
        cand_out->b[0] = b0;
      }
    }

    const cu_info_t *b1 = kvz_cu_array_at_const(cua, x + width - 1, y - 1);
    // The block to the left is always coded before the current one.
    if (b1->type == CU_INTER) {
      cand_out->b[1] = b1;
    }

    if (x != 0) {
      const cu_info_t *b2 = kvz_cu_array_at_const(cua, x - 1, y - 1);
      // The block above and to the left is always coded before the current
      // one.
      if (b2->type == CU_INTER) {
        cand_out->b[2] = b2;
      }
    }
  }
}

static INLINE int16_t get_scaled_mv(int16_t mv, int scale)
{
  int32_t scaled = scale * mv;
  return CLIP(-32768, 32767, (scaled + 127 + (scaled < 0)) >> 8);
}

static void apply_mv_scaling_pocs(int32_t current_poc,
                                  int32_t current_ref_poc,
                                  int32_t neighbor_poc,
                                  int32_t neighbor_ref_poc,
                                  int16_t mv_cand[2])
{
  int32_t diff_current  = current_poc  - current_ref_poc;
  int32_t diff_neighbor = neighbor_poc - neighbor_ref_poc;

  if (diff_current == diff_neighbor) return;

  diff_current  = CLIP(-128, 127, diff_current);
  diff_neighbor = CLIP(-128, 127, diff_neighbor);

  int scale = CLIP(-4096, 4095,
    (diff_current * ((0x4000 + (abs(diff_neighbor) >> 1)) / diff_neighbor) + 32) >> 6);

  mv_cand[0] = get_scaled_mv(mv_cand[0], scale);
  mv_cand[1] = get_scaled_mv(mv_cand[1], scale);
}

static INLINE void apply_mv_scaling(const encoder_state_t *state,
                                    const cu_info_t *current_cu,
                                    const cu_info_t *neighbor_cu,
                                    int8_t current_reflist,
                                    int8_t neighbor_reflist,
                                    int16_t mv_cand[2])
{
  apply_mv_scaling_pocs(state->frame->poc,
                        state->frame->ref->pocs[
                          state->frame->ref_LX[current_reflist][
                            current_cu->inter.mv_ref[current_reflist]]],
                        state->frame->poc,
                        state->frame->ref->pocs[
                          state->frame->ref_LX[neighbor_reflist][
                            neighbor_cu->inter.mv_ref[neighbor_reflist]]],
                        mv_cand);
}

/**
 * \brief Try to add a temporal MVP or merge candidate.
 *
 * \param state         encoder state
 * \param current_ref   index of the picture referenced by the current CU
 * \param colocated     colocated CU
 * \param reflist       either 0 (for L0) or 1 (for L1)
 * \param[out] mv_out   Returns the motion vector
 *
 * \return Whether a temporal candidate was added or not.
 */
static bool add_temporal_candidate(const encoder_state_t *state,
                                   uint8_t current_ref,
                                   const cu_info_t *colocated,
                                   int32_t reflist,
                                   int16_t mv_out[2])
{
  if (!colocated) return false;

  int colocated_ref;
  if (state->frame->ref_LX_size[0] > 0) {
    // get the first reference from L0 if it exists
    colocated_ref = state->frame->ref_LX[0][0];
  } else {
    // otherwise no candidate added
    return false;
  }

  // When there are reference pictures from the future (POC > current POC)
  // in L0 or L1, the primary list for the colocated PU is the inverse of
  // collocated_from_l0_flag. Otherwise it is equal to reflist.
  //
  // Kvazaar always sets collocated_from_l0_flag so the list is L1 when
  // there are future references.
  int col_list = reflist;
  for (int i = 0; i < state->frame->ref->used_size; i++) {
    if (state->frame->ref->pocs[i] > state->frame->poc) {
      col_list = 1;
      break;
    }
  }

  if ((colocated->inter.mv_dir & (col_list + 1)) == 0) {
    // Use the other list if the colocated PU does not have a MV for the
    // primary list.
    col_list = 1 - col_list;
  }

  mv_out[0] = colocated->inter.mv[col_list][0];
  mv_out[1] = colocated->inter.mv[col_list][1];
  apply_mv_scaling_pocs(
    state->frame->poc,
    state->frame->ref->pocs[current_ref],
    state->frame->ref->pocs[colocated_ref],
    state->frame->ref->images[colocated_ref]->ref_pocs[
      state->frame->ref->ref_LXs[colocated_ref]
        [col_list][colocated->inter.mv_ref[col_list]]],
    mv_out
  );

  return true;
}

static INLINE bool add_mvp_candidate(const encoder_state_t *state,
                                     const cu_info_t *cur_cu,
                                     const cu_info_t *cand,
                                     int8_t reflist,
                                     bool scaling,
                                     int16_t mv_cand_out[2])
{
  if (!cand) return false;

  assert(cand->inter.mv_dir != 0);

  for (int i = 0; i < 2; i++) {
    const int cand_list = i == 0 ? reflist : !reflist;

    if ((cand->inter.mv_dir & (1 << cand_list)) == 0) continue;

    if (scaling) {
      mv_cand_out[0] = cand->inter.mv[cand_list][0];
      mv_cand_out[1] = cand->inter.mv[cand_list][1];
      apply_mv_scaling(state, cur_cu, cand, reflist, cand_list, mv_cand_out);
      return true;
    }

    if (cand->inter.mv_dir & (1 << cand_list) &&
        state->frame->ref_LX[cand_list][cand->inter.mv_ref[cand_list]] ==
        state->frame->ref_LX[reflist][cur_cu->inter.mv_ref[reflist]])
    {
      mv_cand_out[0] = cand->inter.mv[cand_list][0];
      mv_cand_out[1] = cand->inter.mv[cand_list][1];
      return true;
    }
  }

  return false;
}

/**
 * \brief Pick two mv candidates from the spatial and temporal candidates.
 */
static void get_mv_cand_from_candidates(const encoder_state_t * const state,
                                        int32_t x,
                                        int32_t y,
                                        int32_t width,
                                        int32_t height,
                                        const merge_candidates_t *merge_cand,
                                        const cu_info_t * const cur_cu,
                                        int8_t reflist,
                                        int16_t mv_cand[2][2])
{
  const cu_info_t *const *a = merge_cand->a;
  const cu_info_t *const *b = merge_cand->b;
  const cu_info_t *c3 = merge_cand->c3;
  const cu_info_t *h  = merge_cand->h;

  uint8_t candidates = 0;
  uint8_t b_candidates = 0;

  // Left predictors without scaling
  for (int i = 0; i < 2; i++) {
    if (add_mvp_candidate(state, cur_cu, a[i], reflist, false, mv_cand[candidates])) {
      candidates++;
      break;
    }
  }

  // Left predictors with scaling
  if (candidates == 0) {
    for (int i = 0; i < 2; i++) {
      if (add_mvp_candidate(state, cur_cu, a[i], reflist, true, mv_cand[candidates])) {
        candidates++;
        break;
      }
    }
  }

  // Top predictors without scaling
  for (int i = 0; i < 3; i++) {
    if (add_mvp_candidate(state, cur_cu, b[i], reflist, false, mv_cand[candidates])) {
      b_candidates++;
      break;
    }
  }

  candidates += b_candidates;

  // When a1 or a0 is available, we dont check for secondary B candidates.
  if (a[0] || a[1]) {
    b_candidates = 1;
  } else if (candidates != 2) {
    b_candidates = 0;
  }

  if (!b_candidates) {
    // Top predictors with scaling
    for (int i = 0; i < 3; i++) {
      if (add_mvp_candidate(state, cur_cu, b[i], reflist, true, mv_cand[candidates])) {
        candidates++;
        break;
      }
    }
  }

  // Remove identical candidate
  if (candidates == 2 && mv_cand[0][0] == mv_cand[1][0] && mv_cand[0][1] == mv_cand[1][1]) {
    candidates = 1;
  }

  // Use Temporal Motion Vector Prediction when enabled.
  // TMVP required at least two sequential P/B-frames.
  bool can_use_tmvp =
    state->encoder_control->cfg.tmvp_enable &&
    state->frame->poc > 1 &&
    state->frame->ref->used_size &&
    candidates < AMVP_MAX_NUM_CANDS &&
    (h != NULL || c3 != NULL);

  if (can_use_tmvp && add_temporal_candidate(state,
                                             state->frame->ref_LX[reflist][cur_cu->inter.mv_ref[reflist]],
                                             (h != NULL) ? h : c3,
                                             reflist,
                                             mv_cand[candidates]))
  {
    candidates++;
  }

  // Fill with (0,0)
  while (candidates < AMVP_MAX_NUM_CANDS) {
    mv_cand[candidates][0] = 0;
    mv_cand[candidates][1] = 0;
    candidates++;
  }
}

/**
 * \brief Get MV prediction for current block.
 *
 * \param state     encoder state
 * \param x         block x position in pixels
 * \param y         block y position in pixels
 * \param width     block width in pixels
 * \param height    block height in pixels
 * \param mv_cand   Return the motion vector candidates.
 * \param cur_cu    current CU
 * \param lcu       current LCU
 * \param reflist   reflist index (either 0 or 1)
 */
void kvz_inter_get_mv_cand(const encoder_state_t * const state,
                           int32_t x,
                           int32_t y,
                           int32_t width,
                           int32_t height,
                           int16_t mv_cand[2][2],
                           const cu_info_t  * const cur_cu,
                           lcu_t *lcu,
                           int8_t reflist)
{
  merge_candidates_t merge_cand = { {0, 0}, {0, 0, 0}, 0, 0 };

  get_spatial_merge_candidates(x, y, width, height,
                               state->tile->frame->width,
                               state->tile->frame->height,
                               lcu,
                               &merge_cand);
  get_temporal_merge_candidates(state, x, y, width, height, 1, 0, &merge_cand);
  get_mv_cand_from_candidates(state, x, y, width, height, &merge_cand, cur_cu, reflist, mv_cand);
}

/**
 * \brief Get MV prediction for current block using state->tile->frame->cu_array.
 *
 * \param state     encoder state
 * \param x         block x position in pixels
 * \param y         block y position in pixels
 * \param width     block width in pixels
 * \param height    block height in pixels
 * \param mv_cand   Return the motion vector candidates.
 * \param cur_cu    current CU
 * \param reflist   reflist index (either 0 or 1)
 */
void kvz_inter_get_mv_cand_cua(const encoder_state_t * const state,
                               int32_t x,
                               int32_t y,
                               int32_t width,
                               int32_t height,
                               int16_t mv_cand[2][2],
                               const cu_info_t* cur_cu,
                               int8_t reflist)
{
  merge_candidates_t merge_cand = { {0, 0}, {0, 0, 0}, 0, 0 };

  const cu_array_t *cua = state->tile->frame->cu_array;
  get_spatial_merge_candidates_cua(cua,
                                   x, y, width, height,
                                   state->tile->frame->width, state->tile->frame->height,
                                   &merge_cand);
  get_temporal_merge_candidates(state, x, y, width, height, 1, 0, &merge_cand);
  get_mv_cand_from_candidates(state, x, y, width, height, &merge_cand, cur_cu, reflist, mv_cand);
}

static bool is_duplicate_candidate(const cu_info_t* cu1, const cu_info_t* cu2)
{
  if (!cu2) return false;
  if (cu1->inter.mv_dir != cu2->inter.mv_dir) return false;

  for (int reflist = 0; reflist < 2; reflist++) {
    if (cu1->inter.mv_dir & (1 << reflist)) {
      if (cu1->inter.mv[reflist][0]  != cu2->inter.mv[reflist][0]  ||
          cu1->inter.mv[reflist][1]  != cu2->inter.mv[reflist][1]  ||
          cu1->inter.mv_ref[reflist] != cu2->inter.mv_ref[reflist]) {
        return false;
      }
    }
  }

  return true;
}

static bool add_merge_candidate(const cu_info_t *cand,
                                const cu_info_t *possible_duplicate1,
                                const cu_info_t *possible_duplicate2,
                                inter_merge_cand_t *merge_cand_out,
                                uint8_t candidates,
                                uint8_t max_num_cands)
{
  if (!cand ||
      is_duplicate_candidate(cand, possible_duplicate1) ||
      is_duplicate_candidate(cand, possible_duplicate2) ||
      candidates >= max_num_cands) {
    return false;
  }

  merge_cand_out->mv[0][0] = cand->inter.mv[0][0];
  merge_cand_out->mv[0][1] = cand->inter.mv[0][1];
  merge_cand_out->mv[1][0] = cand->inter.mv[1][0];
  merge_cand_out->mv[1][1] = cand->inter.mv[1][1];
  merge_cand_out->ref[0]   = cand->inter.mv_ref[0]; // L0/L1 references
  merge_cand_out->ref[1]   = cand->inter.mv_ref[1];
  merge_cand_out->dir      = cand->inter.mv_dir;
  return true;
}

/**
 * \brief Get merge predictions for current block
 * \param state     the encoder state
 * \param x         block x position in SCU
 * \param y         block y position in SCU
 * \param width     block width
 * \param height    block height
 * \param use_a1    true, if candidate a1 can be used
 * \param use_b1    true, if candidate b1 can be used
 * \param mv_cand   Returns the merge candidates.
 * \param lcu       lcu containing the block
 * \return          number of merge candidates
 */
uint8_t kvz_inter_get_merge_cand(const encoder_state_t * const state,
                                 int32_t x, int32_t y,
                                 int32_t width, int32_t height,
                                 bool use_a1, bool use_b1,
                                 inter_merge_cand_t mv_cand[MRG_MAX_NUM_CANDS],
                                 lcu_t *lcu)
{
  uint8_t candidates = 0;
  int8_t zero_idx = 0;

  merge_candidates_t merge_cand = { {0, 0}, {0, 0, 0}, 0, 0 };
  const uint8_t max_num_cands = state->encoder_control->cfg.max_merge;
  get_spatial_merge_candidates(x, y, width, height,
                               state->tile->frame->width,
                               state->tile->frame->height,
                               lcu,
                               &merge_cand);

  const cu_info_t **a = merge_cand.a;
  const cu_info_t **b = merge_cand.b;

  if (!use_a1) a[1] = NULL;
  if (!use_b1) b[1] = NULL;

  if (add_merge_candidate(a[1], NULL, NULL, &mv_cand[candidates], candidates, max_num_cands)) candidates++;
  if (add_merge_candidate(b[1], a[1], NULL, &mv_cand[candidates], candidates, max_num_cands)) candidates++;
  if (add_merge_candidate(b[0], b[1], NULL, &mv_cand[candidates], candidates, max_num_cands)) candidates++;
  if (add_merge_candidate(a[0], a[1], NULL, &mv_cand[candidates], candidates, max_num_cands)) candidates++;
  if (candidates < 4 &&
      add_merge_candidate(b[2], a[1], b[1], &mv_cand[candidates], candidates, max_num_cands)) candidates++;

  bool can_use_tmvp =
    state->encoder_control->cfg.tmvp_enable &&
    candidates < max_num_cands &&
    state->frame->ref->used_size;

  if (can_use_tmvp) {
    mv_cand[candidates].dir = 0;

    const int max_reflist = (state->frame->slicetype == KVZ_SLICE_B ? 1 : 0);
    for (int reflist = 0; reflist <= max_reflist; reflist++) {
      // Fetch temporal candidates for the current CU
      get_temporal_merge_candidates(state, x, y, width, height, 1, 0, &merge_cand);
      // TODO: enable L1 TMVP candidate
      // get_temporal_merge_candidates(state, x, y, width, height, 2, 0, &merge_cand);

      const cu_info_t *temporal_cand =
        (merge_cand.h != NULL) ? merge_cand.h : merge_cand.c3;

      if (add_temporal_candidate(state,
                                 // Reference index 0 is always used for
                                 // the temporal merge candidate.
                                 state->frame->ref_LX[reflist][0],
                                 temporal_cand,
                                 reflist,
                                 mv_cand[candidates].mv[reflist])) {
        mv_cand[candidates].ref[reflist] = 0;
        mv_cand[candidates].dir |= (1 << reflist);
      }
    }

    if (mv_cand[candidates].dir != 0) candidates++;
  }

  if (candidates < max_num_cands && state->frame->slicetype == KVZ_SLICE_B) {
    #define NUM_PRIORITY_LIST 12;
    static const uint8_t priorityList0[] = { 0, 1, 0, 2, 1, 2, 0, 3, 1, 3, 2, 3 };
    static const uint8_t priorityList1[] = { 1, 0, 2, 0, 2, 1, 3, 0, 3, 1, 3, 2 };
    uint8_t cutoff = candidates;
    for (int32_t idx = 0; idx<cutoff*(cutoff - 1) && candidates != max_num_cands; idx++) {
      uint8_t i = priorityList0[idx];
      uint8_t j = priorityList1[idx];
      if (i >= candidates || j >= candidates) break;

      // Find one L0 and L1 candidate according to the priority list
      if ((mv_cand[i].dir & 0x1) && (mv_cand[j].dir & 0x2)) {
        mv_cand[candidates].dir = 3;

        // get Mv from cand[i] and cand[j]
        mv_cand[candidates].mv[0][0]  = mv_cand[i].mv[0][0];
        mv_cand[candidates].mv[0][1]  = mv_cand[i].mv[0][1];
        mv_cand[candidates].mv[1][0]  = mv_cand[j].mv[1][0];
        mv_cand[candidates].mv[1][1]  = mv_cand[j].mv[1][1];
        mv_cand[candidates].ref[0]    = mv_cand[i].ref[0];
        mv_cand[candidates].ref[1]    = mv_cand[j].ref[1];

        if (state->frame->ref_LX[0][mv_cand[i].ref[0]] ==
            state->frame->ref_LX[1][mv_cand[j].ref[1]]
            &&
            mv_cand[i].mv[0][0] == mv_cand[j].mv[1][0] && 
            mv_cand[i].mv[0][1] == mv_cand[j].mv[1][1]) {
          // Not a candidate
        } else {
          candidates++;
        }
      }
    }
  }

  int num_ref = state->frame->ref->used_size;

  if (candidates < max_num_cands && state->frame->slicetype == KVZ_SLICE_B) {
    int j;
    int ref_negative = 0;
    int ref_positive = 0;
    for (j = 0; j < state->frame->ref->used_size; j++) {
      if (state->frame->ref->pocs[j] < state->frame->poc) {
        ref_negative++;
      } else {
        ref_positive++;
      }
    }
    num_ref = MIN(ref_negative, ref_positive);
  }

  // Add (0,0) prediction
  while (candidates != max_num_cands) {
    mv_cand[candidates].mv[0][0] = 0;
    mv_cand[candidates].mv[0][1] = 0;
    mv_cand[candidates].ref[0] = (zero_idx >= num_ref - 1) ? 0 : zero_idx;
    mv_cand[candidates].ref[1] = mv_cand[candidates].ref[0];
    mv_cand[candidates].dir = 1;
    if (state->frame->slicetype == KVZ_SLICE_B) {
      mv_cand[candidates].mv[1][0] = 0;
      mv_cand[candidates].mv[1][1] = 0;
      mv_cand[candidates].dir = 3;
    }
    zero_idx++;
    candidates++;
  }

  return candidates;
}
