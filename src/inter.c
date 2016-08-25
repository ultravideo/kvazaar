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

#include "inter.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "encoder.h"
#include "imagelist.h"
#include "strategies/generic/picture-generic.h"
#include "strategies/strategies-ipol.h"
#include "videoframe.h"


static void inter_recon_frac_luma(const encoder_state_t * const state,
                                  const kvz_picture * const ref,
                                  int32_t xpos,
                                  int32_t ypos,
                                  int32_t block_width,
                                  int32_t block_height,
                                  const int16_t mv_param[2],
                                  lcu_t *lcu)
{
  int mv_frac_x = (mv_param[0] & 3);
  int mv_frac_y = (mv_param[1] & 3);

 #define FILTER_SIZE_Y 8 //Luma filter size

  // Fractional luma 1/4-pel
  kvz_extended_block src = {0, 0, 0, 0};

  // Fractional luma
  kvz_get_extended_block(xpos,
                         ypos,
                         mv_param[0] >> 2,
                         mv_param[1] >> 2,
                         state->tile->lcu_offset_x * LCU_WIDTH,
                         state->tile->lcu_offset_y * LCU_WIDTH,
                         ref->y,
                         ref->width,
                         ref->height,
                         FILTER_SIZE_Y,
                         block_width,
                         block_height,
                         &src);
  kvz_sample_quarterpel_luma(state->encoder_control,
                                     src.orig_topleft,
                                     src.stride,
                                     block_width,
                                     block_height,
                                     lcu->rec.y + (ypos%LCU_WIDTH)*LCU_WIDTH + (xpos%LCU_WIDTH),
                                     LCU_WIDTH,
                                     mv_frac_x,
                                     mv_frac_y,
                                     mv_param);

  if (src.malloc_used) free(src.buffer);
}

static void inter_recon_14bit_frac_luma(const encoder_state_t * const state,
                                        const kvz_picture * const ref,
                                        int32_t xpos,
                                        int32_t ypos,
                                        int32_t block_width,
                                        int32_t block_height,
                                        const int16_t mv_param[2],
                                        hi_prec_buf_t *hi_prec_out)
{
  int mv_frac_x = (mv_param[0] & 3);
  int mv_frac_y = (mv_param[1] & 3);

#define FILTER_SIZE_Y 8 //Luma filter size

  // Fractional luma 1/4-pel
  kvz_extended_block src = { 0, 0, 0, 0 };

  // Fractional luma
  kvz_get_extended_block(xpos,
                         ypos,
                         mv_param[0] >> 2,
                         mv_param[1] >> 2,
                         state->tile->lcu_offset_x * LCU_WIDTH,
                         state->tile->lcu_offset_y * LCU_WIDTH,
                         ref->y,
                         ref->width,
                         ref->height,
                         FILTER_SIZE_Y,
                         block_width,
                         block_height,
                         &src);
  kvz_sample_14bit_quarterpel_luma(state->encoder_control,
                                           src.orig_topleft,
                                           src.stride,
                                           block_width,
                                           block_height,
                                           hi_prec_out->y + (ypos%LCU_WIDTH)*LCU_WIDTH + (xpos%LCU_WIDTH),
                                           LCU_WIDTH,
                                           mv_frac_x,
                                           mv_frac_y,
                                           mv_param);

  if (src.malloc_used) free(src.buffer);
}

static void inter_recon_frac_chroma(const encoder_state_t * const state,
                                    const kvz_picture * const ref,
                                    int32_t xpos,
                                    int32_t ypos,
                                    int32_t block_width,
                                    int32_t block_height,
                                    const int16_t mv_param[2],
                                    lcu_t *lcu)
{
  int mv_frac_x = (mv_param[0] & 7);
  int mv_frac_y = (mv_param[1] & 7);

  // Translate to chroma
  xpos >>= 1;
  ypos >>= 1;
  block_width >>= 1;
  block_height >>= 1;

#define FILTER_SIZE_C 4 //Chroma filter size

  // Fractional chroma 1/8-pel
  kvz_extended_block src_u = { 0, 0, 0, 0 };
  kvz_extended_block src_v = { 0, 0, 0, 0 };

  //Fractional chroma U
  kvz_get_extended_block(xpos, ypos, (mv_param[0] >> 2) >> 1, (mv_param[1] >> 2) >> 1, state->tile->lcu_offset_x * LCU_WIDTH_C, state->tile->lcu_offset_y * LCU_WIDTH_C,
    ref->u, ref->width >> 1, ref->height >> 1, FILTER_SIZE_C, block_width, block_height, &src_u);
  kvz_sample_octpel_chroma(state->encoder_control, src_u.orig_topleft, src_u.stride, block_width,
    block_height, lcu->rec.u + (ypos % LCU_WIDTH_C)*LCU_WIDTH_C + (xpos % LCU_WIDTH_C), LCU_WIDTH_C, mv_frac_x, mv_frac_y, mv_param);

  //Fractional chroma V
  kvz_get_extended_block(xpos, ypos, (mv_param[0] >> 2) >> 1, (mv_param[1] >> 2) >> 1, state->tile->lcu_offset_x * LCU_WIDTH_C, state->tile->lcu_offset_y * LCU_WIDTH_C,
    ref->v, ref->width >> 1, ref->height >> 1, FILTER_SIZE_C, block_width, block_height, &src_v);
  kvz_sample_octpel_chroma(state->encoder_control, src_v.orig_topleft, src_v.stride, block_width,
    block_height, lcu->rec.v + (ypos  % LCU_WIDTH_C)*LCU_WIDTH_C + (xpos % LCU_WIDTH_C), LCU_WIDTH_C, mv_frac_x, mv_frac_y, mv_param);

  if (src_u.malloc_used) free(src_u.buffer);
  if (src_v.malloc_used) free(src_v.buffer);
}

static void inter_recon_14bit_frac_chroma(const encoder_state_t * const state,
                                          const kvz_picture * const ref,
                                          int32_t xpos,
                                          int32_t ypos,
                                          int32_t block_width,
                                          int32_t block_height,
                                          const int16_t mv_param[2],
                                          hi_prec_buf_t *hi_prec_out)
{
  int mv_frac_x = (mv_param[0] & 7);
  int mv_frac_y = (mv_param[1] & 7);

  // Translate to chroma
  xpos >>= 1;
  ypos >>= 1;
  block_width >>= 1;
  block_height >>= 1;

#define FILTER_SIZE_C 4 //Chroma filter size

  // Fractional chroma 1/8-pel
  kvz_extended_block src_u = { 0, 0, 0, 0 };
  kvz_extended_block src_v = { 0, 0, 0, 0 };

  //Fractional chroma U
  kvz_get_extended_block(xpos,
                         ypos,
                         (mv_param[0] >> 2) >> 1,
                         (mv_param[1] >> 2) >> 1,
                         state->tile->lcu_offset_x * LCU_WIDTH_C,
                         state->tile->lcu_offset_y * LCU_WIDTH_C,
                         ref->u,
                         ref->width >> 1,
                         ref->height >> 1,
                         FILTER_SIZE_C,
                         block_width,
                         block_height,
                         &src_u);
  kvz_sample_14bit_octpel_chroma(state->encoder_control,
                                         src_u.orig_topleft,
                                         src_u.stride,
                                         block_width,
                                         block_height,
                                         hi_prec_out->u + (ypos % LCU_WIDTH_C)*LCU_WIDTH_C + (xpos % LCU_WIDTH_C),
                                         LCU_WIDTH_C,
                                         mv_frac_x,
                                         mv_frac_y,
                                         mv_param);

  //Fractional chroma V
  kvz_get_extended_block(xpos,
                         ypos,
                         (mv_param[0] >> 2) >> 1,
                         (mv_param[1] >> 2) >> 1,
                         state->tile->lcu_offset_x * LCU_WIDTH_C,
                         state->tile->lcu_offset_y * LCU_WIDTH_C,
                         ref->v,
                         ref->width >> 1,
                         ref->height >> 1,
                         FILTER_SIZE_C,
                         block_width,
                         block_height,
                         &src_v);
  kvz_sample_14bit_octpel_chroma(state->encoder_control,
                                         src_v.orig_topleft,
                                         src_v.stride,
                                         block_width,
                                         block_height,
                                         hi_prec_out->v + (ypos  % LCU_WIDTH_C)*LCU_WIDTH_C + (xpos % LCU_WIDTH_C),
                                         LCU_WIDTH_C,
                                         mv_frac_x,
                                         mv_frac_y,
                                         mv_param);

  if (src_u.malloc_used) free(src_u.buffer);
  if (src_v.malloc_used) free(src_v.buffer);
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
 * \brief Reconstruct inter block
 *
 * \param state         encoder state
 * \param ref           picture to copy the data from
 * \param xpos          block x position
 * \param ypos          block y position
 * \param width         block width
 * \param height        block height
 * \param mv_param      motion vector
 * \param lcu           destination lcu
 * \param hi_prec_out   destination of high precision output (null if not needed)
*/
void kvz_inter_recon_lcu(const encoder_state_t * const state,
                         const kvz_picture * const ref,
                         int32_t xpos,
                         int32_t ypos,
                         int32_t width,
                         int32_t height,
                         const int16_t mv_param[2],
                         lcu_t *lcu,
                         hi_prec_buf_t *hi_prec_out)
{
  const vector2d_t tile_in_frame = {
    state->tile->lcu_offset_x * LCU_WIDTH,
    state->tile->lcu_offset_y * LCU_WIDTH
  };
  const vector2d_t pu_in_tile = { xpos, ypos };
  const vector2d_t pu_in_lcu = { xpos % LCU_WIDTH, ypos % LCU_WIDTH };

  const vector2d_t mv_in_pu = { mv_param[0] >> 2, mv_param[1] >> 2 };
  const vector2d_t mv_in_frame = {
    mv_in_pu.x + pu_in_tile.x + tile_in_frame.x,
    mv_in_pu.y + pu_in_tile.y + tile_in_frame.y
  };

  const bool mv_is_outside_frame = mv_in_frame.x < 0 ||
      mv_in_frame.y < 0 ||
      mv_in_frame.x + width > ref->width ||
      mv_in_frame.y + height > ref->height;

  // With 420, odd coordinates need interpolation.
  const int8_t fractional_chroma = (mv_in_pu.x & 1) || (mv_in_pu.y & 1);
  const int8_t fractional_luma = ((mv_param[0] & 3) || (mv_param[1] & 3));

  // Generate prediction for luma.
  if (fractional_luma) {
    // With a fractional MV, do interpolation.
    if (state->encoder_control->cfg->bipred && hi_prec_out) {
      inter_recon_14bit_frac_luma(state, ref,
                                  pu_in_tile.x, pu_in_tile.y,
                                  width, height,
                                  mv_param, hi_prec_out);
    } else {
      inter_recon_frac_luma(state, ref,
                            pu_in_tile.x, pu_in_tile.y,
                            width, height,
                            mv_param, lcu);
    }
  } else {
    // With an integer MV, copy pixels directly from the reference.
    const int lcu_pu_index = pu_in_lcu.y * LCU_WIDTH + pu_in_lcu.x;
    if (mv_is_outside_frame) {
      inter_cp_with_ext_border(ref->y, ref->width,
                               ref->width, ref->height,
                               &lcu->rec.y[lcu_pu_index], LCU_WIDTH,
                               width, height,
                               &mv_in_frame);
    } else {
      const int frame_mv_index = mv_in_frame.y * ref->width + mv_in_frame.x;
      kvz_pixels_blit(&ref->y[frame_mv_index],
                      &lcu->rec.y[lcu_pu_index],
                      width, height,
                      ref->width, LCU_WIDTH);
    }
  }

  if (state->encoder_control->chroma_format == KVZ_CSP_400) {
    return;
  }

  // Generate prediction for chroma.
  if (fractional_luma || fractional_chroma) {
    // With a fractional MV, do interpolation.
    if (state->encoder_control->cfg->bipred && hi_prec_out) {
      inter_recon_14bit_frac_chroma(state, ref,
                                    pu_in_tile.x, pu_in_tile.y,
                                    width, height,
                                    mv_param, hi_prec_out);
    } else {
      inter_recon_frac_chroma(state, ref,
                              pu_in_tile.x, pu_in_tile.y,
                              width, height,
                              mv_param, lcu);
    }
  } else {
    // With an integer MV, copy pixels directly from the reference.
    const int lcu_pu_index_c = pu_in_lcu.y / 2 * LCU_WIDTH_C + pu_in_lcu.x / 2;
    const vector2d_t mv_in_frame_c = { mv_in_frame.x / 2, mv_in_frame.y / 2 };

    if (mv_is_outside_frame) {
      inter_cp_with_ext_border(ref->u, ref->width / 2,
                               ref->width / 2, ref->height / 2,
                               &lcu->rec.u[lcu_pu_index_c], LCU_WIDTH_C,
                               width / 2, height / 2,
                               &mv_in_frame_c);
      inter_cp_with_ext_border(ref->v, ref->width / 2,
                               ref->width / 2, ref->height / 2,
                               &lcu->rec.v[lcu_pu_index_c], LCU_WIDTH_C,
                               width / 2, height / 2,
                               &mv_in_frame_c);
    } else {
      const int frame_mv_index = mv_in_frame_c.y * ref->width / 2 + mv_in_frame_c.x;

      kvz_pixels_blit(&ref->u[frame_mv_index],
                      &lcu->rec.u[lcu_pu_index_c],
                      width / 2, height / 2,
                      ref->width / 2, LCU_WIDTH_C);
      kvz_pixels_blit(&ref->v[frame_mv_index],
                      &lcu->rec.v[lcu_pu_index_c],
                      width / 2, height / 2,
                      ref->width / 2, LCU_WIDTH_C);
    }
  }
}

/**
 * \brief Reconstruct bi-pred inter block
 *
 * \param state     encoder state
 * \param ref1      reference picture to copy the data from
 * \param ref2      other reference picture to copy the data from
 * \param xpos      block x position
 * \param ypos      block y position
 * \param width     block width
 * \param height    block height
 * \param mv_param  motion vectors
 * \param lcu       destination lcu
 */
void kvz_inter_recon_lcu_bipred(const encoder_state_t * const state,
                                const kvz_picture * ref1,
                                const kvz_picture * ref2,
                                int32_t xpos,
                                int32_t ypos,
                                int32_t width,
                                int32_t height,
                                int16_t mv_param[2][2],
                                lcu_t* lcu)
{
  kvz_pixel temp_lcu_y[LCU_WIDTH*LCU_WIDTH];
  kvz_pixel temp_lcu_u[LCU_WIDTH_C*LCU_WIDTH_C];
  kvz_pixel temp_lcu_v[LCU_WIDTH_C*LCU_WIDTH_C];
  int temp_x, temp_y;
  int shift = 15 - KVZ_BIT_DEPTH;
  int offset = 1 << (shift - 1);

  const int hi_prec_luma_rec0 = mv_param[0][0] & 3 || mv_param[0][1] & 3;
  const int hi_prec_luma_rec1 = mv_param[1][0] & 3 || mv_param[1][1] & 3;

  const int hi_prec_chroma_rec0 = mv_param[0][0] & 7 || mv_param[0][1] & 7;
  const int hi_prec_chroma_rec1 = mv_param[1][0] & 7 || mv_param[1][1] & 7;

  hi_prec_buf_t* high_precision_rec0 = 0;
  hi_prec_buf_t* high_precision_rec1 = 0;
  if (hi_prec_chroma_rec0) high_precision_rec0 = kvz_hi_prec_buf_t_alloc(LCU_WIDTH*LCU_WIDTH);
  if (hi_prec_chroma_rec1) high_precision_rec1 = kvz_hi_prec_buf_t_alloc(LCU_WIDTH*LCU_WIDTH);
  //Reconstruct both predictors
  kvz_inter_recon_lcu(state, ref1, xpos, ypos, width, height, mv_param[0], lcu, high_precision_rec0);
  if (!hi_prec_luma_rec0){
    memcpy(temp_lcu_y, lcu->rec.y, sizeof(kvz_pixel) * 64 * 64);
  }
  if (!hi_prec_chroma_rec0){
    memcpy(temp_lcu_u, lcu->rec.u, sizeof(kvz_pixel) * 32 * 32);
    memcpy(temp_lcu_v, lcu->rec.v, sizeof(kvz_pixel) * 32 * 32);
  }
  kvz_inter_recon_lcu(state, ref2, xpos, ypos, width, height, mv_param[1], lcu, high_precision_rec1);

  // After reconstruction, merge the predictors by taking an average of each pixel
  for (temp_y = 0; temp_y < height; ++temp_y) {
    int y_in_lcu = ((ypos + temp_y) & ((LCU_WIDTH)-1));
    for (temp_x = 0; temp_x < width; ++temp_x) {
      int x_in_lcu = ((xpos + temp_x) & ((LCU_WIDTH)-1));
      int16_t sample0_y = (hi_prec_luma_rec0 ? high_precision_rec0->y[y_in_lcu * LCU_WIDTH + x_in_lcu] : (temp_lcu_y[y_in_lcu * LCU_WIDTH + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      int16_t sample1_y = (hi_prec_luma_rec1 ? high_precision_rec1->y[y_in_lcu * LCU_WIDTH + x_in_lcu] : (lcu->rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      lcu->rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] = (kvz_pixel)kvz_fast_clip_32bit_to_pixel((sample0_y + sample1_y + offset) >> shift);
    }

  }
  for (temp_y = 0; temp_y < height >> 1; ++temp_y) {
    int y_in_lcu = (((ypos >> 1) + temp_y) & (LCU_WIDTH_C - 1));
    for (temp_x = 0; temp_x < width >> 1; ++temp_x) {
      int x_in_lcu = (((xpos >> 1) + temp_x) & (LCU_WIDTH_C - 1));
      int16_t sample0_u = (hi_prec_chroma_rec0 ? high_precision_rec0->u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (temp_lcu_u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      int16_t sample1_u = (hi_prec_chroma_rec1 ? high_precision_rec1->u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (lcu->rec.u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      lcu->rec.u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] = (kvz_pixel)kvz_fast_clip_32bit_to_pixel((sample0_u + sample1_u + offset) >> shift);

      int16_t sample0_v = (hi_prec_chroma_rec0 ? high_precision_rec0->v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (temp_lcu_v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      int16_t sample1_v = (hi_prec_chroma_rec1 ? high_precision_rec1->v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (lcu->rec.v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      lcu->rec.v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] = (kvz_pixel)kvz_fast_clip_32bit_to_pixel((sample0_v + sample1_v + offset) >> shift);
    }
  }
  if (high_precision_rec0 != 0) kvz_hi_prec_buf_t_free(high_precision_rec0);
  if (high_precision_rec1 != 0) kvz_hi_prec_buf_t_free(high_precision_rec1);
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
* \param encoder encoder control struct to use
* \param x block x position in SCU
* \param y block y position in SCU
* \param width current block width
* \param height current block height
* \param H candidate H
* \param C1 candidate C1
*/
static void kvz_inter_get_temporal_merge_candidates(const encoder_state_t * const state,
                                             int32_t x,
                                             int32_t y,
                                             int32_t width,
                                             int32_t height,
                                             cu_info_t **C3,
                                             cu_info_t **H) {
  /*
  Predictor block locations
  _________
  |CurrentPU|
  | |C0|__  |
  |    |C3| |
  |_________|_
            |H|
  */

  *C3 = NULL;
  *H  = NULL;

  // Find temporal reference, closest POC
  if (state->frame->ref->used_size) {
    uint32_t poc_diff = UINT_MAX;
    int32_t closest_ref = 0;

    for (int temporal_cand = 0; temporal_cand < state->frame->ref->used_size; temporal_cand++) {
      int td = state->frame->poc - state->frame->ref->pocs[temporal_cand];

      td = td < 0 ? -td : td;
      if (td < poc_diff) {
        closest_ref = temporal_cand;
        poc_diff = td;
      }
    }

    cu_array_t *ref_cu_array = state->frame->ref->cu_arrays[closest_ref];
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
          *H = &ref_cu_array->data[H_offset];
        }
      }
    }
    uint32_t xColCtr = x + (width / 2);
    uint32_t yColCtr = y + (height / 2);

    // C3 must be inside the LCU, in the center position of current CU
    if (xColCtr < state->encoder_control->in.width && yColCtr < state->encoder_control->in.height) {
      uint32_t C3_offset = ((xColCtr >> 4) << 4) / SCU_WIDTH + ((((yColCtr >> 4) << 4) / SCU_WIDTH) * cu_per_width);
      if (ref_cu_array->data[C3_offset].type == CU_INTER) {
        *C3 = &ref_cu_array->data[C3_offset];
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
 * \param b0              Returns the b0 candidate.
 * \param b1              Returns the b1 candidate.
 * \param b2              Returns the b2 candidate.
 * \param a0              Returns the a0 candidate.
 * \param a1              Returns the a1 candidate.
 * \param lcu             current LCU
 */
static void get_spatial_merge_candidates(int32_t x,
                                         int32_t y,
                                         int32_t width,
                                         int32_t height,
                                         int32_t picture_width,
                                         int32_t picture_height,
                                         cu_info_t **b0,
                                         cu_info_t **b1,
                                         cu_info_t **b2,
                                         cu_info_t **a0,
                                         cu_info_t **a1,
                                         lcu_t *lcu)
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
    *a1 = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local + height - 1);
    // Do not check (*a1)->coded because the block above is always coded before
    // the current one and the flag is not set when searching an SMP block.
    if ((*a1)->type == CU_INTER) {
      inter_clear_cu_unused(*a1);
    } else {
      *a1 = NULL;
    }

    if (y_local + height < LCU_WIDTH && y + height < picture_height) {
      *a0 = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local + height);
      if ((*a0)->type == CU_INTER && is_a0_cand_coded(x, y, width, height)) {
        inter_clear_cu_unused(*a0);
      } else {
        *a0 = NULL;
      }
    }
  }

  // B0, B1 and B2 availability testing
  if (y != 0) {
    if (x + width < picture_width) {
      if (x_local + width < LCU_WIDTH) {
        *b0 = LCU_GET_CU_AT_PX(lcu, x_local + width, y_local - 1);
      } else if (y_local == 0) {
        // Special case, top-right CU
        *b0 = LCU_GET_TOP_RIGHT_CU(lcu);
      }
    }
    if ((*b0) && (*b0)->type == CU_INTER && is_b0_cand_coded(x, y, width, height)) {
      inter_clear_cu_unused(*b0);
    } else {
      *b0 = NULL;
    }

    *b1 = LCU_GET_CU_AT_PX(lcu, x_local + width - 1, y_local - 1);
    // Do not check (*b1)->coded because the block to the left is always coded
    // before the current one and the flag is not set when searching an SMP
    // block.
    if ((*b1)->type == CU_INTER) {
      inter_clear_cu_unused(*b1);
    } else {
      *b1 = NULL;
    }

    if (x != 0) {
      *b2 = LCU_GET_CU_AT_PX(lcu, x_local - 1, y_local - 1);
      // Do not check (*b2)->coded because the block above and to the left is
      // always coded before the current one.
      if ((*b2)->type == CU_INTER) {
        inter_clear_cu_unused(*b2);
      } else {
        *b2 = NULL;
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
 * \param b0              Returns the b0 candidate.
 * \param b1              Returns the b1 candidate.
 * \param b2              Returns the b2 candidate.
 * \param a0              Returns the a0 candidate.
 * \param a1              Returns the a1 candidate.
 */
static void get_spatial_merge_candidates_cua(const cu_array_t *cua,
                                             int32_t x,
                                             int32_t y,
                                             int32_t width,
                                             int32_t height,
                                             int32_t picture_width,
                                             int32_t picture_height,
                                             const cu_info_t **b0,
                                             const cu_info_t **b1,
                                             const cu_info_t **b2,
                                             const cu_info_t **a0,
                                             const cu_info_t **a1)
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
    *a1 = kvz_cu_array_at_const(cua, x - 1, y + height - 1);
    // The block above is always coded before the current one.
    if ((*a1)->type != CU_INTER) {
      *a1 = NULL;
    }

    if (y_local + height < LCU_WIDTH && y + height < picture_height) {
      *a0 = kvz_cu_array_at_const(cua, x - 1, y + height);
      if ((*a0)->type != CU_INTER || !is_a0_cand_coded(x, y, width, height)) {
        *a0 = NULL;
      }
    }
  }

  // B0, B1 and B2 availability testing
  if (y != 0) {
    if (x + width < picture_width && (x_local + width < LCU_WIDTH || y_local == 0)) {
      *b0 = kvz_cu_array_at_const(cua, x + width, y - 1);
      if ((*b0)->type != CU_INTER || !is_b0_cand_coded(x, y, width, height)) {
        *b0 = NULL;
      }
    }

    *b1 = kvz_cu_array_at_const(cua, x + width - 1, y - 1);
    // The block to the left is always coded before the current one.
    if ((*b1)->type != CU_INTER) {
      *b1 = NULL;
    }

    if (x != 0) {
      *b2 = kvz_cu_array_at_const(cua, x - 1, y - 1);
      // The block above and to the left is always coded before the current
      // one.
      if ((*b2)->type != CU_INTER) {
        *b2 = NULL;
      }
    }
  }
}

/**
 * \brief Pick two mv candidates from the spatial candidates.
 */
static void get_mv_cand_from_spatial(const encoder_state_t * const state,
                                     int32_t x,
                                     int32_t y,
                                     int32_t width,
                                     int32_t height,
                                     const cu_info_t *b0,
                                     const cu_info_t *b1,
                                     const cu_info_t *b2,
                                     const cu_info_t *a0,
                                     const cu_info_t *a1,
                                     const cu_info_t *c3,
                                     const cu_info_t *h,
                                     const cu_info_t *cur_cu,
                                     int8_t reflist,
                                     int16_t mv_cand[2][2])
{
  uint8_t candidates = 0;
  uint8_t b_candidates = 0;
  int8_t reflist2nd = !reflist;

 #define CALCULATE_SCALE(cu,tb,td) ((tb * ((0x4000 + (abs(td)>>1))/td) + 32) >> 6)
#define APPLY_MV_SCALING(cu, cand, list) {int td = state->frame->poc - state->frame->ref->pocs[(cu)->inter.mv_ref[list]];\
                                   int tb = state->frame->poc - state->frame->ref->pocs[cur_cu->inter.mv_ref[reflist]];\
                                   if (td != tb) { \
                                      int scale = CALCULATE_SCALE(cu,tb,td); \
                                       mv_cand[cand][0] = ((scale * (cu)->inter.mv[list][0] + 127 + (scale * (cu)->inter.mv[list][0] < 0)) >> 8 ); \
                                       mv_cand[cand][1] = ((scale * (cu)->inter.mv[list][1] + 127 + (scale * (cu)->inter.mv[list][1] < 0)) >> 8 ); }}

  // Left predictors
  if (a0 && (
    ((a0->inter.mv_dir & 1) && a0->inter.mv_ref[0] == cur_cu->inter.mv_ref[reflist]) ||
    ((a0->inter.mv_dir & 2) && a0->inter.mv_ref[1] == cur_cu->inter.mv_ref[reflist]))) {
    if (a0->inter.mv_dir & (1 << reflist) && a0->inter.mv_ref[reflist] == cur_cu->inter.mv_ref[reflist]) {
      mv_cand[candidates][0] = a0->inter.mv[reflist][0];
      mv_cand[candidates][1] = a0->inter.mv[reflist][1];
    } else {
      mv_cand[candidates][0] = a0->inter.mv[reflist2nd][0];
      mv_cand[candidates][1] = a0->inter.mv[reflist2nd][1];
    }
    candidates++;
  } else if (a1 && (
    ((a1->inter.mv_dir & 1) && a1->inter.mv_ref[0] == cur_cu->inter.mv_ref[reflist]) ||
    ((a1->inter.mv_dir & 2) && a1->inter.mv_ref[1] == cur_cu->inter.mv_ref[reflist]))) {
    if (a1->inter.mv_dir & (1 << reflist) && a1->inter.mv_ref[reflist] == cur_cu->inter.mv_ref[reflist]) {
      mv_cand[candidates][0] = a1->inter.mv[reflist][0];
      mv_cand[candidates][1] = a1->inter.mv[reflist][1];
    } else {
      mv_cand[candidates][0] = a1->inter.mv[reflist2nd][0];
      mv_cand[candidates][1] = a1->inter.mv[reflist2nd][1];
    }
    candidates++;
  }

  if(!candidates) {
      // Left predictors
    if (a0) {
      if (a0->inter.mv_dir & (1 << reflist)) {
        mv_cand[candidates][0] = a0->inter.mv[reflist][0];
        mv_cand[candidates][1] = a0->inter.mv[reflist][1];
        APPLY_MV_SCALING(a0, candidates, reflist);
      } else {
        mv_cand[candidates][0] = a0->inter.mv[reflist2nd][0];
        mv_cand[candidates][1] = a0->inter.mv[reflist2nd][1];
        APPLY_MV_SCALING(a0, candidates, reflist2nd);
      }
      candidates++;
    } else if (a1) {
      if (a1->inter.mv_dir & (1 << reflist)) {
        mv_cand[candidates][0] = a1->inter.mv[reflist][0];
        mv_cand[candidates][1] = a1->inter.mv[reflist][1];
        APPLY_MV_SCALING(a1, candidates, reflist);
      } else {
        mv_cand[candidates][0] = a1->inter.mv[reflist2nd][0];
        mv_cand[candidates][1] = a1->inter.mv[reflist2nd][1];
        APPLY_MV_SCALING(a1, candidates, reflist2nd);
      }
      candidates++;
    }
  }

  // Top predictors
  if (b0 && (
    ((b0->inter.mv_dir & 1) && b0->inter.mv_ref[0] == cur_cu->inter.mv_ref[reflist]) ||
    ((b0->inter.mv_dir & 2) && b0->inter.mv_ref[1] == cur_cu->inter.mv_ref[reflist]))) {
    if (b0->inter.mv_dir & (1 << reflist) && b0->inter.mv_ref[reflist] == cur_cu->inter.mv_ref[reflist]) {
      mv_cand[candidates][0] = b0->inter.mv[reflist][0];
      mv_cand[candidates][1] = b0->inter.mv[reflist][1];
    } else {
      mv_cand[candidates][0] = b0->inter.mv[reflist2nd][0];
      mv_cand[candidates][1] = b0->inter.mv[reflist2nd][1];
    }
    b_candidates++;
  } else if (b1 && (
    ((b1->inter.mv_dir & 1) && b1->inter.mv_ref[0] == cur_cu->inter.mv_ref[reflist]) ||
    ((b1->inter.mv_dir & 2) && b1->inter.mv_ref[1] == cur_cu->inter.mv_ref[reflist]))) {
    if (b1->inter.mv_dir & (1 << reflist) && b1->inter.mv_ref[reflist] == cur_cu->inter.mv_ref[reflist]) {
      mv_cand[candidates][0] = b1->inter.mv[reflist][0];
      mv_cand[candidates][1] = b1->inter.mv[reflist][1];
    } else {
      mv_cand[candidates][0] = b1->inter.mv[reflist2nd][0];
      mv_cand[candidates][1] = b1->inter.mv[reflist2nd][1];
    }
    b_candidates++;
  } else if (b2 && (
    ((b2->inter.mv_dir & 1) && b2->inter.mv_ref[0] == cur_cu->inter.mv_ref[reflist]) ||
    ((b2->inter.mv_dir & 2) && b2->inter.mv_ref[1] == cur_cu->inter.mv_ref[reflist]))) {
    if (b2->inter.mv_dir & (1 << reflist) && b2->inter.mv_ref[reflist] == cur_cu->inter.mv_ref[reflist]) {
      mv_cand[candidates][0] = b2->inter.mv[reflist][0];
      mv_cand[candidates][1] = b2->inter.mv[reflist][1];
    } else {
      mv_cand[candidates][0] = b2->inter.mv[reflist2nd][0];
      mv_cand[candidates][1] = b2->inter.mv[reflist2nd][1];
    }
    b_candidates++;
  }
  candidates += b_candidates;

  // When a1 or a0 is available, we dont check for secondary B candidates
  if (a1 || a0) {
    b_candidates = 1;
  } else if(candidates != 2) {
    b_candidates = 0;
  }

  if(!b_candidates) {
    // Top predictors
    if (b0) {
      if (b0->inter.mv_dir & (1 << reflist)) {
        mv_cand[candidates][0] = b0->inter.mv[reflist][0];
        mv_cand[candidates][1] = b0->inter.mv[reflist][1];
        APPLY_MV_SCALING(b0, candidates, reflist);
      } else {
        mv_cand[candidates][0] = b0->inter.mv[reflist2nd][0];
        mv_cand[candidates][1] = b0->inter.mv[reflist2nd][1];
        APPLY_MV_SCALING(b0, candidates, reflist2nd);
      }
      candidates++;
    } else if (b1) {
      if (b1->inter.mv_dir & (1 << reflist)) {
        mv_cand[candidates][0] = b1->inter.mv[reflist][0];
        mv_cand[candidates][1] = b1->inter.mv[reflist][1];
        APPLY_MV_SCALING(b1, candidates, reflist);
      } else {
        mv_cand[candidates][0] = b1->inter.mv[reflist2nd][0];
        mv_cand[candidates][1] = b1->inter.mv[reflist2nd][1];
        APPLY_MV_SCALING(b1, candidates, reflist2nd);
      }
      candidates++;
    } else if (b2) {
      if (b2->inter.mv_dir & (1 << reflist)) {
        mv_cand[candidates][0] = b2->inter.mv[reflist][0];
        mv_cand[candidates][1] = b2->inter.mv[reflist][1];
        APPLY_MV_SCALING(b2, candidates, reflist);
      } else {
        mv_cand[candidates][0] = b2->inter.mv[reflist2nd][0];
        mv_cand[candidates][1] = b2->inter.mv[reflist2nd][1];
        APPLY_MV_SCALING(b2, candidates, reflist2nd);
      }
      candidates++;
    }
  }

  // Remove identical candidate
  if(candidates == 2 && mv_cand[0][0] == mv_cand[1][0] && mv_cand[0][1] == mv_cand[1][1]) {
    candidates = 1;
  }

  if (state->encoder_control->cfg->tmvp_enable) {
    /*
    Predictor block locations
    _________
    |CurrentPU|
    | |C0|__  |
    |    |C3| |
    |_________|_
    |H|
    */

    // Find temporal reference, closest POC
    if (state->frame->poc > 1 && state->frame->ref->used_size && candidates < AMVP_MAX_NUM_CANDS) {
      uint32_t poc_diff = UINT_MAX;

      for (int temporal_cand = 0; temporal_cand < state->frame->ref->used_size; temporal_cand++) {
        int td = state->frame->poc - state->frame->ref->pocs[temporal_cand];
        td = td < 0 ? -td : td;
        if (td < poc_diff) {
          poc_diff = td;
        }
      }

      const cu_info_t *selected_CU = (h != NULL) ? h : (c3 != NULL) ? c3 : NULL;

      if (selected_CU) {
        int td = selected_CU->inter.mv_ref[reflist] + 1;
        int tb = cur_cu->inter.mv_ref[reflist] + 1;

        int scale = CALCULATE_SCALE(NULL, tb, td);
        mv_cand[candidates][0] = ((scale * selected_CU->inter.mv[0][0] + 127 + (scale * selected_CU->inter.mv[0][0] < 0)) >> 8);
        mv_cand[candidates][1] = ((scale * selected_CU->inter.mv[0][1] + 127 + (scale * selected_CU->inter.mv[0][1] < 0)) >> 8);

        candidates++;
      }
#undef CALCULATE_SCALE
    }
  }

  // Fill with (0,0)
  while (candidates < AMVP_MAX_NUM_CANDS) {
    mv_cand[candidates][0] = 0;
    mv_cand[candidates][1] = 0;
    candidates++;
  }
#undef CALCULATE_SCALE
#undef APPLY_MV_SCALING
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
                           cu_info_t* cur_cu,
                           lcu_t *lcu,
                           int8_t reflist)
{
  cu_info_t *b0, *b1, *b2, *a0, *a1, *c3, *h;
  b0 = b1 = b2 = a0 = a1 = c3 = h = NULL;
  get_spatial_merge_candidates(x, y, width, height,
                               state->tile->frame->width, state->tile->frame->height,
                               &b0, &b1, &b2, &a0, &a1, lcu);
  kvz_inter_get_temporal_merge_candidates(state, x, y, width, height, &c3, &h);
  get_mv_cand_from_spatial(state, x, y, width, height, b0, b1, b2, a0, a1, c3, h, cur_cu, reflist, mv_cand);
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
  const cu_info_t *b0, *b1, *b2, *a0, *a1;
  cu_info_t *c3, *h;
  b0 = b1 = b2 = a0 = a1 = c3 = h = NULL;
  
  const cu_array_t *cua = state->tile->frame->cu_array;
  get_spatial_merge_candidates_cua(cua,
                                   x, y, width, height,
                                   state->tile->frame->width, state->tile->frame->height,
                                   &b0, &b1, &b2, &a0, &a1);
  kvz_inter_get_temporal_merge_candidates(state, x, y, width, height, &c3, &h);
  get_mv_cand_from_spatial(state, x, y, width, height, b0, b1, b2, a0, a1, c3, h, cur_cu, reflist, mv_cand);
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
  int8_t duplicate = 0;

  cu_info_t *b0, *b1, *b2, *a0, *a1;
  int8_t zero_idx = 0;
  b0 = b1 = b2 = a0 = a1 = NULL;
  get_spatial_merge_candidates(x, y, width, height,
                               state->tile->frame->width, state->tile->frame->height,
                               &b0, &b1, &b2, &a0, &a1, lcu);

  if (!use_a1) a1 = NULL;
  if (!use_b1) b1 = NULL;

#define CHECK_DUPLICATE(CU1,CU2) {duplicate = 0; if ((CU2) && \
                                                     (CU1)->inter.mv_dir == (CU2)->inter.mv_dir && \
                                                    (!(((CU1)->inter.mv_dir & 1) && ((CU2)->inter.mv_dir & 1)) || \
                                                      ((CU1)->inter.mv[0][0] == (CU2)->inter.mv[0][0] && \
                                                       (CU1)->inter.mv[0][1] ==  (CU2)->inter.mv[0][1] && \
                                                       (CU1)->inter.mv_ref[0] == (CU2)->inter.mv_ref[0]) ) && \
                                                    (!(((CU1)->inter.mv_dir & 2) && ((CU2)->inter.mv_dir & 2) )  || \
                                                      ((CU1)->inter.mv[1][0] == (CU2)->inter.mv[1][0] && \
                                                       (CU1)->inter.mv[1][1] == (CU2)->inter.mv[1][1] && \
                                                       (CU1)->inter.mv_ref[1] == (CU2)->inter.mv_ref[1]) ) \
                                                      ) duplicate = 1; }

  if (a1) {
    mv_cand[candidates].mv[0][0] = a1->inter.mv[0][0];
    mv_cand[candidates].mv[0][1] = a1->inter.mv[0][1];
    mv_cand[candidates].mv[1][0] = a1->inter.mv[1][0];
    mv_cand[candidates].mv[1][1] = a1->inter.mv[1][1];
    mv_cand[candidates].ref[0] = a1->inter.mv_ref[0];
    mv_cand[candidates].ref[1] = a1->inter.mv_ref[1];
    mv_cand[candidates].dir = a1->inter.mv_dir;
    candidates++;
  }

  if (b1) {
    if(candidates) CHECK_DUPLICATE(b1, a1);
    if(!duplicate) {
      mv_cand[candidates].mv[0][0] = b1->inter.mv[0][0];
      mv_cand[candidates].mv[0][1] = b1->inter.mv[0][1];
      mv_cand[candidates].mv[1][0] = b1->inter.mv[1][0];
      mv_cand[candidates].mv[1][1] = b1->inter.mv[1][1];
      mv_cand[candidates].ref[0] = b1->inter.mv_ref[0];
      mv_cand[candidates].ref[1] = b1->inter.mv_ref[1];
      mv_cand[candidates].dir = b1->inter.mv_dir;
      candidates++;
    }
  }

  if (b0) {
    if(candidates) CHECK_DUPLICATE(b0,b1);
    if(!duplicate) {
      mv_cand[candidates].mv[0][0] = b0->inter.mv[0][0];
      mv_cand[candidates].mv[0][1] = b0->inter.mv[0][1];
      mv_cand[candidates].mv[1][0] = b0->inter.mv[1][0];
      mv_cand[candidates].mv[1][1] = b0->inter.mv[1][1];
      mv_cand[candidates].ref[0] = b0->inter.mv_ref[0];
      mv_cand[candidates].ref[1] = b0->inter.mv_ref[1];
      mv_cand[candidates].dir = b0->inter.mv_dir;
      candidates++;
    }
  }

  if (a0) {
    if(candidates) CHECK_DUPLICATE(a0,a1);
    if(!duplicate) {
      mv_cand[candidates].mv[0][0] = a0->inter.mv[0][0];
      mv_cand[candidates].mv[0][1] = a0->inter.mv[0][1];
      mv_cand[candidates].mv[1][0] = a0->inter.mv[1][0];
      mv_cand[candidates].mv[1][1] = a0->inter.mv[1][1];
      mv_cand[candidates].ref[0] = a0->inter.mv_ref[0];
      mv_cand[candidates].ref[1] = a0->inter.mv_ref[1];
      mv_cand[candidates].dir = a0->inter.mv_dir;
      candidates++;
    }
  }

  if (candidates != 4) {
    if (b2) {
      CHECK_DUPLICATE(b2,a1);
      if(!duplicate) {
        CHECK_DUPLICATE(b2,b1);
        if(!duplicate) {
          mv_cand[candidates].mv[0][0] = b2->inter.mv[0][0];
          mv_cand[candidates].mv[0][1] = b2->inter.mv[0][1];
          mv_cand[candidates].mv[1][0] = b2->inter.mv[1][0];
          mv_cand[candidates].mv[1][1] = b2->inter.mv[1][1];
          mv_cand[candidates].ref[0] = b2->inter.mv_ref[0];
          mv_cand[candidates].ref[1] = b2->inter.mv_ref[1];
          mv_cand[candidates].dir = b2->inter.mv_dir;
          candidates++;
        }
      }
    }
  }
  
  if (state->encoder_control->cfg->tmvp_enable) {
#define CALCULATE_SCALE(cu,tb,td) ((tb * ((0x4000 + (abs(td)>>1))/td) + 32) >> 6)

    if (candidates < MRG_MAX_NUM_CANDS && state->frame->ref->used_size) {

      cu_info_t *c3 = NULL;
      cu_info_t *h = NULL;

      kvz_inter_get_temporal_merge_candidates(state, x, y, width, height, &c3, &h);

      const cu_info_t *selected_CU = (h != NULL) ? h : (c3 != NULL) ? c3 : NULL;

      if (selected_CU) {
        int td = selected_CU->inter.mv_ref[0] + 1;
        int tb = 1;

        int scale = CALCULATE_SCALE(NULL, tb, td);
        mv_cand[candidates].mv[0][0] = ((scale * selected_CU->inter.mv[0][0] + 127 + (scale * selected_CU->inter.mv[0][0] < 0)) >> 8);
        mv_cand[candidates].mv[0][1] = ((scale * selected_CU->inter.mv[0][1] + 127 + (scale * selected_CU->inter.mv[0][1] < 0)) >> 8);

        /*
        ToDo: temporal prediction in B-pictures
        mv_cand[candidates].mv[1][0] = selected_CU->inter.mv[1][0];
        mv_cand[candidates].mv[1][1] = selected_CU->inter.mv[1][1];
        */
        mv_cand[candidates].dir = selected_CU->inter.mv_dir;
        mv_cand[candidates].ref[0] = 0;
        candidates++;
      }
    }
#undef CALCULATE_SCALE
  }

  if (candidates < MRG_MAX_NUM_CANDS && state->frame->slicetype == KVZ_SLICE_B) {
    #define NUM_PRIORITY_LIST 12;
    static const uint8_t priorityList0[] = { 0, 1, 0, 2, 1, 2, 0, 3, 1, 3, 2, 3 };
    static const uint8_t priorityList1[] = { 1, 0, 2, 0, 2, 1, 3, 0, 3, 1, 3, 2 };
    uint8_t cutoff = candidates;
    for (int32_t idx = 0; idx<cutoff*(cutoff - 1) && candidates != MRG_MAX_NUM_CANDS; idx++) {
      uint8_t i = priorityList0[idx];
      uint8_t j = priorityList1[idx];
      if (i >= candidates || j >= candidates) break;

      // Find one L0 and L1 candidate according to the priority list
      if ((mv_cand[i].dir & 0x1) && (mv_cand[j].dir & 0x2)) {
        mv_cand[candidates].dir = 3;

        // get Mv from cand[i] and cand[j]
        mv_cand[candidates].mv[0][0] = mv_cand[i].mv[0][0];
        mv_cand[candidates].mv[0][1] = mv_cand[i].mv[0][1];
        mv_cand[candidates].mv[1][0] = mv_cand[j].mv[1][0];
        mv_cand[candidates].mv[1][1] = mv_cand[j].mv[1][1];
        mv_cand[candidates].ref[0]   = mv_cand[i].ref[0];
        mv_cand[candidates].ref[1]   = mv_cand[j].ref[1];

        if (mv_cand[i].ref[0] == mv_cand[j].ref[1] &&
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

  if (candidates < MRG_MAX_NUM_CANDS && state->frame->slicetype == KVZ_SLICE_B) {
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
  while (candidates != MRG_MAX_NUM_CANDS) {
    mv_cand[candidates].mv[0][0] = 0;
    mv_cand[candidates].mv[0][1] = 0;
    mv_cand[candidates].ref[0] = (zero_idx>=num_ref-1)?0:zero_idx;
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
