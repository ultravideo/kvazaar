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

/*
 * \file
 */

#include "inter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "filter.h"
#include "strategies/strategies-ipol.h"
#include "strategies/generic/ipol-generic.h"
#include "strategies/generic/picture-generic.h"

/**
 * \brief Set block info to the CU structure
 * \param pic picture to use
 * \param x_cu x CU position (smallest CU)
 * \param y_cu y CU position (smallest CU)
 * \param depth current CU depth
 * \param cur_cu CU to take the settings from
 * \returns Void
*/
void inter_set_block(videoframe_t* frame, uint32_t x_cu, uint32_t y_cu, uint8_t depth, cu_info_t* cur_cu)
{
  uint32_t x, y;
  // Width in smallest CU
  int block_scu_width = (LCU_WIDTH>>depth)/(LCU_WIDTH>>MAX_DEPTH);
  int tr_depth = (depth == 0 ? 1 : depth);
  // Loop through all the block in the area of cur_cu
  for (y = y_cu; y < y_cu + block_scu_width; y++) {
    for (x = x_cu; x < x_cu + block_scu_width; x++) {
      cu_info_t * const cu = videoframe_get_cu(frame, x, y);
      // Set all SCU's to this blocks values at the bottom most depth.
      cu->depth = depth;
      cu->type  = CU_INTER;
      cu->part_size = SIZE_2Nx2N;
      memcpy(&cu->inter, &cur_cu->inter, sizeof(cur_cu->inter));
      
      cu->tr_depth = tr_depth;
    }
  }
}

void inter_recon_frac_luma(const encoder_state_t * const state, const kvz_picture * const ref, int32_t xpos, int32_t ypos, int32_t block_width, const int16_t mv_param[2], lcu_t *lcu)
{
  int mv_frac_x = (mv_param[0] & 3);
  int mv_frac_y = (mv_param[1] & 3);

  int y, x;

  #define FILTER_SIZE_Y 8 //Luma filter size

  // Fractional luma 1/4-pel
  kvz_pixel qpel_src_y[(LCU_WIDTH + FILTER_SIZE_Y) * (LCU_WIDTH + FILTER_SIZE_Y)];
  kvz_pixel* qpel_src_off_y = &qpel_src_y[(block_width + FILTER_SIZE_Y)*(FILTER_SIZE_Y >> 1) + (FILTER_SIZE_Y >> 1)];

  // Fractional luma
  extend_borders(xpos, ypos, mv_param[0] >> 2, mv_param[1] >> 2, state->tile->lcu_offset_x * LCU_WIDTH, state->tile->lcu_offset_y * LCU_WIDTH,
    ref->y, ref->width, ref->height, FILTER_SIZE_Y, block_width, block_width, qpel_src_y);
  sample_quarterpel_luma_generic(state->encoder_control, qpel_src_off_y, block_width + FILTER_SIZE_Y, block_width,
    block_width, lcu->rec.y + (ypos%LCU_WIDTH)*LCU_WIDTH + (xpos%LCU_WIDTH), LCU_WIDTH, mv_frac_x, mv_frac_y, mv_param);
}

void inter_recon_14bit_frac_luma(const encoder_state_t * const state, const kvz_picture * const ref, int32_t xpos, int32_t ypos, int32_t block_width, const int16_t mv_param[2], hi_prec_buf_t *hi_prec_out)
{
  int mv_frac_x = (mv_param[0] & 3);
  int mv_frac_y = (mv_param[1] & 3);

  int y, x;

#define FILTER_SIZE_Y 8 //Luma filter size

  // Fractional luma 1/4-pel
  kvz_pixel qpel_src_y[(LCU_WIDTH + FILTER_SIZE_Y) * (LCU_WIDTH + FILTER_SIZE_Y)];
  kvz_pixel* qpel_src_off_y = &qpel_src_y[(block_width + FILTER_SIZE_Y)*(FILTER_SIZE_Y >> 1) + (FILTER_SIZE_Y >> 1)];

  // Fractional luma
  extend_borders(xpos, ypos, mv_param[0] >> 2, mv_param[1] >> 2, state->tile->lcu_offset_x * LCU_WIDTH, state->tile->lcu_offset_y * LCU_WIDTH,
    ref->y, ref->width, ref->height, FILTER_SIZE_Y, block_width, block_width, qpel_src_y);
  sample_14bit_quarterpel_luma_generic(state->encoder_control, qpel_src_off_y, block_width + FILTER_SIZE_Y, block_width,
    block_width, hi_prec_out->y + (ypos%LCU_WIDTH)*LCU_WIDTH + (xpos%LCU_WIDTH), LCU_WIDTH, mv_frac_x, mv_frac_y, mv_param);
}

void inter_recon_frac_chroma(const encoder_state_t * const state, const kvz_picture * const ref, int32_t xpos, int32_t ypos, int32_t block_width, const int16_t mv_param[2], lcu_t *lcu)
{
  int mv_frac_x = (mv_param[0] & 7);
  int mv_frac_y = (mv_param[1] & 7);

  // Translate to chroma
  xpos >>= 1;
  ypos >>= 1;
  block_width >>= 1;

#define FILTER_SIZE_C 4 //Chroma filter size

  // Fractional chroma 1/8-pel
  kvz_pixel octpel_src[((LCU_WIDTH_C) + FILTER_SIZE_C) * ((LCU_WIDTH_C) + FILTER_SIZE_C)];
  kvz_pixel* octpel_src_off = &octpel_src[(block_width + FILTER_SIZE_C)*(FILTER_SIZE_C >> 1) + (FILTER_SIZE_C >> 1)];

  //Fractional chroma U
  extend_borders(xpos, ypos, (mv_param[0] >> 2) >> 1, (mv_param[1] >> 2) >> 1, state->tile->lcu_offset_x * LCU_WIDTH_C, state->tile->lcu_offset_y * LCU_WIDTH_C,
    ref->u, ref->width >> 1, ref->height >> 1, FILTER_SIZE_C, block_width, block_width, octpel_src);
  sample_octpel_chroma_generic(state->encoder_control, octpel_src_off, block_width + FILTER_SIZE_C, block_width,
    block_width, lcu->rec.u + (ypos % LCU_WIDTH_C)*LCU_WIDTH_C + (xpos % LCU_WIDTH_C), LCU_WIDTH_C, mv_frac_x, mv_frac_y, mv_param);

  //Fractional chroma V
  extend_borders(xpos, ypos, (mv_param[0] >> 2) >> 1, (mv_param[1] >> 2) >> 1, state->tile->lcu_offset_x * LCU_WIDTH_C, state->tile->lcu_offset_y * LCU_WIDTH_C,
    ref->v, ref->width >> 1, ref->height >> 1, FILTER_SIZE_C, block_width, block_width, octpel_src);
  sample_octpel_chroma_generic(state->encoder_control, octpel_src_off, block_width + FILTER_SIZE_C, block_width,
    block_width, lcu->rec.v + (ypos  % LCU_WIDTH_C)*LCU_WIDTH_C + (xpos % LCU_WIDTH_C), LCU_WIDTH_C, mv_frac_x, mv_frac_y, mv_param);
}

void inter_recon_14bit_frac_chroma(const encoder_state_t * const state, const kvz_picture * const ref, int32_t xpos, int32_t ypos, int32_t block_width, const int16_t mv_param[2], hi_prec_buf_t *hi_prec_out)
{
  int mv_frac_x = (mv_param[0] & 7);
  int mv_frac_y = (mv_param[1] & 7);

  // Translate to chroma
  xpos >>= 1;
  ypos >>= 1;
  block_width >>= 1;

#define FILTER_SIZE_C 4 //Chroma filter size

  // Fractional chroma 1/8-pel
  kvz_pixel octpel_src[((LCU_WIDTH_C)+FILTER_SIZE_C) * ((LCU_WIDTH_C)+FILTER_SIZE_C)];
  kvz_pixel* octpel_src_off = &octpel_src[(block_width + FILTER_SIZE_C)*(FILTER_SIZE_C >> 1) + (FILTER_SIZE_C >> 1)];

  //Fractional chroma U
  extend_borders(xpos, ypos, (mv_param[0] >> 2) >> 1, (mv_param[1] >> 2) >> 1, state->tile->lcu_offset_x * LCU_WIDTH_C, state->tile->lcu_offset_y * LCU_WIDTH_C,
    ref->u, ref->width >> 1, ref->height >> 1, FILTER_SIZE_C, block_width, block_width, octpel_src);
  sample_14bit_octpel_chroma_generic(state->encoder_control, octpel_src_off, block_width + FILTER_SIZE_C, block_width,
    block_width, hi_prec_out->u + (ypos % LCU_WIDTH_C)*LCU_WIDTH_C + (xpos % LCU_WIDTH_C), LCU_WIDTH_C, mv_frac_x, mv_frac_y, mv_param);

  //Fractional chroma V
  extend_borders(xpos, ypos, (mv_param[0] >> 2) >> 1, (mv_param[1] >> 2) >> 1, state->tile->lcu_offset_x * LCU_WIDTH_C, state->tile->lcu_offset_y * LCU_WIDTH_C,
    ref->v, ref->width >> 1, ref->height >> 1, FILTER_SIZE_C, block_width, block_width, octpel_src);
  sample_14bit_octpel_chroma_generic(state->encoder_control, octpel_src_off, block_width + FILTER_SIZE_C, block_width,
    block_width, hi_prec_out->v + (ypos  % LCU_WIDTH_C)*LCU_WIDTH_C + (xpos % LCU_WIDTH_C), LCU_WIDTH_C, mv_frac_x, mv_frac_y, mv_param);
}

/**
 * \brief Reconstruct inter block
 * \param ref picture to copy the data from
 * \param xpos block x position
 * \param ypos block y position
 * \param width block width
 * \param mv[2] motion vector
 * \param lcu destination lcu
 * \param hi_prec destination of high precision output (null if not needed)
 * \returns Void
*/
void inter_recon_lcu(const encoder_state_t * const state, const kvz_picture * const ref, int32_t xpos, int32_t ypos,int32_t width, const int16_t mv_param[2], lcu_t *lcu, hi_prec_buf_t *hi_prec_out)
{
  int x,y,coord_x,coord_y;
  int16_t mv[2] = { mv_param[0], mv_param[1] };

  int32_t dst_width_c = LCU_WIDTH>>1; //!< Destination picture width in chroma pixels
  int32_t ref_width_c = ref->width>>1; //!< Reference picture width in chroma pixels

  // negative overflow flag
  int8_t overflow_neg_x = (state->tile->lcu_offset_x * LCU_WIDTH + xpos + (mv[0]>>2) < 0)?1:0;
  int8_t overflow_neg_y = (state->tile->lcu_offset_y * LCU_WIDTH + ypos + (mv[1]>>2) < 0)?1:0;

  // positive overflow flag
  int8_t overflow_pos_x = (state->tile->lcu_offset_x * LCU_WIDTH + xpos + (mv[0]>>2) + width > ref->width )?1:0;
  int8_t overflow_pos_y = (state->tile->lcu_offset_y * LCU_WIDTH + ypos + (mv[1]>>2) + width > ref->height)?1:0;

  int8_t chroma_halfpel = ((mv[0]>>2)&1) || ((mv[1]>>2)&1); //!< (luma integer mv) lsb is set -> chroma is half-pel
  // Luma quarter-pel
  int8_t fractional_mv = (mv[0]&1) || (mv[1]&1) || (mv[0]&2) || (mv[1]&2); // either of 2 lowest bits of mv set -> mv is fractional

  if(fractional_mv) {
    if (state->encoder_control->cfg->bipred && hi_prec_out){
      inter_recon_14bit_frac_luma(state, ref, xpos, ypos, width, mv_param, hi_prec_out);
      inter_recon_14bit_frac_chroma(state, ref, xpos, ypos, width, mv_param, hi_prec_out);
    } else {
      inter_recon_frac_luma(state, ref, xpos, ypos, width, mv_param, lcu);
      inter_recon_frac_chroma(state, ref, xpos, ypos, width, mv_param, lcu);
    }
  }

  mv[0] >>= 2;
  mv[1] >>= 2;

  // Chroma half-pel
  // get half-pel interpolated block and push it to output
  if(!fractional_mv) {
    if(chroma_halfpel) {
      if (state->encoder_control->cfg->bipred && hi_prec_out){
        inter_recon_14bit_frac_chroma(state, ref, xpos, ypos, width, mv_param, hi_prec_out);
      } else {
        inter_recon_frac_chroma(state, ref, xpos, ypos, width, mv_param, lcu);
      }
    }

    // With overflow present, more checking
    if (overflow_neg_x || overflow_neg_y || overflow_pos_x || overflow_pos_y) {
      // Copy Luma with boundary checking
      for (y = ypos; y < ypos + width; y++) {
        for (x = xpos; x < xpos + width; x++) {
          int x_in_lcu = (x & ((LCU_WIDTH)-1));
          int y_in_lcu = (y & ((LCU_WIDTH)-1));

          coord_x = (x + state->tile->lcu_offset_x * LCU_WIDTH) + mv[0];
          coord_y = (y + state->tile->lcu_offset_y * LCU_WIDTH) + mv[1];
          overflow_neg_x = (coord_x < 0)?1:0;
          overflow_neg_y = (coord_y < 0)?1:0;

          overflow_pos_x = (coord_x >= ref->width )?1:0;
          overflow_pos_y = (coord_y >= ref->height)?1:0;

          // On x-overflow set coord_x accordingly
          if (overflow_neg_x) {
            coord_x = 0;
          } else if (overflow_pos_x) {
            coord_x = ref->width - 1;
          }

          // On y-overflow set coord_y accordingly
          if (overflow_neg_y) {
            coord_y = 0;
          } else if (overflow_pos_y) {
            coord_y = ref->height - 1;
          }

          // set destination to (corrected) pixel value from the reference
          lcu->rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] = ref->y[coord_y*ref->width + coord_x];
        }
      }

      if(!chroma_halfpel) {
        // Copy Chroma with boundary checking
        for (y = ypos>>1; y < (ypos + width)>>1; y++) {
          for (x = xpos>>1; x < (xpos + width)>>1; x++) {
            int x_in_lcu = (x & ((LCU_WIDTH>>1)-1));
            int y_in_lcu = (y & ((LCU_WIDTH>>1)-1));

            coord_x = (x + state->tile->lcu_offset_x * (LCU_WIDTH >> 1)) + (mv[0]>>1);
            coord_y = (y + state->tile->lcu_offset_y * (LCU_WIDTH >> 1)) + (mv[1]>>1);

            overflow_neg_x = (coord_x < 0)?1:0;
            overflow_neg_y = (y + (mv[1]>>1) < 0)?1:0;

            overflow_pos_x = (coord_x >= ref->width>>1 )?1:0;
            overflow_pos_y = (coord_y >= ref->height>>1)?1:0;

            // On x-overflow set coord_x accordingly
            if(overflow_neg_x) {
              coord_x = 0;
            } else if(overflow_pos_x) {
              coord_x = (ref->width>>1) - 1;
            }

            // On y-overflow set coord_y accordingly
            if(overflow_neg_y) {
              coord_y = 0;
            } else if(overflow_pos_y) {
              coord_y = (ref->height>>1) - 1;
            }

            // set destinations to (corrected) pixel value from the reference
            lcu->rec.u[y_in_lcu*dst_width_c + x_in_lcu] = ref->u[coord_y * ref_width_c + coord_x];
            lcu->rec.v[y_in_lcu*dst_width_c + x_in_lcu] = ref->v[coord_y * ref_width_c + coord_x];
          }
        }
      }
    } else { //If no overflow, we can copy without checking boundaries
      // Copy Luma
      for (y = ypos; y < ypos + width; y++) {
        int y_in_lcu = (y & ((LCU_WIDTH)-1));
        coord_y = ((y + state->tile->lcu_offset_y * LCU_WIDTH) + mv[1]) * ref->width; // pre-calculate
        for (x = xpos; x < xpos + width; x++) {
          int x_in_lcu = (x & ((LCU_WIDTH)-1));

          lcu->rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] = ref->y[coord_y + (x + state->tile->lcu_offset_x * LCU_WIDTH) + mv[0]];
        }
      }

      if(!chroma_halfpel) {
        // Copy Chroma
        // TODO: chroma fractional pixel interpolation
        for (y = ypos>>1; y < (ypos + width)>>1; y++) {
          int y_in_lcu = (y & ((LCU_WIDTH>>1)-1));
          coord_y = ((y + state->tile->lcu_offset_y * (LCU_WIDTH>>1)) + (mv[1]>>1)) * ref_width_c; // pre-calculate
          for (x = xpos>>1; x < (xpos + width)>>1; x++) {
            int x_in_lcu = (x & ((LCU_WIDTH>>1)-1));
            lcu->rec.u[y_in_lcu*dst_width_c + x_in_lcu] = ref->u[coord_y + (x + state->tile->lcu_offset_x * (LCU_WIDTH>>1)) + (mv[0]>>1)];
            lcu->rec.v[y_in_lcu*dst_width_c + x_in_lcu] = ref->v[coord_y + (x + state->tile->lcu_offset_x * (LCU_WIDTH>>1)) + (mv[0]>>1)];
          }
        }
      }
    }
  }
}

/**
* \brief Reconstruct bi-pred inter block
* \param ref1 reference picture to copy the data from
* \param ref2 other reference picture to copy the data from
* \param xpos block x position
* \param ypos block y position
* \param width block width
* \param mv[2][2] motion vectors
* \param lcu destination lcu
* \returns Void
*/

void inter_recon_lcu_bipred(const encoder_state_t * const state, const kvz_picture * ref1, const kvz_picture * ref2, int32_t xpos, int32_t ypos, int32_t width, int16_t mv_param[2][2], lcu_t* lcu) {
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
  if (hi_prec_chroma_rec0) high_precision_rec0 = hi_prec_buf_t_alloc(LCU_WIDTH*LCU_WIDTH);
  if (hi_prec_chroma_rec1) high_precision_rec1 = hi_prec_buf_t_alloc(LCU_WIDTH*LCU_WIDTH);
  //Reconstruct both predictors
  inter_recon_lcu(state, ref1, xpos, ypos, width, mv_param[0], lcu, high_precision_rec0);
  if (!hi_prec_luma_rec0){
    memcpy(temp_lcu_y, lcu->rec.y, sizeof(kvz_pixel) * 64 * 64);
  }
  if (!hi_prec_chroma_rec0){
    memcpy(temp_lcu_u, lcu->rec.u, sizeof(kvz_pixel) * 32 * 32);
    memcpy(temp_lcu_v, lcu->rec.v, sizeof(kvz_pixel) * 32 * 32);
  }
  inter_recon_lcu(state, ref2, xpos, ypos, width, mv_param[1], lcu, high_precision_rec1);

  // After reconstruction, merge the predictors by taking an average of each pixel
  for (temp_y = 0; temp_y < width; ++temp_y) {
    int y_in_lcu = ((ypos + temp_y) & ((LCU_WIDTH)-1));
    for (temp_x = 0; temp_x < width; ++temp_x) {
      int x_in_lcu = ((xpos + temp_x) & ((LCU_WIDTH)-1));
      int16_t sample0_y = (hi_prec_luma_rec0 ? high_precision_rec0->y[y_in_lcu * LCU_WIDTH + x_in_lcu] : (temp_lcu_y[y_in_lcu * LCU_WIDTH + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      int16_t sample1_y = (hi_prec_luma_rec1 ? high_precision_rec1->y[y_in_lcu * LCU_WIDTH + x_in_lcu] : (lcu->rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      lcu->rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] = (kvz_pixel)fast_clip_32bit_to_pixel((sample0_y + sample1_y + offset) >> shift);
    }

  }
  for (temp_y = 0; temp_y < width >> 1; ++temp_y) {
    int y_in_lcu = (((ypos >> 1) + temp_y) & (LCU_WIDTH_C - 1));
    for (temp_x = 0; temp_x < width >> 1; ++temp_x) {
      int x_in_lcu = (((xpos >> 1) + temp_x) & (LCU_WIDTH_C - 1));
      int16_t sample0_u = (hi_prec_chroma_rec0 ? high_precision_rec0->u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (temp_lcu_u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      int16_t sample1_u = (hi_prec_chroma_rec1 ? high_precision_rec1->u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (lcu->rec.u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      lcu->rec.u[y_in_lcu * LCU_WIDTH_C + x_in_lcu] = (kvz_pixel)fast_clip_32bit_to_pixel((sample0_u + sample1_u + offset) >> shift);

      int16_t sample0_v = (hi_prec_chroma_rec0 ? high_precision_rec0->v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (temp_lcu_v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      int16_t sample1_v = (hi_prec_chroma_rec1 ? high_precision_rec1->v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] : (lcu->rec.v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] << (14 - KVZ_BIT_DEPTH)));
      lcu->rec.v[y_in_lcu * LCU_WIDTH_C + x_in_lcu] = (kvz_pixel)fast_clip_32bit_to_pixel((sample0_v + sample1_v + offset) >> shift);
    }
  }
  if (high_precision_rec0 != 0) hi_prec_buf_t_free(high_precision_rec0);
  if (high_precision_rec1 != 0) hi_prec_buf_t_free(high_precision_rec1);
}

/**
 * \brief Set unused L0/L1 motion vectors and reference
 * \param cu coding unit to clear
 */
static void inter_clear_cu_unused(cu_info_t* cu) {
  if(!(cu->inter.mv_dir & 1)) {
    cu->inter.mv[0][0] = 0;
    cu->inter.mv[0][1] = 0;
    cu->inter.mv_ref[0] = 255;
  }
  if(!(cu->inter.mv_dir & 2)) {
    cu->inter.mv[1][0] = 0;
    cu->inter.mv[1][1] = 0;
    cu->inter.mv_ref[1] = 255;
  }
}

/**
 * \brief Get merge candidates for current block
 * \param encoder encoder control struct to use
 * \param x_cu block x position in SCU
 * \param y_cu block y position in SCU
 * \param depth current block depth
 * \param b0 candidate b0
 * \param b1 candidate b1
 * \param b2 candidate b2
 * \param a0 candidate a0
 * \param a1 candidate a1
 */
void inter_get_spatial_merge_candidates(int32_t x, int32_t y, int8_t depth, cu_info_t **b0, cu_info_t **b1,
                                        cu_info_t **b2,cu_info_t **a0,cu_info_t **a1, lcu_t *lcu)
{
  uint8_t cur_block_in_scu = (LCU_WIDTH>>depth) / CU_MIN_SIZE_PIXELS; //!< the width of the current block on SCU
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
  int32_t x_cu = (x & (LCU_WIDTH - 1)) >> MAX_DEPTH; //!< coordinates from top-left of this LCU
  int32_t y_cu = (y & (LCU_WIDTH - 1)) >> MAX_DEPTH;
  cu_info_t* cu = &lcu->cu[LCU_CU_OFFSET];
  // A0 and A1 availability testing
  if (x != 0) {
    *a1 = &cu[x_cu - 1 + (y_cu + cur_block_in_scu - 1) * LCU_T_CU_WIDTH];
    if (!(*a1)->coded) *a1 = NULL;
    if(*a1) inter_clear_cu_unused(*a1);

    if (y_cu + cur_block_in_scu < LCU_WIDTH>>3) {
      *a0 = &cu[x_cu - 1 + (y_cu + cur_block_in_scu) * LCU_T_CU_WIDTH];
      if (!(*a0)->coded) *a0 = NULL;
    }
    if(*a0) inter_clear_cu_unused(*a0);
  }

  // B0, B1 and B2 availability testing
  if (y != 0) {
    if (x_cu + cur_block_in_scu < LCU_WIDTH>>3) {
      *b0 = &cu[x_cu + cur_block_in_scu + (y_cu - 1) * LCU_T_CU_WIDTH];
      if (!(*b0)->coded) *b0 = NULL;
    } else if(y_cu == 0) {
      // Special case, top-right cu from LCU is the last in lcu->cu array
      *b0 = &lcu->cu[LCU_T_CU_WIDTH*LCU_T_CU_WIDTH];
      if (!(*b0)->coded) *b0 = NULL;
    }
    if(*b0) inter_clear_cu_unused(*b0);

    *b1 = &cu[x_cu + cur_block_in_scu - 1 + (y_cu - 1) * LCU_T_CU_WIDTH];
    if (!(*b1)->coded) *b1 = NULL;
    if(*b1) inter_clear_cu_unused(*b1);

    if (x != 0) {
      *b2 = &cu[x_cu - 1 + (y_cu - 1) * LCU_T_CU_WIDTH];
      if(!(*b2)->coded) *b2 = NULL;
    }
    if(*b2) inter_clear_cu_unused(*b2);
  }
}

/**
 * \brief Get MV prediction for current block
 * \param encoder encoder control struct to use
 * \param x_cu block x position in SCU
 * \param y_cu block y position in SCU
 * \param depth current block depth
 * \param mv_pred[2][2] 2x motion vector prediction
 */
void inter_get_mv_cand(const encoder_state_t * const state, int32_t x, int32_t y, int8_t depth, int16_t mv_cand[2][2], cu_info_t* cur_cu, lcu_t *lcu, int8_t reflist)
{
  uint8_t candidates = 0;
  uint8_t b_candidates = 0;
  int8_t reflist2nd = !reflist;

  cu_info_t *b0, *b1, *b2, *a0, *a1;
  b0 = b1 = b2 = a0 = a1 = NULL;
  inter_get_spatial_merge_candidates(x, y, depth, &b0, &b1, &b2, &a0, &a1, lcu);

 #define CALCULATE_SCALE(cu,tb,td) ((tb * ((0x4000 + (abs(td)>>1))/td) + 32) >> 6)
#define APPLY_MV_SCALING(cu, cand, list) {int td = state->global->poc - state->global->ref->pocs[(cu)->inter.mv_ref[list]];\
                                   int tb = state->global->poc - state->global->ref->pocs[cur_cu->inter.mv_ref[reflist]];\
                                   if (td != tb) { \
                                      int scale = CALCULATE_SCALE(cu,tb,td); \
                                       mv_cand[cand][0] = ((scale * (cu)->inter.mv[list][0] + 127 + (scale * (cu)->inter.mv[list][0] < 0)) >> 8 ); \
                                       mv_cand[cand][1] = ((scale * (cu)->inter.mv[list][1] + 127 + (scale * (cu)->inter.mv[list][1] < 0)) >> 8 ); }}

  // Left predictors
  if (a0 && a0->type == CU_INTER && (
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
  } else if (a1 && a1->type == CU_INTER && (
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
    if (a0 && a0->type == CU_INTER) {
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
    } else if (a1 && a1->type == CU_INTER) {
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
  if (b0 && b0->type == CU_INTER && (
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
  } else if (b1 && b1->type == CU_INTER && (
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
  } else if (b2 && b2->type == CU_INTER && (
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
  if((a1 && a1->type == CU_INTER) || (a0 && a0->type == CU_INTER)) {
    b_candidates = 1;
  } else if(candidates != 2) {
    b_candidates = 0;
  }

  if(!b_candidates) {
    // Top predictors
    if (b0 && b0->type == CU_INTER) {
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
    } else if (b1 && b1->type == CU_INTER) {
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
    } else if(b2 && b2->type == CU_INTER) {
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

#if ENABLE_TEMPORAL_MVP
  if(candidates < AMVP_MAX_NUM_CANDS) {
    //TODO: add temporal mv predictor
  }
#endif

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
 * \brief Get merge predictions for current block
 * \param encoder encoder control struct to use
 * \param x_cu block x position in SCU
 * \param y_cu block y position in SCU
 * \param depth current block depth
 * \param mv_pred[MRG_MAX_NUM_CANDS][2] MRG_MAX_NUM_CANDS motion vector prediction
 */
uint8_t inter_get_merge_cand(const encoder_state_t * const state, int32_t x, int32_t y, int8_t depth, inter_merge_cand_t mv_cand[MRG_MAX_NUM_CANDS], lcu_t *lcu)
{
  uint8_t candidates = 0;
  int8_t duplicate = 0;

  cu_info_t *b0, *b1, *b2, *a0, *a1;
  int8_t zero_idx = 0;
  b0 = b1 = b2 = a0 = a1 = NULL;
  inter_get_spatial_merge_candidates(x, y, depth, &b0, &b1, &b2, &a0, &a1, lcu);


#define CHECK_DUPLICATE(CU1,CU2) {duplicate = 0; if ((CU2) && (CU2)->type == CU_INTER && \
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

  if (a1 && a1->type == CU_INTER) {
    mv_cand[candidates].mv[0][0] = a1->inter.mv[0][0];
    mv_cand[candidates].mv[0][1] = a1->inter.mv[0][1];
    mv_cand[candidates].mv[1][0] = a1->inter.mv[1][0];
    mv_cand[candidates].mv[1][1] = a1->inter.mv[1][1];
    mv_cand[candidates].ref[0] = a1->inter.mv_ref[0];
    mv_cand[candidates].ref[1] = a1->inter.mv_ref[1];
    mv_cand[candidates].dir = a1->inter.mv_dir;
    candidates++;
  }

  if (b1 && b1->type == CU_INTER) {
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

  if (b0 && b0->type == CU_INTER) {
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

  if (a0 && a0->type == CU_INTER) {
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
    if(b2 && b2->type == CU_INTER) {
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

#if ENABLE_TEMPORAL_MVP
  if(candidates < AMVP_MAX_NUM_CANDS) {
    //TODO: add temporal mv predictor
  }
#endif

  if (candidates < MRG_MAX_NUM_CANDS && state->global->slicetype == SLICE_B) {
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

  int num_ref = state->global->ref->used_size;

  if (candidates < MRG_MAX_NUM_CANDS && state->global->slicetype == SLICE_B) {
    int j;
    int ref_negative = 0;
    int ref_positive = 0;
    for (j = 0; j < state->global->ref->used_size; j++) {
      if (state->global->ref->pocs[j] < state->global->poc) {
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
    if (state->global->slicetype == SLICE_B) {
      mv_cand[candidates].mv[1][0] = 0;
      mv_cand[candidates].mv[1][1] = 0;
      mv_cand[candidates].dir = 3;
    }
    zero_idx++;
    candidates++;
  }

  return candidates;
}
