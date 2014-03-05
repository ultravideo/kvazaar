/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2014 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
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

/**
 * \brief Set block info to the CU structure
 * \param pic picture to use
 * \param x_cu x CU position (smallest CU)
 * \param y_cu y CU position (smallest CU)
 * \param depth current CU depth
 * \param cur_cu CU to take the settings from
 * \returns Void
*/
void inter_set_block(picture* pic, uint32_t x_cu, uint32_t y_cu, uint8_t depth, cu_info* cur_cu)
{
  uint32_t x, y;
  // Width in smallest CU
  int width_in_scu = pic->width_in_lcu<<MAX_DEPTH;
  int block_scu_width = (LCU_WIDTH>>depth)/(LCU_WIDTH>>MAX_DEPTH);
  int tr_depth = (depth == 0 ? 1 : depth);
  // Loop through all the block in the area of cur_cu
  for (y = y_cu; y < y_cu + block_scu_width; y++) {
    int cu_pos = y * width_in_scu; //!< calculate y-position once, use with every x
    for (x = x_cu; x < x_cu + block_scu_width; x++) {
      // Set all SCU's to this blocks values at the bottom most depth.
      pic->cu_array[MAX_DEPTH][cu_pos + x].depth = depth;
      pic->cu_array[MAX_DEPTH][cu_pos + x].type  = CU_INTER;
      pic->cu_array[MAX_DEPTH][cu_pos + x].part_size = SIZE_2Nx2N;
      pic->cu_array[MAX_DEPTH][cu_pos + x].inter.mode   = cur_cu->inter.mode;
      pic->cu_array[MAX_DEPTH][cu_pos + x].inter.mv[0]  = cur_cu->inter.mv[0];
      pic->cu_array[MAX_DEPTH][cu_pos + x].inter.mv[1]  = cur_cu->inter.mv[1];
      pic->cu_array[MAX_DEPTH][cu_pos + x].inter.mv_dir = cur_cu->inter.mv_dir;
      pic->cu_array[MAX_DEPTH][cu_pos + x].inter.mv_ref = cur_cu->inter.mv_ref;
      pic->cu_array[MAX_DEPTH][cu_pos + x].tr_depth = tr_depth;
    }
  }
}

/**
 * \brief Reconstruct inter block
 * \param ref picture to copy the data from
 * \param xpos block x position
 * \param ypos block y position
 * \param width block width
 * \param mv[2] motion vector
 * \param lcu destination lcu
 * \returns Void
*/
void inter_recon_lcu(picture* ref,int32_t xpos, int32_t ypos,int32_t width, const int16_t mv_param[2], lcu_t *lcu)
{
  int x,y,coord_x,coord_y;
  int16_t mv[2] = { mv_param[0], mv_param[1] };

  int32_t dst_width_c = LCU_WIDTH>>1; //!< Destination picture width in chroma pixels
  int32_t ref_width_c = ref->width>>1; //!< Reference picture width in chroma pixels

  // negative overflow flag
  int8_t overflow_neg_x = (xpos + (mv[0]>>2) < 0)?1:0;
  int8_t overflow_neg_y = (ypos + (mv[1]>>2) < 0)?1:0;

  // positive overflow flag
  int8_t overflow_pos_x = (xpos + (mv[0]>>2) + width > ref->width )?1:0;
  int8_t overflow_pos_y = (ypos + (mv[1]>>2) + width > ref->height)?1:0;

  // Chroma half-pel
  #define HALFPEL_CHROMA_WIDTH ((LCU_WIDTH>>1) + 8)
  int8_t chroma_halfpel = ((mv[0]>>2)&1) || ((mv[1]>>2)&1); //!< (luma integer mv) lsb is set -> chroma is half-pel
  int16_t halfpel_src_u[HALFPEL_CHROMA_WIDTH * HALFPEL_CHROMA_WIDTH]; //!< U source block for interpolation
  int16_t halfpel_src_v[HALFPEL_CHROMA_WIDTH * HALFPEL_CHROMA_WIDTH]; //!< V source block for interpolation
  int16_t *halfpel_src_off_u = &halfpel_src_u[HALFPEL_CHROMA_WIDTH*4 + 4]; //!< halfpel_src_u with offset (4,4)
  int16_t *halfpel_src_off_v = &halfpel_src_v[HALFPEL_CHROMA_WIDTH*4 + 4]; //!< halfpel_src_v with offset (4,4)
  int16_t halfpel_u[LCU_WIDTH * LCU_WIDTH]; //!< interpolated 2W x 2H block (u)
  int16_t halfpel_v[LCU_WIDTH * LCU_WIDTH]; //!< interpolated 2W x 2H block (v)

  // TODO: Fractional pixel support
  mv[0] >>= 2;
  mv[1] >>= 2;

  // Chroma half-pel
  // get half-pel interpolated block and push it to output
  if(chroma_halfpel) {

    int halfpel_y, halfpel_x;
    int abs_mv_x = mv[0]&1;
    int abs_mv_y = mv[1]&1;
    int8_t overflow_neg_y_temp,overflow_pos_y_temp,overflow_neg_x_temp,overflow_pos_x_temp;
    // Fill source blocks with data from reference, -4...width+4
    for (halfpel_y = 0, y = (ypos>>1) - 4; y < ((ypos + width)>>1) + 4; halfpel_y++, y++) {
      // calculate y-pixel offset
      coord_y = y + (mv[1]>>1);

      // On y-overflow set coord_y accordingly
      overflow_neg_y_temp = (coord_y < 0) ? 1 : 0;
      overflow_pos_y_temp = (coord_y >= ref->height>>1) ? 1 : 0;
      if (overflow_neg_y_temp)      coord_y = 0;
      else if (overflow_pos_y_temp) coord_y = (ref->height>>1) - 1;
      coord_y *= ref_width_c;

      for (halfpel_x = 0, x = (xpos>>1) - 4; x < ((xpos + width)>>1) + 4; halfpel_x++, x++) {
        coord_x = x + (mv[0]>>1);

        // On x-overflow set coord_x accordingly
        overflow_neg_x_temp = (coord_x < 0) ? 1 : 0;
        overflow_pos_x_temp = (coord_x >= ref_width_c) ? 1 : 0;
        if (overflow_neg_x_temp)      coord_x = 0;
        else if (overflow_pos_x_temp) coord_x = ref_width_c - 1;

        // Store source block data (with extended borders)
        halfpel_src_u[halfpel_y*HALFPEL_CHROMA_WIDTH + halfpel_x] = ref->u_recdata[coord_y + coord_x];
        halfpel_src_v[halfpel_y*HALFPEL_CHROMA_WIDTH + halfpel_x] = ref->v_recdata[coord_y + coord_x];
      }
    }

    // Filter the block to half-pel resolution
    filter_inter_halfpel_chroma(halfpel_src_off_u, HALFPEL_CHROMA_WIDTH, width>>1, width>>1, halfpel_u, LCU_WIDTH, abs_mv_x, abs_mv_y);
    filter_inter_halfpel_chroma(halfpel_src_off_v, HALFPEL_CHROMA_WIDTH, width>>1, width>>1, halfpel_v, LCU_WIDTH, abs_mv_x, abs_mv_y);

    // Assign filtered pixels to output, take every second half-pel sample with offset of abs_mv_y/x
    for (halfpel_y = abs_mv_y, y = ypos>>1; y < (ypos + width)>>1; halfpel_y += 2, y++) {
      for (halfpel_x = abs_mv_x, x = xpos>>1; x < (xpos + width)>>1; halfpel_x += 2, x++) {
        int x_in_lcu = (x & ((LCU_WIDTH>>1)-1));
        int y_in_lcu = (y & ((LCU_WIDTH>>1)-1));
        lcu->rec.u[y_in_lcu*dst_width_c + x_in_lcu] = (uint8_t)halfpel_u[halfpel_y*LCU_WIDTH + halfpel_x];
        lcu->rec.v[y_in_lcu*dst_width_c + x_in_lcu] = (uint8_t)halfpel_v[halfpel_y*LCU_WIDTH + halfpel_x];
      }
    }
  }

  // With overflow present, more checking
  if (overflow_neg_x || overflow_neg_y || overflow_pos_x || overflow_pos_y) {
    // Copy Luma with boundary checking
    for (y = ypos; y < ypos + width; y++) {
      for (x = xpos; x < xpos + width; x++) {
        int x_in_lcu = (x & ((LCU_WIDTH)-1));
        int y_in_lcu = (y & ((LCU_WIDTH)-1));

        coord_x = x + mv[0];
        coord_y = y + mv[1];
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
        lcu->rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] = ref->y_recdata[coord_y*ref->width + coord_x];
      }
    }

    if(!chroma_halfpel) {
      // Copy Chroma with boundary checking
      // TODO: chroma fractional pixel interpolation
      for (y = ypos>>1; y < (ypos + width)>>1; y++) {
        for (x = xpos>>1; x < (xpos + width)>>1; x++) {
          int x_in_lcu = (x & ((LCU_WIDTH>>1)-1));
          int y_in_lcu = (y & ((LCU_WIDTH>>1)-1));

          coord_x = x + (mv[0]>>1);
          coord_y = y + (mv[1]>>1);

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
          lcu->rec.u[y_in_lcu*dst_width_c + x_in_lcu] = ref->u_recdata[coord_y*ref_width_c + coord_x];
          lcu->rec.v[y_in_lcu*dst_width_c + x_in_lcu] = ref->v_recdata[coord_y*ref_width_c + coord_x];

        }
      }
    }
  } else { //If no overflow, we can copy without checking boundaries
    // Copy Luma
    for (y = ypos; y < ypos + width; y++) {
      int y_in_lcu = (y & ((LCU_WIDTH)-1));
      coord_y = (y + mv[1]) * ref->width; // pre-calculate
      for (x = xpos; x < xpos + width; x++) {
        int x_in_lcu = (x & ((LCU_WIDTH)-1));
        
        lcu->rec.y[y_in_lcu * LCU_WIDTH + x_in_lcu] = ref->y_recdata[coord_y + x + mv[0]];
      }
    }

    if(!chroma_halfpel) {
      // Copy Chroma
      // TODO: chroma fractional pixel interpolation
      for (y = ypos>>1; y < (ypos + width)>>1; y++) {
        int y_in_lcu = (y & ((LCU_WIDTH>>1)-1));
        coord_y = (y + (mv[1]>>1)) * ref_width_c; // pre-calculate
        for (x = xpos>>1; x < (xpos + width)>>1; x++) {
          int x_in_lcu = (x & ((LCU_WIDTH>>1)-1));
          lcu->rec.u[y_in_lcu*dst_width_c + x_in_lcu] = ref->u_recdata[coord_y + x + (mv[0]>>1)];
          lcu->rec.v[y_in_lcu*dst_width_c + x_in_lcu] = ref->v_recdata[coord_y + x + (mv[0]>>1)];
        }
      }
    }
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
void inter_get_spatial_merge_candidates(int32_t x, int32_t y, int8_t depth, cu_info **b0, cu_info **b1,
                                        cu_info **b2,cu_info **a0,cu_info **a1, lcu_t *lcu)
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
  int32_t x_cu = (x & LCU_WIDTH-1) >> MAX_DEPTH; //!< coordinates from top-left of this LCU
  int32_t y_cu = (y & LCU_WIDTH-1) >> MAX_DEPTH;
  cu_info* cu = &lcu->cu[LCU_CU_OFFSET];
  // A0 and A1 availability testing
  if (x != 0) {
    *a1 = &cu[x_cu - 1 + (y_cu + cur_block_in_scu - 1) * LCU_T_CU_WIDTH];
    if (!(*a1)->coded) *a1 = NULL;

    if (y_cu + cur_block_in_scu < LCU_WIDTH>>3) {
      *a0 = &cu[x_cu - 1 + (y_cu + cur_block_in_scu) * LCU_T_CU_WIDTH];
      if (!(*a0)->coded) *a0 = NULL;
    }
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
    
    *b1 = &cu[x_cu + cur_block_in_scu - 1 + (y_cu - 1) * LCU_T_CU_WIDTH];
    if (!(*b1)->coded) *b1 = NULL;

    if (x != 0) {
      *b2 = &cu[x_cu - 1 + (y_cu - 1) * LCU_T_CU_WIDTH];
      if(!(*b2)->coded) *b2 = NULL;
    }
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
void inter_get_mv_cand(encoder_control *encoder, int32_t x, int32_t y, int8_t depth, int16_t mv_cand[2][2], cu_info* cur_cu, lcu_t *lcu)
{
  uint8_t candidates = 0;
  uint8_t b_candidates = 0;

  cu_info *b0, *b1, *b2, *a0, *a1;
  b0 = b1 = b2 = a0 = a1 = NULL;
  inter_get_spatial_merge_candidates(x, y, depth, &b0, &b1, &b2, &a0, &a1, lcu);

 #define CALCULATE_SCALE(cu,tb,td) ((tb * ((0x4000 + (abs(td)>>1))/td) + 32) >> 6)
#define APPLY_MV_SCALING(cu, cand) {int td = encoder->poc - encoder->ref->pics[(cu)->inter.mv_ref]->poc;\
                                   int tb = encoder->poc - encoder->ref->pics[cur_cu->inter.mv_ref]->poc;\
                                   if (td != tb) { \
                                      int scale = CALCULATE_SCALE(cu,tb,td); \
                                       mv_cand[cand][0] = ((scale * (cu)->inter.mv[0] + 127 + (scale * (cu)->inter.mv[0] < 0)) >> 8 ); \
                                       mv_cand[cand][1] = ((scale * (cu)->inter.mv[1] + 127 + (scale * (cu)->inter.mv[1] < 0)) >> 8 ); }}

  // Left predictors
  if (a0 && a0->type == CU_INTER && a0->inter.mv_ref == cur_cu->inter.mv_ref) {
    mv_cand[candidates][0] = a0->inter.mv[0];
    mv_cand[candidates][1] = a0->inter.mv[1];
    candidates++;
  } else if (a1 && a1->type == CU_INTER && a1->inter.mv_ref == cur_cu->inter.mv_ref) {
    mv_cand[candidates][0] = a1->inter.mv[0];
    mv_cand[candidates][1] = a1->inter.mv[1];
    candidates++;
  }

  if(!candidates) {
      // Left predictors
    if (a0 && a0->type == CU_INTER) {
      mv_cand[candidates][0] = a0->inter.mv[0];
      mv_cand[candidates][1] = a0->inter.mv[1];
      APPLY_MV_SCALING(a0, candidates);
      candidates++;
    } else if (a1 && a1->type == CU_INTER) {
      mv_cand[candidates][0] = a1->inter.mv[0];
      mv_cand[candidates][1] = a1->inter.mv[1];
      APPLY_MV_SCALING(a1, candidates);
      candidates++;
    }
  }

  // Top predictors
  if (b0 && b0->type == CU_INTER && b0->inter.mv_ref == cur_cu->inter.mv_ref) {
    mv_cand[candidates][0] = b0->inter.mv[0];
    mv_cand[candidates][1] = b0->inter.mv[1];
    b_candidates++;
  } else if (b1 && b1->type == CU_INTER && b1->inter.mv_ref == cur_cu->inter.mv_ref) {
    mv_cand[candidates][0] = b1->inter.mv[0];
    mv_cand[candidates][1] = b1->inter.mv[1];
    b_candidates++;
  } else if(b2 && b2->type == CU_INTER && b2->inter.mv_ref == cur_cu->inter.mv_ref) {
    mv_cand[candidates][0] = b2->inter.mv[0];
    mv_cand[candidates][1] = b2->inter.mv[1];
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
      mv_cand[candidates][0] = b0->inter.mv[0];
      mv_cand[candidates][1] = b0->inter.mv[1];
      APPLY_MV_SCALING(b0, candidates);
      candidates++;
    } else if (b1 && b1->type == CU_INTER) {
      mv_cand[candidates][0] = b1->inter.mv[0];
      mv_cand[candidates][1] = b1->inter.mv[1];
      APPLY_MV_SCALING(b1, candidates);
      candidates++;
    } else if(b2 && b2->type == CU_INTER) {
      mv_cand[candidates][0] = b2->inter.mv[0];
      mv_cand[candidates][1] = b2->inter.mv[1];
      APPLY_MV_SCALING(b2, candidates);
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
uint8_t inter_get_merge_cand(int32_t x, int32_t y, int8_t depth, int16_t mv_cand[MRG_MAX_NUM_CANDS][3], cu_info* cur_cu, lcu_t *lcu)
{
  uint8_t candidates = 0;
  int8_t duplicate = 0;

  cu_info *b0, *b1, *b2, *a0, *a1;
  int8_t zero_idx = 0;
  b0 = b1 = b2 = a0 = a1 = NULL;
  inter_get_spatial_merge_candidates(x, y, depth, &b0, &b1, &b2, &a0, &a1, lcu);


#define CHECK_DUPLICATE(CU1,CU2) {duplicate = 0; if ((CU2) && (CU2)->type == CU_INTER && \
                                                     (CU1)->inter.mv[0] == (CU2)->inter.mv[0] && \
                                                     (CU1)->inter.mv[1] == (CU2)->inter.mv[1] && \
                                                     (CU1)->inter.mv_ref == (CU2)->inter.mv_ref) duplicate = 1; }

  if (a1 && a1->type == CU_INTER) {
      mv_cand[candidates][0] = a1->inter.mv[0];
      mv_cand[candidates][1] = a1->inter.mv[1];
      mv_cand[candidates][2] = a1->inter.mv_ref;
      candidates++;
  }

  if (b1 && b1->type == CU_INTER) {
    if(candidates) CHECK_DUPLICATE(b1, a1);
    if(!duplicate) {
      mv_cand[candidates][0] = b1->inter.mv[0];
      mv_cand[candidates][1] = b1->inter.mv[1];
      mv_cand[candidates][2] = b1->inter.mv_ref;
      candidates++;
    }
  }

  if (b0 && b0->type == CU_INTER) {
    if(candidates) CHECK_DUPLICATE(b0,b1);
    if(!duplicate) {
      mv_cand[candidates][0] = b0->inter.mv[0];
      mv_cand[candidates][1] = b0->inter.mv[1];
      mv_cand[candidates][2] = b0->inter.mv_ref;
      candidates++;
    }
  }

  if (a0 && a0->type == CU_INTER) {
    if(candidates) CHECK_DUPLICATE(a0,a1);
    if(!duplicate) {
      mv_cand[candidates][0] = a0->inter.mv[0];
      mv_cand[candidates][1] = a0->inter.mv[1];
      mv_cand[candidates][2] = a0->inter.mv_ref;
      candidates++;
    }
  }

  if (candidates != 4) {
    if(b2 && b2->type == CU_INTER) {
      CHECK_DUPLICATE(b2,a1);
      if(!duplicate) {
        CHECK_DUPLICATE(b2,b1);
        if(!duplicate) {
          mv_cand[candidates][0] = b2->inter.mv[0];
          mv_cand[candidates][1] = b2->inter.mv[1];
          mv_cand[candidates][2] = b2->inter.mv_ref;
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

  // Add (0,0) prediction
  if (candidates != 5) {
    mv_cand[candidates][0] = 0;
    mv_cand[candidates][1] = 0;
    mv_cand[candidates][2] = zero_idx;
    zero_idx++;
    candidates++;
  }

  return candidates;
}
