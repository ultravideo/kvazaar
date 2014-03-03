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

/**
 * \file
 * \brief Functions for handling intra frames.
 */

#include "intra.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "encoder.h"


const uint8_t intra_hor_ver_dist_thres[5] = {0,7,1,0,0};


/**
 * \brief Set intrablock mode (and init typedata)
 * \param pic picture to use
 * \param xCtb x CU position (smallest CU)
 * \param yCtb y CU position (smallest CU)
 * \param depth current CU depth
 * \param mode mode to set
 * \returns Void
 */
void intra_set_block_mode(picture *pic,uint32_t x_cu, uint32_t y_cu, uint8_t depth, uint8_t mode, uint8_t part_mode)
{
  uint32_t x, y;
  int width_in_scu = pic->width_in_lcu<<MAX_DEPTH; //!< Width in smallest CU
  int block_scu_width = (LCU_WIDTH>>depth)/(LCU_WIDTH>>MAX_DEPTH);

  if (part_mode == SIZE_NxN) {
    cu_info *cur_cu = &pic->cu_array[MAX_DEPTH][x_cu + y_cu * width_in_scu];
    // Modes are already set.
    cur_cu->depth = depth;
    cur_cu->type = CU_INTRA;
    cur_cu->tr_depth = depth + 1;
    return;
  }

  // Loop through all the blocks in the area of cur_cu
  for (y = y_cu; y < y_cu + block_scu_width; y++) {
    for (x = x_cu; x < x_cu + block_scu_width; x++) {
      cu_info *cur_cu = &pic->cu_array[MAX_DEPTH][x + y * width_in_scu];
      cur_cu->depth = depth;
      cur_cu->type = CU_INTRA;
      cur_cu->intra[0].mode = mode;
      cur_cu->intra[1].mode = mode;
      cur_cu->intra[2].mode = mode;
      cur_cu->intra[3].mode = mode;
      cur_cu->part_size = part_mode;
      cur_cu->tr_depth = depth;
    }
  }
}

/**
 * \brief get intrablock mode
 * \param pic picture data to use
 * \param picwidth width of the picture data
 * \param xpos x-position
 * \param ypos y-position
 * \param width block width
 * \returns DC prediction
*/
pixel intra_get_dc_pred(pixel *pic, uint16_t picwidth, uint8_t width)
{
  int32_t i, sum = 0;

  // pixels on top and left
  for (i = -picwidth; i < width - picwidth; i++) {
    sum += pic[i];
  }
  for (i = -1; i < width * picwidth - 1; i += picwidth) {
    sum += pic[i];
  }

  // return the average
  return (pixel)((sum + width) / (width + width));
}

#define PU_INDEX(x_pu, y_pu) (((x_pu) % 2)  + 2 * ((y_pu) % 2))

/**
 * \brief Function for deriving intra luma predictions
 * \param pic picture to use
 * \param x_cu x CU position (smallest CU)
 * \param y_cu y CU position (smallest CU)
 * \param preds output buffer for 3 predictions
 * \returns (predictions are found)?1:0
 */
int8_t intra_get_dir_luma_predictor(uint32_t x, uint32_t y, int8_t* preds,
                                    cu_info* cur_cu, cu_info* left_cu, cu_info* above_cu)
{
  int x_cu = x>>3;
  int y_cu = y>>3;

  // The default mode if block is not coded yet is INTRA_DC.
  int8_t left_intra_dir  = 1;
  int8_t above_intra_dir = 1;
    
  if (cur_cu->part_size == SIZE_NxN && (x & 7) == 1) {
    // If current CU is NxN and PU is on the right half, take mode from the
    // left half of the same CU.
    left_intra_dir = cur_cu->intra[PU_INDEX(0, y_cu<<1)].mode;
  } else if (left_cu && left_cu->type == CU_INTRA) {
    // Otherwise take the mode from the right side of the CU on the left.
    left_intra_dir = left_cu->intra[PU_INDEX(1, y_cu<<1)].mode;
  }

  if (cur_cu->part_size == SIZE_NxN && (y & 7) == 1) {
    // If current CU is NxN and PU is on the bottom half, take mode from the
    // top half of the same CU.
    above_intra_dir = cur_cu->intra[PU_INDEX(x_cu<<1, 0)].mode;
  } else if (above_cu && above_cu->type == CU_INTRA &&
             (y_cu * (LCU_WIDTH>>MAX_DEPTH)) % LCU_WIDTH != 0)
  {
    // Otherwise take the mode from the bottom half of the CU above.
    above_intra_dir = above_cu->intra[PU_INDEX(x_cu<<1, 1)].mode;
  }

  // If the predictions are the same, add new predictions
  if (left_intra_dir == above_intra_dir) {
    if (left_intra_dir > 1) { // angular modes
      preds[0] = left_intra_dir;
      preds[1] = ((left_intra_dir + 29) % 32) + 2;
      preds[2] = ((left_intra_dir - 1 ) % 32) + 2;
    } else { //non-angular
      preds[0] = 0;//PLANAR_IDX;
      preds[1] = 1;//DC_IDX;
      preds[2] = 26;//VER_IDX;
    }
  } else { // If we have two distinct predictions
    preds[0] = left_intra_dir;
    preds[1] = above_intra_dir;

    // add planar mode if it's not yet present
    if (left_intra_dir && above_intra_dir ) {
      preds[2] = 0; // PLANAR_IDX;
    } else { // else we add 26 or 1
      preds[2] =  (left_intra_dir+above_intra_dir)<2? 26 : 1;
    }
  }

  return 1;
}

/**
 * \brief Intra filtering of the border samples
 * \param ref reference picture data
 * \param x_cu x CU position (smallest CU)
 * \param y_cu y CU position (smallest CU)
 * \param depth current CU depth
 * \param preds output buffer for 3 predictions
 * \returns (predictions are found)?1:0
 */
void intra_filter(pixel *ref, int32_t stride,int32_t width, int8_t mode)
{
  #define FWIDTH (LCU_WIDTH*2+1)
  pixel filtered[FWIDTH * FWIDTH]; //!< temporary buffer for filtered samples
  pixel *filteredShift = &filtered[FWIDTH+1]; //!< pointer to temporary buffer with offset (1,1)
  int x,y;

  if (!mode) {
    // pF[ -1 ][ -1 ] = ( p[ -1 ][ 0 ] + 2*p[ -1 ][ -1 ] + p[ 0 ][ -1 ] + 2 )  >>  2	(8 35)
    filteredShift[-FWIDTH-1] = (ref[-1] + 2*ref[-(int32_t)stride-1] + ref[-(int32_t)stride] + 2) >> 2;

    // pF[ -1 ][ y ] = ( p[ -1 ][ y + 1 ] + 2*p[ -1 ][ y ] + p[ -1 ][ y - 1 ] + 2 )  >>  2 for y = 0..nTbS * 2 - 2	(8 36)
    for (y = 0; y < (int32_t)width * 2 - 1; y++) {
      filteredShift[y*FWIDTH-1] = (ref[(y + 1) * stride - 1] + 2*ref[y * stride - 1] + ref[(y - 1) * stride - 1] + 2) >> 2;
    }

    // pF[ -1 ][ nTbS * 2 - 1 ] = p[ -1 ][ nTbS * 2 - 1 ]		(8 37)
    filteredShift[(width * 2 - 1) * FWIDTH - 1] = ref[(width * 2 - 1) * stride - 1];

    // pF[ x ][ -1 ] = ( p[ x - 1 ][ -1 ] + 2*p[ x ][ -1 ] + p[ x + 1 ][ -1 ] + 2 )  >>  2 for x = 0..nTbS * 2 - 2	(8 38)
    for(x = 0; x < (int32_t)width*2-1; x++) {
      filteredShift[x - FWIDTH] = (ref[x - 1 - stride] + 2*ref[x - stride] + ref[x + 1 - stride] + 2) >> 2;
    }

    // pF[ nTbS * 2 - 1 ][ -1 ] = p[ nTbS * 2 - 1 ][ -1 ]
    filteredShift[(width * 2 - 1) - FWIDTH] = ref[(width * 2 - 1) - stride];

    // Copy filtered samples to the input array
    for (x = -1; x < (int32_t)width * 2; x++) {
      ref[x - stride] = filtered[x + 1];
    }
    for(y = 0; y < (int32_t)width * 2; y++)  {
      ref[y * stride - 1] = filtered[(y + 1) * FWIDTH];
    }

  } else  {
    printf("UNHANDLED: %s: %d\r\n", __FILE__, __LINE__);
    exit(1);
  }
  #undef FWIDTH
}

/** 
 * \brief Helper function to find intra merge costs
 * \returns intra mode coding cost in bits
 */
static uint32_t intra_pred_ratecost(int16_t mode, int8_t *intra_preds)
{
   // merge mode -1 means they are not used -> cost 0
   if(intra_preds[0] == -1) return 0;

   // First candidate needs only one bit and two other need two
   if(intra_preds[0] == mode) {
     return 1; 
   } else if(intra_preds[1] == mode || intra_preds[2] == mode) {
     return 2;
   }
   // Without merging the cost is 5 bits
   return 5;
}

/**
 * \brief Function to test best intra prediction mode
 * \param orig original picture data
 * \param origstride original picture stride
 * \param rec reconstructed picture data
 * \param recstride reconstructed picture stride
 * \param xpos source x-position
 * \param ypos source y-position
 * \param width block size to predict
 * \param dst destination buffer for best prediction
 * \param dststride destination width
 * \param sad_out sad value of best mode
 * \returns best intra mode

 This function derives the prediction samples for planar mode (intra coding).
*/
int16_t intra_prediction(pixel *orig, int32_t origstride, pixel *rec, int16_t recstride, uint16_t xpos,
                         uint16_t ypos, uint8_t width, pixel *dst, int32_t dststride, uint32_t *sad_out, int8_t *intra_preds)
{
  uint32_t best_sad = 0xffffffff;
  uint32_t sad = 0;
  int16_t best_mode = 1;
  int32_t x,y;
  int16_t i;
  uint32_t ratecost = 0;

  cost_16bit_nxn_func cost_func = get_sad_16bit_nxn_func(width);

  // Temporary block arrays
  // TODO: alloc with alignment
  pixel pred[LCU_WIDTH * LCU_WIDTH + 1];
  pixel orig_block[LCU_WIDTH * LCU_WIDTH + 1];
  pixel rec_filtered_temp[(LCU_WIDTH * 2 + 8) * (LCU_WIDTH * 2 + 8) + 1];

  pixel* rec_filtered = &rec_filtered_temp[recstride + 1]; //!< pointer to rec_filtered_temp with offset of (1,1)
  int8_t filter = (width<32); // TODO: chroma support

  uint8_t threshold = intra_hor_ver_dist_thres[g_to_bits[width]]; //!< Intra filtering threshold

  #define COPY_PRED_TO_DST() for (y = 0; y < (int32_t)width; y++)  { for (x = 0; x < (int32_t)width; x++) { dst[x + y*dststride] = pred[x + y*width]; } }
  #define CHECK_FOR_BEST(mode, additional_sad)  sad = cost_func(pred, orig_block); \
                                                sad += additional_sad;\
                                                if(sad < best_sad)\
                                                {\
                                                  best_sad = sad;\
                                                  best_mode = mode;\
                                                  COPY_PRED_TO_DST();\
                                                }

  // Store original block for SAD computation
  i = 0;
  for(y = 0; y < (int32_t)width; y++) {
    for(x = 0; x < (int32_t)width; x++) {
      orig_block[i++] = orig[x + y*origstride];
    }
  }

  // Filtered only needs the borders
  for (y = -1; y < (int32_t)recstride; y++) {
    rec_filtered[y*recstride - 1] = rec[y*recstride - 1];
  }
  for (x = 0; x < (int32_t)recstride; x++) {
    rec_filtered[x - recstride] = rec[x - recstride];
  }
  // Apply filter
  intra_filter(rec_filtered,recstride,width,0);

  // Test DC mode (never filtered)
  {
    pixel val = intra_get_dc_pred(rec, recstride, width);
    for (i = 0; i < (int32_t)(width*width); i++) {
      pred[i] = val;
    }
    ratecost = intra_pred_ratecost(1,intra_preds)*(int)(g_cur_lambda_cost+0.5);
    CHECK_FOR_BEST(1,ratecost);
  }

  // Check angular not requiring filtering
  for (i = 2; i < 35; i++) {
    int distance = MIN(abs(i - 26),abs(i - 10)); //!< Distance from top and left predictions
    if(distance <= threshold) {
      intra_get_angular_pred(rec, recstride, pred, width, width, i, filter);
      ratecost = intra_pred_ratecost(i,intra_preds)*(int)(g_cur_lambda_cost+0.5);
      CHECK_FOR_BEST(i,ratecost);
    }
  }

  // FROM THIS POINT FORWARD, USING FILTERED PREDICTION

  // Test planar mode (always filtered)
  intra_get_planar_pred(rec_filtered, recstride, width, pred, width);
  ratecost = intra_pred_ratecost(0,intra_preds)*(int)(g_cur_lambda_cost+0.5);
  CHECK_FOR_BEST(0,ratecost);

  // Check angular predictions which require filtered samples
  // TODO: add conditions to skip some modes on borders
  // chroma can use only 26 and 10 (if not using luma-prediction)
  for (i = 2; i < 35; i++) {
    int distance = MIN(abs(i-26),abs(i-10)); //!< Distance from top and left predictions
    if(distance > threshold) {
      intra_get_angular_pred(rec_filtered, recstride, pred, width, width, i, filter);
      ratecost = intra_pred_ratecost(i,intra_preds)*(int)(g_cur_lambda_cost+0.5);
      CHECK_FOR_BEST(i,ratecost);
    }
  }

  // assign final sad to output
  *sad_out = best_sad;
  #undef COPY_PRED_TO_DST
  #undef CHECK_FOR_BEST

  return best_mode;
}

/**
 * \brief Reconstruct intra block according to prediction
 * \param rec reconstructed picture data
 * \param recstride reconstructed picture stride
 * \param width block size to predict
 * \param dst destination buffer for best prediction
 * \param dststride destination width
 * \param mode intra mode to use
 * \param chroma chroma-block flag

*/
void intra_recon(pixel* rec, uint32_t recstride, uint32_t width, pixel* dst, int32_t dststride, int8_t mode, int8_t chroma)
{
  int32_t x,y;
  pixel pred[LCU_WIDTH * LCU_WIDTH];
  int8_t filter = !chroma && width < 32;


  // Filtering apply if luma and not DC
  if (!chroma && mode != 1 && width > 4) {
    uint8_t threshold = intra_hor_ver_dist_thres[g_to_bits[width]];
    if(MIN(abs(mode-26),abs(mode-10)) > threshold) {
      intra_filter(rec,recstride,width,0);
    }
  }

  // planar
  if (mode == 0)  {
    intra_get_planar_pred(rec, recstride, width, pred, width);
  } else if (mode == 1) { // DC
    pixel val = intra_get_dc_pred(rec, (uint16_t)recstride, (uint8_t)width);
    for (y = 0; y < (int32_t)width; y++) {
      for (x = 0; x < (int32_t)width; x++) {
        dst[x + y*dststride] = val;
      }
    }
    // Assigned value directly to output, no need to stay here
    return;
  } else {  // directional predictions
    intra_get_angular_pred(rec, recstride,pred, width, width, mode, filter);
  }

  for(y = 0; y < (int32_t)width; y++)  {
    for(x = 0; x < (int32_t)width; x++) {
      dst[x+y*dststride] = pred[x+y*width];
    }
  }

}

/**
 * \brief Build top and left borders for a reference block.
 * \param pic picture to use as a source
 * \param outwidth width of the prediction block
 * \param chroma signaling if chroma is used, 0 = luma, 1 = U and 2 = V
 *
 * The end result is 2*width+8 x 2*width+8 array, with only the top and left
 * edge pixels filled with the reconstructed pixels.
 */
void intra_build_reference_border(int32_t x_luma, int32_t y_luma, int16_t out_width,
                                      pixel *dst, int32_t dst_stride, int8_t chroma,
                                      int32_t pic_width, int32_t pic_height,
                                      lcu_t *lcu)
{
  // Some other function might make use of the arrays num_ref_pixels_top and
  // num_ref_pixels_left in the future, but until that happens lets leave
  // them here.

  /**
   * \brief Table for looking up the number of intra reference pixels based on
   *        prediction units coordinate within an LCU.
   *
   * This table was generated by "tools/generate_ref_pixel_tables.py".
   */
  static const uint8_t num_ref_pixels_top[16][16] = {
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 32, 28, 24, 20, 16, 12,  8,  4, 32, 28, 24, 20, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 64, 60, 56, 52, 48, 44, 40, 36, 32, 28, 24, 20, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 32, 28, 24, 20, 16, 12,  8,  4, 32, 28, 24, 20, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4, 16, 12,  8,  4 },
    {  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 }
  };

  /**
   * \brief Table for looking up the number of intra reference pixels based on
   *        prediction units coordinate within an LCU.
   *
   * This table was generated by "tools/generate_ref_pixel_tables.py".
   */
  static const uint8_t num_ref_pixels_left[16][16] = {
    { 64,  4,  8,  4, 16,  4,  8,  4, 32,  4,  8,  4, 16,  4,  8,  4 },
    { 64,  4,  4,  4, 12,  4,  4,  4, 28,  4,  4,  4, 12,  4,  4,  4 },
    { 64,  4,  8,  4,  8,  4,  8,  4, 24,  4,  8,  4,  8,  4,  8,  4 },
    { 64,  4,  4,  4,  4,  4,  4,  4, 20,  4,  4,  4,  4,  4,  4,  4 },
    { 64,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4 },
    { 64,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4 },
    { 64,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 64,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4 },
    { 64,  4,  8,  4, 16,  4,  8,  4, 32,  4,  8,  4, 16,  4,  8,  4 },
    { 64,  4,  4,  4, 12,  4,  4,  4, 28,  4,  4,  4, 12,  4,  4,  4 },
    { 64,  4,  8,  4,  8,  4,  8,  4, 24,  4,  8,  4,  8,  4,  8,  4 },
    { 64,  4,  4,  4,  4,  4,  4,  4, 20,  4,  4,  4,  4,  4,  4,  4 },
    { 64,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4, 16,  4,  8,  4 },
    { 64,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4, 12,  4,  4,  4 },
    { 64,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4,  8,  4 },
    { 64,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4 }
  };

  const pixel dc_val = 1 << (g_bitdepth - 1);
  const int is_chroma = chroma ? 1 : 0;

  // input picture pointer
  //const pixel * const src = (!chroma) ? pic->y_recdata : ((chroma == 1) ? pic->u_recdata : pic->v_recdata);

  // Convert luma coordinates to chroma coordinates for chroma.
  const int x = chroma ? x_luma / 2 : x_luma;
  const int y = chroma ? y_luma / 2 : y_luma;

  const int y_in_lcu = y_luma % LCU_WIDTH;
  const int x_in_lcu = x_luma % LCU_WIDTH;

  int x_local = (x_luma&0x3f)>>is_chroma, y_local = (y_luma&0x3f)>>is_chroma;

  pixel *left_ref = !chroma ? &lcu->left_ref.y[1] : (chroma == 1) ? &lcu->left_ref.u[1] : &lcu->left_ref.v[1];
  pixel *top_ref  = !chroma ? &lcu->top_ref.y[1]  : (chroma == 1) ? &lcu->top_ref.u[1]  : &lcu->top_ref.v[1];
  pixel *rec_ref  = !chroma ? lcu->rec.y : (chroma == 1) ? lcu->rec.u : lcu->rec.v;

  pixel *left_border = &left_ref[y_local];
  pixel *top_border = &top_ref[x_local];
  uint32_t left_stride = 1;

  if(x_local) {
    left_border = &rec_ref[x_local - 1 + y_local * (LCU_WIDTH>>is_chroma)];
    left_stride = LCU_WIDTH>>is_chroma;
  }

  if(y_local) {
    top_border = &rec_ref[x_local + (y_local - 1) * (LCU_WIDTH>>is_chroma)];
  }

  // Copy pixels for left edge.
  if (x > 0) {
    // Get the number of reference pixels based on the PU coordinate within the LCU.
    int num_ref_pixels = num_ref_pixels_left[y_in_lcu / 4][x_in_lcu / 4] >> is_chroma;
    int i;
    pixel nearest_pixel;

    // Max pixel we can copy from src is yy + outwidth - 1 because the dst
    // extends one pixel to the left.
    num_ref_pixels = MIN(num_ref_pixels, out_width - 1);
    // There are no coded pixels below the frame.
    num_ref_pixels = MIN(num_ref_pixels, pic_height - y);
    // There are no coded pixels below the bottom of the LCU due to raster
    // scan order.
    num_ref_pixels = MIN(num_ref_pixels, (LCU_WIDTH - y_in_lcu) >> is_chroma);

    // Copy pixels from coded CUs.
    for (i = 0; i < num_ref_pixels; ++i) {
      dst[(i + 1) * dst_stride] = left_border[i*left_stride];
    }
    // Extend the last pixel for the rest of the reference values.
    nearest_pixel = dst[i * dst_stride];
    for (i = num_ref_pixels; i < out_width - 1; ++i) {
      dst[i * dst_stride] = nearest_pixel;
    }
  } else {
    // If we are on the left edge, extend the first pixel of the top row.
    pixel nearest_pixel = y > 0 ? top_border[0] : dc_val;
    int i;
    for (i = 1; i < out_width - 1; i++) {
      dst[i * dst_stride] = nearest_pixel;
    }
  }

  // Copy pixels for top edge.
  if (y > 0) {
    // Get the number of reference pixels based on the PU coordinate within the LCU.
    int num_ref_pixels = num_ref_pixels_top[y_in_lcu / 4][x_in_lcu / 4] >> is_chroma;
    int i;
    pixel nearest_pixel;

    // Max pixel we can copy from src is yy + outwidth - 1 because the dst
    // extends one pixel to the left.
    num_ref_pixels = MIN(num_ref_pixels, out_width - 1);
    // All LCUs in the row above have been coded.
    num_ref_pixels = MIN(num_ref_pixels, pic_width - x);

    // Copy pixels from coded CUs.
    for (i = 0; i < num_ref_pixels; ++i) {
      dst[i + 1] = top_border[i];
    }
    // Extend the last pixel for the rest of the reference values.
    nearest_pixel = top_border[num_ref_pixels - 1];
    for (; i < out_width - 1; ++i) {
      dst[i + 1] = nearest_pixel;
    }
  } else {
    // Extend nearest pixel.
    pixel nearest_pixel = x > 0 ? left_border[0] : dc_val;
    int i;
    for(i = 1; i < out_width; i++)
    {
      dst[i] = nearest_pixel;
    }
  }

  // If top-left corner sample doesn't exist, use the sample from below.
  // Unavailable samples on the left boundary are copied from below if
  // available. This is the only place they are available because we don't
  // support constrained intra prediction.
  if (x > 0 && y > 0) {
    // Make sure we always take the top-left pixel from the LCU reference
    // pixel arrays if they are available.
    if (x_local == 0) {
      dst[0] = left_border[-1];
    } else {
      dst[0] = top_border[-1];
    }
  } else {
    dst[0] = dst[dst_stride];
  }
}

const int32_t ang_table[9]     = {0,    2,    5,   9,  13,  17,  21,  26,  32};
const int32_t inv_ang_table[9] = {0, 4096, 1638, 910, 630, 482, 390, 315, 256}; // (256 * 32) / Angle

/**
 * \brief this functions constructs the angular intra prediction from border samples
 *
 */
void intra_get_angular_pred(pixel* src, int32_t src_stride, pixel* dst, int32_t dst_stride, int32_t width, int32_t dir_mode, int8_t filter)
{
  int32_t k,l;
  int32_t blk_size        = width;

  // Map the mode index to main prediction direction and angle
  int8_t mode_hor       = dir_mode < 18;
  int8_t mode_ver       = !mode_hor;
  int32_t intra_pred_angle = mode_ver ? (int32_t)dir_mode - 26 : mode_hor ? -((int32_t)dir_mode - 10) : 0;
  int32_t abs_ang       = abs(intra_pred_angle);
  int32_t sign_ang      = intra_pred_angle < 0 ? -1 : 1;

  // Set bitshifts and scale the angle parameter to block size
  int32_t inv_angle       = inv_ang_table[abs_ang];

  // Do angular predictions
  pixel *ref_main;
  pixel *ref_side;
  pixel  ref_above[2 * LCU_WIDTH + 1];
  pixel  ref_left[2 * LCU_WIDTH + 1];

  abs_ang           = ang_table[abs_ang];
  intra_pred_angle  = sign_ang * abs_ang;

  // Initialise the Main and Left reference array.
  if (intra_pred_angle < 0) {
    int32_t invAngleSum = 128; // rounding for (shift by 8)
    for (k = 0; k < blk_size + 1; k++) {
      ref_above[k + blk_size - 1] = src[k - src_stride - 1];
      ref_left[k + blk_size - 1]  = src[(k - 1) * src_stride - 1];
    }

    ref_main = (mode_ver ? ref_above : ref_left) + (blk_size - 1);
    ref_side = (mode_ver ? ref_left : ref_above) + (blk_size - 1);

    // Extend the Main reference to the left.
    for (k = -1; k > blk_size * intra_pred_angle>>5; k--) {
      invAngleSum += inv_angle;
      ref_main[k] = ref_side[invAngleSum>>8];
    }
  } else {
    for (k = 0; k < 2 * blk_size + 1; k++) {
      ref_above[k] = src[k - src_stride - 1];
      ref_left[k]  = src[(k - 1) * src_stride - 1];
    }
    ref_main = mode_ver ? ref_above : ref_left;
    ref_side = mode_ver ? ref_left  : ref_above;
  }

  if (intra_pred_angle == 0) {
    for (k = 0; k < blk_size; k++) {
      for (l = 0; l < blk_size; l++) {
        dst[k * dst_stride + l] = ref_main[l + 1];
      }
    }

    if (filter) {
      for (k=0;k<blk_size;k++) {
        dst[k * dst_stride] = CLIP(0, (1<<g_bitdepth) - 1, dst[k * dst_stride] + (( ref_side[k + 1] - ref_side[0]) >> 1));
      }
    }
  } else {
    int32_t delta_pos=0;
    int32_t delta_int;
    int32_t delta_fract;
    int32_t minus_delta_fract;
    int32_t ref_main_index;
    for (k = 0; k < blk_size; k++) {
      delta_pos += intra_pred_angle;
      delta_int   = delta_pos >> 5;
      delta_fract = delta_pos & (32 - 1);


      if (delta_fract) {
        minus_delta_fract = (32 - delta_fract);
        // Do linear filtering
        for (l = 0; l < blk_size; l++) {
          ref_main_index        = l + delta_int + 1;
          dst[k * dst_stride + l] = (pixel) ( (minus_delta_fract * ref_main[ref_main_index]
                                                 + delta_fract * ref_main[ref_main_index + 1] + 16) >> 5);
        }
      } else {
        // Just copy the integer samples
        for (l = 0; l < blk_size; l++) {
          dst[k * dst_stride + l] = ref_main[l + delta_int + 1];
        }
      }
    }
  }

  // Flip the block if this is the horizontal mode
  if (mode_hor) {
    pixel tmp;
    for (k=0;k<blk_size-1;k++) {
      for (l=k+1;l<blk_size;l++) {
        tmp                 = dst[k * dst_stride + l];
        dst[k * dst_stride + l] = dst[l * dst_stride + k];
        dst[l * dst_stride + k] = tmp;
      }
    }
  }

}




void intra_dc_pred_filtering(pixel *src, int32_t src_stride, pixel *dst, int32_t dst_stride, int32_t width, int32_t height )
{
  int32_t x, y, dst_stride2, src_stride2;

  // boundary pixels processing
  dst[0] = ((src[-src_stride] + src[-1] + 2 * dst[0] + 2) >> 2);

  for (x = 1; x < width; x++) {
    dst[x] = ((src[x - src_stride] +  3 * dst[x] + 2) >> 2);
  }
  for ( y = 1, dst_stride2 = dst_stride, src_stride2 = src_stride-1;
        y < height; y++, dst_stride2+=dst_stride, src_stride2+=src_stride ) {
    dst[dst_stride2] = ((src[src_stride2] + 3 * dst[dst_stride2] + 2) >> 2);
  }
  return;
}

/**
 * \brief Function for deriving planar intra prediction.
 * \param src source pixel array
 * \param srcstride source width
 * \param width block size to predict
 * \param dst destination buffer for prediction
 * \param dststride destination width

  This function derives the prediction samples for planar mode (intra coding).
*/
void intra_get_planar_pred(pixel* src, int32_t srcstride, uint32_t width, pixel* dst, int32_t dststride)
{
  int32_t k, l, bottom_left, top_right;
  int32_t hor_pred;
  int32_t left_column[LCU_WIDTH+1], top_row[LCU_WIDTH+1], bottom_row[LCU_WIDTH+1], right_column[LCU_WIDTH+1];
  uint32_t blk_size = width;
  uint32_t offset_2d = width;
  uint32_t shift_1d = g_convert_to_bit[ width ] + 2;
  uint32_t shift_2d = shift_1d + 1;

  // Get left and above reference column and row
  for (k = 0; k < (int32_t)blk_size + 1; k++) {
    top_row[k] = src[k - srcstride];
    left_column[k] = src[k * srcstride - 1];
  }

  // Prepare intermediate variables used in interpolation
  bottom_left = left_column[blk_size];
  top_right   = top_row[blk_size];
  for (k = 0; k < (int32_t)blk_size; k++) {
    bottom_row[k]   = bottom_left - top_row[k];
    right_column[k] = top_right   - left_column[k];
    top_row[k]      <<= shift_1d;
    left_column[k]  <<= shift_1d;
  }

  // Generate prediction signal
  for (k = 0; k < (int32_t)blk_size; k++) {
    hor_pred = left_column[k] + offset_2d;
    for (l = 0; l < (int32_t)blk_size; l++) {
      hor_pred += right_column[k];
      top_row[l] += bottom_row[l];
      dst[k * dststride + l] = (pixel)((hor_pred + top_row[l]) >> shift_2d);
    }
  }
}

void intra_recon_lcu(encoder_control* encoder, int x, int y, int depth, lcu_t *lcu, uint32_t pic_width, uint32_t pic_height)
{
  int x_local = (x&0x3f), y_local = (y&0x3f);
  cu_info *cur_cu = &lcu->cu[LCU_CU_OFFSET + (x_local>>3) + (y_local>>3)*LCU_T_CU_WIDTH];
  
  // Pointers to reconstruction arrays
  pixel *recbase_y = &lcu->rec.y[x_local + y_local * LCU_WIDTH];
  pixel *recbase_u = &lcu->rec.u[x_local/2 + (y_local * LCU_WIDTH)/4];
  pixel *recbase_v = &lcu->rec.v[x_local/2 + (y_local * LCU_WIDTH)/4];
  int32_t rec_stride = LCU_WIDTH;


  pixel rec[(LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)];
  pixel *rec_shift  = &rec[(LCU_WIDTH >> (depth)) * 2 + 8 + 1];

  int8_t width = LCU_WIDTH >> depth;
  int8_t width_c = LCU_WIDTH >> (depth + 1);
  static vector2d offsets[4] = {{0,0},{1,0},{0,1},{1,1}};
  int num_pu = (cur_cu->part_size == SIZE_2Nx2N ? 1 : 4);
  int i;

  if (cur_cu->part_size == SIZE_NxN) {
    width = width_c;
  }

  cur_cu->intra[0].mode_chroma = 36; // TODO: Chroma intra prediction

  // Reconstruct chroma
  rec_shift  = &rec[width_c * 2 + 8 + 1];
  intra_build_reference_border(x, y,(int16_t)width_c * 2 + 8, rec, (int16_t)width_c * 2 + 8, 1,
                                   pic_width/2, pic_height/2, lcu);
  intra_recon(rec_shift,
              width_c * 2 + 8,
              width_c,
              recbase_u,
              rec_stride >> 1,
              cur_cu->intra[0].mode_chroma != 36 ? cur_cu->intra[0].mode_chroma : cur_cu->intra[0].mode,
              1);

  intra_build_reference_border(x, y,(int16_t)width_c * 2 + 8, rec, (int16_t)width_c * 2 + 8, 2,
                                   pic_width/2, pic_height/2, lcu);
  intra_recon(rec_shift,
              width_c * 2 + 8,
              width_c,
              recbase_v,
              rec_stride >> 1,
              cur_cu->intra[0].mode_chroma != 36 ? cur_cu->intra[0].mode_chroma : cur_cu->intra[0].mode,
              2);

  for (i = 0; i < num_pu; ++i) {
    // Build reconstructed block to use in prediction with extrapolated borders
    int x_off = offsets[i].x * width;
    int y_off = offsets[i].y * width;
    recbase_y = &lcu->rec.y[x_local + x_off + (y_local+y_off) * LCU_WIDTH];

    rec_shift  = &rec[width * 2 + 8 + 1];
    intra_build_reference_border(x+x_off, y+y_off,(int16_t)width * 2 + 8, rec, (int16_t)width * 2 + 8, 0,
                                 pic_width, pic_height, lcu);
    intra_recon(rec_shift, width * 2 + 8,
                width, recbase_y, rec_stride, cur_cu->intra[i].mode, 0);

    // Filter DC-prediction
    if (cur_cu->intra[i].mode == 1 && width < 32) {
      intra_dc_pred_filtering(rec_shift, width * 2 + 8, recbase_y,
                              rec_stride, width, width);
    }

    // Handle NxN mode by doing quant/transform and inverses for the next NxN block
    if (cur_cu->part_size == SIZE_NxN) {
      encode_transform_tree(encoder, x + x_off, y + y_off, depth+1, lcu);
    }
  }

  // If we coded NxN block, fetch the coded block flags to this level
  if (cur_cu->part_size == SIZE_NxN) {
    cur_cu->coeff_top_y[depth] = cur_cu->coeff_top_y[depth+1] | cur_cu->coeff_top_y[depth+2] | cur_cu->coeff_top_y[depth+3] | cur_cu->coeff_top_y[depth+4];
    cur_cu->coeff_top_u[depth] = cur_cu->coeff_top_u[depth+1];
    cur_cu->coeff_top_v[depth] = cur_cu->coeff_top_v[depth+1];
    return;
  }

  encode_transform_tree(encoder, x, y, depth, lcu);
}