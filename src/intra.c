/**
 * \file
 * \brief Functions for handling intra frames.
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
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
void intra_set_block_mode(picture *pic,uint32_t x_cu, uint32_t y_cu, uint8_t depth, uint8_t mode)
{
  uint32_t x, y;  
  int width_in_scu = pic->width_in_lcu<<MAX_DEPTH; //!< Width in smallest CU
  int block_scu_width = (LCU_WIDTH>>depth)/(LCU_WIDTH>>MAX_DEPTH);

  // Loop through all the blocks in the area of cur_cu
  for (y = y_cu; y < y_cu + block_scu_width; y++) {
    int cu_pos = y * width_in_scu;
    for (x = x_cu; x < x_cu + block_scu_width; x++) {
      pic->cu_array[MAX_DEPTH][cu_pos + x].depth = depth;
      pic->cu_array[MAX_DEPTH][cu_pos + x].type  = CU_INTRA;
      pic->cu_array[MAX_DEPTH][cu_pos + x].intra.mode = mode;
    }
  }
}

/**
 * \brief get intrablock mode
 * \param pic picture to use
 * \param xCtb x CU position (smallest CU)
 * \param yCtb y CU position (smallest CU)
 * \param depth current CU depth
 * \returns mode if it's present, otherwise -1
*/
int8_t intra_get_block_mode(picture *pic, uint32_t x_cu, uint32_t y_cu, uint8_t depth)
{ 
  int width_in_scu = pic->width_in_lcu<<MAX_DEPTH; //!< width in smallest CU
  int cu_pos = y_cu * width_in_scu + x_cu;
  if (pic->cu_array[MAX_DEPTH][cu_pos].type == CU_INTRA) {
    return pic->cu_array[MAX_DEPTH][cu_pos].intra.mode;
  }
  return -1;
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
int16_t intra_get_dc_pred(int16_t *pic, uint16_t picwidth, uint32_t xpos, uint32_t ypos, uint8_t width)
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
  return (sum + width) / (width + width);
}

/** 
 * \brief Function for deriving intra luma predictions
 * \param pic picture to use
 * \param x_cu x CU position (smallest CU)
 * \param y_cu y CU position (smallest CU)
 * \param depth current CU depth
 * \param preds output buffer for 3 predictions 
 * \returns (predictions are found)?1:0
 */
int8_t intra_get_dir_luma_predictor(picture* pic, uint32_t x_cu, uint32_t y_cu, uint8_t depth, int8_t* preds)
{
  int32_t left_intra_dir  = 1; // reset to DC_IDX
  int32_t above_intra_dir = 1; // reset to DC_IDX
  int width_in_scu = pic->width_in_lcu<<MAX_DEPTH;
  int32_t cu_pos = y_cu * width_in_scu + x_cu;
  
  // Left PU predictor
  if(x_cu && pic->cu_array[MAX_DEPTH][cu_pos - 1].type == CU_INTRA && pic->cu_array[MAX_DEPTH][cu_pos - 1].coded) {
    left_intra_dir = pic->cu_array[MAX_DEPTH][cu_pos - 1].intra.mode;
  }

  // Top PU predictor
  if(y_cu && ((y_cu * (LCU_WIDTH>>MAX_DEPTH)) % LCU_WIDTH) != 0
     && pic->cu_array[MAX_DEPTH][cu_pos - width_in_scu].type == CU_INTRA && pic->cu_array[MAX_DEPTH][cu_pos - width_in_scu].coded) {
    above_intra_dir = pic->cu_array[MAX_DEPTH][cu_pos - width_in_scu].intra.mode;
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
void intra_filter(int16_t *ref, int32_t stride,int32_t width, int8_t mode)
{
  #define FWIDTH (LCU_WIDTH*2+1)
  int16_t filtered[FWIDTH * FWIDTH]; //!< temporary buffer for filtered samples
  int16_t *filteredShift = &filtered[FWIDTH+1]; //!< pointer to temporary buffer with offset (1,1)
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
int16_t intra_prediction(pixel *orig, int32_t origstride, int16_t *rec, int32_t recstride, uint32_t xpos,
                         uint32_t ypos, uint32_t width, int16_t *dst, int32_t dststride, uint32_t *sad_out)
{
  typedef uint32_t (*sad_function)(int16_t *block,uint32_t stride1,int16_t *block2, uint32_t stride2);
  uint32_t best_sad = 0xffffffff;
  uint32_t sad = 0;
  int16_t best_mode = 1;
  int32_t x,y,i;
  sad_function calc_sad;

  // Temporary block arrays
  // TODO: alloc with alignment
  int16_t pred[LCU_WIDTH * LCU_WIDTH + 1];  
  int16_t orig_block[LCU_WIDTH * LCU_WIDTH + 1];  
  int16_t rec_filtered_temp[(LCU_WIDTH * 2 + 8) * (LCU_WIDTH * 2 + 8) + 1];
  
  int16_t* rec_filtered = &rec_filtered_temp[recstride + 1]; //!< pointer to rec_filtered_temp with offset of (1,1)
  pixel *orig_shift = &orig[xpos + ypos*origstride];  //!< pointer to orig with offset of (1,1)
  int8_t filter = (width<32); // TODO: chroma support

  sad_function sad_array[5] = {&sad4x4,&sad8x8,&sad16x16,&sad32x32,&sad64x64}; //TODO: get SAD functions from parameters
  uint8_t threshold = intra_hor_ver_dist_thres[g_to_bits[width]]; //!< Intra filtering threshold

  #define COPY_PRED_TO_DST() for (y = 0; y < (int32_t)width; y++)  { for (x = 0; x < (int32_t)width; x++) { dst[x + y*dststride] = pred[x + y*width]; } }
  #define CHECK_FOR_BEST(mode, additional_sad)  sad = calc_sad(pred,width,orig_block,width); \
                                                sad += additional_sad;\
                                                if(sad < best_sad)\
                                                {\
                                                  best_sad = sad;\
                                                  best_mode = mode;\
                                                  COPY_PRED_TO_DST();\
                                                }

  // Choose SAD function according to width
  calc_sad = sad_array[g_to_bits[width]];

  // Store original block for SAD computation
  i = 0;
  for(y = 0; y < (int32_t)width; y++) {
    for(x = 0; x < (int32_t)width; x++) {
      orig_block[i++] = orig_shift[x + y*origstride];
    }
  }

  // Filtered only needs the borders
  for (y = -1; y < (int32_t)recstride; y++) {
    rec_filtered[y*recstride - 1] = rec[y*recstride - 1];
  }
  for (x = 0; x < (int32_t)recstride; x++) {
    rec_filtered[y - recstride] = rec[y - recstride];
  }    
  // Apply filter
  intra_filter(rec_filtered,recstride,width,0);
  

  // Test DC mode (never filtered)
  x = intra_get_dc_pred(rec, recstride, xpos, ypos, width);
  for (i = 0; i < (int32_t)(width*width); i++) {
    pred[i] = x;
  }
  CHECK_FOR_BEST(1,0);
  
  // Check angular not requiring filtering
  for (i = 2; i < 35; i++) {
    int distance = MIN(abs(i - 26),abs(i - 10)); //!< Distance from top and left predictions
    if(distance <= threshold) {
      intra_get_angular_pred(rec, recstride, pred, width, width, width, i, xpos?1:0, ypos?1:0, filter);
      CHECK_FOR_BEST(i,0);
    }
  }
  
  // FROM THIS POINT FORWARD, USING FILTERED PREDICTION

  // Test planar mode (always filtered)
  intra_get_planar_pred(rec_filtered, recstride, xpos, ypos, width, pred, width);
  CHECK_FOR_BEST(0,0);  
  
  // Check angular predictions which require filtered samples
  // TODO: add conditions to skip some modes on borders  
  // chroma can use only 26 and 10 (if not using luma-prediction)  
  for (i = 2; i < 35; i++) {
    int distance = MIN(abs(i-26),abs(i-10)); //!< Distance from top and left predictions
    if(distance > threshold) {
      intra_get_angular_pred(rec_filtered, recstride, pred, width, width, width, i, xpos?1:0, ypos?1:0, filter);
      CHECK_FOR_BEST(i,0);
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
 * \param xpos source x-position
 * \param ypos source y-position
 * \param width block size to predict
 * \param dst destination buffer for best prediction
 * \param dststride destination width
 * \param mode intra mode to use
 * \param chroma chroma-block flag

*/
void intra_recon(int16_t* rec,uint32_t recstride, uint32_t xpos, uint32_t ypos,uint32_t width, int16_t* dst,int32_t dststride, int8_t mode, int8_t chroma)
{
  int32_t x,y,i;
  int16_t pred[LCU_WIDTH * LCU_WIDTH];
  int8_t filter = !chroma&&(width<32);
  #define COPY_PRED_TO_DST() for(y = 0; y < (int32_t)width; y++)  { for(x = 0; x < (int32_t)width; x++) { dst[x+y*dststride] = pred[x+y*width]; } }

  // Filtering apply if luma and not DC
  if (!chroma && mode != 1) {
    uint8_t threshold = intra_hor_ver_dist_thres[g_to_bits[width]];
    if(MIN(abs(mode-26),abs(mode-10)) > threshold) {
      intra_filter(rec,recstride,width,0);
    }
  }

  // planar
  if (mode == 0)  {
    intra_get_planar_pred(rec, recstride, xpos, ypos, width, pred, width); 
  } else if (mode == 1) { // DC
    i = intra_get_dc_pred(rec, recstride, xpos, ypos, width);
    for (y = 0; y < (int32_t)width; y++) {
      for (x = 0; x < (int32_t)width; x++) {
        dst[x + y*dststride] = i;
      }
    }
    // Assigned value directly to output, no need to stay here
    return;
  } else {  // directional predictions
    intra_get_angular_pred(rec, recstride,pred, width, width, width, mode, xpos?1:0, ypos?1:0, filter);
  }

  COPY_PRED_TO_DST();

  #undef COPY_PRED_TO_DST
}


/** 
 * \brief this functions build a reference block (only borders) used for intra predictions
 * \param pic picture to use as a source, should contain full CU-data
 * \param outwidth width of the prediction block
 * \param chroma signaling if chroma is used, 0 = luma, 1 = U and 2 = V    
 *
 */
void intra_build_reference_border(picture *pic, int32_t x_cu, int32_t y_cu,int16_t outwidth, int16_t *dst, int32_t dststride, int8_t chroma)
{
  int32_t left_column; //!< left column iterator
  int16_t val;         //!< variable to store extrapolated value
  int32_t i;           //!< index iterator
  int16_t dc_val       = 1<<(g_bitdepth-1); //!< default predictor value
  int32_t top_row;     //!< top row iterator
  int32_t src_width    = (pic->width>>(chroma?1:0)); //!< source picture width
  int32_t src_height   = (pic->height>>(chroma?1:0));//!< source picture height
  pixel *src         = (!chroma) ? pic->y_recdata : ((chroma == 1) ? pic->u_recdata : pic->v_recdata); //!< input picture pointer
  int16_t scu_width    = LCU_WIDTH>>(MAX_DEPTH+(chroma?1:0)); //!< Smallest Coding Unit width
  pixel *src_shifted = &src[x_cu * scu_width + (y_cu * scu_width) * src_width];  //!< input picture pointer shifted to start from the left-top corner of the current block
  int width_in_scu = pic->width_in_lcu<<MAX_DEPTH;     //!< picture width in smallest CU

  // Fill left column when not on the border
  if (x_cu) {
    // loop SCU's
    for (left_column = 1; left_column < outwidth / scu_width; left_column++) {
      // If over the picture height or block not yet coded, stop
      if ((y_cu + left_column) * scu_width >= src_height || !pic->cu_array[MAX_DEPTH][x_cu - 1 + (y_cu + left_column) * width_in_scu].coded) {
        break;
      }
    }
    // Copy the pixels to output
    for (i = 0; i < left_column*scu_width - 1; i ++) {
      dst[(i + 1) * dststride] = src_shifted[i*src_width - 1];
    }

    // if the loop was not completed, extrapolate the last pixel pushed to output
    if (left_column != outwidth / scu_width) {
      val = src_shifted[(left_column * scu_width - 1) * src_width - 1];
      for(i = (left_column * scu_width); i < outwidth; i++) {
        dst[i * dststride] = val;
      }
    }    
  } else { // If left column not available, copy from toprow or use the default predictor
    val = y_cu ? src_shifted[-src_width] : dc_val;
    for (i = 0; i < outwidth; i++) {
      dst[i * dststride] = val;
    }
  }

  if(y_cu) {
    // Loop top SCU's
    for(top_row = 1; top_row < outwidth / scu_width; top_row++)  {
      // If over the picture width or block not yet coded, stop
      if ((x_cu + top_row) * scu_width >= src_width || !pic->cu_array[MAX_DEPTH][x_cu + top_row+(y_cu - 1) * width_in_scu].coded) {
        break;
      }
    }

    // Copy the pixels to output
    for(i = 0; i < top_row * scu_width - 1; i++) {
      dst[i + 1] = src_shifted[i - src_width];
    }

    if(top_row != outwidth/scu_width) {
      val = src_shifted[(top_row * scu_width) - src_width - 1];
      for(i = (top_row * scu_width); i < outwidth; i++) {
        dst[i] = val;
      }
    }
  } else {
    val = x_cu ? src_shifted[-1] : dc_val;
    for(i = 1; i < outwidth; i++)
    {
      dst[i] = val;
    }
  }
  // Topleft corner sample
  dst[0] = (x_cu && y_cu) ? src_shifted[-src_width - 1] : dst[dststride];

}

const int32_t ang_table[9]     = {0,    2,    5,   9,  13,  17,  21,  26,  32};
const int32_t inv_ang_table[9] = {0, 4096, 1638, 910, 630, 482, 390, 315, 256}; // (256 * 32) / Angle

/** 
 * \brief this functions constructs the angular intra prediction from border samples
 *
 */
void intra_get_angular_pred(int16_t* src, int32_t src_stride, int16_t* dst, int32_t dst_stride, int32_t width,
                           int32_t height, int32_t dir_mode, int8_t left_avail,int8_t top_avail, int8_t filter)
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
  int16_t *ref_main;
  int16_t *ref_side;
  int16_t  ref_above[2 * LCU_WIDTH + 1];
  int16_t  ref_left[2 * LCU_WIDTH + 1];

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
    for (k =- 1; k > blk_size * intra_pred_angle>>5; k--) {
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
          dst[k * dst_stride + l] = (int16_t) ( (minus_delta_fract * ref_main[ref_main_index]
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
    int16_t tmp;
    for (k=0;k<blk_size-1;k++) {
      for (l=k+1;l<blk_size;l++) {
        tmp                 = dst[k * dst_stride + l];
        dst[k * dst_stride + l] = dst[l * dst_stride + k];
        dst[l * dst_stride + k] = tmp;
      }
    }
  }

}




void intra_dc_pred_filtering(int16_t *src, int32_t src_stride, int16_t *dst, int32_t dst_stride, int32_t width, int32_t height )
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
 * \param xpos source x-position
 * \param ypos source y-position
 * \param width block size to predict
 * \param dst destination buffer for prediction
 * \param dststride destination width
 
  This function derives the prediction samples for planar mode (intra coding).
*/
void intra_get_planar_pred(int16_t* src,int32_t srcstride, uint32_t xpos, uint32_t ypos,uint32_t width, int16_t* dst,int32_t dststride)
{
  int16_t dc_val = 1<<(g_bitdepth-1);
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
      dst[k * dststride + l] = ( (hor_pred + top_row[l]) >> shift_2d );
    }
  }
}
