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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rdo.h"
#include "transform.h"
#include "context.h"
#include "cabac.h"
#include "transform.h"


#define QUANT_SHIFT          14
#define SCAN_SET_SIZE        16
#define LOG2_SCAN_SET_SIZE    4
#define SBH_THRESHOLD         4

const uint32_t g_go_rice_range[5] = { 7, 14, 26, 46, 78 };
const uint32_t g_go_rice_prefix_len[5] = { 8, 7, 6, 5, 4 };


#define CTX_ENTROPY_BITS(ctx,val) entropy_bits[(ctx)->uc_state ^ val]
/**
 * Entropy bits to estimate coded bits in RDO / RDOQ (From HM 12.0)
 */
const uint32_t entropy_bits[128] =
{
  0x08000, 0x08000, 0x076da, 0x089a0, 0x06e92, 0x09340, 0x0670a, 0x09cdf, 0x06029, 0x0a67f, 0x059dd, 0x0b01f, 0x05413, 0x0b9bf, 0x04ebf, 0x0c35f,
  0x049d3, 0x0ccff, 0x04546, 0x0d69e, 0x0410d, 0x0e03e, 0x03d22, 0x0e9de, 0x0397d, 0x0f37e, 0x03619, 0x0fd1e, 0x032ee, 0x106be, 0x02ffa, 0x1105d,
  0x02d37, 0x119fd, 0x02aa2, 0x1239d, 0x02836, 0x12d3d, 0x025f2, 0x136dd, 0x023d1, 0x1407c, 0x021d2, 0x14a1c, 0x01ff2, 0x153bc, 0x01e2f, 0x15d5c,
  0x01c87, 0x166fc, 0x01af7, 0x1709b, 0x0197f, 0x17a3b, 0x0181d, 0x183db, 0x016d0, 0x18d7b, 0x01595, 0x1971b, 0x0146c, 0x1a0bb, 0x01354, 0x1aa5a,
  0x0124c, 0x1b3fa, 0x01153, 0x1bd9a, 0x01067, 0x1c73a, 0x00f89, 0x1d0da, 0x00eb7, 0x1da79, 0x00df0, 0x1e419, 0x00d34, 0x1edb9, 0x00c82, 0x1f759,
  0x00bda, 0x200f9, 0x00b3c, 0x20a99, 0x00aa5, 0x21438, 0x00a17, 0x21dd8, 0x00990, 0x22778, 0x00911, 0x23118, 0x00898, 0x23ab8, 0x00826, 0x24458,
  0x007ba, 0x24df7, 0x00753, 0x25797, 0x006f2, 0x26137, 0x00696, 0x26ad7, 0x0063f, 0x27477, 0x005ed, 0x27e17, 0x0059f, 0x287b6, 0x00554, 0x29156,
  0x0050e, 0x29af6, 0x004cc, 0x2a497, 0x0048d, 0x2ae35, 0x00451, 0x2b7d6, 0x00418, 0x2c176, 0x003e2, 0x2cb15, 0x003af, 0x2d4b5, 0x0037f, 0x2de55
};

// Entropy bits scaled so that 50% probability yields 1 bit.
const float f_entropy_bits[128] =
{
  1.0, 1.0,
  0.92852783203125, 1.0751953125,
  0.86383056640625, 1.150390625,
  0.80499267578125, 1.225555419921875,
  0.751251220703125, 1.300750732421875,
  0.702056884765625, 1.375946044921875,
  0.656829833984375, 1.451141357421875,
  0.615203857421875, 1.526336669921875,
  0.576751708984375, 1.601531982421875,
  0.54119873046875, 1.67669677734375,
  0.508209228515625, 1.75189208984375,
  0.47760009765625, 1.82708740234375,
  0.449127197265625, 1.90228271484375,
  0.422637939453125, 1.97747802734375,
  0.39788818359375, 2.05267333984375,
  0.37481689453125, 2.127838134765625,
  0.353240966796875, 2.203033447265625,
  0.33306884765625, 2.278228759765625,
  0.31414794921875, 2.353424072265625,
  0.29644775390625, 2.428619384765625,
  0.279815673828125, 2.5037841796875,
  0.26422119140625, 2.5789794921875,
  0.24957275390625, 2.6541748046875,
  0.235809326171875, 2.7293701171875,
  0.222869873046875, 2.8045654296875,
  0.210662841796875, 2.879730224609375,
  0.199188232421875, 2.954925537109375,
  0.188385009765625, 3.030120849609375,
  0.17822265625, 3.105316162109375,
  0.168609619140625, 3.180511474609375,
  0.1595458984375, 3.255706787109375,
  0.1510009765625, 3.33087158203125,
  0.1429443359375, 3.40606689453125,
  0.135345458984375, 3.48126220703125,
  0.128143310546875, 3.55645751953125,
  0.121368408203125, 3.63165283203125,
  0.114959716796875, 3.706817626953125,
  0.10888671875, 3.782012939453125,
  0.1031494140625, 3.857208251953125,
  0.09771728515625, 3.932403564453125,
  0.09259033203125, 4.007598876953125,
  0.0877685546875, 4.082794189453125,
  0.083160400390625, 4.157958984375,
  0.078826904296875, 4.233154296875,
  0.07470703125, 4.308349609375,
  0.070831298828125, 4.383544921875,
  0.067138671875, 4.458740234375,
  0.06365966796875, 4.533935546875,
  0.06036376953125, 4.609100341796875,
  0.057220458984375, 4.684295654296875,
  0.05426025390625, 4.759490966796875,
  0.05145263671875, 4.834686279296875,
  0.048797607421875, 4.909881591796875,
  0.046295166015625, 4.985076904296875,
  0.043914794921875, 5.06024169921875,
  0.0416259765625, 5.13543701171875,
  0.03948974609375, 5.21063232421875,
  0.0374755859375, 5.285858154296875,
  0.035552978515625, 5.360992431640625,
  0.033721923828125, 5.43621826171875,
  0.031982421875, 5.51141357421875,
  0.03033447265625, 5.586578369140625,
  0.028778076171875, 5.661773681640625,
  0.027313232421875, 5.736968994140625,
};


/**
 * \brief Function to compare RDO costs
 * \param rdo_costs array of current costs
 * \param cost new cost to check
 * \returns -1 if cost is worse than the one in the array or array position for worst cost

 This function derives the prediction samples for planar mode (intra coding).
*/
int intra_rdo_cost_compare(uint32_t *rdo_costs,int8_t rdo_modes_to_check, uint32_t cost)
{
  int i;
  int found = 0;

  for(i = 0; i < rdo_modes_to_check; i++) {
    if(rdo_costs[i] > cost) {
      found = 1;
      break;
    }
  }
  // Search for worst cost
  if(found) {
    uint32_t worst_cost = 0;
    int worst_mode = -1;
    for(i = 0; i < rdo_modes_to_check; i++) {
      if(rdo_costs[i] > worst_cost) {
        worst_cost = rdo_costs[i];
        worst_mode = i;
      }
    }
    return worst_mode;
  }

  return -1;
}


/**
 * \brief RDO function to calculate cost for intra
 * \returns cost to code pred block

 ** Only for luma
 */
uint32_t rdo_cost_intra(encoder_state * const encoder_state, pixel *pred, pixel *orig_block, int width, int8_t mode, int tr_depth)
{
    const encoder_control * const encoder = encoder_state->encoder_control;
    coefficient pre_quant_coeff[LCU_WIDTH*LCU_WIDTH>>2];
    int16_t block[LCU_WIDTH*LCU_WIDTH>>2];
    int16_t temp_block[LCU_WIDTH*LCU_WIDTH>>2];
    coefficient temp_coeff[LCU_WIDTH*LCU_WIDTH>>2];
    int8_t luma_scan_mode = SCAN_DIAG;

    int i = 0,x,y;
    for (y = 0; y < width; y++) {
      for (x = 0; x < width; x++) {
        block[i++] = orig_block[x + y*width]- pred[x + y*width];
      }
    }
    // Scan mode is diagonal, except for 4x4 and 8x8, where:
    // - angular 6-14 = vertical
    // - angular 22-30 = horizontal
    if (width <= 8) {
      if (mode >= 6 && mode <= 14) {
        luma_scan_mode = SCAN_VER;
      } else if (mode >= 22 && mode <= 30) {
        luma_scan_mode = SCAN_HOR;
      }
    }
    transform2d(encoder, block,pre_quant_coeff,width,0);
    if(encoder->rdoq_enable) {
      rdoq(encoder_state, pre_quant_coeff, temp_coeff, width, width, 0, luma_scan_mode, CU_INTRA, tr_depth);
    } else {
      quant(encoder_state, pre_quant_coeff, temp_coeff, width, width, 0, luma_scan_mode, CU_INTRA);
    }
    dequant(encoder_state, temp_coeff, pre_quant_coeff, width, width, 0, CU_INTRA);
    itransform2d(encoder, temp_block,pre_quant_coeff,width,0);

    unsigned ssd = 0;
    // SSD between original and reconstructed
    for (i = 0; i < width*width; i++) {
      //int diff = temp_block[i]-block[i];
      int diff = orig_block[i] - CLIP(0, 255, pred[i] + temp_block[i]);

      ssd += diff*diff;
    }

    double coeff_bits = get_coeff_cost(encoder_state, temp_coeff, width, 0, luma_scan_mode);

    return (uint32_t)(0.5 + ssd + coeff_bits * encoder_state->global->cur_lambda_cost);
}


/** Calculate actual (or really close to actual) bitcost for coding coefficients
 * \param coeff coefficient array
 * \param width coeff block width
 * \param type data type (0 == luma)
 * \returns bits needed to code input coefficients
 */
int32_t get_coeff_cost(const encoder_state * const current_encoder_state, coefficient *coeff, int32_t width, int32_t type, int8_t scan_mode)
{
  int32_t cost = 0;
  int i;
  int found = 0;
  encoder_state temp_encoder;

  // Make sure there are coeffs present
  for(i = 0; i < width*width; i++) {
    if (coeff[i] != 0) {
      found = 1;
      break;
    }
  }

  if(!found) return 0;

  // Store cabac state and contexts
  memcpy(&temp_encoder,current_encoder_state,sizeof(encoder_state));

  // Clear bytes and bits and set mode to "count"
  temp_encoder.cabac.only_count = 1;
  temp_encoder.cabac.num_buffered_bytes = 0;
  temp_encoder.cabac.bits_left = 23;

  // Execute the coding function
  encode_coeff_nxn(&temp_encoder, coeff, width, type, scan_mode, 0);

  // Store bitcost before restoring cabac
  cost = (23-temp_encoder.cabac.bits_left) + (temp_encoder.cabac.num_buffered_bytes << 3);

  return cost;
}



#define COEF_REMAIN_BIN_REDUCTION 3
/** Calculates the cost for specific absolute transform level
 * \param abs_level scaled quantized level
 * \param ctx_num_one current ctxInc for coeff_abs_level_greater1 (1st bin of coeff_abs_level_minus1 in AVC)
 * \param ctx_num_abs current ctxInc for coeff_abs_level_greater2 (remaining bins of coeff_abs_level_minus1 in AVC)
 * \param abs_go_rice Rice parameter for coeff_abs_level_minus3
 * \returns cost of given absolute transform level
 * From HM 12.0
 */
int32_t get_ic_rate  (encoder_state * const encoder_state,
                          uint32_t abs_level,
                          uint16_t ctx_num_one,
                          uint16_t ctx_num_abs,
                          uint16_t abs_go_rice,
                          uint32_t c1_idx,
                          uint32_t c2_idx,
                          int8_t type
                          )
{
  cabac_data * const cabac = &encoder_state->cabac;
  int32_t rate = 32768;
  uint32_t base_level  =  (c1_idx < C1FLAG_NUMBER)? (2 + (c2_idx < C2FLAG_NUMBER)) : 1;
  cabac_ctx *base_one_ctx = (type == 0) ? &(cabac->ctx.cu_one_model_luma[0]) : &(cabac->ctx.cu_one_model_chroma[0]);
  cabac_ctx *base_abs_ctx = (type == 0) ? &(cabac->ctx.cu_abs_model_luma[0]) : &(cabac->ctx.cu_abs_model_chroma[0]);

  if ( abs_level >= base_level ) {
    int32_t symbol     = abs_level - base_level;
    int32_t length;
    if (symbol < (COEF_REMAIN_BIN_REDUCTION << abs_go_rice)) {
      length = symbol>>abs_go_rice;
      rate += (length+1+abs_go_rice)<< 15;
    } else {
      length = abs_go_rice;
      symbol  = symbol - ( COEF_REMAIN_BIN_REDUCTION << abs_go_rice);
      while (symbol >= (1<<length)) {
        symbol -=  (1<<(length++));
      }
      rate += (COEF_REMAIN_BIN_REDUCTION+length+1-abs_go_rice+length)<< 15;
    }
    if (c1_idx < C1FLAG_NUMBER) {
      rate += CTX_ENTROPY_BITS(&base_one_ctx[ctx_num_one],1);

      if (c2_idx < C2FLAG_NUMBER) {
        rate += CTX_ENTROPY_BITS(&base_abs_ctx[ctx_num_abs],1);
      }
    }
  }
  else if( abs_level == 1 ) {
    rate += CTX_ENTROPY_BITS(&base_one_ctx[ctx_num_one],0);
  } else if( abs_level == 2 ) {
    rate += CTX_ENTROPY_BITS(&base_one_ctx[ctx_num_one],1);
    rate += CTX_ENTROPY_BITS(&base_abs_ctx[ctx_num_abs],0);
  }

  return rate;
}

/** Get the best level in RD sense
 * \param coded_cost reference to coded cost
 * \param coded_cost0 reference to cost when coefficient is 0
 * \param coded_cost_sig reference to cost of significant coefficient
 * \param level_double reference to unscaled quantized level
 * \param max_abs_level scaled quantized level
 * \param ctx_num_sig current ctxInc for coeff_abs_significant_flag
 * \param ctx_num_one current ctxInc for coeff_abs_level_greater1 (1st bin of coeff_abs_level_minus1 in AVC)
 * \param ctx_num_abs current ctxInc for coeff_abs_level_greater2 (remaining bins of coeff_abs_level_minus1 in AVC)
 * \param abs_go_rice current Rice parameter for coeff_abs_level_minus3
 * \param q_bits quantization step size
 * \param temp correction factor
 * \param last indicates if the coefficient is the last significant
 * \returns best quantized transform level for given scan position
 * This method calculates the best quantized transform level for a given scan position.
 * From HM 12.0
 */
uint32_t get_coded_level ( encoder_state * const encoder_state, double *coded_cost, double *coded_cost0, double *coded_cost_sig,
                           int32_t level_double, uint32_t max_abs_level,
                           uint16_t ctx_num_sig, uint16_t ctx_num_one, uint16_t ctx_num_abs,
                           uint16_t abs_go_rice,
                           uint32_t c1_idx, uint32_t c2_idx,
                           int32_t q_bits,double temp, int8_t last, int8_t type)
{
  cabac_data * const cabac = &encoder_state->cabac;
  double cur_cost_sig   = 0;
  uint32_t best_abs_level = 0;
  int32_t abs_level;
  int32_t min_abs_level;
  cabac_ctx* base_sig_model = type?(cabac->ctx.cu_sig_model_chroma):(cabac->ctx.cu_sig_model_luma);

  if( !last && max_abs_level < 3 ) {
    *coded_cost_sig = encoder_state->global->cur_lambda_cost * CTX_ENTROPY_BITS(&base_sig_model[ctx_num_sig], 0);
    *coded_cost     = *coded_cost0 + *coded_cost_sig;
    if (max_abs_level == 0) return best_abs_level;
  } else {
    *coded_cost = MAX_DOUBLE;
  }

  if( !last ) {
    cur_cost_sig = encoder_state->global->cur_lambda_cost * CTX_ENTROPY_BITS(&base_sig_model[ctx_num_sig], 1);
  }

  min_abs_level    = ( max_abs_level > 1 ? max_abs_level - 1 : 1 );
  for (abs_level = max_abs_level; abs_level >= min_abs_level ; abs_level-- ) {
    double err       = (double)(level_double - ( abs_level << q_bits ) );
    double cur_cost  = err * err * temp + encoder_state->global->cur_lambda_cost *
                       get_ic_rate( encoder_state, abs_level, ctx_num_one, ctx_num_abs,
                                    abs_go_rice, c1_idx, c2_idx, type);
    cur_cost        += cur_cost_sig;

    if( cur_cost < *coded_cost ) {
      best_abs_level  = abs_level;
      *coded_cost     = cur_cost;
      *coded_cost_sig = cur_cost_sig;
    }
  }

  return best_abs_level;
}


/** Calculates the cost of signaling the last significant coefficient in the block
 * \param pos_x X coordinate of the last significant coefficient
 * \param pos_y Y coordinate of the last significant coefficient
 * \returns cost of last significant coefficient
 * \param uiWidth width of the transform unit (TU)
 *
 * From HM 12.0
*/
static double get_rate_last(const encoder_state * const encoder_state,
                            const uint32_t  pos_x, const uint32_t pos_y,
                            int32_t* last_x_bits, int32_t* last_y_bits)
{
  uint32_t ctx_x   = g_group_idx[pos_x];
  uint32_t ctx_y   = g_group_idx[pos_y];
  double uiCost = last_x_bits[ ctx_x ] + last_y_bits[ ctx_y ];
  if( ctx_x > 3 ) {
    uiCost += 32768.0 * ((ctx_x-2)>>1);
  }
  if( ctx_y > 3 ) {
    uiCost += 32768.0 * ((ctx_y-2)>>1);
  }
  return encoder_state->global->cur_lambda_cost*uiCost;
}

static void calc_last_bits(encoder_state * const encoder_state, int32_t width, int32_t height, int8_t type,
                           int32_t* last_x_bits, int32_t* last_y_bits)
{
  cabac_data * const cabac = &encoder_state->cabac;
  int32_t bits_x = 0, bits_y = 0;
  int32_t blk_size_offset_x, blk_size_offset_y, shiftX, shiftY;
  int32_t ctx;

  cabac_ctx *base_ctx_x = (type ? cabac->ctx.cu_ctx_last_x_chroma : cabac->ctx.cu_ctx_last_x_luma);
  cabac_ctx *base_ctx_y = (type ? cabac->ctx.cu_ctx_last_y_chroma : cabac->ctx.cu_ctx_last_y_luma);

  blk_size_offset_x = type ? 0: (g_convert_to_bit[ width ] *3 + ((g_convert_to_bit[ width ] +1)>>2));
  blk_size_offset_y = type ? 0: (g_convert_to_bit[ height ]*3 + ((g_convert_to_bit[ height ]+1)>>2));
  shiftX = type ? g_convert_to_bit[ width  ] :((g_convert_to_bit[ width  ]+3)>>2);
  shiftY = type ? g_convert_to_bit[ height ] :((g_convert_to_bit[ height ]+3)>>2);


  for (ctx = 0; ctx < g_group_idx[ width - 1 ]; ctx++) {
    int32_t ctx_offset = blk_size_offset_x + (ctx >>shiftX);
    last_x_bits[ ctx ] = bits_x + CTX_ENTROPY_BITS(&base_ctx_x[ ctx_offset ],0);
    bits_x += CTX_ENTROPY_BITS(&base_ctx_x[ ctx_offset ],1);
  }
  last_x_bits[ctx] = bits_x;
  for (ctx = 0; ctx < g_group_idx[ height - 1 ]; ctx++) {
    int32_t ctx_offset = blk_size_offset_y + (ctx >>shiftY);
    last_y_bits[ ctx ] = bits_y + CTX_ENTROPY_BITS(&base_ctx_y[ ctx_offset ],0);
    bits_y +=  CTX_ENTROPY_BITS(&base_ctx_y[ ctx_offset ],1);
  }
  last_y_bits[ctx] = bits_y;
}

void rdoq_sign_hiding(const encoder_state *const encoder_state,
                      const int32_t qp_scaled,
                      const uint32_t *const scan,
                      const int32_t delta_u[32 * 32],
                      const int32_t rate_inc_up[32 * 32],
                      const int32_t rate_inc_down[32 * 32],
                      const int32_t sig_rate_delta[32 * 32],
                      const int32_t width,
                      const coefficient *const coef,
                      coefficient *const dest_coeff)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  const int32_t size = width * width;
  
  int64_t rd_factor = (int64_t)(
    g_inv_quant_scales[qp_scaled % 6] * g_inv_quant_scales[qp_scaled % 6] * (1 << (2 * (qp_scaled / 6)))
    / encoder_state->global->cur_lambda_cost / 16 / (1 << (2 * (encoder->bitdepth - 8)))
    + 0.5);
  int32_t lastCG = -1;
  int32_t absSum = 0;
  int32_t n, subset;

  for (subset = (size - 1) >> LOG2_SCAN_SET_SIZE; subset >= 0; subset--) {
    int32_t  subPos = subset << LOG2_SCAN_SET_SIZE;
    int32_t  firstNZPosInCG = SCAN_SET_SIZE, lastNZPosInCG = -1;
    absSum = 0;

    for (n = SCAN_SET_SIZE - 1; n >= 0; --n) {
      if (dest_coeff[scan[n + subPos]]) {
        lastNZPosInCG = n;
        break;
      }
    }

    for (n = 0; n <SCAN_SET_SIZE; n++) {
      if (dest_coeff[scan[n + subPos]]) {
        firstNZPosInCG = n;
        break;
      }
    }

    for (n = firstNZPosInCG; n <= lastNZPosInCG; n++) {
      absSum += dest_coeff[scan[n + subPos]];
    }

    if (lastNZPosInCG >= 0 && lastCG == -1) lastCG = 1;

    if (lastNZPosInCG - firstNZPosInCG >= SBH_THRESHOLD) {
      int32_t signbit = (dest_coeff[scan[subPos + firstNZPosInCG]]>0 ? 0 : 1);
      if (signbit != (absSum & 0x1)) {  // hide but need tune
        // calculate the cost
        int64_t minCostInc = MAX_INT64, curCost = MAX_INT64;
        int32_t minPos = -1, finalChange = 0, curChange = 0;

        for (n = (lastCG == 1 ? lastNZPosInCG : SCAN_SET_SIZE - 1); n >= 0; --n) {
          uint32_t blkpos = scan[n + subPos];
          if (dest_coeff[blkpos] != 0) {
            int64_t costUp = rd_factor * (-delta_u[blkpos]) + rate_inc_up[blkpos];
            int64_t costDown = rd_factor * (delta_u[blkpos]) + rate_inc_down[blkpos]
              - (abs(dest_coeff[blkpos]) == 1 ? ((1 << 15) + sig_rate_delta[blkpos]) : 0);

            if (lastCG == 1 && lastNZPosInCG == n && abs(dest_coeff[blkpos]) == 1) {
              costDown -= (4 << 15);
            }

            if (costUp<costDown) {
              curCost = costUp;
              curChange = 1;
            } else {
              curChange = -1;
              if (n == firstNZPosInCG && abs(dest_coeff[blkpos]) == 1) {
                curCost = MAX_INT64;
              } else {
                curCost = costDown;
              }
            }
          } else {
            curCost = rd_factor * (-(abs(delta_u[blkpos]))) + (1 << 15) + rate_inc_up[blkpos] + sig_rate_delta[blkpos];
            curChange = 1;

            if (n<firstNZPosInCG) {
              if (((coef[blkpos] >= 0) ? 0 : 1) != signbit) curCost = MAX_INT64;
            }
          }

          if (curCost<minCostInc) {
            minCostInc = curCost;
            finalChange = curChange;
            minPos = blkpos;
          }
        }

        if (dest_coeff[minPos] == 32767 || dest_coeff[minPos] == -32768) {
          finalChange = -1;
        }

        if (coef[minPos] >= 0) {
          dest_coeff[minPos] += (coefficient)finalChange;
        } else {
          dest_coeff[minPos] -= (coefficient)finalChange;
        }
      }
    }
    if (lastCG == 1) lastCG = 0;
  }
}

/** RDOQ with CABAC
 * \returns void
 * Rate distortion optimized quantization for entropy
 * coding engines using probability models like CABAC
 * From HM 12.0
 */
void  rdoq(encoder_state * const encoder_state, coefficient *coef, coefficient *dest_coeff, int32_t width,
           int32_t height, int8_t type, int8_t scan_mode, int8_t block_type, int8_t tr_depth)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  cabac_data * const cabac = &encoder_state->cabac;
  uint32_t log2_tr_size    = g_convert_to_bit[ width ] + 2;
  int32_t  transform_shift = MAX_TR_DYNAMIC_RANGE - encoder->bitdepth - log2_tr_size;  // Represents scaling through forward transform
  uint16_t go_rice_param   = 0;
  uint32_t log2_block_size = g_convert_to_bit[ width ] + 2;
  uint32_t max_num_coeff   = width * height;
  int32_t  scalinglist_type= (block_type == CU_INTRA ? 0 : 3) + (int8_t)("\0\3\1\2"[type]);

  int32_t qp_scaled = get_scaled_qp(type, encoder_state->global->QP, 0);
  uint32_t abs_sum = 0;

  
  int32_t q_bits = QUANT_SHIFT + qp_scaled/6 + transform_shift;

  const int32_t *quant_coeff  = encoder->scaling_list.quant_coeff[log2_tr_size-2][scalinglist_type][qp_scaled%6];
  const double *err_scale     = encoder->scaling_list.error_scale[log2_tr_size-2][scalinglist_type][qp_scaled%6];

  double block_uncoded_cost = 0;

  double cost_coeff [ 32 * 32 ];
  double cost_sig   [ 32 * 32 ];
  double cost_coeff0[ 32 * 32 ];

  int32_t rate_inc_up   [ 32 * 32 ];
  int32_t rate_inc_down [ 32 * 32 ];
  int32_t sig_rate_delta[ 32 * 32 ];
  int32_t delta_u       [ 32 * 32 ];


  const uint32_t *scan_cg = g_sig_last_scan_cg[log2_block_size - 2][scan_mode];
  const int32_t  shift   = 4>>1;
  const uint32_t cg_size = 16;
  const uint32_t num_blk_side    = width >> shift;
  double   cost_coeffgroup_sig[ 64 ];
  uint32_t sig_coeffgroup_flag[ 64 ];

  int32_t  cg_last_scanpos = -1;

  uint16_t    ctx_set        = 0;
  int16_t     c1             = 1;
  int16_t     c2             = 0;
  double      base_cost      = 0;
  int32_t     last_scanpos   = -1;

  uint32_t    c1_idx     = 0;
  uint32_t    c2_idx     = 0;
  int32_t     base_level;

  const uint32_t *scan = g_sig_last_scan[ scan_mode ][ log2_block_size - 1 ];


  uint32_t cg_num = width * height >> 4;
  int32_t  scanpos;

  cabac_ctx *base_coeff_group_ctx = &(cabac->ctx.cu_sig_coeff_group_model[type]);
  cabac_ctx *baseCtx              = (type == 0) ? &(cabac->ctx.cu_sig_model_luma[0]) : &(cabac->ctx.cu_sig_model_chroma[0]);
  cabac_ctx *base_one_ctx = (type == 0) ? &(cabac->ctx.cu_one_model_luma[0]) : &(cabac->ctx.cu_one_model_chroma[0]);

  double  best_cost        = 0;
  int32_t ctx_cbf          = 0;
  int32_t best_last_idx_p1 = 0;
  int8_t found_last        = 0;
  int32_t cg_scanpos, scanpos_in_cg;

  coeffgroup_rd_stats rd_stats;

  int32_t last_x_bits[32],last_y_bits[32];
  calc_last_bits(encoder_state, width, height, type,last_x_bits, last_y_bits);

  FILL_ARRAY(cost_coeff, 0, max_num_coeff);
  FILL_ARRAY(cost_sig, 0, max_num_coeff);
  

  if (encoder->sign_hiding) {
    FILL_ARRAY(rate_inc_up, 0, max_num_coeff);
    FILL_ARRAY(rate_inc_down, 0, max_num_coeff);
    FILL_ARRAY(sig_rate_delta, 0, max_num_coeff);
    FILL_ARRAY(delta_u, 0, max_num_coeff);
  }

  FILL(cost_coeffgroup_sig, 0);
  FILL(sig_coeffgroup_flag, 0);

  for (cg_scanpos = cg_num-1; cg_scanpos >= 0; cg_scanpos--) {
    uint32_t cg_blkpos = scan_cg[ cg_scanpos ];
    uint32_t cg_pos_y   = cg_blkpos / num_blk_side;
    uint32_t cg_pos_x   = cg_blkpos - (cg_pos_y * num_blk_side);
    int32_t  scanpos_in_cg;

    int32_t pattern_sig_ctx = context_calc_pattern_sig_ctx(sig_coeffgroup_flag,
                                                           cg_pos_x, cg_pos_y, width);

    FILL(rd_stats, 0);
    for (scanpos_in_cg = cg_size-1; scanpos_in_cg >= 0; scanpos_in_cg--)  {
      uint32_t blkpos;
      int32_t q;
      double temp, err;
      int32_t level_double;
      uint32_t max_abs_level;

      scanpos = cg_scanpos*cg_size + scanpos_in_cg;
      blkpos          = scan[scanpos];
      q  = quant_coeff[blkpos];
      temp = err_scale[blkpos];
      level_double        = coef[blkpos];
      level_double        = MIN(abs(level_double) * q , MAX_INT - (1 << (q_bits - 1)));
      max_abs_level       = (level_double + (1 << (q_bits - 1))) >> q_bits;

      err               = (double)level_double;
      cost_coeff0[ scanpos ]  = err * err * temp;
      block_uncoded_cost      += cost_coeff0[ scanpos ];
      dest_coeff[ blkpos ] = (coefficient)max_abs_level;

      if ( max_abs_level > 0 && last_scanpos < 0 ) {
        last_scanpos             = scanpos;
        ctx_set                  = (scanpos > 0 && type == 0) ? 2 : 0;
        cg_last_scanpos          = cg_scanpos;
      }

      if ( last_scanpos >= 0 ) {
        //===== coefficient level estimation =====
        int32_t  level;
        uint16_t  one_ctx = 4 * ctx_set + c1;
        uint16_t  abs_ctx = ctx_set + c2;

        if( scanpos == last_scanpos ) {
          level            = get_coded_level(encoder_state, &cost_coeff[ scanpos ], &cost_coeff0[ scanpos ], &cost_sig[ scanpos ],
                                               level_double, max_abs_level, 0, one_ctx, abs_ctx, go_rice_param,
                                               c1_idx, c2_idx, q_bits, temp, 1, type );
        } else {
          uint32_t  pos_y    = blkpos >> log2_block_size;
          uint32_t  pos_x    = blkpos - ( pos_y << log2_block_size );
          uint16_t  ctx_sig  = (uint16_t)context_get_sig_ctx_inc(pattern_sig_ctx, scan_mode, pos_x, pos_y,
                                                       log2_block_size, type);
          level              = get_coded_level(encoder_state, &cost_coeff[ scanpos ], &cost_coeff0[ scanpos ], &cost_sig[ scanpos ],
                                               level_double, max_abs_level, ctx_sig, one_ctx, abs_ctx, go_rice_param,
                                               c1_idx, c2_idx, q_bits, temp, 0, type );
          if (encoder->sign_hiding) {
            sig_rate_delta[blkpos] = CTX_ENTROPY_BITS(&baseCtx[ctx_sig], 1) - CTX_ENTROPY_BITS(&baseCtx[ctx_sig], 0);
          }
        }

        if (encoder->sign_hiding) {
          delta_u[blkpos] = (level_double - ((int32_t)level << q_bits)) >> (q_bits - 8);
          if (level > 0) {
            int32_t rate_now = get_ic_rate(encoder_state, level, one_ctx, abs_ctx, go_rice_param, c1_idx, c2_idx, type);
            rate_inc_up[blkpos] = get_ic_rate(encoder_state, level + 1, one_ctx, abs_ctx, go_rice_param, c1_idx, c2_idx, type) - rate_now;
            rate_inc_down[blkpos] = get_ic_rate(encoder_state, level - 1, one_ctx, abs_ctx, go_rice_param, c1_idx, c2_idx, type) - rate_now;
          } else { // level == 0
            rate_inc_up[blkpos] = CTX_ENTROPY_BITS(&base_one_ctx[one_ctx], 0);
          }
        }

        dest_coeff[blkpos] = (coefficient)level;
        base_cost         += cost_coeff[scanpos];

        base_level = (c1_idx < C1FLAG_NUMBER) ? (2 + (c2_idx < C2FLAG_NUMBER)) : 1;
        if( level >= base_level ) {
          if(level  > 3*(1<<go_rice_param)) {
            go_rice_param = MIN(go_rice_param + 1, 4);
          }
        }
        if (level >= 1) c1_idx ++;

        //===== update bin model =====
        if (level > 1) {
          c1 = 0;
          c2 += (c2 < 2);
          c2_idx ++;
        } else if( (c1 < 3) && (c1 > 0) && level) {
          c1++;
        }

        //===== context set update =====
        if ((scanpos % SCAN_SET_SIZE == 0) && scanpos > 0) {
          c2                = 0;
          go_rice_param     = 0;

          c1_idx   = 0;
          c2_idx   = 0;
          ctx_set = (scanpos == SCAN_SET_SIZE || type!=0) ? 0 : 2;
          if( c1 == 0 ) {
            ctx_set++;
          }
          c1 = 1;
        }
      } else {
        base_cost += cost_coeff0[scanpos];
      }
      rd_stats.sig_cost += cost_sig[scanpos];
      if (scanpos_in_cg == 0 ) {
        rd_stats.sig_cost_0 = cost_sig[scanpos];
      }
      if (dest_coeff[ blkpos ] )  {
        sig_coeffgroup_flag[ cg_blkpos ] = 1;
        rd_stats.coded_level_and_dist += cost_coeff[scanpos] - cost_sig[scanpos];
        rd_stats.uncoded_dist += cost_coeff0[scanpos];
        if ( scanpos_in_cg != 0 ) {
          rd_stats.nnz_before_pos0++;
        }
      }
    } //end for (scanpos_in_cg)

    if (cg_last_scanpos >= 0) {
      if( cg_scanpos ) {
        if (sig_coeffgroup_flag[ cg_blkpos ] == 0) {
          uint32_t ctx_sig  = context_get_sig_coeff_group(sig_coeffgroup_flag, cg_pos_x,
                                                          cg_pos_y, width);
          cost_coeffgroup_sig[ cg_scanpos ] = encoder_state->global->cur_lambda_cost*CTX_ENTROPY_BITS(&base_coeff_group_ctx[ctx_sig],0);
          base_cost += cost_coeffgroup_sig[ cg_scanpos ]  - rd_stats.sig_cost;
        } else {
          if (cg_scanpos < cg_last_scanpos) {//skip the last coefficient group, which will be handled together with last position below.
            double cost_zero_cg;
            uint32_t ctx_sig;
            if (rd_stats.nnz_before_pos0 == 0) {
              base_cost -= rd_stats.sig_cost_0;
              rd_stats.sig_cost -= rd_stats.sig_cost_0;
            }
            // rd-cost if SigCoeffGroupFlag = 0, initialization
            cost_zero_cg = base_cost;

            // add SigCoeffGroupFlag cost to total cost
            ctx_sig  = context_get_sig_coeff_group(sig_coeffgroup_flag, cg_pos_x,
                                                            cg_pos_y, width);
            if (cg_scanpos < cg_last_scanpos) {
              cost_coeffgroup_sig[cg_scanpos] = encoder_state->global->cur_lambda_cost*CTX_ENTROPY_BITS(&base_coeff_group_ctx[ctx_sig],1);
              base_cost    += cost_coeffgroup_sig[cg_scanpos];
              cost_zero_cg += encoder_state->global->cur_lambda_cost*CTX_ENTROPY_BITS(&base_coeff_group_ctx[ctx_sig],0);
            }

            // try to convert the current coeff group from non-zero to all-zero
            cost_zero_cg += rd_stats.uncoded_dist;          // distortion for resetting non-zero levels to zero levels
            cost_zero_cg -= rd_stats.coded_level_and_dist;  // distortion and level cost for keeping all non-zero levels
            cost_zero_cg -= rd_stats.sig_cost;              // sig cost for all coeffs, including zero levels and non-zerl levels

            // if we can save cost, change this block to all-zero block
            if (cost_zero_cg < base_cost) {
              int32_t scanpos_in_cg;
              sig_coeffgroup_flag[ cg_blkpos ] = 0;
              base_cost = cost_zero_cg;
              if (cg_scanpos < cg_last_scanpos) {
                cost_coeffgroup_sig[ cg_scanpos ] = encoder_state->global->cur_lambda_cost*CTX_ENTROPY_BITS(&base_coeff_group_ctx[ctx_sig],0);
              }
              // reset coeffs to 0 in this block
              for (scanpos_in_cg = cg_size-1; scanpos_in_cg >= 0; scanpos_in_cg--) {
                uint32_t blkpos;
                scanpos      = cg_scanpos*cg_size + scanpos_in_cg;
                blkpos = scan[ scanpos ];

                if (dest_coeff[ blkpos ]) {
                  dest_coeff[ blkpos ]  = 0;
                  cost_coeff[ scanpos ] = cost_coeff0[ scanpos ];
                  cost_sig  [ scanpos ] = 0;
                }
              }
            } // end if ( cost_all_zeros < base_cost )
          }
        } // end if if (sig_coeffgroup_flag[ cg_blkpos ] == 0)
      } else {
        sig_coeffgroup_flag[ cg_blkpos ] = 1;
      }
    }
  } //end for (cg_scanpos)

  //===== estimate last position =====
  if (last_scanpos < 0) return;


  if( block_type != CU_INTRA && !type/* && pcCU->getTransformIdx( uiAbsPartIdx ) == 0*/ ) {
    best_cost  = block_uncoded_cost +   encoder_state->global->cur_lambda_cost*CTX_ENTROPY_BITS(&(cabac->ctx.cu_qt_root_cbf_model),0);
    base_cost +=   encoder_state->global->cur_lambda_cost*CTX_ENTROPY_BITS(&(cabac->ctx.cu_qt_root_cbf_model),1);
  } else {
    cabac_ctx* base_cbf_model = type?(cabac->ctx.qt_cbf_model_chroma):(cabac->ctx.qt_cbf_model_luma);
    ctx_cbf   = ( type ? tr_depth : !tr_depth);
    best_cost  = block_uncoded_cost +  encoder_state->global->cur_lambda_cost*CTX_ENTROPY_BITS(&base_cbf_model[ctx_cbf],0);
    base_cost +=   encoder_state->global->cur_lambda_cost*CTX_ENTROPY_BITS(&base_cbf_model[ctx_cbf],1);
  }

  for (cg_scanpos = cg_last_scanpos; cg_scanpos >= 0; cg_scanpos--) {
    uint32_t cg_blkpos = scan_cg[cg_scanpos];

    base_cost -= cost_coeffgroup_sig[cg_scanpos];
    if (sig_coeffgroup_flag[ cg_blkpos ]) {
      for (scanpos_in_cg = cg_size-1; scanpos_in_cg >= 0; scanpos_in_cg--) {
        uint32_t   blkpos;
        scanpos = cg_scanpos*cg_size + scanpos_in_cg;
        if (scanpos > last_scanpos) continue;
        blkpos  = scan[scanpos];

        if( dest_coeff[ blkpos ] ) {
          uint32_t   pos_y       = blkpos >> log2_block_size;
          uint32_t   pos_x       = blkpos - ( pos_y << log2_block_size );

          double cost_last = (scan_mode == SCAN_VER) ? get_rate_last(encoder_state, pos_y, pos_x,last_x_bits,last_y_bits) : get_rate_last(encoder_state, pos_x, pos_y, last_x_bits,last_y_bits );
          double totalCost = base_cost + cost_last - cost_sig[ scanpos ];

          if( totalCost < best_cost ) {
            best_last_idx_p1  = scanpos + 1;
            best_cost         = totalCost;
          }
          if( dest_coeff[ blkpos ] > 1 ) {
            found_last = 1;
            break;
          }
          base_cost  -= cost_coeff[ scanpos ];
          base_cost  += cost_coeff0[ scanpos ];
        } else {
          base_cost  -= cost_sig[ scanpos ];
        }
      } //end for
      if (found_last) break;
    } // end if (sig_coeffgroup_flag[ cg_blkpos ])
  } // end for

  for ( scanpos = 0; scanpos < best_last_idx_p1; scanpos++ ) {
    int32_t blkPos = scan[ scanpos ];
    int32_t level  = dest_coeff[ blkPos ];
    abs_sum += level;
    dest_coeff[ blkPos ] = (coefficient)(( coef[ blkPos ] < 0 ) ? -level : level);
  }

  //===== clean uncoded coefficients =====
  for ( scanpos = best_last_idx_p1; scanpos <= last_scanpos; scanpos++ ) {
    dest_coeff[ scan[ scanpos ] ] = 0;
  }

  if (encoder->sign_hiding && abs_sum >= 2) {
    rdoq_sign_hiding(encoder_state, qp_scaled, scan,
                     delta_u, rate_inc_up, rate_inc_down, sig_rate_delta,
                     width, coef, dest_coeff);
  }
}
