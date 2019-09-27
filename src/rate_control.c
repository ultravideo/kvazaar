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

#include "rate_control.h"

#include <math.h>

#include "encoder.h"
#include "kvazaar.h"


static const int SMOOTHING_WINDOW = 40;
static const double MIN_LAMBDA    = 0.1;
static const double MAX_LAMBDA    = 10000;

/**
 * \brief Clip lambda value to a valid range.
 */
static double clip_lambda(double lambda) {
  if (isnan(lambda)) return MAX_LAMBDA;
  return CLIP(MIN_LAMBDA, MAX_LAMBDA, lambda);
}

/**
 * \brief Update alpha and beta parameters.
 *
 * \param         bits        number of bits spent for coding the area
 * \param         pixels      size of the area in pixels
 * \param         lambda_real lambda used for coding the area
 * \param[in,out] alpha       alpha parameter to update
 * \param[in,out] beta        beta parameter to update
 */
static void update_parameters(uint32_t bits,
                              uint32_t pixels,
                              double lambda_real,
                              double *alpha,
                              double *beta)
{
  const double bpp              = bits / (double)pixels;
  const double lambda_comp      = clip_lambda(*alpha * pow(bpp, *beta));
  const double lambda_log_ratio = log(lambda_real) - log(lambda_comp);

  *alpha += 0.10 * lambda_log_ratio * (*alpha);
  *alpha = CLIP(0.05, 20, *alpha);

  *beta  += 0.05 * lambda_log_ratio * CLIP(-5.0, -1.0, log(bpp));
  *beta  = CLIP(-3, -0.1, *beta);
}

/**
 * \brief Allocate bits for the current GOP.
 * \param state   the main encoder state
 * \return        target number of bits
 */
static double gop_allocate_bits(encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;

  // At this point, total_bits_coded of the current state contains the
  // number of bits written encoder->owf frames before the current frame.
  uint64_t bits_coded = state->frame->total_bits_coded;
  int pictures_coded = MAX(0, state->frame->num - encoder->cfg.owf);

  int gop_offset = (state->frame->gop_offset - encoder->cfg.owf) % MAX(1, encoder->cfg.gop_len);
  
  if (encoder->cfg.gop_len > 0 && gop_offset != encoder->cfg.gop_len - 1 && encoder->cfg.gop_lp_definition.d == 0) {
    // Subtract number of bits in the partially coded GOP.
    bits_coded -= state->frame->cur_gop_bits_coded;
    // Subtract number of pictures in the partially coded GOP.
    pictures_coded -= gop_offset + 1;
  }

  // Equation 12 from https://doi.org/10.1109/TIP.2014.2336550
  double gop_target_bits =
    (encoder->target_avg_bppic * (pictures_coded + SMOOTHING_WINDOW) - bits_coded)
    * MAX(1, encoder->cfg.gop_len) / SMOOTHING_WINDOW;
  // Allocate at least 200 bits for each GOP like HM does.
  return MAX(200, gop_target_bits);
}

/**
 * Estimate number of bits used for headers of the current picture.
 * \param state   the main encoder state
 * \return        number of header bits
 */
static uint64_t pic_header_bits(encoder_state_t * const state)
{
  const kvz_config* cfg = &state->encoder_control->cfg;

  // nal type and slice header
  uint64_t bits = 48 + 24;

  // entry points
  bits += 12 * state->encoder_control->in.height_in_lcu;

  switch (cfg->hash) {
    case KVZ_HASH_CHECKSUM:
      bits += 168;
      break;

    case KVZ_HASH_MD5:
      bits += 456;
      break;

    case KVZ_HASH_NONE:
      break;
  }

  if (encoder_state_must_write_vps(state)) {
    bits += 613;
  }

  if (state->frame->num == 0 && cfg->add_encoder_info) {
    bits += 1392;
  }

  return bits;
}

/**
 * Allocate bits for the current picture.
 * \param state   the main encoder state
 * \return        target number of bits, excluding headers
 */
static double pic_allocate_bits(encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;

  if (encoder->cfg.gop_len == 0 ||
      state->frame->gop_offset == 0 ||
      state->frame->num == 0)
  {
    // A new GOP starts at this frame.
    state->frame->cur_gop_target_bits = gop_allocate_bits(state);
    state->frame->cur_gop_bits_coded  = 0;
  } else {
    state->frame->cur_gop_target_bits =
      state->previous_encoder_state->frame->cur_gop_target_bits;
  }

  if (encoder->cfg.gop_len <= 0) {
    return state->frame->cur_gop_target_bits;
  }

  const double pic_weight = encoder->gop_layer_weights[
    encoder->cfg.gop[state->frame->gop_offset].layer - 1];
  const double pic_target_bits =
    state->frame->cur_gop_target_bits * pic_weight - pic_header_bits(state);
  // Allocate at least 100 bits for each picture like HM does.
  return MAX(100, pic_target_bits);
}

static int8_t lambda_to_qp(const double lambda)
{
  const int8_t qp = 4.2005 * log(lambda) + 13.7223 + 0.5;
  return CLIP_TO_QP(qp);
}

static double solve_cubic_equation(const encoder_state_config_frame_t * const state,
                            int ctu_index,
                            int last_ctu,
                            int layer,
                            double est_lambda,
                            double target_bits) 
{
  double best_lambda = 0.0;
  double para_a = 0.0;
  double para_b = 0.0;
  double para_c = 0.0;
  double para_d = 0.0;
  double delta = 0.0;
  double para_aa = 0.0;
  double para_bb = 0.0;
  double para_cc = 0.0;
  for (int i = ctu_index; i < last_ctu; i++)
  {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double d = 0.0;
    assert(!((state->new_ratecontrol.c_para[layer][i] <= 0) || (state->new_ratecontrol.k_para[layer][i] >= 0))); //Check C and K during each solution 

    double CLCU = state->new_ratecontrol.c_para[layer][i];
    double KLCU = state->new_ratecontrol.k_para[layer][i];
    a = -CLCU * KLCU / pow(state->lcu_stats[i].pixels, KLCU - 1.0);
    b = -1.0 / (KLCU - 1.0);
    d = est_lambda;
    c = pow(a / d, b);
    para_a = para_a - c * pow(b, 3.0) / 6.0;
    para_b = para_b + (pow(b, 2.0) / 2.0 + pow(b, 3.0)*log(d) / 2.0)*c;
    para_c = para_c - (pow(b, 3.0) / 2.0*pow(log(d), 2.0) + pow(b, 2.0)*log(d) + b)*c;
    para_d = para_d + c * (1 + b * log(d) + pow(b, 2.0) / 2 * pow(log(d), 2.0) + pow(b, 3.0) / 6 * pow(log(d), 3.0));
  }

  para_d = para_d - target_bits;
  para_aa = para_b * para_b - 3 * para_a*para_c;
  para_bb = para_b * para_c - 9 * para_a*para_d;
  para_cc = para_c * para_c - 3 * para_b*para_d;

  delta = para_bb * para_bb - 4 * para_aa*para_cc;

  if (delta > 0.0)	//Check whether delta is right
  {
    double temp_x = 0.0;
    double part1 = 0.0;
    double part2 = 0.0;
    double flag1 = 0.0;
    double flag2 = 0.0;
    part1 = para_aa * para_b + 3 * para_a*(-para_bb - pow(delta, 0.5)) / 2.0;
    part2 = para_aa * para_b + 3 * para_a*(-para_bb + pow(delta, 0.5)) / 2.0;
    if (part1 < 0.0) {
      part1 = -part1;
      flag1 = -1.0;
    }
    else {
      flag1 = 1.0;
    }
    if (part2 < 0.0) {
      part2 = -part2;
      flag2 = -1.0;
    }
    else {
      flag2 = 1.0;
    }
    temp_x = (-para_b - flag1 * pow(part1, 1.0 / 3.0) - flag2 * pow(part2, 1.0 / 3.0)) / 3 / para_a;
    best_lambda = exp(temp_x);
  }
  else {
    best_lambda = est_lambda;		//Use the original picture estimated lambda for the current CTU
  }
  best_lambda = CLIP(0.001, 100000000.0, best_lambda);

  return best_lambda;
}

static INLINE double calculate_weights(encoder_state_t* const state, const int layer, const int ctu_count, double estLambda) {
  double total_weight = 0;
  for(int i = 0; i < ctu_count; i++) {
    double c_lcu = state->frame->new_ratecontrol.c_para[layer][i];
    double k_lcu = state->frame->new_ratecontrol.k_para[layer][i];
    double a = -c_lcu * k_lcu / pow(state->frame->lcu_stats[i].pixels, k_lcu - 1.0);
    double b = -1.0 / (k_lcu - 1.0);
    state->frame->lcu_stats[i].weight = pow(a / estLambda, b);
    if (state->frame->lcu_stats[i].weight < 0.01) {
      state->frame->lcu_stats[i].weight = 0.01;
    }
    total_weight += state->frame->lcu_stats[i].weight;
  }
  return total_weight;
}


void kvz_estimate_pic_lambda(encoder_state_t * const state) {
  const encoder_control_t * const encoder = state->encoder_control;

  double bits = pic_allocate_bits(state);
  state->frame->cur_pic_target_bits = bits;
  const int layer = encoder->gop_layer_weights[
    encoder->cfg.gop[state->frame->gop_offset].layer] - (state->frame->is_irap ? 1 : 0);
  const int ctu_count = state->tile->frame->height_in_lcu * state->tile->frame->width_in_lcu;

  double alpha;
  double beta;  
  if(state->frame->poc == 0) {
    alpha = state->frame->rc_alpha;
    beta = state->frame->rc_beta;
  }
  else {
    alpha = -state->frame->new_ratecontrol.pic_c_para[state->frame->gop_offset] *
      state->frame->new_ratecontrol.pic_k_para[state->frame->gop_offset];
    beta = state->frame->new_ratecontrol.pic_k_para[state->frame->gop_offset] - 1;
  }
  double est_lambda;
  double bpp = bits / (state->encoder_control->cfg.width * state->encoder_control->cfg.height);
  if (state->frame->is_irap) {
    // TODO: Intra
    est_lambda = alpha * pow(bpp, beta) * 0.5;
  }
  else {
    est_lambda = alpha * pow(bpp, beta);
  }

  double temp_lambda;
  if ((temp_lambda = state->frame->new_ratecontrol.previous_lambdas[layer]) > 0.0) {
    est_lambda = CLIP(temp_lambda * pow(2.0, -1), temp_lambda * 2, est_lambda);
  }

  if((temp_lambda = state->frame->new_ratecontrol.previous_frame_lambda) > 0.0) {
    est_lambda = CLIP(temp_lambda * pow(2.0, -10.0 / 3.0), temp_lambda * pow(2.0, 10.0 / 3.0), est_lambda);
  }

  est_lambda = MAX(est_lambda, 0.1);

  double total_weight = 0;

  if(!state->frame->is_irap) {
    if(!state->encoder_control->cfg.frame_allocation) {
      double best_lambda = 0.0;
      temp_lambda = est_lambda;
      double taylor_e3;
      int iteration_number = 0;
      do {
        taylor_e3 = 0.0;
        best_lambda = temp_lambda = solve_cubic_equation(state->frame, 0, ctu_count, layer, temp_lambda, bits);
        for (int i = 0; i < ctu_count; ++i) {
          double CLCU = state->frame->new_ratecontrol.c_para[layer][i];
          double KLCU = state->frame->new_ratecontrol.k_para[layer][i];
          double a = -CLCU * KLCU / pow(state->frame->lcu_stats[i].pixels, KLCU - 1.0);
          double b = -1.0 / (KLCU - 1.0);
          taylor_e3 += pow(a / best_lambda, b);
        }
      }
      while (fabs(taylor_e3 - bits) > 0.01 && iteration_number <= 11);
    }
    total_weight = calculate_weights(state, layer, ctu_count, est_lambda);
  }
  else {
    for (int i = 0; i < ctu_count; ++i) {
      state->frame->lcu_stats[i].weight = MAX(0.01,
        state->frame->lcu_stats[i].pixels * pow(est_lambda / state->frame->rc_alpha,
                                                1.0 / state->frame->rc_beta));
      total_weight += state->frame->lcu_stats[i].weight;
    }
  }

  for(int i = 0; i < ctu_count; ++i) {
    state->frame->lcu_stats[i].weight = bits * state->frame->lcu_stats[i].weight / total_weight;
  }

  state->frame->lambda = est_lambda;
  state->frame->QP = lambda_to_qp(est_lambda);
}


static double get_ctu_bits(encoder_state_t * const state, vector2d_t pos) {
  int avg_bits;
  const encoder_control_t * const encoder = state->encoder_control;

  const int layer = encoder->gop_layer_weights[
    encoder->cfg.gop[state->frame->gop_offset].layer] - (state->frame->is_irap ? 1 : 0);

  const int num_ctu = state->encoder_control->in.width_in_lcu * state->encoder_control->in.height_in_lcu;
  const int index = pos.x + pos.y * state->tile->frame->width_in_lcu;

  if (state->frame->is_irap) {
    // TODO: intra
    avg_bits = state->frame->cur_pic_target_bits * ((double)state->frame->lcu_stats[index].pixels / 
      (state->encoder_control->in.height * state->encoder_control->in.width));
  }
  else {
    double total_weight = 0;
    const int used_ctu_count = MIN(4, num_ctu - index); //g_RCLCUSmoothWindowSize, the same as the original RC scheme
    int target_bits = 0;
    double best_lambda = 0.0;
    double temp_lambda = state->frame->lambda;
    double taylor_e3 = 0.0;
    int iter = 0;
    double est_lambda = temp_lambda;

    for (int i = index; i < num_ctu; i++) {
      total_weight += state->frame->lcu_stats[i].weight;
    }

    int last_ctu = index + used_ctu_count;
    for (int i = index; i < last_ctu; i++) {
      target_bits += state->frame->lcu_stats[i].weight;
    }

    target_bits = MAX(target_bits + state->frame->total_bits_coded - (int)total_weight, 10); //obtain the total bit-rate for the realInfluenceLCU (=4) CTUs

    //just similar with the process at frame level, details can refer to the function TEncRCPic::kvz_estimate_pic_lambda
    do {
      taylor_e3 = 0.0;
      best_lambda = solve_cubic_equation(state->frame, index, last_ctu, layer, temp_lambda, target_bits);
      temp_lambda = best_lambda;
      for (int i = index; i < last_ctu; i++) {

        double CLCU = state->frame->new_ratecontrol.c_para[layer][i];
        double KLCU = state->frame->new_ratecontrol.k_para[layer][i];
        double a = -CLCU * KLCU / pow((double)state->frame->lcu_stats[i].pixels, KLCU - 1.0);
        double b = -1.0 / (KLCU - 1.0);
        taylor_e3 += pow(a / best_lambda, b);
      }
      iter++;
    } while (fabs(taylor_e3 - target_bits) > 0.01 && iter < 5);

    double c_ctu = state->frame->new_ratecontrol.c_para[layer][index];
    double k_ctu = state->frame->new_ratecontrol.k_para[layer][index];
    double a = -c_ctu * k_ctu / pow(((double)state->frame->lcu_stats[index].pixels), k_ctu - 1.0);
    double b = -1.0 / (k_ctu - 1.0);

    state->frame->lcu_stats[index].weight = MAX(pow(a / best_lambda, b), 0.01);

    avg_bits = (int)(state->frame->lcu_stats[index].weight + 0.5);
  }

  if (avg_bits < 1) {
    avg_bits = 1;
  }
  return avg_bits;
}


void kvz_set_ctu_qp_lambda(encoder_state_t * const state, vector2d_t pos) {
  double bits = get_ctu_bits(state, pos);

  const encoder_control_t * const encoder = state->encoder_control;
  const int frame_allocation = state->encoder_control->cfg.frame_allocation;
  const int layer = encoder->gop_layer_weights[
    encoder->cfg.gop[state->frame->gop_offset].layer] - (state->frame->is_irap ? 1 : 0);

  int index = pos.x + pos.y * state->encoder_control->in.width_in_lcu;
  double bpp = bits / state->frame->lcu_stats[index].pixels;

  double alpha;
  double beta;
  if (state->frame->poc == 0) {
    alpha = state->frame->rc_alpha;
    beta = state->frame->rc_beta;
  }
  else {
    alpha = -state->frame->new_ratecontrol.c_para[layer][index] *
      state->frame->new_ratecontrol.k_para[layer][index];
    beta = state->frame->new_ratecontrol.k_para[layer][index] - 1;
  }

  double est_lambda = alpha * pow(bpp, beta);
  double clip_lambda = state->frame->lambda;

  double clip_neighbor_lambda = -1;
  for(int temp_index = index - 1; temp_index >= 0; --temp_index) {
    if(state->frame->lcu_stats[index].lambda > 0) {
      clip_neighbor_lambda = state->frame->lcu_stats[index].lambda;
      break;
    }
  }

  if (clip_neighbor_lambda > 0) {
    est_lambda = CLIP(clip_neighbor_lambda * pow(2, -(1.0 + frame_allocation) / 3.0),
      clip_neighbor_lambda * pow(2.0, (1.0 + frame_allocation) / 3.0),
      est_lambda);
  }

  if (clip_lambda > 0) {
    est_lambda = CLIP(clip_lambda * pow(2, -(2.0 + frame_allocation) / 3.0),
      clip_lambda * pow(2.0, (1.0 + frame_allocation) / 3.0),
      est_lambda);
  }
  else {
    est_lambda = CLIP(10.0, 1000.0, est_lambda);
  }

  if (est_lambda < 0.1) {
    est_lambda = 0.1;
  }

  int est_qp = lambda_to_qp(est_lambda);

  int clip_qp = -1;
  for (int temp_index = index - 1; temp_index >= 0; --temp_index) {
    if (state->frame->lcu_stats[index].qp > -1) {
      clip_qp = state->frame->lcu_stats[index].qp;
      break;
    }
  }

  if( clip_qp > -1) {
    est_qp = CLIP(clip_qp - 1 - frame_allocation,
      clip_qp + 1 + frame_allocation,
      clip_qp);
  }

  est_qp = CLIP(state->frame->QP - 2 - frame_allocation,
    state->frame->QP + 2 + frame_allocation,
    est_qp);

  state->lambda = est_lambda;
  state->lambda_sqrt = sqrt(est_lambda);
  state->qp = est_qp;
  state->frame->lcu_stats[index].qp = est_qp;
  state->frame->lcu_stats[index].lambda = est_lambda;
}


static void update_pic_ck(encoder_state_t * const state, double bpp, double distortion, double lambda, int layer) {
  double new_k, new_c;
  if(state->frame->num == 1) {
    new_k = log(distortion / state->frame->new_ratecontrol.intra_pic_distortion) /
      log(bpp / state->frame->new_ratecontrol.intra_pic_bpp);
    new_c = distortion / pow(bpp, new_k);
  }
  else {
    new_k = -bpp * lambda / distortion;
    new_c = distortion / pow(bpp, new_k);
  }
  new_c = CLIP(+.1, 100.0, new_c);
  new_k = CLIP(-3.0, -0.001, new_k);

  if(state->frame->is_irap || state->frame->num <= (4 - state->encoder_control->cfg.frame_allocation)) {
    for(int i = 1; i < 5; i++) {
      state->frame->new_ratecontrol.pic_c_para[i] = new_c;
      state->frame->new_ratecontrol.pic_k_para[i] = new_k;
    }
  }
  else {
    state->frame->new_ratecontrol.pic_c_para[layer] = new_c;
    state->frame->new_ratecontrol.pic_k_para[layer] = new_k;
  }
}


static void update_ck(encoder_state_t * const state, int ctu_index, int layer)
{
  double bpp = (double)state->frame->lcu_stats[ctu_index].bits / state->frame->lcu_stats[ctu_index].pixels;
  double distortion = state->frame->lcu_stats[ctu_index].distortion;
  double lambda = state->frame->lcu_stats[ctu_index].lambda;

  double new_k, new_c;
  if (!state->frame->lcu_stats[ctu_index].skipped) {
    if (state->frame->num == 1) {
      new_k = log(distortion / state->frame->new_ratecontrol.intra_pic_distortion) /
        log(bpp / state->frame->new_ratecontrol.intra_pic_bpp);
      new_c = distortion / pow(bpp, new_k);
    }
    else {
      bpp = CLIP(0.0001, 10.0, bpp);
      new_k = -bpp * lambda / distortion;
      new_c = distortion / pow(bpp, new_k);
    }
    new_c = CLIP(+.1, 100.0, new_c);
    new_k = CLIP(-3.0, -0.001, new_k);

    if (state->frame->is_irap || state->frame->num <= (4 - state->encoder_control->cfg.frame_allocation)) {
      for (int i = 1; i < 5; i++) {
        state->frame->new_ratecontrol.c_para[i][ctu_index] = new_c;
        state->frame->new_ratecontrol.k_para[i][ctu_index] = new_k;
      }
    }
    else {
      state->frame->new_ratecontrol.c_para[layer][ctu_index] = new_c;
      state->frame->new_ratecontrol.k_para[layer][ctu_index] = new_k;
    }
  }
}


void kvz_update_after_picture(encoder_state_t * const state) {
  double total_distortion = 0;
  double lambda = 0;
  double pic_bpp = (double)state->frame->cur_frame_bits_coded / (state->encoder_control->in.width * state->encoder_control->in.height);

  const encoder_control_t * const encoder = state->encoder_control;
  const int layer = encoder->gop_layer_weights[
    encoder->cfg.gop[state->frame->gop_offset].layer] - (state->frame->is_irap ? 1 : 0);
  for(int y_ctu = 0; y_ctu < state->encoder_control->in.height_in_lcu; y_ctu++) {
    for (int x_ctu = 0; x_ctu < state->encoder_control->in.width_in_lcu; x_ctu++) {
      int ctu_distortion = 0;
      lcu_stats_t *ctu = kvz_get_lcu_stats(state, x_ctu, y_ctu);
      for (int y = y_ctu * 64; y < MIN((y_ctu + 1) * 64, state->tile->frame->height); y++) {
        for (int x = x_ctu * 64; x < MIN((x_ctu + 1) * 64, state->tile->frame->width); x++) {
          int temp = (int)state->tile->frame->source->y[x + y * state->encoder_control->in.width] -
            state->tile->frame->rec->y[x + y * state->encoder_control->in.width];
          ctu_distortion += temp * temp;
        }        
      }
      ctu->distortion = (double)ctu_distortion / ctu->pixels;
      total_distortion += (double)ctu_distortion / ctu->pixels;
      lambda += ctu->lambda / (state->encoder_control->in.width_in_lcu * state->encoder_control->in.height_in_lcu);
    }    
  }

  if (state->frame->is_irap) {
    for (int y_ctu = 0; y_ctu < state->encoder_control->in.height_in_lcu; y_ctu++) {
      for (int x_ctu = 0; x_ctu < state->encoder_control->in.width_in_lcu; x_ctu++) {
        lcu_stats_t *ctu = kvz_get_lcu_stats(state, x_ctu, y_ctu);
        state->frame->new_ratecontrol.intra_dis[x_ctu + y_ctu * state->encoder_control->in.width_in_lcu] =
          ctu->distortion;
        state->frame->new_ratecontrol.intra_bpp[x_ctu + y_ctu * state->encoder_control->in.width_in_lcu] =
          ctu->bits / ctu->pixels;
      }
    }
    state->frame->new_ratecontrol.intra_pic_distortion = total_distortion;
    state->frame->new_ratecontrol.intra_pic_bpp = pic_bpp;
  }

  state->frame->new_ratecontrol.previous_frame_lambda = lambda;
  state->frame->new_ratecontrol.previous_lambdas[layer] = lambda;
 
  update_pic_ck(state, pic_bpp, total_distortion, lambda, layer);
  for(int i = 0; i < state->encoder_control->in.width_in_lcu * state->encoder_control->in.height_in_lcu; i++) {
    update_ck(state, i, layer);
  }
}


static double qp_to_lamba(encoder_state_t * const state, int qp)
{
  const encoder_control_t * const ctrl = state->encoder_control;
  const int gop_len = ctrl->cfg.gop_len;
  const int period = gop_len > 0 ? gop_len : ctrl->cfg.intra_period;

  kvz_gop_config const * const gop = &ctrl->cfg.gop[state->frame->gop_offset];

  double lambda = pow(2.0, (qp - 12) / 3.0);

  if (state->frame->slicetype == KVZ_SLICE_I) {
    lambda *= 0.57;

    // Reduce lambda for I-frames according to the number of references.
    if (period == 0) {
      lambda *= 0.5;
    } else {
      lambda *= 1.0 - CLIP(0.0, 0.5, 0.05 * (period - 1));
    }
  } else if (gop_len > 0) {
    lambda *= gop->qp_factor;
  } else {
    lambda *= 0.4624;
  }

  // Increase lambda if not key-frame.
  if (period > 0 && state->frame->poc % period != 0) {
    lambda *= CLIP(2.0, 4.0, (state->frame->QP - 12) / 6.0);
  }

  return lambda;
}

/**
 * \brief Allocate bits and set lambda and QP for the current picture.
 * \param state the main encoder state
 */
void kvz_set_picture_lambda_and_qp(encoder_state_t * const state)
{
  const encoder_control_t * const ctrl = state->encoder_control;

  if (ctrl->cfg.target_bitrate > 0) {
    // Rate control enabled

    if (state->frame->num > ctrl->cfg.owf) {
      // At least one frame has been written.
      update_parameters(state->stats_bitstream_length * 8,
                        ctrl->in.pixels_per_pic,
                        state->frame->lambda,
                        &state->frame->rc_alpha,
                        &state->frame->rc_beta);
    }

    const double pic_target_bits = pic_allocate_bits(state);
    const double target_bpp = pic_target_bits / ctrl->in.pixels_per_pic;
    double lambda = state->frame->rc_alpha * pow(target_bpp, state->frame->rc_beta);
    lambda = clip_lambda(lambda);

    state->frame->lambda              = lambda;
    state->frame->QP                  = lambda_to_qp(lambda);
    state->frame->cur_pic_target_bits = pic_target_bits;

  } else {
    // Rate control disabled
    kvz_gop_config const * const gop = &ctrl->cfg.gop[state->frame->gop_offset];
    const int gop_len = ctrl->cfg.gop_len;

    if (gop_len > 0 && state->frame->slicetype != KVZ_SLICE_I) {
      state->frame->QP = CLIP_TO_QP(ctrl->cfg.qp + gop->qp_offset);
    } else {
      state->frame->QP = ctrl->cfg.qp;
    }

    state->frame->lambda = qp_to_lamba(state, state->frame->QP);
  }
}

/**
 * \brief Allocate bits for a LCU.
 * \param state   the main encoder state
 * \param pos     location of the LCU as number of LCUs from top left
 * \return number of bits allocated for the LCU
 */
static double lcu_allocate_bits(encoder_state_t * const state,
                                vector2d_t pos)
{
  double lcu_weight;
  if (state->frame->num > state->encoder_control->cfg.owf) {
    lcu_weight = kvz_get_lcu_stats(state, pos.x, pos.y)->weight;
  } else {
    const uint32_t num_lcus = state->encoder_control->in.width_in_lcu *
                              state->encoder_control->in.height_in_lcu;
    lcu_weight = 1.0 / num_lcus;
  }

  // Target number of bits for the current LCU.
  const double lcu_target_bits = state->frame->cur_pic_target_bits * lcu_weight;

  // Allocate at least one bit for each LCU.
  return MAX(1, lcu_target_bits);
}

void kvz_set_lcu_lambda_and_qp(encoder_state_t * const state,
                               vector2d_t pos)
{
  const encoder_control_t * const ctrl = state->encoder_control;

  if (ctrl->cfg.roi.dqps != NULL) {
    vector2d_t lcu = {
      pos.x + state->tile->lcu_offset_x,
      pos.y + state->tile->lcu_offset_y
    };
    vector2d_t roi = {
      lcu.x * ctrl->cfg.roi.width / ctrl->in.width_in_lcu,
      lcu.y * ctrl->cfg.roi.height / ctrl->in.height_in_lcu
    };
    int roi_index = roi.x + roi.y * ctrl->cfg.roi.width;
    int dqp = ctrl->cfg.roi.dqps[roi_index];
    state->qp = CLIP_TO_QP(state->frame->QP + dqp);
    state->lambda = qp_to_lamba(state, state->qp);
    state->lambda_sqrt = sqrt(state->lambda);

  }
  else if (ctrl->cfg.target_bitrate > 0) {
    lcu_stats_t *lcu         = kvz_get_lcu_stats(state, pos.x, pos.y);
    const uint32_t pixels    = MIN(LCU_WIDTH, state->tile->frame->width  - LCU_WIDTH * pos.x) *
                               MIN(LCU_WIDTH, state->tile->frame->height - LCU_WIDTH * pos.y);

    if (state->frame->num > ctrl->cfg.owf) {
      update_parameters(lcu->bits,
                        pixels,
                        lcu->lambda,
                        &lcu->rc_alpha,
                        &lcu->rc_beta);
    } else {
      lcu->rc_alpha = state->frame->rc_alpha;
      lcu->rc_beta  = state->frame->rc_beta;
    }

    const double target_bits = lcu_allocate_bits(state, pos);
    const double target_bpp  = target_bits / pixels;

    double lambda = clip_lambda(lcu->rc_alpha * pow(target_bpp, lcu->rc_beta));
    // Clip lambda according to the equations 24 and 26 in
    // https://doi.org/10.1109/TIP.2014.2336550
    if (state->frame->num > ctrl->cfg.owf) {
      const double bpp         = lcu->bits / (double)pixels;
      const double lambda_comp = clip_lambda(lcu->rc_alpha * pow(bpp, lcu->rc_beta));
      lambda = CLIP(lambda_comp * 0.7937005259840998,
                    lambda_comp * 1.2599210498948732,
                    lambda);
    }
    lambda = CLIP(state->frame->lambda * 0.6299605249474366,
                  state->frame->lambda * 1.5874010519681994,
                  lambda);
    lambda = clip_lambda(lambda);

    lcu->lambda        = lambda;
    state->lambda      = lambda;
    state->lambda_sqrt = sqrt(lambda);
    state->qp          = lambda_to_qp(lambda);

  } else {
    state->qp          = state->frame->QP;
    state->lambda      = state->frame->lambda;
    state->lambda_sqrt = sqrt(state->frame->lambda);
  }
}
