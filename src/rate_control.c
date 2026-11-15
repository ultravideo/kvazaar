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

#include "rate_control.h"

#include <math.h>

#include "encoder.h"
#include "kvazaar.h"
#include "pthread.h"
#include "lookahead.h"


static const int MIN_SMOOTHING_WINDOW = 40;
static const double MIN_LAMBDA    = 0.1;
static const double MAX_LAMBDA    = 10000;
#define BETA1 1.2517


/**
 * \brief Clip lambda value to a valid range.
 */
static double clip_lambda(double lambda) {
  if (isnan(lambda)) return MAX_LAMBDA;
  return CLIP(MIN_LAMBDA, MAX_LAMBDA, lambda);
}

kvz_rc_data * kvz_get_rc_data(const encoder_control_t * const encoder) {
  kvz_rc_data* data = calloc(1, sizeof(kvz_rc_data));

  if (data == NULL) return NULL;
  if (pthread_mutex_init(&data->ck_frame_lock, NULL) != 0) return NULL;
  if (pthread_mutex_init(&data->lambda_lock, NULL) != 0) return NULL;
  if (pthread_mutex_init(&data->intra_lock, NULL) != 0) return NULL;
  for (int (i) = 0; (i) < KVZ_MAX_GOP_LAYERS; ++(i)) {
    if (pthread_rwlock_init(&data->ck_ctu_lock[i], NULL) != 0) return NULL;
  }

  const int num_lcus = encoder->in.width_in_lcu * encoder->in.height_in_lcu;

  for (int i = 0; i < KVZ_MAX_GOP_LAYERS; i++) {
    data->c_para[i] = malloc(sizeof(double) * num_lcus);
    if (data->c_para[i] == NULL) return NULL;

    data->k_para[i] = malloc(sizeof(double) * num_lcus);
    if (data->k_para[i] == NULL) return NULL;

    data->pic_c_para[i] = 5.0;
    data->pic_k_para[i] = -0.1;

    for (int j = 0; j < num_lcus; j++) {
      data->c_para[i][j] = 5.0;
      data->k_para[i][j] = -0.1;
    }
  }
  data->intra_bpp = calloc(num_lcus, sizeof(double));
  if (data->intra_bpp == NULL) return NULL;
  data->intra_dis = calloc(num_lcus, sizeof(double));
  if (data->intra_dis == NULL) return NULL;

  memset(data->previous_lambdas, 0, sizeof(data->previous_lambdas));

  data->previous_frame_lambda = 0.0;

  data->intra_pic_bpp = 0.0;
  data->intra_pic_distortion = 0.0;

  data->intra_alpha = 6.7542000000000000;
  data->intra_beta = 1.7860000000000000;
  if(encoder->cfg.stats_file_prefix) {
    char buff[128];
    sprintf(buff, "%sbits.txt", encoder->cfg.stats_file_prefix);
    data->bits_file = fopen(buff, "w");
    sprintf(buff, "%sdist.txt", encoder->cfg.stats_file_prefix);
    data->dist_file = fopen(buff, "w");
    sprintf(buff, "%sqp.txt", encoder->cfg.stats_file_prefix);
    data->qp_file = fopen(buff, "w");
    sprintf(buff, "%slambda.txt", encoder->cfg.stats_file_prefix);
    data->lambda_file = fopen(buff, "w");
  }

  data->smoothing_window = MIN_SMOOTHING_WINDOW;

  return data;
}

void kvz_free_rc_data(kvz_rc_data *data) {
  if (data == NULL) return;

  pthread_mutex_destroy(&data->ck_frame_lock);
  pthread_mutex_destroy(&data->lambda_lock);
  pthread_mutex_destroy(&data->intra_lock);
  for (int i = 0; i < KVZ_MAX_GOP_LAYERS; ++i) {
    pthread_rwlock_destroy(&data->ck_ctu_lock[i]);
  }

  if (data->intra_bpp) FREE_POINTER(data->intra_bpp);
  if (data->intra_dis) FREE_POINTER(data->intra_dis);
  for (int i = 0; i < KVZ_MAX_GOP_LAYERS; i++) {
    if (data->c_para[i]) FREE_POINTER(data->c_para[i]);
    if (data->k_para[i]) FREE_POINTER(data->k_para[i]);
  }
  FREE_POINTER(data);
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
  int pictures_coded = 0;
  if(encoder->cfg.gop_len) {
    pictures_coded = MAX(0, state->frame->num - CEILDIV(encoder->cfg.owf, encoder->cfg.gop_len)*encoder->cfg.gop_len);
  }
  else {
    pictures_coded = MAX(0, state->frame->num - encoder->cfg.owf);
  }
    
  if (encoder->cfg.gop_len > 0 && encoder->cfg.owf > 0) {
    // Subtract number of bits in the partially coded GOP.
    bits_coded -= state->frame->cur_gop_bits_coded;
  }

  state->frame->new_ratecontrol->smoothing_window = MAX(MIN_SMOOTHING_WINDOW, state->frame->new_ratecontrol->smoothing_window - MAX(encoder->cfg.gop_len / 2, 1));
  double gop_target_bits = -1;

  while( gop_target_bits < 0 && state->frame->new_ratecontrol->smoothing_window < 150) {
    // Equation 12 from https://doi.org/10.1109/TIP.2014.2336550
    gop_target_bits =
      (encoder->target_avg_bppic * (pictures_coded + state->frame->new_ratecontrol->smoothing_window) - bits_coded)
      * MAX(1, encoder->cfg.gop_len) / state->frame->new_ratecontrol->smoothing_window;
    if(gop_target_bits < 0) {
      state->frame->new_ratecontrol->smoothing_window += 10;
    }
  }
  // Allocate at least 200 bits for each GOP like HM does.
  return MAX(200, gop_target_bits);
}

static int xCalcHADs8x8_ISlice(kvz_pixel * piOrg, int y, int iStrideOrg)
{
  piOrg += y * iStrideOrg;
  int i, j;
  int diff[64], m1[8][8], m2[8][8], m3[8][8], iSumHad = 0;

  for (int k = 0; k < 64; k += 8) {
    diff[k + 0] = piOrg[0];
    diff[k + 1] = piOrg[1];
    diff[k + 2] = piOrg[2];
    diff[k + 3] = piOrg[3];
    diff[k + 4] = piOrg[4];
    diff[k + 5] = piOrg[5];
    diff[k + 6] = piOrg[6];
    diff[k + 7] = piOrg[7];

    piOrg += iStrideOrg;
  }

  //horizontal
  for (j = 0; j < 8; j++) {
    int jj = j << 3;
    m2[j][0] = diff[jj] + diff[jj + 4];
    m2[j][1] = diff[jj + 1] + diff[jj + 5];
    m2[j][2] = diff[jj + 2] + diff[jj + 6];
    m2[j][3] = diff[jj + 3] + diff[jj + 7];
    m2[j][4] = diff[jj] - diff[jj + 4];
    m2[j][5] = diff[jj + 1] - diff[jj + 5];
    m2[j][6] = diff[jj + 2] - diff[jj + 6];
    m2[j][7] = diff[jj + 3] - diff[jj + 7];

    m1[j][0] = m2[j][0] + m2[j][2];
    m1[j][1] = m2[j][1] + m2[j][3];
    m1[j][2] = m2[j][0] - m2[j][2];
    m1[j][3] = m2[j][1] - m2[j][3];
    m1[j][4] = m2[j][4] + m2[j][6];
    m1[j][5] = m2[j][5] + m2[j][7];
    m1[j][6] = m2[j][4] - m2[j][6];
    m1[j][7] = m2[j][5] - m2[j][7];

    m2[j][0] = m1[j][0] + m1[j][1];
    m2[j][1] = m1[j][0] - m1[j][1];
    m2[j][2] = m1[j][2] + m1[j][3];
    m2[j][3] = m1[j][2] - m1[j][3];
    m2[j][4] = m1[j][4] + m1[j][5];
    m2[j][5] = m1[j][4] - m1[j][5];
    m2[j][6] = m1[j][6] + m1[j][7];
    m2[j][7] = m1[j][6] - m1[j][7];
  }

  //vertical
  for (i = 0; i < 8; i++) {
    m3[0][i] = m2[0][i] + m2[4][i];
    m3[1][i] = m2[1][i] + m2[5][i];
    m3[2][i] = m2[2][i] + m2[6][i];
    m3[3][i] = m2[3][i] + m2[7][i];
    m3[4][i] = m2[0][i] - m2[4][i];
    m3[5][i] = m2[1][i] - m2[5][i];
    m3[6][i] = m2[2][i] - m2[6][i];
    m3[7][i] = m2[3][i] - m2[7][i];

    m1[0][i] = m3[0][i] + m3[2][i];
    m1[1][i] = m3[1][i] + m3[3][i];
    m1[2][i] = m3[0][i] - m3[2][i];
    m1[3][i] = m3[1][i] - m3[3][i];
    m1[4][i] = m3[4][i] + m3[6][i];
    m1[5][i] = m3[5][i] + m3[7][i];
    m1[6][i] = m3[4][i] - m3[6][i];
    m1[7][i] = m3[5][i] - m3[7][i];

    m2[0][i] = m1[0][i] + m1[1][i];
    m2[1][i] = m1[0][i] - m1[1][i];
    m2[2][i] = m1[2][i] + m1[3][i];
    m2[3][i] = m1[2][i] - m1[3][i];
    m2[4][i] = m1[4][i] + m1[5][i];
    m2[5][i] = m1[4][i] - m1[5][i];
    m2[6][i] = m1[6][i] + m1[7][i];
    m2[7][i] = m1[6][i] - m1[7][i];
  }

  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      iSumHad += abs(m2[i][j]);
    }
  }
  iSumHad -= abs(m2[0][0]);
  iSumHad = (iSumHad + 2) >> 2;
  return(iSumHad);
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
    // A new GOP starts at this frame. cur_gop_bits_coded is reset when the
    // bitstream of this frame is written.
    state->frame->cur_gop_target_bits = gop_allocate_bits(state);
  } else {
    state->frame->cur_gop_target_bits =
      state->previous_encoder_state->frame->cur_gop_target_bits;
  }

  if (state->frame->is_irap && encoder->cfg.intra_bit_allocation) {
    int total_cost = 0;
    for (int y = 0; y < encoder->cfg.height; y += 8) {
      for (int x = 0; x < encoder->cfg.width; x += 8) {
        int cost = xCalcHADs8x8_ISlice(state->tile->frame->source->y + x, y, state->tile->frame->source->stride);
        total_cost += cost;
        kvz_get_lcu_stats(state, x / 64, y / 64)->i_cost += cost;
      }
    }
    state->frame->icost = total_cost;
    state->frame->remaining_weight = total_cost;

    double bits = state->frame->cur_gop_target_bits / MAX(encoder->cfg.gop_len, 1);
    double alpha, beta = 0.5582;
    if (bits * 40 < encoder->cfg.width * encoder->cfg.height) {
      alpha = 0.25;
    }
    else {
      alpha = 0.3;
    }

    double low_limit = (encoder->cfg.rc_algorithm == KVZ_LAMBDA && encoder->cfg.rdo < 2 ? 1.0 : 1.2) * pow(state->encoder_control->cfg.framerate, -0.873) * state->encoder_control->cfg.target_bitrate;
    double high_limit = (encoder->cfg.rdo < 2 ? (encoder->cfg.rc_algorithm == KVZ_LAMBDA ? 1.1 : 3.5) : 2.25) * pow(state->encoder_control->cfg.framerate, -0.61) * state->encoder_control->cfg.target_bitrate;
    double original_bits = alpha * pow(state->frame->icost * 4 / bits, beta) * bits;
    double limited = MIN(MAX(low_limit, original_bits), high_limit);
    return limited;
  }

  if (encoder->cfg.gop_len <= 0) {
    return state->frame->cur_gop_target_bits;
  }

  const double pic_weight = encoder->gop_layer_weights[
    encoder->cfg.gop[state->frame->gop_offset].layer - 1];
  const double pic_target_bits =
    state->frame->cur_gop_target_bits * pic_weight - pic_header_bits(state);
  // Allocate at least 100 bits for each picture like HM does.
  const double intra_bits = state->frame->is_irap ? (encoder->cfg.rdo < 2 ? 4 : 6) * state->encoder_control->target_avg_bppic : pic_target_bits;
  return MAX(MAX(100, pic_target_bits), intra_bits);
}

static int8_t lambda_to_qp(const double lambda)
{
  const int8_t qp = 4.2005 * log(lambda) + 13.7223 + 0.5;
  return CLIP_TO_QP(qp);
}

static double solve_cubic_equation(const encoder_state_config_frame_t * const state,
                            int ctu_index,
                            int last_ctu,
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
    assert(!((state->c_para[i] <= 0) || (state->k_para[i] >= 0))); //Check C and K during each solution 

    double CLCU = state->c_para[i];
    double KLCU = state->k_para[i];
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

  if (delta > 0.0)  //Check whether delta is right
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
    best_lambda = est_lambda;  //Use the original picture estimated lambda for the current CTU
  }
  best_lambda = CLIP(0.001, 100000000.0, best_lambda);

  return best_lambda;
}

static INLINE double calculate_weights(encoder_state_t* const state, const int ctu_count, double est_lambda) {
  double total_weight = 0;
  for(int i = 0; i < ctu_count; i++) {
    double c_lcu = state->frame->c_para[i];
    double k_lcu = state->frame->k_para[i];
    double a = -c_lcu * k_lcu / pow(state->frame->lcu_stats[i].pixels, k_lcu - 1.0);
    double b = -1.0 / (k_lcu - 1.0);
    state->frame->lcu_stats[i].original_weight = state->frame->lcu_stats[i].weight = pow(a / est_lambda, b);
    if (state->frame->lcu_stats[i].weight < 0.01) {
      state->frame->lcu_stats[i].weight = 0.01;
    }
    total_weight += state->frame->lcu_stats[i].weight;
  }
  return total_weight;
}


void kvz_estimate_pic_lambda(encoder_state_t * const state) {
  const encoder_control_t * const encoder = state->encoder_control;

  const int layer = encoder->cfg.gop[state->frame->gop_offset].layer - (state->frame->is_irap ? 1 : 0);
  const int ctu_count = state->tile->frame->height_in_lcu * state->tile->frame->width_in_lcu;

  double alpha;
  double beta;
  if(state->frame->is_irap && encoder->cfg.intra_bit_allocation) {
    pthread_mutex_lock(&state->frame->new_ratecontrol->intra_lock);
    alpha = state->frame->new_ratecontrol->intra_alpha;
    beta = state->frame->new_ratecontrol->intra_beta;
    pthread_mutex_unlock(&state->frame->new_ratecontrol->intra_lock);
  }
  else if(state->frame->poc == 0) {
    alpha = state->frame->rc_alpha;
    beta = state->frame->rc_beta;
  }
  else {
    pthread_mutex_lock(&state->frame->new_ratecontrol->ck_frame_lock);
    alpha = -state->frame->new_ratecontrol->pic_c_para[layer] *
      state->frame->new_ratecontrol->pic_k_para[layer];
    beta = state->frame->new_ratecontrol->pic_k_para[layer] - 1;
    pthread_mutex_unlock(&state->frame->new_ratecontrol->ck_frame_lock);
  }
  double bits = pic_allocate_bits(state);
  state->frame->cur_pic_target_bits = bits;

  double est_lambda;
  int32_t num_pixels = state->encoder_control->cfg.width * state->encoder_control->cfg.height;
  double bpp = bits / num_pixels;
  if (state->frame->is_irap) {
    if(encoder->cfg.intra_bit_allocation) {
      state->frame->i_bits_left = bits;
      double temp = pow(state->frame->icost / num_pixels, BETA1);
      est_lambda = alpha / 256 * pow(temp/bpp, beta);
    }
    else {
      // arbitrary reduction to the lambda for intra frames
      est_lambda = alpha * pow(bpp, beta) * 0.5;
    }
  }
  else {
    est_lambda = alpha * pow(bpp, beta);
  }

  double temp_lambda;
  pthread_mutex_lock(&state->frame->new_ratecontrol->lambda_lock);
  if ((temp_lambda = state->frame->new_ratecontrol->previous_lambdas[layer]) > 0.0) {
    temp_lambda = CLIP(0.1, 10000.0, temp_lambda);
    est_lambda = CLIP(temp_lambda * pow(2.0, -1), temp_lambda * 2, est_lambda);
  }

  if((temp_lambda = state->frame->new_ratecontrol->previous_frame_lambda) > 0.0) {
    temp_lambda = CLIP(0.1, 2000.0, temp_lambda);
    est_lambda = CLIP(temp_lambda * pow(2.0, -10.0 / 3.0), temp_lambda * pow(2.0, 10.0 / 3.0), est_lambda);
  }
  pthread_mutex_unlock(&state->frame->new_ratecontrol->lambda_lock);

  est_lambda = CLIP(0.1, 10000.0, est_lambda);

  double total_weight = 0;

  if(!state->frame->is_irap) {
    double best_lambda = est_lambda;
    if(!state->encoder_control->cfg.frame_allocation) {
      pthread_rwlock_rdlock(&state->frame->new_ratecontrol->ck_ctu_lock[layer]);
      memcpy(state->frame->c_para, state->frame->new_ratecontrol->c_para[layer], ctu_count * sizeof(double));
      memcpy(state->frame->k_para, state->frame->new_ratecontrol->k_para[layer], ctu_count * sizeof(double));
      pthread_rwlock_unlock(&state->frame->new_ratecontrol->ck_ctu_lock[layer]);
      temp_lambda = est_lambda;
      double taylor_e3;
      int iteration_number = 0;
      do {
        taylor_e3 = 0.0;
        best_lambda = temp_lambda = solve_cubic_equation(state->frame, 0, ctu_count, temp_lambda, bits);
        for (int i = 0; i < ctu_count; ++i) {
          double CLCU = state->frame->c_para[i];
          double KLCU = state->frame->k_para[i];
          double a = -CLCU * KLCU / pow(state->frame->lcu_stats[i].pixels, KLCU - 1.0);
          double b = -1.0 / (KLCU - 1.0);
          taylor_e3 += pow(a / best_lambda, b);
        }
        iteration_number++;
      }
      while (fabs(taylor_e3 - bits) > 0.01 && iteration_number <= 11);
    }
    total_weight = calculate_weights(state, ctu_count, best_lambda);
    state->frame->remaining_weight = bits;
  }
  else {
    for (int i = 0; i < ctu_count; ++i) {
      state->frame->lcu_stats[i].weight = MAX(0.01,
        state->frame->lcu_stats[i].pixels * pow(est_lambda / alpha,
                                                1.0 / beta));
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
  
  int num_ctu = state->encoder_control->in.width_in_lcu * state->encoder_control->in.height_in_lcu;
  const int index = pos.x + pos.y * state->tile->frame->width_in_lcu;

  if (state->frame->is_irap) {
    if(encoder->cfg.intra_bit_allocation) {
      int cus_left = num_ctu - index + 1;
      int window = MIN(4, cus_left);
      double mad = kvz_get_lcu_stats(state, pos.x, pos.y)->i_cost;

      pthread_mutex_lock(&state->frame->rc_lock);
      double bits_left = state->frame->cur_pic_target_bits - state->frame->cur_frame_bits_coded;
      double weighted_bits_left = (bits_left * window + (bits_left - state->frame->i_bits_left)*cus_left) / window;
      avg_bits = mad * weighted_bits_left / state->frame->remaining_weight;
      state->frame->remaining_weight -= mad;
      state->frame->i_bits_left -= state->frame->cur_pic_target_bits * mad / state->frame->icost;
      pthread_mutex_unlock(&state->frame->rc_lock);
    }
    else {
      avg_bits = state->frame->cur_pic_target_bits * ((double)state->frame->lcu_stats[index].pixels /
        (state->encoder_control->in.height * state->encoder_control->in.width));
    }
  }
  else {
    double total_weight = 0;
    // In case wpp is used only the ctus of the current frame are safe to use
    const int used_ctu_count = MIN(4, (encoder->cfg.wpp ? (pos.y + 1) * encoder->in.width_in_lcu : num_ctu) - index);
    int target_bits = 0;
    double best_lambda = 0.0;
    double temp_lambda = state->frame->lambda;
    double taylor_e3 = 0.0;
    int iter = 0;

    int last_ctu = index + used_ctu_count;
    for (int i = index; i < last_ctu; i++) {
      target_bits += state->frame->lcu_stats[i].weight;
    }

    pthread_mutex_lock(&state->frame->rc_lock);
    total_weight = state->frame->remaining_weight;
    target_bits = MAX(target_bits + state->frame->cur_pic_target_bits - state->frame->cur_frame_bits_coded - (int)total_weight, 10);
    pthread_mutex_unlock(&state->frame->rc_lock);

    //just similar with the process at frame level, details can refer to the function kvz_estimate_pic_lambda
    do {
      taylor_e3 = 0.0;
      best_lambda = solve_cubic_equation(state->frame, index, last_ctu, temp_lambda, target_bits);
      temp_lambda = best_lambda;
      for (int i = index; i < last_ctu; i++) {
        double CLCU = state->frame->c_para[i];
        double KLCU = state->frame->k_para[i];
        double a = -CLCU * KLCU / pow((double)state->frame->lcu_stats[i].pixels, KLCU - 1.0);
        double b = -1.0 / (KLCU - 1.0);
        taylor_e3 += pow(a / best_lambda, b);
      }
      iter++;
    } while (fabs(taylor_e3 - target_bits) > 0.01 && iter < 5);

    double c_ctu = state->frame->c_para[index];
    double k_ctu = state->frame->k_para[index];
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

static double qp_to_lambda(encoder_state_t* const state, int qp)
{
  const int shift_qp = 12;
  double lambda = 0.57 * pow(2.0, (qp - shift_qp) / 3.0);

  // NOTE: HM adjusts lambda for inter according to Hadamard usage in ME.
  //       SATD is currently always enabled for ME, so this has no effect.
  // bool hadamard_me = true;
  // if (!hadamard_me && state->frame->slicetype != KVZ_SLICE_I) {
  //   lambda *= 0.95;
  // }

  return lambda;
}

 void kvz_set_ctu_qp_lambda(encoder_state_t * const state, vector2d_t pos) {
  double bits = get_ctu_bits(state, pos);

  const encoder_control_t * const encoder = state->encoder_control;
  const int frame_allocation = state->encoder_control->cfg.frame_allocation;

  int index = pos.x + pos.y * state->encoder_control->in.width_in_lcu;
  lcu_stats_t* ctu = &state->frame->lcu_stats[index];
  double bpp = bits / ctu->pixels;

  double alpha;
  double beta;
  if (state->frame->is_irap && encoder->cfg.intra_bit_allocation) {
    pthread_mutex_lock(&state->frame->new_ratecontrol->intra_lock);
    alpha = state->frame->new_ratecontrol->intra_alpha;
    beta = state->frame->new_ratecontrol->intra_beta;
    pthread_mutex_unlock(&state->frame->new_ratecontrol->intra_lock);
  }
  else if(state->frame->num == 0) {
    alpha = state->frame->rc_alpha;
    beta = state->frame->rc_beta;
  }
  else {
    alpha = -state->frame->c_para[index] * state->frame->k_para[index];
    beta = state->frame->k_para[index] - 1;
  }
  alpha *= (double)(1 << (KVZ_BIT_DEPTH - 8));
  double est_lambda;
  int est_qp;
  if (state->frame->is_irap && encoder->cfg.intra_bit_allocation) {
    double cost_per_pixel = (double)ctu->i_cost / ctu->pixels;
    cost_per_pixel = pow(cost_per_pixel, BETA1);
    est_lambda = alpha / 256.0 * pow(cost_per_pixel / bpp, beta);
    est_qp = state->frame->QP;
    double max_lambda = exp(((double)est_qp + 2.49 - 13.7122) / 4.2005);
    double min_lambda = exp(((double)est_qp - 2.49 - 13.7122) / 4.2005);
    est_lambda = CLIP(min_lambda, max_lambda, est_lambda);

    est_qp = lambda_to_qp(est_lambda);
  }
  else {
    // In case wpp is used the previous ctus may not be ready from above rows
    const int ctu_limit = encoder->cfg.wpp ? pos.y * encoder->in.width_in_lcu : 0;
    
    est_lambda = alpha * pow(bpp, beta) * (state->frame->is_irap ? 0.5 : 1);
    const double clip_lambda = state->frame->lambda;

    double clip_neighbor_lambda = -1;
    int clip_qp = -1;
    if (encoder->cfg.clip_neighbour || state->frame->num == 0) {
      for (int temp_index = index - 1; temp_index >= ctu_limit; --temp_index) {
        if (state->frame->lcu_stats[temp_index].lambda > 0) {
          clip_neighbor_lambda = state->frame->lcu_stats[temp_index].lambda;
          break;
        }
      }
      for (int temp_index = index - 1; temp_index >= ctu_limit; --temp_index) {
        if (state->frame->lcu_stats[temp_index].qp > -1) {
          clip_qp = state->frame->lcu_stats[temp_index].qp;
          break;
        }
      }
    }
    else {
      
      if (state->frame->lcu_stats[index].lambda > 0) {
        clip_neighbor_lambda = state->frame->previous_layer_state->frame->lcu_stats[index].lambda;
      }
      if (state->frame->lcu_stats[index].qp > 0) {
        clip_qp = state->frame->previous_layer_state->frame->lcu_stats[index].qp;
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

    est_qp = lambda_to_qp(est_lambda);

    if( clip_qp > -1) {
      est_qp = CLIP(clip_qp - 1 - frame_allocation,
        clip_qp + 1 + frame_allocation,
        est_qp);
    }

    est_qp = CLIP(state->frame->QP - 2 - frame_allocation,
      state->frame->QP + 2 + frame_allocation,
      est_qp);
  }

  state->lambda = est_lambda;
  state->lambda_sqrt = sqrt(est_lambda);
  state->qp = est_qp;
  ctu->qp = est_qp;
  ctu->lambda = est_lambda;
  ctu->i_cost = 0;

  // Apply variance adaptive quantization
  if (encoder->cfg.vaq || encoder->cfg.cutree > 0) {
    vector2d_t lcu = {
      pos.x + state->tile->lcu_offset_x,
      pos.y + state->tile->lcu_offset_y
    };
    int id = lcu.x + lcu.y * state->tile->frame->width_in_lcu;
    int aq_offset = round(state->frame->aq_offsets[id]);
    state->qp += aq_offset;
    // Maximum delta QP is clipped according to ITU T-REC-H.265 specification chapter 7.4.9.10 Transform unit semantics
    // Clipping range is a function of bit depth
    // Since this value will be later combined with qp_pred, clip to half of that instead to be safe
    state->qp = CLIP(state->frame->QP + KVZ_QP_DELTA_MIN / 2, state->frame->QP + KVZ_QP_DELTA_MAX / 2, state->qp);
    state->qp = CLIP_TO_QP(state->qp);
    state->lambda = qp_to_lambda(state, state->qp);
    state->lambda_sqrt = sqrt(state->lambda);

    ctu->adjust_lambda = state->lambda;
    ctu->adjust_qp = state->qp;
    //ctu->qp = state->qp;
    //ctu->lambda = state->lambda;
  }
}


static void update_pic_ck(encoder_state_t * const state, double bpp, double distortion, double lambda, int layer) {
  double new_k = 0, new_c;
  if(state->frame->num == 1) {
    new_k = log(distortion / state->frame->new_ratecontrol->intra_pic_distortion) /
      log(bpp / state->frame->new_ratecontrol->intra_pic_bpp);
    new_c = distortion / pow(bpp, new_k);
  }
  new_k = -bpp * lambda / distortion;
  new_c = distortion / pow(bpp, new_k);

  new_c = CLIP(+.1, 100.0, new_c);
  new_k = CLIP(-3.0, -0.001, new_k);

  // The model is read by kvz_estimate_pic_lambda of a later frame, which can
  // run in the main thread at the same time as this (bitstream writing of an
  // earlier frame in a worker), so update it under the same lock.
  pthread_mutex_lock(&state->frame->new_ratecontrol->ck_frame_lock);
  if(state->frame->is_irap || state->frame->num <= (4 - state->encoder_control->cfg.frame_allocation)) {
    for(int i = 1; i < 5; i++) {
      state->frame->new_ratecontrol->pic_c_para[i] = new_c;
      state->frame->new_ratecontrol->pic_k_para[i] = new_k;
    }
  }
  else {
    state->frame->new_ratecontrol->pic_c_para[layer] = new_c;
    state->frame->new_ratecontrol->pic_k_para[layer] = new_k;
  }
  pthread_mutex_unlock(&state->frame->new_ratecontrol->ck_frame_lock);
}


static void update_ck(encoder_state_t * const state, int ctu_index, int layer)
{
  double bpp = (double)state->frame->lcu_stats[ctu_index].bits / state->frame->lcu_stats[ctu_index].pixels;
  double distortion = state->frame->lcu_stats[ctu_index].distortion;
  double lambda = state->frame->lcu_stats[ctu_index].lambda;

  double new_k = 0, new_c = -1;
  if (!state->frame->lcu_stats[ctu_index].skipped) {
    distortion = MAX(distortion, 0.0001);

    bpp = CLIP(0.0001, 10.0, bpp);
    new_k = -bpp * lambda / distortion;
    new_k = CLIP(-3.0, -0.001, new_k);
    new_c = distortion / pow(bpp, new_k);
    
    new_c = CLIP(+.1, 100.0, new_c);

    if (state->frame->is_irap || state->frame->num <= (4 - state->encoder_control->cfg.frame_allocation)) {
      for (int i = 1; i < 5; i++) {
        state->frame->new_ratecontrol->c_para[i][ctu_index] = new_c;
        state->frame->new_ratecontrol->k_para[i][ctu_index] = new_k;
      }
    }
    else {
      state->frame->new_ratecontrol->c_para[layer][ctu_index] = new_c;
      state->frame->new_ratecontrol->k_para[layer][ctu_index] = new_k;
    }
  }
}


static int calc_poc(encoder_state_t * const state) {
  const encoder_control_t * const encoder = state->encoder_control;
  if((encoder->cfg.open_gop && !encoder->cfg.gop_lowdelay) || !encoder->cfg.intra_period) {
    return state->frame->poc;
  }
  if(!encoder->cfg.gop_len || encoder->cfg.open_gop || encoder->cfg.intra_period == 1 || encoder->cfg.gop_lowdelay) {
    return state->frame->poc + state->frame->num / encoder->cfg.intra_period * encoder->cfg.intra_period;
  }
  if (!encoder->cfg.gop_lowdelay && !encoder->cfg.open_gop) {
    return state->frame->poc + state->frame->num / (encoder->cfg.intra_period + 1) * (encoder->cfg.intra_period + 1);
  }
  assert(0);
  return -1;
}


void kvz_update_after_picture(encoder_state_t * const state) {
  double total_distortion = 0;
  double lambda = 0;
  int32_t pixels = (state->encoder_control->in.width * state->encoder_control->in.height);
  double pic_bpp = (double)state->frame->cur_frame_bits_coded / pixels;

  const encoder_control_t * const encoder = state->encoder_control;
  const int layer = encoder->cfg.gop[state->frame->gop_offset].layer - (state->frame->is_irap ? 1 : 0);

  if (state->frame->is_irap && encoder->cfg.intra_bit_allocation) {
    double lnbpp = log(pow(state->frame->icost / pixels, BETA1));
    pthread_mutex_lock(&state->frame->new_ratecontrol->intra_lock);
    double diff_lambda = state->frame->new_ratecontrol->intra_beta * log(state->frame->cur_frame_bits_coded) - log(state->frame->cur_pic_target_bits);

    diff_lambda = CLIP(-0.125, 0.125, 0.25*diff_lambda);

    state->frame->new_ratecontrol->intra_alpha *= exp(diff_lambda);
    state->frame->new_ratecontrol->intra_beta += diff_lambda / lnbpp;
    pthread_mutex_unlock(&state->frame->new_ratecontrol->intra_lock);
  }

  if (encoder->cfg.stats_file_prefix) {
    int poc = calc_poc(state);
    fprintf(state->frame->new_ratecontrol->dist_file, "%d %d %d\n", poc, encoder->in.width_in_lcu, encoder->in.height_in_lcu);
    fprintf(state->frame->new_ratecontrol->bits_file, "%d %d %d\n", poc, encoder->in.width_in_lcu, encoder->in.height_in_lcu);
    fprintf(state->frame->new_ratecontrol->qp_file, "%d %d %d\n", poc, encoder->in.width_in_lcu, encoder->in.height_in_lcu);
    fprintf(state->frame->new_ratecontrol->lambda_file, "%d %d %d\n", poc, encoder->in.width_in_lcu, encoder->in.height_in_lcu);
  }

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
      if(encoder->cfg.stats_file_prefix) {
        fprintf(state->frame->new_ratecontrol->dist_file, "%f ", ctu->distortion);
        fprintf(state->frame->new_ratecontrol->bits_file, "%d ", ctu->bits);
        fprintf(state->frame->new_ratecontrol->qp_file, "%d ", ctu->adjust_qp ? ctu->adjust_qp : ctu->qp);
        fprintf(state->frame->new_ratecontrol->lambda_file, "%f ", ctu->adjust_lambda ? ctu->adjust_lambda : ctu->lambda);
      }
    }
    if (encoder->cfg.stats_file_prefix) {
      fprintf(state->frame->new_ratecontrol->dist_file, "\n");
      fprintf(state->frame->new_ratecontrol->bits_file, "\n");
      fprintf(state->frame->new_ratecontrol->qp_file, "\n");
      fprintf(state->frame->new_ratecontrol->lambda_file, "\n");
    }
  }

  if(encoder->cfg.stats_file_prefix && encoder->cfg.rc_algorithm != KVZ_OBA) return;

  total_distortion /= (state->encoder_control->in.height_in_lcu * state->encoder_control->in.width_in_lcu);
  if (state->frame->is_irap) {
    pthread_mutex_lock(&state->frame->new_ratecontrol->intra_lock);
    for (int y_ctu = 0; y_ctu < state->encoder_control->in.height_in_lcu; y_ctu++) {
      for (int x_ctu = 0; x_ctu < state->encoder_control->in.width_in_lcu; x_ctu++) {
        lcu_stats_t *ctu = kvz_get_lcu_stats(state, x_ctu, y_ctu);
        state->frame->new_ratecontrol->intra_dis[x_ctu + y_ctu * state->encoder_control->in.width_in_lcu] =
          ctu->distortion;
        state->frame->new_ratecontrol->intra_bpp[x_ctu + y_ctu * state->encoder_control->in.width_in_lcu] =
          ctu->bits / ctu->pixels;
      }
    }
    state->frame->new_ratecontrol->intra_pic_distortion = total_distortion;
    state->frame->new_ratecontrol->intra_pic_bpp = pic_bpp;
    pthread_mutex_unlock(&state->frame->new_ratecontrol->intra_lock);
  }

  pthread_mutex_lock(&state->frame->new_ratecontrol->lambda_lock);
  state->frame->new_ratecontrol->previous_frame_lambda = lambda;
  state->frame->new_ratecontrol->previous_lambdas[layer] = lambda;
  pthread_mutex_unlock(&state->frame->new_ratecontrol->lambda_lock);

  update_pic_ck(state, pic_bpp, total_distortion, lambda, layer);
  if (state->frame->num <= 4 || state->frame->is_irap){
    for (int i = 1; i < 5; ++i) {
      pthread_rwlock_wrlock(&state->frame->new_ratecontrol->ck_ctu_lock[i]);
    }
  }
  else{
    pthread_rwlock_wrlock(&state->frame->new_ratecontrol->ck_ctu_lock[layer]);
  }
  for(int i = 0; i < state->encoder_control->in.width_in_lcu * state->encoder_control->in.height_in_lcu; i++) {
    update_ck(state, i, layer);
  }
  if (state->frame->num <= 4 || state->frame->is_irap){
    for (int i = 1; i < 5; ++i) {
      pthread_rwlock_unlock(&state->frame->new_ratecontrol->ck_ctu_lock[i]);
    }
  }
  else{
    pthread_rwlock_unlock(&state->frame->new_ratecontrol->ck_ctu_lock[layer]);
  }
}

/**
 * \brief Allocate bits and set lambda and QP for the current picture.
 * \param state the main encoder state
 */
/**
 * \brief Frame QP from a base QP using the GOP QP offsets as in CQP mode.
 */
static int8_t gop_frame_qp(const encoder_state_t * const state, double base_qp)
{
  const encoder_control_t * const ctrl = state->encoder_control;
  kvz_gop_config const * const gop = &ctrl->cfg.gop[state->frame->gop_offset];

  if (ctrl->cfg.gop_len > 0 && state->frame->slicetype != KVZ_SLICE_I) {
    double qp = base_qp;
    qp += gop->qp_offset;
    qp += CLIP(0.0, 3.0, qp * gop->qp_model_scale + gop->qp_model_offset);
    return CLIP_TO_QP((int)(qp + 0.5));
  }
  return CLIP_TO_QP((int)(base_qp + ctrl->cfg.intra_qp_offset + 0.5));
}


static INLINE double qp_to_qscale(double qp)
{
  return 0.85 * pow(2.0, (qp - 12.0) / 6.0);
}


static INLINE double qscale_to_qp(double qscale)
{
  return 12.0 + 6.0 * log2(qscale / 0.85);
}


/** \brief Prior of the ABR bits model for inter pictures, per pixel. */
#define ABR_PRIOR_K_INTER 0.0015
/** \brief Prior of the ABR bits model for intra pictures, per pixel. */
#define ABR_PRIOR_K_INTRA 0.005
/** \brief Adaptation speed of the ABR bits model. */
#define ABR_K_ALPHA 0.3
/** \brief Time in seconds over which the bitrate error is corrected. */
#define ABR_ERROR_SECONDS 2.0
/** \brief Maximum change of the base QP between GOPs. */
#define ABR_MAX_QP_STEP 2.0
/** \brief Exponent of the quantizer scale in the bits model. */
#define ABR_GAMMA 1.3
/** \brief Typical base QP at ABR_BPP_REF bits per pixel. */
#define ABR_BPP_QP 27.0
#define ABR_BPP_REF 0.02
/** \brief The first picture is not coded more than this below the QP
 *         estimated from the bits per pixel. */
#define ABR_FIRST_QP_MARGIN 4.0


static double abr_class_k(const encoder_control_t *const ctrl, const kvz_rc_data *rc, int cls)
{
  const int intra_class = KVZ_MAX_GOP_LENGTH;
  if (rc->abr_k_count[cls] > 0) return rc->abr_k[cls];
  if (cls == intra_class) return ABR_PRIOR_K_INTRA * ctrl->in.pixels_per_pic;

  // Use the average of the observed inter classes.
  double sum = 0.0;
  int n = 0;
  for (int i = 0; i < MAX(1, ctrl->cfg.gop_len); i++) {
    if (rc->abr_k_count[i] > 0) {
      sum += rc->abr_k[i];
      n++;
    }
  }
  return n ? sum / n : ABR_PRIOR_K_INTER * ctrl->in.pixels_per_pic;
}


/**
 * \brief Set the frame QP for the lookahead based average bitrate control.
 *
 * One pass ABR in the spirit of x264. The bits of a frame are modelled as
 *   bits = k * complexity / qscale^gamma
 * with a separate k for each GOP position and for intra pictures. The base
 * quantizer scale of a frame is
 *   qscale = complexity^(1 - qcomp) / rate_factor
 * where the rate factor is chosen so that the frames from the current GOP to
 * the end of the lookahead would hit the target bitrate, after correcting
 * the bitrate error accumulated so far. The GOP QP offsets are applied on
 * top of the base QP like in the CQP mode.
 */
static double abr_class_k(const encoder_control_t *const ctrl, const kvz_rc_data *rc, int cls);


static INLINE double get_framerate(const encoder_control_t *const ctrl)
{
  return ctrl->cfg.framerate_num != 0 ?
    ctrl->cfg.framerate_num / (double)ctrl->cfg.framerate_denom : ctrl->cfg.framerate;
}


/**
 * \brief Update the bits model with the frame previously coded by the
 * encoder state. Its bitstream has been written since the state is being
 * reused, and the states are reused in coding order.
 */
static void abr_account_previous(encoder_state_t * const state)
{
  const encoder_control_t * const ctrl = state->encoder_control;
  kvz_rc_data *rc = state->frame->new_ratecontrol;

  if (ctrl->cfg.vbv_maxrate > 0 && !rc->vbv_init_done) {
    rc->vbv_fill = ctrl->cfg.vbv_init * ctrl->cfg.vbv_bufsize;
    rc->vbv_init_done = true;
  }

  if (!state->frame->abr_pending) return;

  const double bits = (double)state->stats_bitstream_length * 8;
  const int cls = state->frame->abr_class;
  const double k = bits * pow(state->frame->abr_qscale, ABR_GAMMA) / state->frame->abr_cost;
  if (rc->abr_k_count[cls] == 0) {
    rc->abr_k[cls] = k;
  } else {
    rc->abr_k[cls] += ABR_K_ALPHA * (k - rc->abr_k[cls]);
  }
  rc->abr_k_count[cls]++;
  rc->abr_total_bits += bits;
  rc->abr_frames++;
  rc->abr_inflight_bits -= state->frame->abr_expected_bits;
  rc->abr_inflight_frames--;
  state->frame->abr_pending = false;

  if (ctrl->cfg.vbv_maxrate > 0) {
    rc->vbv_fill -= bits;
    if (rc->vbv_fill < 0 && ctrl->cfg.enable_logging_output) {
      fprintf(stderr, "VBV underflow by %.0f bits (frame %d)\n", -rc->vbv_fill, state->frame->num);
    }
    rc->vbv_fill = MIN((double)ctrl->cfg.vbv_bufsize,
                       rc->vbv_fill + ctrl->cfg.vbv_maxrate / get_framerate(ctrl));
  }
}


/**
 * \brief Average bits model parameter of the inter pictures of a GOP.
 *
 * \param min_count  set to the smallest number of observations of the
 *                   classes, or NULL
 */
static double abr_inter_k_avg(const encoder_control_t *const ctrl, const kvz_rc_data *rc,
                              int *min_count)
{
  double k_avg = 0.0;
  int n = INT32_MAX;
  int classes = 0;
  if (ctrl->cfg.adaptive_gop && ctrl->agop_max_len > ctrl->cfg.gop_len) {
    // The classes of the pictures of a mini-GOP of the maximum length.
    const int len = ctrl->agop_max_len;
    for (int pos = 1; pos <= len; pos++) {
      const int cls = kvz_lookahead_gop_template(&ctrl->cfg, len, pos, pos & -pos);
      k_avg += abr_class_k(ctrl, rc, cls);
      n = MIN(n, rc->abr_k_count[cls]);
      classes++;
    }
  } else {
    classes = MAX(1, ctrl->cfg.gop_len);
    for (int i = 0; i < classes; i++) {
      k_avg += abr_class_k(ctrl, rc, i);
      n = MIN(n, rc->abr_k_count[i]);
    }
  }
  if (min_count) *min_count = n;
  return k_avg / classes;
}


/** \brief Store the bits model information of the frame being started. */
static void abr_register_frame(encoder_state_t * const state, double base_qp)
{
  const encoder_control_t * const ctrl = state->encoder_control;
  kvz_rc_data *rc = state->frame->new_ratecontrol;
  const kvz_la_result *la = state->frame->la;
  const int intra_class = KVZ_MAX_GOP_LENGTH;

  state->frame->abr_pending = true;
  state->frame->abr_qscale = qp_to_qscale(base_qp);
  state->frame->abr_class = state->frame->is_irap ? intra_class :
                            (ctrl->cfg.gop_len > 0 ? state->frame->gop_offset : 0);
  if (la) {
    state->frame->abr_cost = state->frame->is_irap ? MAX(1.0, la->intra_cost) : MAX(1.0, la->cost);
  } else {
    state->frame->abr_cost = 1.0;
  }
  state->frame->abr_expected_bits = abr_class_k(ctrl, rc, state->frame->abr_class) *
                                    state->frame->abr_cost /
                                    pow(state->frame->abr_qscale, ABR_GAMMA);
  rc->abr_inflight_bits += state->frame->abr_expected_bits;
  rc->abr_inflight_frames++;
}


/** \brief Lowest allowed VBV buffer fullness in the simulation. */
#define VBV_MIN_FILL 0.2
/** \brief Maximum QP increase by the VBV. */
#define VBV_MAX_QP_INCREASE 20.0
/** \brief Maximum decrease of the VBV QP increase from one frame to the next. */
#define VBV_MAX_RELEASE 1.0


/**
 * \brief Increase the base quantizer scale of the frame until the VBV
 * buffer simulated over the lookahead does not underflow.
 *
 * The bits of the frame, the frames being coded and the frames in the
 * lookahead are predicted with the bits model. The future frames are
 * assumed to use the same rate factor as the current frame.
 */
static double vbv_limit_qscale(encoder_state_t * const state, double qscale)
{
  const encoder_control_t * const ctrl = state->encoder_control;
  kvz_rc_data *rc = state->frame->new_ratecontrol;
  const kvz_la_result *la = state->frame->la;
  const int intra_class = KVZ_MAX_GOP_LENGTH;
  const double qcomp = ctrl->cfg.qcomp;
  const double bufsize = ctrl->cfg.vbv_bufsize;
  const double refill = ctrl->cfg.vbv_maxrate / get_framerate(ctrl);
  if (!la) return qscale;

  const int cls = state->frame->is_irap ? intra_class :
                  (ctrl->cfg.gop_len > 0 ? state->frame->gop_offset : 0);
  // Be conservative while the bits model has few observations.
#define VBV_UNCERTAINTY(n) (1.0 + 1.0 / (1.0 + (n)))
  const double k_cur = abr_class_k(ctrl, rc, cls) * VBV_UNCERTAINTY(rc->abr_k_count[cls]);
  const double c_cur = state->frame->is_irap ? MAX(1.0, la->intra_cost) : MAX(1.0, la->cost);
  const double r_cur = MAX(1.0, la->blurred_cost);
  double k_avg = 0.0;
  {
    int n;
    k_avg = abr_inter_k_avg(ctrl, rc, &n);
    k_avg = k_avg * VBV_UNCERTAINTY(n);
  }
  const double k_intra = abr_class_k(ctrl, rc, intra_class) * VBV_UNCERTAINTY(rc->abr_k_count[intra_class]);

  // Buffer fullness after the frames being coded.
  double fill_start = rc->vbv_fill - rc->abr_inflight_bits * VBV_UNCERTAINTY(rc->abr_frames / MAX(1, ctrl->cfg.gop_len));
  fill_start = MIN(bufsize, fill_start + rc->abr_inflight_frames * refill);

  const int64_t last = kvz_lookahead_newest(ctrl->lookahead);
  // The increase is released gradually since the bits model is less
  // accurate far from the recently used quantizers.
  const double step = pow(2.0, 0.5 / 6.0);
  double mult = pow(2.0, MAX(0.0, rc->vbv_qp_increase - VBV_MAX_RELEASE) / 6.0);
  for (; mult < pow(2.0, VBV_MAX_QP_INCREASE / 6.0); mult *= step) {
    const double qs = qscale * mult;
    double fill = fill_start - k_cur * c_cur / pow(qs, ABR_GAMMA);
    bool ok = fill >= VBV_MIN_FILL * bufsize;
    fill = MIN(bufsize, fill + refill);
    for (int64_t i = la->index + 1; ok && i <= last; i++) {
      double cost, intra;
      bool key;
      if (!kvz_lookahead_frame_info(ctrl->lookahead, i, &cost, &intra, &key)) continue;
      const double qs_i = qs * pow(cost / r_cur, 1.0 - qcomp);
      const double bits = key ? k_intra * intra / pow(qs_i, ABR_GAMMA) :
                                k_avg * cost / pow(qs_i, ABR_GAMMA);
      fill -= bits;
      if (fill < VBV_MIN_FILL * bufsize) ok = false;
      fill = MIN(bufsize, fill + refill);
    }
    if (ok) break;
  }
  rc->vbv_qp_increase = 6.0 * log2(mult);
  return qscale * mult;
}


void kvz_set_picture_qp_abr(encoder_state_t * const state)
{
  const encoder_control_t * const ctrl = state->encoder_control;
  kvz_rc_data *rc = state->frame->new_ratecontrol;
  const int intra_class = KVZ_MAX_GOP_LENGTH;
  const double qcomp = ctrl->cfg.qcomp;

  abr_account_previous(state);

  const kvz_la_result *la = state->frame->la;
  const double framerate = ctrl->cfg.framerate_num != 0 ?
    ctrl->cfg.framerate_num / (double)ctrl->cfg.framerate_denom : ctrl->cfg.framerate;
  const double target = ctrl->target_avg_bppic;

  // Expected bits of the window of frames at rate factor 1.
  int64_t first = la ? la->index : 0;
  if (la && ctrl->cfg.gop_len > 0 && !ctrl->cfg.gop_lowdelay && !state->frame->new_segment) {
    first = MAX(0, la->index - ctrl->cfg.gop_len + 1);
  }
  const int64_t last = kvz_lookahead_newest(ctrl->lookahead);
  const double k_avg = abr_inter_k_avg(ctrl, rc, NULL);
  const double k_intra = abr_class_k(ctrl, rc, intra_class);

  double expected = 0.0;
  int n = 0;
  double key_sum = 0.0;
  int key_n = 0;
  for (int64_t i = first; i <= last; i++) {
    double cost, intra;
    bool key;
    if (!kvz_lookahead_frame_info(ctrl->lookahead, i, &cost, &intra, &key)) continue;
    const double rceq = pow(cost, (1.0 - qcomp) * ABR_GAMMA);
    if (key) {
      key_sum += k_intra * intra / rceq;
      key_n++;
    } else {
      expected += k_avg * cost / rceq;
    }
    n++;
  }
  if (ctrl->cfg.intra_period > 1) {
    // The bits of the periodic intra pictures are spread over the intra
    // period. Otherwise the rate factor would swing depending on whether
    // the window happens to contain an intra picture.
    if (key_n > 0) rc->abr_key_est = key_sum / key_n;
    expected += rc->abr_key_est * n / MAX(n, ctrl->cfg.intra_period);
  } else {
    expected += key_sum;
  }

  // The rate factor is updated once per GOP so that all the pictures of a
  // GOP use the same rate factor.
  // Adaptive GOP: at the anchor and, in mini-GOPs longer than the GOP
  // table, also at the pictures that start the halves of the hierarchy so
  // that the rate factor is updated as often as with the GOP table.
  const bool agop_anchor = la && la->gop_size > 0 &&
                           (la->gop_pos == la->gop_size ||
                            (la->gop_size > ctrl->cfg.gop_len && la->gop_span >= ctrl->cfg.gop_len / 2));
  const bool gop_start = ctrl->cfg.gop_len == 0 || ctrl->cfg.gop_lowdelay ||
                         state->frame->gop_offset == 0 || state->frame->new_segment ||
                         agop_anchor || rc->abr_rate_factor <= 0.0;
  if (gop_start && n > 0 && expected > 0.0) {
    // Spread the accumulated error over ABR_ERROR_SECONDS.
    // Include the expected bits of the frames that are still being coded.
    const double error = rc->abr_total_bits + rc->abr_inflight_bits -
                         target * (rc->abr_frames + rc->abr_inflight_frames);
    // Early in the sequence the error is corrected faster.
    const double horizon = MAX((double)n, MIN(ABR_ERROR_SECONDS * framerate,
                                              (double)(rc->abr_frames + rc->abr_inflight_frames)));
    double budget = target * n - error * n / horizon;
    budget = MAX(budget, 0.3 * target * n);
    double rate_factor = pow(budget / expected, 1.0 / ABR_GAMMA);
    if (rc->abr_rate_factor > 0.0) {
      // Limit the change of the base QP between GOPs, also at scene cuts
      // where the bits model of the new scene is not known yet.
      const double max_step = ABR_MAX_QP_STEP;
      rate_factor = CLIP(rc->abr_rate_factor * pow(2.0, -max_step / 6.0),
                         rc->abr_rate_factor * pow(2.0, max_step / 6.0),
                         rate_factor);
    }
    rc->abr_rate_factor = rate_factor;
  }

  double qscale;
  if (rc->abr_rate_factor > 0.0) {
    const double cost = la ? MAX(1.0, la->blurred_cost) : 1.0;
    qscale = pow(cost, 1.0 - qcomp) / rc->abr_rate_factor;
  } else {
    qscale = qp_to_qscale(ctrl->cfg.qp);
  }

  if (ctrl->cfg.vbv_maxrate > 0) {
    qscale = vbv_limit_qscale(state, qscale);
  }

  double base_qp = CLIP(0.0, 51.0, qscale_to_qp(qscale));
  if (rc->abr_frames + rc->abr_inflight_frames == 0) {
    // Nothing has been coded yet and the bits model is only a prior. Do not
    // code the first picture at a much lower QP than the bits per pixel
    // suggest, which could spend the budget of many seconds on it.
    const double bpp = target / ctrl->in.pixels_per_pic;
    const double qp_bpp = ABR_BPP_QP + 6.0 / ABR_GAMMA * log2(ABR_BPP_REF / MAX(bpp, 1e-6));
    base_qp = CLIP(0.0, 51.0, MAX(base_qp, qp_bpp - ABR_FIRST_QP_MARGIN));
  }

  abr_register_frame(state, base_qp);

  state->frame->QP = gop_frame_qp(state, base_qp);
  state->frame->lambda = qp_to_lambda(state, state->frame->QP);
}


void kvz_set_qg_qp_lambda(encoder_state_t * const state, int x, int y)
{
  const encoder_control_t * const ctrl = state->encoder_control;
  if (!ctrl->cfg.cutree_qg32 || state->frame->max_qp_delta_depth < 1 ||
      state->tile->frame->source->roi.roi_array) {
    return;
  }
  const int w32 = CEILDIV(ctrl->in.width, 32);
  const int h32 = CEILDIV(ctrl->in.height, 32);
  // CUs of a CTU partially outside the picture may start below or to the
  // right of the last quantization group of the offset array.
  const int bx = MIN((x + state->tile->offset_x) / 32, w32 - 1);
  const int by = MIN((y + state->tile->offset_y) / 32, h32 - 1);
  int qp = state->qp_base + (int)round(state->frame->aq_offsets32[bx + by * w32]);
  qp = CLIP(state->frame->QP + KVZ_QP_DELTA_MIN / 2, state->frame->QP + KVZ_QP_DELTA_MAX / 2, qp);
  state->qp = CLIP_TO_QP(qp);
  state->lambda = qp_to_lambda(state, state->qp);
  state->lambda_sqrt = sqrt(state->lambda);
}


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

    // With CRF the base QP comes from the lookahead complexity estimate.
    double base_qp = ctrl->cfg.qp;
    const bool crf = ctrl->cfg.crf >= 0 && state->frame->la;
    if (crf) {
      base_qp = CLIP(0.0, 51.0, state->frame->la->crf_qp);
      if (ctrl->cfg.vbv_maxrate > 0) {
        abr_account_previous(state);
        base_qp = CLIP(0.0, 51.0, qscale_to_qp(vbv_limit_qscale(state, qp_to_qscale(base_qp))));
        abr_register_frame(state, base_qp);
      }
    }

    if (gop_len > 0 && state->frame->slicetype != KVZ_SLICE_I) {
      double qp = base_qp;
      qp += gop->qp_offset;
      qp += CLIP(0.0, 3.0, qp * gop->qp_model_scale + gop->qp_model_offset);
      state->frame->QP = CLIP_TO_QP((int)(qp + 0.5));

    }
    else if (crf) {
      state->frame->QP = CLIP_TO_QP((int)(base_qp + ctrl->cfg.intra_qp_offset + 0.5));
    }
    else {
      state->frame->QP = CLIP_TO_QP(ctrl->cfg.qp + ctrl->cfg.intra_qp_offset);
    }

    state->frame->lambda = qp_to_lambda(state, state->frame->QP);
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
  lcu_stats_t *lcu = kvz_get_lcu_stats(state, pos.x, pos.y);

  if (state->tile->frame->source->roi.roi_array) {
    vector2d_t lcu_vec = {
      pos.x + state->tile->lcu_offset_x,
      pos.y + state->tile->lcu_offset_y
    };
    vector2d_t roi = {
      lcu_vec.x * state->tile->frame->source->roi.width / ctrl->in.width_in_lcu,
      lcu_vec.y * state->tile->frame->source->roi.height / ctrl->in.height_in_lcu
    };
    int roi_index = roi.x + roi.y * state->tile->frame->source->roi.width;
    int dqp = state->tile->frame->source->roi.roi_array[roi_index];
    if(dqp != 0) {
      pos.x = 0;
    }
    state->qp = CLIP_TO_QP(state->frame->QP + dqp);
    state->lambda = qp_to_lambda(state, state->qp);
    state->lambda_sqrt = sqrt(state->lambda);
  }
  else if (ctrl->cfg.target_bitrate > 0 && ctrl->cfg.rc_algorithm == KVZ_LAMBDA) {
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

    state->lambda      = lambda;
    state->lambda_sqrt = sqrt(lambda);
    state->qp          = lambda_to_qp(lambda);

  } else {
    state->qp          = state->frame->QP;
    state->lambda      = state->frame->lambda;
    state->lambda_sqrt = sqrt(state->frame->lambda);
  }

  lcu->lambda = state->lambda;
  lcu->qp = state->qp;
  state->qp_base = state->qp;

  // Apply variance adaptive quantization
  if (ctrl->cfg.vaq || ctrl->cfg.cutree > 0) {
    vector2d_t lcu_pos = {
      pos.x + state->tile->lcu_offset_x,
      pos.y + state->tile->lcu_offset_y
    };
    int id = lcu_pos.x + lcu_pos.y * state->tile->frame->width_in_lcu;
    int aq_offset = round(state->frame->aq_offsets[id]);
    state->qp += aq_offset;    
    // Maximum delta QP is clipped according to ITU T-REC-H.265 specification chapter 7.4.9.10 Transform unit semantics
    // Clipping range is a function of bit depth
    // Since this value will be later combined with qp_pred, clip to half of that instead to be safe
    state->qp = CLIP(state->frame->QP + KVZ_QP_DELTA_MIN / 2, state->frame->QP + KVZ_QP_DELTA_MAX / 2, state->qp);
    state->qp = CLIP_TO_QP(state->qp);
    state->lambda = qp_to_lambda(state, state->qp);
    state->lambda_sqrt = sqrt(state->lambda);

    lcu->adjust_lambda = state->lambda;
    lcu->adjust_qp = state->qp;
  }
}
