#ifndef RDO_H_
#define RDO_H_
/**
 * \file
 * \brief Handling Rate-Distortion Optimization related functionality
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "global.h"

#include "encoder.h"


typedef struct
{  
  double coded_level_and_dist;
  double uncoded_dist;
  double sig_cost;
  double sig_cost_0;
  int32_t nnz_before_pos0;
} coeffgroup_rd_stats;

const uint32_t g_go_rice_range[5] = { 7, 14, 26, 46, 78 };
const uint32_t g_go_rice_prefix_len[5] = { 8, 7, 6, 5, 4 };


void  rdoq(encoder_control *encoder, coefficient *coef, coefficient *dest_coeff, int32_t width,
           int32_t height, uint32_t *abs_sum, int8_t type, int8_t scan_idx, int8_t block_type, int8_t scan_mode, int8_t tr_depth);


int32_t get_ic_rate( uint32_t abs_level, uint16_t ctx_num_one,uint16_t ctx_num_abs,
                     uint16_t abs_go_rice, uint32_t c1_idx, uint32_t c2_idx, int8_t type);
uint32_t get_coded_level ( encoder_control* encoder, double* coded_cost, double* coded_cost0, double* coded_cost_sig,
                           int32_t level_double, uint32_t max_abs_level,
                           uint16_t ctx_num_sig, uint16_t ctx_num_one, uint16_t ctx_num_abs,
                           uint16_t abs_go_rice,
                           uint32_t c1_idx, uint32_t c2_idx,
                           int32_t q_bits,double temp, int8_t last, int8_t type);


#endif
