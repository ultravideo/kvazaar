#ifndef SEARCH_H_
#define SEARCH_H_
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

/**
 * \ingroup Compression
 * \file
 * \brief Compression of a single coding tree unit (CTU).
 */

#include "cu.h"
#include "encoderstate.h"
#include "global.h" // IWYU pragma: keep
#include "image.h"
#include "constraint.h"

#define MAX_UNIT_STATS_MAP_SIZE MAX(MAX_REF_PIC_COUNT, MRG_MAX_NUM_CANDS)

 // Modify weight of luma SSD.
#ifndef KVZ_LUMA_MULT
# define KVZ_LUMA_MULT 0.8
#endif
// Modify weight of chroma SSD.
#ifndef KVZ_CHROMA_MULT
# define KVZ_CHROMA_MULT 1.5
#endif

 /**
  *  \brief Data collected during search processes.
  * 
  *         The intended use is to collect statistics of the
  *         searched coding/prediction units. Data related to
  *         a specific unit is found at index i. The arrays
  *         should be indexed by elements of the "keys" array
  *         that will be sorted by the RD costs of the units.         
  */
typedef struct unit_stats_map_t {

  cu_info_t unit[MAX_UNIT_STATS_MAP_SIZE]; //!< list of searched units
  double    cost[MAX_UNIT_STATS_MAP_SIZE]; //!< list of matching RD costs
  double    bits[MAX_UNIT_STATS_MAP_SIZE]; //!< list of matching bit costs  
  int8_t    keys[MAX_UNIT_STATS_MAP_SIZE]; //!< list of keys (indices) to elements in the other arrays
  int       size;                    //!< number of active elements in the lists
} unit_stats_map_t;

void kvz_sort_modes(int8_t *__restrict modes, double *__restrict costs, uint8_t length);
void kvz_sort_keys_by_cost(unit_stats_map_t *__restrict map);

void kvz_search_lcu(encoder_state_t *state, int x, int y, const yuv_t *hor_buf, const yuv_t *ver_buf);

double kvz_cu_rd_cost_luma(const encoder_state_t *const state,
                           const int x_px, const int y_px, const int depth,
                           const cu_info_t *const pred_cu,
                           lcu_t *const lcu);
double kvz_cu_rd_cost_chroma(const encoder_state_t *const state,
                             const int x_px, const int y_px, const int depth,
                             const cu_info_t *const pred_cu,
                             lcu_t *const lcu);
void kvz_lcu_fill_trdepth(lcu_t *lcu, int x_px, int y_px, int depth, int tr_depth);

void kvz_intra_recon_lcu_luma(encoder_state_t * const state, int x, int y, int depth, int8_t intra_mode, cu_info_t *cur_cu, lcu_t *lcu);
void kvz_intra_recon_lcu_chroma(encoder_state_t * const state, int x, int y, int depth, int8_t intra_mode, cu_info_t *cur_cu, lcu_t *lcu);

#endif
