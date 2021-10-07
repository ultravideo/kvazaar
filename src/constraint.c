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

#include "constraint.h"

 /**
  * \brief Allocate the constraint_t structure.
  *
  * \param state   encoder state
  * \return the pointer of constraint_t structure
  */
void * kvz_init_constraint(encoder_state_t* state, const encoder_control_t * const encoder) {
  constraint_t* constr = NULL;
  // Allocate the constraint_t strucutre
  constr = MALLOC(constraint_t, 1);
  if (!constr) {
    fprintf(stderr, "Memory allocation failed!\n");
    assert(0);
  }

  // Allocate the ml_intra_ctu_pred_t structure
  constr->ml_intra_depth_ctu = NULL;
  if (encoder->cfg.ml_pu_depth_intra) // TODO: Change this by a new param !!
  {
    constr->ml_intra_depth_ctu = kvz_init_ml_intra_depth_const();
  }
  return constr;
}

/**
 * \brief Deallocate the constraint_t structure.
 *
 * \param state   encoder state
 */
void kvz_constraint_free(encoder_state_t* state) {
  constraint_t* constr = state->constraint;
  if (constr->ml_intra_depth_ctu) 
  {
    kvz_end_ml_intra_depth_const(constr->ml_intra_depth_ctu);
  }
  FREE_POINTER(constr);
}
