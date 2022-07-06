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

#ifndef FAST_COEFF_COST_H_
#define FAST_COEFF_COST_H_

#include <stdio.h>
#include "kvazaar.h"
// #include "encoderstate.h"

#define MAX_FAST_COEFF_COST_QP 50

typedef struct {
  uint64_t wts_by_qp[MAX_FAST_COEFF_COST_QP];
} fast_coeff_table_t;

// Weights for 4 buckets (coeff 0, coeff 1, coeff 2, coeff >= 3), for QPs from
// 0 to MAX_FAST_COEFF_COST_QP
static const double default_fast_coeff_cost_wts[][4] = {
{0.162000, 4.126087, 3.499517, 6.969847},
{0.162000, 4.126087, 3.499517, 6.969847},
{0.162000, 4.126087, 3.499517, 6.969847},
{0.162000, 4.126087, 3.499517, 6.969847},
{0.162000, 4.126087, 3.499517, 6.969847},
{0.162000, 4.126087, 3.499517, 6.969847},
{0.162000, 4.126087, 3.499517, 6.969847},
{0.162000, 4.126087, 3.499517, 6.969847},
{0.162000, 4.126087, 3.499517, 6.969847},
{0.162000, 4.126087, 3.499517, 6.969847},
{0.162000, 4.126087, 3.499517, 6.969847},
{0.157760, 4.037673, 3.558663, 6.895640},
{0.127943, 4.308060, 3.916680, 6.962907},
{0.110555, 4.422860, 3.944640, 6.898343},
{0.094532, 4.479287, 4.161790, 6.804273},
{0.074032, 4.629857, 4.042727, 6.722910},
{0.051644, 4.960970, 4.001523, 6.556783},
{0.039513, 5.133963, 3.951247, 6.472487},
{0.034188, 5.185183, 3.805350, 6.418810},
{0.028981, 5.203517, 3.785043, 6.351090},
{0.022543, 5.315690, 3.796553, 6.347457},
{0.020300, 5.221910, 3.817927, 6.322733},
{0.015400, 5.170127, 3.937963, 6.326643},
{0.010147, 5.088577, 4.143093, 6.293030},
{0.008239, 5.017160, 4.204780, 6.267220},
{0.006386, 4.956723, 4.303120, 6.208533},
{0.004876, 4.912990, 4.400863, 6.175370},
{0.003707, 4.905997, 4.388617, 6.134007},
{0.003089, 4.872320, 4.521937, 6.153827},
{0.002479, 4.864330, 4.591423, 6.152587},
{0.002180, 4.864427, 4.607133, 6.141223},
{0.002556, 4.771863, 4.793583, 6.232397},
{0.001316, 4.793543, 4.787927, 6.272543},
{0.001169, 4.845383, 4.787190, 6.235333},
{0.001000, 4.849327, 4.805003, 6.273347},
{0.000830, 4.839947, 4.866000, 6.346927},
{0.001131, 4.772140, 4.969497, 6.448050},
{0.000553, 4.743423, 5.050670, 6.663760},
{0.000466, 4.800883, 5.034373, 6.601250},
{0.000400, 4.797313, 5.079183, 6.743547},
{0.000333, 4.783170, 5.142737, 6.869933},
{0.000355, 4.915657, 5.217510, 7.225673},
{0.000186, 4.973477, 5.151287, 7.280497},
{0.000113, 5.316010, 4.509893, 6.585287},
{0.000091, 5.304703, 4.553107, 6.773803},
{0.000076, 5.263460, 4.689990, 6.962153},
{0.000064, 5.190947, 4.733550, 7.100820},
{0.000053, 5.180677, 4.833283, 7.340667},
{0.000047, 5.182963, 4.829380, 7.338863},
{0.000032, 5.389257, 4.518127, 7.265003},
{0.000020, 5.970297, 3.981997, 7.201180},
{0.000000, 0.000000, 0.000000, 0.000000},


};

typedef struct encoder_state_t encoder_state_t;

int kvz_fast_coeff_table_parse(fast_coeff_table_t *fast_coeff_table, FILE *fast_coeff_table_f);
void kvz_fast_coeff_use_default_table(fast_coeff_table_t *fast_coeff_table);
uint64_t kvz_fast_coeff_get_weights(const encoder_state_t *state);

#endif // FAST_COEFF_COST_H_
