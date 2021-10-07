#ifndef ML_INTRA_CU_DEPTH_PRED_H_
#define ML_INTRA_CU_DEPTH_PRED_H_
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

#include <stdio.h>
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"




#define LCU_DEPTH_MAT_SIZE 64
#define RESTRAINED_FLAG 1

#define pow2(x) ((x)*(x))
#define CR_XMAX(x_px, block_size, width)       (MIN((x_px) + (block_size), (width))  - (x_px))
#define CR_YMAX(y_px, block_size, height)      (MIN((y_px) + (block_size), (height)) - (y_px))
#define CR_GET_X_LCU(lcu_id, nb_lcu_width)     (((lcu_id) % (nb_lcu_width)) << 6)
#define CR_GET_Y_LCU(lcu_id, nb_lcu_width)     (((lcu_id) / (nb_lcu_width)) << 6)
#define CR_GET_CU_D3(x, y, depth) ((x)*(1 << (3-depth)) + ((y) << (6 - depth)))
#define CR_GET_CU_D4(x, y, depth) ((x)*(1 << (4-depth)) + ((y) << (8 - depth)))
#define CR_GET_DEPTH_MIN(x, y, depth_min_mat) *(depth_min_mat + (x >> 3) + ((y >> 3) << 3))
#define CR_GET_DEPTH_MAX(x, y, depth_max_mat) *(depth_max_mat + (x >> 3) + ((y >> 3) << 3))

typedef struct {
	int32_t x;
	int32_t y;
}vect_2D;


 // Structure used for the CTU depth prediction using Machine Learning 
 // in All Intra 
typedef struct {
	/*!< Number of depth to add to the QT prediction in ''one-shot'' */
	int8_t   i_nb_addDepth;
	/*!< Apply an extra Upper Expansion in the upper_depth */
	bool	 b_extra_up_exp;
	/*!< Matrix used to store the upper and lower QT prediction*/
	uint8_t* _mat_upper_depth; 
	uint8_t* _mat_lower_depth;
} ml_intra_ctu_pred_t;



/*
 * brief generic structure used for the features
 *
 */
typedef struct {
	double variance;
	double merge_variance;
	double sub_variance_0;
	double sub_variance_1;
	double sub_variance_2;
	double sub_variance_3;
	double neigh_variance_A;
	double neigh_variance_B;
	double neigh_variance_C;
	double var_of_sub_mean;
	int 	qp;
	//int   NB_pixels;
	double var_of_sub_var;
}features_s;


typedef int (*tree_predict)(features_s*, double*, double*);

ml_intra_ctu_pred_t* kvz_init_ml_intra_depth_const(void);
void kvz_end_ml_intra_depth_const(ml_intra_ctu_pred_t * ml_intra_depth_ctu);

void kvz_lcu_luma_depth_pred(ml_intra_ctu_pred_t* ml_intra_depth_ctu, kvz_pixel* luma_px, int8_t qp);

#endif
