#ifndef ML_CLASSIFIER_INTRA_DEPTH_PRED
#define ML_CLASSIFIER_INTRA_DEPTH_PRED

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

#include "ml_intra_cu_depth_pred.h"


int tree_predict_merge_depth_1(features_s* p_features, double* p_nb_iter, double* p_nb_bad);
int tree_predict_merge_depth_2(features_s* p_features, double* p_nb_iter, double* p_nb_bad);
int tree_predict_merge_depth_3(features_s* p_features, double* p_nb_iter, double* p_nb_bad);
int tree_predict_merge_depth_4(features_s* p_features, double* p_nb_iter, double* p_nb_bad);


int tree_predict_split_depth_0(features_s* p_features, double* p_nb_iter, double* p_nb_bad);
int tree_predict_split_depth_1(features_s* p_features, double* p_nb_iter, double* p_nb_bad);
int tree_predict_split_depth_2(features_s* p_features, double* p_nb_iter, double* p_nb_bad);
int tree_predict_split_depth_3(features_s* p_features, double* p_nb_iter, double* p_nb_bad);

#endif