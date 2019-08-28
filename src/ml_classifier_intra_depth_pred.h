#ifndef ML_CLASSIFIER_INTRA_DEPTH_PRED
#define ML_CLASSIFIER_INTRA_DEPTH_PRED

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