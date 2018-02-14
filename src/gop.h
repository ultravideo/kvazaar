#ifndef GOP_H_
#define GOP_H_
/*****************************************************************************
* This file is part of Kvazaar HEVC encoder.
*
* Copyright (C) 2018 Tampere University of Technology and others (see
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

#include <kvazaar.h>


static const kvz_gop_config kvz_gop_ra8[8] = {
  {
    .poc_offset = 8,
    .layer      = 1,
    .qp_offset  = 1,
    .qp_factor  = 0.442,
    .is_ref     = 1,
    .ref_neg_count = 3,
    .ref_neg = { 8, 12, 16 },
    .ref_pos_count = 0,
    .ref_pos = { 0 },
  },
  {
    .poc_offset = 4,
    .layer      = 2,
    .qp_offset  = 2,
    .qp_factor  = 0.3536,
    .is_ref     = 1,
    .ref_neg_count = 2,
    .ref_neg = { 4, 8 },
    .ref_pos_count = 1,
    .ref_pos = { 4 },
  },
  {
    .poc_offset = 2,
    .layer      = 3,
    .qp_offset  = 3,
    .qp_factor  = 0.3536,
    .is_ref     = 1,
    .ref_neg_count = 2,
    .ref_neg = { 2, 6 },
    .ref_pos_count = 2,
    .ref_pos = { 2, 6 },
  },
  {
    .poc_offset = 1,
    .layer      = 4,
    .qp_offset  = 4,
    .qp_factor  = 0.68,
    .is_ref     = 0,
    .ref_neg_count = 1,
    .ref_neg = { 1 },
    .ref_pos_count = 3,
    .ref_pos = { 1, 3, 7 },
  },
  {
    .poc_offset = 3,
    .layer      = 4,
    .qp_offset  = 4,
    .qp_factor  = 0.68,
    .is_ref     = 0,
    .ref_neg_count = 2,
    .ref_neg = { 1, 3 },
    .ref_pos_count = 2,
    .ref_pos = { 1, 5 },
  },
  {
    .poc_offset = 6,
    .layer      = 3,
    .qp_offset  = 3,
    .qp_factor  = 0.3536,
    .is_ref     = 1,
    .ref_neg_count = 2,
    .ref_neg = { 2, 6 },
    .ref_pos_count = 1,
    .ref_pos = { 2 },
  },
  {
    .poc_offset = 5,
    .layer      = 4,
    .qp_offset  = 4,
    .qp_factor  = 0.68,
    .is_ref     = 0,
    .ref_neg_count = 2,
    .ref_neg = { 1, 5 },
    .ref_pos_count = 2,
    .ref_pos = { 1, 3 },
  },
  {
    .poc_offset = 7,
    .layer      = 4,
    .qp_offset  = 4,
    .qp_factor  = 0.68,
    .is_ref     = 0,
    .ref_neg_count = 3,
    .ref_neg = { 1, 3, 7 },
    .ref_pos_count = 1,
    .ref_pos = { 1 },
  },
};

#endif // GOP_H_
