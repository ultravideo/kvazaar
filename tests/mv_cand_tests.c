/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License version 2.1 as
 * published by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include <string.h>

#include "greatest/greatest.h"

#include "src/cu.h"
#include "src/inter.h"

TEST test_get_spatial_merge_cand(void)
{
  lcu_t lcu;
  memset(&lcu, 0, sizeof(lcu));
  lcu.cu[20].coded = 1;
  lcu.cu[22].coded = 1;
  lcu.cu[23].coded = 1;
  lcu.cu[56].coded = 1;
  lcu.cu[65].coded = 1;

  cu_info_t *mv_cand[5] = { NULL };
  kvz_inter_get_spatial_merge_candidates(16, 16, // x, y
                                         16, 32, // width, height
                                         &mv_cand[0], // b0
                                         &mv_cand[1], // b1
                                         &mv_cand[2], // b2
                                         &mv_cand[3], // a0
                                         &mv_cand[4], // a1
                                         &lcu);

  ASSERT_EQ(mv_cand[0], &lcu.cu[23]); // b0
  ASSERT_EQ(mv_cand[1], &lcu.cu[22]); // b1
  ASSERT_EQ(mv_cand[2], &lcu.cu[20]); // b2
  ASSERT_EQ(mv_cand[3], &lcu.cu[65]); // a0
  ASSERT_EQ(mv_cand[4], &lcu.cu[56]); // a1

  PASS();
}

SUITE(mv_cand_tests) {
  RUN_TEST(test_get_spatial_merge_cand);
}
