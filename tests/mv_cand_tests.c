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

#include "src/inter.c"

#include <string.h>

#include "greatest/greatest.h"

TEST test_get_spatial_merge_cand(void)
{
  lcu_t lcu;
  memset(&lcu, 0, sizeof(lcu));
  for (int i = 0; i < sizeof(lcu.cu) / sizeof(cu_info_t); i++) {
    lcu.cu[i].type = CU_INTER;
  }

  merge_candidates_t cand = { {0, 0}, {0, 0, 0}, 0, 0 };

  get_spatial_merge_candidates(64 + 32, 64, // x, y
                               32, 24,      // width, height
                               1920, 1080,  // picture size
                               &lcu,
                               &cand);

  ASSERT_EQ(cand.b[0], &lcu.cu[289]);
  ASSERT_EQ(cand.b[1], &lcu.cu[ 16]);
  ASSERT_EQ(cand.b[2], &lcu.cu[  8]);
  ASSERT_EQ(cand.a[0], &lcu.cu[127]);
  ASSERT_EQ(cand.a[1], &lcu.cu[110]);

  PASS();
}

TEST test_is_a0_cand_coded()
{
  // +--+--+
  // |##|  |
  // +--+--+
  // |  |  |
  // +--+--+
  ASSERT_EQ(is_a0_cand_coded(32, 64, 16, 16), true);
  // Same as above with a 2NxN block
  ASSERT_EQ(is_a0_cand_coded(32, 64, 32, 16), true);
  // Same as above with a 2NxnU block
  ASSERT_EQ(is_a0_cand_coded(32, 64, 32, 8), true);
  // Same as above with a 2NxnD block
  ASSERT_EQ(is_a0_cand_coded(32, 64, 32, 24), true);

  // +--+--+
  // |  |##|
  // +--+--+
  // |  |  |
  // +--+--+
  ASSERT_EQ(is_a0_cand_coded(16, 0, 16, 16), false);

  // +--+--+
  // |  |  |
  // +--+--+
  // |  |##|
  // +--+--+
  ASSERT_EQ(is_a0_cand_coded(48, 16, 16, 16), false);
  // Same as above with a Nx2N block
  ASSERT_EQ(is_a0_cand_coded(48, 0, 16, 32), false);
  // Same as above with a nLx2N block
  ASSERT_EQ(is_a0_cand_coded(40, 0, 24, 32), false);
  // Same as above with a nRx2N block
  ASSERT_EQ(is_a0_cand_coded(56, 0, 8, 32), false);

  // +-----+--+--+
  // |     |  |  |
  // |     +--+--+
  // |     |##|  |
  // +-----+--+--+
  // |     |     |
  // |     |     |
  // |     |     |
  // +-----+-----+
  ASSERT_EQ(is_a0_cand_coded(32, 16, 16, 16), false);

  // Same as above with a 2NxnU block
  ASSERT_EQ(is_a0_cand_coded(32, 8, 32, 24), false);
  // Same as above with a 2NxnD block
  ASSERT_EQ(is_a0_cand_coded(32, 24, 32, 8), false);

  // Same as above with a Nx2N block
  ASSERT_EQ(is_a0_cand_coded(32, 0, 16, 32), false);
  // Same as above with a nLx2N block
  ASSERT_EQ(is_a0_cand_coded(32, 0, 8, 32), false);
  // Same as above with a nRx2N block
  ASSERT_EQ(is_a0_cand_coded(32, 0, 24, 32), false);

  // +--+--+-----+
  // |  |  |     |
  // +--+--+     |
  // |##|  |     |
  // +--+--+-----+
  // |     |     |
  // |     |     |
  // |     |     |
  // +-----+-----+
  ASSERT_EQ(is_a0_cand_coded(32, 8, 8, 8), true);

  // Same as above with a 2NxnU block
  ASSERT_EQ(is_a0_cand_coded(32, 4, 16, 12), true);
  // Same as above with a 2NxnD block
  ASSERT_EQ(is_a0_cand_coded(32, 12, 16, 4), true);

  // Same as above with a Nx2N block
  ASSERT_EQ(is_a0_cand_coded(32, 0, 8, 16), true);
  // Same as above with a nLx2N block
  ASSERT_EQ(is_a0_cand_coded(32, 0, 4, 16), true);
  // Same as above with a nRx2N block
  ASSERT_EQ(is_a0_cand_coded(32, 0, 12, 16), true);

  PASS();
}

TEST test_is_b0_cand_coded()
{
  // +--+--+
  // |##|  |
  // +--+--+
  // |  |  |
  // +--+--+
  ASSERT_EQ(is_b0_cand_coded(32, 64, 16, 16), true);
  // Same as above with a Nx2N block
  ASSERT_EQ(is_b0_cand_coded(32, 64, 16, 32), true);
  // Same as above with a nLx2N block
  ASSERT_EQ(is_b0_cand_coded(32, 64, 24, 32), true);
  // Same as above with a nRx2N block
  ASSERT_EQ(is_b0_cand_coded(32, 64, 8, 32), true);

  // +--+--+
  // |  |  |
  // +--+--+
  // |##|  |
  // +--+--+
  ASSERT_EQ(is_b0_cand_coded(32, 16, 16, 16), true);

  // +--+--+
  // |  |  |
  // +--+--+
  // |  |##|
  // +--+--+
  ASSERT_EQ(is_b0_cand_coded(48, 16, 16, 16), false);
  // Same as above with a 2NxN block
  ASSERT_EQ(is_b0_cand_coded(32, 16, 32, 16), false);
  // Same as above with a 2NxnU block
  ASSERT_EQ(is_b0_cand_coded(32, 8, 32, 24), false);
  // Same as above with a 2NxnD block
  ASSERT_EQ(is_b0_cand_coded(32, 24, 32, 8), false);

  // +-----+-----+
  // |     |     |
  // |     |     |
  // |     |     |
  // +-----+--+--+
  // |     |  |##|
  // |     +--+--+
  // |     |  |  |
  // +-----+--+--+
  ASSERT_EQ(is_b0_cand_coded(48, 32, 16, 16), false);

  // Same as above with a 2NxnU block
  ASSERT_EQ(is_b0_cand_coded(32, 32, 32, 8), false);
  // Same as above with a 2NxnD block
  ASSERT_EQ(is_b0_cand_coded(32, 32, 32, 24), false);

  // Same as above with a nLx2N block
  ASSERT_EQ(is_b0_cand_coded(56, 32, 8, 32), false);
  // Same as above with a nRx2N block
  ASSERT_EQ(is_b0_cand_coded(40, 32, 24, 32), false);

  // +--+--+-----+
  // |  |##|     |
  // +--+--+     |
  // |  |  |     |
  // +--+--+-----+
  // |     |     |
  // |     |     |
  // |     |     |
  // +-----+-----+
  ASSERT_EQ(is_b0_cand_coded(16, 0, 16, 16), true);

  // Same as above with a 2NxnU block
  ASSERT_EQ(is_b0_cand_coded(0, 0, 32, 8), true);
  // Same as above with a 2NxnD block
  ASSERT_EQ(is_b0_cand_coded(0, 0, 32, 24), true);

  // Same as above with a nLx2N block
  ASSERT_EQ(is_b0_cand_coded(8, 0, 24, 32), true);
  // Same as above with a nRx2N block
  ASSERT_EQ(is_b0_cand_coded(24, 0, 8, 32), true);

  PASS();
}

SUITE(mv_cand_tests) {
  RUN_TEST(test_get_spatial_merge_cand);
  RUN_TEST(test_is_a0_cand_coded);
  RUN_TEST(test_is_b0_cand_coded);
}
