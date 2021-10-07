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

#include "greatest/greatest.h"

#include "test_strategies.h"

GREATEST_MAIN_DEFS();
#if KVZ_BIT_DEPTH == 8
extern SUITE(sad_tests);
extern SUITE(intra_sad_tests);
extern SUITE(satd_tests);
extern SUITE(speed_tests);
extern SUITE(dct_tests);
#endif //KVZ_BIT_DEPTH == 8

extern SUITE(coeff_sum_tests);
extern SUITE(mv_cand_tests);
extern SUITE(inter_recon_bipred_tests);

int main(int argc, char **argv)
{
  GREATEST_MAIN_BEGIN();

  init_test_strategies(1);
#if KVZ_BIT_DEPTH == 8
  RUN_SUITE(sad_tests);
  RUN_SUITE(intra_sad_tests);
  RUN_SUITE(satd_tests);
  RUN_SUITE(dct_tests);

  if (greatest_info.suite_filter &&
      greatest_name_match("speed", greatest_info.suite_filter))
  {
    RUN_SUITE(speed_tests);
  }
#else
  printf("10-bit tests are not yet supported\n");
#endif //KVZ_BIT_DEPTH == 8

  RUN_SUITE(coeff_sum_tests);

  RUN_SUITE(mv_cand_tests);

  // Doesn't work in git
  //RUN_SUITE(inter_recon_bipred_tests);

  GREATEST_MAIN_END();
}
