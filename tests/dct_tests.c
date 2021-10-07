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
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "greatest/greatest.h"

#include "test_strategies.h"

#include "src/image.h"

#include <math.h>
#include <stdlib.h>


//////////////////////////////////////////////////////////////////////////
// MACROS
#define NUM_TESTS 12
#define NUM_SIZES 5
#define LCU_MAX_LOG_W 5
#define LCU_MIN_LOG_W 2

//////////////////////////////////////////////////////////////////////////
// GLOBALS
static int16_t * dct_bufs[NUM_TESTS] = { 0 }; // SIMD aligned pointers.
static int16_t * dct_actual_bufs[NUM_TESTS] = { 0 }; // pointers returned by malloc.

static int16_t dct_result[NUM_SIZES][LCU_WIDTH*LCU_WIDTH] = { { 0 } };
static int16_t idct_result[NUM_SIZES][LCU_WIDTH*LCU_WIDTH] = { { 0 } };

static struct test_env_t {
  int log_width; // for selecting dim from bufs
  dct_func * tested_func;
  const strategy_t * strategy;
  char msg[1024];
} test_env;


//////////////////////////////////////////////////////////////////////////
// SETUP, TEARDOWN AND HELPER FUNCTIONS
static void init_gradient(int x_px, int y_px, int width, int slope, int16_t *buf)
{
  for (int y = 0; y < width; ++y) {
    for (int x = 0; x < width; ++x) {
      int diff_x = x_px - x;
      int diff_y = y_px - y;
      int val = slope * sqrt(diff_x * diff_x + diff_y * diff_y) + 0.5;
      buf[y * width + x] = CLIP(0, 255, val);
    }
  }
}


static void setup_tests()
{
  for (int test = 0; test < NUM_TESTS; ++test) {

    dct_actual_bufs[test] = malloc(LCU_WIDTH*LCU_WIDTH*sizeof(int16_t) + SIMD_ALIGNMENT);
    dct_bufs[test] = ALIGNED_POINTER(dct_actual_bufs[test], SIMD_ALIGNMENT);
  }

  for (int test = 0; test < NUM_TESTS; ++test) {
      const int width = LCU_WIDTH;
      init_gradient(width, width, width, 255 / width, dct_bufs[test]);
  }

   

  // Select buffer width according to function name for dct function.
  int block = 0;
  for (int s = 0; s < strategies.count && block < NUM_SIZES; ++s)
  {
    strategy_t *strat = &strategies.strategies[s];
    dct_func* dct_generic = 0;
    if (
      (

      strcmp(strat->type, "fast_forward_dst_4x4") == 0 ||
      strcmp(strat->type, "dct_4x4") == 0 ||
      strcmp(strat->type, "dct_8x8") == 0 || 
      strcmp(strat->type, "dct_16x16") == 0 || 
      strcmp(strat->type, "dct_32x32") == 0
      )
      &&
      strcmp(strat->strategy_name, "generic") == 0
      )
    {
      dct_generic = strat->fptr;
      dct_generic(KVZ_BIT_DEPTH, dct_bufs[block], dct_result[block]);
      ++block;
    }
  }

  block = 0;
  for (int s = 0; s < strategies.count && block < NUM_SIZES; ++s)
  {
    strategy_t *strat = &strategies.strategies[s];
    dct_func* idct_generic = 0;
    if (
      (

      strcmp(strat->type, "fast_inverse_dst_4x4") == 0 ||
      strcmp(strat->type, "idct_4x4") == 0 ||
      strcmp(strat->type, "idct_8x8") == 0 ||
      strcmp(strat->type, "idct_16x16") == 0 ||
      strcmp(strat->type, "idct_32x32") == 0
      )
      &&
      strcmp(strat->strategy_name, "generic") == 0
      )
    {
      idct_generic = strat->fptr;
      idct_generic(KVZ_BIT_DEPTH, dct_bufs[block], idct_result[block]);
      ++block;
    }
  }
}

static void tear_down_tests()
{
  for (int test = 0; test < NUM_TESTS; ++test) {
    free(dct_actual_bufs[test]);
  }
}


//////////////////////////////////////////////////////////////////////////
// TESTS
TEST dct(void)
{
  int index = test_env.log_width - 1;
  if (strcmp(test_env.strategy->type, "fast_forward_dst_4x4") == 0) index = 0;

  int16_t *buf = dct_bufs[index];
  ALIGNED(32) int16_t test_result[LCU_WIDTH*LCU_WIDTH] = { 0 };

  test_env.tested_func(KVZ_BIT_DEPTH, buf, test_result);

  for (int i = 0; i < LCU_WIDTH*LCU_WIDTH; ++i){
    ASSERT_EQ(test_result[i], dct_result[index][i]);
  }
  
  PASS();
}

TEST idct(void)
{
  int index = test_env.log_width - 1;
  if (strcmp(test_env.strategy->type, "fast_inverse_dst_4x4") == 0) index = 0;

  int16_t *buf = dct_bufs[index];
  ALIGNED(32) int16_t test_result[LCU_WIDTH*LCU_WIDTH] = { 0 };

  test_env.tested_func(KVZ_BIT_DEPTH, buf, test_result);

  for (int i = 0; i < LCU_WIDTH*LCU_WIDTH; ++i){
    ASSERT_EQ(test_result[i], idct_result[index][i]);
  }

  PASS();
}


//////////////////////////////////////////////////////////////////////////
// TEST FIXTURES
SUITE(dct_tests)
{
  //SET_SETUP(sad_setup);
  //SET_TEARDOWN(sad_teardown);

  setup_tests();

  // Loop through all strategies picking out the intra sad ones and run
  // select strategies though all tests
  for (volatile unsigned i = 0; i < strategies.count; ++i) {
    const strategy_t * strategy = &strategies.strategies[i];

    // Select buffer width according to function name for dct function.
    if (strcmp(strategy->type, "fast_forward_dst_4x4") == 0) {
      test_env.log_width = 2;
    }
    else if (strcmp(strategy->type, "dct_4x4") == 0) {
      test_env.log_width = 2;
    }
    else if (strcmp(strategy->type, "dct_8x8") == 0) {
      test_env.log_width = 3;
    }
    else if (strcmp(strategy->type, "dct_16x16") == 0) {
      test_env.log_width = 4;
    }
    else if (strcmp(strategy->type, "dct_32x32") == 0) {
      test_env.log_width = 5;
    }
    else if (strcmp(strategy->type, "fast_inverse_dst_4x4") == 0) {
      test_env.log_width = 2;
    }
    else if (strcmp(strategy->type, "idct_4x4") == 0) {
      test_env.log_width = 2;
    }
    else if (strcmp(strategy->type, "idct_8x8") == 0) {
      test_env.log_width = 3;
    }
    else if (strcmp(strategy->type, "idct_16x16") == 0) {
      test_env.log_width = 4;
    }
    else if (strcmp(strategy->type, "idct_32x32") == 0) {
      test_env.log_width = 5;
    }
    else {
      test_env.log_width = 0;
    }

    test_env.tested_func = strategies.strategies[i].fptr;
    test_env.strategy = strategy;

    // Call different tests depending on type of function.
    // This allows for selecting a subset of tests with -t parameter.
   if (strncmp(strategy->type, "dct_", 4) == 0 ||
      strcmp(strategy->type, "fast_forward_dst_4x4") == 0)
    {
      RUN_TEST(dct);
    }
    else if (strncmp(strategy->type, "idct_", 4) == 0 ||
      strcmp(strategy->type, "fast_inverse_dst_4x4") == 0)
    {
      RUN_TEST(idct);
    }
  }

  tear_down_tests();
}
