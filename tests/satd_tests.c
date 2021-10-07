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

#include "src/image.h"

#include <math.h>

//////////////////////////////////////////////////////////////////////////
// MACROS
#define NUM_TESTS 3
#define LCU_MAX_LOG_W 6
#define LCU_MIN_LOG_W 2

//////////////////////////////////////////////////////////////////////////
// GLOBALS
static kvz_pixel * satd_bufs[NUM_TESTS][7][2];

static struct {
  int log_width; // for selecting dim from satd_bufs
  cost_pixel_nxn_func * tested_func;
} satd_test_env;


//////////////////////////////////////////////////////////////////////////
// SETUP, TEARDOWN AND HELPER FUNCTIONS
static void setup_tests()
{
  for (int test = 0; test < NUM_TESTS; ++test) {
    for (int w = 0; w <= LCU_MIN_LOG_W; ++w) {
      satd_bufs[test][w][0] = NULL;
      satd_bufs[test][w][1] = NULL;
    }

    for (int w = LCU_MIN_LOG_W; w <= LCU_MAX_LOG_W; ++w) {
      unsigned size = 1 << (w * 2);
      satd_bufs[test][w][0] = malloc(size * sizeof(kvz_pixel));
      satd_bufs[test][w][1] = malloc(size * sizeof(kvz_pixel));
    }
  }

  //Black and white buffers
  int test = 0;
  for (int w = LCU_MIN_LOG_W; w <= LCU_MAX_LOG_W; ++w) {
    unsigned size = 1 << (w * 2);
    FILL_ARRAY(satd_bufs[test][w][0], 0, size);
    FILL_ARRAY(satd_bufs[test][w][1], 255, size);
  }

  //Checker patterns, buffer 1 is negative of buffer 2
  test = 1;
  for (int w = LCU_MIN_LOG_W; w <= LCU_MAX_LOG_W; ++w) {
    unsigned size = 1 << (w * 2);
    for (int i = 0; i < size; ++i){
      satd_bufs[test][w][0][i] = 255 * ( ( ((i >> w)%2) + (i % 2) ) % 2);
      satd_bufs[test][w][1][i] = (satd_bufs[test][w][0][i] + 1) % 2 ;
    }
  }

  //Gradient test pattern
  test = 2;
  for (int w = LCU_MIN_LOG_W; w <= LCU_MAX_LOG_W; ++w) {
    unsigned size = 1 << (w * 2);
    for (int i = 0; i < size; ++i){
      int column = (i % (1 << w) );
      int row = (i / (1 << w) );
      int r = sqrt(row * row + column * column);
      satd_bufs[test][w][0][i] = 255 / (r + 1);
      satd_bufs[test][w][1][i] = 255 - 255 / (r + 1);
    }
  }
}

static void satd_tear_down_tests()
{
  for (int test = 0; test < NUM_TESTS; ++test) {
    for (int log_width = 2; log_width <= 6; ++log_width) {
      free(satd_bufs[test][log_width][0]);
      free(satd_bufs[test][log_width][1]);
    }
  }
}


//////////////////////////////////////////////////////////////////////////
// TESTS

TEST satd_test_black_and_white(void)
{
  const int satd_results[5] = {2040, 4080, 16320, 65280, 261120};
  
  const int test = 0;

  kvz_pixel * buf1 = satd_bufs[test][satd_test_env.log_width][0];
  kvz_pixel * buf2 = satd_bufs[test][satd_test_env.log_width][1];

  unsigned result1 = satd_test_env.tested_func(buf1, buf2);
  unsigned result2 = satd_test_env.tested_func(buf2, buf1);

  ASSERT_EQ(result1, result2);
  ASSERT_EQ(result1, satd_results[satd_test_env.log_width - 2]);

  PASS();
}

TEST satd_test_checkers(void)
{
  const int satd_checkers_results[5] = { 2040, 4080, 16320, 65280, 261120 };

  const int test = 1;

  kvz_pixel * buf1 = satd_bufs[test][satd_test_env.log_width][0];
  kvz_pixel * buf2 = satd_bufs[test][satd_test_env.log_width][1];
  
  unsigned result1 = satd_test_env.tested_func(buf1, buf2);
  unsigned result2 = satd_test_env.tested_func(buf2, buf1);

  ASSERT_EQ(result1, result2);
  ASSERT_EQ(result1, satd_checkers_results[satd_test_env.log_width - 2]);

  PASS();
}


TEST satd_test_gradient(void)
{
  const int satd_gradient_results[5] = {3140,9004,20481,67262,258672};

  const int test = 2;

  kvz_pixel * buf1 = satd_bufs[test][satd_test_env.log_width][0];
  kvz_pixel * buf2 = satd_bufs[test][satd_test_env.log_width][1];
  
  unsigned result1 = satd_test_env.tested_func(buf1, buf2);
  unsigned result2 = satd_test_env.tested_func(buf2, buf1);

  ASSERT_EQ(result1, result2);
  ASSERT_EQ(result1, satd_gradient_results[satd_test_env.log_width - 2]);

  PASS();
}

//////////////////////////////////////////////////////////////////////////
// TEST FIXTURES
SUITE(satd_tests)
{
  setup_tests();

  // Loop through all strategies picking out the intra sad ones and run
  // selectec strategies though all tests.
  for (volatile unsigned i = 0; i < strategies.count; ++i) {
    const char * type = strategies.strategies[i].type;
    
    if (strcmp(type, "satd_4x4") == 0) {
      satd_test_env.log_width = 2;
    }
    else if (strcmp(type, "satd_8x8") == 0) {
      satd_test_env.log_width = 3;
    }
    else if (strcmp(type, "satd_16x16") == 0) {
      satd_test_env.log_width = 4;
    }
    else if (strcmp(type, "satd_32x32") == 0) {
      satd_test_env.log_width = 5;
    }
    else if (strcmp(type, "satd_64x64") == 0) {
      satd_test_env.log_width = 6;
    }
    else {
      continue;
    }

    satd_test_env.tested_func = strategies.strategies[i].fptr;

    // Tests
    RUN_TEST(satd_test_black_and_white);
    RUN_TEST(satd_test_checkers);
    RUN_TEST(satd_test_gradient);
  }

  satd_tear_down_tests();
}
