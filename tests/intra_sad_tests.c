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
#include "src/strategyselector.h"

#include <math.h>


//////////////////////////////////////////////////////////////////////////
// MACROS
#define NUM_TESTS 2
#define LCU_MAX_LOG_W 6
#define LCU_MIN_LOG_W 2

//////////////////////////////////////////////////////////////////////////
// GLOBALS
static kvz_pixel * bufs[NUM_TESTS][7][2];

static struct {
  int log_width; // for selecting dim from bufs
  cost_pixel_nxn_func * tested_func;
} test_env;


//////////////////////////////////////////////////////////////////////////
// SETUP, TEARDOWN AND HELPER FUNCTIONS
static void init_gradient(int x_px, int y_px, int width, int slope, kvz_pixel *buf)
{
  for (int y = 0; y < width; ++y) {
    for (int x = 0; x < width; ++x) {
      int diff_x = x_px - x;
      int diff_y = y_px - y;
      int val = sqrt(diff_x * diff_x + diff_y * diff_y) + 0.5 + slope;
      buf[y * width + x] = CLIP(0, 255, val);
    }
  }
}


static void setup_tests()
{
  for (int test = 0; test < NUM_TESTS; ++test) {
    for (int w = LCU_MIN_LOG_W; w <= LCU_MAX_LOG_W; ++w) {
      bufs[test][w][0] = 0;
      bufs[test][w][1] = 0;
    }

    for (int w = LCU_MIN_LOG_W; w <= LCU_MAX_LOG_W; ++w) {
      unsigned size = 1 << (w * 2);
      bufs[test][w][0] = malloc(size * sizeof(kvz_pixel) + SIMD_ALIGNMENT);
      bufs[test][w][0] = ALIGNED_POINTER(bufs[test][w][0], SIMD_ALIGNMENT);

      bufs[test][w][1] = malloc(size * sizeof(kvz_pixel) + SIMD_ALIGNMENT);
      bufs[test][w][1] = ALIGNED_POINTER(bufs[test][w][1], SIMD_ALIGNMENT);
    }
  }

  int test = 0;
  for (int w = LCU_MIN_LOG_W; w <= LCU_MAX_LOG_W; ++w) {
    unsigned size = 1 << (w * 2);
    FILL_ARRAY(bufs[test][w][0], 0, size);
    FILL_ARRAY(bufs[test][w][1], 255, size);
  }

  test = 1;
  for (int w = LCU_MIN_LOG_W; w <= LCU_MAX_LOG_W; ++w) {
    unsigned width = 1 << w;
    unsigned size = 1 << (w * 2);
    init_gradient(3, 1, width, 1, bufs[test][w][0]);
    //init_gradient(width / 2, 0, width, 1, bufs[test][w][1]);
    FILL_ARRAY(bufs[test][w][1], 128, size);
  }
}

static void tear_down_tests()
{
  for (int test = 0; test < NUM_TESTS; ++test) {
    for (int log_width = 2; log_width <= 6; ++log_width) {
      //free(bufs[test][log_width][0]);
      //free(bufs[test][log_width][1]);
    }
  }
}


static unsigned test_calc_sad(const kvz_pixel * buf1, const kvz_pixel * buf2, int dim)
{
  unsigned result = 0;
  for (int i = 0; i < dim * dim; ++i) {
    result += abs(buf1[i] - buf2[i]);
  }
  return result;
}


//////////////////////////////////////////////////////////////////////////
// TESTS

/**
 * Test that the maximum SAD value for a given buffer size doesn't overflow.
 */
TEST test_black_and_white(void)
{
  const int test = 0;
  const int width = 1 << test_env.log_width;

  kvz_pixel * buf1 = bufs[test][test_env.log_width][0];
  kvz_pixel * buf2 = bufs[test][test_env.log_width][1];

  unsigned result1 = test_env.tested_func(buf1, buf2);
  unsigned result2 = test_env.tested_func(buf2, buf1);

  // Order of parameters must not matter.
  ASSERT_EQ(result1, result2);

  // Result matches trivial implementation.
  ASSERT_EQ(result1, 255 * width * width);

  PASS();
}


/**
* Test that the maximum SAD value for a given buffer size doesn't overflow.
*/
TEST test_gradient(void)
{
  const int test = 1;
  const int width = 1 << test_env.log_width;

  kvz_pixel * buf1 = bufs[test][test_env.log_width][0];
  kvz_pixel * buf2 = bufs[test][test_env.log_width][1];

  unsigned result = test_calc_sad(buf1, buf2, width);
  unsigned result1 = test_env.tested_func(buf1, buf2);
  unsigned result2 = test_env.tested_func(buf2, buf1);
  
  // Order of parameters must not matter.
  ASSERT_EQ(result1, result2);

  // Result matches trivial implementation.
  ASSERT_EQ(result1, result);

  PASS();
}


//////////////////////////////////////////////////////////////////////////
// TEST FIXTURES
SUITE(intra_sad_tests)
{
  //SET_SETUP(sad_setup);
  //SET_TEARDOWN(sad_teardown);

  setup_tests();

  // Loop through all strategies picking out the intra sad ones and run
  // selectec strategies though all tests.
  for (volatile unsigned i = 0; i < strategies.count; ++i) {
    const char * type = strategies.strategies[i].type;

    if (strcmp(type, "sad_4x4") == 0) {
      test_env.log_width = 2;
    } else if (strcmp(type, "sad_8x8") == 0) {
      test_env.log_width = 3;
    } else if (strcmp(type, "sad_16x16") == 0) {
      test_env.log_width = 4;
    } else if (strcmp(type, "sad_32x32") == 0) {
      test_env.log_width = 5;
    } else if (strcmp(type, "sad_64x64") == 0) {
      test_env.log_width = 6;
    }  else {
      continue;
    }

    test_env.tested_func = strategies.strategies[i].fptr;

    // Tests
    RUN_TEST(test_black_and_white);
    RUN_TEST(test_gradient);
  }

  tear_down_tests();
}
