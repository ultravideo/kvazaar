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

#include "greatest/greatest.h"

#include "test_strategies.h"

#include "src/image.h"

#include <string.h>


//////////////////////////////////////////////////////////////////////////
// EXTERNAL FUNCTIONS

//////////////////////////////////////////////////////////////////////////
// DEFINES
#define TEST_SAD(X, Y) image_calc_sad(g_pic, g_ref, 0, 0, (X), (Y), 8, 8, -1)

//////////////////////////////////////////////////////////////////////////
// GLOBALS
static const kvz_pixel ref_data[64] = {
  1,2,2,2,2,2,2,3,
  4,5,5,5,5,5,5,6,
  4,5,5,5,5,5,5,6,
  4,5,5,5,5,5,5,6,
  4,5,5,5,5,5,5,6,
  4,5,5,5,5,5,5,6,
  4,5,5,5,5,5,5,6,
  7,8,8,8,8,8,8,9
};

static const kvz_pixel pic_data[64] = {
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1
};

static kvz_picture *g_pic = 0;
static kvz_picture *g_ref = 0;

//////////////////////////////////////////////////////////////////////////
// SETUP, TEARDOWN AND HELPER FUNCTIONS
static void setup_tests()
{
  g_pic = image_alloc(8, 8);
  for (int i = 0; i < 64; ++i) {
    g_pic->y[i] = pic_data[i] + 48;
  }

  g_ref = image_alloc(8, 8);
  for (int i = 0; i < 64; ++i) {
    g_ref->y[i] = ref_data[i] + 48;
  }
}

static void tear_down_tests()
{
  image_free(g_pic);
  image_free(g_ref);
}


//////////////////////////////////////////////////////////////////////////
// OVERLAPPING BOUNDARY TESTS
TEST test_topleft(void)
{
  ASSERT_EQ(
    1*(4*4) + (2+4)*(4*4) + 5*(4*4) - 64,
    TEST_SAD(-3, -3));
  PASS();
}

TEST test_top(void)
{
  ASSERT_EQ(
    (1+3)*4 + 2*(6*4) + (4+6)*4 + 5*(6*4) - 64,
    TEST_SAD(0, -3));
  PASS();
}

TEST test_topright(void)
{
  ASSERT_EQ(
    3*(4*4) + (2+6)*(4*4) + 5*(4*4) - 64,
    TEST_SAD(3, -3));
  PASS();
}

TEST test_left(void)
{
  ASSERT_EQ(
    (1+7)*4 + 4*(6*4) + (2+8)*4 + 5*(6*4) - 64,
    TEST_SAD(-3, 0));
  PASS();
}

TEST test_no_offset(void)
{
  ASSERT_EQ(
    (1+3+7+9) + (2+4+6+8)*6 + 5*(6*6) - 64,
    TEST_SAD(0, 0));
  PASS();
}

TEST test_right(void)
{
  ASSERT_EQ(
    (3+9)*4 + 6*(4*6) + (2+8)*4 + 5*(6*4) - 64,
    TEST_SAD(3, 0));
  PASS();
}

TEST test_bottomleft(void)
{
  ASSERT_EQ(
    7*(4*4) + (4+8)*(4*4) + 5*(4*4) - 64,
    TEST_SAD(-3, 3));
  PASS();
}

TEST test_bottom(void)
{
  ASSERT_EQ(
    (7+9)*4 + 8*(6*4) + (4+6)*4 + 5*(6*4) - 64,
    TEST_SAD(0, 3));
  PASS();
}

TEST test_bottomright(void)
{
  ASSERT_EQ(
    9*(4*4) + (6+8)*(4*4) + 5*(4*4) - 64,
    TEST_SAD(3, 3));
  PASS();
}

//////////////////////////////////////////////////////////////////////////
// OUT OF FRAME TESTS

#define DIST 10
TEST test_topleft_out(void)
{
  ASSERT_EQ(
    1*(8*8) - 64,
    TEST_SAD(-DIST, -DIST));
  PASS();
}

TEST test_top_out(void)
{
  ASSERT_EQ(
    (1+3)*8 + 2*(6*8) - 64,
    TEST_SAD(0, -DIST));
  PASS();
}

TEST test_topright_out(void)
{
  ASSERT_EQ(
    3*(8*8) - 64,
    TEST_SAD(DIST, -DIST));
  PASS();
}

TEST test_left_out(void)
{
  ASSERT_EQ(
    (1+7)*8 + 4*(6*8) - 64,
    TEST_SAD(-DIST, 0));
  PASS();
}

TEST test_right_out(void)
{
  ASSERT_EQ(
    (3+9)*8 + 6*(6*8) - 64,
    TEST_SAD(DIST, 0));
  PASS();
}

TEST test_bottomleft_out(void)
{
  ASSERT_EQ(
    7*(8*8) - 64,
    TEST_SAD(-DIST, DIST));
  PASS();
}

TEST test_bottom_out(void)
{
  ASSERT_EQ(
    (7+9)*8 + 8*(6*8) - 64,
    TEST_SAD(0, DIST));
  PASS();
}

TEST test_bottomright_out(void)
{
  ASSERT_EQ(
    9*(8*8) - 64,
    TEST_SAD(DIST, DIST));
  PASS();
}


struct sad_test_env_t {
  kvz_picture *g_pic;
  kvz_picture *g_ref;
};


//////////////////////////////////////////////////////////////////////////
// TEST FIXTURES
SUITE(sad_tests)
{
  //SET_SETUP(sad_setup);
  //SET_TEARDOWN(sad_teardown);

  setup_tests();

  for (unsigned i = 0; i < strategies.count; ++i) {
    if (strcmp(strategies.strategies[i].type, "reg_sad") != 0) {
      continue;
    }

    // Change the global reg_sad function pointer.
    reg_sad = strategies.strategies[i].fptr;

    // Tests for movement vectors that overlap frame.
    RUN_TEST(test_topleft);
    RUN_TEST(test_top);
    RUN_TEST(test_topright);

    RUN_TEST(test_left);
    RUN_TEST(test_no_offset);
    RUN_TEST(test_right);

    RUN_TEST(test_bottomleft);
    RUN_TEST(test_bottom);
    RUN_TEST(test_bottomright);

    // Tests for movement vectors that are outside the frame.
    RUN_TEST(test_topleft_out);
    RUN_TEST(test_top_out);
    RUN_TEST(test_topright_out);

    RUN_TEST(test_left_out);
    RUN_TEST(test_right_out);

    RUN_TEST(test_bottomleft_out);
    RUN_TEST(test_bottom_out);
    RUN_TEST(test_bottomright_out);
  }
  
  tear_down_tests();
}
