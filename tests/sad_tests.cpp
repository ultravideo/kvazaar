#include "seatest.h"

#include <stdlib.h>

#include "picture.h"

//////////////////////////////////////////////////////////////////////////
// EXTERNAL FUNCTIONS
unsigned calc_sad(picture *pic, picture *ref, 
                  int pic_x, int pic_y, int ref_x, int ref_y, 
                  int block_width, int block_height);

//////////////////////////////////////////////////////////////////////////
// DEFINES
#define TEST_SAD(X, Y) calc_sad(g_pic, g_ref, 0, 0, (X), (Y), 8, 8)

//////////////////////////////////////////////////////////////////////////
// GLOBALS
const uint8_t ref_data[64] = {
  1,2,2,2,2,2,2,3,
  4,5,5,5,5,5,5,6,
  4,5,5,5,5,5,5,6,
  4,5,5,5,5,5,5,6,
  4,5,5,5,5,5,5,6,
  4,5,5,5,5,5,5,6,
  4,5,5,5,5,5,5,6,
  7,8,8,8,8,8,8,9
};

const uint8_t pic_data[64] = {
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1,
  1,1,1,1,1,1,1,1
};

picture *g_pic = 0;
picture *g_ref = 0;

//////////////////////////////////////////////////////////////////////////
// SETUP, TEARDOWN AND HELPER FUNCTIONS
void sad_setup(void)
{
  unsigned i;
  g_pic = picture_init(8, 8, 1, 1);
  for (i = 0; i < 64; ++i) {
    g_pic->y_data[i] = pic_data[i] + 48;
  }

  g_ref = picture_init(8, 8, 1, 1);
  for (i = 0; i < 64; ++i) {
    g_ref->y_data[i] = ref_data[i] + 48;
  }
}

void sad_teardown(void)
{
  free(g_pic); g_pic = 0;
  free(g_ref); g_ref = 0;
}

//////////////////////////////////////////////////////////////////////////
// OVERLAPPING BOUNDARY TESTS
void test_topleft(void)
{
  assert_ulong_equal(
    1*(4*4) + (2+4)*(4*4) + 5*(4*4) - 64,
    TEST_SAD(-3, -3));
}

void test_top(void)
{
  assert_ulong_equal(
    (1+3)*4 + 2*(6*4) + (4+6)*4 + 5*(6*4) - 64,
    TEST_SAD(0, -3));
}

void test_topright(void)
{
  assert_ulong_equal(
    3*(4*4) + (2+6)*(4*4) + 5*(4*4) - 64,
    TEST_SAD(3, -3));
}

void test_left(void)
{
  assert_ulong_equal(
    (1+7)*4 + 4*(6*4) + (2+8)*4 + 5*(6*4) - 64,
    TEST_SAD(-3, 0));
}

void test_no_offset(void)
{
  assert_ulong_equal(
    (1+3+7+9) + (2+4+6+8)*6 + 5*(6*6) - 64,
    TEST_SAD(0, 0));
}

void test_right(void)
{
  assert_ulong_equal(
    (3+9)*4 + 6*(4*6) + (2+8)*4 + 5*(6*4) - 64,
    TEST_SAD(3, 0));
}

void test_bottomleft(void)
{
  assert_ulong_equal(
    7*(4*4) + (4+8)*(4*4) + 5*(4*4) - 64,
    TEST_SAD(-3, 3));
}

void test_bottom(void)
{
  assert_ulong_equal(
    (7+9)*4 + 8*(6*4) + (4+6)*4 + 5*(6*4) - 64,
    TEST_SAD(0, 3));
}

void test_bottomright(void)
{
  assert_ulong_equal(
    9*(4*4) + (6+8)*(4*4) + 5*(4*4) - 64,
    TEST_SAD(3, 3));
}

//////////////////////////////////////////////////////////////////////////
// OUT OF FRAME TESTS

void test_topleft_out(void)
{
  assert_ulong_equal(
    1*(8*8) - 64,
    TEST_SAD(-8, -8));
}

void test_top_out(void)
{
  assert_ulong_equal(
    (1+3)*8 + 2*(6*8) - 64,
    TEST_SAD(0, -8));
}

void test_topright_out(void)
{
  assert_ulong_equal(
    3*(8*8) - 64,
    TEST_SAD(8, -8));
}

void test_left_out(void)
{
  assert_ulong_equal(
    (1+7)*8 + 4*(6*8) - 64,
    TEST_SAD(-8, 0));
}

void test_right_out(void)
{
  assert_ulong_equal(
    (3+9)*8 + 6*(6*8) - 64,
    TEST_SAD(8, 0));
}

void test_bottomleft_out(void)
{
  assert_ulong_equal(
    7*(8*8) - 64,
    TEST_SAD(-8, 8));
}

void test_bottom_out(void)
{
  assert_ulong_equal(
    (7+9)*8 + 8*(6*8) - 64,
    TEST_SAD(0, 8));
}

void test_bottomright_out(void)
{
  assert_ulong_equal(
    9*(8*8) - 64,
    TEST_SAD(8, 8));
}

//////////////////////////////////////////////////////////////////////////
// TEST FIXTURES
void sad_tests(void)
{
  test_fixture_start();
  fixture_setup(sad_setup);


  // Tests for movement vectors that overlap frame.
  run_test(test_topleft);
  run_test(test_top);
  run_test(test_topright);

  run_test(test_left);
  run_test(test_no_offset);
  run_test(test_right);

  run_test(test_bottomleft);
  run_test(test_bottom);
  run_test(test_bottomright);

  // Tests for movement vectors that are outside the frame.
  run_test(test_topleft_out);
  run_test(test_top_out);
  run_test(test_topright_out);

  run_test(test_left_out);
  run_test(test_right_out);

  run_test(test_bottomleft_out);
  run_test(test_bottom_out);
  run_test(test_bottomright_out);


  fixture_teardown(sad_teardown);
  test_fixture_end();
}

