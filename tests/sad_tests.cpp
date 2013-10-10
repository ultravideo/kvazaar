#include "seatest.h"

#include <stdlib.h>

#include "picture.h"

//////////////////////////////////////////////////////////////////////////
// EXTERNAL FUNCTIONS
unsigned get_block_sad(picture *pic, picture *ref, 
                       int pic_x, int pic_y, int ref_x, int ref_y, 
                       int block_width, int block_height);

//////////////////////////////////////////////////////////////////////////
// DEFINES
#define TEST_SAD(X, Y) get_block_sad(g_pic, g_ref, 0, 0, (X), (Y), 8, 8)

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
  g_pic = picture_init(8, 8, 1, 1);
  memcpy(g_pic->y_data, pic_data, 64);

  g_ref = picture_init(8, 8, 1, 1);
  memcpy(g_ref->y_data, ref_data, 64);
}

void sad_teardown(void)
{
  free(g_pic); g_pic = 0;
  free(g_ref); g_ref = 0;
}

//////////////////////////////////////////////////////////////////////////
// TESTS
void test_topleft(void)
{
  assert_ulong_equal(
    1*(5*5) + (2+4)*(3*5) + 5*(3*3) - 64,
    TEST_SAD(-4, -4));
}

void test_top(void)
{
  assert_ulong_equal(
    1*5 + 2*(6*5) + 3*5 + (4+6)*3 + 5*(6*3) - 64,
    TEST_SAD(0, -4));
}

void test_topright(void)
{
  assert_ulong_equal(
    3*(5*5) + (2+6)*(3*5) + 5*(3*3) - 64,
    TEST_SAD(4, -4));
}

void test_left(void)
{
  assert_ulong_equal(
    1*5 + 4*(6*5) + 7*5 + (2+8)*3 + 5*(6*3) - 64,
    TEST_SAD(-4, 0));
}

void test_no_offset(void)
{
  assert_ulong_equal(
    (1+3+7+9) + (2+4+6+8)*6 + 6*(6*5) - 64,
    TEST_SAD(0, 0));
}

void test_right(void)
{
  assert_ulong_equal(
    (3+9)*5 + 6*(5*6) + (2+8)*3 + 5*(6*3) - 64,
    TEST_SAD(4, 0));
}

void test_bottomleft(void)
{
  assert_ulong_equal(
    7*(5*5) + (4+8)*(3*5) + 5*(3*3) - 64,
    TEST_SAD(-4, 4));
}

void test_bottom(void)
{
  assert_ulong_equal(
    (7+9)*5 + 8*(6*5) + (4+6)*3 + 5*(6*3) - 64,
    TEST_SAD(0, 4));
}

void test_bottomright(void)
{
  assert_ulong_equal(
    9*(5*5) + (6+8)*(3*5) + 5*(3*3) - 64,
    TEST_SAD(-4, 4));
}

//////////////////////////////////////////////////////////////////////////
// TEST FIXTURES
void sad_tests(void)
{
  test_fixture_start();
  fixture_setup(sad_setup);

  run_test(test_topleft);
  run_test(test_top);
  run_test(test_topright);

  run_test(test_left);
  run_test(test_no_offset);
  run_test(test_right);

  run_test(test_bottomleft);
  run_test(test_bottom);
  run_test(test_bottomright);

  fixture_teardown(sad_teardown);
  test_fixture_end();
}

