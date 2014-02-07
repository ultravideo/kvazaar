#include "greatest/greatest.h"

#include "src/picture.h"


//////////////////////////////////////////////////////////////////////////
// EXTERNAL FUNCTIONS

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
void sad_setup(void *environment)
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

void sad_teardown(void *environment)
{
  free(g_pic); g_pic = 0;
  free(g_ref); g_ref = 0;
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

//////////////////////////////////////////////////////////////////////////
// TEST FIXTURES
SUITE(sad_tests)
{
  //SET_SETUP(sad_setup);
  //SET_TEARDOWN(sad_teardown);

  sad_setup(0);
  
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

  sad_setup(0);
}
