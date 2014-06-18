#include "greatest/greatest.h"

#include "src/image.h"
#include "src/strategyselector.h"


//////////////////////////////////////////////////////////////////////////
// MACROS
#define NUM_TESTS 1
#define LCU_MAX_LOG_W 6
#define LCU_MIN_LOG_W 2

//////////////////////////////////////////////////////////////////////////
// GLOBALS
static strategy_list strategies;
pixel * bufs[NUM_TESTS][6][2];

static struct {
  int log_width; // for selecting dim from bufs
  cost_pixel_nxn_func * tested_func;
} test_env;


//////////////////////////////////////////////////////////////////////////
// SETUP, TEARDOWN AND HELPER FUNCTIONS
static void init_strategies()
{
  strategies.allocated = 0;
  strategies.count = 0;
  strategies.strategies = NULL;

  // Init strategyselector because it sets hardware flags.
  strategyselector_init();

  // Collect all strategies.
  if (!strategy_register_picture(&strategies)) {
    fprintf(stderr, "strategy_register_picture failed!\n");
    return;
  }
}


static void setup_tests()
{
  init_strategies();

  for (int test = 0; test < NUM_TESTS; ++test) {
    for (int w = LCU_MIN_LOG_W; w <= LCU_MAX_LOG_W; ++w) {
      bufs[test][w][0] = 0;
      bufs[test][w][1] = 0;
    }

    for (int w = LCU_MIN_LOG_W; w <= LCU_MAX_LOG_W; ++w) {
      unsigned size = 1 << (w * 2);
      bufs[test][w][0] = malloc(size * sizeof(pixel));
      bufs[test][w][1] = malloc(size * sizeof(pixel));
    }
  }

  int test = 0;
  for (int w = LCU_MIN_LOG_W; w <= LCU_MAX_LOG_W; ++w) {
    unsigned size = 1 << (w * 2);
    memset(bufs[test][w][0], 0, size);
    memset(bufs[test][w][1], 255, size);
  }
}

static void tear_down_tests()
{
  for (int test = 0; test < NUM_TESTS; ++test) {
    for (int log_width = 2; log_width <= 6; ++log_width) {
      free(bufs[test][log_width][0]);
      free(bufs[test][log_width][1]);
    }
  }
}


static unsigned test_calc_sad(const pixel * buf1, const pixel * buf2, int dim)
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

  pixel * buf1 = bufs[test][test_env.log_width][0];
  pixel * buf2 = bufs[test][test_env.log_width][1];

  unsigned result1 = test_env.tested_func(buf1, buf2);
  unsigned result2 = test_env.tested_func(buf2, buf1);

  // Order of parameters must not matter.
  ASSERT_EQ(result1, result2);

  // Result matches trivial implementation.
  ASSERT_EQ(result1, 255 * width * width);

  PASS();
}

void sad_intra_test_performance(void)
{
  const int test = 0;
  const int width = 1 << test_env.log_width;

  pixel * buf1 = bufs[test][test_env.log_width][0];
  pixel * buf2 = bufs[test][test_env.log_width][1];

  unsigned result1 = test_env.tested_func(buf1, buf2);
  
  return;
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
  for (unsigned i = 0; i < strategies.count; ++i) {
    const char * type = strategies.strategies[i].type;

    if (strcmp(type, "sad_8bit_4x4") == 0) {
      test_env.log_width = 2;
    } else if (strcmp(type, "sad_8bit_8x8") == 0) {
      test_env.log_width = 3;
    } else if (strcmp(type, "sad_8bit_16x16") == 0) {
      test_env.log_width = 4;
    } else if (strcmp(type, "sad_8bit_32x32") == 0) {
      test_env.log_width = 5;
    } else if (strcmp(type, "sad_8bit_64x64") == 0) {
      test_env.log_width = 6;
    }  else {
      continue;
    }

    test_env.tested_func = strategies.strategies[i].fptr;

    // Tests
    RUN_TEST(test_black_and_white);
    for (int i = 0; i < 100000; ++i){
      sad_intra_test_performance();
    }
  }

  tear_down_tests();
}
