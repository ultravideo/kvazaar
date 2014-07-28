#include "greatest/greatest.h"

#include "src/image.h"
#include "src/strategyselector.h"
#include "src/threads.h"

#include <math.h>
#include <stdlib.h>


//////////////////////////////////////////////////////////////////////////
// MACROS
#define NUM_TESTS 113
#define NUM_CHUNKS 36
#define LCU_MAX_LOG_W 6
#define LCU_MIN_LOG_W 2

// Time per tested function, in seconds.
#define TIME_PER_TEST 1.0

//////////////////////////////////////////////////////////////////////////
// GLOBALS
static strategy_list strategies;
pixel * bufs[NUM_TESTS]; // SIMD aligned pointers.
pixel * actual_bufs[NUM_TESTS]; // pointers returned by malloc.

static struct test_env_t {
  int log_width; // for selecting dim from bufs
  cost_pixel_nxn_func * tested_func;
  const strategy * strategy;
  char msg[1024];
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


static void init_gradient(int x_px, int y_px, int width, int slope, pixel *buf)
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
  init_strategies();

  for (int test = 0; test < NUM_TESTS; ++test) {
    bufs[test] = 0;
    bufs[test] = 0;
    unsigned size = 36*64*64;
    actual_bufs[test] = malloc(size * sizeof(pixel) + SIMD_ALIGNMENT);
    bufs[test] = ALIGNED_POINTER(actual_bufs[test], SIMD_ALIGNMENT);
  }

  for (int test = 0; test < NUM_TESTS; ++test) {
    for (int chunk = 0; chunk < NUM_CHUNKS; ++chunk) {
      const int width = 64;
      int x = (test + chunk) % width;
      int y = (test + chunk) / width;
      init_gradient(width - x, y, width, 255 / width, &bufs[test][chunk * 64*64]);
    }
  }
}

static void tear_down_tests()
{
  for (int test = 0; test < NUM_TESTS; ++test) {
    free(actual_bufs[test]);
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

TEST test_intra_speed(const int width)
{
  const int size = width * width;
  uint64_t call_cnt = 0;
  CLOCK_T clock_now;
  GET_TIME(&clock_now);
  double test_end = CLOCK_T_AS_DOUBLE(clock_now) + TIME_PER_TEST;

  // Loop until time allocated for test has passed.
  for (unsigned i = 0; 
      test_end > CLOCK_T_AS_DOUBLE(clock_now);
      ++i, GET_TIME(&clock_now))
  {
    int test = i % NUM_TESTS;
    uint64_t sum = 0;
    for (int offset = 0; offset < NUM_CHUNKS * 64 * 64; offset += NUM_CHUNKS * size) {
      // Compare the first chunk against the 35 other chunks to simulate real usage.
      pixel * buf1 = &bufs[test][offset];
      for (int chunk = 1; chunk < NUM_CHUNKS; ++chunk) {
        pixel * buf2 = &bufs[test][chunk * size + offset];
        sum += test_env.tested_func(buf1, buf2);
        ++call_cnt;
      }
    }

    ASSERT(sum > 0);
  }

  sprintf(test_env.msg, "%.3fM x %s:%s",
    (double)call_cnt / 1000000.0,
    test_env.strategy->type,
    test_env.strategy->strategy_name);
  PASSm(test_env.msg);
}


TEST test_inter_speed(const int width)
{
  const int size = width * width;
  unsigned call_cnt = 0;
  CLOCK_T clock_now;
  GET_TIME(&clock_now);
  double test_end = CLOCK_T_AS_DOUBLE(clock_now) + TIME_PER_TEST;

  // Loop until time allocated for test has passed.
  for (unsigned i = 0;
      test_end > CLOCK_T_AS_DOUBLE(clock_now);
      ++i, GET_TIME(&clock_now))
  {
    int test = i % NUM_TESTS;
    uint64_t sum = 0;
    for (int offset = 0; offset < NUM_CHUNKS * 64 * 64; offset += NUM_CHUNKS * size) {
      // Treat 4 consecutive chunks as one chunk with double width and height,
      // and do a 8x8 grid search against the first chunk to simulate real usage.
      pixel * buf1 = &bufs[test][offset];
      for (int chunk = 0; chunk < NUM_CHUNKS; chunk += 4) {
        pixel * buf2 = &bufs[test][chunk * size + offset];
        for (int y = 0; y < 8; ++y) {
          for (int x = 0; x < 8; ++x) {
            const int stride1 = 2 * 64;
            const int stride2 = 2 * 64;
            sum += test_env.tested_func(buf1, &buf2[y * stride2 + x], width, width, stride1, stride2);
            ++call_cnt;
          }
        }
      }
    }
    ASSERT(sum > 0);
  }

  sprintf(test_env.msg, "%.3fM x %s(%ix%i):%s",
    (double)call_cnt / 1000000.0,
    test_env.strategy->type,
    width,
    width,
    test_env.strategy->strategy_name);
  PASSm(test_env.msg);
}


TEST intra_sad(void)
{
  const int width = 1 << test_env.log_width;
  return test_intra_speed(width);
}


TEST intra_satd(void)
{
  const int width = 1 << test_env.log_width;
  return test_intra_speed(width);
}


TEST inter_sad(void)
{
  const int width = 1 << test_env.log_width;
  return test_inter_speed(width);
}


//////////////////////////////////////////////////////////////////////////
// TEST FIXTURES
SUITE(speed_tests)
{
  //SET_SETUP(sad_setup);
  //SET_TEARDOWN(sad_teardown);

  setup_tests();

  // Loop through all strategies picking out the intra sad ones and run
  // selectec strategies though all tests
  for (unsigned i = 0; i < strategies.count; ++i) {
    const strategy * strategy = &strategies.strategies[i];

    // Select buffer width according to function name for intra cost functions.
    if (strcmp(strategy->type, "sad_8bit_4x4") == 0) {
      test_env.log_width = 2;
    } else if (strcmp(strategy->type, "sad_8bit_8x8") == 0) {
      test_env.log_width = 3;
    } else if (strcmp(strategy->type, "sad_8bit_16x16") == 0) {
      test_env.log_width = 4;
    } else if (strcmp(strategy->type, "sad_8bit_32x32") == 0) {
      test_env.log_width = 5;
    } else if (strcmp(strategy->type, "sad_8bit_64x64") == 0) {
      test_env.log_width = 6;
    } else if (strcmp(strategy->type, "satd_8bit_4x4") == 0) {
      test_env.log_width = 2;
    } else if (strcmp(strategy->type, "satd_8bit_8x8") == 0) {
      test_env.log_width = 3;
    } else if (strcmp(strategy->type, "satd_8bit_16x16") == 0) {
      test_env.log_width = 4;
    } else if (strcmp(strategy->type, "satd_8bit_32x32") == 0) {
      test_env.log_width = 5;
    } else if (strcmp(strategy->type, "satd_8bit_64x64") == 0) {
      test_env.log_width = 6;
    }

    test_env.tested_func = strategies.strategies[i].fptr;
    test_env.strategy = strategy;

    // Call different tests depending on type of function.
    // This allows for selecting a subset of tests with -t parameter.
    if (strncmp(strategy->type, "satd_8bit_", 10) == 0) {
      RUN_TEST(intra_satd);
    } else if (strncmp(strategy->type, "sad_8bit_", 9) == 0) {
      RUN_TEST(intra_sad);
    } else if (strcmp(strategy->type, "reg_sad") == 0) {
      // Call reg_sad with all the sizes it is actually called with.
      for (int width = 3; width <= 6; ++width) {
        test_env.log_width = width;
        RUN_TEST(inter_sad);
      }
    }
  }

  tear_down_tests();
}
