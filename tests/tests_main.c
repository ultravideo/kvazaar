#include "greatest/greatest.h"

#include "test_strategies.h"

GREATEST_MAIN_DEFS();
extern SUITE(sad_tests);
extern SUITE(intra_sad_tests);
extern SUITE(satd_tests);
extern SUITE(speed_tests);
extern SUITE(dct_tests);

int main(int argc, char **argv)
{
  GREATEST_MAIN_BEGIN();

  init_test_strategies();

  RUN_SUITE(sad_tests);
  RUN_SUITE(intra_sad_tests);
  RUN_SUITE(satd_tests);
  RUN_SUITE(dct_tests);

  if (greatest_info.suite_filter &&
      greatest_name_match("speed", greatest_info.suite_filter))
  {
    RUN_SUITE(speed_tests);
  }

  GREATEST_MAIN_END();
}
