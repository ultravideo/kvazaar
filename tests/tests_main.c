#include "greatest/greatest.h"


GREATEST_MAIN_DEFS();
extern SUITE(sad_tests);
extern SUITE(intra_sad_tests);
//extern SUITE(satd_tests);

int main(int argc, char **argv)
{
  GREATEST_MAIN_BEGIN();
  RUN_SUITE(sad_tests);
  RUN_SUITE(intra_sad_tests);
  //RUN_SUITE(satd_tests);
  GREATEST_MAIN_END();
}
