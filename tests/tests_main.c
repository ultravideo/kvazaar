#include "greatest/greatest.h"


GREATEST_MAIN_DEFS();
extern SUITE(sad_tests);

int main(int argc, char **argv)
{
  GREATEST_MAIN_BEGIN();
  RUN_SUITE(sad_tests);
  GREATEST_MAIN_END();
}
