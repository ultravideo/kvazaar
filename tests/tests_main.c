#include "greatest/greatest.h"


GREATEST_MAIN_DEFS();
extern SUITE(sad_tests);
extern SUITE(picture_list_tests);

int main(int argc, char **argv)
{
  GREATEST_MAIN_BEGIN();
  RUN_SUITE(sad_tests);
  RUN_SUITE(picture_list_tests);
  GREATEST_MAIN_BEGIN();
}
