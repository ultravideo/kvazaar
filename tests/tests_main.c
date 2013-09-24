#include "seatest.h"

#include "picture_list_tests.h"


void all_tests(void)
{
  picture_list_tests();
}

int main()
{
  run_tests(all_tests);
}
