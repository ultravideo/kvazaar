#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>

#include "seatest.h"

void test1(void)
{
  assert_bit_set(0, 0x01);
  assert_bit_set(2, 0x04);
}

void test2(void)
{
  assert_bit_set(0, 0x00);
}

void fixture1(void)
{
  test_fixture_start();
  run_test(test1);
  run_test(test2);
  test_fixture_end();
}

void all_tests(void)
{
  fixture1();
}

int main()
{
  run_tests(all_tests);
}
