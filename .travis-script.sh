#!/bin/sh
set -ev

if [ -n "$VALGRIND_TEST" ]; then
  cd src
  make debug
  valgrind --leak-check=full --error-exitcode=1 ./kvazaar_debug -i ../mandelbrot_${TEST_DIM}.yuv --input-res=${TEST_DIM} -o /dev/null $VALGRIND_TEST
elif [ -n "$EXPECTED_STATUS" ]; then
  cd src
  make cli
  set +e
  ./kvazaar $PARAMS
  EXIT_STATUS=$?
  set -e
  [ "$EXIT_STATUS" = "$EXPECTED_STATUS" ]
else
  cd src
  make
  make tests
fi
