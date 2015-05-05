#!/bin/sh
set -ev

if [ -z "$VALGRIND_TEST" ]; then
  cd src
  make
  make tests
else
  cd src
  make debug
  valgrind --leak-check=full --error-exitcode=1 ./kvazaar_debug -i ../mandelbrot_${TEST_DIM}.yuv --input-res=${TEST_DIM} -o /dev/null $VALGRIND_TEST
fi
