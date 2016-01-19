#!/bin/sh
set -ev

./autogen.sh
./configure
make --jobs=2

if [ -n "$VALGRIND_TEST" ]; then
  libtool execute valgrind --leak-check=full --error-exitcode=1 -- \
    src/kvazaar -i mandelbrot_${TEST_DIM}.yuv --input-res=${TEST_DIM} -o /dev/null \
    $VALGRIND_TEST
elif [ -n "$EXPECTED_STATUS" ]; then
  set +e
  libtool execute src/kvazaar $PARAMS
  EXIT_STATUS=$?
  set -e
  [ "$EXIT_STATUS" = "$EXPECTED_STATUS" ]
else
  make check
fi
