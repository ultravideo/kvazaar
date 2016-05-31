#!/bin/sh
set -ev

./autogen.sh
./configure $KVZ_CONFIGURE_ARGS
make --jobs=2 V=1

if [ -n "$VALGRIND_TEST" ]; then
  libtool execute valgrind --leak-check=full --error-exitcode=1 -- \
    src/kvazaar -i mandelbrot_${TEST_DIM}.yuv --input-res=${TEST_DIM} \
    -o test.265 $VALGRIND_TEST
    ./hmdec-16.10 -b test.265
elif [ -n "$EXPECTED_STATUS" ]; then
  set +e
  libtool execute src/kvazaar $PARAMS
  EXIT_STATUS=$?
  set -e
  [ "$EXIT_STATUS" = "$EXPECTED_STATUS" ]
else
  make check
fi
