#!/bin/sh
set -ev

ver_major=$(grep -e '^VER_MAJOR' src/Makefile | sed -e 's/VER_MAJOR//g; s/=//g; s/ //g;')
ver_minor=$(grep -e '^VER_MINOR' src/Makefile | sed -e 's/VER_MINOR//g; s/=//g; s/ //g;')
ver_release=$(grep -e '^VER_RELEASE' src/Makefile | sed -e 's/VER_RELEASE//g; s/=//g; s/ //g;')
libkvazaar=libkvazaar.so.${ver_major}.${ver_minor}.${ver_release}

if [ -n "$VALGRIND_TEST" ]; then
  cd src
  make debug
  valgrind --leak-check=full --error-exitcode=1 ./kvazaar_debug -i ../mandelbrot_${TEST_DIM}.yuv --input-res=${TEST_DIM} -o /dev/null $VALGRIND_TEST
elif [ -n "$EXPECTED_STATUS" ]; then
  cd src
  make cli
  set +e
  LD_PRELOAD=./$libkvazaar ./kvazaar $PARAMS
  EXIT_STATUS=$?
  set -e
  [ "$EXIT_STATUS" = "$EXPECTED_STATUS" ]
else
  cd src
  make
  make tests
fi
