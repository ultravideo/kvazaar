#!/bin/sh

# Test scalability. TODO: Add more tests

set -eu
. "${0%/*}/util.sh"

# Basic scalable encoding
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r0 --ilr=1 --gop=0
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0

# Test encoder control initing failures when using layers. TODO: strictly necessary?
valgrind_expected_test 264x130 10 1 --ilr=1
valgrind_expected_test 264x130 10 1 --no-wpp --gop=0 --layer --ilr=2
valgrind_expected_test 264x130 10 1 --no-wpp --gop=0 --layer --ilr=2 --layer --ilr=3