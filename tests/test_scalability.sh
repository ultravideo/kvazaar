#!/bin/sh

# Test scalability. TODO: Add more tests

set -eu
. "${0%/*}/util.sh"

# Basic scalable encoding
valgrind_test 512x264 20 --preset=ultrafast -p=12 --layer-res=256x132 -q=30 -r=3 --gop=0 --layer --preset=ultrafast -q=28 -r=0 --ilr=1 --gop=0 --no-wpp
valgrind_test 512x264 20 --preset=ultrafast -p=12 --layer-res=256x132 -q=30 -r=3 --gop=0 --layer --preset=ultrafast -q=28 -r=2 --ilr=1 --gop=0 --no-wpp

# Test encoder control initing failures when using layers. TODO: strictly necessary?
valgrind_expected_test 264x130 10 1 --ilr=1
valgrind_expected_test 264x130 10 1 --no-wpp --gop=0 --layer --ilr=2
valgrind_expected_test 264x130 10 1 --no-wpp --gop=0 --layer --ilr=2 --layer --ilr=3