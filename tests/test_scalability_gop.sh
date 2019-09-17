#!/bin/sh

# Test scalability.

set -eu
. "${0%/*}/util.sh"

# Test with GOP
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1
valgrind_test 512x264 20 --preset=ultrafast -p16 --layer-res=256x132 -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=8

#   Test without threads
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --threads=0 --owf=0
valgrind_test 512x264 20 --preset=ultrafast -p16 --layer-res=256x132 -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=8 --threads=0 --owf=0