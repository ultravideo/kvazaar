#!/bin/sh

# Test scalability.

set -eu
. "${0%/*}/util.sh"

# Basic scalable encoding
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --no-deblock --layer --preset=ultrafast -q28 -r0 --ilr=1 --gop=0
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --no-deblock --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --no-deblock --no-wpp --layer --preset=ultrafast -q28 -r0 --ilr=1 --gop=0 --no-wpp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --no-deblock --no-wpp --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --no-wpp

#   Test without threads
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r0 --ilr=1 --gop=0 --threads=0 --owf=0
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --threads=0 --owf=0
