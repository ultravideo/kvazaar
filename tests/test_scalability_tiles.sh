#!/bin/sh

# Test scalability.

set -eu
. "${0%/*}/util.sh"

#Test tiles
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --tiles 128x66 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --tiles 256x132
valgrind_test 512x264 20 --preset=ultrafast -p12 -q30 -r3 --gop=0 --tiles 256x132 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --tiles 256x132

# Test without threads
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --tiles 128x66 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --threads=0 --owf=0 --tiles 256x132
valgrind_test 512x264 20 --preset=ultrafast -p12 -q30 -r3 --gop=0 --tiles 256x132 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --threads=0 --owf=0 --tiles 256x132