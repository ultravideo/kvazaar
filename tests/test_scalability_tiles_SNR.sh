#!/bin/sh

# Test scalability.

set -eu
. "${0%/*}/util.sh"

#Test tiles
valgrind_test 512x512 20 --preset=ultrafast -p12 -q30 -r3 --gop=0 --tiles=2x2 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --tiles=2x2

# Test without threads
valgrind_test 512x512 20 --preset=ultrafast -p12 -q30 -r3 --gop=0 --tiles=2x2 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --threads=0 --owf=0 --tiles=2x2