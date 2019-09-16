#!/bin/sh

# Test scalability.

set -eu
. "${0%/*}/util.sh"

# Test SAO/deblock
valgrind_test 512x264 20 --preset=ultrafast -p12  --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1  --sao=off 
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1  --sao=off
valgrind_test 512x264 20 --preset=ultrafast -p12  --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1  --sao=full --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1  --sao=full
valgrind_test 512x264 20 --preset=ultrafast -p12   --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=off --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132  --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=off
valgrind_test 512x264 20 --preset=ultrafast -p12   --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=full --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132  --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=full
valgrind_test 512x264 20 --preset=ultrafast -p12  --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=off
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=off
valgrind_test 512x264 20 --preset=ultrafast -p12  --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=full
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=full
