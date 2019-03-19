#!/bin/sh

# Test scalability. TODO: Add more tests

set -eu
. "${0%/*}/util.sh"

# Basic scalable encoding
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r0 --ilr=1 --gop=0 --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --no-wpp --layer --preset=ultrafast -q28 -r0 --ilr=1 --gop=0 --no-wpp --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --no-wpp --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --no-wpp --no-tmvp

#   Test without threads
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r0 --ilr=1 --gop=0 --threads=0 --owf=0 --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --threads=0 --owf=0 --no-tmvp

# Test SNR
valgrind_test 512x264 20 --preset=ultrafast -p12 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r0 --ilr=1 --gop=0 --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 -q30 -r3 --gop=0 --no-wpp --layer --preset=ultrafast -q28 -r0 --ilr=1 --gop=0 --no-wpp --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 -q30 -r3 --gop=0 --no-wpp --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --no-wpp --no-tmvp

#   Test without threads
valgrind_test 512x264 20 --preset=ultrafast -p12 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r0 --ilr=1 --gop=0 --threads=0 --owf=0 --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 -q30 -r3 --gop=0 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=0 --threads=0 --owf=0 --no-tmvp

# Test with GOP
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p16 --layer-res=256x132 -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=8 --no-tmvp

#   Test without threads
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --threads=0 --owf=0 --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p16 --layer-res=256x132 -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --gop=8 --threads=0 --owf=0 --no-tmvp

# Test SAO/deblock
valgrind_test 512x264 20 --preset=ultrafast -p12  --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1  --sao=off --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1  --sao=off --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12  --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1  --sao=full --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1  --sao=full --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12   --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=off --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132  --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=off --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12   --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=full --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132  --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=full --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12  --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=off --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=off --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12  --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=full --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast -p12 --layer-res=256x132 --no-deblock --sao=full -q30 -r3 --layer --preset=ultrafast -q28 -r2 --ilr=1 --no-deblock --sao=full --no-tmvp


# Test encoder control initing failures when using layers. TODO: strictly necessary?
valgrind_expected_test 264x130 10 1 --ilr=1
valgrind_expected_test 264x130 10 1 --no-wpp --gop=0 --layer --ilr=2
valgrind_expected_test 264x130 10 1 --no-wpp --gop=0 --layer --ilr=2 --layer --ilr=3