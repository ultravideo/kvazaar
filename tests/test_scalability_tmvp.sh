#!/bin/sh

# Test scalability.

set -eu
. "${0%/*}/util.sh"

# Test tmvp
valgrind_test 512x264 20 --preset=ultrafast --layer-res=256x132 --layer --preset=ultrafast --no-tmvp
valgrind_test 512x264 20 --preset=ultrafast --layer-res=256x132 --no-tmvp --layer --preset=ultrafast
valgrind_test 512x264 20 --preset=ultrafast --layer-res=256x132 --no-tmvp --layer --preset=ultrafast --no-tmvp

#   Test without owf
valgrind_test 512x264 20 --preset=ultrafast --layer-res=256x132 --layer --preset=ultrafast --no-tmvp --owf=0
valgrind_test 512x264 20 --preset=ultrafast --layer-res=256x132 --no-tmvp --layer --preset=ultrafast --owf=0
valgrind_test 512x264 20 --preset=ultrafast --layer-res=256x132 --no-tmvp --layer --preset=ultrafast --no-tmvp --owf=0
