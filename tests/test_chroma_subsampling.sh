#!/bin/sh

# Test GOP, with and without OWF.

set -eu
. "${0%/*}/util.sh"

common_args='--preset=veryslow -p0 --threads=2 --wpp --no-smp --pu-depth-intra=1-3'
valgrind_test_444 264x130 10 $common_args --gop=8 -p0 --owf=1
valgrind_test_444 264x130 20 $common_args --gop=8 -p16 --owf=0

valgrind_test_444 264x130 10 $common_args --gop=16 -p0 --owf=1

# Do more extensive tests in a private gitlab CI runner
if [ ! -z ${GITLAB_CI+x} ];then valgrind_test_444 264x130 40 $common_args --gop=8 -p32 --owf=4 --no-open-gop; fi

