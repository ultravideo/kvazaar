#!/bin/bash

# Test GOP, with and without OWF.

set -eu
source util.bash

common_args='--gop=8 -p0 --threads=2 --wpp --rd=0 --no-rdoq --no-deblock --no-sao --no-signhide --subme=0 --pu-depth-inter=1-3 --pu-depth-intra=2-3'
valgrind_test 264x130 10 $common_args --owf=1
valgrind_test 264x130 10 $common_args --owf=4
valgrind_test 264x130 20 $common_args --owf=0
