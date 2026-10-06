#!/bin/sh

# Test the lookahead: scene cuts, CU-tree, CRF and ABR rate control.

set -eu
. "${0%/*}/util.sh"

common_args='--threads=2 --rd=0 --no-rdoq --no-deblock --no-sao --no-signhide --subme=0 --pu-depth-inter=1-3 --pu-depth-intra=2-3'
valgrind_test 264x130 20 $common_args --gop=8 -p16 --owf=1 --lookahead=8
valgrind_test 264x130 20 $common_args --gop=8 -p0 --owf=4 --crf=30
valgrind_test 264x130 20 $common_args --gop=0 -p8 --owf=1 --crf=30 --lookahead=4
valgrind_test 264x130 20 $common_args --gop=lp-g4d3t1 --owf=2 --crf=30
valgrind_test 264x130 20 $common_args --gop=8 --owf=1 --bitrate=500000 --rc-algorithm=abr
valgrind_test 264x130 20 $common_args --gop=8 --owf=1 --bitrate=500000 --lookahead=10
valgrind_test 264x130 5 $common_args --gop=8 --owf=1 --crf=30 --lookahead=20
valgrind_test 264x130 20 $common_args --gop=8 --owf=2 --crf=26 --vbv-maxrate=300000 --vbv-bufsize=150000
valgrind_test 264x130 20 $common_args --gop=0 --owf=1 --bitrate=300000 --rc-algorithm=abr --vbv-maxrate=400000 --vbv-bufsize=200000
if [ ! -z ${GITLAB_CI+x} ];then valgrind_test 264x130 40 $common_args --gop=16 -p32 --owf=4 --crf=28 --no-open-gop; fi
if [ ! -z ${GITLAB_CI+x} ];then valgrind_test 264x130 40 $common_args --gop=8 -p0 --owf=0 --bitrate=200000 --rc-algorithm=abr --vaq=5; fi
