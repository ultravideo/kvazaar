#!/bin/bash

set -eu
source util.bash

valgrind_test 512x256 10 --threads=2 --owf=1 --preset=ultrafast --tiles=2x2 --slices=tiles
valgrind_test 264x130 10 --threads=2 --owf=1 --preset=ultrafast --slices=wpp
