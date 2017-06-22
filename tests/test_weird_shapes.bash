#!/bin/bash

set -eu
source util.bash

valgrind_test  16x16  10 --threads=2 --owf=1 --preset=veryslow
valgrind_test 256x16  10 --threads=2 --owf=1 --preset=veryslow
valgrind_test  16x256 10 --threads=2 --owf=1 --preset=veryslow
