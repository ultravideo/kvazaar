#!/bin/sh

# Test scalability.

set -eu
. "${0%/*}/util.sh"

# Test encoder control initing failures when using layers. TODO: strictly necessary?
valgrind_expected_test 264x130 10 1 --ilr=1
valgrind_expected_test 264x130 10 1 --no-wpp --gop=0 --layer --ilr=2
valgrind_expected_test 264x130 10 1 --no-wpp --gop=0 --layer --ilr=2 --layer --ilr=3