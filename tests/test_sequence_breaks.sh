#!/bin/sh

# Test correct behaviour/no memory leaks when exiting from different points in the cli main func

set -eu
. "${0%/*}/util.sh"

valgrind_expected_test 264x130 10 1 --not-a-parameter
valgrind_expected_test 264x130 10 0 --version
valgrind_expected_test 264x130 10 0 --help
valgrind_expected_test 511x261 10 1
valgrind_expected_test 264x130 1 0 --seek=2