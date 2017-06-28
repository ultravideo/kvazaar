#!/bin/sh

# Test trying to use invalid input dimensions.

set -eu
. util.sh

encode_test 1x65 1 1
