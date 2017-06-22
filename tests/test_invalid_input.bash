#!/bin/bash

# Test trying to use invalid input dimensions.

set -eu
source util.bash

encode_test 1x65 1 1
