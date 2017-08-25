#!/bin/sh

# Test temporal scalability

set -eu
. "${0%/*}/util.sh"

# With gop 8
temporal_test 264x130 30 0 --preset=ultrafast --gop=t8
temporal_test 264x130 30 1 --preset=ultrafast --gop=t8
temporal_test 264x130 30 2 --preset=ultrafast --gop=t8
temporal_test 264x130 30 3 --preset=ultrafast --gop=t8
temporal_test 264x130 30 4 --preset=ultrafast --gop=t8

# With lt gop
temporal_test 264x130 30 0 --preset=ultrafast --gop=lpt-d1t5
temporal_test 264x130 30 1 --preset=ultrafast --gop=lpt-d1t5
temporal_test 264x130 30 2 --preset=ultrafast --gop=lpt-d1t5
temporal_test 264x130 30 3 --preset=ultrafast --gop=lpt-d1t5
temporal_test 264x130 30 4 --preset=ultrafast --gop=lpt-d1t5
temporal_test 264x130 30 5 --preset=ultrafast --gop=lpt-d1t5

temporal_test 264x130 30 0 --preset=ultrafast --gop=lpt-d1t4
temporal_test 264x130 30 1 --preset=ultrafast --gop=lpt-d1t4
temporal_test 264x130 30 2 --preset=ultrafast --gop=lpt-d1t4
temporal_test 264x130 30 3 --preset=ultrafast --gop=lpt-d1t4
temporal_test 264x130 30 4 --preset=ultrafast --gop=lpt-d1t4

temporal_test 264x130 30 0 --preset=ultrafast --gop=lpt-d1t3
temporal_test 264x130 30 1 --preset=ultrafast --gop=lpt-d1t3
temporal_test 264x130 30 2 --preset=ultrafast --gop=lpt-d1t3
temporal_test 264x130 30 3 --preset=ultrafast --gop=lpt-d1t3

temporal_test 264x130 30 0 --preset=ultrafast --gop=lpt-d1t2
temporal_test 264x130 30 1 --preset=ultrafast --gop=lpt-d1t2
temporal_test 264x130 30 2 --preset=ultrafast --gop=lpt-d1t2

temporal_test 264x130 30 0 --preset=ultrafast --gop=lpt-d1t1
temporal_test 264x130 30 1 --preset=ultrafast --gop=lpt-d1t1