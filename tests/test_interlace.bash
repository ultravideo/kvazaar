#!/bin/bash

set -eu
source util.bash

valgrind_test 264x130 10 --source-scan-type=tff -p0 --preset=ultrafast --threads=2 --owf=1 --wpp
