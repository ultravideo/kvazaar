#!/bin/sh

# Test GOP, with and without OWF.

set -eu
. "${0%/*}/util.sh"

common_args='--preset=veryslow -p0 --threads=2 --wpp'

# Under valgrind use a lighter preset: valgrind is for memory checking, and
# veryslow under valgrind is ~50x slower (the valgrind job was taking 1.5 h
# with the veryslow tests).
if [ "${KVZ_TEST_VALGRIND:-0}" = '1' ]; then
  common_args='--preset=fast -p0 --threads=2 --wpp'
fi

# Skip valgrind tests if KVZ_TEST_VALGRIND is not set to 1
if [ "${KVZ_TEST_CHROMASUBSAMPLING:-0}" != '0' ]; then

  valgrind_test_422 264x130 10 $common_args --gop=8 -p0 --owf=1
  valgrind_test_422 264x130 10 $common_args --gop=16 -p0 --owf=1

  # Per-CTU QP changes (cu_qp_delta) with 4:2:2 chroma: intra-only with WPP,
  # variance adaptive QP and rate control.
  valgrind_test_422 264x130 10 $common_args -p1 --owf=1 -q34 --vaq=5
  valgrind_test_422 264x130 10 $common_args --gop=8 --owf=1 -q34 --vaq=5
  valgrind_test_422 264x130 10 $common_args --gop=8 --owf=1 --bitrate=300000

  # Do more extensive tests in a private gitlab CI runner. The 40-frame
  # veryslow encodes dominate the valgrind job's runtime, so skip them there.
  if [ ! -z ${GITLAB_CI+x} ] && [ "${KVZ_TEST_VALGRIND:-0}" != '1' ];then valgrind_test_422 264x130 40 $common_args --gop=8 -p32 --owf=4 --no-open-gop; fi

fi