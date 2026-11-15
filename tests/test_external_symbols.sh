#!/bin/sh

# Check for external symbols without kvz_ prefix.

set -eu${BASH+o pipefail}

# The static library of the autotools build by default. The CMake build sets
# KVZ_STATIC_LIBRARY to a static library of the same objects.
library="${KVZ_STATIC_LIBRARY:-../src/.libs/libkvazaar.a}"

if nm -go --defined-only "${library}" | grep -Ev ' _?(kvz_|__[a-z0-9]+(_|\.)get_pc_thunk\.)'; then
    printf '%s\n' 'Only symbols prefixed with "kvz_" should be exported from libkvazaar.'
    false
fi
