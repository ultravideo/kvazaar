#!/bin/bash

# Helper functions for test scripts.

set -euo pipefail

# Temporary files for encoder input and output.
yuvfile="$(mktemp --tmpdir tmp.XXXXXXXXXX.yuv)"
hevcfile="$(mktemp --tmpdir tmp.XXXXXXXXXX.hevc)"

cleanup() {
    rm -rf ${yuvfile} ${hevcfile}
}
trap cleanup EXIT

print_and_run() {
    printf '\n\n$ %s\n' "$*"
    "$@"
}

prepare() {
    cleanup
    print_and_run \
        ffmpeg -f lavfi -i "mandelbrot=size=${1}" \
            -vframes "${2}" -pix_fmt yuv420p \
            "${yuvfile}"
}

valgrind_test() {
    dimensions="$1"
    shift
    frames="$1"
    shift

    prepare "${dimensions}" "${frames}"

    print_and_run \
        libtool execute \
            valgrind --leak-check=full --error-exitcode=1 -- \
            ../src/kvazaar -i "${yuvfile}" "--input-res=${dimensions}" -o "${hevcfile}" "$@"

    print_and_run \
        TAppDecoderStatic -b "${hevcfile}"

    cleanup
}

encode_test() {
    dimensions="$1"
    shift
    frames="$1"
    shift
    expected_status="$1"
    shift

    prepare "${dimensions}" "${frames}"

    set +e
    print_and_run \
        libtool execute \
            ../src/kvazaar -i "${yuvfile}" "--input-res=${dimensions}" -o "${hevcfile}" "$@"
    actual_status="$?"
    set -e
    [[ ${actual_status} = ${expected_status} ]]
}
