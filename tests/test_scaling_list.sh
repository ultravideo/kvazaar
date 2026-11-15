#!/bin/sh

# Test custom scaling lists (--cqmfile) and the default scaling lists.

set -eu

. "${0%/*}/util.sh"

cqmfile="$(mktemp)"

# Writes a list of $2 entries named $1. The entries go from 1 to 40, so some
# of the quantization scales (16 * 26214 / entry) do not fit in 16 bits.
# Lists with the same seed $3 are equal, which codes them as predicted from
# an earlier list.
write_list() {
    printf '%s =\n' "$1"
    i=0
    while [ "$i" -lt "$2" ]; do
        printf '%d,' $(( 1 + (i * 7 + $3 * 13) % 40 ))
        i=$(( i + 1 ))
        [ $(( i % 8 )) -ne 0 ] || printf '\n'
    done
}

{
    for type in INTRA INTER; do
        for comp in LUMA CHROMAU CHROMAV; do
            write_list "${type}4X4_${comp}" 16 1
            write_list "${type}8X8_${comp}" 64 2
            write_list "${type}16X16_${comp}" 64 3
            printf '%s =\n%d\n' "${type}16X16_${comp}_DC" 5
        done
        write_list "${type}32X32_LUMA" 64 4
        printf '%s =\n%d\n' "${type}32X32_LUMA_DC" 30
    done
} > "${cqmfile}"

common_args='--gop=lp-g4d3t1 --threads=2 --owf=1 --pu-depth-intra=0-3 --pu-depth-inter=0-3'
valgrind_test 264x130 4 $common_args --scaling-list=custom --cqmfile="${cqmfile}"
valgrind_test 264x130 4 $common_args --scaling-list=custom --cqmfile="${cqmfile}" --rd=2 --signhide
valgrind_test 264x130 4 $common_args --scaling-list=default
# 4:4:4 uses the 32x32 chroma lists, which are those of 16x16.
valgrind_test_444 264x130 4 $common_args --scaling-list=custom --cqmfile="${cqmfile}"
valgrind_test_444 264x130 4 $common_args --scaling-list=default

rm -f "${cqmfile}"
