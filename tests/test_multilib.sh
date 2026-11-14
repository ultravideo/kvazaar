#!/bin/sh

# Test every bit depth compiled into kvazaar (see --with-bit-depths and
# KVZ_BIT_DEPTHS).

set -eu

. "${0%/*}/util.sh"

bit_depths="$(../libtool execute ../src/kvazaar --version | sed -n 's/^Bit depths://p')"
printf 'Bit depths:%s\n' "${bit_depths}"
[ -n "${bit_depths}" ]

for bit_depth in ${bit_depths}; do
    valgrind_test 264x130 4 --output-bitdepth=${bit_depth} -p0 --threads=2 --owf=1 --gop=lp-g4d3t1
done

# A bit depth that is not compiled in is an error.
encode_test 264x130 1 1 --output-bitdepth=9

# By default, 10-bit input is coded with the lowest bit depth that has at
# least 10 bits.
expected=''
for bit_depth in ${bit_depths}; do
    expected="${bit_depth}"
    [ "${bit_depth}" -lt 10 ] || break
done
cleanup
print_and_run \
    ffmpeg -f lavfi -i "mandelbrot=size=264x130" -vframes 2 -pix_fmt yuv420p10le \
        -f rawvideo "${yuvfile}"
output="$(print_and_run \
    ../libtool execute \
        ../src/kvazaar -i "${yuvfile}" --input-res=264x130 --input-bitdepth=10 \
            -o "${hevcfile}" 2>&1)"
printf '%s\n' "${output}"
printf '%s\n' "${output}" | grep "Bit depth: ${expected} "
print_and_run \
    TAppDecoderStatic -b "${hevcfile}"
