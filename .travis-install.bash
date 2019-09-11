#!/bin/bash

# Download FFmpeg and HM decoder and place them in $PATH.

set -euvo pipefail

mkdir -p "${HOME}/bin"

wget http://ultravideo.cs.tut.fi/ffmpeg-release-4.2.1-32bit-static.tar.xz
sha256sum -c - << EOF
226f55f8a94d71f3d231a20fe59fcbb7f6100cabf663f9bcb887d17b332a91c5  ffmpeg-release-4.2.1-32bit-static.tar.xz
EOF
tar xf ffmpeg-release-4.2.1-32bit-static.tar.xz
cp ffmpeg-4.2.1-i686-static/ffmpeg "${HOME}/bin/ffmpeg"
chmod +x "${HOME}/bin/ffmpeg"

wget http://ultravideo.cs.tut.fi/ubuntu-12.04-hmdec-16.10.tgz
sha256sum -c - << EOF
e00d61dd031a14aab1a03c0b23df315b8f6ec3fab66a0e2ae2162496153ccf92  ubuntu-12.04-hmdec-16.10.tgz
EOF
tar xf ubuntu-12.04-hmdec-16.10.tgz
cp hmdec-16.10 "${HOME}/bin/TAppDecoderStatic"
chmod +x "${HOME}/bin/TAppDecoderStatic"

export PATH="${HOME}/bin:${PATH}"
