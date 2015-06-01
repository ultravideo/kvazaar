#!/bin/sh
set -ev

if [ -n "$USE_NEW_GCC" ]; then
  sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y > /dev/null
  sudo apt-get update -qq
  sudo apt-get install -qq gcc-4.8
  export CC=gcc-4.8
else
  sudo apt-get update -qq
fi

if [ -n "$VALGRIND_TEST" ]; then
  sudo apt-get install -qq valgrind
  sudo apt-get install -qq p7zip-full
  wget http://ultravideo.cs.tut.fi/ffmpeg-release-32bit-static.tar.xz
  7z x ffmpeg-release-32bit-static.tar.xz
  7z x ffmpeg-release-32bit-static.tar
  chmod +x ./ffmpeg-2.6.3-32bit-static/ffmpeg
  ./ffmpeg-2.6.3-32bit-static/ffmpeg -f lavfi -i "mandelbrot=size=${TEST_DIM}:end_pts=10" -vframes $TEST_FRAMES -pix_fmt yuv420p mandelbrot_${TEST_DIM}.yuv
fi

if [ -n "$USE_YASM" ]; then
  sudo apt-get install -qq yasm
fi
