#!/bin/sh
set -ev

if [ -n "$VALGRIND_TEST" ]; then
  wget http://ultravideo.cs.tut.fi/ffmpeg-release-32bit-static.tar.xz
  7z x ffmpeg-release-32bit-static.tar.xz
  7z x ffmpeg-release-32bit-static.tar
  chmod +x ./ffmpeg-2.6.3-32bit-static/ffmpeg
  ./ffmpeg-2.6.3-32bit-static/ffmpeg -f lavfi -i "mandelbrot=size=${TEST_DIM}:end_pts=10" -vframes $TEST_FRAMES -pix_fmt yuv420p mandelbrot_${TEST_DIM}.yuv
fi
