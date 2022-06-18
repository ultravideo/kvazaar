#!/usr/bin/bash
set -e

# Install build dependencies for kvazaar
pacman -S --noconfirm --noprogressbar --needed \
    $MINGW_PACKAGE_PREFIX-gcc \
    $MINGW_PACKAGE_PREFIX-yasm \
	$MINGW_PACKAGE_PREFIX-autotools

# Delete unused packages to reduce space used in the Appveyor cache
pacman -Sc --noconfirm
