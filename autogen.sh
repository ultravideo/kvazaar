#!/bin/sh

git submodule init
git submodule update
mkdir -p m4
autoreconf -if
