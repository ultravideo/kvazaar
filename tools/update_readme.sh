#!/bin/bash

# This file is part of Kvazaar HEVC encoder.
#
# Copyright (C) 2013-2016 Tampere University of Technology and others (see
# COPYING file).
#
# Kvazaar is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License version 2.1 as
# published by the Free Software Foundation.
#
# Kvazaar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.

if [[ $# != 1 ]]; then
    printf "Usage: $0 README.md\n"
    exit 1
fi

tmpfile="$(mktemp)"

{
    sed '/BEGIN KVAZAAR HELP MESSAGE/q' -- "$1";
    printf '```\n';
    kvazaar --help;
    printf '```\n';
    sed -n '/END KVAZAAR HELP MESSAGE/{:a;p;n;ba}' -- "$1";
} >> "$tmpfile"

mv -- "$tmpfile" "$1"
