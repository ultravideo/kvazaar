"""This is a script that generates tables for Kvazaar HEVC encoder.

/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

This script is provided as reference, in case we ever need to change the tables
or generate more similar tables.

Because the CUs are coded in Z-order a particular index in the LCU will always
have the same number of coded reference pixels, except if the PU is on the very
top or left edge of the LCU.

"""

import numpy


def make_z_order_table(width, coord = None, zid=0, min_width=4, result=None):
    """Return a table with the quadtree z-order.

    Args:
        width: width of the area (LCU)
        coord: numpy.array with index 0 as x and 1 as y
        min_width: width at which the recursion is stopped
        result: numpy.array with the current table

    Returns: numpy.array with the quadtree z-order.
    """
    if coord is None:
        coord = numpy.array([0, 0])

    if result is None:
        num_pu = width / min_width
        result = numpy.zeros((num_pu, num_pu), numpy.int16)
    offset = width / 2

    if offset >= min_width:
        # Recurse in quadtree z-order.
        offsets = map(numpy.array, [[0,0],[1,0],[0,1],[1,1]])
        for num, os in enumerate(offsets):
            num_pu = offset**2 / min_width**2
            result = make_z_order_table(offset, coord + os * offset,
                    zid + num * num_pu, min_width, result)
    else:
        pu = coord / min_width
        result[pu[1]][pu[0]] = zid

    return result


def num_lessed_zid_on_left(table, x, y):
    """Z-order table + coord -> number of ref PUs on the left."""
    i = 0
    while True:
        if x == 0:
            return 16
        if y + i >= 16 or table[y + i][x - 1] > table[y][x]:
            return i
        i = i + 1


def num_lessed_zid_on_top(table, x, y):
    """Z-order table + coord -> number of ref PUs on the top."""
    i = 0
    while True:
        if y == 0:
            return 16
        if x + i >= 16 or table[y - 1][x + i] > table[y][x]:
            return i
        i = i + 1


def matrix_to_initializer_list(table):
    """Output a list of lists as an initializer list in C syntax.

    Args:
        table: list(list(int)) representing 2d array
    Returns:
        str
    """
    # Convert the numbers into strings and pad them to be 2-chars wide to make
    # the table look nicer.
    str_nums = (("{0: >2}".format(x) for x in line) for line in table)

    # Get the lines with all the numbers divided by commas.
    lines = (", ".join(line) for line in str_nums)

    # Join the lines with commas and newlines in between.
    result = "{ %s }" % (" },\n{ ".join(lines))

    return result


def main():
    zid_table = make_z_order_table(64)

    num_pu = 16
    left_table = numpy.zeros((num_pu, num_pu), numpy.int16)
    top_table = numpy.zeros((num_pu, num_pu), numpy.int16)

    for y in range(16):
        for x in range(16):
            left_table[y][x] = num_lessed_zid_on_left(zid_table, x, y)
            top_table[y][x] = num_lessed_zid_on_top(zid_table, x, y)

    print zid_table
    print left_table
    print top_table

    # Multiply by number of pixels in a PU
    left_table = left_table * 4
    top_table = top_table * 4

    print
    print "left"
    print matrix_to_initializer_list(left_table)
    print
    print "top"
    print matrix_to_initializer_list(top_table)
    print


if __name__ == '__main__':
    main()
