/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

/*
 * \file
 */

#include <string.h>
#include <stdio.h>

#include "yuv_input.h"

static void fill_after_frame(unsigned height, unsigned array_width,
                             unsigned array_height, pixel_t *data)
{
  pixel_t* p = data + height * array_width;
  pixel_t* end = data + array_width * array_height;

  while (p < end) {
    // Fill the line by copying the line above.
    memcpy(p, p - array_width, array_width);
    p += array_width;
  }
}


static int read_and_fill_frame_data(FILE *file,
                                    unsigned width, unsigned height,
                                    unsigned array_width, pixel_t *data)
{
  pixel_t* p = data;
  pixel_t* end = data + array_width * height;
  pixel_t fill_char;
  unsigned i;

  while (p < end) {
    // Read the beginning of the line from input.
    if (width != fread(p, sizeof(unsigned char), width, file))
      return 0;

    // Fill the rest with the last pixel value.
    fill_char = p[width - 1];

    for (i = width; i < array_width; ++i) {
      p[i] = fill_char;
    }

    p += array_width;
  }
  return 1;
}


/**
 * \brief Read a single frame from a file.
 *
 * Read luma and chroma values from file. Extend pixels if the image buffer
 * is larger than the input image.
 *
 * \param file          input file
 * \param input_width   width of the input video in pixels
 * \param input_height  height of the input video in pixels
 * \param array_width   width of the image buffer in pixels
 * \param array_height  height of the image buffer in pixels
 * \param img_out       image buffer
 *
 * \return              1 on success, 0 on failure
 */
int yuv_input_read(FILE* file,
                   unsigned input_width, unsigned input_height,
                   unsigned array_width, unsigned array_height,
                   image_t *img_out)
{
  const unsigned y_size = input_width * input_height;
  const unsigned uv_input_width  = input_width  / 2;
  const unsigned uv_input_height = input_height / 2;
  const unsigned uv_size = uv_input_width * uv_input_height;

  const unsigned uv_array_width  = array_width  / 2;
  const unsigned uv_array_height = array_height  / 2;

  if (input_width == array_width) {
    // No need to extend pixels.
    const size_t pixel_size = sizeof(unsigned char);
    if (fread(img_out->y, pixel_size, y_size,  file) != y_size)  return 0;
    if (fread(img_out->u, pixel_size, uv_size, file) != uv_size) return 0;
    if (fread(img_out->v, pixel_size, uv_size, file) != uv_size) return 0;
  } else {
    // Need to copy pixels to fill the image in horizontal direction.
    if (!read_and_fill_frame_data(file, input_width,    input_height,    array_width,    img_out->y)) return 0;
    if (!read_and_fill_frame_data(file, uv_input_width, uv_input_height, uv_array_width, img_out->u)) return 0;
    if (!read_and_fill_frame_data(file, uv_input_width, uv_input_height, uv_array_width, img_out->v)) return 0;
  }

  if (input_height != array_height) {
    // Need to copy pixels to fill the image in vertical direction.
    fill_after_frame(input_height,    array_width,    array_height,    img_out->y);
    fill_after_frame(uv_input_height, uv_array_width, uv_array_height, img_out->u);
    fill_after_frame(uv_input_height, uv_array_width, uv_array_height, img_out->v);
  }

  return 1;
}


/**
 * \brief Seek forward in a YUV input file.
 *
 * \param file          the input file
 * \param frames        number of frames to seek
 * \param input_width   width of the input video in pixels
 * \param input_height  height of the input video in pixels
 *
 * \return              1 on success, 0 on failure
 */
int yuv_input_seek(FILE* file, unsigned frames,
                   unsigned input_width, unsigned input_height)
{
    const size_t frame_bytes = input_width * input_height * 3 / 2;
    const size_t skip_bytes = frames * frame_bytes;

    // Attempt to seek normally.
    int error = fseek(file, skip_bytes, SEEK_CUR);
    if (!error) return 1;

    // Seek failed. Skip data by reading.
    unsigned char* tmp[4096];
    size_t bytes_left = skip_bytes;
    while (bytes_left > 0 && !error) {
      const size_t skip = MIN(4096, bytes_left);
      error = fread(tmp, sizeof(unsigned char), skip, file) != skip;
      bytes_left -= skip;
    }

    return !error || feof(file);
}
