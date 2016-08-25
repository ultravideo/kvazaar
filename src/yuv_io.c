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

#include "yuv_io.h"

static void fill_after_frame(unsigned height, unsigned array_width,
                             unsigned array_height, kvz_pixel *data)
{
  kvz_pixel* p = data + height * array_width;
  kvz_pixel* end = data + array_width * array_height;

  while (p < end) {
    // Fill the line by copying the line above.
    memcpy(p, p - array_width, array_width);
    p += array_width;
  }
}


static int read_and_fill_frame_data(FILE *file,
                                    unsigned width, unsigned height, unsigned bytes_per_sample,
                                    unsigned array_width, kvz_pixel *data)
{
  kvz_pixel* p = data;
  kvz_pixel* end = data + array_width * height;
  kvz_pixel fill_char;
  unsigned i;

  while (p < end) {
    // Read the beginning of the line from input.
    if (width != fread(p, bytes_per_sample, width, file))
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


static void swap_16b_buffer_bytes(kvz_pixel* input, int size)
{
  for (int i = 0; i < size; ++i) {
    input[i] = ((input[i] & 0xff) << 8) + ((input[i] & 0xff00) >> 8);
  }
}


static void shift_to_bitdepth(kvz_pixel* input, int size, int from_bitdepth, int to_bitdepth)
{
  int shift = from_bitdepth - to_bitdepth;
  for (int i = 0; i < size; ++i) {
    // Shifting by a negative number is undefined.
    if (shift > 0) {
      input[i] <<= shift;
    } else {
      input[i] >>= shift;
    }
  }
}


bool machine_is_big_endian()
{
  uint16_t number = 1;
  char first_byte = *(char*)&number;

  return (first_byte != 0);
}


static int yuv_io_read_plane(
    FILE* file,
    unsigned in_width, unsigned in_height, unsigned in_bitdepth,
    unsigned out_width, unsigned out_height, unsigned out_bitdepth,
    kvz_pixel *out_buf)
{
  unsigned bytes_per_sample = in_bitdepth > 8 ? 2 : 1;
  unsigned buf_length = in_width * in_height;
  unsigned buf_bytes = buf_length * bytes_per_sample;

  if (in_width == out_width) {
    // No need to extend pixels.
    const size_t pixel_size = sizeof(unsigned char);
    if (fread(out_buf, pixel_size, buf_bytes, file) != buf_bytes)  return 0;
  } else {
    // Need to copy pixels to fill the image in horizontal direction.
    if (!read_and_fill_frame_data(file, in_width, in_height, bytes_per_sample, out_width, out_buf)) return 0;
  }

  if (in_height != out_height) {
    // Need to copy pixels to fill the image in vertical direction.
    fill_after_frame(in_height, out_width, out_height, out_buf);
  }

  if (in_bitdepth > 8) {
    if (machine_is_big_endian()) {
      swap_16b_buffer_bytes(out_buf, buf_length);
    }

    if (in_bitdepth != out_bitdepth) {
      shift_to_bitdepth(out_buf, buf_length, in_bitdepth, out_bitdepth);
    }
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
 * \param img_out       image buffer
 *
 * \return              1 on success, 0 on failure
 */
int yuv_io_read(FILE* file,
                unsigned in_width, unsigned out_width,
                unsigned in_bitdepth, unsigned out_bitdepth,
                kvz_picture *img_out)
{
  assert(in_width % 2 == 0);
  assert(out_width % 2 == 0);

  int ok;

  ok = yuv_io_read_plane(
      file, 
      in_width, out_width, in_bitdepth,
      img_out->width, img_out->height, out_bitdepth,
      img_out->y);
  if (!ok) return 0;

  if (img_out->chroma_format != KVZ_CSP_400) {
    unsigned uv_width_in = in_width / 2;
    unsigned uv_height_in = out_width / 2;
    unsigned uv_width_out = img_out->width / 2;
    unsigned uv_height_out = img_out->height / 2;

    ok = yuv_io_read_plane(
        file,
        uv_width_in, uv_height_in, in_bitdepth,
        uv_width_out, uv_height_out, out_bitdepth,
        img_out->u);
    if (!ok) return 0;

    ok = yuv_io_read_plane(
        file, 
        uv_width_in, uv_height_in, in_bitdepth,
        uv_width_out, uv_height_out, out_bitdepth,
        img_out->v);
    if (!ok) return 0;
  }
  
  return 1;
}


/**
 * \brief Seek forward in a YUV file.
 *
 * \param file          the input file
 * \param frames        number of frames to seek
 * \param input_width   width of the input video in pixels
 * \param input_height  height of the input video in pixels
 *
 * \return              1 on success, 0 on failure
 */
int yuv_io_seek(FILE* file, unsigned frames,
                unsigned input_width, unsigned input_height)
{
    const size_t frame_bytes = input_width * input_height * 3 / 2;
    const int64_t skip_bytes = (int64_t)(frames * frame_bytes);

    // Attempt to seek normally.
    size_t error = fseek(file, skip_bytes, SEEK_CUR);
    if (!error) return 1;

    // Seek failed. Skip data by reading.
    error = 0;
    unsigned char* tmp[4096];
    size_t bytes_left = skip_bytes;
    while (bytes_left > 0 && !error) {
      const size_t skip = MIN(4096, bytes_left);
      error = fread(tmp, sizeof(unsigned char), skip, file) != skip;
      bytes_left -= skip;
    }

    return !error || feof(file);
}


/**
 * \brief Write a single frame to a file.
 *
 * \param file           output file
 * \param img            image to output
 * \param output_width   width of the output in pixels
 * \param output_height  height of the output in pixels
 *
 * \return              1 on success, 0 on failure
 */
int yuv_io_write(FILE* file,
                const kvz_picture *img,
                unsigned output_width, unsigned output_height)
{
  const int width = img->width;
  for (int y = 0; y < output_height; ++y) {
    fwrite(&img->y[y * width], sizeof(*img->y), output_width, file);
    // TODO: Check that fwrite succeeded.
  }

  if (img->chroma_format != KVZ_CSP_400) {
    for (int y = 0; y < output_height / 2; ++y) {
      fwrite(&img->u[y * width / 2], sizeof(*img->u), output_width / 2, file);
    }
    for (int y = 0; y < output_height / 2; ++y) {
      fwrite(&img->v[y * width / 2], sizeof(*img->v), output_width / 2, file);
    }
  }

  return 1;
}
