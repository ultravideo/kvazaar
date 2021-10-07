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
  int shift = to_bitdepth - from_bitdepth;
  kvz_pixel bitdepth_mask = (1 << from_bitdepth) - 1;

  for (int i = 0; i < size; ++i) {
    // Shifting by a negative number is undefined.
    if (shift > 0) {
      input[i] = (input[i] & bitdepth_mask) << shift;
    } else {
      input[i] = (input[i] & bitdepth_mask) >> shift;
    }
  }
}


// Shift and copy 1-byte aligned samples to 2-byte aligned array
static void shift_to_bitdepth_and_spread(kvz_pixel *input,
                                         int size,
                                         int from_bitdepth,
                                         int to_bitdepth)
{
  assert(sizeof(kvz_pixel) > 1);
  int shift = to_bitdepth - from_bitdepth;
  unsigned char *byte_buf = (unsigned char *)input;
  kvz_pixel bitdepth_mask = (1 << from_bitdepth) - 1;
  
  // Starting from the back of the 1-byte samples, copy each sample to it's
  // place in the 2-byte per sample array, overwriting the bytes that have
  // already been copied in the process.
  // Even though the two pointers are aliased, this should work because the
  // future values read through byte_buf poiner never change as a result of
  // writing through input pointer.
  for (int i = size - 1; i >= 0; --i) {
    // Shifting by a negative number is undefined.
    if (shift > 0) {
      input[i] = (byte_buf[i] & bitdepth_mask) << shift;
    } else {
      input[i] = (byte_buf[i] & bitdepth_mask) >> shift;
    }
  }
}


static bool machine_is_big_endian()
{
  // Big and little endianess refers to which end of the egg you prefer to eat
  // first. Therefore in big endian system, the most significant bits are in
  // the first address.

  uint16_t number = 1;
  char first_byte = *(char*)&number;

  return (first_byte == 0);
}


static void mask_to_bitdepth(kvz_pixel *buf, unsigned length, unsigned bitdepth)
{
  kvz_pixel bitdepth_mask = (1 << bitdepth) - 1;
  for (int i = 0; i < length; ++i) {
    buf[i] = buf[i] & bitdepth_mask;
  }
}


static int yuv_io_read_plane(
    FILE* file,
    unsigned in_width, unsigned in_height, unsigned in_bitdepth,
    unsigned out_width, unsigned out_height, unsigned out_bitdepth,
    kvz_pixel *out_buf)
{
  unsigned bytes_per_sample = in_bitdepth > 8 ? 2 : 1;
  unsigned buf_bytes = in_width * in_height * bytes_per_sample;
  unsigned out_length = out_width * out_height;

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
    // Assume little endian input.
    if (machine_is_big_endian()) {
      swap_16b_buffer_bytes(out_buf, out_length);
    }
  }

  // Shift the data to the correct bitdepth.
  // Ignore any bits larger than in_bitdepth to guarantee ouput data will be
  // in the correct range.
  if (in_bitdepth <= 8 && out_bitdepth > 8) {
    shift_to_bitdepth_and_spread(out_buf, out_length, in_bitdepth, out_bitdepth);
  } else if (in_bitdepth != out_bitdepth) {
    shift_to_bitdepth(out_buf, out_length, in_bitdepth, out_bitdepth);
  } else if (in_bitdepth % 8 != 0) {
    mask_to_bitdepth(out_buf, out_length, out_bitdepth);
  }

  return 1;
}


static int read_frame_header(FILE* input) {
  char buffer[256];
  bool frame_start = false;

  while (!frame_start) {
    for (int i = 0; i < 256; i++) {
      buffer[i] = getc(input);
      if (buffer[i] == EOF) return 0;
      // ToDo: frame headers can have some information structured same as start headers
      // This info is just skipped for now, since it's not clear what it could be.
      if (buffer[i] == 0x0A) {
        frame_start = true;
        break;
      }
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
                kvz_picture *img_out, unsigned file_format)
{
  assert(in_width % 2 == 0);
  assert(out_width % 2 == 0);

  int ok;

  if (file_format == KVZ_FORMAT_Y4M) {
    ok = read_frame_header(file);
    if (!ok) return 0;
  }

  

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
                unsigned input_width, unsigned input_height,
                unsigned file_format)
{
    const size_t frame_bytes = input_width * input_height * 3 / 2;

    if (file_format == KVZ_FORMAT_Y4M) {
      for (unsigned i = 0; i < frames; i++) {
        if (!read_frame_header(file)) return 0;
        if (fseek(file, frame_bytes, SEEK_CUR)) return 0;
      }
      return 1;
    }

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
