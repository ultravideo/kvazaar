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

/*
 * The CLI is compiled once for all the bit depths that the library supports,
 * so the samples of the pictures are accessed with the sample size of the
 * bit depth in use (1 byte for 8-bit, 2 bytes otherwise) instead of
 * kvz_pixel.
 */

static unsigned bytes_per_sample(unsigned bitdepth)
{
  return bitdepth > 8 ? 2 : 1;
}


static void fill_after_frame(unsigned height, unsigned array_width,
                             unsigned array_height, unsigned sample_bytes,
                             void *data)
{
  const size_t line_bytes = (size_t)array_width * sample_bytes;
  uint8_t* p = (uint8_t*)data + height * line_bytes;
  uint8_t* end = (uint8_t*)data + array_height * line_bytes;

  while (p < end) {
    // Fill the line by copying the line above.
    memcpy(p, p - line_bytes, line_bytes);
    p += line_bytes;
  }
}


static int read_and_fill_frame_data(FILE *file,
                                    unsigned width, unsigned height, unsigned in_sample_bytes,
                                    unsigned array_width, unsigned out_sample_bytes,
                                    void *data)
{
  unsigned i;
  // If the input has 1 byte per sample and the picture 2 (8-bit input coded
  // with a higher bit depth), the lines are read with 1 byte per sample and
  // spread out to 2 bytes per sample later on.
  const unsigned sample_bytes = MIN(in_sample_bytes, out_sample_bytes);

  if (sample_bytes == 1) {
    uint8_t* p = (uint8_t*)data;
    uint8_t* end = (uint8_t*)data + array_width * height;
    while (p < end) {
      // Read the beginning of the line from input.
      if (width != fread(p, 1, width, file)) return 0;
      // Fill the rest with the last pixel value.
      const uint8_t fill_char = p[width - 1];
      for (i = width; i < array_width; ++i) {
        p[i] = fill_char;
      }
      p += array_width;
    }
  } else {
    uint16_t* p = (uint16_t*)data;
    uint16_t* end = (uint16_t*)data + array_width * height;
    while (p < end) {
      // Read the beginning of the line from input.
      if (width != fread(p, 2, width, file)) return 0;
      // Fill the rest with the last pixel value.
      const uint16_t fill_char = p[width - 1];
      for (i = width; i < array_width; ++i) {
        p[i] = fill_char;
      }
      p += array_width;
    }
  }
  return 1;
}


static void swap_16b_buffer_bytes(uint16_t* input, int size)
{
  for (int i = 0; i < size; ++i) {
    input[i] = ((input[i] & 0xff) << 8) + ((input[i] & 0xff00) >> 8);
  }
}


static void shift_to_bitdepth(uint16_t* input, int size, int from_bitdepth, int to_bitdepth)
{
  int shift = to_bitdepth - from_bitdepth;
  uint16_t bitdepth_mask = (1 << from_bitdepth) - 1;

  for (int i = 0; i < size; ++i) {
    // Shifting by a negative number is undefined.
    if (shift > 0) {
      input[i] = (input[i] & bitdepth_mask) << shift;
    } else {
      input[i] = (input[i] & bitdepth_mask) >> -shift;
    }
  }
}


// Shift and copy 1-byte aligned samples to 2-byte aligned array
static void shift_to_bitdepth_and_spread(uint16_t *input,
                                         int size,
                                         int from_bitdepth,
                                         int to_bitdepth)
{
  int shift = to_bitdepth - from_bitdepth;
  unsigned char *byte_buf = (unsigned char *)input;
  uint16_t bitdepth_mask = (1 << from_bitdepth) - 1;
  
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
      input[i] = (byte_buf[i] & bitdepth_mask) >> -shift;
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


static void mask_to_bitdepth(uint16_t *buf, unsigned length, unsigned bitdepth)
{
  uint16_t bitdepth_mask = (1 << bitdepth) - 1;
  for (int i = 0; i < length; ++i) {
    buf[i] = buf[i] & bitdepth_mask;
  }
}


static int yuv_io_read_plane(
    FILE* file,
    unsigned in_width, unsigned in_height, unsigned in_bitdepth,
    unsigned out_width, unsigned out_height, unsigned out_bitdepth,
    void *out_buf)
{
  const unsigned in_sample_bytes = bytes_per_sample(in_bitdepth);
  const unsigned out_sample_bytes = bytes_per_sample(out_bitdepth);
  unsigned buf_bytes = in_width * in_height * in_sample_bytes;
  unsigned out_length = out_width * out_height;

  // The library refuses input with more bits than the output bit depth if
  // the output is 8-bit, so there is always room for the input samples.
  assert(in_sample_bytes <= out_sample_bytes);

  if (in_width == out_width) {
    // No need to extend pixels.
    const size_t pixel_size = sizeof(unsigned char);
    if (fread(out_buf, pixel_size, buf_bytes, file) != buf_bytes)  return 0;
  } else {
    // Need to copy pixels to fill the image in horizontal direction.
    if (!read_and_fill_frame_data(file, in_width, in_height, in_sample_bytes,
                                  out_width, out_sample_bytes, out_buf)) return 0;
  }

  if (in_height != out_height) {
    // Need to copy pixels to fill the image in vertical direction. The
    // lines have the input sample size at this point.
    fill_after_frame(in_height, out_width, out_height, in_sample_bytes, out_buf);
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
  char header[4096];
  int count = 0;
  const int max_scan = (int)sizeof(header) - 1; // prevent infinite loops

  while (count < max_scan) {
    const int c = getc(input);
    if (c == EOF) {
      return 0; // EOF: no frame header present.
    }
    // ToDo: frame headers can have some information structured same as start headers
    // This info is just skipped for now, since it's not clear what it could be.
    if (c == 0x0A) {
      // A y4m frame header must start with the "FRAME" tag. Rejecting
      // anything else keeps a malformed file from silently misaligning
      // the frame data that follows.
      return count >= 5 && strncmp(header, "FRAME", 5) == 0;
    }
    header[count++] = (char)c;
  }

  return 0; // EOF or scan limit reached
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
                unsigned in_width, unsigned in_height,
                unsigned in_bitdepth, unsigned out_bitdepth,
                kvz_picture *img_out, unsigned file_format)
{
  assert(in_width % 2 == 0);
  assert(in_height % 2 == 0);

  int ok;

  if (file_format == KVZ_FORMAT_Y4M) {
    ok = read_frame_header(file);
    if (!ok) return 0;
  }

  

  ok = yuv_io_read_plane(
      file, 
      in_width, in_height, in_bitdepth,
      img_out->width, img_out->height, out_bitdepth,
      img_out->y);
  if (!ok) return 0;

  const uint8_t chroma_shift_w = (img_out->chroma_format == 0 || img_out->chroma_format == 3) ? 0 : 1;
  const uint8_t chroma_shift_h = img_out->chroma_format == 1 ? 1 : 0;

  if (img_out->chroma_format != KVZ_CSP_400) {
    unsigned uv_width_in = in_width >> chroma_shift_w;
    unsigned uv_height_in = in_height >> chroma_shift_h;
    unsigned uv_width_out = img_out->width >> chroma_shift_w;
    unsigned uv_height_out = img_out->height >> chroma_shift_h;

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
    // ToDo: allow other chroma subsamplings
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
    unsigned char tmp[4096];
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
 * \param bitdepth       bit depth of the image
 *
 * \return              1 on success, 0 on failure
 */
int yuv_io_write(FILE* file,
                const kvz_picture *img,
                unsigned output_width, unsigned output_height,
                unsigned bitdepth)
{
  const unsigned sample_bytes = bytes_per_sample(bitdepth);
  const int width = img->width;
  const uint8_t *y_buf = (const uint8_t *)img->y;
  for (int y = 0; y < output_height; ++y) {
    fwrite(&y_buf[(size_t)y * width * sample_bytes], sample_bytes, output_width, file);
    // TODO: Check that fwrite succeeded.
  }

  if (img->chroma_format != KVZ_CSP_400) {
    const uint8_t chroma_shift_w = img->chroma_format == 0 || img->chroma_format == 3 ? 0 : 1;
    const uint8_t chroma_shift_h = img->chroma_format == 1 ? 1 : 0;
    const size_t stride_bytes = (size_t)(width >> chroma_shift_w) * sample_bytes;
    const uint8_t *u_buf = (const uint8_t *)img->u;
    const uint8_t *v_buf = (const uint8_t *)img->v;
    for (int y = 0; y < (output_height >> chroma_shift_h); ++y) {
      fwrite(&u_buf[y * stride_bytes], sample_bytes, output_width >> chroma_shift_w, file);
    }
    for (int y = 0; y < (output_height >> chroma_shift_h); ++y) {
      fwrite(&v_buf[y * stride_bytes], sample_bytes, output_width >> chroma_shift_w, file);
    }
  }

  return 1;
}
