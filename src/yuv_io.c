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
                                    unsigned width, unsigned height,
                                    unsigned array_width, kvz_pixel *data)
{
  kvz_pixel* p = data;
  kvz_pixel* end = data + array_width * height;
  kvz_pixel fill_char;
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
* \brief Convert 8 bit (single byte per pixel) to 10bit (two bytes per pixel) array
*
* \param input   input/output buffer
* \return        1
*/
int frame_8bit_to_10bit(kvz_pixel* input, int width, int height) {
  uint8_t* temp_buffer = (uint8_t*)input;
  const uint32_t pixels = width*height;
  for (int i = pixels - 1; i >= 0; i--) {
    input[i] = temp_buffer[i] << 2;
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
                unsigned input_width, unsigned input_height,
                kvz_picture *img_out)
{
  assert(input_width % 2 == 0);
  assert(input_height % 2 == 0);

  const unsigned y_size = input_width * input_height;
  const unsigned uv_input_width  = input_width  / 2;
  const unsigned uv_input_height = input_height / 2;
  const unsigned uv_size = uv_input_width * uv_input_height;

  const unsigned uv_array_width  = img_out->width  / 2;
  const unsigned uv_array_height = img_out->height  / 2;

  if (input_width == img_out->width) {
    // No need to extend pixels.
    const size_t pixel_size = sizeof(unsigned char);
    if (fread(img_out->y, pixel_size, y_size,  file) != y_size)  return 0;
    if (fread(img_out->u, pixel_size, uv_size, file) != uv_size) return 0;
    if (fread(img_out->v, pixel_size, uv_size, file) != uv_size) return 0;
  } else {
    // Need to copy pixels to fill the image in horizontal direction.
    if (!read_and_fill_frame_data(file, input_width,    input_height,    img_out->width, img_out->y)) return 0;
    if (!read_and_fill_frame_data(file, uv_input_width, uv_input_height, uv_array_width, img_out->u)) return 0;
    if (!read_and_fill_frame_data(file, uv_input_width, uv_input_height, uv_array_width, img_out->v)) return 0;
  }

  if (input_height != img_out->height) {
    // Need to copy pixels to fill the image in vertical direction.
    fill_after_frame(input_height,    img_out->width, img_out->height,    img_out->y);
    fill_after_frame(uv_input_height, uv_array_width, uv_array_height, img_out->u);
    fill_after_frame(uv_input_height, uv_array_width, uv_array_height, img_out->v);
  }

#if KVZ_BIT_DEPTH == 10
  frame_8bit_to_10bit(img_out->y, img_out->width, img_out->height);
	frame_8bit_to_10bit(img_out->u, img_out->width >> 1, img_out->height >> 1);
	frame_8bit_to_10bit(img_out->v, img_out->width >> 1, img_out->height >> 1);
#endif

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
  for (int y = 0; y < output_height / 2; ++y) {
    fwrite(&img->u[y * width / 2], sizeof(*img->u), output_width / 2, file);
  }
  for (int y = 0; y < output_height / 2; ++y) {
    fwrite(&img->v[y * width / 2], sizeof(*img->v), output_width / 2, file);
  }

  return 1;
}

/**
* \brief Separate a single field from a frame.
*
* \param frame_in           input frame to extract field from
* \param source_scan_type   scan type of input material (0: progressive, 1:top field first, 2:bottom field first)
* \param field parity   
* \param field_out
*
* \return              1 on success, 0 on failure
*/
int yuv_io_extract_field(const kvz_picture *frame_in, unsigned source_scan_type, unsigned field_parity, kvz_picture *field_out)
{
  if ((source_scan_type != 1) && (source_scan_type != 2)) return 0;
  if ((field_parity != 0)     && (field_parity != 1))     return 0;

  unsigned offset = 0;
  if (source_scan_type == 1) offset = field_parity ? 1 : 0;
  else if (source_scan_type == 2) offset = field_parity ? 0 : 1;  

  //Luma
  for (int i = 0; i < field_out->height; ++i){
    kvz_pixel *row_in  = frame_in->y + MIN(frame_in->height - 1, 2 * i + offset) * frame_in->stride;
    kvz_pixel *row_out = field_out->y + i * field_out->stride;
    memcpy(row_out, row_in, sizeof(kvz_pixel) * frame_in->width);
  }

  //Chroma
  for (int i = 0; i < field_out->height / 2; ++i){
    kvz_pixel *row_in = frame_in->u + MIN(frame_in->height / 2 - 1, 2 * i + offset) * frame_in->stride / 2;
    kvz_pixel *row_out = field_out->u + i * field_out->stride / 2;
    memcpy(row_out, row_in, sizeof(kvz_pixel) * frame_in->width / 2);
  }

  for (int i = 0; i < field_out->height / 2; ++i){
    kvz_pixel *row_in = frame_in->v + MIN(frame_in->height / 2 - 1, 2 * i + offset) * frame_in->stride / 2;
    kvz_pixel *row_out = field_out->v + i * field_out->stride / 2;
    memcpy(row_out, row_in, sizeof(kvz_pixel) * frame_in->width / 2);
  }

  return 1;
}