#ifndef SCALER_H_
#define SCALER_H_
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

typedef int pic_data_t; //Use some other type?

/**
 * \brief Picture buffer type for operating on image data.
 */
typedef struct{
  pic_data_t* data; //Contain main data
  pic_data_t* tmp_row; //A temporary buffer row that may be used to hold data when operating on buffer

  int width;
  int height;

} pic_buffer_t;

/**
 * \brief Create/Initialize a Picture buffer. The caller is responsible for deallocation
 */
pic_buffer_t* newPictureBuffer(int width, int height, int has_tmp_row);

/**
 * \brief Deallocate a picture buffer.
 */
void deallocatePictureBuffer(pic_buffer_t* buffer);

/**
 * \brief Copies data from one buffer to the other.
 * \param src is the source buffer
 * \param dst is the destination buffer
 * \param fill signals if the inds in dst not overlapped by src should be filled
*    with values adjacent to the said index.
 */
void copyPictureBuffer(pic_buffer_t* src, pic_buffer_t* dst, int fill);

typedef struct{
  //Original parameters
  int src_width;
  int src_height;

  int trgt_width;
  int trgt_height;

  //Resampling parameters
  int rnd_trgt_width;
  int rnd_trgt_height;

  int right_offset;
  int bottom_offset;

  int shift_x;
  int shift_y;

  int scale_x;
  int scale_y;

  int add_x;
  int add_y;

} scaling_parameter_t;

#endif