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
* \brief Picture buffer type for yuv frames.
*/
typedef struct{
  pic_buffer_t* y;
  pic_buffer_t* u;
  pic_buffer_t* v;
} yuv_buffer_t;

/**
 * \brief Create a Picture buffer. The caller is responsible for deallocation
 */
pic_buffer_t* newPictureBuffer(int width, int height, int has_tmp_row);

/**
* \brief Create/Initialize a Picture buffer. Widht/height should be the width/height of the data. The caller is responsible for deallocation.
*/
pic_buffer_t* newPictureBuffer_double(double* data, int width, int height, int has_tmp_row);

/**
* \brief Create/Initialize a yuv buffer. Widht/height should be the width/height of the data. The caller is responsible for deallocation
*/
yuv_buffer_t* newYuvBuffer_double(double* y_data, double* u_data, double* v_data, int width, int height, int is_420);

/**
* \brief Create/Initialize a Picture buffer. The caller is responsible for deallocation
*/
void deallocateYuvBuffer(yuv_buffer_t* yuv);

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

//TODO: Move to .c?
//TODO: Add offsets/cropping
typedef struct{
  //Original parameters
  int src_width;
  int src_height;

  int trgt_width;
  int trgt_height;

  //Resampling parameters
  int rnd_trgt_width;
  int rnd_trgt_height;

  int rnd_src_width;
  int rnd_src_height;

  //Sample positional parameters
  int right_offset;
  int bottom_offset;

  int shift_x;
  int shift_y;

  int scale_x;
  int scale_y;

  int add_x;
  int add_y;

} scaling_parameter_t;


/**
* \brief Function for getting initial scaling parameters given src and trgt size parameters.
*/
scaling_parameter_t newScalingParameters(int src_width, int src_height, int trgt_width, int trgt_height);


//TODO: Return/recycle the same buffer for the scaled yuv
/**
* \brief Function for scaling a yuv picture.
*/
yuv_buffer_t* yuvDownscaling(const yuv_buffer_t* const yuv, const scaling_parameter_t* const base_param, int is_420);


#endif