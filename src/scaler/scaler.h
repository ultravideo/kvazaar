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

typedef pic_data_t int; //Use some other type?

/**
 * \brief Picture buffer type for operating on image data.
 */
typedef struct{
  pic_data_t* data;
  int width;
  int height;
} pic_buffer_t;

/**
* \brief Create/Initialize a Picture buffer. The caller is responsible for deallocation
*/
pic_buffer_t newPictureBuffer(int width, int height);

/**
* \brief Deallocate a picture buffer.
*/
void deallocatePictureBuffer(pic_buffer_t buffer);

/**
* \brief Struct for passing scaling parameters.
*/
typedef struct{
  int src_width;
  int src_height;

  int trgt_width;
  int trgt_height;
}

#endif