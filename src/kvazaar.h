#ifndef KVAZAAR_H_
#define KVAZAAR_H_
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

/**
 * \file
 * \brief This file defines the public API of Kvazaar when used as a library.
 */

#include <stddef.h>
#include <stdint.h>

#include "kvazaar_version.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct config_t kvz_cfg;
typedef struct encoder_state_t encoder_state_t;
typedef struct kvz_payload kvz_payload;
typedef struct encoder_control_t encoder_control_t;

/**
* A payload unit containing at most a single frame.
* If next is not NULL, the bytestream continues in that payload unit.
*/
typedef struct kvz_payload {
  uint32_t     type;
  uint32_t     size_bytes;
  uint8_t     *payload;
  kvz_payload *next;
} kvz_payload;

typedef struct kvz_picture {
  void*    planes[3];
  uint32_t stride[3];
  uint8_t  bit_depth;
  uint8_t  poc;
} kvz_picture;

/**
 * Main datastructure representing one instance of the encoder.
 * - encoder_open
 * - encoder_close
 */
typedef struct kvz_encoder {
  encoder_control_t* control;
  encoder_state_t* states;
  int num_encoder_states;
  int cur_state_num;

  size_t bitstream_length;
} kvz_encoder;


typedef struct kvz_api {
  kvz_cfg *     (*config_alloc)(void);
  int           (*config_destroy)(kvz_cfg *);
  int           (*config_init)(kvz_cfg *);
  int           (*config_read)(kvz_cfg *, int argc, char *argv[]);

  kvz_picture * (*picture_create)(void);
  void          (*picture_destroy)(kvz_picture *);

  kvz_encoder * (*encoder_open)(kvz_cfg *);
  void          (*encoder_close)(kvz_encoder *);

  // \brief Encode one picture.
  // \param encoder  
  // \param pic_in    Picture containing the encoded data.
  // \param pic_out   Picture containing the reconstructed data.
  // \param nals_out  The first NAL containing bitstream generated, or NULL.
  // \return 1 on success, negative on error.
  int           (*encoder_encode)(kvz_encoder *encoder, kvz_picture *pic_in, kvz_picture *pic_out, kvz_payload **payload);
} kvz_api;

// Append API version to the getters name to prevent linking against incompatible versions.
#define KVZ_API_CONCAT(func, version) func ## _apiv ## version
#define KVZ_API_EXPAND_VERSION(func, version) KVZ_API_CONCAT(func, version)
#define kvz_api_get KVZ_API_EXPAND_VERSION(kvz_api_get, KVZ_API_VERSION)
const kvz_api* kvz_api_get(int bit_depth);

#ifdef __cplusplus
}
#endif

#endif // KVAZAAR_H_
