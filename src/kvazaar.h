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

/**
 * Maximum length of a GoP structure.
 */
#define KVZ_MAX_GOP_LENGTH 32

#define KVZ_BIT_DEPTH 8
#if KVZ_BIT_DEPTH == 8
typedef uint8_t pixel_t;
#else
typedef uint16_t pixel_t;
#endif

/**
 * \brief GoP picture configuration.
 */
typedef struct kvz_gop_config {
  double qp_factor;
  int8_t qp_offset;    /*!< \brief QP offset */
  int8_t poc_offset;   /*!< \brief POC offset */
  int8_t layer;        /*!< \brief Current layer */
  int8_t is_ref;       /*!< \brief Flag if this picture is used as a reference */
  int8_t ref_pos_count;/*!< \brief Reference picture count */
  int8_t ref_pos[16];  /*!< \brief reference picture offset list */
  int8_t ref_neg_count;/*!< \brief Reference picture count */
  int8_t ref_neg[16];  /*!< \brief reference picture offset list */
} kvz_gop_config;

/**
 * \brief Struct which contains all configuration data
 */
typedef struct kvz_config
{
  char *input;      /*!< \brief Pointer to input filename  */
  char *output;     /*!< \brief Pointer to output filename */
  char *debug;      /*!< \brief Pointer to debug output    */
  int32_t qp;        /*!< \brief Quantization parameter */
  int32_t intra_period; /*!< \brief the period of intra frames in stream */
  int32_t vps_period; /*!< \brief how often the vps is re-sent */
  int32_t frames;  /*!< \brief Number of frames to decode */
  int32_t width;   /*!< \brief frame width */
  int32_t height;  /*!< \brief frame height */
  double framerate; /*!< \brief Input framerate */
  int32_t deblock_enable; /*!< \brief Flag to enable deblocking filter */
  int32_t sao_enable;     /*!< \brief Flag to enable sample adaptive offset filter */
  int32_t rdoq_enable;    /*!< \brief Flag to enable RD optimized quantization. */
  int32_t signhide_enable;   /*!< \brief Flag to enable sign hiding. */
  int32_t rdo;            /*!< \brief RD-calculation level (0..2) */
  int32_t full_intra_search; /*!< \brief If true, don't skip modes in intra search. */
  int32_t trskip_enable;    /*!< \brief Flag to enable transform skip (for 4x4 blocks). */
  int32_t tr_depth_intra; /*!< \brief Maximum transform depth for intra. */
  int8_t  ime_algorithm;  /*!< \brief Integer motion estimation algorithm. */
  int32_t fme_level;      /*!< \brief Fractional pixel motion estimation level (0: disabled, 1: enabled). */
  int32_t bipred;         /*!< \brief Bi-prediction (0: disabled, 1: enabled). */
  int32_t deblock_beta;   /*!< \brief (deblocking) beta offset (div 2), range -6...6 */
  int32_t deblock_tc;     /*!< \brief (deblocking) tc offset (div 2), range -6...6 */
  struct
  {
    int32_t sar_width;   /*!< \brief the horizontal size of the sample aspect ratio (in arbitrary units) */
    int32_t sar_height;  /*!< \brief the vertical size of the sample aspect ratio (in the same arbitrary units as sar_width). */
    int8_t overscan;     /*!< \brief Crop overscan setting */
    int8_t videoformat;  /*!< \brief Video format */
    int8_t fullrange;    /*!< \brief Flag to indicate full-range */
    int8_t colorprim;    /*!< \brief Color primaries */
    int8_t transfer;     /*!< \brief Transfer characteristics */
    int8_t colormatrix;  /*!< \brief Color matrix coefficients */
    int32_t chroma_loc;   /*!< \brief Chroma sample location */
  } vui;
  int32_t aud_enable;     /*!< \brief Flag to use access unit delimiters */
  int32_t ref_frames;     /*!< \brief number of reference frames to use */
  char * cqmfile;        /*!< \brief Pointer to custom quantization matrices filename */
  int32_t seek;           /*!< \brief Number of frames to skip in the beginning of input. */
  
  int32_t tiles_width_count;      /*!< \brief number of tiles separation in x direction */
  int32_t tiles_height_count;      /*!< \brief number of tiles separation in y direction */
  int32_t* tiles_width_split;      /*!< \brief tiles split x coordinates (dimension: tiles_width_count) */
  int32_t* tiles_height_split;      /*!< \brief tiles split y coordinates (dimension: tiles_height_count) */
  
  int wpp;
  int owf;
  
  int32_t slice_count;
  int32_t* slice_addresses_in_ts;
  
  int32_t threads;
  int32_t cpuid;

  struct {
    int32_t min;
    int32_t max;
  } pu_depth_inter, pu_depth_intra;

  int32_t add_encoder_info;
  int8_t gop_len;            /*!< \brief length of GOP for the video sequence */
  kvz_gop_config gop[KVZ_MAX_GOP_LENGTH];  /*!< \brief Array of GOP settings */

  int32_t target_bitrate;
} kvz_config;

typedef struct encoder_state_t encoder_state_t;
typedef struct encoder_control_t encoder_control_t;
typedef struct bitstream_chunk_t kvz_payload;

/**
* \brief Struct which contains all picture data
*/
typedef struct kvz_picture {
  pixel_t *fulldata;         //!< \brief Allocated buffer (only used in the base_image)

  pixel_t *y;                //!< \brief Pointer to luma pixel array.
  pixel_t *u;                //!< \brief Pointer to chroma U pixel array.
  pixel_t *v;                //!< \brief Pointer to chroma V pixel array.
  pixel_t *data[3]; //!< \brief Alternate access method to same data.

  int32_t width;           //!< \brief Luma pixel array width.
  int32_t height;          //!< \brief Luma pixel array height.

  int32_t stride;          //!< \brief Luma pixel array width for the full picture (should be used as stride)

  struct kvz_picture *base_image; //!< \brief Pointer to the picture which owns the pixels
  int32_t refcount;        //!< \brief Number of references to the picture
} kvz_picture;

/**
 * Main datastructure representing one instance of the encoder.
 * - encoder_open
 * - encoder_close
 */
typedef struct kvz_encoder {
  encoder_control_t* control;
  encoder_state_t* states;
  unsigned num_encoder_states;
  unsigned cur_state_num;
  unsigned frames_started;
  unsigned frames_done;

  size_t bitstream_length;
} kvz_encoder;

typedef struct kvz_api {
  kvz_config *  (*config_alloc)(void);
  int           (*config_destroy)(kvz_config *);
  int           (*config_init)(kvz_config *);
  int           (*config_parse)(kvz_config *, const char *name, const char *value);

  kvz_encoder * (*encoder_open)(kvz_config *);
  void          (*encoder_close)(kvz_encoder *);

  // \brief Encode one picture.
  // \param encoder   Encoder
  // \param pic_in    Input frame
  // \param pic_out   Returns the reconstructed picture.
  // \param payload   Returns the encoded data.
  // \return 1 on success, 0 on error.
  int           (*encoder_encode)(kvz_encoder *encoder, kvz_picture *pic_in, kvz_picture **pic_out, kvz_payload **payload);
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
