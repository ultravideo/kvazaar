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

#if defined(KVZ_STATIC_LIB)
  // Using or building kvazaar as a static library.
  #define KVZ_PUBLIC
#elif defined(_WIN32) || defined(__CYGWIN__)
  #ifdef KVZ_DLL_EXPORTS
    // Building kvazaar on windows.
    #define KVZ_PUBLIC __declspec(dllexport)
  #else
    // Using kvazaar as a DLL on windows.
    #define KVZ_PUBLIC __declspec(dllimport)
  #endif
#elif defined(__GNUC__)
  // Using GCC and not on windows.
  #define KVZ_PUBLIC __attribute__ ((visibility ("default")))
#else
  #define KVZ_PUBLIC
#endif

/**
 * Maximum length of a GoP structure.
 */
#define KVZ_MAX_GOP_LENGTH 32

/**
 * Size of data chunks.
 */
#define KVZ_DATA_CHUNK_SIZE 4096

#define KVZ_BIT_DEPTH 8
#if KVZ_BIT_DEPTH == 8
typedef uint8_t kvz_pixel;
#else
typedef uint16_t kvz_pixel;
#endif

/**
 * \brief Opaque data structure representing one instance of the encoder.
 */
typedef struct kvz_encoder kvz_encoder;

/**
 * \brief Integer motion estimation algorithms.
 */
enum kvz_ime_algorithm {
  KVZ_IME_HEXBS = 0,
  KVZ_IME_TZ = 1,
};

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
  int32_t qp;        /*!< \brief Quantization parameter */
  int32_t intra_period; /*!< \brief the period of intra frames in stream */
  int32_t vps_period; /*!< \brief how often the vps is re-sent */
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
  enum kvz_ime_algorithm ime_algorithm;  /*!< \brief Integer motion estimation algorithm. */
  int32_t fme_level;      /*!< \brief Fractional pixel motion estimation level (0: disabled, 1: enabled). */
  int8_t source_scan_type; /*!< \brief Source scan type (0: progressive, 1: top field first, 2: bottom field first).*/
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

/**
* \brief Struct which contains all picture data
*/
typedef struct kvz_picture {
  kvz_pixel *fulldata;         //!< \brief Allocated buffer (only used in the base_image)

  kvz_pixel *y;                //!< \brief Pointer to luma pixel array.
  kvz_pixel *u;                //!< \brief Pointer to chroma U pixel array.
  kvz_pixel *v;                //!< \brief Pointer to chroma V pixel array.
  kvz_pixel *data[3]; //!< \brief Alternate access method to same data.

  int32_t width;           //!< \brief Luma pixel array width.
  int32_t height;          //!< \brief Luma pixel array height.

  int32_t stride;          //!< \brief Luma pixel array width for the full picture (should be used as stride)

  struct kvz_picture *base_image; //!< \brief Pointer to the picture which owns the pixels
  int32_t refcount;        //!< \brief Number of references to the picture
} kvz_picture;

/**
 * \brief A linked list of chunks of data.
 *
 * Used for returning the encoded data.
 */
typedef struct kvz_data_chunk {
  /// \brief Buffer for the data.
  uint8_t data[KVZ_DATA_CHUNK_SIZE];

  /// \brief Number of bytes filled in this chunk.
  uint32_t len;

  /// \brief Next chunk in the list.
  struct kvz_data_chunk *next;
} kvz_data_chunk;

typedef struct kvz_api {
  kvz_config *  (*config_alloc)(void);
  int           (*config_destroy)(kvz_config *);
  int           (*config_init)(kvz_config *);
  int           (*config_parse)(kvz_config *, const char *name, const char *value);

  kvz_picture * (*picture_alloc)(int32_t width, int32_t height);
  void          (*picture_free)(kvz_picture *pic);

  /**
   * \brief Free a list of data chunks.
   */
  void          (*chunk_free)(kvz_data_chunk *chunk);

  kvz_encoder * (*encoder_open)(const kvz_config *);
  void          (*encoder_close)(kvz_encoder *);

  /**
   * \brief Encode one picture.
   *
   * The caller must not modify pic_in after passing it to this function.
   *
   * If pic_out and data_out are set to non-NULL values, the caller is
   * responsible for calling picture_free and chunk_free on them.
   *
   * \param encoder   Encoder
   * \param pic_in    Input frame
   * \param data_out  Returns the encoded data.
   * \param len_out   Returns number of bytes in the encoded data.
   * \param pic_out   Returns the reconstructed picture.
   * \return 1 on success, 0 on error.
   */
  int           (*encoder_encode)(kvz_encoder *encoder,
                                  kvz_picture *pic_in,
                                  kvz_data_chunk **data_out,
                                  uint32_t *len_out,
                                  kvz_picture **pic_out);
} kvz_api;

// Append API version to the getters name to prevent linking against incompatible versions.
#define KVZ_API_CONCAT(func, version) func ## _apiv ## version
#define KVZ_API_EXPAND_VERSION(func, version) KVZ_API_CONCAT(func, version)
#define kvz_api_get KVZ_API_EXPAND_VERSION(kvz_api_get, KVZ_API_VERSION)

KVZ_PUBLIC const kvz_api * kvz_api_get(int bit_depth);

#ifdef __cplusplus
}
#endif

#endif // KVAZAAR_H_
