#ifndef CONFIG_H_
#define CONFIG_H_
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
 * \brief Handles parsing and storing of configuration of the encoder.
 */

#include "global.h"


/*!
    \brief Struct which contains all configuration data
*/
typedef struct
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
  int32_t deblock_enable; /*!< \brief Flag to enable deblocking filter */
  int32_t sao_enable;     /*!< \brief Flag to enable sample adaptive offset filter */
  int32_t rdoq_enable;    /*!< \brief Flag to enable RD optimized quantization. */
  bool signhide_enable;
  int32_t rdo;            /*!< \brief RD-calculation level (0..2) */
  bool full_intra_search; /*!< \brief Don't skip modes in intra search.e */
  int32_t trskip_enable;    /*!< \brief Flag to enable transform skip (for 4x4 blocks). */
  int32_t tr_depth_intra; /*!< \brief Maximum transform depth for intra. */
  int32_t fme_level;      /*!< \brief Fractional pixel motion estimation level (0: disabled, 1: enabled). */
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
} config;

/* Function definitions */
config *config_alloc(void);
int config_init(config *cfg);
int config_destroy(config *cfg);
int config_read(config *cfg,int argc, char *argv[]);
int config_validate(config *cfg);
int config_set_owf_auto(config *cfg);

#endif
