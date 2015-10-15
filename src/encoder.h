#ifndef ENCODER_H_
#define ENCODER_H_
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
 * \brief The highest level of the encoder.
 */

#include "global.h"

#include "image.h"
#include "bitstream.h"
#include "cabac.h"
#include "config.h"
#include "tables.h"
#include "scalinglist.h"
#include "threadqueue.h"


enum { FORMAT_400 = 0, FORMAT_420, FORMAT_422, FORMAT_444 };

/* Encoder control options, the main struct */
typedef struct encoder_control_t
{
  /* Configuration */
  const kvz_config *cfg;
  
  /* Input */
  struct {
    int32_t width;
    int32_t height;
    int32_t width_in_lcu;
    int32_t height_in_lcu;
    int32_t real_width;  /*!< \brief real input picture width */
    int32_t real_height; /*!< \brief real input picture width */
    int8_t video_format;
    int8_t bitdepth;  /*!< \brief input bit depth (8,10) */
    int64_t pixels_per_pic;
    int8_t source_scan_type;
  } in;
  
  /* TODO: add ME data */
  struct {
    void(*IME)();
    void(*FME)();
    int range;
  } me;
  
  int8_t bitdepth;
  int8_t tr_depth_intra;

  int8_t fme_level;

  /* Filtering */
  int8_t deblock_enable; // \brief Flag to enable deblocking filter
  int8_t sao_enable;     // \brief Flag to enable sample adaptive offset filter
  int8_t rdoq_enable;    // \brief Whether RDOQ is enabled or not.
  int8_t rdo;            // \brief RDO level
  int8_t full_intra_search; // \brief Whether to skip intra modes during search.
  int8_t trskip_enable;    // \brief Flag to enable transform skipping (4x4 intra)
  int8_t beta_offset_div2; // \brief (deblocking) beta offset (div 2), range -6...6
  int8_t tc_offset_div2;   // \brief (deblocking)tc offset (div 2), range -6...6

  /* VUI */
  struct
  {
    int16_t sar_width;
    int16_t sar_height;
    int8_t overscan;
    int8_t videoformat;
    int8_t fullrange;
    int8_t colorprim;
    int8_t transfer;
    int8_t colormatrix;
    int8_t chroma_loc;

    int8_t field_seq_flag;
    int8_t frame_field_info_present_flag;
  } vui;

  int8_t aud_enable;

  //scaling list
  scaling_list_t scaling_list;
  
  //spec: references to variables defined in Rec. ITU-T H.265 (04/2013)
  int8_t tiles_enable; /*!<spec: tiles_enabled */
  
  int8_t tiles_uniform_spacing_flag; /*!<spec: uniform_spacing_flag */
  
  uint8_t tiles_num_tile_columns; /*!<spec: num_tile_columns_minus1 + 1 */
  uint8_t tiles_num_tile_rows; /*!<spec: num_tile_rows_minus1 + 1*/
  
  const int32_t *tiles_col_width; /*!<spec: colWidth (6.5.1); dimension: tiles_num_tile_columns */
  const int32_t *tiles_row_height; /*!<spec: rowHeight (6.5.1); dimension: tiles_num_tile_rows */
  
  const int32_t *tiles_col_bd; /*!<spec: colBd (6.5.1); dimension: tiles_num_tile_columns + 1 */
  const int32_t *tiles_row_bd; /*!<spec: rowBd (6.5.1); dimension: tiles_num_tile_rows + 1  */
  
  //PicSizeInCtbsY = height_in_lcu * width_in_lcu
  const int32_t *tiles_ctb_addr_rs_to_ts; /*!<spec:  CtbAddrRsToTs (6.5.1); dimension: PicSizeInCtbsY */
  const int32_t *tiles_ctb_addr_ts_to_rs; /*!<spec:  CtbAddrTsToRs (6.5.1); dimension: PicSizeInCtbsY */
  
  const int32_t *tiles_tile_id; /*!<spec:  TileId (6.5.1); dimension: PicSizeInCtbsY */
  
  //WPP
  int wpp;
  
  //OWF 0 = no owf, 1 = 1 frame, 2 = 2 frames, etc.
  int owf;
  
  //Slices
  int slice_count;
  const int* slice_addresses_in_ts;
  
  threadqueue_queue_t *threadqueue;

  struct {
    uint8_t min;
    uint8_t max;
  } pu_depth_inter, pu_depth_intra;
  
  // How often Video Parameter Set is re-sent.
  int32_t vps_period;

  bool sign_hiding;

  //! Target average bits per picture.
  double target_avg_bppic;

  //! Target average bits per pixel.
  double target_avg_bpp;

  //! Picture weights when GOP is used.
  double gop_layer_weights[MAX_GOP_LAYERS];

} encoder_control_t;

encoder_control_t* kvz_encoder_control_init(const kvz_config *cfg);
void kvz_encoder_control_free(encoder_control_t *encoder);

void kvz_encoder_control_input_init(encoder_control_t *encoder, int32_t width, int32_t height);
unsigned kvz_get_padding(unsigned width_or_height);
#endif
