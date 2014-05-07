#ifndef ENCODER_H_
#define ENCODER_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2014 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

/*
 * \file
 * \brief The highest level of the encoder.
 */

#include "global.h"

#include "picture.h"
#include "bitstream.h"
#include "cabac.h"
#include "config.h"
#include "tables.h"
#include "scalinglist.h"


/* TODO: add ME data */
typedef struct
{
  void (*IME)();
  void (*FME)();
  int range;
} encoder_me;

enum { FORMAT_400 = 0, FORMAT_420, FORMAT_422, FORMAT_444 };

/* Encoder control options, the main struct */
typedef struct
{
  /* Configuration */
  const config *cfg;
  
  /* Input */
  struct {
    FILE *file;
    int32_t width;
    int32_t height;
    int32_t width_in_lcu;
    int32_t height_in_lcu;
    int32_t real_width;  /*!< \brief real input picture width */
    int32_t real_height; /*!< \brief real input picture width */
    int8_t video_format;
    int8_t bitdepth;  /*!< \brief input bit depth (8,10) */
  } in;
  
  /* Output */
  struct {
    FILE *file;
  } out;
  
  encoder_me me;
  
  int8_t bitdepth;

  /* Filtering */
  int8_t deblock_enable; // \brief Flag to enable deblocking filter
  int8_t sao_enable;     // \brief Flag to enable sample adaptive offset filter
  int8_t rdoq_enable;    // \brief Whether RDOQ is enabled or not.
  int8_t rdo;            // \brief RDO level
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
  } vui;

  int8_t aud_enable;

  //scaling list
  scaling_list scaling_list;
  
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
  
  //Slices
  int slice_count;
  const int* slice_addresses_in_ts;
  
} encoder_control;

typedef enum {
  ENCODER_STATE_TYPE_INVALID = 'i',
  ENCODER_STATE_TYPE_MAIN = 'M',
  ENCODER_STATE_TYPE_SLICE = 'S',
  ENCODER_STATE_TYPE_TILE = 'T',
  ENCODER_STATE_TYPE_WAVEFRONT_ROW = 'W',
} encoder_state_type;



typedef struct {
  double cur_lambda_cost;
  
  int32_t frame;
  int32_t poc; /*!< \brief picture order count */
  
  int8_t QP;   //!< \brief Quantization parameter
  
  //Current picture available references
  picture_list *ref;
  int8_t ref_list;
  //int8_t ref_idx_num[2];
  
} encoder_state_config_global;

typedef struct {
  //Current picture to encode
  picture *cur_pic;
  
  int32_t id;
  
  //Tile: offset in LCU for current encoder_state in global coordinates
  int32_t lcu_offset_x;
  int32_t lcu_offset_y;
  
  //Position of the first element in tile scan in global coordinates
  int32_t lcu_offset_in_ts;
} encoder_state_config_tile;

typedef struct {
  int32_t id;
  
  //Global coordinates
  int32_t start_in_ts;
  int32_t end_in_ts;
  
  //Global coordinates
  int32_t start_in_rs;
  int32_t end_in_rs;
} encoder_state_config_slice;

typedef struct {
  //Row in image coordinates of the wavefront
  int32_t lcu_offset_y;
} encoder_state_config_wfrow;

typedef struct encoder_state {
  const encoder_control *encoder_control;
  encoder_state_type type;

  //List of children, the last item of this list is a pseudo-encoder with encoder_control = NULL
  //Use for (i = 0; encoder_state->children[i].encoder_control; ++i) {
  struct encoder_state *children;
  struct encoder_state *parent;
  
  encoder_state_config_global *global;
  encoder_state_config_tile   *tile;
  encoder_state_config_slice  *slice;
  encoder_state_config_wfrow  *wfrow;
  
  bitstream stream;
  cabac_data cabac;
} encoder_state;

int encoder_control_init(encoder_control *encoder, const config *cfg);
int encoder_control_finalize(encoder_control *encoder);

void encoder_control_input_init(encoder_control *encoder, int32_t width, int32_t height);

int encoder_state_init(encoder_state * child_state, encoder_state * parent_state);
void encoder_state_finalize(encoder_state *encoder_state);
void encoder_state_init_lambda(encoder_state *encoder_state);

void encode_one_frame(encoder_state *encoder_state);
int read_one_frame(FILE* file, const encoder_state *encoder);

void encoder_next_frame(encoder_state *encoder_state);

void encode_seq_parameter_set(encoder_state *encoder);
void encode_pic_parameter_set(encoder_state *encoder);
void encode_vid_parameter_set(encoder_state *encoder);
void encode_slice_header(encoder_state * encoder);
void encode_access_unit_delimiter(encoder_state *encoder);
void encode_prefix_sei_version(encoder_state *encoder);
void encode_coding_tree(encoder_state *encoder, uint16_t x_ctb,
                        uint16_t y_ctb, uint8_t depth);

void encode_last_significant_xy(encoder_state *encoder,
                                uint8_t lastpos_x, uint8_t lastpos_y,
                                uint8_t width, uint8_t height,
                                uint8_t type, uint8_t scan);
void encode_coeff_nxn(encoder_state *encoder, int16_t *coeff, uint8_t width,
                      uint8_t type, int8_t scan_mode, int8_t tr_skip);
void encode_transform_tree(encoder_state *encoder_state, int32_t x, int32_t y, uint8_t depth, lcu_t* lcu );
void encode_transform_coeff(encoder_state *encoder_state, int32_t x_cu, int32_t y_cu,
                            int8_t depth, int8_t tr_depth, uint8_t parent_coeff_u, uint8_t parent_coeff_v);
void encode_block_residual(const encoder_control * const encoder,
                           uint16_t x_ctb, uint16_t y_ctb, uint8_t depth);

static const uint8_t g_group_idx[32] = {
  0, 1, 2, 3, 4, 4, 5, 5, 6, 6,
  6, 6, 7, 7, 7, 7, 8, 8, 8, 8,
  8, 8, 8, 8, 9, 9, 9, 9, 9, 9,
  9, 9 };

static const uint8_t g_min_in_group[10] = {
  0, 1, 2, 3, 4, 6, 8, 12, 16, 24 };




#define C1FLAG_NUMBER 8 // maximum number of largerThan1 flag coded in one chunk
#define C2FLAG_NUMBER 1 // maximum number of largerThan2 flag coded in one chunk



#endif
