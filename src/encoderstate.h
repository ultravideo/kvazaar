#ifndef ENCODERSTATE_H_
#define ENCODERSTATE_H_
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
 * \brief
 */

#include "global.h"

#include "videoframe.h"
#include "encoder.h"
#include "image.h"
#include "bitstream.h"
#include "cabac.h"
#include "config.h"
#include "tables.h"
#include "scalinglist.h"
#include "threadqueue.h"
#include "imagelist.h"


// Submodules
// Functions to obtain geometry information from LCU
#include "encoder_state-geometry.h"
// Constructors/destructors
#include "encoder_state-ctors_dtors.h"
// Functions writing bitstream parts
#include "encoder_state-bitstream.h"


typedef enum {
  ENCODER_STATE_TYPE_INVALID = 'i',
  ENCODER_STATE_TYPE_MAIN = 'M',
  ENCODER_STATE_TYPE_SLICE = 'S',
  ENCODER_STATE_TYPE_TILE = 'T',
  ENCODER_STATE_TYPE_WAVEFRONT_ROW = 'W',
} encoder_state_type;



typedef struct {
  double cur_lambda_cost; //!< \brief Lambda for SSE
  double cur_lambda_cost_sqrt; //!< \brief Lambda for SAD and SATD
  
  int32_t frame;
  int32_t poc; /*!< \brief picture order count */
  int8_t gop_offset; /*!< \brief offset in the gop structure */
  
  int8_t QP;   //!< \brief Quantization parameter
  double QP_factor; //!< \brief Quantization factor
  
  //Current picture available references
  image_list_t *ref;
  int8_t ref_list;

  struct {
    int32_t poc;
    int8_t list;
    int8_t idx;
  } refmap[16];
  
  bool is_idr_frame;
  uint8_t pictype;
  uint8_t slicetype;

  //! Total number of bits written.
  uint64_t total_bits_coded;

  //! Number of bits written in the current GOP.
  uint64_t cur_gop_bits_coded;

  //! Number of bits targeted for the current GOP.
  double cur_gop_target_bits;

  // Parameters used in rate control
  double rc_alpha;
  double rc_beta;

} encoder_state_config_global_t;

typedef struct {
  //Current sub-frame
  videoframe_t *frame;
  
  int32_t id;
  
  //Tile: offset in LCU for current encoder_state in global coordinates
  int32_t lcu_offset_x;
  int32_t lcu_offset_y;
  
  //Position of the first element in tile scan in global coordinates
  int32_t lcu_offset_in_ts;
  
  //Buffer for search
  //order by row of (LCU_WIDTH * cur_pic->width_in_lcu) pixels
  yuv_t *hor_buf_search;
  //order by column of (LCU_WIDTH * encoder_state->height_in_lcu) pixels (there is no more extra pixel, since we can use a negative index)
  yuv_t *ver_buf_search;
  
  yuv_t *hor_buf_before_sao;
  yuv_t *ver_buf_before_sao;
  
  //Jobs for each individual LCU of a wavefront row.
  threadqueue_job_t **wf_jobs;
} encoder_state_config_tile_t;

typedef struct {
  int32_t id;
  
  //Global coordinates
  int32_t start_in_ts;
  int32_t end_in_ts;
  
  //Global coordinates
  int32_t start_in_rs;
  int32_t end_in_rs;
} encoder_state_config_slice_t;

typedef struct {
  //Row in tile coordinates of the wavefront
  int32_t lcu_offset_y;
} encoder_state_config_wfrow_t;

typedef struct lcu_order_element {
  //This it used for leaf of the encoding tree. All is relative to the tile.
  int id;
  int index;
  struct encoder_state_t *encoder_state;
  vector2d_t position;
  vector2d_t position_px; //Top-left
  vector2d_t size;
  int first_column;
  int first_row;
  int last_column;
  int last_row;
  
  struct lcu_order_element *above;
  struct lcu_order_element *below;
  struct lcu_order_element *left;
  struct lcu_order_element *right;
} lcu_order_element_t;

typedef struct encoder_state_t {
  const encoder_control_t *encoder_control;
  encoder_state_type type;

  //List of children, the last item of this list is a pseudo-encoder with encoder_control = NULL
  //Use for (i = 0; encoder_state->children[i].encoder_control; ++i) {
  struct encoder_state_t *children;
  struct encoder_state_t *parent;
  
  //Pointer to the encoder_state of the previous frame
  struct encoder_state_t *previous_encoder_state;
  
  encoder_state_config_global_t *global;
  encoder_state_config_tile_t   *tile;
  encoder_state_config_slice_t  *slice;
  encoder_state_config_wfrow_t  *wfrow;
  
  int is_leaf; //A leaf encoder state is one which should encode LCUs...
  lcu_order_element_t *lcu_order;
  uint32_t lcu_order_count;
  
  bitstream_t stream;
  cabac_data_t cabac;

  /**
   * \brief Indicates that this encoder state is ready for encoding the
   * next frame i.e. encoder_next_frame has been called.
   */
  int prepared;

  /**
   * \brief Indicates that the previous frame has been encoded and the
   * encoded data written and the encoding the next frame has not been
   * started yet.
   */
  int frame_done;

  uint32_t stats_bitstream_length; //Bitstream length written in bytes
  
  //Jobs to wait for
  threadqueue_job_t * tqj_recon_done; //Reconstruction is done
  threadqueue_job_t * tqj_bitstream_written; //Bitstream is written
} encoder_state_t;



void encode_one_frame(encoder_state_t *state);
int encoder_feed_frame(encoder_state_t* const state, image_t* const img_in);

void encoder_compute_stats(encoder_state_t *state, FILE * const recout, double psnr[3], uint64_t *bitstream_length);
void encoder_next_frame(encoder_state_t *state);


void encode_coding_tree(encoder_state_t *state, uint16_t x_ctb,
                        uint16_t y_ctb, uint8_t depth);

void encode_last_significant_xy(encoder_state_t *state,
                                uint8_t lastpos_x, uint8_t lastpos_y,
                                uint8_t width, uint8_t height,
                                uint8_t type, uint8_t scan);
void encode_coeff_nxn(encoder_state_t *state, int16_t *coeff, uint8_t width,
                      uint8_t type, int8_t scan_mode, int8_t tr_skip);
void encode_transform_coeff(encoder_state_t *state, int32_t x_cu, int32_t y_cu,
                            int8_t depth, int8_t tr_depth, uint8_t parent_coeff_u, uint8_t parent_coeff_v);
void encode_block_residual(const encoder_control_t * const encoder,
                           uint16_t x_ctb, uint16_t y_ctb, uint8_t depth);

int encoder_state_match_children_of_previous_frame(encoder_state_t * const state);

coeff_scan_order_t get_scan_order(int8_t cu_type, int intra_mode, int depth);

void encoder_ref_insertion_sort(int reflist[16], int length);

static const uint8_t g_group_idx[32] = {
  0, 1, 2, 3, 4, 4, 5, 5, 6, 6,
  6, 6, 7, 7, 7, 7, 8, 8, 8, 8,
  8, 8, 8, 8, 9, 9, 9, 9, 9, 9,
  9, 9 };

static const uint8_t g_min_in_group[10] = {
  0, 1, 2, 3, 4, 6, 8, 12, 16, 24 };


#define C1FLAG_NUMBER 8 // maximum number of largerThan1 flag coded in one chunk
#define C2FLAG_NUMBER 1 // maximum number of largerThan2 flag coded in one chunk

//Get the data for vertical buffer position at the left of LCU identified by the position in pixel
#define OFFSET_VER_BUF(position_x, position_y, cur_pic, i) ((position_y) + i + ((position_x)/LCU_WIDTH - 1) * (cur_pic)->height)
#define OFFSET_VER_BUF_C(position_x, position_y, cur_pic, i) ((position_y/2) + i + ((position_x)/LCU_WIDTH - 1) * (cur_pic)->height / 2)

//Get the data for horizontal buffer position at the top of LCU identified by the position in pixel
#define OFFSET_HOR_BUF(position_x, position_y, cur_pic, i) ((position_x) + i + ((position_y)/LCU_WIDTH - 1) * (cur_pic)->width)
#define OFFSET_HOR_BUF_C(position_x, position_y, cur_pic, i) ((position_x/2) + i + ((position_y)/LCU_WIDTH - 1) * (cur_pic)->width / 2)
  

#endif //ENCODERSTATE_H_
