#ifndef ENCODER_H_
#define ENCODER_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

/**
 * \ingroup Control
 * \file
 * Initialization of encoder_control_t.
 */

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "scalinglist.h"
#include "threadqueue.h"
#include "fast_coeff_cost.h"

/* Encoder control options, the main struct */
typedef struct encoder_control_t
{
  /**
   * \brief Configuration.
   *
   * NOTE: The following fields are not copied from the config passed to
   * kvz_encoder_control_init and must not be accessed:
   *    - cqmfile
   *    - tiles_width_split
   *    - tiles_height_split
   *    - slice_addresses_in_ts
   * Use appropriate fields in encoder_control_t instead.
   */
  kvz_config cfg;

  /* Input */
  struct {
    int32_t width;
    int32_t height;
    int32_t width_in_lcu;
    int32_t height_in_lcu;
    int32_t real_width;  /*!< \brief real input picture width */
    int32_t real_height; /*!< \brief real input picture height */
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
  enum kvz_chroma_format chroma_format;

  /* VUI */
  struct
  {
    /* Timing info */
    int32_t num_units_in_tick; /*!< \brief Timing scale numerator */
    int32_t time_scale; /*!< \brief Timing scale denominator */

    int8_t field_seq_flag;
    int8_t frame_field_info_present_flag;

    int8_t timing_info_present_flag;
  } vui;

  //scaling list
  scaling_list_t scaling_list;

  //spec: references to variables defined in Rec. ITU-T H.265 (04/2013)
  int8_t tiles_enable; /*!<spec: tiles_enabled */

  int8_t tiles_uniform_spacing_flag; /*!<spec: uniform_spacing_flag */

  const int32_t *tiles_col_width; /*!<spec: colWidth (6.5.1); dimension: tiles_num_tile_columns */
  const int32_t *tiles_row_height; /*!<spec: rowHeight (6.5.1); dimension: tiles_num_tile_rows */

  const int32_t *tiles_col_bd; /*!<spec: colBd (6.5.1); dimension: tiles_num_tile_columns + 1 */
  const int32_t *tiles_row_bd; /*!<spec: rowBd (6.5.1); dimension: tiles_num_tile_rows + 1  */

  //PicSizeInCtbsY = height_in_lcu * width_in_lcu
  const int32_t *tiles_ctb_addr_rs_to_ts; /*!<spec:  CtbAddrRsToTs (6.5.1); dimension: PicSizeInCtbsY */
  const int32_t *tiles_ctb_addr_ts_to_rs; /*!<spec:  CtbAddrTsToRs (6.5.1); dimension: PicSizeInCtbsY */

  const int32_t *tiles_tile_id; /*!<spec:  TileId (6.5.1); dimension: PicSizeInCtbsY */

  //Slices
  int slice_count;
  const int* slice_addresses_in_ts;

  threadqueue_queue_t *threadqueue;

  //! Target average bits per picture.
  double target_avg_bppic;

  //! Target average bits per pixel.
  double target_avg_bpp;

  //! Picture weights when GOP is used.
  double gop_layer_weights[MAX_GOP_LAYERS];

  FILE *roi_file;

  int tr_depth_inter;

  //! pic_parameter_set
  struct {
    uint8_t dependent_slice_segments_enabled_flag;
  } pps;

  //! Maximum motion vector distance as number of LCUs.
  struct {
    int right;
    int down;
  } max_inter_ref_lcu;

  int32_t poc_lsb_bits;

  fast_coeff_table_t fast_coeff_table;

} encoder_control_t;

encoder_control_t* kvz_encoder_control_init(const kvz_config *cfg);
void kvz_encoder_control_free(encoder_control_t *encoder);

void kvz_encoder_control_input_init(encoder_control_t *encoder, int32_t width, int32_t height);
#endif
