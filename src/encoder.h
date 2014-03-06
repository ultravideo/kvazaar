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
#include "config.h"


/* TODO: add ME data */
typedef struct
{
  void (*IME)();
  void (*FME)();
  int range;
} encoder_me;

enum { FORMAT_400 = 0, FORMAT_420, FORMAT_422, FORMAT_444 };

/* Input info struct */
typedef struct
{
  FILE *file;
  int32_t width;  /*!< \brief input picture width */
  int32_t height; /*!< \brief input picture height */
  int32_t real_width;  /*!< \brief real input picture width */
  int32_t real_height; /*!< \brief real input picture width */
  int32_t height_in_lcu; /*!< \brief input picture width in LCU*/
  int32_t width_in_lcu;  /*!< \brief input picture height in LCU */
  picture *cur_pic;
  int8_t video_format;
  int8_t bitdepth;  /*!< \brief input bit depth (8,10) */
} encoder_input;

/* Encoder control options, the main struct */
typedef struct
{
  int32_t frame;
  int32_t poc; /*!< \brief picture order count */
  config *cfg;
  encoder_input in;
  encoder_me me;
  bitstream *stream;
  FILE *output;
  picture_list *ref;
  int8_t ref_list;
  int8_t ref_idx_num[2];
  int8_t QP;             // \brief Quantization parameter
  int8_t bitdepth;

  /* Filtering */
  int8_t deblock_enable; // \brief Flag to enable deblocking filter
  int8_t sao_enable;     // \brief Flag to enable sample adaptive offset filter
  int8_t rdoq_enable;    // \brief Whether RDOQ is enabled or not.
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
  FILE *cqmfile; // \brief Costum Quantization Matrices
} encoder_control;

void init_tables(void);
void init_lambda(encoder_control* encoder);
void free_tables(void);
encoder_control *init_encoder_control(config *cfg);
void init_encoder_input(encoder_input *input, FILE* inputfile,
                        int32_t width, int32_t height);
void encode_one_frame(encoder_control *encoder);
int read_one_frame(FILE *file, encoder_control *encoder);

void encode_seq_parameter_set(encoder_control *encoder);
void encode_pic_parameter_set(encoder_control *encoder);
void encode_vid_parameter_set(encoder_control *encoder);
void encode_slice_data(encoder_control *encoder);
void encode_slice_header(encoder_control *encoder);
void encode_access_unit_delimiter(encoder_control* encoder);
void encode_prefix_sei_version(encoder_control* encoder);
void encode_coding_tree(encoder_control *encoder, uint16_t x_ctb,
                        uint16_t y_ctb, uint8_t depth);
void encode_last_significant_xy(encoder_control *encoder, uint8_t lastpos_x,
                                uint8_t lastpos_y, uint8_t width, uint8_t height,
                                uint8_t type, uint8_t scan);
void encode_coeff_nxn(encoder_control *encoder, int16_t *coeff, uint8_t width,
                      uint8_t type, int8_t scan_mode);
void encode_transform_tree(encoder_control* encoder, int32_t x, int32_t y, uint8_t depth, lcu_t *lcu);
void encode_transform_coeff(encoder_control *encoder, int32_t x_cu, int32_t y_cu,
                            int8_t depth, int8_t tr_depth, uint8_t parent_coeff_u, uint8_t parent_coeff_v);
void encode_block_residual(encoder_control *encoder,
                           uint16_t x_ctb, uint16_t y_ctb, uint8_t depth);

extern double g_lambda_cost[55];
extern double g_cur_lambda_cost;
extern uint32_t* g_sig_last_scan[3][7];
extern int8_t g_convert_to_bit[LCU_WIDTH + 1];
static int8_t g_bitdepth     = 8;
static int8_t g_bit_increment = 0;

#define MAX_NUM_SPU_W ((1 << (MAX_DEPTH)) / 4)
static uint32_t g_z_scan_to_raster[MAX_NUM_SPU_W * MAX_NUM_SPU_W] = { 0, };
static uint32_t g_raster_to_z_scan[MAX_NUM_SPU_W * MAX_NUM_SPU_W] = { 0, };

static const uint8_t g_group_idx[32] = {
  0, 1, 2, 3, 4, 4, 5, 5, 6, 6,
  6, 6, 7, 7, 7, 7, 8, 8, 8, 8,
  8, 8, 8, 8, 9, 9, 9, 9, 9, 9,
  9, 9 };

static const uint8_t g_min_in_group[10] = {
  0, 1, 2, 3, 4, 6, 8, 12, 16, 24 };

static uint32_t g_sig_last_scan_32x32[64] = {
  0,  8,  1,  16, 9,  2,  24, 17, 10, 3,
  32, 25, 18, 11, 4,  40, 33, 26, 19, 12,
  5,  48, 41, 34, 27, 20, 13, 6,  56, 49,
  42, 35, 28, 21, 14, 7,  57, 50, 43, 36,
  29, 22, 15, 58, 51, 44, 37, 30, 23, 59,
  52, 45, 38, 31, 60, 53, 46, 39, 61, 54,
  47, 62, 55, 63 };

static const uint32_t g_sig_last_scan_8x8[3][4] =
{ {0, 2, 1, 3},
  {0, 1, 2, 3},
  {0, 2, 1, 3}
};


//4 8 16 32 64 128
//0 1  2  3  4   5
static const uint8_t g_to_bits[129] =
{
  0,
  0,0,0,0,
  0,0,0,1,
  0,0,0,0,0,0,0,2,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
};
#define TOBITS(len) g_to_bits[len]


#define C1FLAG_NUMBER 8 // maximum number of largerThan1 flag coded in one chunk
#define C2FLAG_NUMBER 1 // maximum number of largerThan2 flag coded in one chunk

enum COEFF_SCAN_TYPE
{
  SCAN_DIAG = 0, // up-right diagonal scan
  SCAN_HOR,      // horizontal first scan
  SCAN_VER       // vertical first scan
};


#endif
