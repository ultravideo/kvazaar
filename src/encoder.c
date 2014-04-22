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
 */

#include "encoder.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "tables.h"
#include "config.h"
#include "cabac.h"
#include "picture.h"
#include "nal.h"
#include "context.h"
#include "transform.h"
#include "intra.h"
#include "inter.h"
#include "filter.h"
#include "search.h"
#include "sao.h"
#include "rdo.h"

/* Local functions. */
static void add_checksum(encoder_state *encoder);
static void encode_VUI(encoder_state *encoder);
static void encode_sao(encoder_state *encoder,
                       cabac_data *cabac,
                       unsigned x_lcu, uint16_t y_lcu,
                       sao_info *sao_luma, sao_info *sao_chroma);

/*!
  \brief Initializes lambda-value for current QP

  Implementation closer to HM (Used HM12 as reference)
   - Still missing functionality when GOP and B-pictures are used
 */
void encoder_state_init_lambda(encoder_state * const encoder_state)
{
  const picture * const cur_pic = encoder_state->cur_pic;
  double qp = encoder_state->QP;
  double lambda_scale = 1.0;
  double qp_temp      = qp - 12;
  double lambda;

  // Default QP-factor from HM config
  double qp_factor = 0.4624;

  if (cur_pic->slicetype == SLICE_I) {
    qp_factor=0.57*lambda_scale;
  }

  lambda = qp_factor*pow( 2.0, qp_temp/3.0 );

  if (cur_pic->slicetype != SLICE_I ) {
    lambda *= 0.95;
  }

  encoder_state->cur_lambda_cost = lambda;
}

int encoder_control_init(encoder_control * const encoder, const config * const cfg) {
  if (!cfg) {
    fprintf(stderr, "Config object must not be null!\n");
    return 0;
  }
  
  // Config pointer to config struct
  encoder->cfg = cfg;
  encoder->bitdepth = 8;
  
  // deblocking filter
  encoder->deblock_enable    = 1;
  encoder->beta_offset_div2  = 0;
  encoder->tc_offset_div2    = 0;
  // SAO
  encoder->sao_enable = 1;
  // Rate-distortion optimization level
  encoder->rdo        = 1;
  
  // Initialize the scaling list
  scalinglist_init(&encoder->scaling_list);
  
  // CQM
  {
    FILE* cqmfile;
    cqmfile = cfg->cqmfile ? fopen(cfg->cqmfile, "rb") : NULL;
    if (cqmfile) {
      scalinglist_parse(&encoder->scaling_list, cqmfile);
      fclose(cqmfile);
    }
  }
  scalinglist_process(&encoder->scaling_list, encoder->bitdepth);
  
  return 1;
}

int encoder_control_finalize(encoder_control * const encoder) {
  scalinglist_destroy(&encoder->scaling_list);
  
  return 1;
}

int encoder_state_init(encoder_state * const encoder_state, const encoder_control * const encoder) {
  picture_list    *pic_list = NULL;
  
  encoder_state->encoder_control = encoder;
  
  // Allocate the bitstream struct
  if (!bitstream_init(&encoder_state->stream, BITSTREAM_TYPE_FILE)) {
    fprintf(stderr, "Could not initialize stream!\n");
    return 0;
  }
  
  pic_list = picture_list_init(MAX_REF_PIC_COUNT);
  if(!pic_list) {
    fprintf(stderr, "Failed to allocate the picture list!\n");
    return 0;
  }
  
  encoder_state->ref = pic_list;
  encoder_state->ref_list = REF_PIC_LIST_0;
  
  encoder_state->frame = 0;
  encoder_state->poc = 0;
  
  // Allocate the picture and CU array
  encoder_state->cur_pic = picture_init(encoder->in.width, encoder->in.height,
                                encoder->in.width_in_lcu, encoder->in.height_in_lcu);

  if (!encoder_state->cur_pic) {
    printf("Error allocating picture!\r\n");
    return 0;
  }
  
  // Init coeff data table
  encoder_state->cur_pic->coeff_y = MALLOC(coefficient, encoder->in.width * encoder->in.height);
  encoder_state->cur_pic->coeff_u = MALLOC(coefficient, (encoder->in.width * encoder->in.height) >> 2);
  encoder_state->cur_pic->coeff_v = MALLOC(coefficient, (encoder->in.width * encoder->in.height) >> 2);

  // Init predicted data table
  encoder_state->cur_pic->pred_y = MALLOC(pixel, encoder->in.width * encoder->in.height);
  encoder_state->cur_pic->pred_u = MALLOC(pixel, (encoder->in.width * encoder->in.height) >> 2);
  encoder_state->cur_pic->pred_v = MALLOC(pixel, (encoder->in.width * encoder->in.height) >> 2);
  
  encoder_state->children = NULL;
  
  encoder_state->stream.file.output = encoder->out.file;
  
  // Set CABAC output bitstream
  encoder_state->cabac.stream = &encoder_state->stream;
  
  return 1;
}

int encoder_state_finalize(encoder_state * const encoder_state) {
  picture_destroy(encoder_state->cur_pic);
  FREE_POINTER(encoder_state->cur_pic);
  
  picture_list_destroy(encoder_state->ref);
  bitstream_finalize(&encoder_state->stream);
  return 1;
}

void encoder_control_input_init(encoder_control * const encoder, FILE *inputfile,
                        const int32_t width, const int32_t height)
{
  encoder->in.file = inputfile;
  encoder->in.width = width;
  encoder->in.height = height;
  encoder->in.real_width = width;
  encoder->in.real_height = height;

  // If input dimensions are not divisible by the smallest block size, add
  // pixels to the dimensions, so that they are. These extra pixels will be
  // compressed along with the real ones but they will be cropped out before
  // rendering.
  if (encoder->in.width % CU_MIN_SIZE_PIXELS) {
    encoder->in.width += CU_MIN_SIZE_PIXELS - (width % CU_MIN_SIZE_PIXELS);
  }

  if (encoder->in.height % CU_MIN_SIZE_PIXELS) {
    encoder->in.height += CU_MIN_SIZE_PIXELS - (height % CU_MIN_SIZE_PIXELS);
  }

  encoder->in.height_in_lcu = encoder->in.height / LCU_WIDTH;
  encoder->in.width_in_lcu  = encoder->in.width / LCU_WIDTH;

  // Add one extra LCU when image not divisible by LCU_WIDTH
  if (encoder->in.height_in_lcu * LCU_WIDTH < height) {
    encoder->in.height_in_lcu++;
  }

  if (encoder->in.width_in_lcu * LCU_WIDTH < width) {
    encoder->in.width_in_lcu++;
  }



  #ifdef _DEBUG
  if (width != i_width || height != i_height) {
    printf("Picture buffer has been extended to be a multiple of the smallest block size:\r\n");
    printf("  Width = %d (%d), Height = %d (%d)\r\n", width, encoder->in.width, height,
           encoder->in.height);
  }
  #endif
}

static void write_aud(encoder_state * const encoder_state)
{
  bitstream * const stream = &encoder_state->stream;
  encode_access_unit_delimiter(encoder_state);
  nal_write(stream, AUD_NUT, 0, 1);
  bitstream_align(stream);
}

void encode_one_frame(encoder_state * const encoder_state)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  bitstream * const stream = &encoder_state->stream;
  
  yuv_t *hor_buf = alloc_yuv_t(encoder_state->cur_pic->width);
  // Allocate 2 extra luma pixels so we get 1 extra chroma pixel for the
  // for the extra pixel on the top right.
  yuv_t *ver_buf = alloc_yuv_t(LCU_WIDTH + 2);

  const int is_first_frame = (encoder_state->frame == 0);
  const int is_i_radl = (encoder->cfg->intra_period == 1 && encoder_state->frame % 2 == 0);
  const int is_p_radl = (encoder->cfg->intra_period > 1 && (encoder_state->frame % encoder->cfg->intra_period) == 0);
  const int is_radl_frame = is_first_frame || is_i_radl || is_p_radl;


  /** IDR picture when: period == 0 and frame == 0
   *                    period == 1 && frame%2 == 0
   *                    period != 0 && frame%period == 0
   **/
  if (is_radl_frame) {
    // Clear the reference list
    while (encoder_state->ref->used_size) {
      picture_list_rem(encoder_state->ref, encoder_state->ref->used_size - 1, 1);
    }

    encoder_state->poc = 0;

    encoder_state->cur_pic->slicetype = SLICE_I;
    encoder_state->cur_pic->type = NAL_IDR_W_RADL;

    // Access Unit Delimiter (AUD)
    if (encoder->aud_enable)
      write_aud(encoder_state);

    // Video Parameter Set (VPS)
    nal_write(stream, NAL_VPS_NUT, 0, 1);
    encode_vid_parameter_set(encoder_state);
    bitstream_align(stream);

    // Sequence Parameter Set (SPS)
    nal_write(stream, NAL_SPS_NUT, 0, 1);
    encode_seq_parameter_set(encoder_state);
    bitstream_align(stream);

    // Picture Parameter Set (PPS)
    nal_write(stream, NAL_PPS_NUT, 0, 1);
    encode_pic_parameter_set(encoder_state);
    bitstream_align(stream);

    if (encoder_state->frame == 0) {
      // Prefix SEI
      nal_write(stream, PREFIX_SEI_NUT, 0, 0);
      encode_prefix_sei_version(encoder_state);
      bitstream_align(stream);
    }
  } else {
    // When intra period == 1, all pictures are intra
    encoder_state->cur_pic->slicetype = encoder->cfg->intra_period==1 ? SLICE_I : SLICE_P;
    encoder_state->cur_pic->type = NAL_TRAIL_R;

    // Access Unit Delimiter (AUD)
    if (encoder->aud_enable)
      write_aud(encoder_state);
  }

  {
    // Not quite sure if this is correct, but it seems to have worked so far
    // so I tried to not change it's behavior.
    int long_start_code = is_radl_frame || encoder->aud_enable ? 0 : 1;

    nal_write(stream,
              is_radl_frame ? NAL_IDR_W_RADL : NAL_TRAIL_R, 0, long_start_code);
  }

  cabac_start(&encoder_state->cabac);
  init_contexts(&encoder_state->cabac, encoder_state->QP, encoder_state->cur_pic->slicetype);
  encode_slice_header(encoder_state);
  bitstream_align(stream);

  // Initialize lambda value(s) to use in search
  encoder_state_init_lambda(encoder_state);

  {
    picture* const cur_pic = encoder_state->cur_pic;
    vector2d lcu;
    const vector2d size = { cur_pic->width, cur_pic->height };
    const vector2d size_lcu = { cur_pic->width_in_lcu, cur_pic->height_in_lcu };

    for (lcu.y = 0; lcu.y < cur_pic->height_in_lcu; lcu.y++) {
      for (lcu.x = 0; lcu.x < cur_pic->width_in_lcu; lcu.x++) {
        const vector2d px = { lcu.x * LCU_WIDTH, lcu.y * LCU_WIDTH };

        // Handle partial LCUs on the right and bottom.
        const vector2d lcu_dim = {
          MIN(LCU_WIDTH, size.x - px.x), MIN(LCU_WIDTH, size.y - px.y)
        };
        const int right = px.x + lcu_dim.x;
        const int bottom = px.y + lcu_dim.y;

        search_lcu(encoder_state, &encoder_state->cabac, px.x, px.y, hor_buf, ver_buf);

        // Take the bottom right pixel from the LCU above and put it as the
        // first pixel in this LCUs rightmost pixels.
        if (lcu.y > 0) {
          ver_buf->y[0] = hor_buf->y[right - 1];
          ver_buf->u[0] = hor_buf->u[right / 2 - 1];
          ver_buf->v[0] = hor_buf->v[right / 2 - 1];
        }

        // Take bottom and right pixels from this LCU to be used on the search of next LCU.
        picture_blit_pixels(&cur_pic->y_recdata[(bottom - 1) * size.x + px.x],
                            &hor_buf->y[px.x],
                            lcu_dim.x, 1, size.x, size.x);
        picture_blit_pixels(&cur_pic->u_recdata[(bottom / 2 - 1) * size.x / 2 + px.x / 2],
                            &hor_buf->u[px.x / 2],
                            lcu_dim.x / 2, 1, size.x / 2, size.x / 2);
        picture_blit_pixels(&cur_pic->v_recdata[(bottom / 2 - 1) * size.x / 2 + px.x / 2],
                            &hor_buf->v[px.x / 2],
                            lcu_dim.x / 2, 1, size.x / 2, size.x / 2);

        picture_blit_pixels(&cur_pic->y_recdata[px.y * size.x + right - 1],
                            &ver_buf->y[1],
                            1, lcu_dim.y, size.x, 1);
        picture_blit_pixels(&cur_pic->u_recdata[px.y * size.x / 4 + (right / 2) - 1],
                            &ver_buf->u[1],
                            1, lcu_dim.y / 2, size.x / 2, 1);
        picture_blit_pixels(&cur_pic->v_recdata[px.y * size.x / 4 + (right / 2) - 1],
                            &ver_buf->v[1],
                            1, lcu_dim.y / 2, size.x / 2, 1);

        if (encoder->deblock_enable) {
          filter_deblock_lcu(encoder_state, px.x, px.y);
        }

        if (encoder->sao_enable) {
          const int stride = cur_pic->width_in_lcu;
          sao_info *sao_luma = &cur_pic->sao_luma[lcu.y * stride + lcu.x];
          sao_info *sao_chroma = &cur_pic->sao_chroma[lcu.y * stride + lcu.x];
          init_sao_info(sao_luma);
          init_sao_info(sao_chroma);

          {
            sao_info *sao_top = lcu. y != 0 ? &cur_pic->sao_luma[(lcu.y - 1) * stride + lcu.x] : NULL;
            sao_info *sao_left = lcu.x != 0 ? &cur_pic->sao_luma[lcu.y * stride + lcu.x -1] : NULL;
            sao_search_luma(encoder_state, cur_pic, lcu.x, lcu.y, sao_luma, sao_top, sao_left);
          }

          {
            sao_info *sao_top = lcu.y != 0 ? &cur_pic->sao_chroma[(lcu.y - 1) * stride + lcu.x] : NULL;
            sao_info *sao_left = lcu.x != 0 ? &cur_pic->sao_chroma[lcu.y * stride + lcu.x - 1] : NULL;
            sao_search_chroma(encoder_state, cur_pic, lcu.x, lcu.y, sao_chroma, sao_top, sao_left);
          }

          // Merge only if both luma and chroma can be merged
          sao_luma->merge_left_flag = sao_luma->merge_left_flag & sao_chroma->merge_left_flag;
          sao_luma->merge_up_flag = sao_luma->merge_up_flag & sao_chroma->merge_up_flag;

          encode_sao(encoder_state, &encoder_state->cabac, lcu.x, lcu.y, sao_luma, sao_chroma);
        }
        
        encode_coding_tree(encoder_state, &encoder_state->cabac, lcu.x << MAX_DEPTH, lcu.y << MAX_DEPTH, 0);

        {
          const int last_lcu = (lcu.x == size_lcu.x - 1 && lcu.y == size_lcu.y - 1);
          cabac_encode_bin_trm(&encoder_state->cabac, last_lcu ? 1 : 0);  // end_of_slice_segment_flag
        }
      }
    }
  }

  cabac_flush(&encoder_state->cabac);
  bitstream_align(stream);

  if (encoder->sao_enable) {
    sao_reconstruct_frame(encoder_state);
  }

  // Calculate checksum
  add_checksum(encoder_state);

  encoder_state->cur_pic->poc = encoder_state->poc;

  dealloc_yuv_t(hor_buf);
  dealloc_yuv_t(ver_buf);
}

static void fill_after_frame(unsigned height, unsigned array_width,
                             unsigned array_height, pixel *data)
{
  pixel* p = data + height * array_width;
  pixel* end = data + array_width * array_height;

  while (p < end) {
    // Fill the line by copying the line above.
    memcpy(p, p - array_width, array_width);
    p += array_width;
  }
}

static int read_and_fill_frame_data(FILE *file,
                                    unsigned width, unsigned height,
                                    unsigned array_width, pixel *data)
{
  pixel* p = data;
  pixel* end = data + array_width * height;
  pixel fill_char;
  unsigned i;

  while (p < end) {
    // Read the beginning of the line from input.
    if (width != fread(p, sizeof(unsigned char), width, file))
      return 0;

    // Fill the rest with the last pixel value.
    fill_char = p[width - 1];

    for (i = width; i < array_width; ++i) {
      p[i] = fill_char;
    }

    p += array_width;
  }
  return 1;
}

int read_one_frame(FILE* file, const encoder_state * const encoder_state)
{
  unsigned width = encoder_state->encoder_control->in.real_width;
  unsigned height = encoder_state->encoder_control->in.real_height;
  unsigned array_width = encoder_state->cur_pic->width;
  unsigned array_height = encoder_state->cur_pic->height;

  if (width != array_width) {
    // In the case of frames not being aligned on 8 bit borders, bits need to be copied to fill them in.
    if (!read_and_fill_frame_data(file, width, height, array_width,
                                  encoder_state->cur_pic->y_data) ||
        !read_and_fill_frame_data(file, width >> 1, height >> 1, array_width >> 1,
                                  encoder_state->cur_pic->u_data) ||
        !read_and_fill_frame_data(file, width >> 1, height >> 1, array_width >> 1,
                                  encoder_state->cur_pic->v_data))
      return 0;
  } else {
    // Otherwise the data can be read directly to the array.
    unsigned y_size = width * height;
    unsigned uv_size = (width >> 1) * (height >> 1);
    if (y_size  != fread(encoder_state->cur_pic->y_data, sizeof(unsigned char),
                         y_size, file) ||
        uv_size != fread(encoder_state->cur_pic->u_data, sizeof(unsigned char),
                         uv_size, file) ||
        uv_size != fread(encoder_state->cur_pic->v_data, sizeof(unsigned char),
                         uv_size, file))
      return 0;
  }

  if (height != array_height) {
    fill_after_frame(height, array_width, array_height,
                     encoder_state->cur_pic->y_data);
    fill_after_frame(height >> 1, array_width >> 1, array_height >> 1,
                     encoder_state->cur_pic->u_data);
    fill_after_frame(height >> 1, array_width >> 1, array_height >> 1,
                     encoder_state->cur_pic->v_data);
  }
  return 1;
}

/**
 * \brief Add a checksum SEI message to the bitstream.
 * \param encoder The encoder.
 * \returns Void
 */
static void add_checksum(encoder_state * const encoder_state)
{
  bitstream * const stream = &encoder_state->stream;
  const picture * const cur_pic = encoder_state->cur_pic;
  unsigned char checksum[3][SEI_HASH_MAX_LENGTH];
  uint32_t checksum_val;
  unsigned int i;

  nal_write(stream, NAL_SUFFIT_SEI_NUT, 0, 0);

  picture_checksum(cur_pic, checksum);

  WRITE_U(stream, 132, 8, "sei_type");
  WRITE_U(stream, 13, 8, "size");
  WRITE_U(stream, 2, 8, "hash_type"); // 2 = checksum

  for (i = 0; i < 3; ++i) {
    // Pack bits into a single 32 bit uint instead of pushing them one byte
    // at a time.
    checksum_val = (checksum[i][0] << 24) + (checksum[i][1] << 16) +
                   (checksum[i][2] << 8) + (checksum[i][3]);
    WRITE_U(stream, checksum_val, 32, "picture_checksum");
  }

  bitstream_align(stream);
}

void encode_access_unit_delimiter(encoder_state * const encoder_state)
{
  bitstream * const stream = &encoder_state->stream;
  const picture * const cur_pic = encoder_state->cur_pic;
  uint8_t pic_type = cur_pic->slicetype == SLICE_I ? 0
                   : cur_pic->slicetype == SLICE_P ? 1
                   :                                             2;
  WRITE_U(stream, pic_type, 3, "pic_type");
}

void encode_prefix_sei_version(encoder_state * const encoder_state)
{
#define STR_BUF_LEN 1000
  bitstream * const stream = &encoder_state->stream;
  int i, length;
  char buf[STR_BUF_LEN] = { 0 };
  char *s = buf + 16;
  const config * const cfg = encoder_state->encoder_control->cfg;

  // random uuid_iso_iec_11578 generated with www.famkruithof.net/uuid/uuidgen
  static const uint8_t uuid[16] = {
    0x32, 0xfe, 0x46, 0x6c, 0x98, 0x41, 0x42, 0x69,
    0xae, 0x35, 0x6a, 0x91, 0x54, 0x9e, 0xf3, 0xf1
  };
  memcpy(buf, uuid, 16);

  // user_data_payload_byte
  s += sprintf(s, "Kvazaar HEVC Encoder v. " VERSION_STRING " - "
                  "Copyleft 2012-2014 - http://ultravideo.cs.tut.fi/ - options:");
  s += sprintf(s, " %dx%d", cfg->width, cfg->height);
  s += sprintf(s, " deblock=%d:%d:%d", cfg->deblock_enable,
               cfg->deblock_beta, cfg->deblock_tc);
  s += sprintf(s, " sao=%d", cfg->sao_enable);
  s += sprintf(s, " intra_period=%d", cfg->intra_period);
  s += sprintf(s, " qp=%d", cfg->qp);
  s += sprintf(s, " ref=%d", cfg->ref_frames);

  length = (int)(s - buf + 1);  // length, +1 for \0

  // Assert this so that in the future if the message gets longer, we remember
  // to increase the buf len. Divide by 2 for margin.
  assert(length < STR_BUF_LEN / 2);

  // payloadType = 5 -> user_data_unregistered
  WRITE_U(stream, 5, 8, "last_payload_type_byte");

  // payloadSize
  for (i = 0; i <= length - 255; i += 255)
    WRITE_U(stream, 255, 8, "ff_byte");
  WRITE_U(stream, length - i, 8, "last_payload_size_byte");

  for (i = 0; i < length; i++)
    WRITE_U(stream, ((uint8_t *)buf)[i], 8, "sei_payload");

#undef STR_BUF_LEN
}

void encode_pic_parameter_set(encoder_state * const encoder_state)
{
  bitstream * const stream = &encoder_state->stream;
#ifdef _DEBUG
  printf("=========== Picture Parameter Set ID: 0 ===========\n");
#endif
  WRITE_UE(stream, 0, "pic_parameter_set_id");
  WRITE_UE(stream, 0, "seq_parameter_set_id");
  WRITE_U(stream, 0, 1, "dependent_slice_segments_enabled_flag");
  WRITE_U(stream, 0, 1, "output_flag_present_flag");
  WRITE_U(stream, 0, 3, "num_extra_slice_header_bits");
  WRITE_U(stream, ENABLE_SIGN_HIDING, 1, "sign_data_hiding_flag");
  WRITE_U(stream, 0, 1, "cabac_init_present_flag");

  WRITE_UE(stream, 0, "num_ref_idx_l0_default_active_minus1");
  WRITE_UE(stream, 0, "num_ref_idx_l1_default_active_minus1");
  WRITE_SE(stream, ((int8_t)encoder_state->QP)-26, "pic_init_qp_minus26");
  WRITE_U(stream, 0, 1, "constrained_intra_pred_flag");
  WRITE_U(stream, encoder_state->encoder_control->trskip_enable, 1, "transform_skip_enabled_flag");
  WRITE_U(stream, 0, 1, "cu_qp_delta_enabled_flag");
  //if cu_qp_delta_enabled_flag
  //WRITE_UE(stream, 0, "diff_cu_qp_delta_depth");

  //TODO: add QP offsets
  WRITE_SE(stream, 0, "pps_cb_qp_offset");
  WRITE_SE(stream, 0, "pps_cr_qp_offset");
  WRITE_U(stream, 0, 1, "pps_slice_chroma_qp_offsets_present_flag");
  WRITE_U(stream, 0, 1, "weighted_pred_flag");
  WRITE_U(stream, 0, 1, "weighted_bipred_idc");

  //WRITE_U(stream, 0, 1, "dependent_slices_enabled_flag");
  WRITE_U(stream, 0, 1, "transquant_bypass_enable_flag");
  WRITE_U(stream, 0, 1, "tiles_enabled_flag");
  WRITE_U(stream, 0, 1, "entropy_coding_sync_enabled_flag");
  //TODO: enable tiles for concurrency
  //IF tiles
  //ENDIF
  WRITE_U(stream, 0, 1, "loop_filter_across_slice_flag");
  WRITE_U(stream, 1, 1, "deblocking_filter_control_present_flag");
  //IF deblocking_filter
    WRITE_U(stream, 0, 1, "deblocking_filter_override_enabled_flag");
  WRITE_U(stream, encoder_state->encoder_control->deblock_enable ? 0 : 1, 1,
          "pps_disable_deblocking_filter_flag");

    //IF !disabled
  if (encoder_state->encoder_control->deblock_enable) {
     WRITE_SE(stream, encoder_state->encoder_control->beta_offset_div2, "beta_offset_div2");
     WRITE_SE(stream, encoder_state->encoder_control->tc_offset_div2, "tc_offset_div2");
    }

    //ENDIF
  //ENDIF
  WRITE_U(stream, 0, 1, "pps_scaling_list_data_present_flag");
  //IF scaling_list
  //ENDIF
  WRITE_U(stream, 0, 1, "lists_modification_present_flag");
  WRITE_UE(stream, 0, "log2_parallel_merge_level_minus2");
  WRITE_U(stream, 0, 1, "slice_segment_header_extension_present_flag");
  WRITE_U(stream, 0, 1, "pps_extension_flag");
}

static void encode_PTL(encoder_state * const encoder_state)
{
  bitstream * const stream = &encoder_state->stream;
  int i;
  // PTL
  // Profile Tier
  WRITE_U(stream, 0, 2, "general_profile_space");
  WRITE_U(stream, 0, 1, "general_tier_flag");
  // Main Profile == 1
  WRITE_U(stream, 1, 5, "general_profile_idc");
  /* Compatibility flags should be set at general_profile_idc
   *  (so with general_profile_idc = 1, compatibility_flag[1] should be 1)
   * According to specification, when compatibility_flag[1] is set,
   *  compatibility_flag[2] should be set too.
   */
  WRITE_U(stream, 3<<29, 32, "general_profile_compatibility_flag[]");

  WRITE_U(stream, 1, 1, "general_progressive_source_flag");
  WRITE_U(stream, 0, 1, "general_interlaced_source_flag");
  WRITE_U(stream, 0, 1, "general_non_packed_constraint_flag");
  WRITE_U(stream, 0, 1, "general_frame_only_constraint_flag");

  WRITE_U(stream, 0, 32, "XXX_reserved_zero_44bits[0..31]");
  WRITE_U(stream, 0, 12, "XXX_reserved_zero_44bits[32..43]");

  // end Profile Tier

  // Level 6.2 (general_level_idc is 30 * 6.2)
  WRITE_U(stream, 186, 8, "general_level_idc");

  WRITE_U(stream, 0, 1, "sub_layer_profile_present_flag");
  WRITE_U(stream, 0, 1, "sub_layer_level_present_flag");

  for (i = 1; i < 8; i++) {
    WRITE_U(stream, 0, 2, "reserved_zero_2bits");
  }

  // end PTL
}

static void encode_scaling_list(encoder_state * const encoder_state)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  bitstream * const stream = &encoder_state->stream;
  uint32_t size_id;
  for (size_id = 0; size_id < SCALING_LIST_SIZE_NUM; size_id++) {
    int32_t list_id;
    for (list_id = 0; list_id < g_scaling_list_num[size_id]; list_id++) {
      uint8_t scaling_list_pred_mode_flag = 1;
      int32_t pred_list_idx;
      int32_t i;
      uint32_t ref_matrix_id = UINT32_MAX;

      for (pred_list_idx = list_id; pred_list_idx >= 0; pred_list_idx--) {
        const int32_t * const pred_list  = (list_id == pred_list_idx) ?
                                     scalinglist_get_default(size_id, pred_list_idx) :
                                     encoder->scaling_list.scaling_list_coeff[size_id][pred_list_idx];

        if (!memcmp(encoder->scaling_list.scaling_list_coeff[size_id][list_id], pred_list, sizeof(int32_t) * MIN(8, g_scaling_list_size[size_id])) &&
            ((size_id < SCALING_LIST_16x16) ||
             (encoder->scaling_list.scaling_list_dc[size_id][list_id] == encoder->scaling_list.scaling_list_dc[size_id][pred_list_idx]))) {
          ref_matrix_id = pred_list_idx;
          scaling_list_pred_mode_flag = 0;
          break;
        }
      }
      WRITE_U(stream, scaling_list_pred_mode_flag, 1, "scaling_list_pred_mode_flag" );

      if (!scaling_list_pred_mode_flag) {
        WRITE_UE(stream, list_id - ref_matrix_id, "scaling_list_pred_matrix_id_delta");
      } else {
        int32_t delta;
        const int32_t coef_num = MIN(MAX_MATRIX_COEF_NUM, g_scaling_list_size[size_id]);
        const uint32_t * const scan_cg = (size_id == 0) ? g_sig_last_scan_16x16 : g_sig_last_scan_32x32;
        int32_t next_coef = 8;
        const int32_t * const coef_list = encoder->scaling_list.scaling_list_coeff[size_id][list_id];

        if (size_id >= SCALING_LIST_16x16) {
          WRITE_SE(stream, encoder->scaling_list.scaling_list_dc[size_id][list_id] - 8, "scaling_list_dc_coef_minus8");
          next_coef = encoder->scaling_list.scaling_list_dc[size_id][list_id];
        }

        for (i = 0; i < coef_num; i++) {
          delta     = coef_list[scan_cg[i]] - next_coef;
          next_coef = coef_list[scan_cg[i]];
          if (delta > 127)
            delta -= 256;
          if (delta < -128)
            delta += 256;

          WRITE_SE(stream, delta, "scaling_list_delta_coef");
        }
      }
    }
  }
}

void encode_seq_parameter_set(encoder_state * const encoder_state)
{
  bitstream * const stream = &encoder_state->stream;
  const picture * const cur_pic = encoder_state->cur_pic;

#ifdef _DEBUG
  printf("=========== Sequence Parameter Set ID: 0 ===========\n");
#endif

  // TODO: profile IDC and level IDC should be defined later on
  WRITE_U(stream, 0, 4, "sps_video_parameter_set_id");
  WRITE_U(stream, 1, 3, "sps_max_sub_layers_minus1");
  WRITE_U(stream, 0, 1, "sps_temporal_id_nesting_flag");

  encode_PTL(encoder_state);

  WRITE_UE(stream, 0, "sps_seq_parameter_set_id");
  WRITE_UE(stream, encoder_state->encoder_control->in.video_format,
           "chroma_format_idc");

  if (encoder_state->encoder_control->in.video_format == 3) {
    WRITE_U(stream, 0, 1, "separate_colour_plane_flag");
  }

  WRITE_UE(stream, cur_pic->width, "pic_width_in_luma_samples");
  WRITE_UE(stream, cur_pic->height, "pic_height_in_luma_samples");

  if (cur_pic->width != encoder_state->encoder_control->in.real_width || cur_pic->height != encoder_state->encoder_control->in.real_height) {
    // The standard does not seem to allow setting conf_win values such that
    // the number of luma samples is not a multiple of 2. Options are to either
    // hide one line or show an extra line of non-video. Neither seems like a
    // very good option, so let's not even try.
    assert(!(cur_pic->width % 2));
    WRITE_U(stream, 1, 1, "conformance_window_flag");
    WRITE_UE(stream, 0, "conf_win_left_offset");
    WRITE_UE(stream, (cur_pic->width - encoder_state->encoder_control->in.real_width) >> 1,
             "conf_win_right_offset");
    WRITE_UE(stream, 0, "conf_win_top_offset");
    WRITE_UE(stream, (cur_pic->height - encoder_state->encoder_control->in.real_height) >> 1,
             "conf_win_bottom_offset");
  } else {
    WRITE_U(stream, 0, 1, "conformance_window_flag");
  }

  //IF window flag
  //END IF

  WRITE_UE(stream, encoder_state->encoder_control->bitdepth-8, "bit_depth_luma_minus8");
  WRITE_UE(stream, encoder_state->encoder_control->bitdepth-8, "bit_depth_chroma_minus8");
  WRITE_UE(stream, 0, "log2_max_pic_order_cnt_lsb_minus4");
  WRITE_U(stream, 0, 1, "sps_sub_layer_ordering_info_present_flag");

  //for each layer
  WRITE_UE(stream, 0, "sps_max_dec_pic_buffering");
  WRITE_UE(stream, 0, "sps_num_reorder_pics");
  WRITE_UE(stream, 0, "sps_max_latency_increase");
  //end for

  WRITE_UE(stream, MIN_SIZE-3, "log2_min_coding_block_size_minus3");
  WRITE_UE(stream, MAX_DEPTH, "log2_diff_max_min_coding_block_size");
  WRITE_UE(stream, 0, "log2_min_transform_block_size_minus2");   // 4x4
  WRITE_UE(stream, 3, "log2_diff_max_min_transform_block_size"); // 4x4...32x32
  WRITE_UE(stream, TR_DEPTH_INTER, "max_transform_hierarchy_depth_inter");
  WRITE_UE(stream, TR_DEPTH_INTRA, "max_transform_hierarchy_depth_intra");

  // scaling list
  WRITE_U(stream, encoder_state->encoder_control->scaling_list.enable, 1, "scaling_list_enable_flag");
  if (encoder_state->encoder_control->scaling_list.enable) {
    WRITE_U(stream, 1, 1, "sps_scaling_list_data_present_flag");
    encode_scaling_list(encoder_state);
  }

  WRITE_U(stream, 0, 1, "amp_enabled_flag");
  WRITE_U(stream, encoder_state->encoder_control->sao_enable ? 1 : 0, 1,
          "sample_adaptive_offset_enabled_flag");
  WRITE_U(stream, ENABLE_PCM, 1, "pcm_enabled_flag");
  #if ENABLE_PCM == 1
    WRITE_U(stream, 7, 4, "pcm_sample_bit_depth_luma_minus1");
    WRITE_U(stream, 7, 4, "pcm_sample_bit_depth_chroma_minus1");
    WRITE_UE(stream, 0, "log2_min_pcm_coding_block_size_minus3");
    WRITE_UE(stream, 2, "log2_diff_max_min_pcm_coding_block_size");
    WRITE_U(stream, 1, 1, "pcm_loop_filter_disable_flag");
  #endif

  WRITE_UE(stream, 0, "num_short_term_ref_pic_sets");

  //IF num short term ref pic sets
  //ENDIF

  WRITE_U(stream, 0, 1, "long_term_ref_pics_present_flag");

  //IF long_term_ref_pics_present
  //ENDIF

  WRITE_U(stream, ENABLE_TEMPORAL_MVP, 1,
          "sps_temporal_mvp_enable_flag");
  WRITE_U(stream, 0, 1, "sps_strong_intra_smoothing_enable_flag");
  WRITE_U(stream, 1, 1, "vui_parameters_present_flag");

  encode_VUI(encoder_state);

  WRITE_U(stream, 0, 1, "sps_extension_flag");
}

void encode_vid_parameter_set(encoder_state * const encoder_state)
{
  bitstream * const stream = &encoder_state->stream;
  int i;
#ifdef _DEBUG
  printf("=========== Video Parameter Set ID: 0 ===========\n");
#endif

  WRITE_U(stream, 0, 4, "vps_video_parameter_set_id");
  WRITE_U(stream, 3, 2, "vps_reserved_three_2bits" );
  WRITE_U(stream, 0, 6, "vps_reserved_zero_6bits" );
  WRITE_U(stream, 1, 3, "vps_max_sub_layers_minus1");
  WRITE_U(stream, 0, 1, "vps_temporal_id_nesting_flag");
  WRITE_U(stream, 0xffff, 16, "vps_reserved_ffff_16bits");

  encode_PTL(encoder_state);

  WRITE_U(stream, 0, 1, "vps_sub_layer_ordering_info_present_flag");

  //for each layer
  for (i = 0; i < 1; i++) {
  WRITE_UE(stream, 1, "vps_max_dec_pic_buffering");
  WRITE_UE(stream, 0, "vps_num_reorder_pics");
  WRITE_UE(stream, 0, "vps_max_latency_increase");
  }

  WRITE_U(stream, 0, 6, "vps_max_nuh_reserved_zero_layer_id");
  WRITE_UE(stream, 0, "vps_max_op_sets_minus1");
  WRITE_U(stream, 0, 1, "vps_timing_info_present_flag");

  //IF timing info
  //END IF

  WRITE_U(stream, 0, 1, "vps_extension_flag");
}

static void encode_VUI(encoder_state * const encoder_state)
{
  bitstream * const stream = &encoder_state->stream;
  const encoder_control * const encoder = encoder_state->encoder_control;
#ifdef _DEBUG
  printf("=========== VUI Set ID: 0 ===========\n");
#endif
  if (encoder->vui.sar_width > 0 && encoder->vui.sar_height > 0) {
    int i;
    static const struct
    {
      uint8_t width;
      uint8_t height;
      uint8_t idc;
    } sar[] = {
      // aspect_ratio_idc = 0 -> unspecified
      {  1,  1, 1 }, { 12, 11, 2 }, { 10, 11, 3 }, { 16, 11, 4 },
      { 40, 33, 5 }, { 24, 11, 6 }, { 20, 11, 7 }, { 32, 11, 8 },
      { 80, 33, 9 }, { 18, 11, 10}, { 15, 11, 11}, { 64, 33, 12},
      {160, 99, 13}, {  4,  3, 14}, {  3,  2, 15}, {  2,  1, 16},
      // aspect_ratio_idc = [17..254] -> reserved
      { 0, 0, 255 }
    };

    for (i = 0; sar[i].idc != 255; i++)
      if (sar[i].width  == encoder->vui.sar_width &&
          sar[i].height == encoder->vui.sar_height)
        break;

    WRITE_U(stream, 1, 1, "aspect_ratio_info_present_flag");
    WRITE_U(stream, sar[i].idc, 8, "aspect_ratio_idc");
    if (sar[i].idc == 255) {
      // EXTENDED_SAR
      WRITE_U(stream, encoder->vui.sar_width, 16, "sar_width");
      WRITE_U(stream, encoder->vui.sar_height, 16, "sar_height");
    }
  } else
    WRITE_U(stream, 0, 1, "aspect_ratio_info_present_flag");

  //IF aspect ratio info
  //ENDIF

  if (encoder->vui.overscan > 0) {
    WRITE_U(stream, 1, 1, "overscan_info_present_flag");
    WRITE_U(stream, encoder->vui.overscan - 1, 1, "overscan_appropriate_flag");
  } else
    WRITE_U(stream, 0, 1, "overscan_info_present_flag");

  //IF overscan info
  //ENDIF

  if (encoder->vui.videoformat != 5 || encoder->vui.fullrange ||
      encoder->vui.colorprim != 2 || encoder->vui.transfer != 2 ||
      encoder->vui.colormatrix != 2) {
    WRITE_U(stream, 1, 1, "video_signal_type_present_flag");
    WRITE_U(stream, encoder->vui.videoformat, 3, "video_format");
    WRITE_U(stream, encoder->vui.fullrange, 1, "video_full_range_flag");

    if (encoder->vui.colorprim != 2 || encoder->vui.transfer != 2 ||
        encoder->vui.colormatrix != 2) {
      WRITE_U(stream, 1, 1, "colour_description_present_flag");
      WRITE_U(stream, encoder->vui.colorprim, 8, "colour_primaries");
      WRITE_U(stream, encoder->vui.transfer, 8, "transfer_characteristics");
      WRITE_U(stream, encoder->vui.colormatrix, 8, "matrix_coeffs");
    } else
      WRITE_U(stream, 0, 1, "colour_description_present_flag");
  } else
    WRITE_U(stream, 0, 1, "video_signal_type_present_flag");

  //IF video type
  //ENDIF

  if (encoder->vui.chroma_loc > 0) {
    WRITE_U(stream, 1, 1, "chroma_loc_info_present_flag");
    WRITE_UE(stream, encoder->vui.chroma_loc, "chroma_sample_loc_type_top_field");
    WRITE_UE(stream, encoder->vui.chroma_loc, "chroma_sample_loc_type_bottom_field");
  } else
    WRITE_U(stream, 0, 1, "chroma_loc_info_present_flag");

  //IF chroma loc info
  //ENDIF

  WRITE_U(stream, 0, 1, "neutral_chroma_indication_flag");
  WRITE_U(stream, 0, 1, "field_seq_flag");
  WRITE_U(stream, 0, 1, "frame_field_info_present_flag");
  WRITE_U(stream, 0, 1, "default_display_window_flag");

  //IF default display window
  //ENDIF

  WRITE_U(stream, 0, 1, "vui_timing_info_present_flag");

  //IF timing info
  //ENDIF

  WRITE_U(stream, 0, 1, "bitstream_restriction_flag");

  //IF bitstream restriction
  //ENDIF
}

void encode_slice_header(encoder_state * const encoder_state)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  bitstream * const stream = &encoder_state->stream;
  const picture * const cur_pic = encoder_state->cur_pic;

#ifdef _DEBUG
  printf("=========== Slice ===========\n");
#endif

  WRITE_U(stream, 1, 1, "first_slice_segment_in_pic_flag");

  if (cur_pic->type >= NAL_BLA_W_LP
      && cur_pic->type <= NAL_RSV_IRAP_VCL23) {
    WRITE_U(stream, 1, 1, "no_output_of_prior_pics_flag");
  }

  WRITE_UE(stream, 0, "slice_pic_parameter_set_id");

  //WRITE_U(stream, 0, 1, "dependent_slice_segment_flag");

  WRITE_UE(stream, cur_pic->slicetype, "slice_type");

  // if !entropy_slice_flag

    //if output_flag_present_flag
      //WRITE_U(stream, 1, 1, "pic_output_flag");
    //end if
    //if( IdrPicFlag ) <- nal_unit_type == 5
  if (cur_pic->type != NAL_IDR_W_RADL
      && cur_pic->type != NAL_IDR_N_LP) {
      int j;
      int ref_negative = encoder_state->ref->used_size;
      int ref_positive = 0;
      WRITE_U(stream, encoder_state->poc&0xf, 4, "pic_order_cnt_lsb");
      WRITE_U(stream, 0, 1, "short_term_ref_pic_set_sps_flag");
      WRITE_UE(stream, ref_negative, "num_negative_pics");
      WRITE_UE(stream, ref_positive, "num_positive_pics");

    for (j = 0; j < ref_negative; j++) {
      int32_t delta_poc_minus1 = 0;
      WRITE_UE(stream, delta_poc_minus1, "delta_poc_s0_minus1");
      WRITE_U(stream,1,1, "used_by_curr_pic_s0_flag");
    }

    //WRITE_UE(stream, 0, "short_term_ref_pic_set_idx");
  }

    //end if
  //end if
  if (encoder->sao_enable) {
    WRITE_U(stream, cur_pic->slice_sao_luma_flag, 1, "slice_sao_luma_flag");
    WRITE_U(stream, cur_pic->slice_sao_chroma_flag, 1, "slice_sao_chroma_flag");
  }

  if (cur_pic->slicetype != SLICE_I) {
      WRITE_U(stream, 1, 1, "num_ref_idx_active_override_flag");
        WRITE_UE(stream, encoder_state->ref->used_size-1, "num_ref_idx_l0_active_minus1");
      WRITE_UE(stream, 5-MRG_MAX_NUM_CANDS, "five_minus_max_num_merge_cand");
  }

  if (cur_pic->slicetype == SLICE_B) {
      WRITE_U(stream, 0, 1, "mvd_l1_zero_flag");
  }

  // Skip flags that are not present
  // if !entropy_slice_flag
    WRITE_SE(stream, 0, "slice_qp_delta");
    //WRITE_U(stream, 1, 1, "alignment");
}


static void encode_sao_color(const encoder_state * const encoder_state, cabac_data *cabac, sao_info *sao,
                             color_index color_i)
{
  const picture * const cur_pic = encoder_state->cur_pic;
  sao_eo_cat i;

  // Skip colors with no SAO.
  if (color_i == COLOR_Y && !cur_pic->slice_sao_luma_flag) return;
  if (color_i != COLOR_Y && !cur_pic->slice_sao_chroma_flag) return;

  /// sao_type_idx_luma:   TR, cMax = 2, cRiceParam = 0, bins = {0, bypass}
  /// sao_type_idx_chroma: TR, cMax = 2, cRiceParam = 0, bins = {0, bypass}
  // Encode sao_type_idx for Y and U+V.
  if (color_i != COLOR_V) {
    cabac->ctx = &(cabac->ctx_sao_type_idx_model);;
    CABAC_BIN(cabac, sao->type == SAO_TYPE_NONE ? 0 : 1, "sao_type_idx");
    if (sao->type == SAO_TYPE_BAND) {
      CABAC_BIN_EP(cabac, 0, "sao_type_idx_ep");
    } else if (sao->type == SAO_TYPE_EDGE) {
      CABAC_BIN_EP(cabac, 1, "sao_type_idx_ep");
    }
  }

  if (sao->type == SAO_TYPE_NONE) return;

  /// sao_offset_abs[][][][]: TR, cMax = (1 << (Min(bitDepth, 10) - 5)) - 1,
  ///                         cRiceParam = 0, bins = {bypass x N}
  for (i = SAO_EO_CAT1; i <= SAO_EO_CAT4; ++i) {
    cabac_write_unary_max_symbol_ep(cabac, abs(sao->offsets[i]), SAO_ABS_OFFSET_MAX);
  }

  /// sao_offset_sign[][][][]: FL, cMax = 1, bins = {bypass}
  /// sao_band_position[][][]: FL, cMax = 31, bins = {bypass x N}
  /// sao_eo_class_luma:       FL, cMax = 3, bins = {bypass x 3}
  /// sao_eo_class_chroma:     FL, cMax = 3, bins = {bypass x 3}
  if (sao->type == SAO_TYPE_BAND) {
    for (i = SAO_EO_CAT1; i <= SAO_EO_CAT4; ++i) {
      // Positive sign is coded as 0.
      if(sao->offsets[i] != 0) {
        CABAC_BIN_EP(cabac, sao->offsets[i] < 0 ? 1 : 0, "sao_offset_sign");
      }
    }
    // TODO: sao_band_position
    // FL cMax=31 (5 bits)
    CABAC_BINS_EP(cabac, sao->band_position, 5, "sao_band_position");
  } else if (color_i != COLOR_V) {
    CABAC_BINS_EP(cabac, sao->eo_class, 2, "sao_eo_class");
  }
}

static void encode_sao_merge_flags(sao_info *sao, cabac_data *cabac,
                                   unsigned x_ctb, unsigned y_ctb)
{
  // SAO merge flags are not present for the first row and column.
  if (x_ctb > 0) {
    cabac->ctx = &(cabac->ctx_sao_merge_flag_model);
    CABAC_BIN(cabac, sao->merge_left_flag ? 1 : 0, "sao_merge_left_flag");
  }
  if (y_ctb > 0 && !sao->merge_left_flag) {
    cabac->ctx = &(cabac->ctx_sao_merge_flag_model);
    CABAC_BIN(cabac, sao->merge_up_flag ? 1 : 0, "sao_merge_up_flag");
  }
}

/**
 * \brief Encode SAO information.
 */
static void encode_sao(encoder_state * const encoder_state,
                       cabac_data *cabac,
                       unsigned x_lcu, uint16_t y_lcu,
                       sao_info *sao_luma, sao_info *sao_chroma)
{
  // TODO: transmit merge flags outside sao_info
  encode_sao_merge_flags(sao_luma, cabac, x_lcu, y_lcu);

  // If SAO is merged, nothing else needs to be coded.
  if (!sao_luma->merge_left_flag && !sao_luma->merge_up_flag) {
    encode_sao_color(encoder_state, cabac, sao_luma, COLOR_Y);
    encode_sao_color(encoder_state, cabac, sao_chroma, COLOR_U);
    encode_sao_color(encoder_state, cabac, sao_chroma, COLOR_V);
  }
}


void encode_coding_tree(encoder_state * const encoder_state, cabac_data *cabac,
                        uint16_t x_ctb, uint16_t y_ctb, uint8_t depth)
{
  const picture * const cur_pic = encoder_state->cur_pic;
  cu_info *cur_cu = &cur_pic->cu_array[MAX_DEPTH][x_ctb + y_ctb * (cur_pic->width_in_lcu << MAX_DEPTH)];
  uint8_t split_flag = GET_SPLITDATA(cur_cu, depth);
  uint8_t split_model = 0;

  // Check for slice border
  uint8_t border_x = ((cur_pic->width) < (x_ctb * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth))) ? 1 : 0;
  uint8_t border_y = ((cur_pic->height) < (y_ctb * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth))) ? 1 : 0;
  uint8_t border_split_x = ((cur_pic->width)  < ((x_ctb + 1) * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> (depth + 1)))) ? 0 : 1;
  uint8_t border_split_y = ((cur_pic->height) < ((y_ctb + 1) * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> (depth + 1)))) ? 0 : 1;
  uint8_t border = border_x | border_y; /*!< are we in any border CU */


  // When not in MAX_DEPTH, insert split flag and split the blocks if needed
  if (depth != MAX_DEPTH) {
    // Implisit split flag when on border
    if (!border) {
      // Get left and top block split_flags and if they are present and true, increase model number
      if (x_ctb > 0 && GET_SPLITDATA(&(cur_pic->cu_array[MAX_DEPTH][x_ctb - 1 + y_ctb * (cur_pic->width_in_lcu << MAX_DEPTH)]), depth) == 1) {
        split_model++;
      }

      if (y_ctb > 0 && GET_SPLITDATA(&(cur_pic->cu_array[MAX_DEPTH][x_ctb + (y_ctb - 1) * (cur_pic->width_in_lcu << MAX_DEPTH)]), depth) == 1) {
        split_model++;
      }

      cabac->ctx = &(cabac->ctx_split_flag_model[split_model]);
      CABAC_BIN(cabac, split_flag, "SplitFlag");
    }

    if (split_flag || border) {
      // Split blocks and remember to change x and y block positions
      uint8_t change = 1<<(MAX_DEPTH-1-depth);
      encode_coding_tree(encoder_state, cabac, x_ctb, y_ctb, depth + 1); // x,y

      // TODO: fix when other half of the block would not be completely over the border
      if (!border_x || border_split_x) {
        encode_coding_tree(encoder_state, cabac, x_ctb + change, y_ctb, depth + 1);
      }
      if (!border_y || border_split_y) {
        encode_coding_tree(encoder_state, cabac, x_ctb, y_ctb + change, depth + 1);
      }
      if (!border || (border_split_x && border_split_y)) {
        encode_coding_tree(encoder_state, cabac, x_ctb + change, y_ctb + change, depth + 1);
      }
      return;
    }
  }



    // Encode skip flag
  if (cur_pic->slicetype != SLICE_I) {
    int8_t ctx_skip = 0; // uiCtxSkip = aboveskipped + leftskipped;
    int ui;
    int16_t num_cand = MRG_MAX_NUM_CANDS;
    // Get left and top skipped flags and if they are present and true, increase context number
    if (x_ctb > 0 && (&cur_pic->cu_array[MAX_DEPTH][x_ctb - 1 + y_ctb * (cur_pic->width_in_lcu << MAX_DEPTH)])->skipped) {
      ctx_skip++;
    }

    if (y_ctb > 0 && (&cur_pic->cu_array[MAX_DEPTH][x_ctb + (y_ctb - 1) * (cur_pic->width_in_lcu << MAX_DEPTH)])->skipped) {
      ctx_skip++;
    }

    cabac->ctx = &(cabac->ctx_cu_skip_flag_model[ctx_skip]);
    CABAC_BIN(cabac, cur_cu->skipped, "SkipFlag");

    // IF SKIP
    if (cur_cu->skipped) {
      if (num_cand > 1) {
        for (ui = 0; ui < num_cand - 1; ui++) {
          int32_t symbol = (ui != cur_cu->merge_idx);
          if (ui == 0) {
            cabac->ctx = &(cabac->ctx_cu_merge_idx_ext_model);
            CABAC_BIN(cabac, symbol, "MergeIndex");
          } else {
            CABAC_BIN_EP(cabac,symbol,"MergeIndex");
          }
          if (symbol == 0) {
            break;
          }
        }
      }
      return;
    }
  }

  // ENDIF SKIP

  // Prediction mode
  if (cur_pic->slicetype != SLICE_I) {
    cabac->ctx = &(cabac->ctx_cu_pred_mode_model);
    CABAC_BIN(cabac, (cur_cu->type == CU_INTRA), "PredMode");
  }

  // part_mode
  if (cur_cu->type == CU_INTRA) {
    if (depth == MAX_DEPTH) {
      cabac->ctx = &(cabac->ctx_part_size_model[0]);
      if (cur_cu->part_size == SIZE_2Nx2N) {
        CABAC_BIN(cabac, 1, "part_mode 2Nx2N");
      } else {
        CABAC_BIN(cabac, 0, "part_mode NxN");
      }
    }
  } else {
    // TODO: Handle inter sizes other than 2Nx2N
    cabac->ctx = &(cabac->ctx_part_size_model[0]);
    CABAC_BIN(cabac, 1, "part_mode 2Nx2N");
  }

  //end partsize
  if (cur_cu->type == CU_INTER) {
    // FOR each part
    // Mergeflag
    int16_t num_cand = 0;
    cabac->ctx = &(cabac->ctx_cu_merge_flag_ext_model);
    CABAC_BIN(cabac, cur_cu->merged, "MergeFlag");
    num_cand = MRG_MAX_NUM_CANDS;
    if (cur_cu->merged) { //merge
      if (num_cand > 1) {
        int32_t ui;
        for (ui = 0; ui < num_cand - 1; ui++) {
          int32_t symbol = (ui != cur_cu->merge_idx);
          if (ui == 0) {
            cabac->ctx = &(cabac->ctx_cu_merge_idx_ext_model);
            CABAC_BIN(cabac, symbol, "MergeIndex");
          } else {
            CABAC_BIN_EP(cabac,symbol,"MergeIndex");
          }
          if (symbol == 0) break;
        }
      }
    } else {
      uint32_t ref_list_idx;
      /*
      // Void TEncSbac::codeInterDir( TComDataCU* pcCU, UInt uiAbsPartIdx )
      if(cur_pic->slicetype == SLICE_B)
      {
        // Code Inter Dir
        const UInt uiInterDir = pcCU->getInterDir( uiAbsPartIdx ) - 1;
        const UInt uiCtx      = pcCU->getCtxInterDir( uiAbsPartIdx );
        ContextModel *pCtx    = m_cCUInterDirSCModel.get( 0 );
        if (pcCU->getPartitionSize(uiAbsPartIdx) == SIZE_2Nx2N || pcCU->getHeight(uiAbsPartIdx) != 8 )
        {
          m_pcBinIf->encodeBin( uiInterDir == 2 ? 1 : 0, *( pCtx + uiCtx ) );
        }
        if (uiInterDir < 2)
        {
          m_pcBinIf->encodeBin( uiInterDir, *( pCtx + 4 ) );
        }
      }
      */

      for (ref_list_idx = 0; ref_list_idx < 2; ref_list_idx++) {
            //if(encoder_state->ref_idx_num[uiRefListIdx] > 0)
            {
          if (cur_cu->inter.mv_dir & (1 << ref_list_idx)) {
            if (encoder_state->ref->used_size != 1) { //encoder_state->ref_idx_num[uiRefListIdx] != 1)//NumRefIdx != 1)
              // parseRefFrmIdx
              int32_t ref_frame = cur_cu->inter.mv_ref;

              cabac->ctx = &(cabac->ctx_cu_ref_pic_model[0]);
              CABAC_BIN(cabac, (ref_frame == 0) ? 0 : 1, "ref_frame_flag");

              if (ref_frame > 0) {
                int32_t i;
                int32_t ref_num = encoder_state->ref->used_size - 2;

                cabac->ctx = &(cabac->ctx_cu_ref_pic_model[1]);
                ref_frame--;

                for (i = 0; i < ref_num; ++i) {
                  const uint32_t symbol = (i == ref_frame) ? 0 : 1;

                  if (i == 0) {
                    CABAC_BIN(cabac, symbol, "ref_frame_flag2");
                  } else {
                    CABAC_BIN_EP(cabac, symbol, "ref_frame_flag2");
                  }
                  if (symbol == 0) break;
                }
              }
            }

            if (!(/*pcCU->getSlice()->getMvdL1ZeroFlag() &&*/ encoder_state->ref_list == REF_PIC_LIST_1 && cur_cu->inter.mv_dir == 3)) {
              const int32_t mvd_hor = cur_cu->inter.mvd[0];
              const int32_t mvd_ver = cur_cu->inter.mvd[1];
              const int8_t hor_abs_gr0 = mvd_hor != 0;
              const int8_t ver_abs_gr0 = mvd_ver != 0;
              const uint32_t mvd_hor_abs = abs(mvd_hor);
              const uint32_t mvd_ver_abs = abs(mvd_ver);

              cabac->ctx = &(cabac->ctx_cu_mvd_model[0]);
              CABAC_BIN(cabac, (mvd_hor!=0)?1:0, "abs_mvd_greater0_flag_hor");
              CABAC_BIN(cabac, (mvd_ver!=0)?1:0, "abs_mvd_greater0_flag_ver");

              cabac->ctx = &(cabac->ctx_cu_mvd_model[1]);

              if (hor_abs_gr0) {
                CABAC_BIN(cabac, (mvd_hor_abs>1)?1:0, "abs_mvd_greater1_flag_hor");
              }

              if (ver_abs_gr0) {
                CABAC_BIN(cabac, (mvd_ver_abs>1)?1:0, "abs_mvd_greater1_flag_ver");
              }

              if (hor_abs_gr0) {
                if (mvd_hor_abs > 1) {
                  cabac_write_ep_ex_golomb(cabac,mvd_hor_abs-2, 1);
                }

                CABAC_BIN_EP(cabac, (mvd_hor>0)?0:1, "mvd_sign_flag_hor");
              }

              if (ver_abs_gr0) {
                if (mvd_ver_abs > 1) {
                  cabac_write_ep_ex_golomb(cabac,mvd_ver_abs-2, 1);
                }

                CABAC_BIN_EP(cabac, (mvd_ver>0)?0:1, "mvd_sign_flag_ver");
              }
            }

            // Signal which candidate MV to use
            cabac_write_unary_max_symbol(cabac, cabac->ctx_mvp_idx_model, cur_cu->inter.mv_cand, 1,
                                        AMVP_MAX_NUM_CANDS - 1);
          }
          }
        } // for ref_list
    } // if !merge


    // Only need to signal coded block flag if not skipped or merged
    // skip = no coded residual, merge = coded residual
    if (!cur_cu->merged) {
      cabac->ctx = &(cabac->ctx_cu_qt_root_cbf_model);
      CABAC_BIN(cabac, cur_cu->coeff_top_y[depth] | cur_cu->coeff_top_u[depth] | cur_cu->coeff_top_v[depth], "rqt_root_cbf");
    }
    // Code (possible) coeffs to bitstream

    if(cur_cu->coeff_top_y[depth] | cur_cu->coeff_top_u[depth] | cur_cu->coeff_top_v[depth]) {
      encode_transform_coeff(encoder_state, cabac, x_ctb * 2, y_ctb * 2, depth, 0, 0, 0);
    }


    // END for each part
  } else if (cur_cu->type == CU_INTRA) {
    uint8_t intra_pred_mode[4] = {
      cur_cu->intra[0].mode, cur_cu->intra[1].mode,
      cur_cu->intra[2].mode, cur_cu->intra[3].mode };
    uint8_t intra_pred_mode_chroma = 36; // 36 = Chroma derived from luma
    int8_t intra_preds[4][3] = {{-1, -1, -1},{-1, -1, -1},{-1, -1, -1},{-1, -1, -1}};
    int8_t mpm_preds[4] = {-1, -1, -1, -1};
    int i, j;
    uint32_t flag[4];
    int num_pred_units = (cur_cu->part_size == SIZE_2Nx2N ? 1 : 4);

    #if ENABLE_PCM == 1
    // Code must start after variable initialization
    cabac_encode_bin_trm(cabac, 0); // IPCMFlag == 0
    #endif

    // PREDINFO CODING
    // If intra prediction mode is found from the predictors,
    // it can be signaled with two EP's. Otherwise we can send
    // 5 EP bins with the full predmode
    for (j = 0; j < num_pred_units; ++j) {
      static const vector2d offset[4] = {{0,0},{1,0},{0,1},{1,1}};
      cu_info *left_cu = 0;
      cu_info *above_cu = 0;

      if (x_ctb > 0) {
        left_cu = &cur_pic->cu_array[MAX_DEPTH][x_ctb - 1 + y_ctb * (cur_pic->width_in_lcu << MAX_DEPTH)];
      }
      // Don't take the above CU across the LCU boundary.
      if (y_ctb > 0 && (y_ctb & 7) != 0) {
        above_cu = &cur_pic->cu_array[MAX_DEPTH][x_ctb + (y_ctb - 1) * (cur_pic->width_in_lcu << MAX_DEPTH)];
      }

      intra_get_dir_luma_predictor((x_ctb<<3) + (offset[j].x<<2),
                                   (y_ctb<<3) + (offset[j].y<<2),
                                   intra_preds[j], cur_cu,
                                   left_cu, above_cu);
      for (i = 0; i < 3; i++) {
        if (intra_preds[j][i] == intra_pred_mode[j]) {
          mpm_preds[j] = (int8_t)i;
          break;
        }
      }
      flag[j] = (mpm_preds[j] == -1) ? 0 : 1;
    }

    cabac->ctx = &(cabac->ctx_intra_mode_model);
    for (j = 0; j < num_pred_units; ++j) {
      CABAC_BIN(cabac, flag[j], "prev_intra_luma_pred_flag");
    }

    for (j = 0; j < num_pred_units; ++j) {
      // Signal index of the prediction mode in the prediction list.
      if (flag[j]) {
        CABAC_BIN_EP(cabac, (mpm_preds[j] == 0 ? 0 : 1), "mpm_idx");
        if (mpm_preds[j] != 0) {
          CABAC_BIN_EP(cabac, (mpm_preds[j] == 1 ? 0 : 1), "mpm_idx");
        }
      } else {
        // Signal the actual prediction mode.
        int32_t tmp_pred = intra_pred_mode[j];

        // Sort prediction list from lowest to highest.
        if (intra_preds[j][0] > intra_preds[j][1]) SWAP(intra_preds[j][0], intra_preds[j][1], int8_t);
        if (intra_preds[j][0] > intra_preds[j][2]) SWAP(intra_preds[j][0], intra_preds[j][2], int8_t);
        if (intra_preds[j][1] > intra_preds[j][2]) SWAP(intra_preds[j][1], intra_preds[j][2], int8_t);

        // Reduce the index of the signaled prediction mode according to the
        // prediction list, as it has been already signaled that it's not one
        // of the prediction modes.
        for (i = 2; i >= 0; i--) {
          tmp_pred = (tmp_pred > intra_preds[j][i] ? tmp_pred - 1 : tmp_pred);
        }

        CABAC_BINS_EP(cabac, tmp_pred, 5, "rem_intra_luma_pred_mode");
      }
    }

    {  // start intra chroma pred mode coding
      unsigned pred_mode = 5;
      unsigned chroma_pred_modes[4] = {0, 26, 10, 1};

      if (intra_pred_mode_chroma == 36) {
        pred_mode = 4;
      } else if (intra_pred_mode_chroma == 34) {
        // Angular 34 mode is possible only if intra pred mode is one of the
        // possible chroma pred modes, in which case it is signaled with that
        // duplicate mode.
        for (i = 0; i < 4; ++i) {
          if (intra_pred_mode[0] == chroma_pred_modes[i]) pred_mode = i;
        }
      } else {
        for (i = 0; i < 4; ++i) {
          if (intra_pred_mode_chroma == chroma_pred_modes[i]) pred_mode = i;
        }
      }

      /**
       * Table 9-35 - Binarization for intra_chroma_pred_mode
       *   intra_chroma_pred_mode  bin_string
       *                        4           0
       *                        0         100
       *                        1         101
       *                        2         110
       *                        3         111
       * Table 9-37 - Assignment of ctxInc to syntax elements with context coded bins
       *   intra_chroma_pred_mode[][] = 0, bypass, bypass
       */
      cabac->ctx = &(cabac->ctx_chroma_pred_model[0]);
      if (pred_mode == 4) {
        CABAC_BIN(cabac, 0, "intra_chroma_pred_mode");
      } else {
        CABAC_BIN(cabac, 1, "intra_chroma_pred_mode");
        CABAC_BINS_EP(cabac, pred_mode, 2, "intra_chroma_pred_mode");
      }
    }  // end intra chroma pred mode coding

    encode_transform_coeff(encoder_state, cabac, x_ctb * 2, y_ctb * 2, depth, 0, 0, 0);
  }

    #if ENABLE_PCM == 1
  // Code IPCM block
  if (cur_cu->type == CU_PCM) {
    cabac_encode_bin_trm(cabac, 1); // IPCMFlag == 1
      cabac_finish(cabac);
      bitstream_align(cabac.stream);
    // PCM sample
      {
      unsigned y, x;

      pixel *base_y = &cur_pic->y_data[x_ctb * (LCU_WIDTH >> (MAX_DEPTH))    + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH))) * encoder->in.width];
      pixel *base_u = &cur_pic->u_data[(x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1))) * encoder->in.width / 2)];
      pixel *base_v = &cur_pic->v_data[(x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1))) * encoder->in.width / 2)];

      // Luma
      for (y = 0; y < LCU_WIDTH >> depth; y++) {
        for (x = 0; x < LCU_WIDTH >> depth; x++) {
          bitstream_put(cabac.stream, base_y[x + y * encoder->in.width], 8);
          }
        }

      // Chroma
      if (encoder->in.video_format != FORMAT_400) {
        for (y = 0; y < LCU_WIDTH >> (depth + 1); y++) {
          for (x = 0; x < LCU_WIDTH >> (depth + 1); x++) {
            bitstream_put(cabac.stream, base_u[x + y * (encoder->in.width >> 1)], 8);
          }
        }
        for (y = 0; y < LCU_WIDTH >> (depth + 1); y++) {
          for (x = 0; x < LCU_WIDTH >> (depth + 1); x++) {
            bitstream_put(cabac.stream, base_v[x + y * (encoder->in.width >> 1)], 8);
          }
        }
      }
    }
    // end PCM sample
      cabac_start(cabac);
  } // end Code IPCM block
#endif /* END ENABLE_PCM */
  else { /* Should not happend */
    printf("UNHANDLED TYPE!\r\n");
    exit(1);
  }

   /* end prediction unit */
  /* end coding_unit */
}

static void transform_chroma(encoder_state * const encoder_state, cabac_data *cabac, cu_info *cur_cu,
                             int depth, pixel *base_u, pixel *pred_u,
                             coefficient *coeff_u, int8_t scan_idx_chroma,
                             coefficient *pre_quant_coeff, coefficient *block)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  int base_stride = LCU_WIDTH;
  int pred_stride = LCU_WIDTH;

  int8_t width_c = LCU_WIDTH >> (depth + 1);

  int i = 0;
  unsigned ac_sum = 0;

  int y, x;

  for (y = 0; y < width_c; y++) {
    for (x = 0; x < width_c; x++) {
      block[i] = ((int16_t)base_u[x + y * (base_stride >> 1)]) -
                  pred_u[x + y * (pred_stride >> 1)];
      i++;
    }
  }

  transform2d(encoder, block, pre_quant_coeff, width_c, 65535);
  if (encoder->rdoq_enable) {
    rdoq(encoder_state, cabac, pre_quant_coeff, coeff_u, width_c, width_c, &ac_sum, 2,
         scan_idx_chroma, cur_cu->type, cur_cu->tr_depth-cur_cu->depth);
  } else {
    quant(encoder_state, pre_quant_coeff, coeff_u, width_c, width_c, &ac_sum, 2,
          scan_idx_chroma, cur_cu->type);
  }
}

static void reconstruct_chroma(const encoder_state * const encoder_state, cu_info *cur_cu,
                               int depth, int has_coeffs, coefficient *coeff_u,
                               pixel *recbase_u, pixel *pred_u, int color_type,
                               coefficient *pre_quant_coeff, coefficient *block)
{
  int8_t width_c = LCU_WIDTH >> (depth + 1);
  const int pred_stride = LCU_WIDTH;
  const int recbase_stride = LCU_WIDTH;

  int i, y, x;

  if (has_coeffs) {
    // RECONSTRUCT for predictions
    dequant(encoder_state, coeff_u, pre_quant_coeff, width_c, width_c, (int8_t)color_type, cur_cu->type);
    itransform2d(encoder_state->encoder_control, block, pre_quant_coeff, width_c, 65535);

    i = 0;

    for (y = 0; y < width_c; y++) {
      for (x = 0; x < width_c; x++) {
        int16_t val = block[i++] + pred_u[x + y * (pred_stride >> 1)];
        //TODO: support 10+bits
        recbase_u[x + y * (recbase_stride >> 1)] = (uint8_t)CLIP(0, 255, val);
      }
    }

    // END RECONTRUCTION
  } else {
    // without coeffs, we only use the prediction
    for (y = 0; y < width_c; y++) {
      for (x = 0; x < width_c; x++) {
        recbase_u[x + y * (recbase_stride >> 1)] = (uint8_t)CLIP(0, 255, pred_u[x + y * (pred_stride >> 1)]);
      }
    }
  }
}

void encode_transform_tree(encoder_state * const encoder_state, cabac_data* cabac, int32_t x, int32_t y, const uint8_t depth, lcu_t* lcu)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  // we have 64>>depth transform size
  int x_local = (x&0x3f), y_local = (y&0x3f);
  cu_info *cur_cu = &lcu->cu[LCU_CU_OFFSET + (x_local>>3) + (y_local>>3)*LCU_T_CU_WIDTH];

  int i;
  const int8_t width = LCU_WIDTH>>depth;
  const int8_t width_c = (depth == MAX_DEPTH + 1 ? width : width  / 2);

  // Tell clang-analyzer what is up. For some reason it can't figure out from
  // asserting just depth.
  assert(width == 4 || width == 8 || width == 16 || width == 32 || width == 64);

  // Split transform and increase depth
  if (depth == 0 || cur_cu->tr_depth > depth) {
    int offset = width_c;
    encode_transform_tree(encoder_state, cabac, x,          y,          depth+1, lcu);
    encode_transform_tree(encoder_state, cabac, x + offset, y,          depth+1, lcu);
    encode_transform_tree(encoder_state, cabac, x,          y + offset, depth+1, lcu);
    encode_transform_tree(encoder_state, cabac, x + offset, y + offset, depth+1, lcu);

    // Derive coded coeff flags from the next depth
    if (depth == MAX_DEPTH) {
      cur_cu->coeff_top_y[depth] = cur_cu->coeff_top_y[depth+1] | cur_cu->coeff_top_y[depth+2] | cur_cu->coeff_top_y[depth+3] | cur_cu->coeff_top_y[depth+4];
      cur_cu->coeff_top_u[depth] = cur_cu->coeff_top_u[depth+1];
      cur_cu->coeff_top_v[depth] = cur_cu->coeff_top_v[depth+1];
    } else {
      cu_info *cu_a =  &lcu->cu[LCU_CU_OFFSET + ((x_local + offset)>>3) +  (y_local>>3)        *LCU_T_CU_WIDTH];
      cu_info *cu_b =  &lcu->cu[LCU_CU_OFFSET +  (x_local>>3)           + ((y_local+offset)>>3)*LCU_T_CU_WIDTH];
      cu_info *cu_c =  &lcu->cu[LCU_CU_OFFSET + ((x_local + offset)>>3) + ((y_local+offset)>>3)*LCU_T_CU_WIDTH];
      cur_cu->coeff_top_y[depth] = cur_cu->coeff_top_y[depth+1] | cu_a->coeff_top_y[depth+1] | cu_b->coeff_top_y[depth+1]
                                    | cu_c->coeff_top_y[depth+1];
      cur_cu->coeff_top_u[depth] = cur_cu->coeff_top_u[depth+1] | cu_a->coeff_top_u[depth+1] | cu_b->coeff_top_u[depth+1]
                                    | cu_c->coeff_top_u[depth+1];
      cur_cu->coeff_top_v[depth] = cur_cu->coeff_top_v[depth+1] | cu_a->coeff_top_v[depth+1] | cu_b->coeff_top_v[depth+1]
                                    | cu_c->coeff_top_v[depth+1];
    }

    return;
  }

  {
    // INTRAPREDICTION VARIABLES
      // Pointers to reconstruction arrays
    pixel *recbase_y = &lcu->rec.y[x_local + y_local * LCU_WIDTH];
    pixel *recbase_u = &lcu->rec.u[x_local/2 + (y_local * LCU_WIDTH)/4];
    pixel *recbase_v = &lcu->rec.v[x_local/2 + (y_local * LCU_WIDTH)/4];
    int32_t recbase_stride = LCU_WIDTH;


    pixel *base_y = &lcu->ref.y[x_local + y_local * LCU_WIDTH];
    pixel *base_u = &lcu->ref.u[x_local/2 + (y_local * LCU_WIDTH)/4];
    pixel *base_v = &lcu->ref.v[x_local/2 + (y_local * LCU_WIDTH)/4];
    int32_t base_stride = LCU_WIDTH;

    pixel pred_y[LCU_WIDTH*LCU_WIDTH];
    pixel pred_u[LCU_WIDTH*LCU_WIDTH>>2];
    pixel pred_v[LCU_WIDTH*LCU_WIDTH>>2];
    int32_t pred_stride = LCU_WIDTH;

    coefficient coeff_y[LCU_WIDTH*LCU_WIDTH];
    coefficient coeff_u[LCU_WIDTH*LCU_WIDTH>>2];
    coefficient coeff_v[LCU_WIDTH*LCU_WIDTH>>2];
    coefficient *orig_coeff_y = &lcu->coeff.y[x_local + y_local * LCU_WIDTH];
    coefficient *orig_coeff_u = &lcu->coeff.u[x_local/2 + (y_local * LCU_WIDTH)/4];
    coefficient *orig_coeff_v = &lcu->coeff.v[x_local/2 + (y_local * LCU_WIDTH)/4];
    int32_t coeff_stride = LCU_WIDTH;

    // Quant and transform here...
    int16_t block[LCU_WIDTH*LCU_WIDTH>>2];
    int16_t pre_quant_coeff[LCU_WIDTH*LCU_WIDTH>>2];

    // INTRA PREDICTION

    uint32_t ac_sum = 0;
    uint8_t scan_idx_luma   = SCAN_DIAG;
    uint8_t scan_idx_chroma = SCAN_DIAG;

    int32_t x_pu = x_local >> 2;
    int32_t y_pu = y_local >> 2;

    int cbf_y;

    #if OPTIMIZATION_SKIP_RESIDUAL_ON_THRESHOLD
    uint32_t residual_sum = 0;
    #endif

    // Pick coeff scan mode according to intra prediction mode.
    if (cur_cu->type == CU_INTRA) {
      int pu_index = PU_INDEX(x_pu, y_pu);
      int luma_mode = cur_cu->intra[pu_index].mode;
      int chroma_mode = cur_cu->intra[0].mode_chroma;
      if (chroma_mode == 36) {
        chroma_mode = luma_mode;
      }
      scan_idx_luma = SCAN_DIAG;
      scan_idx_chroma = SCAN_DIAG;

      // Scan mode is diagonal, except for 4x4+8x8 luma and 4x4 chroma, where:
      // - angular 6-14 = vertical
      // - angular 22-30 = horizontal
      if (width <= 8) {
        if (luma_mode >= 6 && luma_mode <= 14) {
          scan_idx_luma = SCAN_VER;
        } else if (luma_mode >= 22 && luma_mode <= 30) {
          scan_idx_luma = SCAN_HOR;
        }

        if (chroma_mode >= 6 && chroma_mode <= 14) {
          scan_idx_chroma = SCAN_VER;
        } else if (chroma_mode >= 22 && chroma_mode <= 30) {
          scan_idx_chroma = SCAN_HOR;
        }
      }
    }


    // Copy Luma and Chroma to the pred-block
    for(y = 0; y < width; y++) {
      for(x = 0; x < width; x++) {
        pred_y[x+y*pred_stride]=recbase_y[x+y*recbase_stride];
      }
    }

    // INTRA PREDICTION ENDS HERE

    // Get residual by subtracting prediction
    i = 0;
    ac_sum = 0;

    for (y = 0; y < width; y++) {
      for (x = 0; x < width; x++) {
        block[i] = ((int16_t)base_y[x + y * base_stride]) -
                   pred_y[x + y * pred_stride];
        #if OPTIMIZATION_SKIP_RESIDUAL_ON_THRESHOLD
        residual_sum += block[i];
        #endif
        i++;
      }
    }
    #if OPTIMIZATION_SKIP_RESIDUAL_ON_THRESHOLD
    #define RESIDUAL_THRESHOLD 500
    if(residual_sum < RESIDUAL_THRESHOLD/(width)) {
      memset(block, 0, sizeof(int16_t)*(width)*(width));
    }
    #endif

    // For 4x4 blocks, check for transform skip
    if(width == 4 && encoder->trskip_enable) {
      int i;
      coefficient temp_block[16];  coefficient temp_coeff[16];
      coefficient temp_block2[16]; coefficient temp_coeff2[16];
      uint32_t cost = 0,cost2 = 0;
      uint32_t coeffcost = 0,coeffcost2 = 0;

      // Test for transform skip
      transformskip(encoder, block,pre_quant_coeff,width);
      if (encoder->rdoq_enable) {
        rdoq(encoder_state, cabac, pre_quant_coeff, temp_coeff, 4, 4, &ac_sum, 0, scan_idx_luma, cur_cu->type,0);
      } else {
        quant(encoder_state, pre_quant_coeff, temp_coeff, 4, 4, &ac_sum, 0, scan_idx_luma, cur_cu->type);
      }
      dequant(encoder_state, temp_coeff, pre_quant_coeff, 4, 4, 0, cur_cu->type);
      itransformskip(encoder, temp_block,pre_quant_coeff,width);

      transform2d(encoder, block,pre_quant_coeff,width,0);
      if (encoder->rdoq_enable) {
        rdoq(encoder_state, cabac, pre_quant_coeff, temp_coeff2, 4, 4, &ac_sum, 0, scan_idx_luma, cur_cu->type,0);
      } else {
        quant(encoder_state, pre_quant_coeff, temp_coeff2, 4, 4, &ac_sum, 0, scan_idx_luma, cur_cu->type);
      }
      dequant(encoder_state, temp_coeff2, pre_quant_coeff, 4, 4, 0, cur_cu->type);
      itransform2d(encoder, temp_block2,pre_quant_coeff,width,0);

      // SSD between original and reconstructed
      for (i = 0; i < 16; i++) {
        int diff = temp_block[i]-block[i];
        cost += diff*diff;

        diff = temp_block2[i] - block[i];
        cost2 += diff*diff;
      }

      // Simple RDO
      if(encoder->rdo == 1) {
        // SSD between reconstruction and original + sum of coeffs
        for (i = 0; i < 16; i++) {
          coeffcost += abs((int)temp_coeff[i]);
          coeffcost2 += abs((int)temp_coeff2[i]);
        }
        cost += (1 + coeffcost + (coeffcost>>1))*((int)encoder_state->cur_lambda_cost+0.5);
        cost2 += (coeffcost2 + (coeffcost2>>1))*((int)encoder_state->cur_lambda_cost+0.5);
        // Full RDO
      } else if(encoder->rdo == 2) {
        coeffcost = get_coeff_cost(encoder, cabac, temp_coeff, 4, 0, scan_idx_luma);
        coeffcost2 = get_coeff_cost(encoder, cabac, temp_coeff2, 4, 0, scan_idx_luma);

        cost  += coeffcost*((int)encoder_state->cur_lambda_cost+0.5);
        cost2 += coeffcost2*((int)encoder_state->cur_lambda_cost+0.5);
      }

      cur_cu->intra[PU_INDEX(x_pu, y_pu)].tr_skip = (cost < cost2);
    }

    // Transform and quant residual to coeffs
    if(width == 4 && cur_cu->intra[PU_INDEX(x_pu, y_pu)].tr_skip) {
      transformskip(encoder, block,pre_quant_coeff,width);
    } else {
      transform2d(encoder, block,pre_quant_coeff,width,0);
    }

    if (encoder->rdoq_enable) {
      rdoq(encoder_state, cabac, pre_quant_coeff, coeff_y, width, width, &ac_sum, 0,
           scan_idx_luma, cur_cu->type, cur_cu->tr_depth-cur_cu->depth);
    } else {
      quant(encoder_state, pre_quant_coeff, coeff_y, width, width, &ac_sum, 0, scan_idx_luma, cur_cu->type);
    }

    // Check for non-zero coeffs
    cbf_y = 0;
    for (i = 0; i < width * width; i++) {
      if (coeff_y[i] != 0) {
        // Found one, we can break here
        cbf_y = 1;
        if (depth <= MAX_DEPTH) {
          int d;
          for (d = 0; d <= depth; ++d) {
            cur_cu->coeff_top_y[d] = 1;
          }
        } else {
          int pu_index = (x_pu & 1) + 2 * (y_pu & 1);
          int d;
          cur_cu->coeff_top_y[depth + pu_index] = 1;
          for (d = 0; d < depth; ++d) {
            cur_cu->coeff_top_y[d] = 1;
          }
        }
        break;
      }
    }

    if (cbf_y) {
      // Combine inverese quantized coefficients with the prediction to get
      // reconstructed image.
      //picture_set_block_residual(cur_pic,x_cu,y_cu,depth,1);
      i = 0;
      for (y = 0; y < width; y++) {
        for (x = 0; x < width; x++) {
          orig_coeff_y[x + y * coeff_stride] = coeff_y[i];
          i++;
        }
      }

      dequant(encoder_state, coeff_y, pre_quant_coeff, width, width, 0, cur_cu->type);
      if(width == 4 && cur_cu->intra[PU_INDEX(x_pu, y_pu)].tr_skip) {
        itransformskip(encoder, block,pre_quant_coeff,width);
      } else {
        itransform2d(encoder, block,pre_quant_coeff,width,0);
      }

      i = 0;

      for (y = 0; y < width; y++) {
        for (x = 0; x < width; x++) {
          int val = block[i++] + pred_y[x + y * pred_stride];
          //TODO: support 10+bits
          recbase_y[x + y * recbase_stride] = (pixel)CLIP(0, 255, val);
        }
      }
    } else {
      // Without coefficients, just copy the prediction as the reconstructed image.
      for (y = 0; y < width; y++) {
        for (x = 0; x < width; x++) {
          recbase_y[x + y * recbase_stride] = (pixel)CLIP(0, 255, pred_y[x + y * pred_stride]);
        }
      }
    }

    // If luma is 4x4, do chroma for the 8x8 luma area when handling the top
    // left PU because the coordinates are correct.
    if (depth <= MAX_DEPTH || (x_pu % 2 == 0 && y_pu % 2 == 0)) {
      int chroma_depth = (depth == MAX_PU_DEPTH ? depth - 1 : depth);
      int chroma_size = LCU_CHROMA_SIZE >> (chroma_depth * 2);

      // These are some weird indices for quant and dequant and should be
      // replaced later with color_index.
      int color_type_u = 2;
      int color_type_v = 3;

      for(y = 0; y < width_c; y++) {
        for(x = 0; x < width_c; x++) {
          pred_u[x+y*(pred_stride>>1)]=recbase_u[x+y*(recbase_stride>>1)];
          pred_v[x+y*(pred_stride>>1)]=recbase_v[x+y*(recbase_stride>>1)];
        }
      }

      transform_chroma(encoder_state, cabac, cur_cu, chroma_depth, base_u, pred_u, coeff_u, scan_idx_chroma, pre_quant_coeff, block);
      for (i = 0; i < chroma_size; i++) {
        if (coeff_u[i] != 0) {
          int d;
          for (d = 0; d <= depth; ++d) {
            cur_cu->coeff_top_u[d] = 1;
          }
          break;
        }
      }
      transform_chroma(encoder_state, cabac, cur_cu, chroma_depth, base_v, pred_v, coeff_v, scan_idx_chroma, pre_quant_coeff, block);
      for (i = 0; i < chroma_size; i++) {
        if (coeff_v[i] != 0) {
          int d;
          for (d = 0; d <= depth; ++d) {
            cur_cu->coeff_top_v[d] = 1;
          }
          break;
        }
      }

      // Save coefficients to cu.
      if (cur_cu->coeff_top_u[depth] || cur_cu->coeff_top_v[depth]) {
        i = 0;
        for (y = 0; y < width_c; y++) {
          for (x = 0; x < width_c; x++) {
            orig_coeff_u[x + y * (coeff_stride>>1)] = coeff_u[i];
            orig_coeff_v[x + y * (coeff_stride>>1)] = coeff_v[i];
            i++;
          }
        }
      }

      reconstruct_chroma(encoder_state, cur_cu, chroma_depth,
                         cur_cu->coeff_top_u[depth],
                         coeff_u, recbase_u, pred_u, color_type_u,
                         pre_quant_coeff, block);
      reconstruct_chroma(encoder_state, cur_cu, chroma_depth,
                         cur_cu->coeff_top_v[depth],
                         coeff_v, recbase_v, pred_v, color_type_v,
                         pre_quant_coeff, block);
    }

    return;
  }

  // end Residual Coding
}

static void encode_transform_unit(encoder_state * const encoder_state, cabac_data *cabac,
                                  int x_pu, int y_pu, int depth, int tr_depth)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  const picture * const cur_pic = encoder_state->cur_pic;
  uint8_t width = LCU_WIDTH >> depth;
  uint8_t width_c = (depth == MAX_PU_DEPTH ? width : width / 2);

  int x_cu = x_pu / 2;
  int y_cu = y_pu / 2;
  cu_info *cur_cu = &cur_pic->cu_array[MAX_DEPTH][x_cu + y_cu * (cur_pic->width_in_lcu << MAX_DEPTH)];

  coefficient coeff_y[LCU_WIDTH*LCU_WIDTH+1];
  coefficient coeff_u[LCU_WIDTH*LCU_WIDTH>>2];
  coefficient coeff_v[LCU_WIDTH*LCU_WIDTH>>2];
  int32_t coeff_stride = cur_pic->width;

  uint32_t ctx_idx;
  int8_t scan_idx = SCAN_DIAG;
  uint32_t dir_mode;

  int cbf_y;
  if (depth <= MAX_DEPTH) {
    cbf_y = cur_cu->coeff_top_y[depth];
  } else {
    int pu_index = x_pu % 2 + 2 * (y_pu % 2);
    cbf_y = cur_cu->coeff_top_y[depth + pu_index];
  }

  if (cbf_y) {
    int x = x_pu * (LCU_WIDTH >> MAX_PU_DEPTH);
    int y = y_pu * (LCU_WIDTH >> MAX_PU_DEPTH);
    coefficient *orig_pos = &cur_pic->coeff_y[x + y * cur_pic->width];
    for (y = 0; y < width; y++) {
      for (x = 0; x < width; x++) {
        coeff_y[x+y*width] = orig_pos[x];
      }
      orig_pos += coeff_stride;
    }
  }

  switch (width) {
    case  2: ctx_idx = 6; break;
    case  4: ctx_idx = 5; break;
    case  8: ctx_idx = 4; break;
    case 16: ctx_idx = 3; break;
    case 32: ctx_idx = 2; break;
    case 64: ctx_idx = 1; break;
    default: ctx_idx = 0; break;
  }

  ctx_idx -= tr_depth;

  // CoeffNxN
  // Residual Coding
  if (cbf_y) {
    scan_idx = SCAN_DIAG;
    if (cur_cu->type == CU_INTRA) {
      // Luma (Intra) scanmode
      if (depth <= MAX_DEPTH) {
        dir_mode = cur_cu->intra[0].mode;
      } else {
        int pu_index = x_pu % 2 + 2 * (y_pu % 2);
        dir_mode = cur_cu->intra[pu_index].mode;
      }

      // Scan mode is diagonal, except for 4x4 and 8x8, where:
      // - angular 6-14 = vertical
      // - angular 22-30 = horizontal
      if (width <= 8) {
        if (dir_mode >= 6 && dir_mode <= 14) {
          scan_idx = SCAN_VER;
        } else if (dir_mode >= 22 && dir_mode <= 30) {
          scan_idx = SCAN_HOR;
        }
      }
    }

    encode_coeff_nxn(encoder, cabac, coeff_y, width, 0, scan_idx, cur_cu->intra[PU_INDEX(x_pu, y_pu)].tr_skip);
  }

  if (depth == MAX_DEPTH + 1 && !(x_pu % 2 && y_pu % 2)) {
    // For size 4x4 luma transform the corresponding chroma transforms are
    // also of size 4x4 covering 8x8 luma pixels. The residual is coded
    // in the last transform unit so for the other ones, don't do anything.
    return;
  }

  if (cur_cu->coeff_top_u[depth] || cur_cu->coeff_top_v[depth]) {
    int x, y;
    coefficient *orig_pos_u, *orig_pos_v;

    if (depth <= MAX_DEPTH) {
      x = x_pu * (LCU_WIDTH >> (MAX_PU_DEPTH + 1));
      y = y_pu * (LCU_WIDTH >> (MAX_PU_DEPTH + 1));
    } else {
      // for 4x4 select top left pixel of the CU.
      x = x_cu * (LCU_WIDTH >> (MAX_DEPTH + 1));
      y = y_cu * (LCU_WIDTH >> (MAX_DEPTH + 1));
    }
    orig_pos_u = &cur_pic->coeff_u[x + y * (cur_pic->width >> 1)];
    orig_pos_v = &cur_pic->coeff_v[x + y * (cur_pic->width >> 1)];
    for (y = 0; y < (width_c); y++) {
      for (x = 0; x < (width_c); x++) {
        coeff_u[x+y*(width_c)] = orig_pos_u[x];
        coeff_v[x+y*(width_c)] = orig_pos_v[x];
      }
      orig_pos_u += coeff_stride>>1;
      orig_pos_v += coeff_stride>>1;
    }

    if(cur_cu->type == CU_INTER) {
      scan_idx = SCAN_DIAG;
    } else {
      // Chroma scanmode
      ctx_idx++;
      dir_mode = cur_cu->intra[0].mode_chroma;

      if (dir_mode == 36) {
        // TODO: support NxN
        dir_mode = cur_cu->intra[0].mode;
      }

      scan_idx = SCAN_DIAG;

      if (ctx_idx > 4 && ctx_idx < 7) { // if multiple scans supported for transform size
        // mode is diagonal, except for 4x4 and 8x8, where:
        // - angular 6-14 = vertical
        // - angular 22-30 = horizontal
        scan_idx = abs((int32_t) dir_mode - 26) < 5 ? 1 : (abs((int32_t)dir_mode - 10) < 5 ? 2 : 0);
      }
    }

    if (cur_cu->coeff_top_u[depth]) {
      encode_coeff_nxn(encoder, cabac, coeff_u, width_c, 2, scan_idx, 0);
    }

    if (cur_cu->coeff_top_v[depth]) {
      encode_coeff_nxn(encoder, cabac, coeff_v, width_c, 2, scan_idx, 0);
    }
  }
}

/**
 * \param encoder
 * \param x_pu            Prediction units' x coordinate.
 * \param y_pu            Prediction units' y coordinate.
 * \param depth           Depth from LCU.
 * \param tr_depth        Depth from last CU.
 * \param parent_coeff_u  What was signaled at previous level for cbf_cb.
 * \param parent_coeff_v  What was signlaed at previous level for cbf_cr.
 */
void encode_transform_coeff(encoder_state * const encoder_state, cabac_data *cabac, int32_t x_pu,int32_t y_pu,
                            int8_t depth, int8_t tr_depth, uint8_t parent_coeff_u, uint8_t parent_coeff_v)
{
  int32_t x_cu = x_pu / 2;
  int32_t y_cu = y_pu / 2;
  const picture * const cur_pic = encoder_state->cur_pic;
  cu_info *cur_cu = &cur_pic->cu_array[MAX_DEPTH][x_cu + y_cu * (cur_pic->width_in_lcu << MAX_DEPTH)];

  // NxN signifies implicit transform split at the first transform level.
  // There is a similar implicit split for inter, but it is only used when
  // transform hierarchy is not in use.
  int intra_split_flag = (cur_cu->type == CU_INTRA && cur_cu->part_size == SIZE_NxN);

  // The implicit split by intra NxN is not counted towards max_tr_depth.
  int max_tr_depth = (cur_cu->type == CU_INTRA ? TR_DEPTH_INTRA + intra_split_flag : TR_DEPTH_INTER);

  int8_t split = (cur_cu->tr_depth > depth);

  int8_t cb_flag_u = cur_cu->coeff_top_u[depth];
  int8_t cb_flag_v = cur_cu->coeff_top_v[depth];
  int cb_flag_y;
  if (depth <= MAX_DEPTH) {
    cb_flag_y = cur_cu->coeff_top_y[depth];
  } else {
    int pu_index = x_pu % 2 + 2 * (y_pu % 2);
    cb_flag_y = cur_cu->coeff_top_y[depth + pu_index];
  }

  // The split_transform_flag is not signaled when:
  // - transform size is greater than 32 (depth == 0)
  // - transform size is 4 (depth == MAX_PU_DEPTH)
  // - transform depth is max
  // - cu is intra NxN and it's the first split
  if (depth > 0 &&
      depth < MAX_PU_DEPTH &&
      tr_depth < max_tr_depth &&
      !(intra_split_flag && tr_depth == 0))
  {
    cabac->ctx = &(cabac->ctx_trans_subdiv_model[5 - ((g_convert_to_bit[LCU_WIDTH] + 2) - depth)]);
    CABAC_BIN(cabac, split, "split_transform_flag");
  }

  // Chroma cb flags are not signaled when one of the following:
  // - transform size is 4 (2x2 chroma transform doesn't exist)
  // - they have already been signaled to 0 previously
  // When they are not present they are inferred to be 0, except for size 4
  // when the flags from previous level are used.
  if (depth < MAX_PU_DEPTH) {
    cabac->ctx = &(cabac->ctx_qt_cbf_model_chroma[tr_depth]);
    if (tr_depth == 0 || parent_coeff_u) {
      CABAC_BIN(cabac, cb_flag_u, "cbf_cb");
    }
    if (tr_depth == 0 || parent_coeff_v) {
      CABAC_BIN(cabac, cb_flag_v, "cbf_cr");
    }
  }

  if (split) {
    uint8_t pu_offset = 1 << (MAX_PU_DEPTH - (depth + 1));
    encode_transform_coeff(encoder_state, cabac, x_pu, y_pu, depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
    encode_transform_coeff(encoder_state, cabac, x_pu + pu_offset, y_pu,  depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
    encode_transform_coeff(encoder_state, cabac, x_pu, y_pu + pu_offset,  depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
    encode_transform_coeff(encoder_state, cabac, x_pu + pu_offset, y_pu + pu_offset,  depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
    return;
  }

  // Luma coded block flag is signaled when one of the following:
  // - prediction mode is intra
  // - transform depth > 0
  // - we have chroma coefficients at this level
  // When it is not present, it is inferred to be 1.
  if(cur_cu->type == CU_INTRA || tr_depth > 0 || cb_flag_u || cb_flag_v) {
      cabac->ctx = &(cabac->ctx_qt_cbf_model_luma[!tr_depth]);
      CABAC_BIN(cabac, cb_flag_y, "cbf_luma");
  }

  if (cb_flag_y | cb_flag_u | cb_flag_v) {
    encode_transform_unit(encoder_state, cabac, x_pu, y_pu, depth, tr_depth);
  }
}

void encode_coeff_nxn(const encoder_control * const encoder, cabac_data *cabac, coefficient *coeff, uint8_t width,
                      uint8_t type, int8_t scan_mode, int8_t tr_skip)
{
  int c1 = 1;
  uint8_t last_coeff_x = 0;
  uint8_t last_coeff_y = 0;
  int32_t i;
  uint32_t sig_coeffgroup_flag[64];

  uint32_t num_nonzero = 0;
  int32_t scan_pos_last = -1;
  int32_t pos_last = 0;
  int32_t shift   = 4>>1;
  int8_t be_valid = ENABLE_SIGN_HIDING;
  int32_t scan_pos_sig;
  int32_t last_scan_set;
  uint32_t go_rice_param = 0;
  uint32_t blk_pos, pos_y, pos_x, sig, ctx_sig;

  // CONSTANTS
  const uint32_t num_blk_side    = width >> shift;
  const uint32_t log2_block_size = g_convert_to_bit[width] + 2;
  const uint32_t *scan           =
    g_sig_last_scan[scan_mode][log2_block_size - 1];
  const uint32_t *scan_cg = g_sig_last_scan_cg[log2_block_size - 2][scan_mode];

  // Init base contexts according to block type
  cabac_ctx *base_coeff_group_ctx = &(cabac->ctx_cu_sig_coeff_group_model[type]);
  cabac_ctx *baseCtx           = (type == 0) ? &(cabac->ctx_cu_sig_model_luma[0]) :
                                 &(cabac->ctx_cu_sig_model_chroma[0]);
  memset(sig_coeffgroup_flag,0,sizeof(uint32_t)*64);

  // transform skip flag
  if(width == 4 && encoder->trskip_enable) {
    cabac->ctx = (type == 0) ? &(cabac->ctx_transform_skip_model_luma) : &(cabac->ctx_transform_skip_model_chroma);
    CABAC_BIN(cabac, tr_skip, "transform_skip_flag");
  }

  // Count non-zero coeffs
  for (i = 0; i < width * width; i++) {
    if (coeff[i] != 0) {
      num_nonzero++;
    }
  }

  scan_pos_last = -1;

  // Significance mapping
  while (num_nonzero > 0) {
    pos_last = scan[++scan_pos_last];
#define POSY (pos_last >> log2_block_size)
#define POSX (pos_last - ( POSY << log2_block_size ))

    if (coeff[pos_last] != 0) {
      sig_coeffgroup_flag[(num_blk_side * (POSY >> shift) + (POSX >> shift))] = 1;
    }

    num_nonzero -= (coeff[pos_last] != 0) ? 1 : 0;
    #undef POSY
    #undef POSX
  }

  last_coeff_x = pos_last & (width - 1);
  last_coeff_y = (uint8_t)(pos_last >> log2_block_size);

  // Code last_coeff_x and last_coeff_y
  encode_last_significant_xy(cabac, last_coeff_x, last_coeff_y, width, width,
                             type, scan_mode);

  scan_pos_sig  = scan_pos_last;
  last_scan_set = (scan_pos_last >> 4);

  // significant_coeff_flag
  for (i = last_scan_set; i >= 0; i--) {
    int32_t sub_pos        = i << 4; // LOG2_SCAN_SET_SIZE;
    int32_t abs_coeff[16];
    int32_t cg_blk_pos     = scan_cg[i];
    int32_t cg_pos_y       = cg_blk_pos / num_blk_side;
    int32_t cg_pos_x       = cg_blk_pos - (cg_pos_y * num_blk_side);

    uint32_t coeff_signs   = 0;
    int32_t last_nz_pos_in_cg = -1;
    int32_t first_nz_pos_in_cg = 16;
    int32_t num_non_zero = 0;
    go_rice_param = 0;

    if (scan_pos_sig == scan_pos_last) {
      abs_coeff[0] = abs(coeff[pos_last]);
      coeff_signs  = (coeff[pos_last] < 0);
      num_non_zero = 1;
      last_nz_pos_in_cg  = scan_pos_sig;
      first_nz_pos_in_cg = scan_pos_sig;
      scan_pos_sig--;
    }

    if (i == last_scan_set || i == 0) {
      sig_coeffgroup_flag[cg_blk_pos] = 1;
    } else {
      uint32_t sig_coeff_group   = (sig_coeffgroup_flag[cg_blk_pos] != 0);
      uint32_t ctx_sig  = context_get_sig_coeff_group(sig_coeffgroup_flag, cg_pos_x,
                                                      cg_pos_y, width);
      cabac->ctx = &base_coeff_group_ctx[ctx_sig];
      CABAC_BIN(cabac, sig_coeff_group, "significant_coeff_group");
    }

    if (sig_coeffgroup_flag[cg_blk_pos]) {
      int32_t pattern_sig_ctx = context_calc_pattern_sig_ctx(sig_coeffgroup_flag,
                                                             cg_pos_x, cg_pos_y, width);

      for (; scan_pos_sig >= sub_pos; scan_pos_sig--) {
        blk_pos = scan[scan_pos_sig];
        pos_y   = blk_pos >> log2_block_size;
        pos_x   = blk_pos - (pos_y << log2_block_size);
        sig    = (coeff[blk_pos] != 0) ? 1 : 0;

        if (scan_pos_sig > sub_pos || i == 0 || num_non_zero) {
          ctx_sig  = context_get_sig_ctx_inc(pattern_sig_ctx, scan_mode, pos_x, pos_y,
                                             log2_block_size, type);
          cabac->ctx = &baseCtx[ctx_sig];
          CABAC_BIN(cabac, sig, "significant_coeff_flag");
        }

        if (sig) {
          abs_coeff[num_non_zero] = abs(coeff[blk_pos]);
          coeff_signs              = 2 * coeff_signs + (coeff[blk_pos] < 0);
          num_non_zero++;

          if (last_nz_pos_in_cg == -1) {
            last_nz_pos_in_cg = scan_pos_sig;
          }

          first_nz_pos_in_cg  = scan_pos_sig;
        }
      }
    } else {
      scan_pos_sig = sub_pos - 1;
    }

    if (num_non_zero > 0) {
      int8_t sign_hidden = (last_nz_pos_in_cg - first_nz_pos_in_cg >=
                            4 /*SBH_THRESHOLD*/) ? 1 : 0;
      uint32_t ctx_set  = (i > 0 && type == 0) ? 2 : 0;
      cabac_ctx *base_ctx_mod;
      int32_t num_c1_flag, first_c2_flag_idx, idx, first_coeff2;

      if (c1 == 0) {
        ctx_set++;
      }

      c1 = 1;

      base_ctx_mod     = (type == 0) ? &(cabac->ctx_cu_one_model_luma[4 * ctx_set]) :
                         &(cabac->ctx_cu_one_model_chroma[4 * ctx_set]);
      num_c1_flag      = MIN(num_non_zero, C1FLAG_NUMBER);
      first_c2_flag_idx = -1;

      for (idx = 0; idx < num_c1_flag; idx++) {
        uint32_t symbol = (abs_coeff[idx] > 1) ? 1 : 0;
        cabac->ctx = &base_ctx_mod[c1];
        CABAC_BIN(cabac, symbol, "significant_coeff2_flag");

        if (symbol) {
          c1 = 0;

          if (first_c2_flag_idx == -1) {
            first_c2_flag_idx = idx;
          }
        } else if ((c1 < 3) && (c1 > 0)) {
          c1++;
        }
      }

      if (c1 == 0) {
        base_ctx_mod = (type == 0) ? &(cabac->ctx_cu_abs_model_luma[ctx_set]) :
                       &(cabac->ctx_cu_abs_model_chroma[ctx_set]);

        if (first_c2_flag_idx != -1) {
          uint8_t symbol = (abs_coeff[first_c2_flag_idx] > 2) ? 1 : 0;
          cabac->ctx      = &base_ctx_mod[0];
          CABAC_BIN(cabac,symbol,"first_c2_flag");
        }
      }

      if (be_valid && sign_hidden) {
        CABAC_BINS_EP(cabac, (coeff_signs >> 1), (num_non_zero - 1), "");
      } else {
        CABAC_BINS_EP(cabac, coeff_signs, num_non_zero, "");
      }

      if (c1 == 0 || num_non_zero > C1FLAG_NUMBER) {
        first_coeff2 = 1;

        for (idx = 0; idx < num_non_zero; idx++) {
          int32_t base_level  = (idx < C1FLAG_NUMBER) ? (2 + first_coeff2) : 1;

          if (abs_coeff[idx] >= base_level) {
            cabac_write_coeff_remain(cabac, abs_coeff[idx] - base_level, go_rice_param);

            if (abs_coeff[idx] > 3 * (1 << go_rice_param)) {
              go_rice_param = MIN(go_rice_param + 1, 4);
            }
          }

          if (abs_coeff[idx] >= 2) {
            first_coeff2 = 0;
          }
        }
      }
    }
  }
}

/*!
 \brief Encode (X,Y) position of the last significant coefficient
 \param lastpos_x X component of last coefficient
 \param lastpos_y Y component of last coefficient
 \param width  Block width
 \param height Block height
 \param type plane type / luminance or chrominance
 \param scan scan type (diag, hor, ver)

 This method encodes the X and Y component within a block of the last significant coefficient.
*/
void encode_last_significant_xy(cabac_data *cabac,
                                uint8_t lastpos_x, uint8_t lastpos_y,
                                uint8_t width, uint8_t height,
                                uint8_t type, uint8_t scan)
{
  uint8_t offset_x  = type?0:((TOBITS(width)*3) + ((TOBITS(width)+1)>>2)),offset_y = offset_x;
  uint8_t shift_x   = type?(TOBITS(width)):((TOBITS(width)+3)>>2), shift_y = shift_x;
  int group_idx_x;
  int group_idx_y;
  int last_x,last_y,i;
  cabac_ctx *base_ctx_x = (type ? cabac->ctx_cu_ctx_last_x_chroma : cabac->ctx_cu_ctx_last_x_luma);
  cabac_ctx *base_ctx_y = (type ? cabac->ctx_cu_ctx_last_y_chroma : cabac->ctx_cu_ctx_last_y_luma);

  if (scan == SCAN_VER) {
    SWAP( lastpos_x, lastpos_y,uint8_t );
  }

  group_idx_x   = g_group_idx[lastpos_x];
  group_idx_y   = g_group_idx[lastpos_y];

  // Last X binarization
  for (last_x = 0; last_x < group_idx_x ; last_x++) {
    cabac->ctx = &base_ctx_x[offset_x + (last_x >> shift_x)];
    CABAC_BIN(cabac,1,"LastSignificantX");
  }

  if (group_idx_x < g_group_idx[width - 1]) {
    cabac->ctx = &base_ctx_x[offset_x + (last_x >> shift_x)];
    CABAC_BIN(cabac,0,"LastSignificantX");
  }

  // Last Y binarization
  for (last_y = 0; last_y < group_idx_y ; last_y++) {
    cabac->ctx = &base_ctx_y[offset_y + (last_y >> shift_y)];
    CABAC_BIN(cabac,1,"LastSignificantY");
  }

  if (group_idx_y < g_group_idx[height - 1]) {
    cabac->ctx = &base_ctx_y[offset_y + (last_y >> shift_y)];
    CABAC_BIN(cabac,0,"LastSignificantY");
  }

  // Last X
  if (group_idx_x > 3) {
    lastpos_x -= g_min_in_group[group_idx_x];

    for (i = ((group_idx_x - 2) >> 1) - 1; i >= 0; i--) {
      CABAC_BIN_EP(cabac,(lastpos_x>>i) & 1,"LastSignificantX");
    }
  }

  // Last Y
  if (group_idx_y > 3) {
    lastpos_y -= g_min_in_group[group_idx_y];

    for (i = ((group_idx_y - 2) >> 1) - 1; i >= 0; i--) {
      CABAC_BIN_EP(cabac,(lastpos_y>>i) & 1,"LastSignificantY");
    }
  }

  // end LastSignificantXY
}
