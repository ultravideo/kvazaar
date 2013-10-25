/**
 * \file
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "encoder.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

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

int16_t g_lambda_cost[55];
uint32_t* g_sig_last_scan[3][7];

/* Local functions. */
static void add_checksum(encoder_control* encoder);

void init_sig_last_scan(uint32_t *buff_d, uint32_t *buff_h, uint32_t *buff_v,
                        int32_t width, int32_t height)
{
  uint32_t num_scan_pos  = width * width;
  uint32_t next_scan_pos = 0;
  int32_t  xx, yy, x, y;
  uint32_t scan_line;
  uint32_t blk_y, blk_x;
  uint32_t blk;
  uint32_t cnt = 0;

  if (width < 16) {
    uint32_t *buff_tmp = buff_d;

    if (width == 8) {
      buff_tmp = (uint32_t *)g_sig_last_scan_32x32;
    }

    for (scan_line = 0; next_scan_pos < num_scan_pos; scan_line++) {
      int    primary_dim  = scan_line;
      int    second_dim  = 0;

      while (primary_dim >= width) {
        second_dim++;
        primary_dim--;
      }

      while (primary_dim >= 0 && second_dim < width) {
        buff_tmp[next_scan_pos] = primary_dim * width + second_dim ;
        next_scan_pos++;
        second_dim++;
        primary_dim--;
      }
    }
  }

  if (width > 4) {
    uint32_t num_blk_side = width >> 2;
    uint32_t num_blks    = num_blk_side * num_blk_side;
    uint32_t log2_blk      = g_convert_to_bit[num_blk_side] + 1;

    for (blk = 0; blk < num_blks; blk++) {
      uint32_t init_blk_pos = g_sig_last_scan[SCAN_DIAG][log2_blk][blk];
      next_scan_pos   = 0;

      if (width == 32) {
        init_blk_pos = g_sig_last_scan_32x32[blk];
      }

      {
        uint32_t offset_y    = init_blk_pos / num_blk_side;
        uint32_t offset_x    = init_blk_pos - offset_y * num_blk_side;
        uint32_t offset_d    = 4 * (offset_x + offset_y * width);
        uint32_t offset_scan = 16 * blk;

        for (scan_line = 0; next_scan_pos < 16; scan_line++) {
          int    primary_dim  = scan_line;
          int    second_dim  = 0;

          //TODO: optimize
          while (primary_dim >= 4) {
            second_dim++;
            primary_dim--;
          }

          while (primary_dim >= 0 && second_dim < 4) {
            buff_d[next_scan_pos + offset_scan] = primary_dim * width +
                                                  second_dim + offset_d;
            next_scan_pos++;
            second_dim++;
            primary_dim--;
          }
        }
      }
    }
  }  
  
  if (width > 2) {
    uint32_t num_blk_side = width >> 2;

    for (blk_y = 0; blk_y < num_blk_side; blk_y++) {
      for (blk_x = 0; blk_x < num_blk_side; blk_x++) {
        uint32_t offset    = blk_y * 4 * width + blk_x * 4;

        for (y = 0; y < 4; y++) {
          for (x = 0; x < 4; x++) {
            buff_h[cnt] = y * width + x + offset;
            cnt ++;
          }
        }
      }
    }

    cnt = 0;

    for (blk_x = 0; blk_x < num_blk_side; blk_x++) {
      for (blk_y = 0; blk_y < num_blk_side; blk_y++) {
        uint32_t offset = blk_y * 4 * width + blk_x * 4;

        for (x = 0; x < 4; x++) {
          for (y = 0; y < 4; y++) {
            buff_v[cnt] = y * width + x + offset;
            cnt ++;
          }
        }
      }
    }
  } else {
    for (yy = 0; yy < height; yy++) {
      for (xx = 0; xx < width; xx++) {
        buff_h[cnt] = yy * width + xx;
        cnt ++;
  }
      }

    cnt = 0;

    for (xx = 0; xx < width; xx++) {
      for (yy = 0; yy < height; yy++) {
        buff_v[cnt] = yy * width + xx;
        cnt ++;
      }
    }
  }
}


void init_tables(void)
{
  int i;
  int c = 0;
  memset( g_convert_to_bit,-1, sizeof( g_convert_to_bit ) );  

  for (i = 4; i < (1 << 7); i *= 2) {
    g_convert_to_bit[i] = c;
    c++;
  }

  g_convert_to_bit[i] = c;
  c = 2;

  for (i = 0; i < 7; i++) {
    g_sig_last_scan[0][i] = (uint32_t*)malloc(c*c*sizeof(uint32_t));
    g_sig_last_scan[1][i] = (uint32_t*)malloc(c*c*sizeof(uint32_t));
    g_sig_last_scan[2][i] = (uint32_t*)malloc(c*c*sizeof(uint32_t));

    init_sig_last_scan(g_sig_last_scan[0][i], g_sig_last_scan[1][i],
                       g_sig_last_scan[2][i], c, c);
    c <<= 1;
  }

  // Lambda cost
  // TODO: cleanup
  //g_lambda_cost = (int16_t*)malloc(sizeof(int16_t)*55);
  for (i = 0; i < 55; i++) {
    if (i < 12) {
      g_lambda_cost[i] = 0;
    } else {
      g_lambda_cost[i] = (int16_t)sqrt(0.57 * pow(2.0, (i - 12) / 3));
    }

    //g_lambda_cost[i] = g_lambda_cost[i]*g_lambda_cost[i];
  }

}

void init_encoder_control(encoder_control* control,bitstream* output)
{
  control->stream = output;  
}

void init_encoder_input(encoder_input *input, FILE *inputfile,
                        int32_t width, int32_t height)
{
  input->file = inputfile;
  input->width = width;
  input->height = height;
  input->real_width = width;
  input->real_height = height;
  
  // If input dimensions are not divisible by the smallest block size, add 
  // pixels to the dimensions, so that they are. These extra pixels will be 
  // compressed along with the real ones but they will be cropped out before 
  // rendering.
  if (width % CU_MIN_SIZE_PIXELS) {
    input->width += CU_MIN_SIZE_PIXELS - (width % CU_MIN_SIZE_PIXELS);
  }

  if (height % CU_MIN_SIZE_PIXELS) {
    input->height += CU_MIN_SIZE_PIXELS - (height % CU_MIN_SIZE_PIXELS);
  }

  input->height_in_lcu = input->height / LCU_WIDTH;
  input->width_in_lcu  = input->width / LCU_WIDTH;

  // Add one extra LCU when image not divisible by LCU_WIDTH
  if (input->height_in_lcu * LCU_WIDTH < height) {
    input->height_in_lcu++;
  }

  if (input->width_in_lcu * LCU_WIDTH < width) {
    input->width_in_lcu++;
  }

  // Allocate the picture and CU array
  input->cur_pic = picture_init(input->width, input->height,
                                input->width_in_lcu,
                                input->height_in_lcu);

  if (!input->cur_pic) {
    printf("Error allocating picture!\r\n");
    exit(1);
  }

  #ifdef _DEBUG
  if (width != input->width || height != input->height) {
    printf("Picture buffer has been extended to be a multiple of the smallest block size:\r\n");
    printf("  Width = %d (%d), Height = %d (%d)\r\n", width, input->width, height,
           input->height);
  }
  #endif
  }  

void encode_one_frame(encoder_control* encoder)
{
  // output parameters before first frame
  if (encoder->frame == 0) {
    // Video Parameter Set (VPS)
    encode_vid_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer,
              encoder->stream->buffer_pos, 0, NAL_VPS_NUT, 0);
    bitstream_clear_buffer(encoder->stream);

    // Sequence Parameter Set (SPS)
    encode_seq_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer,
              encoder->stream->buffer_pos, 0, NAL_SPS_NUT, 0);
    bitstream_clear_buffer(encoder->stream);
        
    // Picture Parameter Set (PPS)
    encode_pic_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer,
              encoder->stream->buffer_pos, 0, NAL_PPS_NUT, 0);
    bitstream_clear_buffer(encoder->stream);

    // First slice is IDR
    cabac_start(&cabac);
    encoder->in.cur_pic->slicetype = SLICE_I;
    encoder->in.cur_pic->type = NAL_IDR_W_RADL;
    search_slice_data(encoder);

    encode_slice_header(encoder);
    bitstream_align(encoder->stream);
    encode_slice_data(encoder);
    cabac_flush(&cabac);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer,
              encoder->stream->buffer_pos, 0, NAL_IDR_W_RADL, 0);
    bitstream_clear_buffer(encoder->stream);
  } else {
    cabac_start(&cabac);
    encoder->in.cur_pic->slicetype = SLICE_P;
    encoder->in.cur_pic->type = NAL_TRAIL_R;
    search_slice_data(encoder);

    encode_slice_header(encoder);
    bitstream_align(encoder->stream);
    encode_slice_data(encoder);
    cabac_flush(&cabac);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer,
              encoder->stream->buffer_pos, 0, NAL_TRAIL_R, 0);
    bitstream_clear_buffer(encoder->stream);
  }  
  
  // Filtering
  if(encoder->deblock_enable) {
    filter_deblock(encoder);
  }
  // Calculate checksum
  add_checksum(encoder);
}

void fill_after_frame(FILE *file, unsigned height, unsigned array_width,
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

void read_and_fill_frame_data(FILE *file, unsigned width, unsigned height,
                              unsigned array_width, pixel *data)
{
  pixel* p = data;
  pixel* end = data + array_width * height;
  pixel fill_char;
  unsigned i;

  while (p < end) {
    // Read the beginning of the line from input.
    fread(p, sizeof(unsigned char), width, file);

    // Fill the rest with the last pixel value.
    fill_char = p[width - 1];

    for (i = width; i < array_width; ++i) {
      p[i] = fill_char;
    }

    p += array_width;
  }
}

void read_one_frame(FILE* file, encoder_control* encoder)
{
  encoder_input* in = &encoder->in;
  unsigned width = in->real_width;
  unsigned height = in->real_height;
  unsigned array_width = in->cur_pic->width;
  unsigned array_height = in->cur_pic->height;

  if (width != array_width) {
    // In the case of frames not being aligned on 8 bit borders, bits need to be copied to fill them in.
    read_and_fill_frame_data(file, width, height, array_width,
                             in->cur_pic->y_data);
    read_and_fill_frame_data(file, width >> 1, height >> 1, array_width >> 1,
                             in->cur_pic->u_data);
    read_and_fill_frame_data(file, width >> 1, height >> 1, array_width >> 1,
                             in->cur_pic->v_data);
  } else {
    // Otherwise the data can be read directly to the array.
    fread(in->cur_pic->y_data, sizeof(unsigned char),
          width * height, file);
    fread(in->cur_pic->u_data, sizeof(unsigned char), 
          (width >> 1) * (height >> 1), file);
    fread(in->cur_pic->v_data, sizeof(unsigned char),
          (width >> 1) * (height >> 1), file);
  }

  if (height != array_height) {
    fill_after_frame(file, height, array_width, array_height,
                     in->cur_pic->y_data);
    fill_after_frame(file, height >> 1, array_width >> 1, array_height >> 1,
                     in->cur_pic->u_data);
    fill_after_frame(file, height >> 1, array_width >> 1, array_height >> 1,
                     in->cur_pic->v_data);
  }
}

/**
 * \brief Add a checksum SEI message to the bitstream.
 * \param encoder The encoder.
 * \returns Void
 */
static void add_checksum(encoder_control* encoder)
{
  unsigned char checksum[3][SEI_HASH_MAX_LENGTH];
  uint32_t checksum_val;
  unsigned int i;

  picture_checksum(encoder->in.cur_pic, checksum);

  WRITE_U(encoder->stream, 132, 8, "sei_type");
  WRITE_U(encoder->stream, 13, 8, "size");
  WRITE_U(encoder->stream, 2, 8, "hash_type"); // 2 = checksum

  for (i = 0; i < 3; ++i) {
    // Pack bits into a single 32 bit uint instead of pushing them one byte 
    // at a time.
    checksum_val = (checksum[i][0] << 24) + (checksum[i][1] << 16) +
                   (checksum[i][2] << 8) + (checksum[i][3]);
    WRITE_U(encoder->stream, checksum_val, 32, "picture_checksum");
  }

  bitstream_align(encoder->stream);
  bitstream_flush(encoder->stream);
  nal_write(encoder->output, encoder->stream->buffer,
            encoder->stream->buffer_pos, 0, NAL_SUFFIT_SEI_NUT, 0);
  bitstream_clear_buffer(encoder->stream);
}

void encode_pic_parameter_set(encoder_control* encoder)
{
#ifdef _DEBUG
  printf("=========== Picture Parameter Set ID: 0 ===========\n");
#endif
  WRITE_UE(encoder->stream, 0, "pic_parameter_set_id");
  WRITE_UE(encoder->stream, 0, "seq_parameter_set_id");
  WRITE_U(encoder->stream, 0, 1, "dependent_slice_segments_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "output_flag_present_flag");
  WRITE_U(encoder->stream, 0, 3, "num_extra_slice_header_bits");
  WRITE_U(encoder->stream, ENABLE_SIGN_HIDING, 1, "sign_data_hiding_flag");
  WRITE_U(encoder->stream, 0, 1, "cabac_init_present_flag");

  WRITE_UE(encoder->stream, 0, "num_ref_idx_l0_default_active_minus1");
  WRITE_UE(encoder->stream, 0, "num_ref_idx_l1_default_active_minus1");
  WRITE_SE(encoder->stream, ((int8_t)encoder->QP)-26, "pic_init_qp_minus26");
  WRITE_U(encoder->stream, 0, 1, "constrained_intra_pred_flag");
  WRITE_U(encoder->stream, 0, 1, "transform_skip_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "cu_qp_delta_enabled_flag");
  //if cu_qp_delta_enabled_flag
  //WRITE_UE(encoder->stream, 0, "diff_cu_qp_delta_depth");

  //TODO: add QP offsets
  WRITE_SE(encoder->stream, 0, "pps_cb_qp_offset");
  WRITE_SE(encoder->stream, 0, "pps_cr_qp_offset");
  WRITE_U(encoder->stream, 0, 1, "pps_slice_chroma_qp_offsets_present_flag");
  WRITE_U(encoder->stream, 0, 1, "weighted_pred_flag");
  WRITE_U(encoder->stream, 0, 1, "weighted_bipred_idc");

  //WRITE_U(encoder->stream, 0, 1, "dependent_slices_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "transquant_bypass_enable_flag");
  WRITE_U(encoder->stream, 0, 1, "tiles_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "entropy_coding_sync_enabled_flag");
  //TODO: enable tiles for concurrency
  //IF tiles
  //ENDIF
  WRITE_U(encoder->stream, 0, 1, "loop_filter_across_slice_flag");
  WRITE_U(encoder->stream, 1, 1, "deblocking_filter_control_present_flag");
  //IF deblocking_filter
    WRITE_U(encoder->stream, 0, 1, "deblocking_filter_override_enabled_flag");
  WRITE_U(encoder->stream, encoder->deblock_enable ? 0 : 1, 1,
          "pps_disable_deblocking_filter_flag");

    //IF !disabled
  if (encoder->deblock_enable) {
     WRITE_SE(encoder->stream, encoder->beta_offset_div2, "beta_offset_div2");
     WRITE_SE(encoder->stream, encoder->tc_offset_div2, "tc_offset_div2");
    }

    //ENDIF
  //ENDIF
  WRITE_U(encoder->stream, 0, 1, "pps_scaling_list_data_present_flag");
  //IF scaling_list
  //ENDIF
  WRITE_U(encoder->stream, 0, 1, "lists_modification_present_flag");
  WRITE_UE(encoder->stream, 0, "log2_parallel_merge_level_minus2");
  WRITE_U(encoder->stream, 0, 1, "slice_segment_header_extension_present_flag");
  WRITE_U(encoder->stream, 0, 1, "pps_extension_flag");
}

void encode_PTL(encoder_control *encoder)
{
  int i;
  // PTL
  // Profile Tier
  WRITE_U(encoder->stream, 0, 2, "XXX_profile_space[]");
  WRITE_U(encoder->stream, 0, 1, "XXX_tier_flag[]");
  WRITE_U(encoder->stream, 0, 5, "XXX_profile_idc[]");
  WRITE_U(encoder->stream, 0, 32, "XXX_profile_compatibility_flag[][j]");

  WRITE_U(encoder->stream, 1, 1, "general_progressive_source_flag");
  WRITE_U(encoder->stream, 0, 1, "general_interlaced_source_flag");
  WRITE_U(encoder->stream, 0, 1, "general_non_packed_constraint_flag");
  WRITE_U(encoder->stream, 0, 1, "general_frame_only_constraint_flag");

  WRITE_U(encoder->stream, 0, 32, "XXX_reserved_zero_44bits[0..31]");
  WRITE_U(encoder->stream, 0, 12, "XXX_reserved_zero_44bits[32..43]");
  
  // end Profile Tier
  
  WRITE_U(encoder->stream, 0, 8, "general_level_idc");
  
  WRITE_U(encoder->stream, 0, 1, "sub_layer_profile_present_flag");
  WRITE_U(encoder->stream, 0, 1, "sub_layer_level_present_flag");

  for (i = 1; i < 8; i++) {
    WRITE_U(encoder->stream, 0, 2, "reserved_zero_2bits");
  }
  
  // end PTL
}

void encode_seq_parameter_set(encoder_control* encoder)
{
  encoder_input* const in = &encoder->in;

#ifdef _DEBUG
  printf("=========== Sequence Parameter Set ID: 0 ===========\n");
#endif

  // TODO: profile IDC and level IDC should be defined later on
  WRITE_U(encoder->stream, 0, 4, "sps_video_parameter_set_id");
  WRITE_U(encoder->stream, 1, 3, "sps_max_sub_layers_minus1");
  WRITE_U(encoder->stream, 0, 1, "sps_temporal_id_nesting_flag");
    
  encode_PTL(encoder);

  WRITE_UE(encoder->stream, 0, "sps_seq_parameter_set_id");
  WRITE_UE(encoder->stream, encoder->in.video_format,
           "chroma_format_idc");

  if (encoder->in.video_format == 3) {
    WRITE_U(encoder->stream, 0, 1, "separate_colour_plane_flag");
  }

  WRITE_UE(encoder->stream, encoder->in.width, "pic_width_in_luma_samples");
  WRITE_UE(encoder->stream, encoder->in.height, "pic_height_in_luma_samples");

  if (in->width != in->real_width || in->height != in->real_height) {
    // The standard does not seem to allow setting conf_win values such that
    // the number of luma samples is not a multiple of 2. Options are to either
    // hide one line or show an extra line of non-video. Neither seems like a
    // very good option, so let's not even try.
    assert(!(in->width % 2));
    WRITE_U(encoder->stream, 1, 1, "conformance_window_flag");
    WRITE_UE(encoder->stream, 0, "conf_win_left_offset");
    WRITE_UE(encoder->stream, (in->width - in->real_width) >> 1,
             "conf_win_right_offset");
    WRITE_UE(encoder->stream, 0, "conf_win_top_offset");
    WRITE_UE(encoder->stream, (in->height - in->real_height) >> 1,
             "conf_win_bottom_offset");
  } else {
    WRITE_U(encoder->stream, 0, 1, "conformance_window_flag");
  }
  
  //IF window flag
  //END IF
  
  WRITE_UE(encoder->stream, encoder->bitdepth-8, "bit_depth_luma_minus8");
  WRITE_UE(encoder->stream, encoder->bitdepth-8, "bit_depth_chroma_minus8");
  WRITE_UE(encoder->stream, 0, "log2_max_pic_order_cnt_lsb_minus4");
  WRITE_U(encoder->stream, 0, 1, "sps_sub_layer_ordering_info_present_flag");

  //for each layer
  WRITE_UE(encoder->stream, 0, "sps_max_dec_pic_buffering");
  WRITE_UE(encoder->stream, 0, "sps_num_reorder_pics");
  WRITE_UE(encoder->stream, 0, "sps_max_latency_increase");
  //end for

  WRITE_UE(encoder->stream, MIN_SIZE-3, "log2_min_coding_block_size_minus3");
  WRITE_UE(encoder->stream, MAX_DEPTH, "log2_diff_max_min_coding_block_size");
  WRITE_UE(encoder->stream, 0, "log2_min_transform_block_size_minus2");   // 4x4
  WRITE_UE(encoder->stream, 3, "log2_diff_max_min_transform_block_size"); // 4x4...32x32
  WRITE_UE(encoder->stream, 2, "max_transform_hierarchy_depth_inter");
  WRITE_UE(encoder->stream, 2, "max_transform_hierarchy_depth_intra");
  
  // Use default scaling list
  WRITE_U(encoder->stream, ENABLE_SCALING_LIST, 1, "scaling_list_enable_flag");  
  #if ENABLE_SCALING_LIST == 1
    WRITE_U(encoder->stream, 0, 1, "sps_scaling_list_data_present_flag");
  #endif
  
  WRITE_U(encoder->stream, 0, 1, "amp_enabled_flag");
  WRITE_U(encoder->stream, encoder->sao_enable ? 1 : 0, 1,
          "sample_adaptive_offset_enabled_flag");
  WRITE_U(encoder->stream, ENABLE_PCM, 1, "pcm_enabled_flag");
  #if ENABLE_PCM == 1
    WRITE_U(encoder->stream, 7, 4, "pcm_sample_bit_depth_luma_minus1");
    WRITE_U(encoder->stream, 7, 4, "pcm_sample_bit_depth_chroma_minus1");
    WRITE_UE(encoder->stream, 0, "log2_min_pcm_coding_block_size_minus3");
    WRITE_UE(encoder->stream, 2, "log2_diff_max_min_pcm_coding_block_size");
    WRITE_U(encoder->stream, 1, 1, "pcm_loop_filter_disable_flag");
  #endif

  WRITE_UE(encoder->stream, 0, "num_short_term_ref_pic_sets"); 

  //IF num short term ref pic sets
  //ENDIF
  
  WRITE_U(encoder->stream, 0, 1, "long_term_ref_pics_present_flag");

  //IF long_term_ref_pics_present
  //ENDIF

  WRITE_U(encoder->stream, ENABLE_TEMPORAL_MVP, 1,
          "sps_temporal_mvp_enable_flag");
  WRITE_U(encoder->stream, 0, 1, "sps_strong_intra_smoothing_enable_flag");
  WRITE_U(encoder->stream, 0, 1, "vui_parameters_present_flag");

  //TODO: VUI?
  //encode_VUI(encoder);
  
	WRITE_U(encoder->stream, 0, 1, "sps_extension_flag");
}

void encode_vid_parameter_set(encoder_control* encoder)
{
  int i;
#ifdef _DEBUG
  printf("=========== Video Parameter Set ID: 0 ===========\n");
#endif

  WRITE_U(encoder->stream, 0, 4, "vps_video_parameter_set_id");
  WRITE_U(encoder->stream, 3, 2, "vps_reserved_three_2bits" );
  WRITE_U(encoder->stream, 0, 6, "vps_reserved_zero_6bits" );
  WRITE_U(encoder->stream, 1, 3, "vps_max_sub_layers_minus1");
  WRITE_U(encoder->stream, 0, 1, "vps_temporal_id_nesting_flag");
  WRITE_U(encoder->stream, 0xffff, 16, "vps_reserved_ffff_16bits");

  encode_PTL(encoder);

  WRITE_U(encoder->stream, 0, 1, "vps_sub_layer_ordering_info_present_flag");

  //for each layer
  for (i = 0; i < 1; i++) {
  WRITE_UE(encoder->stream, 1, "vps_max_dec_pic_buffering");
  WRITE_UE(encoder->stream, 0, "vps_num_reorder_pics");
  WRITE_UE(encoder->stream, 0, "vps_max_latency_increase");
  }

  WRITE_U(encoder->stream, 0, 6, "vps_max_nuh_reserved_zero_layer_id");
  WRITE_UE(encoder->stream, 0, "vps_max_op_sets_minus1");
  WRITE_U(encoder->stream, 0, 1, "vps_timing_info_present_flag");

  //IF timing info
  //END IF

	WRITE_U(encoder->stream, 0, 1, "vps_extension_flag");
}

void encode_VUI(encoder_control* encoder)
{
#ifdef _DEBUG
  printf("=========== VUI Set ID: 0 ===========\n");
#endif
  WRITE_U(encoder->stream, 0, 1, "aspect_ratio_info_present_flag");

  //IF aspect ratio info
  //ENDIF

  WRITE_U(encoder->stream, 0, 1, "overscan_info_present_flag");

  //IF overscan info
  //ENDIF

  WRITE_U(encoder->stream, 0, 1, "video_signal_type_present_flag");

  //IF video type
  //ENDIF

  WRITE_U(encoder->stream, 0, 1, "chroma_loc_info_present_flag");

  //IF chroma loc info
  //ENDIF
  
  WRITE_U(encoder->stream, 0, 1, "neutral_chroma_indication_flag");
  WRITE_U(encoder->stream, 0, 1, "field_seq_flag");
  WRITE_U(encoder->stream, 0, 1, "frame_field_info_present_flag");
  WRITE_U(encoder->stream, 0, 1, "default_display_window_flag");
  
  //IF default display window
  //ENDIF

  WRITE_U(encoder->stream, 0, 1, "vui_timing_info_present_flag");
  
  //IF timing info
  //ENDIF

  WRITE_U(encoder->stream, 0, 1, "bitstream_restriction_flag");
 
  //IF bitstream restriction
  //ENDIF
}

void encode_slice_header(encoder_control* encoder)
{
  picture *cur_pic = encoder->in.cur_pic;

#ifdef _DEBUG
  printf("=========== Slice ===========\n");
#endif

  WRITE_U(encoder->stream, 1, 1, "first_slice_segment_in_pic_flag");

  if (encoder->in.cur_pic->type >= NAL_BLA_W_LP
      && encoder->in.cur_pic->type <= NAL_RSV_IRAP_VCL23) {
    WRITE_U(encoder->stream, 1, 1, "no_output_of_prior_pics_flag");
  }

  WRITE_UE(encoder->stream, 0, "slice_pic_parameter_set_id");

  //WRITE_U(encoder->stream, 0, 1, "dependent_slice_segment_flag");
  
  WRITE_UE(encoder->stream, encoder->in.cur_pic->slicetype, "slice_type");

  // if !entropy_slice_flag
  
    //if output_flag_present_flag
      //WRITE_U(encoder->stream, 1, 1, "pic_output_flag");
    //end if
    //if( IdrPicFlag ) <- nal_unit_type == 5
  if (encoder->in.cur_pic->type != NAL_IDR_W_RADL
      && encoder->in.cur_pic->type != NAL_IDR_N_LP) {
      int j;
      int ref_negative = 1;
      int ref_positive = 0;
      WRITE_U(encoder->stream, encoder->frame&0xf, 4, "pic_order_cnt_lsb");
      WRITE_U(encoder->stream, 0, 1, "short_term_ref_pic_set_sps_flag");
      WRITE_UE(encoder->stream, ref_negative, "num_negative_pics");
      WRITE_UE(encoder->stream, ref_positive, "num_positive_pics");

    for (j = 0; j < ref_negative; j++) {
        WRITE_UE(encoder->stream, 0, "delta_poc_s0_minus1");
        WRITE_U(encoder->stream,1,1, "used_by_curr_pic_s0_flag");
      }

      //WRITE_UE(encoder->stream, 0, "short_term_ref_pic_set_idx");
    }

    //end if
  //end if
  if (encoder->sao_enable) {
    WRITE_U(encoder->stream, cur_pic->slice_sao_luma_flag, 1, "slice_sao_luma_flag");
    WRITE_U(encoder->stream, cur_pic->slice_sao_chroma_flag, 1, "slice_sao_chroma_flag");
  }
    
  if (encoder->in.cur_pic->slicetype != SLICE_I) {
      WRITE_U(encoder->stream, 0, 1, "num_ref_idx_active_override_flag");
      WRITE_UE(encoder->stream, 0, "five_minus_max_num_merge_cand");
    }

  if (encoder->in.cur_pic->slicetype == SLICE_B) {
      WRITE_U(encoder->stream, 0, 1, "mvd_l1_zero_flag");
    }

  // Skip flags that are not present
  // if !entropy_slice_flag
    WRITE_SE(encoder->stream, 0, "slice_qp_delta");
    //WRITE_U(encoder->stream, 1, 1, "alignment");
}


  // TODO: move somewhere else (sao.h?)
#define Y_INDEX 0
#define U_INDEX 1
#define V_INDEX 2
#define YUV_INDEX_END 3
#define NUM_SAO_OFFSETS 4

typedef enum { COLOR_Y = 0, COLOR_U = 1, COLOR_V = 2, NUM_COLORS } color_index;
typedef enum { SAO_TYPE_NONE = 0, SAO_TYPE_BAND, SAO_TYPE_EDGE } sao_type;
typedef enum { SAO_EO0 = 0, SAO_EO1, SAO_EO2, SAO_EO3, SAO_NUM_EO } sao_eo_class;
typedef enum { SAO_EO_CAT0 = 0, SAO_EO_CAT1, SAO_EO_CAT2, SAO_EO_CAT3, SAO_EO_CAT4, NUM_SAO_EDGE_CATEGORIES } sao_eo_cat;

typedef struct {
  sao_type type;
  sao_eo_class eo_class;
  int ddistortion;
  int merge_left_flag;
  int merge_up_flag;
  int offsets[NUM_SAO_OFFSETS];
} sao_info;

//#define SIGN3(x) ((x) > 0) ? +1 : ((x) == 0 ? 0 : -1)
#define SIGN3(x) (((x) > 0) - ((x) < 0))
#define NUM_SAO_EDGE_DIRS 4;

typedef struct {
  int x;
  int y;
} vector2d;

// Offsets of a and b in relation to c.
// dir_offset[dir][a or b]
// |       |   a   | a     |     a |
// | a c b |   c   |   c   |   c   |
// |       |   b   |     b | b     |
static const vector2d g_sao_edge_offsets[4][2] = { 
  { { 0, -1 }, { 0, 1 } },
  { { -1, 0 }, { 1, 0 } },
  { { -1, -1 }, { 1, 1 } },
  { { -1, 1 }, { 1, -1 } }
};
// Mapping of edge_idx values to eo-classes.
static const unsigned g_sao_eo_idx_to_eo_category[] = { 1, 2, 0, 3, 4 };
// Mapping relationships between a, b and c to eo_idx.
#define EO_IDX(a, b, c) (2 + SIGN3((c) - (a)) + SIGN3((c) - (b)))

/**
 * \param orig_data  Original pixel data. 64x64 for luma, 32x32 for chroma.
 * \param rec_data  Reconstructed pixel data. 64x64 for luma, 32x32 for chroma.
 * \param dir_offsets
 * \param is_chroma  0 for luma, 1 for chroma. Indicates 
 */
void calc_sao_edge_dir(const pixel *orig_data, const pixel *rec_data,
                       int eo_class, int block_width,
                       int cat_sum_cnt[2][NUM_SAO_EDGE_CATEGORIES])
{
  int y, x;
  vector2d a_ofs = g_sao_edge_offsets[eo_class][0];
  vector2d b_ofs = g_sao_edge_offsets[eo_class][1];
  // Arrays orig_data and rec_data are quarter size for chroma.

  // Don't sample the edge pixels because this function doesn't have access to
  // their neighbours.
  for (y = 1; y < block_width - 1; ++y) {
    for (x = 1; x < block_width - 1; ++x) {
      const pixel *c_data = &rec_data[y * block_width + x];
      pixel a = c_data[a_ofs.y * block_width + a_ofs.x];
      pixel c = c_data[0];
      pixel b = c_data[b_ofs.y * block_width + b_ofs.x];
      
      int eo_idx = EO_IDX(a, b, c);
      int eo_cat = g_sao_eo_idx_to_eo_category[eo_idx];

      cat_sum_cnt[0][eo_cat] += orig_data[y * block_width + x] - c;
      cat_sum_cnt[1][eo_cat] += 1;
    }
  }
}

void sao_reconstruct_color(pixel *rec_data, const sao_info *sao, color_index color)
{
  unsigned y, x;
  vector2d a_ofs = g_sao_edge_offsets[sao->eo_class][0];
  vector2d b_ofs = g_sao_edge_offsets[sao->eo_class][1];
  // Arrays orig_data and rec_data are quarter size for chroma.
  unsigned block_width = LCU_WIDTH >> !(color == COLOR_Y);

  for (y = 0; y < block_width; ++y) {
    for (x = 0; x < block_width; ++x) {
      pixel *c_data = &rec_data[y * block_width + x];
      pixel a = c_data[a_ofs.y * block_width + a_ofs.x];
      pixel c = c_data[0];
      pixel b = c_data[b_ofs.y * block_width + b_ofs.x];
      
      int eo_idx = EO_IDX(a, b, c);
      int eo_cat = g_sao_eo_idx_to_eo_category[eo_idx];

      c_data[0] += sao->offsets[eo_cat];
    }
  }
}

void sao_reconstruct(picture *pic, unsigned x_ctb, unsigned y_ctb, 
                     const sao_info *sao_luma, const sao_info *sao_chroma)
{
  pixel rec_y[LCU_LUMA_SIZE];
  pixel *y_recdata = &pic->y_recdata[CU_TO_PIXEL(x_ctb, y_ctb, 0, pic->width)];
  // TODO: sao chroma reconstruct

  // Data to tmp buffer.
  picture_blit_pixels(y_recdata, rec_y, LCU_WIDTH, LCU_WIDTH, pic->width, LCU_WIDTH);

  sao_reconstruct_color(rec_y, sao_luma, COLOR_Y);
  //sao_reconstruct_color(rec_u, sao_chroma, COLOR_U);
  //sao_reconstruct_color(rec_v, sao_chroma, COLOR_V);
  
  // Copy reconstructed block from tmp buffer to rec image.
  picture_blit_pixels(rec_y, y_recdata, LCU_WIDTH, LCU_WIDTH, LCU_WIDTH, pic->width);
}

void sao_search_best_mode(const pixel *data, const pixel *recdata, 
                          unsigned block_width, unsigned buf_size, unsigned buf_cnt,
                          sao_info *sao_out)
{
  sao_eo_class edge_class;
  // This array is used to calculate the mean offset used to minimize distortion.
  int cat_sum_cnt[2][NUM_SAO_EDGE_CATEGORIES];
  memset(cat_sum_cnt, 0, 2 * NUM_SAO_EDGE_CATEGORIES);

  sao_out->ddistortion = 0;

  for (edge_class = SAO_EO0; edge_class <= SAO_EO3; ++edge_class) {
    int edge_offset[NUM_SAO_EDGE_CATEGORIES];
    int sum_ddistortion = 0;
    sao_eo_cat edge_cat;

    // Call calc_sao_edge_dir once for luma and twice for chroma.
    while (buf_cnt--) {
      calc_sao_edge_dir(data, recdata, edge_class, block_width, cat_sum_cnt);
      data += buf_size;
      recdata += buf_size;
    }
    
    for (edge_cat = SAO_EO_CAT1; edge_cat <= SAO_EO_CAT4; ++edge_cat) {
      int cat_sum = cat_sum_cnt[0][edge_cat];
      int cat_cnt = cat_sum_cnt[1][edge_cat];
      
      // The optimum offset can be calculated by getting the minima of the
      // fast ddistortion estimation formula. The minima is the mean error
      // and we round that to the nearest integer.
      int offset = (cat_sum + (cat_cnt >> 1)) / cat_cnt;
      edge_offset[edge_cat] = offset;
      // The ddistortion is amount by which the SSE of data changes. It should
      // be negative for all categories, if offset was chosen correctly.
      // ddistortion = N * h^2 - 2 * h * E, where N is the number of samples 
      // and E is the sum of errors.
      // It basically says that all pixels that are not improved by offset
      // increase increase SSE by h^2 and all pixels that are improved by
      // offset decrease SSE by h*E.
      sum_ddistortion += cat_cnt * offset * offset - 2 * offset * cat_sum;
    }
    // SAO is not applied for category 0.
    edge_offset[SAO_EO_CAT0] = 0;

    // Choose the offset class that offers the least error after offset.
    if (sum_ddistortion < sao_out->ddistortion) {
      sao_out->eo_class = edge_class;
      sao_out->ddistortion = sum_ddistortion;
      memcpy(sao_out->offsets, edge_offset, NUM_SAO_EDGE_CATEGORIES);
    }
  }
}

sao_info sao_search_chroma(const picture *pic, unsigned x_ctb, unsigned y_ctb)
{
  sao_info sao;
  sao.merge_left_flag = 0;
  sao.merge_up_flag = 0;
  sao.type = SAO_TYPE_NONE;
  return sao;
}

sao_info sao_search_luma(const picture *pic, unsigned x_ctb, unsigned y_ctb)
{
  // These buffers are needed only until we switch to a LCU based data
  // structure for pixels. Then we can give pointers directly to that structure
  // without making copies.
  // It's 2-dimensional because sao_search_best_mode takes arguments as arrays.
  pixel orig_y[LCU_LUMA_SIZE];
  pixel rec_y[LCU_LUMA_SIZE];
  pixel *y_data = &pic->y_data[CU_TO_PIXEL(x_ctb, y_ctb, 0, pic->width)];
  pixel *y_recdata = &pic->y_recdata[CU_TO_PIXEL(x_ctb, y_ctb, 0, pic->width)];
  sao_info sao_params;

  // Fill temporary buffers with picture data.
  picture_blit_pixels(y_data, orig_y, LCU_WIDTH, LCU_WIDTH, pic->width, LCU_WIDTH);
  picture_blit_pixels(y_recdata, rec_y, LCU_WIDTH, LCU_WIDTH, pic->width, LCU_WIDTH);

  sao_search_best_mode(orig_y, rec_y, LCU_WIDTH, LCU_LUMA_SIZE, 1, &sao_params);

  return sao_params;
}

void encode_sao_offsets(encoder_control *encoder, sao_info *sao)
{
  int i;

  for (i = 0; i < NUM_SAO_OFFSETS; ++i) {
    CABAC_BIN(&cabac, sao->offsets[i] > 0 ? 0 : 1, "sao_offset_sign");
  }

  if (sao->type == SAO_TYPE_EDGE) {
    for (i = 0; i < NUM_SAO_OFFSETS; ++i) {
      if (sao->offsets[i] != 0) {
        // For edge SAO positive sign is encoded as 0.
        CABAC_BIN(&cabac, sao->offsets[i] > 0 ? 0 : 1, "sao_offset_sign");
        // TODO: CABAC_BIN sao_band_position[color_i]
      } else {
        // TODO: CABAC_BIN sao_eo_class[color_i]
      }
    }
  }
}

void encode_sao_color(encoder_control *encoder, sao_info *sao, color_index color_i)
{
  picture *pic = encoder->in.cur_pic;

  // Skip colors with no SAO.
  if (color_i == COLOR_Y && !pic->slice_sao_luma_flag) {
    return;
  } else if (!pic->slice_sao_chroma_flag) {
    return;
  }

  cabac.ctx = &g_sao_type_idx_model;
  if (color_i == COLOR_Y) {
    CABAC_BIN(&cabac, sao->type, "sao_type_idx_luma");
  } else if (color_i == COLOR_U) {
    // SAO type is only coded for the first chroma.
    CABAC_BIN(&cabac, sao->type, "sao_type_idx_chroma");
  }

  if (sao->type != SAO_TYPE_NONE) {
    encode_sao_offsets(encoder, 0);
  }
}

void encode_sao_merge_flags(encoder_control *encoder, sao_info *sao,
                            unsigned x_ctb, unsigned y_ctb)
{
  // SAO merge flags are not present if merge candidate is not in the same
  // slice AND tile, but there isn't any such segmentation right now.
  assert(!USE_SLICES && !USE_TILES);

  // SAO merge flags are not present for the first row and column.
  if (x_ctb > 0) {
    cabac.ctx = &g_sao_merge_flag_model;
    CABAC_BIN(&cabac, sao->merge_left_flag ? 1 : 0, "sao_merge_left_flag");
  }
  if (y_ctb > 0 && !sao->merge_left_flag) {
    cabac.ctx = &g_sao_merge_flag_model;
    CABAC_BIN(&cabac, sao->merge_up_flag ? 1 : 0, "sao_merge_up_flag");
  }
}

/**
 * \brief Stub that encodes all LCU's as none type.
 */
void encode_sao(encoder_control *encoder, unsigned x_lcu, uint16_t y_lcu,
                sao_info *sao_luma, sao_info *sao_chroma)
{
  unsigned sao_type[3] = {SAO_TYPE_NONE, SAO_TYPE_NONE, SAO_TYPE_NONE};
  picture *pic = encoder->in.cur_pic;
  
  // TODO: transmit merge flags outside sao_info
  encode_sao_merge_flags(encoder, sao_luma, x_lcu, y_lcu);

  // If SAO is merged, nothing else needs to be coded.
  if (!sao_luma->merge_left_flag && !sao_luma->merge_up_flag) {
    encode_sao_color(encoder, sao_luma, COLOR_Y);
    encode_sao_color(encoder, sao_chroma, COLOR_U);
    encode_sao_color(encoder, sao_chroma, COLOR_V);
  }
}

void encode_slice_data(encoder_control* encoder)
{
  uint16_t x_ctb, y_ctb;

  scalinglist_process();
  init_contexts(encoder,encoder->in.cur_pic->slicetype);

  // Loop through every LCU in the slice
  for (y_ctb = 0; y_ctb < encoder->in.height_in_lcu; y_ctb++) {
    uint8_t last_cu_y = (y_ctb == (encoder->in.height_in_lcu - 1)) ? 1 : 0;

    for (x_ctb = 0; x_ctb < encoder->in.width_in_lcu; x_ctb++) {
      uint8_t last_cu_x = (x_ctb == (encoder->in.width_in_lcu - 1)) ? 1 : 0;
      uint8_t depth = 0;

      if (encoder->sao_enable) {
        sao_info sao_luma = sao_search_luma(encoder->in.cur_pic, x_ctb, y_ctb);
        sao_info sao_chroma = sao_search_chroma(encoder->in.cur_pic, x_ctb, y_ctb);
        
        // sao_do_merge(encoder, x_ctb, y_ctb, sao_luma, sao_chroma);
        // sao_do_rdo(encoder, x_ctb, y_ctb, sao_luma, sao_chroma);
        encode_sao(encoder, x_ctb, y_ctb, &sao_luma, &sao_chroma);
      }

      // Recursive function for looping through all the sub-blocks
      encode_coding_tree(encoder, x_ctb << MAX_DEPTH, y_ctb << MAX_DEPTH, depth);

      // signal Terminating bit
      if (!last_cu_x || !last_cu_y) {
        cabac_encode_bin_trm(&cabac, 0);
      }
    }
  }
}

void encode_coding_tree(encoder_control *encoder, uint16_t x_ctb,
                        uint16_t y_ctb, uint8_t depth)
{ 
  cu_info *cur_cu = &encoder->in.cur_pic->cu_array[MAX_DEPTH][x_ctb + y_ctb * (encoder->in.width_in_lcu << MAX_DEPTH)];
  uint8_t split_flag = GET_SPLITDATA(cur_cu, depth);
  uint8_t split_model = 0;

  // Check for slice border
  uint8_t border_x = ((encoder->in.width) < (x_ctb * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth))) ? 1 : 0;
  uint8_t border_y = ((encoder->in.height) < (y_ctb * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth))) ? 1 : 0;
  uint8_t border_split_x = ((encoder->in.width)  < ((x_ctb + 1) * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> (depth + 1)))) ? 0 : 1;
  uint8_t border_split_y = ((encoder->in.height) < ((y_ctb + 1) * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> (depth + 1)))) ? 0 : 1;
  uint8_t border = border_x | border_y; /*!< are we in any border CU */
  

  // When not in MAX_DEPTH, insert split flag and split the blocks if needed
  if (depth != MAX_DEPTH) {
    // Implisit split flag when on border
    if (!border) {
      // Get left and top block split_flags and if they are present and true, increase model number
      if (x_ctb > 0 && GET_SPLITDATA(&(encoder->in.cur_pic->cu_array[MAX_DEPTH][x_ctb - 1 + y_ctb * (encoder->in.width_in_lcu << MAX_DEPTH)]), depth) == 1) {
        split_model++;
      }

      if (y_ctb > 0 && GET_SPLITDATA(&(encoder->in.cur_pic->cu_array[MAX_DEPTH][x_ctb + (y_ctb - 1) * (encoder->in.width_in_lcu << MAX_DEPTH)]), depth) == 1) {
        split_model++;
      }

      cabac.ctx = &g_split_flag_model[split_model];
      CABAC_BIN(&cabac, split_flag, "SplitFlag");
    }

    if (split_flag || border) {
      // Split blocks and remember to change x and y block positions
      uint8_t change = 1<<(MAX_DEPTH-1-depth);
      encode_coding_tree(encoder, x_ctb, y_ctb, depth + 1); // x,y

      // TODO: fix when other half of the block would not be completely over the border 
      if (!border_x || border_split_x) { 
        encode_coding_tree(encoder, x_ctb + change, y_ctb, depth + 1);
      }
      if (!border_y || border_split_y) {
        encode_coding_tree(encoder, x_ctb, y_ctb + change, depth + 1);
      }
      if (!border || (border_split_x && border_split_y)) {
        encode_coding_tree(encoder, x_ctb + change, y_ctb + change, depth + 1);
      }
      
      return;
    }
  }
  
  // Encode skip flag
  if (encoder->in.cur_pic->slicetype != SLICE_I) {
    int8_t ctx_skip = 0;
    // uiCtxSkip = aboveskipped + leftskipped;
    cabac.ctx = &g_cu_skip_flag_model[ctx_skip];
    CABAC_BIN(&cabac, (cur_cu->type == CU_SKIP) ? 1 : 0, "SkipFlag");
  }

  // IF SKIP
  if (cur_cu->type == CU_SKIP) {
    // Encode merge index
    //TODO: calculate/fetch merge candidates
    int16_t unary_idx = 0; //pcCU->getMergeIndex( uiAbsPartIdx );
    int16_t num_cand = 0; //pcCU->getSlice()->getMaxNumMergeCand();
    int32_t ui;

    if (num_cand > 1) {
      for (ui = 0; ui < num_cand - 1; ui++) {
        int32_t symbol = (ui == unary_idx) ? 0 : 1;

        if (ui == 0) {
          cabac.ctx = &g_cu_merge_idx_ext_model;
          CABAC_BIN(&cabac, symbol, "MergeIndex");
        } else {
          CABAC_BIN_EP(&cabac,symbol,"MergeIndex");
        }

        if (symbol == 0) {
          break;
        }
      }
    }

    return;
  }

  // ENDIF SKIP

  // Prediction mode
  if (encoder->in.cur_pic->slicetype != SLICE_I) {
    cabac.ctx = &g_cu_pred_mode_model;
    CABAC_BIN(&cabac, (cur_cu->type == CU_INTRA) ? 1 : 0, "PredMode");
  }

  // Signal PartSize on max depth
  if (depth == MAX_DEPTH || cur_cu->type != CU_INTRA) {
    // TODO: Handle inter sizes other than 2Nx2N
    cabac.ctx = &g_part_size_model[0];
    CABAC_BIN(&cabac, 1, "PartSize");
    // TODO: add AMP modes
  }
    
  //end partsize
  if (cur_cu->type == CU_INTER) {
    // FOR each part
    // Mergeflag
    uint8_t merge_flag = 0;
        cabac.ctx = &g_cu_merge_flag_ext_model;
    CABAC_BIN(&cabac, merge_flag, "MergeFlag");

    if (merge_flag) { //merge
      // MergeIndex
      int16_t unary_idx = 0; //pcCU->getMergeIndex( uiAbsPartIdx );
      int16_t num_cand  = 0; //pcCU->getSlice()->getMaxNumMergeCand();
          int32_t ui;

      if (num_cand > 1) {
        for (ui = 0; ui < num_cand - 1; ui++) {
          int32_t symbol = (ui == unary_idx) ? 0 : 1;

          if (ui == 0) {
                cabac.ctx = &g_cu_merge_idx_ext_model;
                CABAC_BIN(&cabac, symbol, "MergeIndex");
          } else {
                CABAC_BIN_EP(&cabac,symbol,"MergeIndex");
              }

          if (symbol == 0) {
                break;
              }
            }
          }
    } else {
      uint32_t ref_list_idx;
          int16_t mv_cand[2][2];

          /*
          // Void TEncSbac::codeInterDir( TComDataCU* pcCU, UInt uiAbsPartIdx )
          if(encoder->in.cur_pic->slicetype == SLICE_B)
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
            //if(encoder->ref_idx_num[uiRefListIdx] > 0)
            {
          if (cur_cu->inter.mv_dir & (1 << ref_list_idx)) {
            if (0) { //encoder->ref_idx_num[uiRefListIdx] != 1)//NumRefIdx != 1)
              // parseRefFrmIdx
              int32_t ref_frame = cur_cu->inter.mv_ref;
                  
                  cabac.ctx = &g_cu_ref_pic_model[0];
              CABAC_BIN(&cabac, (ref_frame == 0) ? 0 : 1, "ref_frame_flag");
    
              if (ref_frame > 0) {
                uint32_t i;
                uint32_t ref_num = encoder->ref_idx_num[ref_list_idx] - 2;

                    cabac.ctx = &g_cu_ref_pic_model[1];
                ref_frame--;

                for (i = 0; i < ref_num; ++i) {
                  const uint32_t symbol = (i == ref_frame) ? 0 : 1;

                  if (i == 0) {
                    CABAC_BIN(&cabac, symbol, "ref_frame_flag2");
                  } else {
                    CABAC_BIN_EP(&cabac, symbol, "ref_frame_flag2");
                      }

                  if (symbol == 0) {
                        break;
                      }
                    }
                  }
                }

            // Get MV candidates
            inter_get_mv_cand(encoder, x_ctb, y_ctb, depth, mv_cand);

            // Select better candidate
            cur_cu->inter.mv_ref = 0; // Default to candidate 0

            // Only check when candidates are different
                if (mv_cand[0][0] != mv_cand[1][0] || mv_cand[0][1] != mv_cand[1][1]) {
              uint16_t cand_1_diff = abs(cur_cu->inter.mv[0] - mv_cand[0][0]) + abs(
                                       cur_cu->inter.mv[1] - mv_cand[0][1]);
              uint16_t cand_2_diff = abs(cur_cu->inter.mv[0] - mv_cand[1][0]) + abs(
                                       cur_cu->inter.mv[1] - mv_cand[1][1]);

              // Select candidate 1 if it's closer
                  if (cand_2_diff < cand_1_diff) {
                cur_cu->inter.mv_ref = 1;
                  }
                }

            if (!(/*pcCU->getSlice()->getMvdL1ZeroFlag() &&*/ encoder->ref_list == REF_PIC_LIST_1 && cur_cu->inter.mv_dir == 3)) {
              const int32_t mvd_hor = cur_cu->inter.mv[0] - mv_cand[cur_cu->inter.mv_ref][0];
              const int32_t mvd_ver = cur_cu->inter.mv[1] - mv_cand[cur_cu->inter.mv_ref][1];
              const int8_t hor_abs_gr0 = mvd_hor != 0;
              const int8_t ver_abs_gr0 = mvd_ver != 0;
                    const uint32_t mvd_hor_abs = abs(mvd_hor);
                    const uint32_t mvd_ver_abs = abs(mvd_ver);

                    cabac.ctx = &g_cu_mvd_model[0];
                    CABAC_BIN(&cabac, (mvd_hor!=0)?1:0, "abs_mvd_greater0_flag_hor");
                    CABAC_BIN(&cabac, (mvd_ver!=0)?1:0, "abs_mvd_greater0_flag_ver");

                    cabac.ctx = &g_cu_mvd_model[1];

              if (hor_abs_gr0) {
                      CABAC_BIN(&cabac, (mvd_hor_abs>1)?1:0, "abs_mvd_greater1_flag_hor");
                    }

              if (ver_abs_gr0) {
                      CABAC_BIN(&cabac, (mvd_ver_abs>1)?1:0, "abs_mvd_greater1_flag_ver");
                    }

              if (hor_abs_gr0) {
                if (mvd_hor_abs > 1) {
                        cabac_write_ep_ex_golomb(&cabac,mvd_hor_abs-2, 1);
                      }

                      CABAC_BIN_EP(&cabac, (mvd_hor>0)?0:1, "mvd_sign_flag_hor");
                    }

              if (ver_abs_gr0) {
                if (mvd_ver_abs > 1) {
                        cabac_write_ep_ex_golomb(&cabac,mvd_ver_abs-2, 1);
                      }

                      CABAC_BIN_EP(&cabac, (mvd_ver>0)?0:1, "mvd_sign_flag_ver");
                    }

              // Inter reconstruction
              inter_recon(encoder->ref->pics[0], x_ctb * CU_MIN_SIZE_PIXELS,
                          y_ctb * CU_MIN_SIZE_PIXELS, LCU_WIDTH >> depth, cur_cu->inter.mv,
                          encoder->in.cur_pic);

              // Mark this block as "coded" (can be used for predictions..)
              picture_set_block_coded(encoder->in.cur_pic, x_ctb, y_ctb, depth, 1);
                }

            // Signal which candidate MV to use
            cabac_write_unary_max_symbol(&cabac, g_mvp_idx_model, cur_cu->inter.mv_ref, 1,
                                         AMVP_MAX_NUM_CANDS - 1);
              }
            }
          }

          
          if (1) {
            pixel *base_y  = &encoder->in.cur_pic->y_data[x_ctb*(LCU_WIDTH>>(MAX_DEPTH))   + (y_ctb*(LCU_WIDTH>>(MAX_DEPTH)))  *encoder->in.width];
            pixel *base_u = &encoder->in.cur_pic->u_data[x_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (y_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
            pixel *base_v = &encoder->in.cur_pic->v_data[x_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (y_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
            uint32_t width = LCU_WIDTH>>depth;

            /* INTRAPREDICTION VARIABLES */
            int16_t pred[LCU_WIDTH*LCU_WIDTH+1];
            int16_t predU[LCU_WIDTH*LCU_WIDTH>>2];
            int16_t predV[LCU_WIDTH*LCU_WIDTH>>2];

            pixel *recbase_y = &encoder->in.cur_pic->y_recdata[x_ctb*(LCU_WIDTH>>(MAX_DEPTH))   + (y_ctb*(LCU_WIDTH>>(MAX_DEPTH)))  *encoder->in.width];
            pixel *recbase_u = &encoder->in.cur_pic->u_recdata[x_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (y_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
            pixel *recbase_v = &encoder->in.cur_pic->v_recdata[x_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (y_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];

            /* TODO: dynamic memory allocation */
            int16_t coeff_y[LCU_WIDTH*LCU_WIDTH*2];
            int16_t coeff_u[LCU_WIDTH*LCU_WIDTH>>1];
            int16_t coeff_v[LCU_WIDTH*LCU_WIDTH>>1];
            int8_t residual = 0;

            /* Initialize helper structure for transform */
            transform_info ti;
            memset(&ti, 0, sizeof(transform_info));

            ti.x_ctb = x_ctb; ti.y_ctb = y_ctb;

            /* Base pointers */
            ti.base =  base_y; ti.base_u = base_u; ti.base_v = base_v;
            ti.base_stride = encoder->in.width;

            // Prediction pointers
            ti.pred =  pred; ti.pred_u = predU; ti.pred_v = predV;
            ti.pred_stride = (LCU_WIDTH>>depth);

            // Reconstruction pointers
            ti.recbase = recbase_y; ti.recbase_u = recbase_u; ti.recbase_v = recbase_v;
            ti.recbase_stride = encoder->in.width;

            // Coeff pointers
            ti.coeff[0] = coeff_y; ti.coeff[1] = coeff_u; ti.coeff[2] = coeff_v;
            ti.block_type = CU_INTER;
       
            // Handle transforms, quant and reconstruction
            ti.idx = 0;
            encode_transform_tree(encoder,&ti, depth);

            // Coded block pattern
            ti.cb_top[0] = (ti.cb[0] & 0x1 || ti.cb[1] & 0x1 || ti.cb[2] & 0x1 || ti.cb[3] & 0x1)?1:0;
            ti.cb_top[1] = (ti.cb[0] & 0x2 || ti.cb[1] & 0x2 || ti.cb[2] & 0x2 || ti.cb[3] & 0x2)?1:0;
            ti.cb_top[2] = (ti.cb[0] & 0x4 || ti.cb[1] & 0x4 || ti.cb[2] & 0x4 || ti.cb[3] & 0x4)?1:0;

            residual = ti.cb_top[0] | ti.cb_top[1] | ti.cb_top[2];
            if(depth == 0)  {
              picture_set_block_residual(encoder->in.cur_pic,x_ctb    ,y_ctb    ,depth+1,ti.cb[0] & 0x1);
              picture_set_block_residual(encoder->in.cur_pic,x_ctb + 4,y_ctb    ,depth+1,ti.cb[1] & 0x1);
              picture_set_block_residual(encoder->in.cur_pic,x_ctb    ,y_ctb + 4,depth+1,ti.cb[2] & 0x1);
              picture_set_block_residual(encoder->in.cur_pic,x_ctb + 4,y_ctb + 4,depth+1,ti.cb[3] & 0x1);
            } else  {
              picture_set_block_residual(encoder->in.cur_pic,x_ctb,y_ctb,depth,ti.cb_top[0]);
            }
            

            cabac.ctx = &g_cu_qt_root_cbf_model;
            CABAC_BIN(&cabac, residual, "rqt_root_cbf");
            // Code (possible) coeffs to bitstream
            ti.idx = 0;
            if(residual) {
              encode_transform_coeff(encoder, &ti,depth, 0);
            }
          }
        }

    // END for each part
  } else if (cur_cu->type == CU_INTRA) {
    uint8_t intra_pred_mode = cur_cu->intra.mode;
    uint8_t intra_pred_mode_chroma = 36; // 36 = Chroma derived from luma
    int8_t intra_preds[3] = { -1, -1, -1};
    int8_t mpm_preds = -1;
    int i;
    uint32_t flag;
    pixel *base_y = &encoder->in.cur_pic->y_data[x_ctb * (LCU_WIDTH >> (MAX_DEPTH))     + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH)))     * encoder->in.width];
    pixel *base_u = &encoder->in.cur_pic->u_data[x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1))) * (encoder->in.width >> 1)];
    pixel *base_v = &encoder->in.cur_pic->v_data[x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1))) * (encoder->in.width >> 1)];
    uint32_t width = LCU_WIDTH>>depth;

    // INTRAPREDICTION VARIABLES
    int16_t pred_y[LCU_WIDTH * LCU_WIDTH + 1];
    int16_t pred_u[LCU_WIDTH * LCU_WIDTH >> 2];
    int16_t pred_v[LCU_WIDTH * LCU_WIDTH >> 2];

    pixel *recbase_y = &encoder->in.cur_pic->y_recdata[x_ctb * (LCU_WIDTH >> (MAX_DEPTH))     + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH)))     * encoder->in.width];
    pixel *recbase_u = &encoder->in.cur_pic->u_recdata[x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1))) * (encoder->in.width >> 1)];
    pixel *recbase_v = &encoder->in.cur_pic->v_recdata[x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1))) * (encoder->in.width >> 1)];

    // SEARCH BEST INTRA MODE (AGAIN)
    int16_t rec[(LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)];
    int16_t *rec_shift = &rec[(LCU_WIDTH >> (depth)) * 2 + 8 + 1];
    intra_build_reference_border(encoder->in.cur_pic, x_ctb, y_ctb,
                                 (LCU_WIDTH >> (depth)) * 2 + 8, rec, 
                                 (LCU_WIDTH >> (depth)) * 2 + 8, 0);
    cur_cu->intra.mode = (int8_t)intra_prediction(encoder->in.cur_pic->y_data,
                                                  encoder->in.width, 
                                                  rec_shift, 
                                                  (LCU_WIDTH >> (depth)) * 2 + 8,
                                                  x_ctb * (LCU_WIDTH >> (MAX_DEPTH)), 
                                                  y_ctb * (LCU_WIDTH >> (MAX_DEPTH)), 
                                                  width, pred_y, width, 
                                                  &cur_cu->intra.cost);
    intra_pred_mode = cur_cu->intra.mode;
    intra_set_block_mode(encoder->in.cur_pic, x_ctb, y_ctb, depth,
                         intra_pred_mode);
      
    #if ENABLE_PCM == 1
    // Code must start after variable initialization
    cabac_encode_bin_trm(&cabac, 0); // IPCMFlag == 0
    #endif
      
    // PREDINFO CODING
    // If intra prediction mode is found from the predictors,
    // it can be signaled with two EP's. Otherwise we can send
    // 5 EP bins with the full predmode
    intra_get_dir_luma_predictor(encoder->in.cur_pic, x_ctb, y_ctb, depth,
                                 intra_preds);
      
    for (i = 0; i < 3; i++) {
      if (intra_preds[i] == intra_pred_mode) {
        mpm_preds = i;
          break;
        }
      }

    // For each part {
    flag = (mpm_preds == -1) ? 0 : 1;
      cabac.ctx = &g_intra_mode_model;
      CABAC_BIN(&cabac,flag,"IntraPred");
    // } End for each part

    // Intrapredmode signaling
    // If found from predictors, we can simplify signaling
    if (flag) {
      flag = (mpm_preds == 0) ? 0 : 1;
        CABAC_BIN_EP(&cabac, flag, "intraPredMode");

      if (mpm_preds != 0) {
        flag = (mpm_preds == 1) ? 0 : 1;
          CABAC_BIN_EP(&cabac, flag, "intraPredMode");
        }
    } else { 
      // we signal the "full" predmode
      int32_t intra_pred_mode_temp = intra_pred_mode;

      if (intra_preds[0] > intra_preds[1]) {
        SWAP(intra_preds[0], intra_preds[1], int8_t);
      }

      if (intra_preds[0] > intra_preds[2]) {
        SWAP(intra_preds[0], intra_preds[2], int8_t);
        }

      if (intra_preds[1] > intra_preds[2]) {
        SWAP(intra_preds[1], intra_preds[2], int8_t);
        }

      for (i = 2; i >= 0; i--) {
        intra_pred_mode_temp = intra_pred_mode_temp > intra_preds[i] ?
                               intra_pred_mode_temp - 1 : intra_pred_mode_temp;
        }

      CABAC_BINS_EP(&cabac, intra_pred_mode_temp, 5, "intraPredMode");
        }

    // If we have chroma, signal it
    if (encoder->in.video_format != FORMAT_400) {
      // Chroma intra prediction
      cabac.ctx = &g_chroma_pred_model[0];
      CABAC_BIN(&cabac, ((intra_pred_mode_chroma != 36) ? 1 : 0), "IntraPredChroma");

      // If not copied from luma, signal it
      if (intra_pred_mode_chroma != 36) {
        int8_t intra_pred_mode_chroma_temp = intra_pred_mode_chroma;
        // Default chroma predictors
        uint32_t allowed_chroma_dir[5] = { 0, 26, 10, 1, 36 };
          
        // If intra is the same as one of the default predictors, replace it
        for (i = 0; i < 4; i++) {
          if (intra_pred_mode == allowed_chroma_dir[i]) {
            allowed_chroma_dir[i] = 34; /* VER+8 mode */
              break;
            }
          }

        for (i = 0; i < 4; i++) {
          if (intra_pred_mode_chroma_temp == allowed_chroma_dir[i]) {
            intra_pred_mode_chroma_temp = i;
            break;
          }
        }

        CABAC_BINS_EP(&cabac, intra_pred_mode_chroma_temp, 2, "intraPredModeChroma");
      }
    }

    // END OF PREDINFO CODING
    
    // Coeff
    // Transform tree
      {
        // TODO: dynamic memory allocation
        int16_t coeff_y[LCU_WIDTH * LCU_WIDTH * 2];
        int16_t coeff_u[LCU_WIDTH * LCU_WIDTH >> 1];
        int16_t coeff_v[LCU_WIDTH * LCU_WIDTH >> 1];

        // Initialize helper structure for transform
        transform_info ti;
        memset(&ti, 0, sizeof(transform_info));

        ti.x_ctb = x_ctb; ti.y_ctb = y_ctb;

        // Base pointers
        ti.base =  base_y; ti.base_u = base_u; ti.base_v = base_v;
        ti.base_stride = encoder->in.width;

        // Prediction pointers
        ti.pred =  pred_y; ti.pred_u = pred_u; ti.pred_v = pred_v;
        ti.pred_stride = (LCU_WIDTH>>depth);

        // Reconstruction pointers
        ti.recbase = recbase_y; ti.recbase_u = recbase_u; ti.recbase_v = recbase_v;
        ti.recbase_stride = encoder->in.width;

        // Coeff pointers
        ti.coeff[0] = coeff_y; ti.coeff[1] = coeff_u; ti.coeff[2] = coeff_v;

        // Prediction info
        ti.intra_pred_mode = intra_pred_mode;
        ti.intra_pred_mode_chroma = intra_pred_mode_chroma;
        
        // Handle transforms, quant and reconstruction
        ti.idx = 0;
        ti.block_type = CU_INTRA;
        encode_transform_tree(encoder,&ti, depth);

        // Coded block pattern
        ti.cb_top[0] = (ti.cb[0] & 0x1 || ti.cb[1] & 0x1 || ti.cb[2] & 0x1
                        || ti.cb[3] & 0x1) ? 1 : 0;
        ti.cb_top[1] = (ti.cb[0] & 0x2 || ti.cb[1] & 0x2 || ti.cb[2] & 0x2
                        || ti.cb[3] & 0x2) ? 1 : 0;
        ti.cb_top[2] = (ti.cb[0] & 0x4 || ti.cb[1] & 0x4 || ti.cb[2] & 0x4
                        || ti.cb[3] & 0x4) ? 1 : 0;
        
        // Code (possible) coeffs to bitstream
        ti.idx = 0;
        encode_transform_coeff(encoder, &ti,depth, 0);
      }
    // end Transform tree
    // end Coeff

    }

    #if ENABLE_PCM == 1
  // Code IPCM block
  if (cur_cu->type == CU_PCM) {
    cabac_encode_bin_trm(&cabac, 1); // IPCMFlag == 1
      cabac_finish(&cabac);
      bitstream_align(cabac.stream);
    // PCM sample
      {
      unsigned y, x;

      pixel *base_y = &encoder->in.cur_pic->y_data[x_ctb * (LCU_WIDTH >> (MAX_DEPTH))    + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH))) * encoder->in.width];
      pixel *base_u = &encoder->in.cur_pic->u_data[(x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1))) * encoder->in.width / 2)];
      pixel *base_v = &encoder->in.cur_pic->v_data[(x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)) + (y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1))) * encoder->in.width / 2)];

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
      cabac_start(&cabac);

  } // end Code IPCM block

    #endif /* END ENABLE_PCM */
  else { /* Should not happend */
      printf("UNHANDLED TYPE!\r\n");
      exit(1);
    }

   /* end prediction unit */
  /* end coding_unit */
  
}

void encode_transform_tree(encoder_control *encoder, transform_info *ti,
                           uint8_t depth)
{
  // we have 64>>depth transform size
  int x,y,i;
  int32_t width = LCU_WIDTH>>depth;

  if (depth == 0) { // Split 64x64
    // Prepare for multi-level splitting
    ti->split[ti->idx] = 1<<depth;
  }

  // Split transform and increase depth
  if (ti->split[ti->idx] & (1 << depth)) {
    uint8_t change = 1<<(MAX_DEPTH-1-depth);    
    ti->idx = 0; encode_transform_tree(encoder,ti,depth+1);
    ti->x_ctb += change;
    ti->idx = 1; encode_transform_tree(encoder,ti,depth+1);
    ti->x_ctb -= change; ti->y_ctb += change;
    ti->idx = 2; encode_transform_tree(encoder,ti,depth+1);
    ti->x_ctb += change;
    ti->idx = 3; encode_transform_tree(encoder,ti,depth+1);
    return;
  }
  
  {
    uint8_t cb_y = 0, cb_u = 0, cb_v = 0;
    int32_t coeff_fourth = ((LCU_WIDTH>>(depth))*(LCU_WIDTH>>(depth)))+1;

    int32_t base_stride    = ti->base_stride;
    int32_t recbase_stride = ti->recbase_stride;
    int32_t pred_stride    = ti->pred_stride;

    int32_t recbase_offset[4] = {
      0, width, ti->recbase_stride * width,
      ti->recbase_stride * width + width
    };
    int32_t base_offset[4]    = {
      0, width, ti->base_stride * width,
      ti->base_stride * width + width 
    };
    int32_t pred_offset[4]    = {
      0, width, ti->pred_stride * width, 
      ti->pred_stride * width + width
    };
    
    int32_t recbase_offset_c[4] = {
      0, width >> 1, (ti->recbase_stride >> 1) * (width >> 1),
      (ti->recbase_stride >> 1) *(width >> 1) + (width >> 1) 
    };
    int32_t base_offset_c[4]    = {
      0, width >> 1, (ti->base_stride >> 1) * (width >> 1),
      (ti->base_stride >> 1) * (width >> 1)   + (width >> 1) 
    };
    int32_t pred_offset_c[4]    = {
      0, width >> 1, (ti->pred_stride >> 1) * (width >> 1), 
      (ti->pred_stride >> 1) * (width >> 1)   + (width >> 1) 
    };
    
    pixel *base_y    = &ti->base[base_offset[ti->idx]];
    pixel *base_u    = &ti->base_u[base_offset_c[ti->idx]];
    pixel *base_v    = &ti->base_v[base_offset_c[ti->idx]];
    
    pixel *recbase_y = &ti->recbase[recbase_offset[ti->idx]];
    pixel *recbase_u = &ti->recbase_u[recbase_offset_c[ti->idx]];
    pixel *recbase_v = &ti->recbase_v[recbase_offset_c[ti->idx]];
    
    int16_t *pred_y    = &ti->pred[pred_offset[ti->idx]];
    int16_t *pred_u    = &ti->pred_u[pred_offset_c[ti->idx]];
    int16_t *pred_v    = &ti->pred_v[pred_offset_c[ti->idx]];
    
    int16_t *coeff_y   = &ti->coeff[0][ti->idx * coeff_fourth];
    int16_t *coeff_u   = &ti->coeff[1][ti->idx * coeff_fourth >> 1];
    int16_t *coeff_v   = &ti->coeff[2][ti->idx * coeff_fourth >> 1];

    // Quant and transform here...
    int16_t block[LCU_WIDTH*LCU_WIDTH>>2];
    int16_t pre_quant_coeff[LCU_WIDTH*LCU_WIDTH>>2];

    // INTRA PREDICTION
    // TODO: split to a function!
    int16_t rec[(LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)];
    int16_t *rec_shift  = &rec[(LCU_WIDTH >> (depth)) * 2 + 8 + 1];
    int16_t *rec_shift_u = &rec[(LCU_WIDTH >> (depth + 1)) * 2 + 8 + 1];

    uint32_t ac_sum = 0;
    uint32_t ctx_idx;
    uint32_t scan_idx_luma = SCAN_DIAG;
    uint32_t scan_idx_chroma = SCAN_DIAG;
    uint8_t dir_mode;
    #if OPTIMIZATION_SKIP_RESIDUAL_ON_THRESHOLD
    uint32_t residual_sum = 0;
    #endif

    switch (width) {
      case  2: ctx_idx = 6; break;
      case  4: ctx_idx = 5; break;
      case  8: ctx_idx = 4; break;
      case 16: ctx_idx = 3; break;
      case 32: ctx_idx = 2; break;
      case 64: ctx_idx = 1; break;
      default: ctx_idx = 0; break;
    }

    if(ti->block_type == CU_INTRA)
    {
      //if multiple scans supported for transform size
      if (ctx_idx > 3 && ctx_idx < 6) {
        scan_idx_luma = abs((int32_t) ti->intra_pred_mode - 26) < 5 ? 1 : (abs((int32_t)ti->intra_pred_mode - 10) < 5 ? 2 : 0);
      }

      // Chroma scanmode
      ctx_idx++;
      dir_mode = ti->intra_pred_mode_chroma;

      if (dir_mode == 36) {
        // TODO: support NxN
        dir_mode = ti->intra_pred_mode;
      }

      if (ctx_idx > 4 && ctx_idx < 7) { // if multiple scans supported for transform size
        scan_idx_chroma = abs((int32_t) dir_mode - 26) < 5 ? 1 : (abs((int32_t)dir_mode - 10) < 5 ? 2 : 0);
      }

      // Build reconstructed block to use in prediction with extrapolated borders
      intra_build_reference_border(encoder->in.cur_pic, ti->x_ctb, ti->y_ctb,
                                   (LCU_WIDTH >> (depth)) * 2 + 8, rec, (LCU_WIDTH >> (depth)) * 2 + 8, 0);
      intra_recon(rec_shift, (LCU_WIDTH >> (depth)) * 2 + 8,
                  ti->x_ctb * (LCU_WIDTH >> (MAX_DEPTH)), ti->y_ctb * (LCU_WIDTH >> (MAX_DEPTH)),
                  width, pred_y, pred_stride, ti->intra_pred_mode, 0);

      // Filter DC-prediction
      if (ti->intra_pred_mode == 1 && width < 32) {
        intra_dc_pred_filtering(rec_shift, (LCU_WIDTH >> (depth)) * 2 + 8, pred_y,
                                width, LCU_WIDTH >> depth, LCU_WIDTH >> depth);
      }

      if (ti->intra_pred_mode_chroma != 36
          && ti->intra_pred_mode_chroma == ti->intra_pred_mode) {
          ti->intra_pred_mode_chroma = 36;
        }
    
      intra_build_reference_border(encoder->in.cur_pic, ti->x_ctb, ti->y_ctb,
                                   (LCU_WIDTH >> (depth + 1)) * 2 + 8, rec,
                                   (LCU_WIDTH >> (depth + 1)) * 2 + 8,
                                   1);
      intra_recon(rec_shift_u, 
                  (LCU_WIDTH >> (depth + 1)) * 2 + 8,
                  ti->x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)),
                  ti->y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)),
                  width >> 1,
                  pred_u,
                  pred_stride >> 1,
                  ti->intra_pred_mode_chroma != 36 ? ti->intra_pred_mode_chroma : ti->intra_pred_mode,
                  1);
      intra_build_reference_border(encoder->in.cur_pic, ti->x_ctb, ti->y_ctb,
                                   (LCU_WIDTH >> (depth + 1)) * 2 + 8,
                                   rec, (LCU_WIDTH >> (depth + 1)) * 2 + 8,
                                   2);
      intra_recon(rec_shift_u, (LCU_WIDTH >> (depth + 1)) * 2 + 8,
                  ti->x_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)),
                  ti->y_ctb * (LCU_WIDTH >> (MAX_DEPTH + 1)),
                  width >> 1,
                  pred_v,
                  pred_stride >> 1,
                  ti->intra_pred_mode_chroma != 36 ? ti->intra_pred_mode_chroma : ti->intra_pred_mode,
                  1);

    // This affects reconstruction, do after that
      picture_set_block_coded(encoder->in.cur_pic, ti->x_ctb, ti->y_ctb, depth, 1);
    } else  { // Inter mode
      for(y = 0; y < LCU_WIDTH>>depth; y++) {
        for(x = 0; x < LCU_WIDTH>>depth; x++) {
          pred_y[x+y*pred_stride]=recbase_y[x+y*base_stride];
        }
      }
      for(y = 0; y < LCU_WIDTH>>(depth+1); y++) {
        for(x = 0; x < LCU_WIDTH>>(depth+1); x++) {
          pred_u[x+y*(pred_stride>>1)]=recbase_u[x+y*(base_stride>>1)];
          pred_v[x+y*(pred_stride>>1)]=recbase_v[x+y*(base_stride>>1)];
        }
      }
    }
    // INTRA PREDICTION ENDS HERE

    // Get residual by subtracting prediction
    i = 0;
    ac_sum = 0;

    for (y = 0; y < LCU_WIDTH >> depth; y++) {
      for (x = 0; x < LCU_WIDTH >> depth; x++) {
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
    if(residual_sum < RESIDUAL_THRESHOLD/(LCU_WIDTH >> depth)) {
      memset(block, 0, sizeof(int16_t)*(LCU_WIDTH >> depth)*(LCU_WIDTH >> depth));
    }
    #endif

    // Transform and quant residual to coeffs
    transform2d(block,pre_quant_coeff,width,0);
    quant(encoder, pre_quant_coeff, coeff_y, width, width, &ac_sum, 0, scan_idx_luma, ti->block_type);

    // Check for non-zero coeffs
    for (i = 0; i < width * width; i++) {
      if (coeff_y[i] != 0) {
        // Found one, we can break here
        cb_y = 1;
        break;
      }
    }
        
    // if non-zero coeffs
    if (cb_y) {
      // RECONSTRUCT for predictions
      dequant(encoder, coeff_y, pre_quant_coeff, width, width, 0, ti->block_type);
      itransform2d(block,pre_quant_coeff,width,0);

      i = 0;

      for (y = 0; y < LCU_WIDTH >> depth; y++) {
        for (x = 0; x < LCU_WIDTH >> depth; x++) {
          int16_t val = block[i++] + pred_y[x + y * pred_stride];
          //TODO: support 10+bits
          recbase_y[x + y * recbase_stride] = (uint8_t)CLIP(0, 255, val);
        }
      }
      // END RECONTRUCTION
    } else {
      // without coeffs, we only use the prediction
      for (y = 0; y < LCU_WIDTH >> depth; y++) {
        for (x = 0; x < LCU_WIDTH >> depth; x++) {
          recbase_y[x + y * recbase_stride] = (uint8_t)CLIP(0, 255, pred_y[x + y * pred_stride]);
    }
        }
      }

    if (encoder->in.video_format != FORMAT_400) {
      // Chroma U
      i = 0;
      ac_sum = 0;

      for (y = 0; y < LCU_WIDTH >> (depth + 1); y++) {
        for (x = 0; x < LCU_WIDTH >> (depth + 1); x++) {
          block[i] = ((int16_t)base_u[x + y * (base_stride >> 1)]) -
                     pred_u[x + y * (pred_stride >> 1)];
          i++;
        }
      }

      transform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),65535);
      quant(encoder, pre_quant_coeff, coeff_u, width >> 1, width >> 1, &ac_sum, 2,
            scan_idx_chroma, ti->block_type);

      for (i = 0; i < width *width >> 2; i++) {
        if (coeff_u[i] != 0) {
          // Found one, we can break here
          cb_u = 1;
          break;
        }
      }

      // Chroma V
      i = 0;
      ac_sum = 0;

      for (y = 0; y < LCU_WIDTH >> (depth + 1); y++) {
        for (x = 0; x < LCU_WIDTH >> (depth + 1); x++) {
          block[i] = ((int16_t)base_v[x + y * (base_stride >> 1)]) -
                     pred_v[x + y * (pred_stride >> 1)];
          i++;
        }
      }

      transform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),65535);
      quant(encoder, pre_quant_coeff, coeff_v, width >> 1, width >> 1, &ac_sum, 3,
            scan_idx_chroma, ti->block_type);

      for (i = 0; i < width *width >> 2; i++) {
        if (coeff_v[i] != 0) {
          // Found one, we can break here
          cb_v = 1;
          break;
        }
      }
          
      if (cb_u) {
        // RECONSTRUCT for predictions
        dequant(encoder, coeff_u, pre_quant_coeff, width >> 1, width >> 1, 2, ti->block_type);
        itransform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),65535);

        i = 0;

        for (y = 0; y < LCU_WIDTH >> (depth + 1); y++) {
          for (x = 0; x < LCU_WIDTH >> (depth + 1); x++) {
            int16_t val = block[i++] + pred_u[x + y * (pred_stride >> 1)];
            //TODO: support 10+bits
            recbase_u[x + y * (recbase_stride >> 1)] = (uint8_t)CLIP(0, 255, val);
          }
        }

        // END RECONTRUCTION
      } else {
        // without coeffs, we only use the prediction
        for (y = 0; y < LCU_WIDTH >> (depth + 1); y++) {
          for (x = 0; x < LCU_WIDTH >> (depth + 1); x++) {
            recbase_u[x + y * (recbase_stride >> 1)] = (uint8_t)CLIP(0, 255,
                                                                     pred_u[x + y * (pred_stride >> 1)]);
      }
          }
        }
      
      if (cb_v) {
        // RECONSTRUCT for predictions
        dequant(encoder, coeff_v, pre_quant_coeff, width >> 1, width >> 1, 3, ti->block_type);
        itransform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),65535);

        i = 0;

        for (y = 0; y < LCU_WIDTH >> (depth + 1); y++) {
          for (x = 0; x < LCU_WIDTH >> (depth + 1); x++) {
            int16_t val = block[i++] + pred_v[x + y * (pred_stride >> 1)];
            //TODO: support 10+bits
            recbase_v[x + y * (recbase_stride >> 1)] = (uint8_t)CLIP(0, 255, val);
          }
        }

        // END RECONTRUCTION
      } else {
        // without coeffs, we only use the prediction
        for (y = 0; y < LCU_WIDTH >> (depth + 1); y++) {
          for (x = 0; x < LCU_WIDTH >> (depth + 1); x++) {
            recbase_v[x + y * (recbase_stride >> 1)] = (uint8_t)CLIP(0, 255,
                                                                     pred_v[x + y * (pred_stride >> 1)]);
      }
          }
        }
      }
    
    // Store coded block pattern
    ti->cb[ti->idx] = cb_y | (cb_u << 1) | (cb_v << 2);
    // END INTRAPREDICTION
    return;
  }

  // end Residual Coding
}

void encode_transform_coeff(encoder_control *encoder, transform_info *ti,
                            int8_t depth, int8_t tr_depth)
{
  int8_t width = LCU_WIDTH>>depth;
  int8_t split = ((ti->split[ti->idx]&(1<<depth))||!depth)?1:0;
  int8_t cb_y, cb_u, cb_v;
  int32_t coeff_fourth = ((LCU_WIDTH>>(depth))*(LCU_WIDTH>>(depth)))+1;
  
  if (depth != 0 && depth != MAX_DEPTH + 1) {
    cabac.ctx = &g_trans_subdiv_model[5 - ((g_convert_to_bit[LCU_WIDTH] + 2) -
                                           depth)];
    CABAC_BIN(&cabac,split,"TransformSubdivFlag");
  }

  // Signal if chroma data is present
  // Chroma data is also signaled BEFORE transform split
  // Chroma data is not signaled if it was set to 0 before split
  if (encoder->in.video_format != FORMAT_400) {
    // Non-zero chroma U Tcoeffs
    int8_t cb_flag = (tr_depth == 0) ? ti->cb_top[1] : ((ti->cb[ti->idx] & 0x2) ? 1
                                                        : 0);
    cabac.ctx = &g_qt_cbf_model_chroma[tr_depth];

    if (tr_depth == 0 || ti->cb_top[1]) {
      CABAC_BIN(&cabac, cb_flag, "cbf_chroma_u");
    }

    // Non-zero chroma V Tcoeffs
    // NOTE: Using the same ctx as before
    cb_flag = (tr_depth == 0) ? ti->cb_top[2] : ((ti->cb[ti->idx] & 0x4) ? 1 : 0);

    if (tr_depth == 0 || ti->cb_top[2]) {
      CABAC_BIN(&cabac, cb_flag, "cbf_chroma_v");
    }
  }
  
  if (split) {
    ti->idx = 0; encode_transform_coeff(encoder, ti, depth + 1, tr_depth + 1);
    ti->idx = 1; encode_transform_coeff(encoder, ti, depth + 1, tr_depth + 1);
    ti->idx = 2; encode_transform_coeff(encoder, ti, depth + 1, tr_depth + 1);
    ti->idx = 3; encode_transform_coeff(encoder, ti, depth + 1, tr_depth + 1);
    return;
  }

  cb_y = ti->cb[ti->idx] & 0x1;
  cb_u = (ti->cb[ti->idx] & 0x2) ? 1 : 0;
  cb_v = (ti->cb[ti->idx] & 0x4) ? 1 : 0;
  if(ti->block_type == CU_INTRA || tr_depth || cb_u || cb_v) {
      // Non-zero luma Tcoeffs
      cabac.ctx = &g_qt_cbf_model_luma[tr_depth ? 0 : 1];
      CABAC_BIN(&cabac, cb_y, "cbf_luma");
  }


  {
    uint32_t ctx_idx;
    uint32_t scan_idx = SCAN_DIAG;
    uint32_t dir_mode;

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
    if (cb_y) {
      if (ti->block_type == CU_INTER) {
        scan_idx = SCAN_DIAG;
      } else {
        // Luma (Intra) scanmode
        dir_mode = ti->intra_pred_mode;

        //if multiple scans supported for transform size
        if (ctx_idx > 3 && ctx_idx < 6) {
          scan_idx = abs((int32_t) dir_mode - 26) < 5 ? 1 : (abs((int32_t)dir_mode - 10) < 5 ? 2 : 0);
        }
      }

      encode_coeff_nxn(encoder, &ti->coeff[0][ti->idx * coeff_fourth], width, 0, scan_idx);
    }

    if (cb_u || cb_v) {
      int8_t chroma_width = width >> 1;
      if(ti->block_type == CU_INTER) {
        scan_idx = SCAN_DIAG;
      } else {
        // Chroma scanmode
        ctx_idx++;
        dir_mode = ti->intra_pred_mode_chroma;

        if (dir_mode == 36) {
          // TODO: support NxN
          dir_mode = ti->intra_pred_mode;
        }

        scan_idx = SCAN_DIAG;

        if (ctx_idx > 4 && ctx_idx < 7) { // if multiple scans supported for transform size
          scan_idx = abs((int32_t) dir_mode - 26) < 5 ? 1 : (abs((int32_t)dir_mode - 10) < 5 ? 2 : 0);
        }
      }

      if (cb_u) {
        encode_coeff_nxn(encoder, &ti->coeff[1][ti->idx * coeff_fourth >> 1],
                         chroma_width, 2, scan_idx);
      }

      if (cb_v) {
        encode_coeff_nxn(encoder, &ti->coeff[2][ti->idx * coeff_fourth >> 1],
                         chroma_width, 2, scan_idx);
      }
    }
  }
}

void encode_coeff_nxn(encoder_control *encoder, int16_t *coeff, uint8_t width,
                      uint8_t type, int8_t scan_mode)
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
  const uint32_t *scan_cg         = NULL;

  // Init base contexts according to block type
  cabac_ctx *base_coeff_group_ctx = &g_cu_sig_coeff_group_model[type];
  cabac_ctx *baseCtx           = (type == 0) ? &g_cu_sig_model_luma[0] :
                                 &g_cu_sig_model_chroma[0];
  memset(sig_coeffgroup_flag,0,sizeof(uint32_t)*64);
  
  // Count non-zero coeffs
  for (i = 0; i < width * width; i++) {
    if (coeff[i] != 0) {
      num_nonzero++;
    }
  }  

  scan_cg = g_sig_last_scan[scan_mode][log2_block_size > 3 ? log2_block_size - 3 : 0];

  if (log2_block_size == 3) {
    scan_cg = g_sig_last_scan_8x8[scan_mode];
  } else if (log2_block_size == 5) {
    scan_cg = g_sig_last_scan_32x32;
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
  last_coeff_y = pos_last >> log2_block_size;

  // Code last_coeff_x and last_coeff_y
  encode_last_significant_xy(encoder, last_coeff_x, last_coeff_y, width, width,
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
      cabac.ctx = &base_coeff_group_ctx[ctx_sig];
      CABAC_BIN(&cabac, sig_coeff_group, "significant_coeff_group");
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
                                             log2_block_size, width, type);
          cabac.ctx = &baseCtx[ctx_sig];
          CABAC_BIN(&cabac, sig, "significant_coeff_flag");
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

      base_ctx_mod     = (type == 0) ? &g_cu_one_model_luma[4 * ctx_set] :
                         &g_cu_one_model_chroma[4 * ctx_set];
      num_c1_flag      = MIN(num_non_zero, C1FLAG_NUMBER);
      first_c2_flag_idx = -1;

      for (idx = 0; idx < num_c1_flag; idx++) {
        uint32_t symbol = (abs_coeff[idx] > 1) ? 1 : 0;
        cabac.ctx = &base_ctx_mod[c1];
        CABAC_BIN(&cabac, symbol, "significant_coeff2_flag");

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
        base_ctx_mod = (type == 0) ? &g_cu_abs_model_luma[ctx_set] :
                       &g_cu_abs_model_chroma[ctx_set];

        if (first_c2_flag_idx != -1) {
          uint8_t symbol = (abs_coeff[first_c2_flag_idx] > 2) ? 1 : 0;
          cabac.ctx      = &base_ctx_mod[0];
          CABAC_BIN(&cabac,symbol,"first_c2_flag");
        }
      }      
      
      if (be_valid && sign_hidden) {
        CABAC_BINS_EP(&cabac, (coeff_signs >> 1), (num_non_zero - 1), "");
      } else {
        CABAC_BINS_EP(&cabac, coeff_signs, num_non_zero, "");
      }
              
      if (c1 == 0 || num_non_zero > C1FLAG_NUMBER) {
        first_coeff2 = 1;

        for (idx = 0; idx < num_non_zero; idx++) {
          int32_t base_level  = (idx < C1FLAG_NUMBER) ? (2 + first_coeff2) : 1;

          if (abs_coeff[idx] >= base_level) {
            cabac_write_coeff_remain(&cabac, abs_coeff[idx] - base_level, go_rice_param);

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
void encode_last_significant_xy(encoder_control *encoder, 
                                uint8_t lastpos_x, uint8_t lastpos_y,
                                uint8_t width, uint8_t height,
                                uint8_t type, uint8_t scan)
{
  uint8_t offset_x  = type?0:((TOBITS(width)*3) + ((TOBITS(width)+1)>>2)),offset_y = offset_x;
  uint8_t shift_x   = type?(TOBITS(width)):((TOBITS(width)+3)>>2), shift_y = shift_x;
  int group_idx_x;
  int group_idx_y;
  int last_x,last_y,i;
  cabac_ctx *base_ctx_x = (type ? g_cu_ctx_last_x_chroma : g_cu_ctx_last_x_luma);
  cabac_ctx *base_ctx_y = (type ? g_cu_ctx_last_y_chroma : g_cu_ctx_last_y_luma);

  if (scan == SCAN_VER) {
    SWAP( lastpos_x, lastpos_y,uint8_t );
  }

  group_idx_x   = g_group_idx[lastpos_x];
  group_idx_y   = g_group_idx[lastpos_y];

  // Last X binarization
  for (last_x = 0; last_x < group_idx_x ; last_x++) {
    cabac.ctx = &base_ctx_x[offset_x + (last_x >> shift_x)];
    CABAC_BIN(&cabac,1,"LastSignificantX");
  }

  if (group_idx_x < g_group_idx[width - 1]) {
    cabac.ctx = &base_ctx_x[offset_x + (last_x >> shift_x)];
    CABAC_BIN(&cabac,0,"LastSignificantX");
  }

  // Last Y binarization
  for (last_y = 0; last_y < group_idx_y ; last_y++) {
    cabac.ctx = &base_ctx_y[offset_y + (last_y >> shift_y)];
    CABAC_BIN(&cabac,1,"LastSignificantY");
  }

  if (group_idx_y < g_group_idx[height - 1]) {
    cabac.ctx = &base_ctx_y[offset_y + (last_y >> shift_y)];
    CABAC_BIN(&cabac,0,"LastSignificantY");
  }

  // Last X
  if (group_idx_x > 3) {
    lastpos_x -= g_min_in_group[group_idx_x];

    for (i = ((group_idx_x - 2) >> 1) - 1; i >= 0; i--) {
      CABAC_BIN_EP(&cabac,(lastpos_x>>i) & 1,"LastSignificantX");
    }
  }          

  // Last Y
  if (group_idx_y > 3) {
    lastpos_y -= g_min_in_group[group_idx_y];

    for (i = ((group_idx_y - 2) >> 1) - 1; i >= 0; i--) {
      CABAC_BIN_EP(&cabac,(lastpos_y>>i) & 1,"LastSignificantY");
    }
  }

  // end LastSignificantXY
}
