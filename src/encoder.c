/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file encoder.c
    \brief Encoding related functions
    \author Marko Viitanen
    \date 2012-06
    
    Encoder main level
*/
/* Suppress some windows warnings */
#ifdef WIN32
  #define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "config.h"
#include "encoder.h"
#include "cabac.h"
#include "picture.h"
#include "nal.h"

void init_encoder_control(encoder_control* control,bitstream* output) {control->stream = output;};
void init_encoder_input(encoder_input* input,FILE* inputfile, uint32_t width, uint32_t height)
{
  input->file = inputfile;
  input->width = width;
  input->height = height;

  input->cur_pic.width = width;
  input->cur_pic.height = height;
  input->cur_pic.referenced = 0;
  /* Allocate buffers */
  input->cur_pic.yData = (uint8_t *)malloc(width*height);
  input->cur_pic.uData = (uint8_t *)malloc((width*height)>>2);
  input->cur_pic.vData = (uint8_t *)malloc((width*height)>>2);
};


void encode_one_frame(encoder_control* encoder)
{
  //output parameters before first frame
  if(encoder->frame == 0)
  {
    encode_seq_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 1, NAL_SEQ_PARAMETER_SET, 0);
        
    bitstream_clear_buffer(encoder->stream);

    encode_pic_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 1, NAL_PIC_PARAMETER_SET, 0);

  }

}

void encode_pic_parameter_set(encoder_control* encoder)
{
#ifdef _DEBUG
  printf("=========== Picture Parameter Set ID: 0 ===========\n");
#endif
  WRITE_UE(encoder->stream, 0, "pic_parameter_set_id");
  WRITE_UE(encoder->stream, 0, "seq_parameter_set_id");

  WRITE_U(encoder->stream, 0, 1, "sign_data_hiding_flag");
  WRITE_U(encoder->stream, 0, 1, "cabac_init_present_flag");

  WRITE_U(encoder->stream, 0, 3, "num_ref_idx_l0_default_active_minus1");
  WRITE_U(encoder->stream, 0, 3, "num_ref_idx_l1_default_active_minus1");
  /*
  WRITE_UE(encoder->stream, 0, "num_ref_idx_l0_default_active_minus1");
  WRITE_UE(encoder->stream, 0, "num_ref_idx_l1_default_active_minus1");
  */
  WRITE_SE(encoder->stream, 0, "pic_init_qp_minus26");

  WRITE_U(encoder->stream, 0, 1, "constrained_intra_pred_flag");
  WRITE_U(encoder->stream, 0, 1, "enable_temporal_mvp_flag");

  WRITE_U(encoder->stream, 0, 2, "slice_granularity");

  WRITE_UE(encoder->stream, 0, "max_cu_qp_delta_depth");

  WRITE_SE(encoder->stream, 0, "cb_qp_offset");
  WRITE_SE(encoder->stream, 0, "cr_qp_offset");

  WRITE_U(encoder->stream, 0, 1, "weighted_pred_flag");
  WRITE_U(encoder->stream, 0, 2, "weighted_bipred_idc");
  
  WRITE_U(encoder->stream, 1, 1, "output_flag_present_flag");
  
  WRITE_U(encoder->stream, 0, 1, "deblocking_filter_control_present_flag");

  WRITE_UE(encoder->stream, 0, "log2_parallel_merge_level_minus2");

  WRITE_U(encoder->stream, 0, 1, "pps_extension_flag"); 
  
}

void encode_seq_parameter_set(encoder_control* encoder)
{
#ifdef _DEBUG
  printf("=========== Sequence Parameter Set ID: 0 ===========\n");
#endif
  WRITE_U(encoder->stream, 0, 8, "profile_idc");
  WRITE_U(encoder->stream, 0, 8, "reserved_zero_8bits");
  WRITE_U(encoder->stream, 0, 8, "level_idc");
  WRITE_UE(encoder->stream, 0, "seq_parameter_set_id");
  WRITE_UE(encoder->stream, 0, "chroma_format_idc"); /* 0 = 4:0:0, 1 = 4:2:0, 2 = 4:2:2, 3 = 4:4:4 */
  WRITE_U(encoder->stream, 0, 3, "max_temporal_layers_minus1");
  WRITE_UE(encoder->stream, encoder->in.width, "pic_width_in_luma_samples");
  WRITE_UE(encoder->stream, encoder->in.height, "pic_height_in_luma_samples");
  WRITE_U(encoder->stream, 0, 1, "pic_cropping_flag");

  WRITE_UE(encoder->stream, 0, "bit_depth_luma_minus8");
  WRITE_UE(encoder->stream, 0, "bit_depth_chroma_minus8");

  WRITE_U(encoder->stream, 0, 1, "pcm_enabled_flag");

  WRITE_U(encoder->stream, 0, 1, "qpprime_y_zero_transquant_bypass_flag");

  WRITE_UE(encoder->stream, 0, "log2_max_pic_order_cnt_lsb_minus4");

  WRITE_UE(encoder->stream, 0, "max_dec_pic_buffering");
  WRITE_UE(encoder->stream, 0, "num_reorder_pics");
  WRITE_UE(encoder->stream, 0, "max_latency_increase");

  WRITE_U(encoder->stream, 0, 1, "restricted_ref_pic_lists_flag");

  WRITE_UE(encoder->stream, 0, "log2_min_coding_block_size_minus3");
  WRITE_UE(encoder->stream, 3, "log2_diff_max_min_coding_block_size");
  WRITE_UE(encoder->stream, 0, "log2_min_transform_block_size_minus2");
  WRITE_UE(encoder->stream, 3, "log2_diff_max_min_transform_block_size");

  WRITE_U(encoder->stream, 0, 1, "unknown_flag");

  WRITE_UE(encoder->stream, 2, "max_transform_hierarchy_depth_inter");
  WRITE_UE(encoder->stream, 2, "max_transform_hierarchy_depth_intra");
  
  WRITE_U(encoder->stream, 0, 1, "scaling_list_enable_flag");
  WRITE_U(encoder->stream, 0, 1, "chroma_pred_from_luma_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "transform_skip_enabled_flag");

  WRITE_U(encoder->stream, 0, 1, "deblocking_filter_in_aps_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "seq_loop_filter_across_slices_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "asymmetric_motion_partitions_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "nsrqt_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "sample_adaptive_offset_enabled_flag");
	WRITE_U(encoder->stream, 0, 1, "adaptive_loop_filter_enabled_flag");
	
  WRITE_U(encoder->stream, 0, 1, "temporal_id_nesting_flag");
	
  WRITE_UE(encoder->stream, 0, "num_short_term_ref_pic_sets");
	
  WRITE_U(encoder->stream, 0, 1, "long_term_ref_pics_present_flag");
	
  WRITE_U(encoder->stream, 0, 2, "tiles_or_entropy_coding_sync_idc");
	
	WRITE_U(encoder->stream, 0, 1, "sps_extension_flag");
  
}
  

