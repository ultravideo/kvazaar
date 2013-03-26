/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file encoder.c
    \brief Encoding related functions
    \author Marko Viitanen
    \date 2013-03
    
    Encoder main level
*/
/* Suppress some visual studio warnings */
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
#include "context.h"
#include "transform.h"
#include "intra.h"

void initSigLastScan(uint32_t* pBuffD, uint32_t* pBuffH, uint32_t* pBuffV, int32_t iWidth, int32_t iHeight)
{
  uint32_t uiNumScanPos  = iWidth * iWidth;
  uint32_t uiNextScanPos = 0;
  int32_t  iX,iY,x,y;
  uint32_t uiScanLine;
  uint32_t blkY,blkX;
  uint32_t uiBlk;
  uint32_t uiCnt = 0;

  if( iWidth < 16 )
  {
    uint32_t* pBuffTemp = pBuffD;
    if( iWidth == 8 )
    {
      pBuffTemp = (uint32_t *)g_sigLastScanCG32x32;
    }
    for( uiScanLine = 0; uiNextScanPos < uiNumScanPos; uiScanLine++ )
    {
      int    iPrimDim  = uiScanLine;
      int    iScndDim  = 0;
      while( iPrimDim >= iWidth )
      {
        iScndDim++;
        iPrimDim--;
      }
      while( iPrimDim >= 0 && iScndDim < iWidth )
      {
        pBuffTemp[ uiNextScanPos ] = iPrimDim * iWidth + iScndDim ;
        uiNextScanPos++;
        iScndDim++;
        iPrimDim--;
      }
    }
  }
  if( iWidth > 4 )
  {
    uint32_t uiNumBlkSide = iWidth >> 2;
    uint32_t uiNumBlks    = uiNumBlkSide * uiNumBlkSide;
    uint32_t log2Blk      = g_aucConvertToBit[ uiNumBlkSide ] + 1;

    for(uiBlk = 0; uiBlk < uiNumBlks; uiBlk++ )
    {      
      uint32_t initBlkPos = g_auiSigLastScan[ SCAN_DIAG ][ log2Blk ][ uiBlk ];
      uiNextScanPos   = 0;
      if( iWidth == 32 )
      {
        initBlkPos = g_sigLastScanCG32x32[ uiBlk ];
      }
      {
        uint32_t offsetY    = initBlkPos / uiNumBlkSide;
        uint32_t offsetX    = initBlkPos - offsetY * uiNumBlkSide;
        uint32_t offsetD    = 4 * ( offsetX + offsetY * iWidth );
        uint32_t offsetScan = 16 * uiBlk;
        for( uiScanLine = 0; uiNextScanPos < 16; uiScanLine++ )
        {
          int    iPrimDim  = uiScanLine;
          int    iScndDim  = 0;
          //ToDo: optimize
          while( iPrimDim >= 4 )
          {
            iScndDim++;
            iPrimDim--;
          }
          while( iPrimDim >= 0 && iScndDim < 4 )
          {
            pBuffD[ uiNextScanPos + offsetScan ] = iPrimDim * iWidth + iScndDim + offsetD;
            uiNextScanPos++;
            iScndDim++;
            iPrimDim--;
          }
        }
      }
    }
  }  
  
  if( iWidth > 2 )
  {
    uint32_t numBlkSide = iWidth >> 2;
    for(blkY=0; blkY < numBlkSide; blkY++)
    {
      for(blkX=0; blkX < numBlkSide; blkX++)
      {
        uint32_t offset    = blkY * 4 * iWidth + blkX * 4;
        for(y=0; y < 4; y++)
        {
          for(x=0; x < 4; x++)
          {
            pBuffH[uiCnt] = y*iWidth + x + offset;
            uiCnt ++;
          }
        }
      }
    }
    uiCnt = 0;
    for(blkX=0; blkX < numBlkSide; blkX++)
    {
      for(blkY=0; blkY < numBlkSide; blkY++)
      {
        uint32_t offset    = blkY * 4 * iWidth + blkX * 4;
        for(x=0; x < 4; x++)
        {
          for(y=0; y < 4; y++)
          {
            pBuffV[uiCnt] = y*iWidth + x + offset;
            uiCnt ++;
          }
        }
      }
    }
  }
  else
  {
    for(iY=0; iY < iHeight; iY++)
    {
      for(iX=0; iX < iWidth; iX++)
      {
        pBuffH[uiCnt] = iY*iWidth + iX;
        uiCnt ++;
      }
    }

    uiCnt = 0;
    for(iX=0; iX < iWidth; iX++)
    {
      for(iY=0; iY < iHeight; iY++)
      {
        pBuffV[uiCnt] = iY*iWidth + iX;
        uiCnt ++;
      }
    }
  }
}

void init_tables(void)
{
  int i;
  int c = 0;
  memset( g_aucConvertToBit,-1, sizeof( g_aucConvertToBit ) );  
  for ( i=4; i<(1<<7); i*=2 )
  {
    g_aucConvertToBit[i] = c;
    c++;
  }
  g_aucConvertToBit[i] = c;

  c = 2;
  for ( i=0; i<7; i++ )
  {
    g_auiSigLastScan[0][i] = (uint32_t*)malloc(c*c*sizeof(uint32_t));
    g_auiSigLastScan[1][i] = (uint32_t*)malloc(c*c*sizeof(uint32_t));
    g_auiSigLastScan[2][i] = (uint32_t*)malloc(c*c*sizeof(uint32_t));

    initSigLastScan( g_auiSigLastScan[0][i], g_auiSigLastScan[1][i], g_auiSigLastScan[2][i], c, c);
    c <<= 1;
  }

}
void init_encoder_control(encoder_control* control,bitstream* output)
{
  control->stream = output;  
}

void init_encoder_input(encoder_input* input,FILE* inputfile, uint32_t width, uint32_t height)
{
  int i;
  input->file = inputfile;
  input->width = width;
  input->height = height;

  input->height_in_LCU = height / LCU_WIDTH;
  input->width_in_LCU  =  width / LCU_WIDTH;

  /* Add one extra LCU when image not divisible by LCU_WIDTH */
  if(input->height_in_LCU * LCU_WIDTH < height)
  {
    input->height_in_LCU++;
  }
  if(input->width_in_LCU * LCU_WIDTH < width)
  {
    input->width_in_LCU++;
  }

  input->cur_pic.width  = width;
  input->cur_pic.height = height;
  input->cur_pic.referenced = 0;
  /* Allocate buffers */
  input->cur_pic.yData = (uint8_t *)malloc(width*height);
  input->cur_pic.uData = (uint8_t *)malloc((width*height)>>2);
  input->cur_pic.vData = (uint8_t *)malloc((width*height)>>2);

  /* Reconstruction buffers */
  input->cur_pic.yRecData = (uint8_t *)malloc(width*height);
  input->cur_pic.uRecData = (uint8_t *)malloc((width*height)>>2);
  input->cur_pic.vRecData = (uint8_t *)malloc((width*height)>>2);

  memset(input->cur_pic.uRecData, 128, (width*height)>>2);
  memset(input->cur_pic.vRecData, 128, (width*height)>>2);

  /* Allocate memory for CU info 2D array */
  //ToDo: we don't need this much space on LCU...MAX_DEPTH-1
  input->cur_pic.CU = (CU_info**)malloc((MAX_DEPTH+1)*sizeof(CU_info*));
  for(i=0; i < MAX_DEPTH+1; i++)
  {
    /* Allocate height_in_SCU x width_in_SCU x sizeof(CU_info) */
    input->cur_pic.CU[i] = (CU_info*)malloc((input->height_in_LCU<<MAX_DEPTH)*(input->width_in_LCU<<MAX_DEPTH)*sizeof(CU_info));
    memset(input->cur_pic.CU[i], 0, (input->height_in_LCU<<MAX_DEPTH)*(input->width_in_LCU<<MAX_DEPTH)*sizeof(CU_info));
  }  
}


void encode_one_frame(encoder_control* encoder)
{
  int i;
  /* output parameters before first frame */
  if(encoder->frame == 0)
  {
    /* Sequence Parameter Set (SPS) */
    encode_seq_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, NAL_SEQ_PARAMETER_SET, 0);
    bitstream_clear_buffer(encoder->stream);

    /* Video Parameter Set (VPS) */    
    encode_vid_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, NAL_VID_PARAMETER_SET, 0);
    bitstream_clear_buffer(encoder->stream);
    
    /* Picture Parameter Set (PPS) */
    encode_pic_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, NAL_PIC_PARAMETER_SET, 0);
    bitstream_clear_buffer(encoder->stream);

    /* First slice is IDR */
    cabac_start(&cabac);
    encoder->in.cur_pic.slicetype = SLICE_I;
    encoder->in.cur_pic.type = NAL_IDR_SLICE;
    encode_slice_header(encoder);
    bitstream_align(encoder->stream);
    encode_slice_data(encoder);
    cabac_flush(&cabac);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, NAL_IDR_SLICE, 0);
    bitstream_clear_buffer(encoder->stream);
  }  
  else
  {
    /* ToDo: add intra/inter search before encoding */

    cabac_start(&cabac);
    encoder->in.cur_pic.slicetype = SLICE_I;
    encoder->in.cur_pic.type = 0;
    encode_slice_header(encoder);
    bitstream_align(encoder->stream);
    encode_slice_data(encoder);
    cabac_flush(&cabac);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, 0, encoder->frame);
    bitstream_clear_buffer(encoder->stream);
  }  
  #ifdef _DEBUG
  /*
  {
    int x,y;
    for(y = 0; y < encoder->in.height_in_LCU*2;y++)
    {
      for(x = 0;x < encoder->in.width_in_LCU*2;x++)
      {
        i = (x<<2)+(y<<2)*(encoder->in.width_in_LCU<<MAX_DEPTH);
        printf("(%d,%d) Intramode: %d\n", x<<2, y<<2,encoder->in.cur_pic.CU[0][i].intra.mode);
      }
    }
  }
  */
  #endif

  /* Clear prediction data */
  /* ToDo: store as reference data */
  for(i=0; i < MAX_DEPTH+1; i++)
  {    
    memset(encoder->in.cur_pic.CU[i], 0, (encoder->in.height_in_LCU<<MAX_DEPTH)*(encoder->in.width_in_LCU<<MAX_DEPTH)*sizeof(CU_info));
  } 

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

  //ToDo: add QP offsets
  WRITE_SE(encoder->stream, 0, "pps_cb_qp_offset");
  WRITE_SE(encoder->stream, 0, "pps_cr_qp_offset");
  WRITE_U(encoder->stream, 0, 1, "pps_slice_chroma_qp_offsets_present_flag");
  WRITE_U(encoder->stream, 0, 1, "weighted_pred_flag");
  WRITE_U(encoder->stream, 0, 1, "weighted_bipred_idc");

  //WRITE_U(encoder->stream, 0, 1, "dependent_slices_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "transquant_bypass_enable_flag");
  WRITE_U(encoder->stream, 0, 1, "tiles_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "entropy_coding_sync_enabled_flag");
  //ToDo: enable tiles for concurrency
  //IF tiles
  //ENDIF
  WRITE_U(encoder->stream, 0, 1, "loop_filter_across_slice_flag");
  WRITE_U(encoder->stream, 1, 1, "deblocking_filter_control_present_flag");
  //IF deblocking_filter
    WRITE_U(encoder->stream, 0, 1, "deblocking_filter_override_enabled_flag");
    WRITE_U(encoder->stream, 1, 1, "pps_disable_deblocking_filter_flag");
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
  /*PTL*/
  /*Profile Tier*/
  WRITE_U(encoder->stream, 0, 2, "XXX_profile_space[]");
  WRITE_U(encoder->stream, 0, 1, "XXX_tier_flag[]");
  WRITE_U(encoder->stream, 0, 5, "XXX_profile_idc[]");
  WRITE_U(encoder->stream, 0, 32, "XXX_profile_compatibility_flag[][j]");

  WRITE_U(encoder->stream, 0, 1, "general_progressive_source_flag");
  WRITE_U(encoder->stream, 0, 1, "general_interlaced_source_flag");
  WRITE_U(encoder->stream, 0, 1, "general_non_packed_constraint_flag");
  WRITE_U(encoder->stream, 0, 1, "general_frame_only_constraint_flag");

  WRITE_U(encoder->stream, 0, 32, "XXX_reserved_zero_44bits[0..31]");
  WRITE_U(encoder->stream, 0, 12, "XXX_reserved_zero_44bits[32..43]");
  
  /*end Profile Tier */

  WRITE_U(encoder->stream, 0, 8, "general_level_idc");

  WRITE_U(encoder->stream, 0, 1, "sub_layer_profile_present_flag");
  WRITE_U(encoder->stream, 0, 1, "sub_layer_level_present_flag");
  for(i = 1; i < 8; i++)
  {
    WRITE_U(encoder->stream, 0, 2, "reserved_zero_2bits");
  }

  /*end PTL*/
}

void encode_seq_parameter_set(encoder_control* encoder)
{
#ifdef _DEBUG
  printf("=========== Sequence Parameter Set ID: 0 ===========\n");
#endif
  /* ToDo: profile IDC and level IDC should be defined later on */
  WRITE_U(encoder->stream, 0, 4, "sps_video_parameter_set_id");
  WRITE_U(encoder->stream, 1, 3, "sps_max_sub_layers_minus1");

  WRITE_U(encoder->stream, 0, 1, "sps_temporal_id_nesting_flag");
    
  encode_PTL(encoder);

  WRITE_UE(encoder->stream, 0, "sps_seq_parameter_set_id");
  WRITE_UE(encoder->stream, encoder->in.video_format, "chroma_format_idc"); /* 0 = 4:0:0, 1 = 4:2:0, 2 = 4:2:2, 3 = 4:4:4 */
  if(encoder->in.video_format == 3)
  {
    WRITE_U(encoder->stream, 0, 1, "separate_colour_plane_flag");
  }
  WRITE_UE(encoder->stream, encoder->in.width, "pic_width_in_luma_samples");
  WRITE_UE(encoder->stream, encoder->in.height, "pic_height_in_luma_samples");
  WRITE_U(encoder->stream, 0, 1, "conformance_window_flag");
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
  WRITE_UE(encoder->stream, 0, "log2_min_transform_block_size_minus2");   /* 4x4 */
  WRITE_UE(encoder->stream, 3, "log2_diff_max_min_transform_block_size"); /* 4x4...32x32 */
  WRITE_UE(encoder->stream, 2, "max_transform_hierarchy_depth_inter");
  WRITE_UE(encoder->stream, 2, "max_transform_hierarchy_depth_intra");

  /* Use default scaling list */
  WRITE_U(encoder->stream, 1, 1, "scaling_list_enable_flag");
  //IF scaling list
    WRITE_U(encoder->stream, 0, 1, "sps_scaling_list_data_present_flag");
  //ENDIF

  WRITE_U(encoder->stream, 0, 1, "amp_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "sample_adaptive_offset_enabled_flag");

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
  WRITE_U(encoder->stream, 0, 1, "sps_temporal_mvp_enable_flag");

  WRITE_U(encoder->stream, 0, 1, "sps_strong_intra_smoothing_enable_flag");

  WRITE_U(encoder->stream, 0, 1, "vui_parameters_present_flag");
  //ToDo: VUI?
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
  for(i = 0; i < 1; i++)
  {
  WRITE_UE(encoder->stream, 1, "vps_max_dec_pic_buffering");
  WRITE_UE(encoder->stream, 0, "vps_num_reorder_pics");
  WRITE_UE(encoder->stream, 0, "vps_max_latency_increase");
  }
  //end for
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
#ifdef _DEBUG
  printf("=========== Slice ===========\n");
#endif

  WRITE_U(encoder->stream, 1, 1, "first_slice_segment_in_pic_flag");
  if(encoder->in.cur_pic.type == NAL_IDR_SLICE)
  {
    WRITE_U(encoder->stream, 0, 1, "no_output_of_prior_pics_flag");
  }
  WRITE_UE(encoder->stream, 0, "slice_pic_parameter_set_id");

  //WRITE_U(encoder->stream, 0, 1, "dependent_slice_segment_flag");
  
  /* ToDo: add more slice types */
  WRITE_UE(encoder->stream, encoder->in.cur_pic.slicetype, "slice_type");

  // if !entropy_slice_flag
  
    //if output_flag_present_flag
      //WRITE_U(encoder->stream, 1, 1, "pic_output_flag");
    //end if
    //if( IdrPicFlag ) <- nal_unit_type == 5
    if(encoder->in.cur_pic.type == NAL_IDR_SLICE)
    {
      //WRITE_UE(encoder->stream, encoder->frame&3, "idr_pic_id");      
    }
    else
    {
      WRITE_U(encoder->stream, encoder->frame, 4, "pic_order_cnt_lsb");
      WRITE_U(encoder->stream, 1, 1, "short_term_ref_pic_set_sps_flag");
      //WRITE_UE(encoder->stream, 0, "short_term_ref_pic_set_idx");
    }
    //end if
  //end if
    //IF sao
    /*
    WRITE_U(encoder->stream, 0,1, "slice_sao_luma_flag" );
    WRITE_U(encoder->stream, 0,1, "slice_sao_chroma_flag" );
    */
    //ENDIF
  /* Skip flags that are not present */
  // if !entropy_slice_flag
    WRITE_SE(encoder->stream, 0, "slice_qp_delta");
    //WRITE_U(encoder->stream, 1, 1, "alignment");
}

void encode_slice_data(encoder_control* encoder)
{
  uint16_t xCtb,yCtb;

  scalinglist_process();
  init_contexts(encoder,encoder->in.cur_pic.slicetype);

  /* Loop through every LCU in the slice */
  for(yCtb = 0; yCtb < encoder->in.height_in_LCU; yCtb++)
  {
    uint8_t lastCUy = (yCtb == (encoder->in.height_in_LCU-1))?1:0;
    for(xCtb = 0; xCtb < encoder->in.width_in_LCU; xCtb++)
    {
      uint8_t lastCUx = (xCtb == (encoder->in.width_in_LCU-1))?1:0;
      uint8_t depth = 0;

      /* Recursive function for looping through all the sub-blocks */
      encode_coding_tree(encoder, xCtb<<MAX_DEPTH,yCtb<<MAX_DEPTH, depth);

      /* signal Terminating bit */
      if(!lastCUx || !lastCUy)
      {
        cabac_encodeBinTrm(&cabac, 0);
      }
    }
  }
}

void encode_coding_tree(encoder_control* encoder,uint16_t xCtb,uint16_t yCtb, uint8_t depth)
{    
  uint8_t split_flag = (depth<3)?1:0; /* ToDo: get from CU data */
  uint8_t split_model = 0;

  /* Check for slice border */
  uint8_t border_x = ((encoder->in.width)<(uint32_t)( xCtb*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>depth) ))?1:0;
  uint8_t border_y = ((encoder->in.height)<(uint32_t)( yCtb*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>depth) ))?1:0;
  uint8_t border = border_x | border_y; /*!< are we in any border CU */
  CU_info *cur_CU = &encoder->in.cur_pic.CU[depth][(xCtb>>(MAX_DEPTH-depth))+(yCtb>>(MAX_DEPTH-depth))*(encoder->in.width_in_LCU<<MAX_DEPTH)];

  /* When not in MAX_DEPTH, insert split flag and split the blocks if needed */
  if(depth != MAX_DEPTH)
  {
    SET_SPLITDATA(cur_CU,split_flag);
    //Implisit split flag when on border
    if(!border)
    {
      /* Get left and top block split_flags and if they are present and true, increase model number */
      if(xCtb > 0 && GET_SPLITDATA(&(encoder->in.cur_pic.CU[depth][(xCtb>>(MAX_DEPTH-depth))-1+(yCtb>>(MAX_DEPTH-depth))*(encoder->in.width_in_LCU<<MAX_DEPTH)])) == 1)
      {
        split_model++;
      }
      if(yCtb > 0 && GET_SPLITDATA(&(encoder->in.cur_pic.CU[depth][(xCtb>>(MAX_DEPTH-depth))+((yCtb>>(MAX_DEPTH-depth))-1)*(encoder->in.width_in_LCU<<MAX_DEPTH)])) == 1)
      {
        split_model++;
      }
      cabac.ctx = &g_SplitFlagSCModel[split_model];    
      CABAC_BIN(&cabac, split_flag, "SplitFlag");
    }
    if(split_flag || border)
    {
      /* Split blocks and remember to change x and y block positions */
      uint8_t change = 1<<(MAX_DEPTH-1-depth);
      encode_coding_tree(encoder,xCtb,yCtb,depth+1);
      if(!border_x)
      {
        encode_coding_tree(encoder,xCtb+change,yCtb,depth+1);
      }
      if(!border_y)
      {
        encode_coding_tree(encoder,xCtb,yCtb+change,depth+1);      
      }
      if(!border)
      {
        encode_coding_tree(encoder,xCtb+change,yCtb+change,depth+1);
      }
      /* We don't need to do anything else here */
      return;
    }
  }
  
  /* Set every block as intra for now */
  {
    cur_CU->type = CU_INTRA;
  }

    /* Signal PartSize on max depth */
    if(depth == MAX_DEPTH)
    {
      cabac.ctx = &g_PartSizeSCModel;
      CABAC_BIN(&cabac, 1, "PartSize");
    }
    /*end partsize*/

    if(cur_CU->type == CU_INTRA)
    {
      uint8_t intraPredMode = 1;
      uint8_t intraPredModeChroma = 36; /* 36 = Chroma derived from luma */
      int8_t intraPreds[3] = {-1, -1, -1};
      int8_t mpmPred = -1;
      int i;
      uint32_t flag;
      uint32_t bestSAD;
      uint8_t *base  = &encoder->in.cur_pic.yData[xCtb*(LCU_WIDTH>>(MAX_DEPTH))   + (yCtb*(LCU_WIDTH>>(MAX_DEPTH)))  *encoder->in.width];
      uint8_t *baseU = &encoder->in.cur_pic.uData[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
      uint8_t *baseV = &encoder->in.cur_pic.vData[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
      uint32_t width = LCU_WIDTH>>depth;

      /* INTRAPREDICTION */
      /* ToDo: split to a function */
      int16_t pred[LCU_WIDTH*LCU_WIDTH];
      int16_t predU[LCU_WIDTH*LCU_WIDTH>>2];
      int16_t predV[LCU_WIDTH*LCU_WIDTH>>2];

      int16_t rec[(LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)];
      int16_t *recShift  = &rec[(LCU_WIDTH>>(depth))*2+8+1];
      int16_t *recShiftU = &rec[(LCU_WIDTH>>(depth+1))*2+8+1];
      uint8_t *recbase   = &encoder->in.cur_pic.yRecData[xCtb*(LCU_WIDTH>>(MAX_DEPTH))   + (yCtb*(LCU_WIDTH>>(MAX_DEPTH)))  *encoder->in.width];
      uint8_t *recbaseU  = &encoder->in.cur_pic.uRecData[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
      uint8_t *recbaseV  = &encoder->in.cur_pic.vRecData[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];

      #if ENABLE_PCM == 1
      /* Code must start after variable initialization */
      cabac_encodeBinTrm(&cabac, 0); /* IPCMFlag == 0 */
      #endif

      /* Build reconstructed block to use in prediction with extrapolated borders */      
      intra_buildReferenceBorder(&encoder->in.cur_pic, xCtb, yCtb,(LCU_WIDTH>>(depth))*2+8, rec, (LCU_WIDTH>>(depth))*2+8, 0);

      intraPredMode = (uint8_t)intra_prediction(encoder->in.cur_pic.yData,encoder->in.width,recShift,(LCU_WIDTH>>(depth))*2+8,xCtb*(LCU_WIDTH>>(MAX_DEPTH)),yCtb*(LCU_WIDTH>>(MAX_DEPTH)),width,pred,width,&bestSAD);
      
      /* Filter DC-prediction */
      if(intraPredMode == 1 && (LCU_WIDTH>>depth) < 32 && xCtb && yCtb)
      {
        intra_DCPredFiltering(recbase,encoder->in.width,recbase,encoder->in.width,LCU_WIDTH>>depth,LCU_WIDTH>>depth);
      }
      
      /* ToDo: separate chroma prediction(?) */
      /* intraPredModeChroma = 1; */
      intra_buildReferenceBorder(&encoder->in.cur_pic, xCtb, yCtb,(LCU_WIDTH>>(depth+1))*2+8, rec, (LCU_WIDTH>>(depth+1))*2+8, 1);
      intra_recon(recShiftU,(LCU_WIDTH>>(depth+1))*2+8,xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)),yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)),width>>1,predU,width>>1,intraPredModeChroma!=36?intraPredModeChroma:intraPredMode,1);
      intra_buildReferenceBorder(&encoder->in.cur_pic, xCtb, yCtb,(LCU_WIDTH>>(depth+1))*2+8, rec, (LCU_WIDTH>>(depth+1))*2+8, 2);
      intra_recon(recShiftU,(LCU_WIDTH>>(depth+1))*2+8,xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)),yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)),width>>1,predV,width>>1,intraPredModeChroma!=36?intraPredModeChroma:intraPredMode,1);
      
      /* This affects reconstruction, do after that */
      intra_setBlockMode(&encoder->in.cur_pic, xCtb, yCtb, depth, intraPredMode);

      /*
        PREDINFO CODING
        If intra prediction mode is found from the predictors,
        it can be signaled with two EP's. Otherwise we can send
        5 EP bins with the full predmode
        ToDo: split to a function
      */
      intra_getDirLumaPredictor(&encoder->in.cur_pic, xCtb, yCtb, depth, intraPreds);
            
      for(i = 0; i < 3; i++)
      {
        if(intraPreds[i] == intraPredMode)
        {
          mpmPred = i;
          break;
        }
      }
      /* For each part { */
      flag = (mpmPred==-1)?0:1;
      cabac.ctx = &g_IntraModeSCModel;
      CABAC_BIN(&cabac,flag,"IntraPred");
      /*} End for each part */

      /*           Intrapredmode signaling 
        If found from predictors, we can simplify signaling
      */
      if(flag)
      {
        flag = (mpmPred==0)?0:1;
        CABAC_BIN_EP(&cabac, flag, "intraPredMode");
        if(mpmPred!=0)
        {
          flag = (mpmPred==1)?0:1;
          CABAC_BIN_EP(&cabac, flag, "intraPredMode");
        }
      }
      else /* Else we signal the "full" predmode */
      {
        int8_t intraPredModeTemp = intraPredMode;
        if (intraPreds[0] > intraPreds[1])
        { 
          SWAP(intraPreds[0], intraPreds[1], int8_t); 
        }
        if (intraPreds[0] > intraPreds[2])
        {
          SWAP(intraPreds[0], intraPreds[2], int8_t); 
        }
        if (intraPreds[1] > intraPreds[2])
        {
          SWAP(intraPreds[1], intraPreds[2], int8_t); 
        }
        for(i = 2; i >= 0; i--)
        {
          intraPredModeTemp = intraPredModeTemp > intraPreds[i] ? intraPredModeTemp - 1 : intraPredModeTemp;
        }
        CABAC_BINS_EP(&cabac, intraPredModeTemp, 5, "intraPredMode");
      }

      /* If we have chroma, signal it */
      if(encoder->in.video_format != FORMAT_400)
      {
        /* Chroma intra prediction */
        cabac.ctx = &g_ChromaPredSCModel[0];
        CABAC_BIN(&cabac,((intraPredModeChroma!=36)?1:0),"IntraPredChroma");

        /* If not copied from luma, signal it */
        if(intraPredModeChroma!=36)
        {
          int8_t intraPredModeChromaTemp = intraPredModeChroma;
          /* Default chroma predictors */
          uint32_t allowedChromaDir[ 5 ] = { 0, 26, 10, 1, 36 };
          
          /* If intra is the same as one of the default predictors, replace it */
          for(i = 0; i < 4; i++ )
          {
            if( intraPredMode == allowedChromaDir[i] )
            {
              allowedChromaDir[i] = 34; /* VER+8 mode */
              break;
            }
          }

          for(i = 0; i < 4; i++ )
          {
            if( intraPredModeChromaTemp == allowedChromaDir[i] )
            {
              intraPredModeChromaTemp = i;
              break;
            }
          }
          CABAC_BINS_EP(&cabac, intraPredModeChromaTemp, 2, "intraPredModeChroma");
        }
      }
      /*
      END OF PREDINFO CODING
      */

      /* Coeff */
      /* Transform tree */
      {
        int16_t coeff[LCU_WIDTH*LCU_WIDTH*2];
        int16_t coeffU[LCU_WIDTH*LCU_WIDTH>>1];
        int16_t coeffV[LCU_WIDTH*LCU_WIDTH>>1];
        int8_t split = 0;
        /* Handle transforms, quant and reconstruction */
        int32_t cb = encode_transform_tree( encoder,base, baseU, baseV, encoder->in.width,
                                            recbase,recbaseU, recbaseV, encoder->in.width,
                                            pred,predU,predV,(LCU_WIDTH>>depth),
                                            coeff,coeffU,coeffV,
                                            depth, &split);

        /* Code coeffs to bitstream */
        encode_transform_coeff(encoder, coeff, coeffU, coeffV, cb, cb, cb,depth, intraPredMode, intraPredModeChroma, split?1:0, 1);
      }
      /* end Transform tree */
      /* end Coeff */

    }
    #if ENABLE_PCM == 1
    /* Code IPCM block */
    else if(cur_CU->type == CU_PCM)
    {
      cabac_encodeBinTrm(&cabac, 1); /* IPCMFlag == 1 */
      cabac_finish(&cabac);
      bitstream_align(cabac.stream);
       /* PCM sample */
      {
        uint8_t *base   = &encoder->in.cur_pic.yData[xCtb*(LCU_WIDTH>>(MAX_DEPTH))    + (yCtb*(LCU_WIDTH>>(MAX_DEPTH)))*encoder->in.width];
        uint8_t *baseCb = &encoder->in.cur_pic.uData[(xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*encoder->in.width/2)];
        uint8_t *baseCr = &encoder->in.cur_pic.vData[(xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*encoder->in.width/2)];
        for(y = 0; y < LCU_WIDTH>>depth; y++)
        {
          for(x = 0; x < LCU_WIDTH>>depth; x++)
          {
            bitstream_put(cabac.stream, base[x+y*encoder->in.width], 8);
          }
        }       
        if(encoder->in.video_format != FORMAT_400)
        {
          /* Cb */
          for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
          {
            for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
            {
              bitstream_put(cabac.stream, baseCb[x+y*(encoder->in.width>>1)], 8);
            }
          }

          /* Cr */
          for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
          {
            for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
            {
              bitstream_put(cabac.stream, baseCr[x+y*(encoder->in.width>>1)], 8);
            }
          }
        }
      }
      /* end PCM sample */
      cabac_start(&cabac);

    } /* end Code IPCM block */ 
    #endif /* END ENABLE_PCM */
    else /* Should not happend */
    {
      printf("UNHANDLED TYPE!\r\n");
      exit(1);
    }
   /* end prediction unit */
  /* end coding_unit */
  
}

int32_t encode_transform_tree(encoder_control* encoder,uint8_t *base, uint8_t *baseU, uint8_t *baseV,int32_t base_stride,
                                                    uint8_t *recbase,uint8_t *recbaseU, uint8_t *recbaseV,int32_t recbase_stride,
                                                    int16_t *pred, int16_t *predU, int16_t *predV,int32_t pred_stride,
                                                    int16_t *coeff, int16_t *coeffU, int16_t *coeffV,
                                                    uint8_t depth, int8_t* split)
{
  /* we have 64>>depth transform size */
  int x,y,i;
  int32_t width = LCU_WIDTH>>depth;
  int32_t half_width = width >> 1;
  uint8_t CbY = 0,CbU = 0,CbV = 0;

  if(depth == 0)
  {
    /* Prepare for multi-level splitting */
    *split |= 1<<depth;
  }

  if((*split)&(1<<depth))
  {
    int32_t recbase_offset_y   = recbase_stride*(half_width);
    int32_t base_offset_y      = base_stride*(half_width);
    int32_t pred_offset_y      = pred_stride*(half_width);
    int32_t recbase_offset_c_y = (recbase_stride>>1)*(half_width>>1);
    int32_t base_offset_c_y    = (base_stride>>1)*(half_width>>1);
    int32_t pred_offset_c_y    = (pred_stride>>1)*(half_width>>1);
    int32_t coeff_fourth = (LCU_WIDTH*LCU_WIDTH>>2);
    int32_t output       = 0;
    int32_t outhelper    = 0;
    output = encode_transform_tree(encoder,base, baseU, baseV, base_stride,
                          recbase,recbaseU, recbaseV, recbase_stride,
                          pred,predU,predV,pred_stride,
                          coeff,coeffU,coeffV,
                          depth+1, split);
    outhelper = encode_transform_tree(encoder,&base[half_width], &baseU[half_width>>1], &baseV[half_width>>1], base_stride,
                          &recbase[half_width],&recbaseU[half_width>>1], &recbaseV[half_width>>1], recbase_stride,
                          &pred[half_width],&predU[half_width>>1],&predV[half_width>>1],pred_stride,
                          &coeff[coeff_fourth],&coeffU[coeff_fourth>>1],&coeffV[coeff_fourth>>1],
                          depth+1, split);
    output |= outhelper<<3;
    outhelper = encode_transform_tree(encoder,&base[base_offset_y], &baseU[base_offset_c_y], &baseV[base_offset_c_y], base_stride,
                          &recbase[recbase_offset_y],&recbaseU[recbase_offset_c_y], &recbaseV[recbase_offset_c_y], recbase_stride,
                          &pred[pred_offset_y],&predU[pred_offset_c_y>>1],&predV[pred_offset_c_y>>1],pred_stride,
                          &coeff[coeff_fourth<<1],&coeffU[coeff_fourth],&coeffV[coeff_fourth],
                          depth+1, split);
    output |= outhelper<<6;
    outhelper = encode_transform_tree(encoder,&base[base_offset_y+half_width], &baseU[base_offset_c_y+(half_width>>1)], &baseV[base_offset_c_y+(half_width>>1)], base_stride,
                          &recbase[recbase_offset_y+half_width],&recbaseU[recbase_offset_c_y+(half_width>>1)], &recbaseV[recbase_offset_c_y+(half_width>>1)], recbase_stride,
                          &pred[pred_offset_y+half_width],&predU[pred_offset_c_y+(half_width>>1)],&predV[pred_offset_c_y+(half_width>>1)],pred_stride,
                          &coeff[3*coeff_fourth],&coeffU[(3*coeff_fourth)>>1],&coeffV[(3*coeff_fourth)>>1],
                          depth+1, split);
    output |= outhelper<<9;
    return output;
  }

  
  {
    
      
    /*
      Quant and transform here...
    */
    int16_t block[LCU_WIDTH*LCU_WIDTH>>2];
    int16_t pre_quant_coeff[LCU_WIDTH*LCU_WIDTH>>2];

    /* Get residual by subtracting prediction */
    i = 0;          
    for(y = 0; y < LCU_WIDTH>>depth; y++)
    {
      for(x = 0; x < LCU_WIDTH>>depth; x++)
      {
        block[i++]=((int16_t)base[x+y*encoder->in.width])-pred[x+y*pred_stride];
      }
    }

    /* Transform and quant residual to coeffs */          
    transform2d(block,pre_quant_coeff,width,0);
    quant(encoder,pre_quant_coeff,coeff,width, width,0, 0, SCAN_DIAG);

    /* Check for non-zero coeffs */
    for(i = 0; i < width*width; i++)
    {
      if(coeff[i] != 0)
      {
        /* Found one, we can break here */
        CbY = 1;
        break;
      }
    }
        
    /* if non-zero coeffs */
    if(CbY)
    {
      /* RECONSTRUCT for predictions */
      dequant(encoder,coeff,pre_quant_coeff,width, width,0);
      itransform2d(block,pre_quant_coeff,width,0);

      i = 0;
      for(y = 0; y < LCU_WIDTH>>depth; y++)
      {
        for(x = 0; x < LCU_WIDTH>>depth; x++)
        {
          int16_t val = block[i++]+pred[x+y*pred_stride];
          //ToDo: support 10+bits
          recbase[x+y*encoder->in.width] = (uint8_t)CLIP(0,255,val);
        }
      }
      /* END RECONTRUCTION */
    }
    /* without coeffs, we only use the prediction */
    else
    {
      for(y = 0; y < LCU_WIDTH>>depth; y++)
      {
        for(x = 0; x < LCU_WIDTH>>depth; x++)
        {
          recbase[x+y*encoder->in.width] = (uint8_t)CLIP(0,255,pred[x+y*pred_stride]);
        }
      }
    }

    if(encoder->in.video_format != FORMAT_400)
    {
      /* U */
      i = 0;
      for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
      {
        for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
        {
          block[i++]=((int16_t)baseU[x+y*(encoder->in.width>>1)])-predU[x+y*(pred_stride>>1)];
        }
      }
      transform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),0);
      quant(encoder,pre_quant_coeff,coeffU, width>>1, width>>1, 0,2,SCAN_DIAG);
      for(i = 0; i < width*width>>2; i++)
      {
        if(coeffU[i] != 0)
        {
          /* Found one, we can break here */
          CbU = 1;
          break;
        }
      }

      /* V */   
      i = 0;
      for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
      {
        for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
        {
          block[i++]=((int16_t)baseV[x+y*(encoder->in.width>>1)])-predV[x+y*(pred_stride>>1)];
        }
      }
      transform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),0);
      quant(encoder,pre_quant_coeff,coeffV, width>>1, width>>1, 0,3,SCAN_DIAG);
      for(i = 0; i < width*width>>2; i++)
      {
        if(coeffV[i] != 0)
        {
          /* Found one, we can break here */
          CbV = 1;
          break;
        }
      }
          
      if(CbU)
      {
        /* RECONSTRUCT for predictions */
        dequant(encoder,coeffU,pre_quant_coeff,width>>1, width>>1,2);
        itransform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),0);

        i = 0;
        for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
        {
          for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
          {
            int16_t val = block[i++]+predU[x+y*(pred_stride>>1)];
            //ToDo: support 10+bits
            recbaseU[x+y*(encoder->in.width>>1)] = (uint8_t)CLIP(0,255,val);
          }
        }
        /* END RECONTRUCTION */
      }
      /* without coeffs, we only use the prediction */
      else
      {
        for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
        {
          for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
          {
            recbaseU[x+y*(encoder->in.width>>1)] = (uint8_t)CLIP(0,255,predU[x+y*(pred_stride>>1)]);
          }
        }
      }
        
      if(CbV)
      {
        /* RECONSTRUCT for predictions */
        dequant(encoder,coeffV,pre_quant_coeff,width>>1, width>>1,3);
        itransform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),0);

        i = 0;
        for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
        {
          for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
          {
            int16_t val = block[i++]+predV[x+y*(pred_stride>>1)];
            //ToDo: support 10+bits
            recbaseV[x+y*(encoder->in.width>>1)] = (uint8_t)CLIP(0,255,val);
          }
        }
        /* END RECONTRUCTION */
      }
      /* without coeffs, we only use the prediction */
      else
      {
        for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
        {
          for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
          {
            recbaseV[x+y*(encoder->in.width>>1)] = (uint8_t)CLIP(0,255,predV[x+y*(pred_stride>>1)]);
          }
        }
      }
    }
        
    /* END INTRAPREDICTION */
  }

    /* end Residual Coding */
  return CbY | (CbU<<1) | (CbV<<2);
}


void encode_transform_coeff(encoder_control* encoder, int16_t *coeff, int16_t *coeffU, int16_t *coeffV, 
                            int32_t CbY, int32_t CbU, int32_t CbV,int8_t depth, int8_t intraPredMode, int8_t intraPredModeChroma, int8_t split, int8_t toplevel)
{
  int8_t width = LCU_WIDTH>>depth;
  if(depth != 0 && depth != MAX_DEPTH)
  {
    cabac.ctx = &g_TransSubdivSCModel[5-(g_aucConvertToBit[LCU_WIDTH]+2-depth)];
    CABAC_BIN(&cabac,split,"TransformSubdivFlag");
  }

  /* Signal if chroma data is present */
  if(encoder->in.video_format != FORMAT_400)
  {
    /* Non-zero chroma U Tcoeffs */
    cabac.ctx = &g_QtCbfSCModelU[0];
    CABAC_BIN(&cabac,(CbU&0x2||toplevel)?1:0,"cbf_chroma_u");
    /* Non-zero chroma V Tcoeffs */
    /* NOTE: Using the same ctx as before */
    CABAC_BIN(&cabac,(CbV&0x4||toplevel)?1:0,"cbf_chroma_v");
  }
  
  if(split)
  {
    encode_transform_coeff(encoder, coeff, coeffU, coeffV,
                           CbY, CbU, CbV,depth+1,
                           intraPredMode, intraPredModeChroma, 0, 0);
    encode_transform_coeff(encoder, coeff, coeffU, coeffV,
                           CbY, CbU, CbV,depth+1,
                           intraPredMode, intraPredModeChroma, 0, 0);
    encode_transform_coeff(encoder, coeff, coeffU, coeffV,
                           CbY, CbU, CbV,depth+1,
                           intraPredMode, intraPredModeChroma, 0, 0);
    encode_transform_coeff(encoder, coeff, coeffU, coeffV,
                           CbY, CbU, CbV,depth+1,
                           intraPredMode, intraPredModeChroma, 0, 0);
    return;
  }
  CbY = (CbY&0x1)?1:0;
  CbU = (CbU&0x2)?1:0;
  CbV = (CbV&0x4)?1:0;

  /* Non-zero luma Tcoeffs */
  cabac.ctx = &g_QtCbfSCModelY[1];
  CABAC_BIN(&cabac,CbY,"cbf_luma");

  {
    uint32_t uiCTXIdx;
    uint32_t uiScanIdx = SCAN_DIAG;
    uint32_t uiDirMode;
    switch(width)
    {
      case  2: uiCTXIdx = 6; break;
      case  4: uiCTXIdx = 5; break;
      case  8: uiCTXIdx = 4; break;
      case 16: uiCTXIdx = 3; break;
      case 32: uiCTXIdx = 2; break;
      case 64: uiCTXIdx = 1; break;
      default: uiCTXIdx = 0; break;
    }
    /* CoeffNxN */
    /* Residual Coding */
    if(CbY)
    {
      /* Luma (Intra) scanmode */
      uiDirMode = intraPredMode;
      if (uiCTXIdx >3 && uiCTXIdx < 6) //if multiple scans supported for transform size
      {
        uiScanIdx = abs((int32_t) uiDirMode - 26) < 5 ? 1 : (abs((int32_t)uiDirMode - 10) < 5 ? 2 : 0);
      }
      encode_CoeffNxN(encoder,coeff, width, 0, uiScanIdx);
    }
    if(CbU||CbV)
    {
      /* Chroma scanmode */
      uiCTXIdx++;
      uiDirMode = intraPredModeChroma;
      if(uiDirMode==36)
      {
        /* ToDo: support NxN */
        uiDirMode = intraPredMode;
      }
      uiScanIdx = SCAN_DIAG;
      if (uiCTXIdx >4 && uiCTXIdx < 7) //if multiple scans supported for transform size
      {
        uiScanIdx = abs((int32_t) uiDirMode - 26) < 5 ? 1 : (abs((int32_t)uiDirMode - 10) < 5 ? 2 : 0);
      }

      if(CbU)
      {
        encode_CoeffNxN(encoder,coeffU, width>>1, 2, uiScanIdx);
      }
      if(CbV)
      {
        encode_CoeffNxN(encoder,coeffV, width>>1, 2, uiScanIdx);
      }
    }
  }
}

void encode_CoeffNxN(encoder_control* encoder,int16_t* coeff, uint8_t width, uint8_t type, int8_t scanMode)
{
  int c1 = 1;
  uint8_t last_coeff_x = 0;
  uint8_t last_coeff_y = 0;
  int32_t i;
  uint32_t sig_coeffgroup_flag[64];
  
  uint32_t num_nonzero = 0;
  int32_t scanPosLast  = -1;
  int32_t posLast = 0;
  int32_t shift   = 4>>1;
  int8_t beValid  = ENABLE_SIGN_HIDING;
  int32_t iScanPosSig;
  int32_t iLastScanSet;
  uint32_t uiGoRiceParam = 0;
  uint32_t uiBlkPos, uiPosY, uiPosX, uiSig, uiCtxSig;  

  /* CONSTANTS */
  const uint32_t uiNumBlkSide    = width >> shift;
  const uint32_t uiLog2BlockSize = g_aucConvertToBit[ width ] + 2;
  const uint32_t* scan           = g_auiSigLastScan[ scanMode ][ uiLog2BlockSize - 1 ];
  const uint32_t* scanCG         = NULL;

  /* Init base contexts according to block type */
  cabac_ctx* baseCoeffGroupCtx = &g_CUSigCoeffGroupSCModel[type];
  cabac_ctx* baseCtx           = (type==0) ? &g_CUSigSCModel_luma[0] :&g_CUSigSCModel_chroma[0];
  memset(sig_coeffgroup_flag,0,sizeof(uint32_t)*64);
  
  /* Count non-zero coeffs */
  for(i = 0; i < width*width; i++)
  {
    if(coeff[i] != 0)
    {
      num_nonzero++;
    }
  }

  scanCG = g_auiSigLastScan[ scanMode ][ uiLog2BlockSize > 3 ? uiLog2BlockSize-2-1 : 0 ];
  if( uiLog2BlockSize == 3 )
  {
    scanCG = g_sigLastScan8x8[ scanMode ];
  }
  else if( uiLog2BlockSize == 5 )
  {
    scanCG = g_sigLastScanCG32x32;
  }

  scanPosLast = -1;
  /* Significance mapping */
  while(num_nonzero > 0)
  {
    posLast = scan[ ++scanPosLast ];
    #define POSY (posLast >> uiLog2BlockSize)
    #define POSX (posLast - ( POSY << uiLog2BlockSize ))
    if( coeff[ posLast ] != 0 )
    {
      sig_coeffgroup_flag[(uiNumBlkSide * (POSY >> shift) + (POSX >> shift))] = 1;
    }
    num_nonzero -= ( coeff[ posLast ] != 0 )?1:0;
    #undef POSY
    #undef POSX
  }
          
  last_coeff_x = posLast & (width-1);
  last_coeff_y = posLast>> uiLog2BlockSize;

  /* Code last_coeff_x and last_coeff_y */
  encode_lastSignificantXY(encoder,last_coeff_x, last_coeff_y, width, width, type, 0);
          
  iScanPosSig  = scanPosLast;
  iLastScanSet = (scanPosLast >> 4);
  /* significant_coeff_flag */
  for(i = iLastScanSet; i >= 0; i-- )
  {
    int32_t iSubPos       = i << 4 /*LOG2_SCAN_SET_SIZE*/;
    int32_t abs_coeff[16];
    int32_t iCGBlkPos     = scanCG[ i ];
    int32_t iCGPosY       = iCGBlkPos / uiNumBlkSide;
    int32_t iCGPosX       = iCGBlkPos - (iCGPosY * uiNumBlkSide);
    uint32_t coeffSigns   = 0;
    int32_t lastNZPosInCG = -1, firstNZPosInCG = 16;
    int32_t numNonZero    = 0;
    uiGoRiceParam         = 0;

    if( iScanPosSig == scanPosLast )
    {
      abs_coeff[ 0 ] = abs( coeff[ posLast ] );
      coeffSigns     = ( coeff[ posLast ] < 0 )?1:0;
      numNonZero     = 1;
      lastNZPosInCG  = iScanPosSig;
      firstNZPosInCG = iScanPosSig;
      iScanPosSig--;
    }
    if( i == iLastScanSet || i == 0)
    {
      sig_coeffgroup_flag[ iCGBlkPos ] = 1;
    }
    else
    {
      uint32_t uiSigCoeffGroup   = (sig_coeffgroup_flag[ iCGBlkPos ] != 0);
      uint32_t uiCtxSig  = context_get_sigCoeffGroup(sig_coeffgroup_flag, iCGPosX, iCGPosY,width);
      cabac.ctx = &baseCoeffGroupCtx[ uiCtxSig ];
      CABAC_BIN(&cabac,uiSigCoeffGroup,"significant_coeff_group");
    }

    if( sig_coeffgroup_flag[ iCGBlkPos ] )
    {
      int32_t patternSigCtx = context_calcPatternSigCtx( sig_coeffgroup_flag, iCGPosX, iCGPosY, width);
      for( ; iScanPosSig >= iSubPos; iScanPosSig-- )
      {
        uiBlkPos = scan[ iScanPosSig ]; 
        uiPosY   = uiBlkPos >> uiLog2BlockSize;
        uiPosX   = uiBlkPos - ( uiPosY << uiLog2BlockSize );
        uiSig    = (coeff[ uiBlkPos ] != 0)?1:0;
        if( iScanPosSig > iSubPos || i == 0 || numNonZero )
        {
          uiCtxSig  = context_getSigCtxInc( patternSigCtx, scanMode, uiPosX, uiPosY, uiLog2BlockSize, width, type );
          cabac.ctx = &baseCtx[ uiCtxSig ];
          CABAC_BIN(&cabac,uiSig,"significant_coeff_flag");
        }
                
        if( uiSig )
        {
          abs_coeff[ numNonZero ] = abs( coeff[ uiBlkPos ] );
          coeffSigns              = 2 * coeffSigns + ( coeff[ uiBlkPos ] < 0 );
          numNonZero++;
          if( lastNZPosInCG == -1 )
          {
            lastNZPosInCG = iScanPosSig;
          }
          firstNZPosInCG  = iScanPosSig;
        }
      }
    }
    else
    {
      iScanPosSig = iSubPos - 1;
    }            

    if( numNonZero > 0 )
    {
      uint8_t signHidden = ( lastNZPosInCG - firstNZPosInCG >= 4 /*SBH_THRESHOLD*/ );
      uint32_t uiCtxSet  = (i > 0 && type==0) ? 2 : 0;
      cabac_ctx* baseCtxMod;
      int32_t numC1Flag,firstC2FlagIdx,idx,iFirstCoeff2;
      if( c1 == 0 )
      {
        uiCtxSet++;
      }
      c1 = 1;

      baseCtxMod     = ( type==0 ) ? &g_CUOneSCModel_luma[4 * uiCtxSet] : &g_CUOneSCModel_chroma[4 * uiCtxSet];      
      numC1Flag      = MIN(numNonZero, C1FLAG_NUMBER);
      firstC2FlagIdx = -1;
      for(idx = 0; idx < numC1Flag; idx++ )
      {
        uint32_t uiSymbol = (abs_coeff[ idx ] > 1)?1:0;
        cabac.ctx = &baseCtxMod[c1];
        CABAC_BIN(&cabac,uiSymbol,"significant_coeff2_flag");
        if( uiSymbol )
        {
          c1 = 0;
          if (firstC2FlagIdx == -1)
          {
            firstC2FlagIdx = idx;
          }
        }
        else if( (c1 < 3) && (c1 > 0) )
        {
          c1++;
        }
      }
      
      if (c1 == 0)
      {
        baseCtxMod = ( type==0 ) ? &g_cCUAbsSCModel_luma[uiCtxSet] : &g_cCUAbsSCModel_chroma[uiCtxSet];
        if ( firstC2FlagIdx != -1)
        {
          uint8_t symbol = (abs_coeff[ firstC2FlagIdx ] > 2)?1:0;
          cabac.ctx      = &baseCtxMod[0];
          CABAC_BIN(&cabac,symbol,"first_c2_flag");
        }
      }      
      
      if( beValid && signHidden )
      {
        CABAC_BINS_EP(&cabac,(coeffSigns >> 1),(numNonZero-1),"");
      }
      else
      {
        CABAC_BINS_EP(&cabac,coeffSigns,numNonZero,"");
      }
              
      if (c1 == 0 || numNonZero > C1FLAG_NUMBER)
      {
        iFirstCoeff2 = 1;
        for (idx = 0; idx < numNonZero; idx++ )
        {
          int32_t baseLevel  = (idx < C1FLAG_NUMBER)? (2 + iFirstCoeff2 ) : 1;

          if( abs_coeff[ idx ] >= baseLevel)
          {
            cabac_writeCoeffRemain(&cabac, abs_coeff[ idx ] - baseLevel, uiGoRiceParam );
            if(abs_coeff[idx] > 3*(1<<uiGoRiceParam))
            {
                uiGoRiceParam = MIN(uiGoRiceParam+1, 4);
            }
          }
          if(abs_coeff[ idx ] >= 2)
          {
            iFirstCoeff2 = 0;
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
void encode_lastSignificantXY(encoder_control* encoder,uint8_t lastpos_x, uint8_t lastpos_y, uint8_t width, uint8_t height, uint8_t type, uint8_t scan)
{
  uint8_t offset_x  = type?0:((TOBITS(width)*3) + ((TOBITS(width)+1)>>2)),offset_y = offset_x;
  uint8_t shift_x   = type?(TOBITS(width)):((TOBITS(width)+3)>>2), shift_y = shift_x;
  int uiGroupIdxX;
  int uiGroupIdxY;
  int last_x,last_y,i;
  cabac_ctx* basectxX = (type?g_CuCtxLastX_chroma:g_CuCtxLastX_luma);
  cabac_ctx* basectxY = (type?g_CuCtxLastY_chroma:g_CuCtxLastY_luma);

  if( scan == SCAN_VER )
  {
    SWAP( lastpos_x, lastpos_y,uint8_t );
  }
  uiGroupIdxX   = g_uiGroupIdx[lastpos_x];
  uiGroupIdxY   = g_uiGroupIdx[lastpos_y];

  /* Last X binarization */
  for(last_x = 0; last_x < uiGroupIdxX ; last_x++)
  {
    cabac.ctx = &basectxX[offset_x+(last_x>>shift_x)];
    CABAC_BIN(&cabac,1,"LastSignificantX");
  }
  if(uiGroupIdxX < g_uiGroupIdx[width-1])
  {
    cabac.ctx = &basectxX[offset_x+(last_x>>shift_x)];
    CABAC_BIN(&cabac,0,"LastSignificantX");
  }

  /* Last Y binarization */
  for(last_y = 0; last_y < uiGroupIdxY ; last_y++)
  {
    cabac.ctx = &basectxY[offset_y+(last_y>>shift_y)];
    CABAC_BIN(&cabac,1,"LastSignificantY");
  }
  if(uiGroupIdxY < g_uiGroupIdx[height-1])
  {
    cabac.ctx = &basectxY[offset_y+(last_y>>shift_y)];
    CABAC_BIN(&cabac,0,"LastSignificantY");
  }

  /* Last X */
  if(uiGroupIdxX > 3)
  {
    lastpos_x -= g_uiMinInGroup[uiGroupIdxX];
    for(i = ((uiGroupIdxX-2)>>1)-1; i>=0; i--) 
    {
      CABAC_BIN_EP(&cabac,(lastpos_x>>i) & 1,"LastSignificantX");
    }
  }          
  /* Last Y */
  if(uiGroupIdxY > 3)
  {
    lastpos_y -= g_uiMinInGroup[uiGroupIdxY];
    for(i = ((uiGroupIdxY-2)>>1)-1; i>=0; i--) 
    {
      CABAC_BIN_EP(&cabac,(lastpos_y>>i) & 1,"LastSignificantY");
    }
  }
  /* end LastSignificantXY */
}