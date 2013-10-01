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
      pBuffTemp = (uint32_t *)g_sig_last_scan_32x32;
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
    uint32_t log2Blk      = g_convert_to_bit[ uiNumBlkSide ] + 1;

    for(uiBlk = 0; uiBlk < uiNumBlks; uiBlk++ )
    {      
      uint32_t initBlkPos = g_sig_last_scan[ SCAN_DIAG ][ log2Blk ][ uiBlk ];
      uiNextScanPos   = 0;
      if( iWidth == 32 )
      {
        initBlkPos = g_sig_last_scan_32x32[ uiBlk ];
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
          //TODO: optimize
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
        uint32_t offset = blkY * 4 * iWidth + blkX * 4;
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
  memset( g_convert_to_bit,-1, sizeof( g_convert_to_bit ) );  
  for ( i=4; i<(1<<7); i*=2 )
  {
    g_convert_to_bit[i] = c;
    c++;
  }
  g_convert_to_bit[i] = c;

  c = 2;
  for ( i=0; i<7; i++ )
  {
    g_sig_last_scan[0][i] = (uint32_t*)malloc(c*c*sizeof(uint32_t));
    g_sig_last_scan[1][i] = (uint32_t*)malloc(c*c*sizeof(uint32_t));
    g_sig_last_scan[2][i] = (uint32_t*)malloc(c*c*sizeof(uint32_t));

    initSigLastScan( g_sig_last_scan[0][i], g_sig_last_scan[1][i], g_sig_last_scan[2][i], c, c);
    c <<= 1;
  }

  /* Lambda cost */
  /* TODO: cleanup */
  //g_lambda_cost = (int16_t*)malloc(sizeof(int16_t)*55);
  for(i = 0; i < 55; i++)
  {
    if(i < 12) g_lambda_cost[i]= 0;
    else g_lambda_cost[i] = (int16_t)sqrt(0.57*pow(2.0,(i-12)/3));
    //g_lambda_cost[i] = g_lambda_cost[i]*g_lambda_cost[i];
  }

}
void init_encoder_control(encoder_control* control,bitstream* output)
{
  control->stream = output;  
}

void init_encoder_input(encoder_input* input,FILE* inputfile, int32_t width, int32_t height)
{
  input->file = inputfile;
  input->width = width;
  input->height = height;
  input->real_width = width;
  input->real_height = height;
  
  // If input dimensions are not divisible by the smallest block size, add pixels to the dimensions, so that they are.
  // These extra pixels will be compressed along with the real ones but they will be cropped out before rendering.
  if (width % CU_MIN_SIZE_PIXELS) {
    input->width += CU_MIN_SIZE_PIXELS - (width % CU_MIN_SIZE_PIXELS);
  }
  if (height % CU_MIN_SIZE_PIXELS) {
    input->height += CU_MIN_SIZE_PIXELS - (height % CU_MIN_SIZE_PIXELS);
  }

  input->height_in_lcu = input->height / LCU_WIDTH;
  input->width_in_lcu  = input->width / LCU_WIDTH;

  /* Add one extra LCU when image not divisible by LCU_WIDTH */
  if(input->height_in_lcu * LCU_WIDTH < height)
  {
    input->height_in_lcu++;
  }
  if(input->width_in_lcu * LCU_WIDTH < width)
  {
    input->width_in_lcu++;
  }

  /* Allocate the picture and CU array */
  input->cur_pic = picture_init(input->width, input->height, input->width_in_lcu,input->height_in_lcu);

  if(!input->cur_pic)
  {
    printf("Error allocating picture!\r\n");
    exit(1);
  }

  #ifdef _DEBUG
  if (width != input->width || height != input->height) {
    printf("Picture buffer has been extended to be a multiple of the smallest block size:\r\n");
    printf("  Width = %d (%d), Height = %d (%d)\r\n", width, input->width, height, input->height);
  }
  #endif



  }  

void encode_one_frame(encoder_control* encoder)
{
  /* output parameters before first frame */
  if(encoder->frame == 0)
  {

   /* Video Parameter Set (VPS) */    
    encode_vid_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, NAL_VPS_NUT, 0);
    bitstream_clear_buffer(encoder->stream);

    /* Sequence Parameter Set (SPS) */
    encode_seq_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, NAL_SPS_NUT, 0);
    bitstream_clear_buffer(encoder->stream);
        
    /* Picture Parameter Set (PPS) */
    encode_pic_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, NAL_PPS_NUT, 0);
    bitstream_clear_buffer(encoder->stream);

    /* First slice is IDR */
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
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, NAL_IDR_W_RADL, 0);
    bitstream_clear_buffer(encoder->stream);
  }
  
  //else if(encoder->frame == 1)
  //{
  //  /*
  //  cabac_start(&cabac);
  //  encoder->in.cur_pic->slicetype = SLICE_P;
  //  encoder->in.cur_pic->type = 1;

  //  encode_slice_header(encoder);
  //  bitstream_align(encoder->stream);
  //  cabac_flush(&cabac);
  //  bitstream_align(encoder->stream);
  //  bitstream_flush(encoder->stream);
  //  nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0,1,0);
  //  bitstream_clear_buffer(encoder->stream);*/


  //  cabac_start(&cabac);
  //  encoder->in.cur_pic->slicetype = SLICE_P;
  //  encoder->in.cur_pic->type = 1;
  //  search_slice_data(encoder);

  //  encode_slice_header(encoder);
  //  bitstream_align(encoder->stream);
  //  encode_slice_data(encoder);
  //  cabac_flush(&cabac);
  //  bitstream_align(encoder->stream);
  //  bitstream_flush(encoder->stream);
  //  nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0,1,1);
  //  bitstream_clear_buffer(encoder->stream);
  //}
  else
  {    
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
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0,NAL_TRAIL_R,0);
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
        printf("(%d,%d) Intramode: %d\n", x<<2, y<<2,encoder->in.cur_pic->CU[0][i].intra.mode);
      }
    }
  }
  */
  #endif
  
  /* Filtering */
  /* TODO: Check for correct deblock condition on inter blocks */
  if(encoder->deblock_enable && encoder->in.cur_pic->slicetype == SLICE_I)
  {
    filter_deblock(encoder);
  }


  /* Calculate checksum */
  add_checksum(encoder);

  }

void fill_after_frame(FILE* file, unsigned height, unsigned array_width, unsigned array_height, unsigned char* data)
{
  unsigned char* p = data + height * array_width;
  unsigned char* end = data + array_width * array_height;

  while (p < end) {
    // Fill the line by copying the line above.
    memcpy(p, p - array_width, array_width);
    p += array_width;
  }
}

void read_and_fill_frame_data(FILE* file, unsigned width, unsigned height, unsigned array_width, unsigned char* data)
{
  unsigned char* p = data;
  unsigned char* end = data + array_width * height;
  unsigned char fill_char;
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
    read_and_fill_frame_data(file, width, height, array_width, in->cur_pic->y_data);
    read_and_fill_frame_data(file, width >> 1, height >> 1, array_width >> 1, in->cur_pic->u_data);
    read_and_fill_frame_data(file, width >> 1, height >> 1, array_width >> 1, in->cur_pic->v_data);
  } else {
    // Otherwise the data can be read directly to the array.
    fread(in->cur_pic->y_data, sizeof(unsigned char), width * height, file);
    fread(in->cur_pic->u_data, sizeof(unsigned char), (width >> 1) * (height >> 1), file);
    fread(in->cur_pic->v_data, sizeof(unsigned char), (width >> 1) * (height >> 1), file);
  }

  if (height != array_height) {
    fill_after_frame(file, height, array_width, array_height, in->cur_pic->y_data);
    fill_after_frame(file, height >> 1, array_width >> 1, array_height >> 1, in->cur_pic->u_data);
    fill_after_frame(file, height >> 1, array_width >> 1, array_height >> 1, in->cur_pic->v_data);
  }
}

/*!
 \brief Add a checksum SEI message to the bitstream.
 \param encoder The encoder.
 \returns Void
 */
static void add_checksum(encoder_control* encoder)
{
  unsigned char checksum[3][SEI_HASH_MAX_LENGTH];
  uint32_t checksum_val;
  unsigned int i;

  picture_checksum(encoder->in.cur_pic, checksum);

  WRITE_U(encoder->stream, 132, 8, "sei_type");
  WRITE_U(encoder->stream, 13, 8, "size");
  WRITE_U(encoder->stream, 2, 8, "hash_type"); /* 2 = checksum*/

  for (i = 0; i < 3; ++i) {
    /* Pack bits into a single 32 bit uint instead of pushing them one byte at a time. */
    checksum_val = (checksum[i][0] << 24) + (checksum[i][1] << 16) + (checksum[i][2] << 8) + (checksum[i][3]);
    WRITE_U(encoder->stream, checksum_val, 32, "picture_checksum");
  }

  bitstream_align(encoder->stream);
  bitstream_flush(encoder->stream);
  nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, NAL_SUFFIT_SEI_NUT, 0);
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
    WRITE_U(encoder->stream, encoder->deblock_enable?0:1, 1, "pps_disable_deblocking_filter_flag");
    //IF !disabled
    if(encoder->deblock_enable)
    {
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
  /*PTL*/
  /*Profile Tier*/
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
  encoder_input* const in = &encoder->in;

#ifdef _DEBUG
  printf("=========== Sequence Parameter Set ID: 0 ===========\n");
#endif
  /* TODO: profile IDC and level IDC should be defined later on */
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

  if (in->width != in->real_width || in->height != in->real_height) {
    // The standard does not seem to allow setting conf_win values such that the number of
    // luma samples is not a multiple of 2. Options are to either hide one line or show an
    // extra line of non-video. Neither seems like a very good option, so let's not even try.
    assert(!(in->width % 2));
    WRITE_U(encoder->stream, 1, 1, "conformance_window_flag");
    WRITE_UE(encoder->stream, 0, "conf_win_left_offset");
    WRITE_UE(encoder->stream, (in->width - in->real_width) >> 1, "conf_win_right_offset");
    WRITE_UE(encoder->stream, 0, "conf_win_top_offset");
    WRITE_UE(encoder->stream, (in->height - in->real_height) >> 1, "conf_win_bottom_offset");
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
  WRITE_UE(encoder->stream, 0, "log2_min_transform_block_size_minus2");   /* 4x4 */
  WRITE_UE(encoder->stream, 3, "log2_diff_max_min_transform_block_size"); /* 4x4...32x32 */
  WRITE_UE(encoder->stream, 2, "max_transform_hierarchy_depth_inter");
  WRITE_UE(encoder->stream, 2, "max_transform_hierarchy_depth_intra");
  
  /* Use default scaling list */
  WRITE_U(encoder->stream, ENABLE_SCALING_LIST, 1, "scaling_list_enable_flag");  
  #if ENABLE_SCALING_LIST == 1
    WRITE_U(encoder->stream, 0, 1, "sps_scaling_list_data_present_flag");
  #endif
  
  WRITE_U(encoder->stream, 0, 1, "amp_enabled_flag");
  WRITE_U(encoder->stream, encoder->sao_enable?1:0, 1, "sample_adaptive_offset_enabled_flag");

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
  WRITE_U(encoder->stream, ENABLE_TEMPORAL_MVP, 1, "sps_temporal_mvp_enable_flag");

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
  if(encoder->in.cur_pic->type >= NAL_BLA_W_LP && encoder->in.cur_pic->type <= NAL_RSV_IRAP_VCL23)
  {
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
    if(encoder->in.cur_pic->type != NAL_IDR_W_RADL && encoder->in.cur_pic->type != NAL_IDR_N_LP)
    {
      int j;
      int ref_negative = 1;
      int ref_positive = 0;
      WRITE_U(encoder->stream, encoder->frame&0xf, 4, "pic_order_cnt_lsb");
      WRITE_U(encoder->stream, 0, 1, "short_term_ref_pic_set_sps_flag");
      WRITE_UE(encoder->stream, ref_negative, "num_negative_pics");
      WRITE_UE(encoder->stream, ref_positive, "num_positive_pics");
      for(j=0; j <ref_negative; j++)
      {
        WRITE_UE(encoder->stream, 0, "delta_poc_s0_minus1");
        WRITE_U(encoder->stream,1,1, "used_by_curr_pic_s0_flag");
      }

      //WRITE_UE(encoder->stream, 0, "short_term_ref_pic_set_idx");
    }
    //end if
  //end if
    if(encoder->sao_enable)
    {    
      WRITE_U(encoder->stream, 1,1, "slice_sao_luma_flag");
      WRITE_U(encoder->stream, 0,1, "slice_sao_chroma_flag");
    }
    
    if(encoder->in.cur_pic->slicetype != SLICE_I)
    {
      WRITE_U(encoder->stream, 0, 1, "num_ref_idx_active_override_flag");
      WRITE_UE(encoder->stream, 0, "five_minus_max_num_merge_cand");
    }

    if(encoder->in.cur_pic->slicetype == SLICE_B)
    {
      WRITE_U(encoder->stream, 0, 1, "mvd_l1_zero_flag");
    }
  /* Skip flags that are not present */
  // if !entropy_slice_flag
    WRITE_SE(encoder->stream, 0, "slice_qp_delta");
    //WRITE_U(encoder->stream, 1, 1, "alignment");
}

void encode_slice_data(encoder_control* encoder)
{
  uint16_t xCtb,yCtb;

  scalinglist_process();
  init_contexts(encoder,encoder->in.cur_pic->slicetype);

  /* Loop through every LCU in the slice */
  for(yCtb = 0; yCtb < encoder->in.height_in_lcu; yCtb++)
  {
    uint8_t lastCUy = (yCtb == (encoder->in.height_in_lcu-1))?1:0;
    for(xCtb = 0; xCtb < encoder->in.width_in_lcu; xCtb++)
    {
      uint8_t lastCUx = (xCtb == (encoder->in.width_in_lcu-1))?1:0;
      uint8_t depth = 0;

      /* Recursive function for looping through all the sub-blocks */
      encode_coding_tree(encoder, xCtb<<MAX_DEPTH,yCtb<<MAX_DEPTH, depth);

      /* signal Terminating bit */
      if(!lastCUx || !lastCUy)
      {
        cabac_encode_bin_trm(&cabac, 0);
      }
    }
  }
}

void encode_coding_tree(encoder_control* encoder,uint16_t xCtb,uint16_t yCtb, uint8_t depth)
{ 
  cu_info *cur_CU = &encoder->in.cur_pic->cu_array[MAX_DEPTH][xCtb+yCtb*(encoder->in.width_in_lcu<<MAX_DEPTH)];
  uint8_t split_flag = GET_SPLITDATA(cur_CU, depth);
  uint8_t split_model = 0;

  /* Check for slice border */
  uint8_t border_x = ((encoder->in.width)<( xCtb*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>depth) ))?1:0;
  uint8_t border_y = ((encoder->in.height)<( yCtb*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>depth) ))?1:0;
  uint8_t border_split_x = ((encoder->in.width)  < ( (xCtb+1)*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>(depth+1)) ))?0:1;
  uint8_t border_split_y = ((encoder->in.height) < ( (yCtb+1)*(LCU_WIDTH>>MAX_DEPTH) + (LCU_WIDTH>>(depth+1)) ))?0:1;
  uint8_t border = border_x | border_y; /*!< are we in any border CU */
  

  /* When not in MAX_DEPTH, insert split flag and split the blocks if needed */
  if(depth != MAX_DEPTH)
  {
    /* Implisit split flag when on border */
    if(!border)
    {
      /* Get left and top block split_flags and if they are present and true, increase model number */
      if(xCtb > 0 && GET_SPLITDATA(&(encoder->in.cur_pic->cu_array[MAX_DEPTH][xCtb-1+yCtb*(encoder->in.width_in_lcu<<MAX_DEPTH)]),depth) == 1)
      {
        split_model++;
      }
      if(yCtb > 0 && GET_SPLITDATA(&(encoder->in.cur_pic->cu_array[MAX_DEPTH][xCtb+(yCtb-1)*(encoder->in.width_in_lcu<<MAX_DEPTH)]),depth) == 1)
      {
        split_model++;
      }
      cabac.ctx = &g_split_flag_model[split_model];
      CABAC_BIN(&cabac, split_flag, "SplitFlag");
    }
    if(split_flag || border)
    {
      /* Split blocks and remember to change x and y block positions */
      uint8_t change = 1<<(MAX_DEPTH-1-depth);
      encode_coding_tree(encoder,xCtb,yCtb,depth+1); /* x,y */
      if(!border_x || border_split_x) /* TODO: fix when other half of the block would not be completely over the border */
      {
        encode_coding_tree(encoder,xCtb+change,yCtb,depth+1); /* x+1,y */
      }
      if(!border_y || border_split_y) /* TODO: fix when other half of the block would not be completely over the border */
      {
        encode_coding_tree(encoder,xCtb,yCtb+change,depth+1); /* x,y+1 */
      }
      if(!border || (border_split_x && border_split_y) ) /* TODO: fix when other half of the block would not be completely over the border */
      {
        encode_coding_tree(encoder,xCtb+change,yCtb+change,depth+1); /* x+1,y+1 */
      }
      /* We don't need to do anything else here */
      return;
    }
  }
  
  /* Encode skip flag */
  if(encoder->in.cur_pic->slicetype != SLICE_I)
  {
    int8_t uiCtxSkip = 0;    
    /* uiCtxSkip = aboveskipped + leftskipped; */
    cabac.ctx = &g_cu_skip_flag_model[uiCtxSkip];
    CABAC_BIN(&cabac, (cur_CU->type == CU_SKIP)?1:0, "SkipFlag");
  }

  /* IF SKIP */
  if(cur_CU->type == CU_SKIP)
  {
    /* Encode merge index */
    //TODO: calculate/fetch merge candidates
    int16_t unaryIdx = 0;//pcCU->getMergeIndex( uiAbsPartIdx );
    int16_t numCand = 0;//pcCU->getSlice()->getMaxNumMergeCand();
    int32_t ui;
    if ( numCand > 1 )
    {
      for(ui = 0; ui < numCand - 1; ui++ )
      {
        int32_t symbol = (ui == unaryIdx) ? 0 : 1;
        if ( ui==0 )
        {
          cabac.ctx = &g_cu_merge_idx_ext_model;
          CABAC_BIN(&cabac, symbol, "MergeIndex");
        }
        else
        {
          CABAC_BIN_EP(&cabac,symbol,"MergeIndex");
        }
        if( symbol == 0 )
        {
          break;
        }
      }
    }
    return;
  }
  /* ENDIF SKIP */

  /* Prediction mode */
  if(encoder->in.cur_pic->slicetype != SLICE_I)
  {
    cabac.ctx = &g_cu_pred_mode_model;
    CABAC_BIN(&cabac, (cur_CU->type == CU_INTRA)?1:0, "PredMode");
  }

  /* Signal PartSize on max depth */
  if(depth == MAX_DEPTH || cur_CU->type != CU_INTRA)
  {
    /* TODO: Handle inter sizes other than 2Nx2N */
    cabac.ctx = &g_part_size_model[0];
    CABAC_BIN(&cabac, 1, "PartSize");
    /* TODO: add AMP modes */
  }
    
    /*end partsize*/
    if(cur_CU->type == CU_INTER)
    {
      /* FOR each part */
        /* Mergeflag */
        uint8_t mergeFlag = 0;
        cabac.ctx = &g_cu_merge_flag_ext_model;
        CABAC_BIN(&cabac, mergeFlag, "MergeFlag");
        if(mergeFlag) //merge
        {
          /* MergeIndex */
          int16_t unaryIdx = 0;//pcCU->getMergeIndex( uiAbsPartIdx );
          int16_t numCand  = 0;//pcCU->getSlice()->getMaxNumMergeCand();
          int32_t ui;
          if (numCand > 1)
          {
            for(ui = 0; ui < numCand - 1; ui++ )
            {
              int32_t symbol = (ui == unaryIdx) ? 0 : 1;
              if (ui == 0)
              {
                cabac.ctx = &g_cu_merge_idx_ext_model;
                CABAC_BIN(&cabac, symbol, "MergeIndex");
              }
              else
              {
                CABAC_BIN_EP(&cabac,symbol,"MergeIndex");
              }
              if(symbol == 0)
              {
                break;
              }
            }
          }
        }
        else
        {
          uint32_t uiRefListIdx;

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

          
          
          for(uiRefListIdx = 0; uiRefListIdx < 2; uiRefListIdx++)
          {
            //if(encoder->ref_idx_num[uiRefListIdx] > 0)
            {
              if(cur_CU->inter.mv_dir & (1 << uiRefListIdx))
              {
                if(0)//encoder->ref_idx_num[uiRefListIdx] != 1)//NumRefIdx != 1)
                {
                  /* parseRefFrmIdx */
                  int32_t iRefFrame = cur_CU->inter.mv_ref;
                  
                  cabac.ctx = &g_cu_ref_pic_model[0];
                  CABAC_BIN(&cabac, (iRefFrame==0)?0:1, "ref_frame_flag"); 
    
                  if(iRefFrame > 0)
                  {
                    uint32_t ui;
                    uint32_t uiRefNum = encoder->ref_idx_num[uiRefListIdx]-2;

                    cabac.ctx = &g_cu_ref_pic_model[1];
                    iRefFrame--;
                    for(ui = 0; ui < uiRefNum; ++ui)
                    {
                      const uint32_t uiSymbol = (ui==iRefFrame)?0:1;
                      if(ui == 0)
                      {
                        CABAC_BIN(&cabac, uiSymbol, "ref_frame_flag2");
                      }
                      else
                      {
                        CABAC_BIN_EP(&cabac,uiSymbol,"ref_frame_flag2");
                      }
                      if(uiSymbol == 0)
                      {
                        break;
                      }
                    }
                  }
                }

                /* Get MV candidates */
                inter_get_mv_cand(encoder, xCtb, yCtb, depth, mv_cand);

                /* Select better candidate */
                cur_CU->inter.mv_ref = 0; /* Default to candidate 0 */
                /* Only check when candidates are different */
                if (mv_cand[0][0] != mv_cand[1][0] || mv_cand[0][1] != mv_cand[1][1]) {
                  uint16_t cand_1_diff = abs(cur_CU->inter.mv[0]-mv_cand[0][0]) + abs(cur_CU->inter.mv[1]-mv_cand[0][1]);
                  uint16_t cand_2_diff = abs(cur_CU->inter.mv[0]-mv_cand[1][0]) + abs(cur_CU->inter.mv[1]-mv_cand[1][1]);
                  /* Select candidate 1 if it's closer */
                  if (cand_2_diff < cand_1_diff) {
                    cur_CU->inter.mv_ref = 1;
                  }
                }

                if (!(/*pcCU->getSlice()->getMvdL1ZeroFlag() &&*/ encoder->ref_list == REF_PIC_LIST_1 && cur_CU->inter.mv_dir==3))
                {
                    const int32_t mvd_hor = cur_CU->inter.mv[0]-mv_cand[cur_CU->inter.mv_ref][0];
                    const int32_t mvd_ver = cur_CU->inter.mv[1]-mv_cand[cur_CU->inter.mv_ref][1];
                    const int8_t bHorAbsGr0 = mvd_hor != 0;
                    const int8_t bVerAbsGr0 = mvd_ver != 0;
                    const uint32_t mvd_hor_abs = abs(mvd_hor);
                    const uint32_t mvd_ver_abs = abs(mvd_ver);

                    cabac.ctx = &g_cu_mvd_model[0];
                    CABAC_BIN(&cabac, (mvd_hor!=0)?1:0, "abs_mvd_greater0_flag_hor");
                    CABAC_BIN(&cabac, (mvd_ver!=0)?1:0, "abs_mvd_greater0_flag_ver");

                    cabac.ctx = &g_cu_mvd_model[1];

                    if(bHorAbsGr0)
                    {
                      CABAC_BIN(&cabac, (mvd_hor_abs>1)?1:0, "abs_mvd_greater1_flag_hor");
                    }

                    if(bVerAbsGr0)
                    {
                      CABAC_BIN(&cabac, (mvd_ver_abs>1)?1:0, "abs_mvd_greater1_flag_ver");
                    }

                    if(bHorAbsGr0)
                    {
                      if(mvd_hor_abs > 1)
                      {
                        cabac_write_ep_ex_golomb(&cabac,mvd_hor_abs-2, 1);
                      }
                      CABAC_BIN_EP(&cabac, (mvd_hor>0)?0:1, "mvd_sign_flag_hor");
                    }

                    if(bVerAbsGr0)
                    {
                      if(mvd_ver_abs > 1)
                      {
                        cabac_write_ep_ex_golomb(&cabac,mvd_ver_abs-2, 1);
                      }
                      CABAC_BIN_EP(&cabac, (mvd_ver>0)?0:1, "mvd_sign_flag_ver");
                    }

                    /* Inter reconstruction */
                    inter_recon(encoder->ref->pics[0],xCtb*CU_MIN_SIZE_PIXELS,yCtb*CU_MIN_SIZE_PIXELS,LCU_WIDTH>>depth,cur_CU->inter.mv,encoder->in.cur_pic);

                    /* Mark this block as "coded" (can be used for predictions..) */
                    picture_set_block_coded(encoder->in.cur_pic,xCtb, yCtb, depth, 1);
                }
                /* Signal which candidate MV to use */
                cabac_write_unary_max_symbol(&cabac,g_mvp_idx_model, cur_CU->inter.mv_ref,1,AMVP_MAX_NUM_CANDS-1);
              }
            }
          }

          cabac.ctx = &g_cu_qt_root_cbf_model;
          CABAC_BIN(&cabac, 1, "rqt_root_cbf");
          if(1)
          {
            uint8_t *base  = &encoder->in.cur_pic->y_data[xCtb*(LCU_WIDTH>>(MAX_DEPTH))   + (yCtb*(LCU_WIDTH>>(MAX_DEPTH)))  *encoder->in.width];
            uint8_t *baseU = &encoder->in.cur_pic->u_data[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
            uint8_t *baseV = &encoder->in.cur_pic->v_data[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
            uint32_t width = LCU_WIDTH>>depth;

            /* INTRAPREDICTION VARIABLES */
            int16_t pred[LCU_WIDTH*LCU_WIDTH+1];
            int16_t predU[LCU_WIDTH*LCU_WIDTH>>2];
            int16_t predV[LCU_WIDTH*LCU_WIDTH>>2];

            uint8_t *recbase   = &encoder->in.cur_pic->y_recdata[xCtb*(LCU_WIDTH>>(MAX_DEPTH))   + (yCtb*(LCU_WIDTH>>(MAX_DEPTH)))  *encoder->in.width];
            uint8_t *recbaseU  = &encoder->in.cur_pic->u_recdata[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
            uint8_t *recbaseV  = &encoder->in.cur_pic->v_recdata[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];

            /* TODO: dynamic memory allocation */
            int16_t coeff[LCU_WIDTH*LCU_WIDTH*2];
            int16_t coeffU[LCU_WIDTH*LCU_WIDTH>>1];
            int16_t coeffV[LCU_WIDTH*LCU_WIDTH>>1];

            /* Initialize helper structure for transform */
            transform_info ti;
            memset(&ti, 0, sizeof(transform_info));

            ti.x_ctb = xCtb; ti.y_ctb = yCtb;

            /* Base pointers */
            ti.base =  base; ti.base_u = baseU; ti.base_v = baseV;
            ti.base_stride = encoder->in.width;

            // Prediction pointers
            ti.pred =  pred; ti.pred_u = predU; ti.pred_v = predV;
            ti.pred_stride = (LCU_WIDTH>>depth);

            // Reconstruction pointers
            ti.recbase = recbase; ti.recbase_u = recbaseU; ti.recbase_v = recbaseV;
            ti.recbase_stride = encoder->in.width;

            // Coeff pointers
            ti.coeff[0] = coeff; ti.coeff[1] = coeffU; ti.coeff[2] = coeffV;
            ti.block_type = CU_INTER;
       
            // Handle transforms, quant and reconstruction
            ti.idx = 0;
            encode_transform_tree(encoder,&ti, depth);

            // Coded block pattern
            ti.cb_top[0] = (ti.cb[0] & 0x1 || ti.cb[1] & 0x1 || ti.cb[2] & 0x1 || ti.cb[3] & 0x1)?1:0;
            ti.cb_top[1] = (ti.cb[0] & 0x2 || ti.cb[1] & 0x2 || ti.cb[2] & 0x2 || ti.cb[3] & 0x2)?1:0;
            ti.cb_top[2] = (ti.cb[0] & 0x4 || ti.cb[1] & 0x4 || ti.cb[2] & 0x4 || ti.cb[3] & 0x4)?1:0;
        
            // Code (possible) coeffs to bitstream
            ti.idx = 0;
            encode_transform_coeff(encoder, &ti,depth, 0);
          }
        }

      /* END for each part */
    }
    else if(cur_CU->type == CU_INTRA)
    {
      uint8_t intraPredMode = cur_CU->intra.mode;
      uint8_t intraPredModeChroma = 36; /* 36 = Chroma derived from luma */
      int8_t intraPreds[3] = {-1, -1, -1};
      int8_t mpmPred = -1;
      int i;
      uint32_t flag;
      uint8_t *base  = &encoder->in.cur_pic->y_data[xCtb*(LCU_WIDTH>>(MAX_DEPTH))   + (yCtb*(LCU_WIDTH>>(MAX_DEPTH)))  *encoder->in.width];
      uint8_t *baseU = &encoder->in.cur_pic->u_data[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
      uint8_t *baseV = &encoder->in.cur_pic->v_data[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
      uint32_t width = LCU_WIDTH>>depth;

      /* INTRAPREDICTION VARIABLES */
      int16_t pred[LCU_WIDTH*LCU_WIDTH+1];
      int16_t predU[LCU_WIDTH*LCU_WIDTH>>2];
      int16_t predV[LCU_WIDTH*LCU_WIDTH>>2];

      uint8_t *recbase   = &encoder->in.cur_pic->y_recdata[xCtb*(LCU_WIDTH>>(MAX_DEPTH))   + (yCtb*(LCU_WIDTH>>(MAX_DEPTH)))  *encoder->in.width];
      uint8_t *recbaseU  = &encoder->in.cur_pic->u_recdata[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];
      uint8_t *recbaseV  = &encoder->in.cur_pic->v_recdata[xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*(encoder->in.width>>1)];


      /* SEARCH BEST INTRA MODE (AGAIN) */  
      
      int16_t rec[(LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)];
      int16_t *recShift = &rec[(LCU_WIDTH>>(depth))*2+8+1];      
      intra_build_reference_border(encoder->in.cur_pic, xCtb, yCtb,(LCU_WIDTH>>(depth))*2+8, rec, (LCU_WIDTH>>(depth))*2+8, 0);
      cur_CU->intra.mode = (int8_t)intra_prediction(encoder->in.cur_pic->y_data,encoder->in.width,recShift,(LCU_WIDTH>>(depth))*2+8,xCtb*(LCU_WIDTH>>(MAX_DEPTH)),yCtb*(LCU_WIDTH>>(MAX_DEPTH)),width,pred,width,&cur_CU->intra.cost);
      intraPredMode = cur_CU->intra.mode;
      intra_set_block_mode(encoder->in.cur_pic,xCtb, yCtb, depth, intraPredMode);
      
      #if ENABLE_PCM == 1
      /* Code must start after variable initialization */
      cabac_encode_bin_trm(&cabac, 0); /* IPCMFlag == 0 */
      #endif
      
      
      /*
        PREDINFO CODING
        If intra prediction mode is found from the predictors,
        it can be signaled with two EP's. Otherwise we can send
        5 EP bins with the full predmode
        TODO: split to a function
      */
      intra_get_dir_luma_predictor(encoder->in.cur_pic, xCtb, yCtb, depth, intraPreds);
      
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
      cabac.ctx = &g_intra_mode_model;
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
        int32_t intraPredModeTemp = intraPredMode;
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
        cabac.ctx = &g_chroma_pred_model[0];
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
        /* TODO: dynamic memory allocation */
        int16_t coeff[LCU_WIDTH*LCU_WIDTH*2];
        int16_t coeffU[LCU_WIDTH*LCU_WIDTH>>1];
        int16_t coeffV[LCU_WIDTH*LCU_WIDTH>>1];

        /* Initialize helper structure for transform */
        transform_info ti;
        memset(&ti, 0, sizeof(transform_info));

        ti.x_ctb = xCtb; ti.y_ctb = yCtb;

        /* Base pointers */
        ti.base =  base; ti.base_u = baseU; ti.base_v = baseV;
        ti.base_stride = encoder->in.width;

        /* Prediction pointers */
        ti.pred =  pred; ti.pred_u = predU; ti.pred_v = predV;
        ti.pred_stride = (LCU_WIDTH>>depth);

        /* Reconstruction pointers */
        ti.recbase = recbase; ti.recbase_u = recbaseU; ti.recbase_v = recbaseV;
        ti.recbase_stride = encoder->in.width;

        /* Coeff pointers */
        ti.coeff[0] = coeff; ti.coeff[1] = coeffU; ti.coeff[2] = coeffV;

        /* Prediction info */
        ti.intra_pred_mode = intraPredMode; ti.intra_pred_mode_chroma = intraPredModeChroma;
        
        /* Handle transforms, quant and reconstruction */
        ti.idx = 0;
        ti.block_type = CU_INTRA;
        encode_transform_tree(encoder,&ti, depth);

        /* Coded block pattern */
        ti.cb_top[0] = (ti.cb[0] & 0x1 || ti.cb[1] & 0x1 || ti.cb[2] & 0x1 || ti.cb[3] & 0x1)?1:0;
        ti.cb_top[1] = (ti.cb[0] & 0x2 || ti.cb[1] & 0x2 || ti.cb[2] & 0x2 || ti.cb[3] & 0x2)?1:0;
        ti.cb_top[2] = (ti.cb[0] & 0x4 || ti.cb[1] & 0x4 || ti.cb[2] & 0x4 || ti.cb[3] & 0x4)?1:0;
        
        /* Code (possible) coeffs to bitstream */
        ti.idx = 0;
        encode_transform_coeff(encoder, &ti,depth, 0);
      }
      /* end Transform tree */
      /* end Coeff */

    }
    #if ENABLE_PCM == 1
    /* Code IPCM block */
    else if(cur_CU->type == CU_PCM)
    {
      cabac_encode_bin_trm(&cabac, 1); /* IPCMFlag == 1 */
      cabac_finish(&cabac);
      bitstream_align(cabac.stream);
       /* PCM sample */
      {
        uint8_t *base   = &encoder->in.cur_pic->y_data[xCtb*(LCU_WIDTH>>(MAX_DEPTH))    + (yCtb*(LCU_WIDTH>>(MAX_DEPTH)))*encoder->in.width];
        uint8_t *baseCb = &encoder->in.cur_pic->u_data[(xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*encoder->in.width/2)];
        uint8_t *baseCr = &encoder->in.cur_pic->v_data[(xCtb*(LCU_WIDTH>>(MAX_DEPTH+1)) + (yCtb*(LCU_WIDTH>>(MAX_DEPTH+1)))*encoder->in.width/2)];
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

void encode_transform_tree(encoder_control* encoder,transform_info* ti,uint8_t depth)
{
  /* we have 64>>depth transform size */
  int x,y,i;
  int32_t width = LCU_WIDTH>>depth;

  if(depth == 0) /* Split 64x64 */
  {
    /* Prepare for multi-level splitting */
    ti->split[ti->idx] = 1<<depth;
  }

  /* Split transform and increase depth */
  if(ti->split[ti->idx] & (1<<depth))
  {
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
    uint8_t CbY = 0,CbU = 0,CbV = 0;
    int32_t coeff_fourth = ((LCU_WIDTH>>(depth))*(LCU_WIDTH>>(depth)))+1;

    int32_t base_stride    = ti->base_stride;
    int32_t recbase_stride = ti->recbase_stride;
    int32_t pred_stride    = ti->pred_stride;

    int32_t recbase_offset[4]   = {0, width   , ti->recbase_stride*(width)        , ti->recbase_stride*(width)        +width     };
    int32_t base_offset[4]      = {0, width   , ti->base_stride*(width)           , ti->base_stride*(width)           +width     };
    int32_t pred_offset[4]      = {0, width   , ti->pred_stride*(width)           , ti->pred_stride*(width)           +width     };
    int32_t recbase_offset_c[4] = {0, width>>1, (ti->recbase_stride>>1)*(width>>1), (ti->recbase_stride>>1)*(width>>1)+(width>>1)};
    int32_t base_offset_c[4]    = {0, width>>1, (ti->base_stride>>1)*(width>>1)   , (ti->base_stride>>1)*(width>>1)   +(width>>1)};
    int32_t pred_offset_c[4]    = {0, width>>1, (ti->pred_stride>>1)*(width>>1)   , (ti->pred_stride>>1)*(width>>1)   +(width>>1)};
    
    uint8_t* base     = &ti->base[base_offset[ti->idx]];
    uint8_t* baseU    = &ti->base_u[base_offset_c[ti->idx]];
    uint8_t* baseV    = &ti->base_v[base_offset_c[ti->idx]];
    
    uint8_t* recbase  = &ti->recbase[recbase_offset[ti->idx]];
    uint8_t* recbaseU = &ti->recbase_u[recbase_offset_c[ti->idx]];
    uint8_t* recbaseV = &ti->recbase_v[recbase_offset_c[ti->idx]];
    
    int16_t* pred     = &ti->pred[pred_offset[ti->idx]];
    int16_t* predU    = &ti->pred_u[pred_offset_c[ti->idx]];
    int16_t* predV    = &ti->pred_v[pred_offset_c[ti->idx]];
    
    int16_t* coeff    = &ti->coeff[0][ti->idx*coeff_fourth];
    int16_t* coeffU   = &ti->coeff[1][ti->idx*coeff_fourth>>1];
    int16_t* coeffV   = &ti->coeff[2][ti->idx*coeff_fourth>>1];
    
    /*
      Quant and transform here...
    */
    int16_t block[LCU_WIDTH*LCU_WIDTH>>2];
    int16_t pre_quant_coeff[LCU_WIDTH*LCU_WIDTH>>2];

    /* INTRA PREDICTION */
    /* TODO: split to a function! */
    int16_t rec[(LCU_WIDTH*2+8)*(LCU_WIDTH*2+8)];
    int16_t *recShift  = &rec[(LCU_WIDTH>>(depth))*2+8+1];
    int16_t *recShiftU = &rec[(LCU_WIDTH>>(depth+1))*2+8+1];

    uint32_t ac_sum = 0;
    if(ti->block_type == CU_INTRA)
    {
      /* Build reconstructed block to use in prediction with extrapolated borders */
      intra_build_reference_border(encoder->in.cur_pic, ti->x_ctb, ti->y_ctb,(LCU_WIDTH>>(depth))*2+8, rec, (LCU_WIDTH>>(depth))*2+8, 0);
      intra_recon(recShift,(LCU_WIDTH>>(depth))*2+8,ti->x_ctb*(LCU_WIDTH>>(MAX_DEPTH)),ti->y_ctb*(LCU_WIDTH>>(MAX_DEPTH)),width,pred,pred_stride,ti->intra_pred_mode,0);

      /* Filter DC-prediction */
      if(ti->intra_pred_mode == 1 && width < 32)
      {
        intra_dc_pred_filtering(recShift,(LCU_WIDTH>>(depth))*2+8,pred,width,LCU_WIDTH>>depth,LCU_WIDTH>>depth);
      }      
      if(ti->intra_pred_mode_chroma != 36 && ti->intra_pred_mode_chroma == ti->intra_pred_mode)
      {
        ti->intra_pred_mode_chroma = 36;
      }
      intra_build_reference_border(encoder->in.cur_pic, ti->x_ctb, ti->y_ctb,(LCU_WIDTH>>(depth+1))*2+8, rec, (LCU_WIDTH>>(depth+1))*2+8, 1);
      intra_recon(recShiftU,(LCU_WIDTH>>(depth+1))*2+8,ti->x_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)),ti->y_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)),width>>1,predU,pred_stride>>1,ti->intra_pred_mode_chroma!=36?ti->intra_pred_mode_chroma:ti->intra_pred_mode,1);
      intra_build_reference_border(encoder->in.cur_pic, ti->x_ctb, ti->y_ctb,(LCU_WIDTH>>(depth+1))*2+8, rec, (LCU_WIDTH>>(depth+1))*2+8, 2);
      intra_recon(recShiftU,(LCU_WIDTH>>(depth+1))*2+8,ti->x_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)),ti->y_ctb*(LCU_WIDTH>>(MAX_DEPTH+1)),width>>1,predV,pred_stride>>1,ti->intra_pred_mode_chroma!=36?ti->intra_pred_mode_chroma:ti->intra_pred_mode,1);
    
      /* This affects reconstruction, do after that */ 
      picture_set_block_coded(encoder->in.cur_pic, ti->x_ctb, ti->y_ctb, depth, 1);
    } else  { // Inter mode
      for(y = 0; y < LCU_WIDTH>>depth; y++) {
        for(x = 0; x < LCU_WIDTH>>depth; x++) {
          pred[x+y*pred_stride]=recbase[x+y*base_stride];
        }
      }
      for(y = 0; y < LCU_WIDTH>>(depth+1); y++) {
        for(x = 0; x < LCU_WIDTH>>(depth+1); x++) {
          predU[x+y*(pred_stride>>1)]=recbaseU[x+y*(base_stride>>1)];
          predV[x+y*(pred_stride>>1)]=recbaseV[x+y*(base_stride>>1)];
        }
      }

    }
    /* INTRA PREDICTION ENDS HERE */

    /* Get residual by subtracting prediction */
    i = 0;
    ac_sum = 0;
    for(y = 0; y < LCU_WIDTH>>depth; y++)
    {
      for(x = 0; x < LCU_WIDTH>>depth; x++)
      {
        block[i]=((int16_t)base[x+y*base_stride])-pred[x+y*pred_stride];
        i++;
      }
    }

    /* Transform and quant residual to coeffs */
    transform2d(block,pre_quant_coeff,width,0);
    quant(encoder,pre_quant_coeff,coeff,width, width,&ac_sum, 0, SCAN_DIAG);

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
          //TODO: support 10+bits
          recbase[x+y*recbase_stride] = (uint8_t)/*(val&0xff);//*/CLIP(0,255,val);
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
          recbase[x+y*recbase_stride] = (uint8_t)CLIP(0,255,pred[x+y*pred_stride]);
        }
      }
    }

    if(encoder->in.video_format != FORMAT_400)
    {
      /* U */
      i = 0;
      ac_sum = 0;
      for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
      {
        for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
        {
          block[i]=((int16_t)baseU[x+y*(base_stride>>1)])-predU[x+y*(pred_stride>>1)];
          i++;
        }
      }
      transform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),65535);
      quant(encoder,pre_quant_coeff,coeffU, width>>1, width>>1, &ac_sum,2,SCAN_DIAG);
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
      ac_sum = 0;
      for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
      {
        for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
        {
          block[i]=((int16_t)baseV[x+y*(base_stride>>1)])-predV[x+y*(pred_stride>>1)];
          i++;
        }
      }
      transform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),65535);
      quant(encoder,pre_quant_coeff,coeffV, width>>1, width>>1, &ac_sum,3,SCAN_DIAG);
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
        itransform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),65535);

        i = 0;
        for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
        {
          for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
          {
            int16_t val = block[i++]+predU[x+y*(pred_stride>>1)];
            //TODO: support 10+bits
            recbaseU[x+y*(recbase_stride>>1)] = (uint8_t)CLIP(0,255,val);
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
            recbaseU[x+y*(recbase_stride>>1)] = (uint8_t)CLIP(0,255,predU[x+y*(pred_stride>>1)]);
          }
        }
      }
      
      if(CbV)
      {
        /* RECONSTRUCT for predictions */
        dequant(encoder,coeffV,pre_quant_coeff,width>>1, width>>1,3);
        itransform2d(block,pre_quant_coeff,LCU_WIDTH>>(depth+1),65535);

        i = 0;
        for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
        {
          for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
          {
            int16_t val = block[i++]+predV[x+y*(pred_stride>>1)];
            //TODO: support 10+bits
            recbaseV[x+y*(recbase_stride>>1)] = (uint8_t)CLIP(0,255,val);
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
            recbaseV[x+y*(recbase_stride>>1)] = (uint8_t)CLIP(0,255,predV[x+y*(pred_stride>>1)]);
          }
        }
      }
    }
    
    /* Store coded block pattern */
    ti->cb[ti->idx] = CbY | (CbU<<1) | (CbV<<2);
    /* END INTRAPREDICTION */
    return;
  }

    /* end Residual Coding */
  
}


void encode_transform_coeff(encoder_control* encoder,transform_info* ti,int8_t depth, int8_t trDepth)
{
  int8_t width = LCU_WIDTH>>depth;
  int8_t split = ((ti->split[ti->idx]&(1<<depth))||!depth)?1:0;
  int8_t CbY,CbU,CbV;
  int32_t coeff_fourth = ((LCU_WIDTH>>(depth))*(LCU_WIDTH>>(depth)))+1;
  
  if(depth != 0 && depth != MAX_DEPTH+1)
  {
    cabac.ctx = &g_trans_subdiv_model[5-((g_convert_to_bit[LCU_WIDTH]+2)-depth)];
    CABAC_BIN(&cabac,split,"TransformSubdivFlag");
  }

  /* Signal if chroma data is present */
  /* Chroma data is also signaled BEFORE transform split */
  /* Chroma data is not signaled if it was set to 0 before split */
  if(encoder->in.video_format != FORMAT_400)
  {
    /* Non-zero chroma U Tcoeffs */    
    int8_t Cb_flag = (trDepth==0)?ti->cb_top[1]:((ti->cb[ti->idx]&0x2)?1:0);
    cabac.ctx = &g_qt_cbf_model_chroma[trDepth];
    if(trDepth == 0 || ti->cb_top[1])
    {
      CABAC_BIN(&cabac,Cb_flag,"cbf_chroma_u");
    }
    /* Non-zero chroma V Tcoeffs */
    /* NOTE: Using the same ctx as before */    
    Cb_flag = (trDepth==0)?ti->cb_top[2]:((ti->cb[ti->idx]&0x4)?1:0);
    if(trDepth == 0 || ti->cb_top[2])
    {
      CABAC_BIN(&cabac,Cb_flag,"cbf_chroma_v");
    }
  }
  
  if(split)
  {    
    ti->idx = 0; encode_transform_coeff(encoder,ti,depth+1,trDepth+1);
    ti->idx = 1; encode_transform_coeff(encoder,ti,depth+1,trDepth+1);
    ti->idx = 2; encode_transform_coeff(encoder,ti,depth+1,trDepth+1);
    ti->idx = 3; encode_transform_coeff(encoder,ti,depth+1,trDepth+1);
    return;
  }
  CbY = ti->cb[ti->idx]&0x1;
  CbU = (ti->cb[ti->idx]&0x2)?1:0;
  CbV = (ti->cb[ti->idx]&0x4)?1:0;

  if(ti->block_type == CU_INTRA || trDepth || CbU || CbV) {
    // Non-zero luma Tcoeffs
    cabac.ctx = &g_qt_cbf_model_luma[trDepth?0:1];
    CABAC_BIN(&cabac,CbY,"cbf_luma");
  }

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
    uiCTXIdx -= trDepth;
    /* CoeffNxN */
    /* Residual Coding */
    if(CbY)
    {
      if(ti->block_type == CU_INTER) {
        uiScanIdx = SCAN_DIAG;
      } else {      
        /* Luma (Intra) scanmode */
        uiDirMode = ti->intra_pred_mode;
        if (uiCTXIdx >3 && uiCTXIdx < 6) { //if multiple scans supported for transform size
          uiScanIdx = abs((int32_t) uiDirMode - 26) < 5 ? 1 : (abs((int32_t)uiDirMode - 10) < 5 ? 2 : 0);
        }
      }
      encode_CoeffNxN(encoder,&ti->coeff[0][ti->idx*coeff_fourth], width, 0, uiScanIdx);
    }
    if(CbU||CbV)
    {
      int8_t chromaWidth = width>>1;

      if(ti->block_type == CU_INTER) {
        uiScanIdx = SCAN_DIAG;
      } else {
        /* Chroma scanmode */
        uiCTXIdx++;
        uiDirMode = ti->intra_pred_mode_chroma;
        if(uiDirMode==36) {
          /* TODO: support NxN */
          uiDirMode = ti->intra_pred_mode;
        }
        uiScanIdx = SCAN_DIAG;
        if (uiCTXIdx >4 && uiCTXIdx < 7) { //if multiple scans supported for transform size
          uiScanIdx = abs((int32_t) uiDirMode - 26) < 5 ? 1 : (abs((int32_t)uiDirMode - 10) < 5 ? 2 : 0);
        }
      }

      if(CbU)
      {
        encode_CoeffNxN(encoder,&ti->coeff[1][ti->idx*coeff_fourth>>1], chromaWidth, 2, uiScanIdx);
      }
      if(CbV)
      {
        encode_CoeffNxN(encoder,&ti->coeff[2][ti->idx*coeff_fourth>>1], chromaWidth, 3, uiScanIdx);
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
  const uint32_t uiLog2BlockSize = g_convert_to_bit[ width ] + 2;
  const uint32_t* scan           = g_sig_last_scan[ scanMode ][ uiLog2BlockSize - 1 ];
  const uint32_t* scanCG         = NULL;

  /* Init base contexts according to block type */
  cabac_ctx* baseCoeffGroupCtx = &g_cu_sig_coeff_group_model[type];
  cabac_ctx* baseCtx           = (type==0) ? &g_cu_sig_model_luma[0] :&g_cu_sig_model_chroma[0];
  memset(sig_coeffgroup_flag,0,sizeof(uint32_t)*64);
  
  /* Count non-zero coeffs */
  for(i = 0; i < width*width; i++)
  {
    if(coeff[i] != 0)
    {
      num_nonzero++;
    }
  }  

  scanCG = g_sig_last_scan[ scanMode ][ uiLog2BlockSize > 3 ? uiLog2BlockSize-3 : 0 ];
  if( uiLog2BlockSize == 3 )
  {
    scanCG = g_sig_last_scan_8x8[ scanMode ];
  }
  else if( uiLog2BlockSize == 5 )
  {
    scanCG = g_sig_last_scan_32x32;
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
  encode_lastSignificantXY(encoder,last_coeff_x, last_coeff_y, width, width, type, scanMode);
  
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
      coeffSigns     = ( coeff[ posLast ] < 0 );
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
      uint32_t uiCtxSig  = context_get_sig_coeff_group(sig_coeffgroup_flag, iCGPosX, iCGPosY,width);
      cabac.ctx = &baseCoeffGroupCtx[ uiCtxSig ];
      CABAC_BIN(&cabac,uiSigCoeffGroup,"significant_coeff_group");
    }

    if( sig_coeffgroup_flag[ iCGBlkPos ] )
    {
      int32_t patternSigCtx = context_calc_pattern_sig_ctx( sig_coeffgroup_flag, iCGPosX, iCGPosY, width);
      for( ; iScanPosSig >= iSubPos; iScanPosSig-- )
      {
        uiBlkPos = scan[ iScanPosSig ]; 
        uiPosY   = uiBlkPos >> uiLog2BlockSize;
        uiPosX   = uiBlkPos - ( uiPosY << uiLog2BlockSize );
        uiSig    = (coeff[ uiBlkPos ] != 0)?1:0;
        if( iScanPosSig > iSubPos || i == 0 || numNonZero )
        {
          uiCtxSig  = context_get_sig_ctx_inc( patternSigCtx, scanMode, uiPosX, uiPosY, uiLog2BlockSize, width, type );
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
      int8_t signHidden = ( lastNZPosInCG - firstNZPosInCG >= 4 /*SBH_THRESHOLD*/ )?1:0;
      uint32_t uiCtxSet  = (i > 0 && type==0) ? 2 : 0;
      cabac_ctx* baseCtxMod;
      int32_t numC1Flag,firstC2FlagIdx,idx,iFirstCoeff2;
      if( c1 == 0 )
      {
        uiCtxSet++;
      }
      c1 = 1;

      baseCtxMod     = ( type==0 ) ? &g_cu_one_model_luma[4 * uiCtxSet] : &g_cu_one_model_chroma[4 * uiCtxSet];      
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
        baseCtxMod = ( type==0 ) ? &g_cu_abs_model_luma[uiCtxSet] : &g_cu_abs_model_chroma[uiCtxSet];
        if (firstC2FlagIdx != -1)
        {
          uint8_t symbol = (abs_coeff[ firstC2FlagIdx ] > 2)?1:0;
          cabac.ctx      = &baseCtxMod[0];
          CABAC_BIN(&cabac,symbol,"first_c2_flag");
        }
      }      
      
      if(beValid && signHidden)
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
            cabac_write_coeff_remain(&cabac, abs_coeff[ idx ] - baseLevel, uiGoRiceParam );
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
  cabac_ctx* basectxX = (type?g_cu_ctx_last_x_chroma:g_cu_ctx_last_x_luma);
  cabac_ctx* basectxY = (type?g_cu_ctx_last_y_chroma:g_cu_ctx_last_y_luma);

  if( scan == SCAN_VER )
  {
    SWAP( lastpos_x, lastpos_y,uint8_t );
  }
  uiGroupIdxX   = g_group_idx[lastpos_x];
  uiGroupIdxY   = g_group_idx[lastpos_y];

  /* Last X binarization */
  for(last_x = 0; last_x < uiGroupIdxX ; last_x++)
  {
    cabac.ctx = &basectxX[offset_x+(last_x>>shift_x)];
    CABAC_BIN(&cabac,1,"LastSignificantX");
  }
  if(uiGroupIdxX < g_group_idx[width-1])
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
  if(uiGroupIdxY < g_group_idx[height-1])
  {
    cabac.ctx = &basectxY[offset_y+(last_y>>shift_y)];
    CABAC_BIN(&cabac,0,"LastSignificantY");
  }

  /* Last X */
  if(uiGroupIdxX > 3)
  {
    lastpos_x -= g_min_in_group[uiGroupIdxX];
    for(i = ((uiGroupIdxX-2)>>1)-1; i>=0; i--) 
    {
      CABAC_BIN_EP(&cabac,(lastpos_x>>i) & 1,"LastSignificantX");
    }
  }          
  /* Last Y */
  if(uiGroupIdxY > 3)
  {
    lastpos_y -= g_min_in_group[uiGroupIdxY];
    for(i = ((uiGroupIdxY-2)>>1)-1; i>=0; i--) 
    {
      CABAC_BIN_EP(&cabac,(lastpos_y>>i) & 1,"LastSignificantY");
    }
  }
  /* end LastSignificantXY */
}