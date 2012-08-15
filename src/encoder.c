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
#include "context.h"

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
  input->width_in_LCU =  width / LCU_WIDTH;
  if(input->height_in_LCU * LCU_WIDTH < height)
    input->height_in_LCU++;
  if(input->width_in_LCU * LCU_WIDTH < width)
    input->width_in_LCU++;

  input->cur_pic.width = width;
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

  /* Allocate memory for CU info 2D array */
  //ToDo: we don't need this much space on LCU...MAX_DEPTH-1
  input->cur_pic.CU = (CU_info**)malloc((MAX_DEPTH+1)*sizeof(CU_info*));
  for(i=0; i < MAX_DEPTH+1; i++)
  {
    input->cur_pic.CU[i] = (CU_info*)malloc((input->height_in_LCU<<2)*(input->width_in_LCU<<2)*sizeof(CU_info));
    memset(input->cur_pic.CU[i], 0, (input->height_in_LCU<<2)*(input->width_in_LCU<<2)*sizeof(CU_info));
  }  
}


void encode_one_frame(encoder_control* encoder)
{
  /* output parameters before first frame */
  if(encoder->frame == 0)
  {
    /* Sequence Parameter Set (SPS) */
    encode_seq_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 1, NAL_SEQ_PARAMETER_SET, 1);
    bitstream_clear_buffer(encoder->stream);

    /* Video Parameter Set (VPS) */    
    encode_vid_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 1, NAL_VID_PARAMETER_SET, 0);
    bitstream_clear_buffer(encoder->stream);
    
    /* Picture Parameter Set (PPS) */
    encode_pic_parameter_set(encoder);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 1, NAL_PIC_PARAMETER_SET, 0);
    bitstream_clear_buffer(encoder->stream);

    /* First slice is IDR */
    cabac_start(&cabac);
    encoder->in.cur_pic.type = NAL_IDR_SLICE;
    encode_slice_header(encoder);
    bitstream_align(encoder->stream);
    encode_slice_data(encoder);
    cabac_flush(&cabac);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, NAL_IDR_SLICE, 1);
    bitstream_clear_buffer(encoder->stream);
  }
  else if(encoder->frame < 3)
  {
    /* Non-IDR slice */
    cabac_start(&cabac);
    encoder->in.cur_pic.type = NAL_NONIDR_SLICE;
    encode_slice_header(encoder);
    bitstream_align(encoder->stream);
    encode_slice_data(encoder);
    cabac_flush(&cabac);
    bitstream_align(encoder->stream);
    bitstream_flush(encoder->stream);
    nal_write(encoder->output, encoder->stream->buffer, encoder->stream->buffer_pos, 0, NAL_NONIDR_SLICE, encoder->frame+1);
    bitstream_clear_buffer(encoder->stream);
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
  WRITE_UE(encoder->stream, 0, "num_ref_idx_l0_default_active_minus1");
  WRITE_UE(encoder->stream, 0, "num_ref_idx_l1_default_active_minus1");
  WRITE_SE(encoder->stream, encoder->QP-26, "pic_init_qp_minus26");
  WRITE_U(encoder->stream, 0, 1, "constrained_intra_pred_flag");
  WRITE_U(encoder->stream, 0, 1, "transform_skip_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "cu_qp_delta_enabled_flag");
  //if cu_qp_delta_enabled_flag
  //WRITE_UE(encoder->stream, 0, "diff_cu_qp_delta_depth");

  WRITE_SE(encoder->stream, 0, "cb_qp_offset");
  WRITE_SE(encoder->stream, 0, "cr_qp_offset");
  WRITE_U(encoder->stream, 0, 1, "slicelevel_chroma_qp_flag");
  WRITE_U(encoder->stream, 0, 1, "weighted_pred_flag");
  WRITE_U(encoder->stream, 0, 1, "weighted_bipred_idc");
  WRITE_U(encoder->stream, 1, 1, "output_flag_present_flag");
  WRITE_U(encoder->stream, 0, 1, "dependent_slices_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "transquant_bypass_enable_flag");
  WRITE_U(encoder->stream, 0, 2, "tiles_or_entropy_coding_sync_idc");
  WRITE_U(encoder->stream, 0, 1, "loop_filter_across_slice_flag");
  WRITE_U(encoder->stream, 0, 1, "deblocking_filter_control_present_flag");
  WRITE_U(encoder->stream, 0, 1, "pps_scaling_list_data_present_flag");
  WRITE_UE(encoder->stream, 0, "log2_parallel_merge_level_minus2");
  WRITE_U(encoder->stream, 0, 1, "slice_header_extension_present_flag");
  WRITE_U(encoder->stream, 0, 1, "pps_extension_flag");
}

void encode_seq_parameter_set(encoder_control* encoder)
{
  int i;
#ifdef _DEBUG
  printf("=========== Sequence Parameter Set ID: 0 ===========\n");
#endif
  /* ToDo: profile IDC and level IDC should be defined later on */
  WRITE_U(encoder->stream, 0, 3, "profile_space");
  WRITE_U(encoder->stream, 0, 5, "profile_idc");
  WRITE_U(encoder->stream, 0, 16, "reserved_indicator_flags");
  WRITE_U(encoder->stream, 0, 8, "level_idc");
  WRITE_U(encoder->stream, 0, 32, "profile_compatibility");
  WRITE_UE(encoder->stream, 0, "seq_parameter_set_id");
  WRITE_UE(encoder->stream, 0, "video_parameter_set_id");
  WRITE_UE(encoder->stream, encoder->in.video_format, "chroma_format_idc"); /* 0 = 4:0:0, 1 = 4:2:0, 2 = 4:2:2, 3 = 4:4:4 */
  WRITE_U(encoder->stream, 0, 3, "max_temporal_layers_minus1");
  WRITE_UE(encoder->stream, encoder->in.width, "pic_width_in_luma_samples");
  WRITE_UE(encoder->stream, encoder->in.height, "pic_height_in_luma_samples");
  WRITE_U(encoder->stream, 0, 1, "pic_cropping_flag");
  /* ToDo: 10bit support? */
  WRITE_UE(encoder->stream, 0, "bit_depth_luma_minus8");
  WRITE_UE(encoder->stream, 0, "bit_depth_chroma_minus8");
  WRITE_U(encoder->stream, ENABLE_PCM, 1, "pcm_enabled_flag");
  #if ENABLE_PCM == 1
    WRITE_U(encoder->stream, 7, 4, "pcm_bit_depth_luma_minus1");
    WRITE_U(encoder->stream, 7, 4, "pcm_bit_depth_chroma_minus1");
  #endif
  WRITE_UE(encoder->stream, 4, "log2_max_pic_order_cnt_lsb_minus4");
  WRITE_UE(encoder->stream, 0, "max_dec_pic_buffering");
  WRITE_UE(encoder->stream, 0, "num_reorder_pics");
  WRITE_UE(encoder->stream, 0, "max_latency_increase");
  WRITE_U(encoder->stream, 0, 1, "restricted_ref_pic_lists_flag");
  WRITE_UE(encoder->stream, 1, "log2_min_coding_block_size_minus3");
  WRITE_UE(encoder->stream, MAX_DEPTH, "log2_diff_max_min_coding_block_size");
  WRITE_UE(encoder->stream, 0, "log2_min_transform_block_size_minus2");
  WRITE_UE(encoder->stream, 3, "log2_diff_max_min_transform_block_size");

  //If log2MinCUSize == 3
  //WRITE_U(encoder->stream, 0, 1, "DisInter4x4");

  #if ENABLE_PCM == 1
    WRITE_UE(encoder->stream, 0, "log2_min_pcm_coding_block_size_minus3");
    WRITE_UE(encoder->stream, 2, "log2_diff_max_min_pcm_coding_block_size");
  #endif
  

  WRITE_UE(encoder->stream, 2, "max_transform_hierarchy_depth_inter");
  WRITE_UE(encoder->stream, 2, "max_transform_hierarchy_depth_intra");

  WRITE_U(encoder->stream, 0, 1, "scaling_list_enable_flag");
  WRITE_U(encoder->stream, 0, 1, "asymmetric_motion_partitions_enabled_flag");
  WRITE_U(encoder->stream, 0, 1, "sample_adaptive_offset_enabled_flag");
	//WRITE_U(encoder->stream, 0, 1, "adaptive_loop_filter_enabled_flag");
  #if ENABLE_PCM == 1
    WRITE_U(encoder->stream, 1, 1, "pcm_loop_filter_disable_flag");
  #endif
  WRITE_U(encoder->stream, 0, 1, "temporal_id_nesting_flag");
  WRITE_UE(encoder->stream, 0, "num_short_term_ref_pic_sets");  
  //WRITE_U(encoder->stream, 0, 1, "inter_ref_pic_set_prediction_flag");
  WRITE_U(encoder->stream, 0, 1, "long_term_ref_pics_present_flag");
  WRITE_U(encoder->stream, 0, 1, "sps_temporal_mvp_enable_flag");
  for(i = 0; i < MAX_DEPTH; i++)
  {
    WRITE_U(encoder->stream, 0, 1, "AMVP modeflag");
  }
	WRITE_U(encoder->stream, 0, 1, "sps_extension_flag");
}

void encode_vid_parameter_set(encoder_control* encoder)
{
#ifdef _DEBUG
  printf("=========== Video Parameter Set ID: 0 ===========\n");
#endif
  WRITE_U(encoder->stream, 0, 3, "vps_max_temporal_layers_minus1");
  WRITE_U(encoder->stream, 0, 5, "vps_max_layers_minus1");
  WRITE_UE(encoder->stream, 0, "video_parameter_set_id");
  WRITE_U(encoder->stream, 0, 1, "vps_temporal_id_nesting_flag");

  WRITE_UE(encoder->stream, 0, "vps_max_dec_pic_buffering");
  WRITE_UE(encoder->stream, 0, "vps_num_reorder_pics");
  WRITE_UE(encoder->stream, 0, "vps_max_latency_increase");

	WRITE_U(encoder->stream, 0, 1, "vps_extension_flag");
}

void encode_slice_header(encoder_control* encoder)
{
#ifdef _DEBUG
  printf("=========== Slice ===========\n");
#endif

  WRITE_U(encoder->stream, 1, 1, "first_slice_in_pic_flag");
  if(encoder->in.cur_pic.type == NAL_IDR_SLICE)
  {
    WRITE_U(encoder->stream, 0, 1, "no_output_of_prior_pics_flag");
  }
  WRITE_UE(encoder->stream, 0, "pic_parameter_set_id");
  
  WRITE_UE(encoder->stream, SLICE_I, "slice_type");

  WRITE_U(encoder->stream, 0, 1, "dependent_slice_flag");

  // if !entropy_slice_flag
  
    //if output_flag_present_flag
      WRITE_U(encoder->stream, 1, 1, "pic_output_flag");
    //end if
    //if( IdrPicFlag ) <- nal_unit_type == 5
    if(encoder->in.cur_pic.type == NAL_IDR_SLICE)
    {
      //WRITE_UE(encoder->stream, encoder->frame&3, "idr_pic_id");      
    }
    else
    {
      WRITE_U(encoder->stream, encoder->frame+1, 8, "pic_order_cnt_lsb");
      WRITE_U(encoder->stream, 1, 1, "short_term_ref_pic_set_sps_flag");
      //WRITE_U(encoder->stream, 1, 1, "inter_ref_pic_set_prediction_flag");
      WRITE_UE(encoder->stream, 0, "short_term_ref_pic_set_idx");
    }
    //end if
  //end if
  /* Skip flags that are not present */
  // if !entropy_slice_flag
    WRITE_SE(encoder->stream, 0, "slice_qp_delta");
    WRITE_UE(encoder->stream, 0, "5_minus_max_num_merge_cand");

    WRITE_U(encoder->stream, 1, 1, "alignment");
}
  



void encode_slice_data(encoder_control* encoder)
{
  uint16_t xCtb,yCtb;

  init_contexts(encoder);

  //encoder->in.cur_pic.CU[1][2].type = CU_INTRA;
  //encoder->in.cur_pic.CU[1][3].type = CU_INTRA;
  
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
  int i,x,y;
  uint8_t split_flag = (depth!=1)?1:0;
  uint8_t split_model = 0;
  uint8_t border_x = (encoder->in.width+1)<((xCtb>>(MAX_DEPTH-depth))+1)*(LCU_WIDTH>>depth);
  uint8_t border_y = (encoder->in.height+1)<((yCtb>>(MAX_DEPTH-depth))+1)*(LCU_WIDTH>>depth);
  uint8_t border = border_x | border_y;

  /* Get left and top block split_flags and if they are present and true, increase model number */
  if(xCtb > 0 && GET_SPLITDATA(&(encoder->in.cur_pic.CU[depth][(xCtb>>(MAX_DEPTH-depth))-1+(yCtb>>(MAX_DEPTH-depth))*(encoder->in.width_in_LCU<<MAX_DEPTH)])) == 1)
  {
    split_model++;
  }
  if(yCtb > 0 && GET_SPLITDATA(&(encoder->in.cur_pic.CU[depth][(xCtb>>(MAX_DEPTH-depth))+((yCtb>>(MAX_DEPTH-depth))-1)*(encoder->in.width_in_LCU<<MAX_DEPTH)])) == 1)
  {
    split_model++;
  }
  
  /* When not in MAX_DEPTH, insert split flag and split the blocks if needed */
  if(depth != MAX_DEPTH || border)
  {
    SET_SPLITDATA(&(encoder->in.cur_pic.CU[depth][(xCtb>>(MAX_DEPTH-depth))+(yCtb>>(MAX_DEPTH-depth))*(encoder->in.width_in_LCU<<MAX_DEPTH)]),split_flag);
    cabac.ctx = &g_SplitFlagSCModel[split_model];
    //Implisit split flag when on border
    if(!border)
    {
      CABAC_BIN(&cabac, split_flag, "SplitFlag");
    }
    if(split_flag)
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
  /* coding_unit( x0, y0, log2CbSize ) */
   /* prediction_unit 2Nx2N*/
    //if !intra PREDMODE
    /* if depth = MAX_DEPTH */
     //PartSize
     if(depth == MAX_DEPTH)
     {
       cabac.ctx = &g_PartSizeSCModel;
       CABAC_BIN(&cabac, 1, "PartSize");
     }
   /*end partsize*/
   //If MODE_INTRA
    //cabac.ctx = &PCMFlagSCModel;
     /* Code IPCM block */
    if(encoder->in.cur_pic.CU[depth][(xCtb>>(MAX_DEPTH-depth))+(yCtb>>(MAX_DEPTH-depth))*(encoder->in.width_in_LCU<<MAX_DEPTH)].type <= CU_PCM)
    {
      cabac_encodeBinTrm(&cabac, 1); /* IPCMFlag == 1 */
      //printf("\tIPCMFlag = 1\n");
      cabac_finish(&cabac);
      WRITE_U(cabac.stream, 1, 1, "stop_bit");
      WRITE_U(cabac.stream, 0, 1, "numSubseqIPCM_flag");
      bitstream_align(cabac.stream);
       /* PCM sample */
      {
        uint8_t *base = &encoder->in.cur_pic.yData[xCtb*(LCU_WIDTH>>(depth+1)) + (yCtb*(LCU_WIDTH>>(depth+1)))*encoder->in.width];
        uint8_t *baseCb = &encoder->in.cur_pic.uData[(xCtb*(LCU_WIDTH>>(depth+2)) + (yCtb*(LCU_WIDTH>>(depth+2)))*encoder->in.width/2)];
        uint8_t *baseCr = &encoder->in.cur_pic.vData[(xCtb*(LCU_WIDTH>>(depth+2)) + (yCtb*(LCU_WIDTH>>(depth+2)))*encoder->in.width/2)];
        for(y = 0; y < LCU_WIDTH>>depth; y++)
        {
          for(x = 0; x < LCU_WIDTH>>depth; x++)
          {
            bitstream_put(cabac.stream, base[x+y*encoder->in.width], 8);
          }
        }
        if(encoder->in.video_format != FORMAT_400)
        {
          //Cb
          for(y = 0; y < LCU_WIDTH>>(depth+1); y++)
          {
            for(x = 0; x < LCU_WIDTH>>(depth+1); x++)
            {
              bitstream_put(cabac.stream, baseCb[x+y*(encoder->in.width>>1)], 8);
            }
          }

          //Cr
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
    else
    {
      cabac_encodeBinTrm(&cabac, 0); /* IPCMFlag == 0 */
       
      cabac.ctx = &g_IntraModeSCModel;
      CABAC_BIN(&cabac,0,"IntraPred");

      /*
      Int preds[3] = {-1, -1, -1};
      Int predNum = pcCU->getIntraDirLumaPredictor(absPartIdx+partOffset*j, preds);  
      */
      CABAC_BINS_EP(&cabac, 0, 5, "intraPredMode");
 
      if(encoder->in.video_format != FORMAT_400)
      {
        cabac.ctx = &g_ChromaPredSCModel[0];
        CABAC_BIN(&cabac,0,"IntraPredChroma");
      }

      /* Coeff */
      /* Transform tree */
      cabac.ctx = &g_TransSubdivSCModel[1]; /* //uiLog2TransformBlockSize */
      CABAC_BIN(&cabac,0,"TransformSubdivFlag");

      /* We don't subdiv and we have 64>>depth transform size */
      /* ToDo: allow other sized */
      {
        uint8_t CbY = 0,CbU = 0,CbV = 0;
      
        /*
         Quant and transform here...
        */
        CbY = 1; /* Let's pretend we have luma coefficients */

        if(encoder->in.video_format != FORMAT_400)
        {
          /* Non-zero chroma U Tcoeffs */
          cabac.ctx = &g_QtCbfSCModelU[0];
          CABAC_BIN(&cabac,CbU,"cbf_chroma_u");

          /* Non-zero chroma V Tcoeffs */
          /* Using the same ctx as before */
          CABAC_BIN(&cabac,CbV,"cbf_chroma_v");
        }

        /* Non-zero luma Tcoeffs */
        cabac.ctx = &g_QtCbfSCModelY[1];
        CABAC_BIN(&cabac,CbY,"cbf_luma");

        /* CoeffNxN */
        if(CbY)
        {
          int c1,c1_num;
          //int patternSigCtx;
          int numNonZero = 16;
          /* scanCG == g_sigLastScanCG32x32 */
          /* Residual Coding */
          /* LastSignificantXY */
          encode_lastSignificantXY(encoder,31/*last_coeff_x */, 31/* last_coeff_y */, 32, 32, 0, 0);
          
          for(i = 15; i >= 0; i-- )
          {
            /* significant_coeff_flag */
            cabac.ctx = &g_CUSigSCModel_luma[21+((i<7)?1:0)]; /* 21 = uiCtxSig =TComTrQuant::getSigCtxInc( patternSigCtx, uiPosX, uiPosY, blockType, uiWidth, uiHeight, eTType );*/
            CABAC_BIN(&cabac,1,"significant_coeff_flag");
          }

          /*UInt uiSigCoeffGroupFlag[ MLS_GRP_NUM ];
            ::memset( uiSigCoeffGroupFlag, 0, sizeof(UInt) * MLS_GRP_NUM );*/
          /* n = 15 .. 0 coeff_abs_level_greater1_flag[ n ] */

          /* coeff_abs_level_greater2_flag[ firstGreater1CoeffIdx] */

          c1 = 1;
          c1_num = MIN(numNonZero,C1FLAG_NUMBER);
          for(i = 0; i < c1_num; i++)
          {
            int val = 0;
            /* significant_coeff_flag */
            cabac.ctx = &g_CUOneSCModel_luma[4*2+c1]; /* 4 * uiCtxSet +c1 */
            CABAC_BIN(&cabac,val,"coeff_abs_level_greater1_flag");
            if(val)
            {
              c1 = 0;
            }
            else if(c1 < 3)
            {
              c1++;
            }
          }

          /* end Residual Coding */
        }
      }
      /* end Transform tree */
      /* end Coeff */
    }
   //endif
   /* end prediction unit */

   //cabac_encodeBin(&cabac, 0); //prev_intra_luma_pred_flag

   //cabac_encodeBin(&cabac, 1); //rem_intra_luma_pred_mode

  /* end coding_unit */
  
}


void encode_lastSignificantXY(encoder_control* encoder,uint8_t lastpos_x, uint8_t lastpos_y, uint8_t width, uint8_t height, uint8_t type, uint8_t scan)
{
  uint8_t offset_x  = type?0:((TOBITS(width)*3) + ((TOBITS(width)+1)>>2)),offset_y = offset_x;
  uint8_t shift_x   = type?(TOBITS(width)):((TOBITS(width)+3)>>2), shift_y = shift_x;
  int uiGroupIdxX   = g_uiGroupIdx[lastpos_x];
  int uiGroupIdxY   = g_uiGroupIdx[lastpos_y];
  int last_x,last_y,i;

  if(width != height)
  {
    shift_y = (TOBITS(height)+3)>>2;
    offset_y = TOBITS(height)*3 + ((TOBITS(height)+1)>>2);
  }

  /* Last X binarization */
  for(last_x = 0; last_x < uiGroupIdxX ; last_x++)
  {
    cabac.ctx = &g_CuCtxLastX_luma[offset_x+(last_x>>shift_x)];
    CABAC_BIN(&cabac,1,"LastSignificantX");
  }
  if(uiGroupIdxX < g_uiGroupIdx[32-1])
  {
    cabac.ctx = &g_CuCtxLastX_luma[offset_x+(last_x>>shift_x)];
    CABAC_BIN(&cabac,0,"LastSignificantX");
  }

  /* Last Y binarization */
  for(last_y = 0; last_y < uiGroupIdxY ; last_y++)
  {
    cabac.ctx = &g_CuCtxLastY_luma[offset_y+(last_y>>shift_y)];
    CABAC_BIN(&cabac,1,"LastSignificantY");
  }
  if(uiGroupIdxY < g_uiGroupIdx[32-1])
  {
    cabac.ctx = &g_CuCtxLastY_luma[offset_y+(last_y>>shift_y)];
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