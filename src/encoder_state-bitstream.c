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

#include "encoder_state-bitstream.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bitstream.h"
#include "cabac.h"
#include "checkpoint.h"
#include "cu.h"
#include "encoder.h"
#include "encoder_state-geometry.h"
#include "encoderstate.h"
#include "imagelist.h"
#include "kvazaar.h"
#include "kvz_math.h"
#include "nal.h"
#include "scalinglist.h"
#include "tables.h"
#include "threadqueue.h"
#include "videoframe.h"


static void encoder_state_write_bitstream_aud(encoder_state_t * const state)
{
  bitstream_t * const stream = &state->stream;
  // ***********************************************
  // Modified for SHVC. TODO: only in base layer?
  kvz_nal_ext_write(stream, KVZ_NAL_AUD_NUT, 0, 1, state->layer->layer_id);
  // ***********************************************

  uint8_t pic_type = state->frame->slicetype == KVZ_SLICE_I ? 0
                   : state->frame->slicetype == KVZ_SLICE_P ? 1
                   :                                       2;
  WRITE_U(stream, pic_type, 3, "pic_type");

  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

static void encoder_state_write_bitstream_PTL(bitstream_t *stream,
                                              encoder_state_t * const state)
{
  // PTL
  // Profile Tier
  WRITE_U(stream, 0, 2, "general_profile_space");
  WRITE_U(stream, 0, 1, "general_tier_flag");
  // Main Profile == 1,  Main 10 profile == 2
  WRITE_U(stream, (state->encoder_control->bitdepth == 8)?1:2, 5, "general_profile_idc");
  /* Compatibility flags should be set at general_profile_idc
   *  (so with general_profile_idc = 1, compatibility_flag[1] should be 1)
   * According to specification, when compatibility_flag[1] is set,
   *  compatibility_flag[2] should be set too.
   */
  WRITE_U(stream, 3<<29, 32, "general_profile_compatibility_flag[]");

  WRITE_U(stream, 1, 1, "general_progressive_source_flag");
  WRITE_U(stream, state->encoder_control->in.source_scan_type!= 0, 1, "general_interlaced_source_flag");
  WRITE_U(stream, 0, 1, "general_non_packed_constraint_flag");
  WRITE_U(stream, 0, 1, "general_frame_only_constraint_flag");

  WRITE_U(stream, 0, 32, "XXX_reserved_zero_44bits[0..31]");
  WRITE_U(stream, 0, 12, "XXX_reserved_zero_44bits[32..43]");

  // end Profile Tier

  // Level 6.2 (general_level_idc is 30 * 6.2)
  WRITE_U(stream, 186, 8, "general_level_idc");

  WRITE_U(stream, 0, 1, "sub_layer_profile_present_flag");
  WRITE_U(stream, 0, 1, "sub_layer_level_present_flag");

  for (int i = 1; i < 8; i++) {
    WRITE_U(stream, 0, 2, "reserved_zero_2bits");
  }

  // end PTL
}

//*******************************************
//For scalability extension. TODO: merge with encoder_state_write_bitstream_PTL? Add asserts
//Handle case when profilePresentFlag is not set
static void encoder_state_write_bitstream_PTL_no_profile(bitstream_t *stream,
  encoder_state_t * const state)
{
  // PTL
  // Level 6.2 (general_level_idc is 30 * 6.2)
  WRITE_U(stream, 186, 8, "general_level_idc");

  WRITE_U(stream, 0, 1, "sub_layer_profile_present_flag");
  WRITE_U(stream, 0, 1, "sub_layer_level_present_flag");

  for (int i = 1; i < 8; i++) {
    WRITE_U(stream, 0, 2, "reserved_zero_2bits");
  }

  // end PTL
}

//TODO: Merger with other ptl / use profile function
static void encoder_state_write_bitstream_PTL_scalable(bitstream_t *stream,
  encoder_state_t * const state)
{
  // PTL
  // Profile Tier
  WRITE_U(stream, 0, 2, "general_profile_space");
  WRITE_U(stream, 0, 1, "general_tier_flag");

  // Main Profile == 1,  Main 10 profile == 2, Scalable Main (10) Profile == 7,
  // Use scalable profile here
  //uint8_t general_profile_idc = (state->encoder_control->bitdepth == 8) ? 7 : 8;

  WRITE_U(stream, 7, 5, "general_profile_idc"); //Should be 7 for 8 and 10 bit version

  /* Compatibility flags should be set at general_profile_idc
  *  (so with general_profile_idc = 1, compatibility_flag[1] should be 1)
  * According to specification, when compatibility_flag[1] is set,
  *  compatibility_flag[2] should be set too.
  */
  WRITE_U(stream, 3 << 29, 32, "general_profile_compatibility_flag[]");

  WRITE_U(stream, 1, 1, "general_progressive_source_flag");
  WRITE_U(stream, state->encoder_control->in.source_scan_type != 0, 1, "general_interlaced_source_flag");
  WRITE_U(stream, 0, 1, "general_non_packed_constraint_flag");
  WRITE_U(stream, 0, 1, "general_frame_only_constraint_flag");

  //Write settings specifyed in the specification
  WRITE_U(stream, 1, 1, "general_max_12bit_constraint_flag");
  WRITE_U(stream, 1, 1, "general_max_10bit_constraint_flag");
  WRITE_U(stream, (state->encoder_control->bitdepth == 8) ? 1 : 0, 1, "general_max_8bit_constraint_flag");
  WRITE_U(stream, 1, 1, "general_max_422chroma_constraint_flag");
  WRITE_U(stream, 1, 1, "general_max_420chroma_constraint_flag");
  WRITE_U(stream, 0, 1, "general_max_monochrome_constraint_flag");
  WRITE_U(stream, 0, 1, "general_intra_constraint_flag");
  WRITE_U(stream, 0, 1, "general_one_picture_only_constraint_flag");
  WRITE_U(stream, 1, 1, "general_lower_bit_rate_constraint_flag");

  WRITE_U(stream, 0, 32, "XXX_reserved_zero_34bits[0..31]");
  WRITE_U(stream, 0, 2, "XXX_reserved_zero_34bits[32..34]");

  WRITE_U(stream, 0, 1, "inbld_flag");
  // end Profile Tier

  // Level 6.2 (general_level_idc is 30 * 6.2)
  WRITE_U(stream, 186, 8, "general_level_idc");

  WRITE_U(stream, 0, 1, "sub_layer_profile_present_flag");
  WRITE_U(stream, 0, 1, "sub_layer_level_present_flag");

  for (int i = 1; i < 8; i++) {
    WRITE_U(stream, 0, 2, "reserved_zero_2bits");
  }

  // end PTL
}

//Write exstension data
static void encoder_state_write_bitsream_vps_extension(bitstream_t* stream,
                                                       encoder_state_t * const state)
{
  encoder_state_write_bitstream_PTL_no_profile(stream, state);
  
  uint8_t splitting_flag = 0; //TODO: implement splitting_flag in configuration?
  WRITE_U(stream, splitting_flag, 1, "splitting_flag");
  
  uint16_t  num_scalability_types = 1; //TODO: calculate from actual mask?

  //for (int i = 0; i < 16; i++){
    //TODO: Implement scalability mask flags in configuration?
    //WRITE_U(stream, 0, 1, "scalability_mask_flag[i]");
    //num_scalability_types += 0;
  //}
  //Write in one operation
  //TODO: implement scalability mask proberly
  WRITE_U(stream, 0x2000, 16, "scalability_mask_flag[i]"); // 0x2000 = 0010000000000000

  //Value from SHM
  uint8_t dimension_id_len[1] = { 1 };
  for (int j = 0; j < (num_scalability_types - splitting_flag); j++){
    WRITE_U(stream, dimension_id_len[j]-1, 3, "dimension_id_len_minus1[j]");
  }

  uint8_t vps_nuh_layer_id_present_flag = 0;
  WRITE_U(stream, vps_nuh_layer_id_present_flag, 1, "vps_nuh_layer_id_present_flag");

  uint16_t dimension_id[2][1] = { 0, 1 }; //Values from SHM
  
  //TODO: implement settings?
  for (int i = 1; i < state->layer->max_layers; i++) {
    if (vps_nuh_layer_id_present_flag){
      WRITE_U(stream, i, 6, "layer_id_in_nuh[i]"); //Sets alternate ids for layers. Defaults to i
    }
    if (!splitting_flag){
      for (int j = 0; j < num_scalability_types; j++){
        WRITE_U(stream, dimension_id[i][j], dimension_id_len[j], "dimension_id[i][j]");
      }
    }
  }

  uint8_t view_id_len = 0;
  WRITE_U(stream, view_id_len, 4, "view_id_len");
  //Not used. TODO: implement?
  //if (view_id_len > 0){
  //  //write "view_id_val[i]"
  //}

  //Generate a default layer reference structure where a layer only depends on the layer directly "below".
  //TODO: implement dependecy specification in config?
  for (int layer = 1; layer < state->layer->max_layers; layer++){
    for (int ref_layer = 0; ref_layer < layer; ref_layer++){
      uint8_t dep_flag = 0;
      if (ref_layer == layer - 1){
        dep_flag = 1;
      }
      WRITE_U(stream, dep_flag, 1, "direct_dependency_flag[i][j]");
    }
  }
  //TODO: implement num independent layers?
  //if numIndependentLayers > 1 write "num_add_layer_sets"

  //TODO: Not really necessary. Can be infered from maxTLayers
  uint8_t vps_sub_layers_max_minus1_present_flag = 1;
  WRITE_U(stream, vps_sub_layers_max_minus1_present_flag, 1, "vps_sub_layers_max_minus1_present_flag");
  //if vps_sub_layers_max_minus1_present_flag write sub_layers_max_minus1[i]=0
  if (vps_sub_layers_max_minus1_present_flag) {
    for (int i = 0; i < state->layer->max_layers; i++) {
      WRITE_U(stream, 0, 3, "sub_layers_vps_max_minus1[i]");
    }
  }

  //Control max number of temporal references
  WRITE_U(stream, 0, 1, "max_tid_ref_present_flag");
  //if max_tid_ref_present_flag write "max_tid_il_ref_pics_plus1[i][j]"

  WRITE_U(stream, 0, 1, "default_ref_layers_active_flag");

  //TODO: Add ptl for scalable main etc.
  uint8_t vps_num_profile_tier_level_minus1 = 2;
  WRITE_UE(stream, vps_num_profile_tier_level_minus1, "vps_num_profile_tier_level_minus1");
  
  // int i = vps_base_layer_internal_flag ? 2 : 1
  for (int i = 2; i <= vps_num_profile_tier_level_minus1; i++) {
    WRITE_U(stream, 1, 1, "vps_profile_present_flag[i]");
    encoder_state_write_bitstream_PTL_scalable(stream, state);
  }
  
  //TODO: Find out proper values? Set in config
  if (state->layer->num_layer_sets > 1) {
    WRITE_UE(stream, state->layer->num_output_layer_sets - state->layer->num_layer_sets, "num_add_olss");
    WRITE_U(stream, 1, 2, "default_output_layer_idc"); //value 1 says the layer with the largest nuh_layer_id is the output layer
  }
  
  //TODO: Move to a better place
  //ptl_idx for layer sets
  uint16_t ptl_idx[2][2] = { 0, 0, 1, 2 };

  //TODO: Add proper conditions
  uint8_t num_layers_in_id_list = 2;
  for (int i = 1; i < state->layer->num_output_layer_sets; i++) {
    //If numLayerSets > 2 && i >= numLayerSets write "layer_set_idx_for_ols_minus1[i]"
    //If i > vps_num_layer_sets_minus1 || defaultOutputLayerIdc == 2 write "output_layer_flag[i][j]
    
    //For j=0;j<NumLayersInIDList[OlsIdxToLsIdx[i]];j++ 
    //OlsIdxToLsIdx[1] == 1; NumLayersInIDList[1] == 2
    for (int j = 0; j < num_layers_in_id_list; j++) {
      //If NecessaryLayerFlag[i][j] && vps_num_prifile_tier_level_minus1 > 0
      //NecessaryLayerFlag[1][0] == true; NecessaryLayerFlag[1][1] == true
      if (vps_num_profile_tier_level_minus1 > 0) {
        WRITE_U(stream, ptl_idx[i][j], kvz_math_ceil_log2(vps_num_profile_tier_level_minus1+1), "profile_tier_level_idx[i][j]");
      }
    }

    //If numOutputLayersInOutputLayerSet[i] == 1 && NumDirectRefLayers[OlsHigestOutputLyerId[i]] > 0
    // numOutputLayersInOutputLayerSet[1] == 1; OlsHighestOutputLayerId[1] == 1; NumDirectRefLayers[1] == 1;
    WRITE_U(stream, 0, 1, "alt_output_layer_flag[i]");

  }


  //Write rep formats here
  //TODO: Implement own settings for alternate rep formats.
  //TODO: Need to add separate rep format when the layers have different sizes
  uint8_t vps_num_rep_formats_minus1 = state->layer->max_layers-1; //TODO: add rep format only for different sizes?
  WRITE_UE(stream, vps_num_rep_formats_minus1, "vps_num_rep_formats_minus1");
  
  {
    const encoder_control_t* encoder = state->encoder_control;
    for (int i = 0; i <= vps_num_rep_formats_minus1; i++) {
    //rep_format(){  
      //unsigned shift = i==0?0:1; //TODO: remove. Only for quick testing.
      //unsigned offset = i==0?0:8; //TODO: remove. Only for quick testing.
      WRITE_U(stream, encoder->in.width, 16, "pic_width_vps_in_luma_samples");
      WRITE_U(stream, encoder->in.height, 16, "pic_height_vps_in_luma_samples");

      uint8_t chroma_and_bit_depth_vps_present_flag = i==0; //Has to be one in the first rep format
      WRITE_U(stream, chroma_and_bit_depth_vps_present_flag, 1, "chroma_and_bit_depth_vps_present_flag");

      if (chroma_and_bit_depth_vps_present_flag) {
        WRITE_U(stream, encoder->chroma_format, 2, "chroma_format_vps_idc");
        if (encoder->chroma_format == KVZ_CSP_444) {
          WRITE_U(stream, 0, 1, "separate_colour_plane_flag");
        }
        WRITE_U(stream, encoder->bitdepth - 8, 4, "bit_depth_luma_minus8");
        WRITE_U(stream, encoder->bitdepth - 8, 4, "bit_depth_chroma_minus8");
      }
      uint8_t conformance_window_flag = encoder->in.width != encoder->in.real_width || encoder->in.height != encoder->in.real_height;
      WRITE_U(stream, conformance_window_flag, 1, "conformance_window_pvs_flag");

      if (conformance_window_flag) {
        WRITE_UE(stream, 0, "conf_win_vps_left_offset");
        WRITE_UE(stream, (encoder->in.width - encoder->in.real_width) >> 1, "conf_win_vps_right_offset");
        WRITE_UE(stream, 0, "conf_win_vps_top_offset");
        WRITE_UE(stream, (encoder->in.height - encoder->in.real_height) >> 1, "conf_win_vps_bottom_offset");
      }
  //}
      encoder = encoder->next_enc_ctrl;
    }
  }

  //if vps_num_rep_formats_minus1 > 0 Write rep_format_idx_present_flag here
  if( vps_num_rep_formats_minus1 > 0) {
    WRITE_U(stream, 0, 1, "rep_format_idx_present_flag");
  }

  //if rep_format_idx_present_flag Write rep_format_idx[i] here

  //Only allow one interlayer ref. TODO: allow multible?
  WRITE_U(stream, 1, 1, "max_one_active_ref_layer_flag");

  WRITE_U(stream, 0, 1, "vps_poc_lsb_aligned_flag");
  
  //TODO: implement numDirectRefLayer and layer_id_in_nuh
  //for (int i = 1; i < state->layer->max_layers; i++){
  //  if (NumDirectRefLayers[layer_id_in_nuh[i]] == 0){
  //    Starts from the first EL
  //    WRITE_U(stream, 0, 1, "poc_lsb_not_present_flag[0]"); 
  //  }
  //}


  //dpb_size(){
  //TODO: Implement properly.
  uint16_t max_dec_pic_buffering_minus1 = state->encoder_control->cfg->ref_frames + state->encoder_control->cfg->gop_len;
  uint16_t max_vps_dec_pic_buffering_minus1[2][2][1] = { 0, 0, max_dec_pic_buffering_minus1,
                                                         max_dec_pic_buffering_minus1 }; //needs to be in line (<=) sps values
  //for (int i = 1; i < state->layer->num_output_layer_sets; i++) {
  //  state->layer->num_output_layer_sets == 2
  WRITE_U(stream, 0, 1, "sub_layer_flag_info_present_flag[i]");
  //  for (int j=0; j <= MaxSubLayersInLayerSetMinus1[OlsIdxToLsIdx[i]]; j++){
  //    OlsIdxToLsIdx[1]==1; MaxSubLayerInLayersSetMinus1[1] == 0
  //    if( j>0 && sub_layer_flag_info_present_flag[i] ) write "sub_layer_dpb_info_present_flag[i][j]
  //    //Always signal for j==0 => sub_layer_dpb_info_present_flag[1][0] == true
  //    if(sub_layer_dpb_info_present_flag[i][j]){
  for(int k=0; k < num_layers_in_id_list; k++){
  //        if (NecessaryLayerFlag[i][k] && (vps_base_layer_internal_flag||(LayerSetLayerIdList[OlsIdxToLsIdx[i]][k] != 0 ) ) ) {
  //        NecessaryLayerFlag[i][k] == true; vps_base_layer_internal_flag == true
    WRITE_UE(stream, max_vps_dec_pic_buffering_minus1[1][k][0], "max_vps_dec_pic_buffering_minus1[i][k][j]"); //Max pics in DPF
  //        }
  }
  WRITE_UE(stream, 3, "max_vps_num_reorder_pics[i][j]"); //value from SHM
  WRITE_UE(stream, 0, "max_vps_latency_increase_plus1[i][j]"); //value from SHM
  //    }
  //  }
  //}
  //}

  uint8_t direct_dep_type_len_minus2 = 0;
  WRITE_UE(stream, direct_dep_type_len_minus2, "direct_dep_type_len_minus2");
  WRITE_U(stream, 1, 1, "direct_dependency_all_layers_flag");
  //if "direct_dependency_all_layers_flag"
  WRITE_U(stream, 2, direct_dep_type_len_minus2+2, "direct_dependency_all_layers_type"); //Value 2 used in SHM. 0: Only use IL sample prediction, 1: Only use IL Motion prediction, 2: Use both
  //Else write separately for each layer

  WRITE_UE(stream, 0, "vps_non_vui_extension_length");
  //Write non vui extension data here

  WRITE_U(stream, 0, 1, "vps_vui_present_flag");
  //write vui stuff here
}
//*********************************************************

static void encoder_state_write_bitstream_vid_parameter_set(bitstream_t* stream,
                                                            encoder_state_t * const state)
{
#ifdef KVZ_DEBUG
  printf("=========== Video Parameter Set ID: 0 ===========\n");
#endif

  WRITE_U(stream, 0, 4, "vps_video_parameter_set_id");
  WRITE_U(stream, 3, 2, "vps_reserved_three_2bits" ); //Vps_base_layer_internal_flag and vps_base_layer_available_flag

  //*********************************************
  //For scalable extension. TODO: Move somewhere else?
   WRITE_U(stream, state->layer->max_layers-1, 6, "vps_max_layers_minus1" );
  //*********************************************

  WRITE_U(stream, 1, 3, "vps_max_sub_layers_minus1");
  WRITE_U(stream, 0, 1, "vps_temporal_id_nesting_flag");
  WRITE_U(stream, 0xffff, 16, "vps_reserved_ffff_16bits");

  encoder_state_write_bitstream_PTL(stream, state);

  WRITE_U(stream, 0, 1, "vps_sub_layer_ordering_info_present_flag");

  //for each layer
  for (int i = 0; i < 1; i++) {
    //*********************************************
    //For scalable extension. TODO: Why was it previously 1?
    WRITE_UE(stream, state->encoder_control->cfg->ref_frames
              + state->encoder_control->cfg->gop_len, "vps_max_dec_pic_buffering_minus1[i]");
    //*********************************************
    WRITE_UE(stream, 0, "vps_num_reorder_pics");
    WRITE_UE(stream, 0, "vps_max_latency_increase");
  }

  //*********************************************
  //For scalable extension. TODO: Move somewhere else?
  WRITE_U(stream, state->layer->max_layers - 1, 6, "vps_max_layer_id");
  
  //TODO: Find out what it does
  WRITE_UE(stream, state->layer->num_layer_sets - 1, "vps_num_layer_sets_minus1");
  for (int i = 1; i < state->layer->num_layer_sets; i++) {
    for (int j = 0; j < state->layer->max_layers; j++) {
      WRITE_U(stream, 1, 1, "layer_id_included_flag[i][j]"); //By default include all.
    }
  }
  //*********************************************

  WRITE_U(stream, 0, 1, "vps_timing_info_present_flag");

  //IF timing info
  //END IF

  //WRITE_U(stream, 0, 1, "vps_extension_flag");
  //*********************************************
  //For scalable extension. TODO: Move somewhere else? Set based on extension used
  WRITE_U(stream, (state->layer->max_layers - 1) > 0 ? 1 : 0 , 1, "vps_extension_flag");

  //Align with ones
  //TODO: a better way?
  if (state->layer->max_layers > 0){
    //while (stream->cur_bit != 0) {
    if (stream->cur_bit != 0){
      //kvz_bitstream_align(stream);
      // while(!aligned) "vbs_extension_alignment_bit_equal_to_one"
      //WRITE_U(stream, 1, 1, "vps_extension_alignment_bit_equal_to_one");
      WRITE_U(stream, 0xffff, 8-(stream->cur_bit & 7), "vps_extension_alignment_bit_equal_to_one");
    }

    //Write vps_extension()
    encoder_state_write_bitsream_vps_extension(stream, state);

    WRITE_U(stream, 0, 1, "vps_extension2_flag");
  }
  //*********************************************

  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

static void encoder_state_write_bitstream_scaling_list(bitstream_t *stream,
                                                       encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;
  uint32_t size_id;
  for (size_id = 0; size_id < SCALING_LIST_SIZE_NUM; size_id++) {
    int32_t list_id;
    for (list_id = 0; list_id < kvz_g_scaling_list_num[size_id]; list_id++) {
      uint8_t scaling_list_pred_mode_flag = 1;
      int32_t pred_list_idx;
      int32_t i;
      uint32_t ref_matrix_id = UINT32_MAX;

      for (pred_list_idx = list_id; pred_list_idx >= 0; pred_list_idx--) {
        const int32_t * const pred_list  = (list_id == pred_list_idx) ?
                                     kvz_scalinglist_get_default(size_id, pred_list_idx) :
                                     encoder->scaling_list.scaling_list_coeff[size_id][pred_list_idx];

        if (!memcmp(encoder->scaling_list.scaling_list_coeff[size_id][list_id], pred_list, sizeof(int32_t) * MIN(8, kvz_g_scaling_list_size[size_id])) &&
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
        const int32_t coef_num = MIN(MAX_MATRIX_COEF_NUM, kvz_g_scaling_list_size[size_id]);
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


static void encoder_state_write_bitstream_VUI(bitstream_t *stream,
                                              encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;
#ifdef KVZ_DEBUG
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
    WRITE_U(stream, encoder->vui.videoformat, 3, "chroma_format");
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
  WRITE_U(stream, encoder->vui.field_seq_flag, 1, "field_seq_flag"); // 0: frames, 1: fields
  WRITE_U(stream, encoder->vui.frame_field_info_present_flag, 1, "frame_field_info_present_flag");
  WRITE_U(stream, 0, 1, "default_display_window_flag");

  //IF default display window
  //ENDIF

  WRITE_U(stream, encoder->vui.timing_info_present_flag, 1, "vui_timing_info_present_flag");
  if (encoder->vui.timing_info_present_flag) {
    WRITE_U(stream, encoder->vui.num_units_in_tick, 32, "vui_num_units_in_tick");
    WRITE_U(stream, encoder->vui.time_scale, 32, "vui_time_scale");

    WRITE_U(stream, 0, 1, "vui_poc_proportional_to_timing_flag");
    WRITE_U(stream, 0, 1, "vui_hrd_parameters_present_flag");    
  }
  
  WRITE_U(stream, 0, 1, "bitstream_restriction_flag");

  //IF bitstream restriction
  //ENDIF
}


static void encoder_state_write_bitstream_SPS_extension(bitstream_t *stream,
                                                        encoder_state_t * const state)
{
  if (state->encoder_control->cfg->implicit_rdpcm &&
      state->encoder_control->cfg->lossless) {
    WRITE_U(stream, 1, 1, "sps_extension_present_flag");

    WRITE_U(stream, 1, 1, "sps_range_extension_flag");
    WRITE_U(stream, 0, 1, "sps_multilayer_extension_flag");
    WRITE_U(stream, 0, 1, "sps_3d_extension_flag");
    WRITE_U(stream, 0, 5, "sps_extension_5bits");

    WRITE_U(stream, 0, 1, "transform_skip_rotation_enabled_flag");
    WRITE_U(stream, 0, 1, "transform_skip_context_enabled_flag");
    WRITE_U(stream, 1, 1, "implicit_rdpcm_enabled_flag");
    WRITE_U(stream, 0, 1, "explicit_rdpcm_enabled_flag");
    WRITE_U(stream, 0, 1, "extended_precision_processing_flag");
    WRITE_U(stream, 0, 1, "intra_smoothing_disabled_flag");
    WRITE_U(stream, 0, 1, "high_precision_offsets_enabled_flag");
    WRITE_U(stream, 0, 1, "persistent_rice_adaptation_enabled_flag");
    WRITE_U(stream, 0, 1, "cabac_bypass_alignment_enabled_flag");
  } else {
    WRITE_U(stream, 0, 1, "sps_extension_present_flag");
  }
}

//*******************************************
//For scalability extension. TODO: remove if pointless
void writeSTermRSet_(encoder_state_t* const state, int neg_ref_minus)
{
  
  const encoder_control_t* const encoder = state->encoder_control;
  bitstream_t* const stream = &state->stream;

  //IF stRpsIdx != 0
  if( neg_ref_minus != 0) {
    WRITE_U(stream, 0, 1, "inter_ref_pic_set_prediction_flag"); //TODO: Use inter rps pred?
  }

  int j;
  int ref_negative = state->encoder_control->cfg->ref_frames - neg_ref_minus;
  ref_negative = ref_negative < 0 ? 0 : ref_negative;
  int ref_positive = 0;
  
  //TODO: Make a better implementation. Use neg_ref_minus?
  if( ref_negative > 0 && state->layer->layer_id > 0) --ref_negative; //One frame needs to be for the ILR ref

  int last_poc = 0;
  int poc_shift = 0;

  WRITE_UE(stream, ref_negative, "num_negative_pics");
  WRITE_UE(stream, ref_positive, "num_positive_pics");
  for (j = 0; j < ref_negative; j++) {
    int8_t delta_poc = 0;

    if (encoder->cfg->gop_len) {
      int8_t found = 0;
      do {
        delta_poc = encoder->cfg->gop[state->frame->gop_offset].ref_neg[j + poc_shift];
        for (int i = 0; i < state->frame->ref->used_size; i++) {
          if (state->frame->ref->pocs[i] == state->frame->poc - delta_poc) {
            found = 1;
            break;
          }
        }
        if (!found) poc_shift++;
        if (j + poc_shift == ref_negative) {
          fprintf(stderr, "Failure, reference not found!");
          exit(EXIT_FAILURE);
        }
      } while (!found);
    }

    WRITE_UE(stream, encoder->cfg->gop_len?delta_poc - last_poc - 1:0, "delta_poc_s0_minus1");
    last_poc = delta_poc;
    WRITE_U(stream,1,1, "used_by_curr_pic_s0_flag");
  }
  last_poc = 0;
  poc_shift = 0;
  for (j = 0; j < ref_positive; j++) {
    int8_t delta_poc = 0;

    if (encoder->cfg->gop_len) {
      int8_t found = 0;
      do {
        delta_poc = encoder->cfg->gop[state->frame->gop_offset].ref_pos[j + poc_shift];
        for (int i = 0; i < state->frame->ref->used_size; i++) {
          if (state->frame->ref->pocs[i] == state->frame->poc + delta_poc) {
            found = 1;
            break;
          }
        }
        if (!found) poc_shift++;
        if (j + poc_shift == ref_positive) {
          fprintf(stderr, "Failure, reference not found!");
          exit(EXIT_FAILURE);
        }
      } while (!found);
    }

    WRITE_UE(stream, encoder->cfg->gop_len ? delta_poc - last_poc - 1 : 0, "delta_poc_s1_minus1");
    last_poc = delta_poc;
    WRITE_U(stream, 1, 1, "used_by_curr_pic_s1_flag");
  }
}
//Use inter rps_set pred. Write ref sets for the case that n prev frames are referenced.
void writeSTermRSet(encoder_state_t* const state, unsigned idx )
{
  
  const encoder_control_t* const encoder = state->encoder_control;
  bitstream_t* const stream = &state->stream;

  //IF stRpsIdx != 0
  if( idx != 0) {
    WRITE_U(stream, 1, 1, "inter_ref_pic_set_prediction_flag"); //TODO: Use inter rps pred?
  }

  int j;
  int ref_negative = idx + 1;//state->encoder_control->cfg->ref_frames;
  int ref_positive = 0;
  
  //TODO: Make a better implementation. Use neg_ref_minus?
  if(state->layer->layer_id > 0 && encoder->cfg->ref_frames <= ref_negative+ref_positive) --ref_negative; //One frame needs to be for the ILR ref

  if( idx != 0) {
    //IF slice header WRITE "delta_idx_minus1" ?
    WRITE_U( stream, 1, 1, "delta_rps_sign" ); //only ref prev frame rps
    WRITE_UE( stream, 0, "abs_delta_rps_minus1" ); //only ref prev frames rps
    for( int i = 0; i < ref_negative + ref_positive; i++ ) {
      WRITE_U( stream, 1, 1, "used_by_curr_pic_flag[j]" );
      //IF !used_by_curr_pic_flag[j] WRITE use_delta_flag[j]
    }
  }
  else {
    int last_poc = 0;
    int poc_shift = 0;

    WRITE_UE(stream, ref_negative, "num_negative_pics");
    WRITE_UE(stream, ref_positive, "num_positive_pics");
    for (j = 0; j < ref_negative; j++) {
      int8_t delta_poc = 0;

      if (encoder->cfg->gop_len) {
        int8_t found = 0;
        do {
          delta_poc = encoder->cfg->gop[state->frame->gop_offset].ref_neg[j + poc_shift];
          for (int i = 0; i < state->frame->ref->used_size; i++) {
            if (state->frame->ref->pocs[i] == state->frame->poc - delta_poc) {
              found = 1;
              break;
            }
          }
          if (!found) poc_shift++;
          if (j + poc_shift == ref_negative) {
            fprintf(stderr, "Failure, reference not found!");
            exit(EXIT_FAILURE);
          }
        } while (!found);
      }

      WRITE_UE(stream, encoder->cfg->gop_len ? delta_poc - last_poc - 1 : 0, "delta_poc_s0_minus1");
      last_poc = delta_poc;
      WRITE_U(stream, 1, 1, "used_by_curr_pic_s0_flag");
    }
    last_poc = 0;
    poc_shift = 0;
    for (j = 0; j < ref_positive; j++) {
      int8_t delta_poc = 0;

      if (encoder->cfg->gop_len) {
        int8_t found = 0;
        do {
          delta_poc = encoder->cfg->gop[state->frame->gop_offset].ref_pos[j + poc_shift];
          for (int i = 0; i < state->frame->ref->used_size; i++) {
            if (state->frame->ref->pocs[i] == state->frame->poc + delta_poc) {
              found = 1;
              break;
            }
          }
          if (!found) poc_shift++;
          if (j + poc_shift == ref_positive) {
            fprintf(stderr, "Failure, reference not found!");
            exit(EXIT_FAILURE);
          }
        } while (!found);
      }

      WRITE_UE(stream, encoder->cfg->gop_len ? delta_poc - last_poc - 1 : 0, "delta_poc_s1_minus1");
      last_poc = delta_poc;
      WRITE_U(stream, 1, 1, "used_by_curr_pic_s1_flag");
    }
  }
}
//*******************************************
  
static void encoder_state_write_bitstream_seq_parameter_set(bitstream_t* stream,
  encoder_state_t * const state)
{
  const encoder_control_t * encoder = state->encoder_control;

#ifdef KVZ_DEBUG
  printf("=========== Sequence Parameter Set ID: 0 ===========\n");
#endif

  // TODO: profile IDC and level IDC should be defined later on
  WRITE_U(stream, 0, 4, "sps_video_parameter_set_id");

  //*******************************************
  //For scalability extension. TODO: add asserts.
  //TODO: set sps_max_sub_layers in cfg?
  if (state->layer->layer_id == 0) {
    WRITE_U(stream, 1, 3, "sps_max_sub_layers_minus1");
  }
  else {
    WRITE_U(stream, 7, 3, "sps_ext_or_max_sub_layers_minus1")
  }
  //TODO: Add sps_ext_of_max_sub_layers_minus1 to cfg?
  uint8_t multi_layer_ext_sps_flag = (state->layer->layer_id != 0); // && sps_ext_or_max_sub_layers_minus1 == 7);

  if (!multi_layer_ext_sps_flag) {
    WRITE_U(stream, 0, 1, "sps_temporal_id_nesting_flag");
    encoder_state_write_bitstream_PTL(stream, state);
  }

  // ***********************************************
  // Modified for SHVC
  // parameter set ids are global (layer id does not matter)
  // so set the correct sps id (use layer id as the pset id)
  // vps is the same for all layers
  WRITE_UE(stream, state->layer->layer_id, "sps_seq_parameter_set_id");
  // ***********************************************
 

  if (multi_layer_ext_sps_flag) {
    WRITE_U(stream, 0, 1, "update_rep_format_flag");
  }
  else {
    WRITE_UE(stream, encoder->chroma_format, "chroma_format_idc");
  
    if (encoder->chroma_format == KVZ_CSP_444) {
      WRITE_U(stream, 0, 1, "separate_colour_plane_flag");
    }



    WRITE_UE(stream, encoder->in.width, "pic_width_in_luma_samples");
    WRITE_UE(stream, encoder->in.height, "pic_height_in_luma_samples");

    if (encoder->in.width != encoder->in.real_width || encoder->in.height != encoder->in.real_height) {
      // The standard does not seem to allow setting conf_win values such that
      // the number of luma samples is not a multiple of 2. Options are to either
      // hide one line or show an extra line of non-video. Neither seems like a
      // very good option, so let's not even try.
      assert(!(encoder->in.width % 2));
      WRITE_U(stream, 1, 1, "conformance_window_flag");
      WRITE_UE(stream, 0, "conf_win_left_offset");
      WRITE_UE(stream, (encoder->in.width - encoder->in.real_width) >> 1,
        "conf_win_right_offset");
      WRITE_UE(stream, 0, "conf_win_top_offset");
      WRITE_UE(stream, (encoder->in.height - encoder->in.real_height) >> 1,
        "conf_win_bottom_offset");
    }
    else {
      WRITE_U(stream, 0, 1, "conformance_window_flag");
    }

    //IF window flag
    //END IF

    WRITE_UE(stream, encoder->bitdepth - 8, "bit_depth_luma_minus8");
    WRITE_UE(stream, encoder->bitdepth - 8, "bit_depth_chroma_minus8");

  }
  
  
  WRITE_UE(stream, 1, "log2_max_pic_order_cnt_lsb_minus4");
  if (!multi_layer_ext_sps_flag) {
    WRITE_U(stream, 0, 1, "sps_sub_layer_ordering_info_present_flag");

    //for each layer
    if (encoder->cfg->gop_lowdelay) {
      WRITE_UE(stream, encoder->cfg->ref_frames, "sps_max_dec_pic_buffering");
      WRITE_UE(stream, 0, "sps_num_reorder_pics");
    }
    else {
      WRITE_UE(stream, encoder->cfg->ref_frames + encoder->cfg->gop_len, "sps_max_dec_pic_buffering");
      WRITE_UE(stream, encoder->cfg->gop_len, "sps_num_reorder_pics");
    }
    WRITE_UE(stream, 0, "sps_max_latency_increase");
    //end for
  }

  WRITE_UE(stream, MIN_SIZE-3, "log2_min_coding_block_size_minus3");
  WRITE_UE(stream, MAX_DEPTH, "log2_diff_max_min_coding_block_size");
  WRITE_UE(stream, 0, "log2_min_transform_block_size_minus2");   // 4x4
  WRITE_UE(stream, 3, "log2_diff_max_min_transform_block_size"); // 4x4...32x32
  WRITE_UE(stream, TR_DEPTH_INTER, "max_transform_hierarchy_depth_inter");
  WRITE_UE(stream, encoder->tr_depth_intra, "max_transform_hierarchy_depth_intra");

  // scaling list
  WRITE_U(stream, encoder->scaling_list.enable, 1, "scaling_list_enable_flag");
  if (encoder->scaling_list.enable) {
    //TODO: Infer scaling list from previous layer?
    uint8_t sps_infer_scaling_list_flag = 0;

    if (multi_layer_ext_sps_flag) {
      WRITE_U(stream, sps_infer_scaling_list_flag, 1, "sps_infer_scaling_list_flag");
    }
    if (sps_infer_scaling_list_flag) {
      WRITE_U(stream, state->layer->layer_id - 1, 6, "sps_scaling_list_ref_layer_id");
    }
    else {
      WRITE_U(stream, 1, 1, "sps_scaling_list_data_present_flag");
      encoder_state_write_bitstream_scaling_list(stream, state);
    }
  }

  WRITE_U(stream, (encoder->cfg->amp_enable ? 1 : 0), 1, "amp_enabled_flag");

  WRITE_U(stream, encoder->sao_enable ? 1 : 0, 1,
          "sample_adaptive_offset_enabled_flag");
  WRITE_U(stream, ENABLE_PCM, 1, "pcm_enabled_flag");
  #if ENABLE_PCM == 1
    WRITE_U(stream, 7, 4, "pcm_sample_bit_depth_luma_minus1");
    WRITE_U(stream, 7, 4, "pcm_sample_bit_depth_chroma_minus1");
    WRITE_UE(stream, 0, "log2_min_pcm_coding_block_size_minus3");
    WRITE_UE(stream, 2, "log2_diff_max_min_pcm_coding_block_size");
    WRITE_U(stream, 1, 1, "pcm_loop_filter_disable_flag");
  #endif

  if (state->layer->max_layers > 1) {
    //Need to make sure the first frames don't reference non-existant pocs
    uint8_t num_short_term_ref_pic_sets = state->encoder_control->cfg->ref_frames; //TODO: a beter implementation?
    if (state->layer->layer_id > 0) num_short_term_ref_pic_sets = num_short_term_ref_pic_sets > 1 ? num_short_term_ref_pic_sets - 1 : 1; //Reserve one "reference" for ILR 
    WRITE_UE(stream, num_short_term_ref_pic_sets, "num_short_term_ref_pic_sets");
         //IF num short term ref pic sets
    for (int i = 0; i < num_short_term_ref_pic_sets; i++) {
      writeSTermRSet(state, i);
    }
  }
  else {
    WRITE_UE(stream, 0, "num_short_term_ref_pic_sets");
  }
  
  //ENDIF


//*******************************************  
  WRITE_U(stream, 0, 1, "long_term_ref_pics_present_flag");

  //IF long_term_ref_pics_present
  //ENDIF

  WRITE_U(stream, state->encoder_control->cfg->tmvp_enable, 1,
          "sps_temporal_mvp_enable_flag");
  WRITE_U(stream, 0, 1, "sps_strong_intra_smoothing_enable_flag");
  WRITE_U(stream, 1, 1, "vui_parameters_present_flag");

  encoder_state_write_bitstream_VUI(stream, state);

  encoder_state_write_bitstream_SPS_extension(stream, state);

  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

static void encoder_state_write_bitstream_pic_parameter_set(bitstream_t* stream,
                                                            encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;
#ifdef KVZ_DEBUG
  printf("=========== Picture Parameter Set ID: 0 ===========\n");
#endif
  // ***********************************************
  // Modified for SHVC
  // parameter set ids are global (layer id does not matter)
  // so set the correct pps/sps id (use layer id as the pset id)
 
  WRITE_UE(stream, state->layer->layer_id, "pic_parameter_set_id");
  WRITE_UE(stream, state->layer->layer_id, "seq_parameter_set_id");
  // ***********************************************
 
  
  WRITE_U(stream, 0, 1, "dependent_slice_segments_enabled_flag");
  WRITE_U(stream, 0, 1, "output_flag_present_flag");
  WRITE_U(stream, 0, 3, "num_extra_slice_header_bits");
  WRITE_U(stream, encoder->sign_hiding, 1, "sign_data_hiding_flag");
  WRITE_U(stream, 0, 1, "cabac_init_present_flag");

  WRITE_UE(stream, 0, "num_ref_idx_l0_default_active_minus1");
  WRITE_UE(stream, 0, "num_ref_idx_l1_default_active_minus1");
  WRITE_SE(stream, ((int8_t)encoder->cfg->qp) - 26, "pic_init_qp_minus26");
  WRITE_U(stream, 0, 1, "constrained_intra_pred_flag");
  WRITE_U(stream, encoder->trskip_enable, 1, "transform_skip_enabled_flag");
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
  WRITE_U(stream, encoder->cfg->lossless, 1, "transquant_bypass_enable_flag");
  WRITE_U(stream, encoder->tiles_enable, 1, "tiles_enabled_flag");
  //wavefronts
  WRITE_U(stream, encoder->wpp, 1, "entropy_coding_sync_enabled_flag");

  if (encoder->tiles_enable) {
    WRITE_UE(stream, encoder->tiles_num_tile_columns - 1, "num_tile_columns_minus1");
    WRITE_UE(stream, encoder->tiles_num_tile_rows - 1, "num_tile_rows_minus1");
    
    WRITE_U(stream, encoder->tiles_uniform_spacing_flag, 1, "uniform_spacing_flag");
    
    if (!encoder->tiles_uniform_spacing_flag) {
      int i;
      for (i = 0; i < encoder->tiles_num_tile_columns - 1; ++i) {
        WRITE_UE(stream, encoder->tiles_col_width[i] - 1, "column_width_minus1[...]");
      }
      for (i = 0; i < encoder->tiles_num_tile_rows - 1; ++i) {
        WRITE_UE(stream, encoder->tiles_row_height[i] - 1, "row_height_minus1[...]");
      }
    }
    WRITE_U(stream, 0, 1, "loop_filter_across_tiles_enabled_flag");
    
  }
  
  WRITE_U(stream, 0, 1, "loop_filter_across_slice_flag");
  WRITE_U(stream, 1, 1, "deblocking_filter_control_present_flag");

  //IF deblocking_filter
    WRITE_U(stream, 0, 1, "deblocking_filter_override_enabled_flag");
  WRITE_U(stream, encoder->deblock_enable ? 0 : 1, 1,
          "pps_disable_deblocking_filter_flag");

    //IF !disabled
  if (encoder->deblock_enable) {
     WRITE_SE(stream, encoder->beta_offset_div2, "beta_offset_div2");
     WRITE_SE(stream, encoder->tc_offset_div2, "tc_offset_div2");
    }

    //ENDIF
  //ENDIF
  WRITE_U(stream, 0, 1, "pps_scaling_list_data_present_flag");
  //IF scaling_list
  //ENDIF

  //*******************************************
  //For scalability extension. TODO: add asserts.
  //TODO: Make a better check? Move to somewhere else?
  state->layer->list_modification_present_flag = (state->layer->layer_id>0)&&(state->encoder_control->cfg->ref_frames>1)&&(state->encoder_control->cfg->intra_period!=1)?1:0;
  WRITE_U(stream, state->layer->list_modification_present_flag, 1, "lists_modification_present_flag");
  WRITE_UE(stream, 0, "log2_parallel_merge_level_minus2");
    
  //TODO: set in cfg?
  uint8_t slice_segment_header_extension_present_flag = 0;
  uint8_t pps_extension_flag = state->layer->layer_id!=0; // state->layer->max_layers > 1 ? 1 : 0;
  uint8_t pps_multilayer_extension_flag = pps_extension_flag; //pps_extension_flag;
  WRITE_U(stream, slice_segment_header_extension_present_flag, 1, "slice_segment_header_extension_present_flag");
  WRITE_U(stream, pps_extension_flag, 1, "pps_extension_flag");

  //Write all possible extension flags here if pps_extension_flag defined
  if (pps_extension_flag) {
    WRITE_U(stream, 0, 1, "pps_range_extension_flag");
    WRITE_U(stream, pps_multilayer_extension_flag, 1, "pps_multilayer_extension_flag");
    WRITE_U(stream, 0, 1, "pps_3d_extension_flag");
    WRITE_U(stream, 0, 5, "pps_extension_5bits");
  }

  if (pps_multilayer_extension_flag) {
    //pps_multilayer_extension()
    WRITE_U(stream, 0, 1, "poc_reset_info_present_flag");

    uint8_t pps_infer_scaling_list_flag = 0; //TODO: Infer list?
    WRITE_U(stream, 0, 1, "pps_infer_scaling_list_flag");
    if (pps_infer_scaling_list_flag) {
      WRITE_U(stream, state->layer->layer_id - 1, 6, "pps_scaling_list_ref_layer_id");
    }

    //TODO: Check actual number of differently sized layers?
    //TODO: Account for more complex ref structure?
    uint8_t num_ref_loc_offsets = 1;//state->layer->max_layers-1; 
    WRITE_UE(stream, num_ref_loc_offsets, "num_ref_loc_offsets");
    for( int i = 0; i < num_ref_loc_offsets; i++ ) {
      WRITE_U(stream, state->layer->layer_id-1-i, 6, "ref_loc_offset_layer_id[i]");
      
      //TODO: For some reason offsets are shifted in the SHM. Why? 
      uint8_t scaled_ref_layer_offset_present_flag = 0; //TODO: get from somewhere
      WRITE_U(stream, scaled_ref_layer_offset_present_flag, 1, "scaled_ref_layer_offset_present_flag" );
      if( scaled_ref_layer_offset_present_flag ) {
        WRITE_SE(stream, 0>>1, "scaled_ref_layer_left_offset");
        WRITE_SE(stream, 0>>1, "scaled_ref_layer_top_offset");
        WRITE_SE(stream, 0>>1, "scaled_ref_layer_right_offset");
        WRITE_SE(stream, 0>>1, "scaled_ref_layer_bottom_offset");
      }

      uint8_t ref_region_layer_offset_present_flag = 1; //TODO: get from somewhere
      WRITE_U(stream, ref_region_layer_offset_present_flag, 1, "ref_region_layer_offset_present_flag" );
      //TODO: Calculate properly.
      if( ref_region_layer_offset_present_flag ) {
        WRITE_SE(stream, 0>>1, "ref_region_layer_left_offset");
        WRITE_SE(stream, 0>>1, "ref_region_layer_top_offset");
        WRITE_SE(stream, 0>>1, "ref_region_layer_right_offset");
        WRITE_SE(stream, 4>>1, "ref_region_layer_bottom_offset");
      }

      uint8_t resample_phase_set_presetn_flag = 0; //TODO: get from somewhere
      WRITE_U(stream, resample_phase_set_presetn_flag, 1, "resample_phase_set_present_flag");
      if( resample_phase_set_presetn_flag ) {
        WRITE_UE(stream, 0, "phase_hor_luma");
        WRITE_UE(stream, 0, "phase_ver_luma");
        WRITE_UE(stream, 8, "phase_hor_chroma_plus8");
        WRITE_UE(stream, 8, "phase_ver_chorma_plus8");
      }
    }

    WRITE_U(stream, 0, 1, "colour_mapping_enabled_flag");
    //if colour_mapping_enabled_flag write colormapping 
  }

  //*******************************************

  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

static void encoder_state_write_bitstream_prefix_sei_version(encoder_state_t * const state)
{
#define STR_BUF_LEN 1000
  bitstream_t * const stream = &state->stream;
  int i, length;
  char buf[STR_BUF_LEN] = { 0 };
  char *s = buf + 16;
  const kvz_config * const cfg = state->encoder_control->cfg;

  // random uuid_iso_iec_11578 generated with www.famkruithof.net/uuid/uuidgen
  static const uint8_t uuid[16] = {
    0x32, 0xfe, 0x46, 0x6c, 0x98, 0x41, 0x42, 0x69,
    0xae, 0x35, 0x6a, 0x91, 0x54, 0x9e, 0xf3, 0xf1
  };
  memcpy(buf, uuid, 16);

  // user_data_payload_byte
  s += sprintf(s, "Kvazaar HEVC Encoder v. " VERSION_STRING " - "
                  "Copyleft 2012-2015 - http://ultravideo.cs.tut.fi/ - options:");
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

  // The bitstream is already aligned, but align it anyway.
  kvz_bitstream_align(stream);

#undef STR_BUF_LEN
}

/*
static void encoder_state_write_active_parameter_sets_sei_message(encoder_state_t * const state) {

  const encoder_control_t * const encoder = state->encoder_control;
  bitstream_t * const stream = &state->stream;

  int i = 0;

  int active_vps_id = 0;
  int self_contained_cvs_flag = 0;
  int no_parameter_set_update_flag = 0;
  int num_sps_ids_minus1 = 0;
  int layer_sps_idx = 0;
  int active_seq_parameter_set_id = 0;
  int vps_base_layer_internal_flag = 0;

  int max_layers_minus1 = 0;

  WRITE_U(stream, 129, 8, "last_payload_type_byte"); //active_parameter_sets
  WRITE_U(stream, 2, 8, "last_payload_size_byte");
  WRITE_U(stream, active_vps_id, 4, "active_video_parameter_set_id");
  WRITE_U(stream, self_contained_cvs_flag, 1, "self_contained_cvs_flag");
  WRITE_U(stream, no_parameter_set_update_flag, 1, "no_parameter_set_update_flag");
  WRITE_UE(stream, num_sps_ids_minus1, "num_sps_ids_minus1");
  //for (i = 0; i <= num_sps_ids_minus1; ++i) {
  WRITE_UE(stream, active_seq_parameter_set_id, "active_seq_parameter_set_id");
  //}
  // for (i = vps_base_layer_internal_flag; i <= max_layers_minus1; ++i){
  WRITE_UE(stream, layer_sps_idx, "layer_sps_idx");
  //}

  kvz_bitstream_rbsp_trailing_bits(stream); //rbsp_trailing_bits
}
*/

static void encoder_state_write_picture_timing_sei_message(encoder_state_t * const state) {

  bitstream_t * const stream = &state->stream;

  if (state->encoder_control->vui.frame_field_info_present_flag){

    int8_t odd_picture = state->frame->num % 2;
    int8_t pic_struct = 0; //0: progressive picture, 1: top field, 2: bottom field, 3...
    int8_t source_scan_type = 1; //0: interlaced, 1: progressive

    switch (state->tile->frame->source->interlacing){
    case 0: //Progressive frame
      pic_struct = 0;
      source_scan_type = 1;
      break;
    case 1: //Top field first
      pic_struct = odd_picture ? 2 : 1;
      source_scan_type = 0;
      break;
    case 2: //Bottom field first
      pic_struct = odd_picture ? 1 : 2;
      source_scan_type = 0;
      break;
    default:
      assert(0); //Should never execute
      break;
    }

    WRITE_U(stream, 1, 8, "last_payload_type_byte"); //pic_timing
    WRITE_U(stream, 1, 8, "last_payload_size_byte");
    WRITE_U(stream, pic_struct, 4, "pic_struct");
    WRITE_U(stream, source_scan_type, 2, "source_scan_type");
    WRITE_U(stream, 0, 1, "duplicate_flag");

    kvz_bitstream_align(stream);
  }
}


static void encoder_state_entry_points_explore(const encoder_state_t * const state, int * const r_count, int * const r_max_length) {
  int i;
  for (i = 0; state->children[i].encoder_control; ++i) {
    if (state->children[i].is_leaf) {
      const int my_length = kvz_bitstream_tell(&state->children[i].stream)/8;
      ++(*r_count);
      if (my_length > *r_max_length) {
        *r_max_length = my_length;
      }
    } else {
      encoder_state_entry_points_explore(&state->children[i], r_count, r_max_length);
    }
  }
}

static void encoder_state_write_bitstream_entry_points_write(bitstream_t * const stream, const encoder_state_t * const state, const int num_entry_points, const int write_length, int * const r_count) {
  int i;
  for (i = 0; state->children[i].encoder_control; ++i) {
    if (state->children[i].is_leaf) {
      const int my_length = kvz_bitstream_tell(&state->children[i].stream)/8;
      ++(*r_count);
      //Don't write the last one
      if (*r_count < num_entry_points) {
        WRITE_U(stream, my_length - 1, write_length, "entry_point_offset-minus1")
      }
    } else {
      encoder_state_write_bitstream_entry_points_write(stream, &state->children[i], num_entry_points, write_length, r_count);
    }
  }
}

void kvz_encoder_state_write_bitstream_slice_header(encoder_state_t * const state)
{
  // ***********************************************
  // Modified for SHVC
  const encoder_control_t * const encoder = state->encoder_control;
  bitstream_t * const stream = &state->stream;
  int j;
  int ref_negative = 0;
  int ref_positive = 0;
  if (encoder->cfg->gop_len) {
    for (j = 0; j < state->frame->ref->used_size; j++) {
      if (state->frame->ref->pocs[j] < state->frame->poc) {
        ref_negative++;
      } else {
        ref_positive++;
      }
    }
  } else ref_negative = state->frame->ref->used_size;
  //TODO: Make a better implementation
  //if( ref_negative > 0 && state->frame->ref->pocs[state->frame->ref->used_size-1] == state->frame->poc) --ref_negative;
  // ***********************************************
  
#ifdef KVZ_DEBUG
  printf("=========== Slice ===========\n");
#endif
  WRITE_U(stream, (state->slice->start_in_rs == 0), 1, "first_slice_segment_in_pic_flag");

  if (state->frame->pictype >= KVZ_NAL_BLA_W_LP
      && state->frame->pictype <= KVZ_NAL_RSV_IRAP_VCL23) {
    WRITE_U(stream, 1, 1, "no_output_of_prior_pics_flag");
  }

  // ***********************************************
  // Modified for SHVC
  // parameter set ids are global (layer id does not matter)
  // so set the correct pps id (use layer id as the pset id)
  WRITE_UE(stream, state->layer->layer_id, "slice_pic_parameter_set_id");
  // ***********************************************
  
  if (state->slice->start_in_rs > 0) {
    //For now, we don't support dependent slice segments
    //WRITE_U(stream, 0, 1, "dependent_slice_segment_flag");
    WRITE_UE(stream, state->slice->start_in_rs, "slice_segment_address");
  }

  WRITE_UE(stream, state->frame->slicetype, "slice_type");

  // if !entropy_slice_flag

    //if output_flag_present_flag
      //WRITE_U(stream, 1, 1, "pic_output_flag");
    //end if
    //if( IdrPicFlag ) <- nal_unit_type == 5

  // ***********************************************
  // Modified for SHVC
  // TODO: Get proper values.
  uint8_t poc_lsb_not_present_flag = 0; //Infered to be 0 if not present
  if ((state->layer->layer_id > 0 && !poc_lsb_not_present_flag) || (state->frame->pictype != KVZ_NAL_IDR_W_RADL
      && state->frame->pictype != KVZ_NAL_IDR_N_LP)) {  
    WRITE_U(stream, state->frame->poc&0x1f, 5, "pic_order_cnt_lsb");
  }
  if (state->frame->pictype != KVZ_NAL_IDR_W_RADL
      && state->frame->pictype != KVZ_NAL_IDR_N_LP) {
    int last_poc = 0;
    int poc_shift = 0;

    //WRITE_U(stream, state->frame->poc & 0x1f, 5, "pic_order_cnt_lsb");
//***********************************************
    
    uint8_t num_short_term_ref_pic_sets = state->encoder_control->cfg->ref_frames; //TODO: get this number from somewhere else
    if( state->layer->layer_id > 0 ) num_short_term_ref_pic_sets = num_short_term_ref_pic_sets > 1 ? num_short_term_ref_pic_sets-1 : 1; //Reserve a "reference" for IRL
    uint8_t ref_sets = state->layer->max_layers > 1 ? 1 : 0; //TODO: Remove if pointless
    WRITE_U(stream, ref_sets, 1, "short_term_ref_pic_set_sps_flag");
    if (ref_sets == 0) {
      WRITE_UE(stream, ref_negative, "num_negative_pics");
      WRITE_UE(stream, ref_positive, "num_positive_pics");
    for (j = 0; j < ref_negative; j++) {      
      int8_t delta_poc = 0;
      
      if (encoder->cfg->gop_len) {
        int8_t found = 0;
        do {
          delta_poc = encoder->cfg->gop[state->frame->gop_offset].ref_neg[j + poc_shift];
          for (int i = 0; i < state->frame->ref->used_size; i++) {
            if (state->frame->ref->pocs[i] == state->frame->poc - delta_poc) {
              found = 1;
              break;
              }
            }
            if (!found) poc_shift++;
            if (j + poc_shift == ref_negative) {
              fprintf(stderr, "Failure, reference not found!");
              exit(EXIT_FAILURE);
            }
          } while (!found);
        }

        WRITE_UE(stream, encoder->cfg->gop_len ? delta_poc - last_poc - 1 : 0, "delta_poc_s0_minus1");
        last_poc = delta_poc;
        WRITE_U(stream, 1, 1, "used_by_curr_pic_s0_flag");
      }
      last_poc = 0;
      poc_shift = 0;
      for (j = 0; j < ref_positive; j++) {
        int8_t delta_poc = 0;



      if (encoder->cfg->gop_len) {
        int8_t found = 0;
        do {
          delta_poc = encoder->cfg->gop[state->frame->gop_offset].ref_pos[j + poc_shift];
          for (int i = 0; i < state->frame->ref->used_size; i++) {
            if (state->frame->ref->pocs[i] == state->frame->poc + delta_poc) {
              found = 1;
              break;
            }
          }
            if (!found) poc_shift++;
            if (j + poc_shift == ref_positive) {
              fprintf(stderr, "Failure, reference not found!");
              exit(EXIT_FAILURE);
            }
          } while (!found);
        }

        WRITE_UE(stream, encoder->cfg->gop_len ? delta_poc - last_poc - 1 : 0, "delta_poc_s1_minus1");
        last_poc = delta_poc;
        WRITE_U(stream, 1, 1, "used_by_curr_pic_s1_flag");
      }
      //WRITE_UE(stream, 0, "short_term_ref_pic_set_idx");
    }
    else if( num_short_term_ref_pic_sets > 1 ){
      uint32_t poc = state->frame->poc;
      //int stRpsIdx = num_short_term_ref_pic_sets <= poc ? 0 : num_short_term_ref_pic_sets - poc; //stRpsIdx==0 should be the one with max ref, so for the first frames use num_short_term_ref_pic_sets - Poc
      int stRpsIdx = (num_short_term_ref_pic_sets <= poc ? num_short_term_ref_pic_sets : poc) - 1; // num of (available) refs increases with idx/poc
      WRITE_U(stream, stRpsIdx, kvz_math_ceil_log2(num_short_term_ref_pic_sets), "short_term_ref_pic_set_idx"); //TODO: Get correct idx from somewhere
    }
    
    if (state->encoder_control->cfg->tmvp_enable) {
      WRITE_U(stream, 1, 1, "slice_temporal_mvp_enabled_flag");
    }
  }

    //end if
  //end if

  // ***********************************************
  // Modified for SHVC
  // TODO: Get proper values.
  if (state->layer->layer_id > 0) {  //&& !default_ref_layers_active_flag && NumDirectRefLayers[nuh_layer_id] > 0 ){
    //default_ref_layers_active_flag == 0; NumDirectRefLayers[1] == 1
    uint8_t inter_layer_pred_enabled_flag = state->frame->slicetype != KVZ_SLICE_I; //TODO: A better way?
    WRITE_U(stream, inter_layer_pred_enabled_flag, 1, "inter_layer_pred_enabled_flag");
    
    if( inter_layer_pred_enabled_flag /* && NumDirectRefLayers[nu_layer_id] > 1 */) {
      //TODO: Write stuff if there are more than one inter layer ref for this slice
    }
  }
  

  if (encoder->sao_enable) {
    WRITE_U(stream, 1, 1, "slice_sao_luma_flag");
    if (encoder->chroma_format != KVZ_CSP_400) {
      WRITE_U(stream, 1, 1, "slice_sao_chroma_flag");
    }
  }

  if (state->frame->slicetype != KVZ_SLICE_I) {
      WRITE_U(stream, 1, 1, "num_ref_idx_active_override_flag");
      WRITE_UE(stream, ref_negative != 0 ? ref_negative - 1: 0, "num_ref_idx_l0_active_minus1");
        if (state->frame->slicetype == KVZ_SLICE_B) {
          WRITE_UE(stream, ref_positive != 0 ? ref_positive - 1 : 0, "num_ref_idx_l1_active_minus1");
    }
    //TODO: Make a better check?
    //if( lists_modification_present_flag && NumPicTotalCurr>1 )
    if( state->layer->list_modification_present_flag && state->frame->ref->used_size > 1) {
    //  ref_pic_lists_modification()
      uint8_t ref_pic_lists_modification_flag_l0 = 1;
      WRITE_U(stream, ref_pic_lists_modification_flag_l0, 1,"ref_pic_list_modification_flag_l0");
      if( ref_pic_lists_modification_flag_l0 ) {
        for (int i = 0; i < ref_negative; ++i) {
          //We want to move the ILR pic first
          WRITE_U(stream, (ref_negative+i-1)%ref_negative, kvz_math_ceil_log2(state->frame->ref->used_size), "list_entry_l0[i]");
        }
      }
      if (state->frame->slicetype == KVZ_SLICE_B) {
        uint8_t ref_pic_lists_modification_flag_l1 = 0;
        WRITE_U(stream, ref_pic_lists_modification_flag_l1, 1, "ref_pic_list_modification_flag_l0");
        if (ref_pic_lists_modification_flag_l1) {
          for (int i = 0; i < ref_positive; ++i) {
           WRITE_U(stream, i, kvz_math_ceil_log2(state->frame->ref->used_size), "list_entry_l1[i]");
          }
        }
      }
    }
    //
    if (state->frame->slicetype == KVZ_SLICE_B) {
      WRITE_U(stream, 0, 1, "mvd_l1_zero_flag");
    }

      // ToDo: handle B-frames with TMVP
      if (state->encoder_control->cfg->tmvp_enable && ref_negative > 1) {
        WRITE_UE(stream, 0, "collocated_ref_idx");
      }

      WRITE_UE(stream, 5-MRG_MAX_NUM_CANDS, "five_minus_max_num_merge_cand");
  }
  //***********************************************
  {
    int slice_qp_delta = state->frame->QP - encoder->cfg->qp;
    WRITE_SE(stream, slice_qp_delta, "slice_qp_delta");
  }
   
  if (encoder->tiles_enable || encoder->wpp) {
    int num_entry_points = 0;
    int max_length_seen = 0;
    
    encoder_state_entry_points_explore(state, &num_entry_points, &max_length_seen);
    
    WRITE_UE(stream, num_entry_points - 1, "num_entry_point_offsets");
    if (num_entry_points > 0) {
      int entry_points_written = 0;
      int offset_len = kvz_math_floor_log2(max_length_seen) + 1;
      WRITE_UE(stream, offset_len - 1, "offset_len_minus1");
      encoder_state_write_bitstream_entry_points_write(stream, state, num_entry_points, offset_len, &entry_points_written); 
    }
  }
}


/**
 * \brief Add a checksum SEI message to the bitstream.
 * \param encoder The encoder.
 * \returns Void
 */
static void add_checksum(encoder_state_t * const state)
{
  bitstream_t * const stream = &state->stream;
  const videoframe_t * const frame = state->tile->frame;
  unsigned char checksum[3][SEI_HASH_MAX_LENGTH];

  // ***********************************************
  // Modified for SHVC
  kvz_nal_ext_write(stream, KVZ_NAL_SUFFIX_SEI_NUT, 0, 0, state->layer->layer_id);
  // ***********************************************

  WRITE_U(stream, 132, 8, "sei_type");

  int num_colors = (state->encoder_control->chroma_format == KVZ_CSP_400 ? 1 : 3);

  switch (state->encoder_control->cfg->hash)
  {
  case KVZ_HASH_CHECKSUM:
    kvz_image_checksum(frame->rec, checksum, state->encoder_control->bitdepth);

    WRITE_U(stream, 1 + num_colors * 4, 8, "size");
    WRITE_U(stream, 2, 8, "hash_type");  // 2 = checksum

    for (int i = 0; i < num_colors; ++i) {
      uint32_t checksum_val = (
        (checksum[i][0] << 24) + (checksum[i][1] << 16) +
        (checksum[i][2] << 8) + (checksum[i][3]));
      WRITE_U(stream, checksum_val, 32, "picture_checksum");
      CHECKPOINT("checksum[%d] = %u", i, checksum_val);
    }

    break;

  case KVZ_HASH_MD5:
    kvz_image_md5(frame->rec, checksum, state->encoder_control->bitdepth);

    WRITE_U(stream, 1 + num_colors * 16, 8, "size");
    WRITE_U(stream, 0, 8, "hash_type");  // 0 = md5

    for (int i = 0; i < num_colors; ++i) {
      for (int b = 0; b < 16; ++b) {
        WRITE_U(stream, checksum[i][b], 8, "picture_md5");
      }
    }

    break;

  case KVZ_HASH_NONE:
    // Means we shouldn't be writing this SEI.
    assert(0);
  }

  kvz_bitstream_align(stream);

  // spec:sei_rbsp() rbsp_trailing_bits
  kvz_bitstream_add_rbsp_trailing_bits(stream);
}

/**
 * \brief Move child state bitstreams to the parent stream.
 */
static void encoder_state_write_bitstream_children(encoder_state_t * const state)
{
  for (int i = 0; state->children[i].encoder_control; ++i) {
    kvz_encoder_state_write_bitstream(&state->children[i]);
    kvz_bitstream_move(&state->stream, &state->children[i].stream);
  }
}

static void encoder_state_write_bitstream_main(encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;
  bitstream_t * const stream = &state->stream;
  uint64_t curpos = kvz_bitstream_tell(stream);

  // The first NAL unit of the access unit must use a long start code.
  bool first_nal_in_au = true;

  // Access Unit Delimiter (AUD)
  if (encoder->aud_enable) {
    first_nal_in_au = false;
    encoder_state_write_bitstream_aud(state);
  }
  
  if ((encoder->vps_period > 0 && state->frame->num % encoder->vps_period == 0)
      || (state->frame->num == 0 && encoder->vps_period >= 0))
  {
    first_nal_in_au = false;
    kvz_encoder_state_write_parameter_sets(&state->stream, state);
  }
  
  // ***********************************************
  // Modified for SHVC. TODO: only in base layer?
  // Send Kvazaar version information only in the first frame.
  if (state->frame->num == 0 && encoder->cfg->add_encoder_info && state->layer->layer_id == 0) {
    kvz_nal_write(stream, KVZ_NAL_PREFIX_SEI_NUT, 0, first_nal_in_au);
    encoder_state_write_bitstream_prefix_sei_version(state);

    // spec:sei_rbsp() rbsp_trailing_bits
    kvz_bitstream_add_rbsp_trailing_bits(stream);
  }
 

  //SEI messages for interlacing
  if (encoder->vui.frame_field_info_present_flag){
    // These should be optional, needed for earlier versions
    // of HM decoder to accept bitstream
    //kvz_nal_write(stream, KVZ_NAL_PREFIX_SEI_NUT, 0, 0);
    //encoder_state_write_active_parameter_sets_sei_message(state);
    //kvz_bitstream_rbsp_trailing_bits(stream);

    kvz_nal_write(stream, KVZ_NAL_PREFIX_SEI_NUT, 0, first_nal_in_au);
    encoder_state_write_picture_timing_sei_message(state);

    // spec:sei_rbsp() rbsp_trailing_bits
    kvz_bitstream_add_rbsp_trailing_bits(stream);
  }

  {
    uint8_t nal_type = (state->frame->is_idr_frame ? KVZ_NAL_IDR_W_RADL : KVZ_NAL_TRAIL_R);
    kvz_nal_ext_write(stream, nal_type, 0, first_nal_in_au, state->layer->layer_id);
  }
 // ***********************************************
  {
    PERFORMANCE_MEASURE_START(KVZ_PERF_FRAME);
    encoder_state_write_bitstream_children(state);
    PERFORMANCE_MEASURE_END(KVZ_PERF_FRAME, encoder->threadqueue, "type=write_bitstream_append,frame=%d,encoder_type=%c", state->frame->num, state->type);
  }
  
  if (state->encoder_control->cfg->hash != KVZ_HASH_NONE) {
    PERFORMANCE_MEASURE_START(KVZ_PERF_FRAME);
    // Calculate checksum
    add_checksum(state);
    PERFORMANCE_MEASURE_END(KVZ_PERF_FRAME, encoder->threadqueue, "type=write_bitstream_checksum,frame=%d,encoder_type=%c", state->frame->num, state->type);
  }
  
  //Get bitstream length for stats
  uint64_t newpos = kvz_bitstream_tell(stream);
  state->stats_bitstream_length = (newpos >> 3) - (curpos >> 3);

  if (state->frame->num > 0) {
    state->frame->total_bits_coded = state->previous_encoder_state->frame->total_bits_coded;
  }
  state->frame->total_bits_coded += newpos - curpos;

  if (encoder->cfg->gop_len > 0 && state->frame->gop_offset > 0) {
    state->frame->cur_gop_bits_coded = state->previous_encoder_state->frame->cur_gop_bits_coded;
  } else {
    state->frame->cur_gop_bits_coded = 0;
  }
  state->frame->cur_gop_bits_coded += newpos - curpos;
}

void kvz_encoder_state_write_bitstream_leaf(encoder_state_t * const state)
{
  const encoder_control_t * const encoder = state->encoder_control;
  //Write terminator of the leaf
  assert(state->is_leaf);
  
  //Last LCU
  {
    const lcu_order_element_t * const lcu = &state->lcu_order[state->lcu_order_count - 1];
    const int lcu_addr_in_ts = lcu->id + state->tile->lcu_offset_in_ts;
    const int end_of_slice_segment_flag = kvz_lcu_at_slice_end(encoder, lcu_addr_in_ts);
  
    kvz_cabac_encode_bin_trm(&state->cabac, end_of_slice_segment_flag);  // end_of_slice_segment_flag
  
    if (!end_of_slice_segment_flag) {
      assert(kvz_lcu_at_tile_end(encoder, lcu_addr_in_ts) || lcu->position.x == (state->tile->frame->width_in_lcu - 1));
      kvz_cabac_encode_bin_trm(&state->cabac, 1); // end_of_sub_stream_one_bit == 1
      kvz_cabac_flush(&state->cabac);
    } else {
      kvz_cabac_flush(&state->cabac);
      kvz_bitstream_align_zero(&state->stream);
    }
  }
}

void kvz_encoder_state_worker_write_bitstream_leaf(void * opaque)
{
  kvz_encoder_state_write_bitstream_leaf((encoder_state_t *) opaque);
}

static void encoder_state_write_bitstream_tile(encoder_state_t * const state)
{
  encoder_state_write_bitstream_children(state);
}

static void encoder_state_write_bitstream_slice(encoder_state_t * const state)
{
  kvz_encoder_state_write_bitstream_slice_header(state);
  kvz_bitstream_add_rbsp_trailing_bits(&state->stream);
  encoder_state_write_bitstream_children(state);
}

void kvz_encoder_state_write_bitstream(encoder_state_t * const state)
{
  if (!state->is_leaf) {
    switch (state->type) {
      case ENCODER_STATE_TYPE_MAIN:
        encoder_state_write_bitstream_main(state);
        break;
      case ENCODER_STATE_TYPE_TILE:
        encoder_state_write_bitstream_tile(state);
        break;
      case ENCODER_STATE_TYPE_SLICE:
        encoder_state_write_bitstream_slice(state);
        break;
      default:
        fprintf(stderr, "Unsupported node type %c!\n", state->type);
        assert(0);
    }
  }
}

void kvz_encoder_state_worker_write_bitstream(void * opaque)
{
  kvz_encoder_state_write_bitstream((encoder_state_t *) opaque);
}

void kvz_encoder_state_write_parameter_sets(bitstream_t *stream,
                                            encoder_state_t * const state)
{
  // ***********************************************
  // Modified for SHVC
  // SPS and PPS need layer id in nalu header. TODO: Write VPS and SPS only when first in sequence?

  // Video Parameter Set (VPS)
  if (state->layer->layer_id == 0) //only call for base layer
  {
    kvz_nal_write(stream, KVZ_NAL_VPS_NUT, 0, 1);
    encoder_state_write_bitstream_vid_parameter_set(stream, state);
  }

  // Sequence Parameter Set (SPS)
  kvz_nal_ext_write(stream, KVZ_NAL_SPS_NUT, 0, 1, state->layer->layer_id);
  encoder_state_write_bitstream_seq_parameter_set(stream, state);

  // Picture Parameter Set (PPS)
  kvz_nal_ext_write(stream, KVZ_NAL_PPS_NUT, 0, 1, state->layer->layer_id);
  encoder_state_write_bitstream_pic_parameter_set(stream, state);
  //***********************************************
}
