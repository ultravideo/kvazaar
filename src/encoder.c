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

#include "encoder.h"

#include <stdio.h>
#include <stdlib.h>

#include "cfg.h"
#include "gop.h"
#include "rdo.h"
#include "strategyselector.h"
#include "kvz_math.h"
#include "fast_coeff_cost.h"

static int encoder_control_init_gop_layer_weights(encoder_control_t * const);

static unsigned cfg_num_threads(void)
{
  if (kvz_g_hardware_flags.logical_cpu_count == 0) {
    // Default to 4 if we don't know the number of CPUs.
    return 4;
  }

  return kvz_g_hardware_flags.logical_cpu_count;
}


static int get_max_parallelism(const encoder_control_t *const encoder)
{
  const int width_lcu  = CEILDIV(encoder->cfg.width, LCU_WIDTH);
  const int height_lcu = CEILDIV(encoder->cfg.height, LCU_WIDTH);
  const int wpp_limit  = MIN(height_lcu, CEILDIV(width_lcu, 2));
  const int par_frames = encoder->cfg.owf + 1;

  int parallelism = 0;

  if (encoder->cfg.intra_period == 1) {
    int threads_per_frame;
    if (encoder->cfg.wpp) {
      // Usually limited by width because starting to code a CTU requires
      // that the next two CTUs in the row above have been completed.
      threads_per_frame = wpp_limit;
    } else {
      // One thread for each tile.
      threads_per_frame = encoder->cfg.tiles_width_count *
                          encoder->cfg.tiles_height_count;
    }
    // Divide by two since all frames cannot achieve the maximum
    // parallelism all the time.
    parallelism = par_frames * threads_per_frame / 2;

  } else {
    if (encoder->cfg.wpp) {
      const int last_diagonal = (width_lcu - 1) + (height_lcu - 1) * 2;

      // Index of a diagonal. The diagonal contains CTUs whose coordinates
      // satisfy x + 2*y == diagonal. We start the sum from the longest
      // diagonal.
      int diagonal = CEILDIV(last_diagonal, 2);

      // Difference between diagonal indices in consecutive frames.
      const int frame_delay = 1 + encoder->max_inter_ref_lcu.right +
                              2 * encoder->max_inter_ref_lcu.down;
      int step = frame_delay;
      int direction = -1;

      // Compute number of threads for each parallel frame.
      for (int num_frames = 0; num_frames < par_frames; num_frames++) {
        if (diagonal < 0 || diagonal > last_diagonal) {
          // No room for more threads.
          break;
        }

        // Count number of CTUs on the diagonal.
        if (diagonal < MIN(2 * height_lcu, width_lcu)) {
          parallelism += 1 + diagonal / 2;
        } else {
          parallelism += MIN(
            wpp_limit,
            height_lcu + CEILDIV(width_lcu, 2) - 1 - CEILDIV(diagonal, 2)
          );
        }
        diagonal += direction * step;
        step += frame_delay;
        direction = -direction;
      }

    } else {
      parallelism = encoder->cfg.tiles_width_count *
                    encoder->cfg.tiles_height_count;
    }
  }

  return parallelism;
}


/**
 * \brief Allocate and initialize an encoder control structure.
 *
 * \param cfg   encoder configuration
 * \return      initialized encoder control or NULL on failure
 */
encoder_control_t* kvz_encoder_control_init(const kvz_config *const cfg)
{
  encoder_control_t *encoder = NULL;

  if (!cfg) {
    fprintf(stderr, "Config object must not be null!\n");
    goto init_failed;
  }

  // Make sure that the parameters make sense.
  if (!kvz_config_validate(cfg)) {
    goto init_failed;
  }

  encoder = calloc(1, sizeof(encoder_control_t));
  if (!encoder) {
    fprintf(stderr, "Failed to allocate encoder control.\n");
    goto init_failed;
  }

  // Take a copy of the config.
  memcpy(&encoder->cfg, cfg, sizeof(encoder->cfg));
  // Set fields that are not copied to NULL.
  encoder->cfg.cqmfile = NULL;
  encoder->cfg.tiles_width_split = NULL;
  encoder->cfg.tiles_height_split = NULL;
  encoder->cfg.slice_addresses_in_ts = NULL;
  encoder->cfg.fast_coeff_table_fn = NULL;

  if (encoder->cfg.gop_len > 0) {
    if (encoder->cfg.gop_lowdelay) {
      if (encoder->cfg.gop_len == 4 && encoder->cfg.ref_frames == 4) {
        memcpy(encoder->cfg.gop, kvz_gop_lowdelay4, sizeof(kvz_gop_lowdelay4));
      } else {
        kvz_config_process_lp_gop(&encoder->cfg);
      }
    }
  } 
  
  if( encoder->cfg.intra_qp_offset_auto ) {
    // Limit offset to -3 since HM/VTM seems to use it even for 32 frame gop
    encoder->cfg.intra_qp_offset = encoder->cfg.gop_len > 1 ? MAX(-(int8_t)kvz_math_ceil_log2( encoder->cfg.gop_len ) + 1, -3) : 0;
  }

  // Disable GOP and QP offset for all-intra coding
  if (encoder->cfg.intra_period == 1) {
    encoder->cfg.gop_len = 0;
    encoder->cfg.intra_qp_offset = 0;
  }

  encoder->poc_lsb_bits = MAX(4, kvz_math_ceil_log2(encoder->cfg.gop_len * 2 + 1));

  encoder->max_inter_ref_lcu.right = 1;
  encoder->max_inter_ref_lcu.down  = 1;

  int max_threads = encoder->cfg.threads;
  if (max_threads < 0) {
    max_threads = cfg_num_threads();
  }
  max_threads = MAX(1, max_threads);

  // Need to set owf before initializing threadqueue.
  if (encoder->cfg.owf < 0) {
    int best_parallelism = 0;

    for (encoder->cfg.owf = 0; true; encoder->cfg.owf++) {
      int parallelism = get_max_parallelism(encoder);

      if (parallelism <= best_parallelism) {
        // No improvement over previous OWF.
        encoder->cfg.owf--;
        break;
      }

      best_parallelism = parallelism;
      if (parallelism >= max_threads) {
        // Cannot have more parallelism than there are threads.
        break;
      }
    }

    // Add two frames so that we have frames ready to be coded when one is
    // completed.
    encoder->cfg.owf += 2;

    if (cfg->enable_logging_output) fprintf(stderr, "--owf=auto value set to %d.\n", encoder->cfg.owf);
  }

  if (encoder->cfg.threads < 0) {
    encoder->cfg.threads = MIN(max_threads, get_max_parallelism(encoder));
    if (cfg->enable_logging_output) fprintf(stderr, "--threads=auto value set to %d.\n", encoder->cfg.threads);
  }

  if (encoder->cfg.source_scan_type != KVZ_INTERLACING_NONE) {
    // If using interlaced coding with OWF, the OWF has to be an even number
    // to ensure that the pair of fields will be output for the same picture.
    if (encoder->cfg.owf % 2 == 1) {
      encoder->cfg.owf += 1;
    }
  }

  encoder->threadqueue = kvz_threadqueue_init(encoder->cfg.threads);
  if (!encoder->threadqueue) {
    fprintf(stderr, "Could not initialize threadqueue.\n");
    goto init_failed;
  }

  encoder->bitdepth = KVZ_BIT_DEPTH;

  encoder->chroma_format = KVZ_FORMAT2CSP(encoder->cfg.input_format);

  // Interlacing
  encoder->in.source_scan_type = (int8_t)encoder->cfg.source_scan_type;
  encoder->vui.field_seq_flag = encoder->cfg.source_scan_type != 0;
  encoder->vui.frame_field_info_present_flag = encoder->cfg.source_scan_type != 0;

  // Initialize the scaling list
  kvz_scalinglist_init(&encoder->scaling_list);

  // CQM
  if (cfg->scaling_list == KVZ_SCALING_LIST_CUSTOM && cfg->cqmfile) {
    FILE* cqmfile = fopen(cfg->cqmfile, "rb");
    if (cqmfile) {
      kvz_scalinglist_parse(&encoder->scaling_list, cqmfile);
      fclose(cqmfile);
    } else {
      fprintf(stderr, "Could not open CQM file.\n");
      goto init_failed;
    }
  } else if (cfg->scaling_list == KVZ_SCALING_LIST_DEFAULT) {
    // Enable scaling lists if default lists are used
    encoder->scaling_list.enable = 1;
    encoder->scaling_list.use_default_list = 1;
  }

  // ROI / delta QP
  if (cfg->roi.file_path) {
    const char *mode[2] = { "r", "rb" };
    encoder->roi_file = fopen(cfg->roi.file_path, mode[cfg->roi.format]);
    if (!encoder->roi_file) {
      fprintf(stderr, "Could not open ROI file.\n");
      goto init_failed;
    }
  }

  if (cfg->fast_coeff_table_fn) {
    FILE *fast_coeff_table_f = fopen(cfg->fast_coeff_table_fn, "rb");
    if (fast_coeff_table_f == NULL) {
      fprintf(stderr, "Could not open fast coeff table file.\n");
      goto init_failed;
    }
    if (kvz_fast_coeff_table_parse(&encoder->fast_coeff_table, fast_coeff_table_f) != 0) {
      fprintf(stderr, "Failed to parse fast coeff table, using default\n");
      kvz_fast_coeff_use_default_table(&encoder->fast_coeff_table);
    }
    fclose(fast_coeff_table_f);
  } else {
    kvz_fast_coeff_use_default_table(&encoder->fast_coeff_table);
  }

  if (cfg->fastrd_sampling_on || cfg->fastrd_accuracy_check_on) {
    if (cfg->fastrd_learning_outdir_fn == NULL) {
      fprintf(stderr, "No output file defined for Fast RD sampling or accuracy check.\n");
      goto init_failed;
    }
    if (kvz_init_rdcost_outfiles(cfg->fastrd_learning_outdir_fn) != 0) {
      goto init_failed;
    }
  }

  kvz_scalinglist_process(&encoder->scaling_list, encoder->bitdepth);

  kvz_encoder_control_input_init(encoder, encoder->cfg.width, encoder->cfg.height);

  if (encoder->cfg.framerate_num != 0) {
    double framerate = encoder->cfg.framerate_num / (double)encoder->cfg.framerate_denom;
    encoder->target_avg_bppic = encoder->cfg.target_bitrate / framerate;
  } else {
    encoder->target_avg_bppic = encoder->cfg.target_bitrate / encoder->cfg.framerate;
  }
  encoder->target_avg_bpp = encoder->target_avg_bppic / encoder->in.pixels_per_pic;

  if (encoder->cfg.target_bitrate > 0 &&
      !encoder_control_init_gop_layer_weights(encoder))
  {
    goto init_failed;
  }

  // NOTE: When tr_depth_inter is equal to 0, the transform is still split
  // for SMP and AMP partition units.
  encoder->tr_depth_inter = 0;

  //Tiles
  encoder->tiles_enable = encoder->cfg.tiles_width_count > 1 ||
                          encoder->cfg.tiles_height_count > 1;

  {
    const int num_ctbs = encoder->in.width_in_lcu * encoder->in.height_in_lcu;

    //Temporary pointers to allow encoder fields to be const
    int32_t *tiles_col_width, *tiles_row_height, *tiles_ctb_addr_rs_to_ts, *tiles_ctb_addr_ts_to_rs, *tiles_tile_id, *tiles_col_bd, *tiles_row_bd;

    if (encoder->cfg.tiles_width_count > encoder->in.width_in_lcu) {
      fprintf(stderr, "Too many tiles (width)!\n");
      goto init_failed;

    } else if (encoder->cfg.tiles_height_count > encoder->in.height_in_lcu) {
      fprintf(stderr, "Too many tiles (height)!\n");
      goto init_failed;
    }

    //Will be (perhaps) changed later
    encoder->tiles_uniform_spacing_flag = 1;

    encoder->tiles_col_width = tiles_col_width =
      MALLOC(int32_t, encoder->cfg.tiles_width_count);
    encoder->tiles_row_height = tiles_row_height =
      MALLOC(int32_t, encoder->cfg.tiles_height_count);

    encoder->tiles_col_bd = tiles_col_bd =
      MALLOC(int32_t, encoder->cfg.tiles_width_count + 1);
    encoder->tiles_row_bd = tiles_row_bd =
      MALLOC(int32_t, encoder->cfg.tiles_height_count + 1);

    encoder->tiles_ctb_addr_rs_to_ts = tiles_ctb_addr_rs_to_ts =
      MALLOC(int32_t, num_ctbs);
    encoder->tiles_ctb_addr_ts_to_rs = tiles_ctb_addr_ts_to_rs =
      MALLOC(int32_t, num_ctbs);
    encoder->tiles_tile_id = tiles_tile_id =
      MALLOC(int32_t, num_ctbs);

    if (!tiles_col_width ||
        !tiles_row_height ||
        !tiles_row_bd ||
        !tiles_col_bd ||
        !tiles_ctb_addr_rs_to_ts ||
        !tiles_ctb_addr_ts_to_rs ||
        !tiles_tile_id) {
      goto init_failed;
    }

    //(6-3) and (6-4) in ITU-T Rec. H.265 (04/2013)
    if (!cfg->tiles_width_split) {
      for (int i = 0; i < encoder->cfg.tiles_width_count; ++i) {
        tiles_col_width[i] =
          (i+1) * encoder->in.width_in_lcu / encoder->cfg.tiles_width_count -
           i    * encoder->in.width_in_lcu / encoder->cfg.tiles_width_count;
      }
    } else {
      int32_t last_pos_in_px = 0;
      tiles_col_width[encoder->cfg.tiles_width_count - 1] = encoder->in.width_in_lcu;
      for (int i = 0; i < encoder->cfg.tiles_width_count - 1; ++i) {
        int32_t column_width_in_lcu = (cfg->tiles_width_split[i] - last_pos_in_px) / LCU_WIDTH;
        last_pos_in_px = cfg->tiles_width_split[i];
        tiles_col_width[i] = column_width_in_lcu;
        tiles_col_width[encoder->cfg.tiles_width_count - 1] -= column_width_in_lcu;
      }
      encoder->tiles_uniform_spacing_flag = 0;
    }

    if (!cfg->tiles_height_split) {
      for (int i = 0; i < encoder->cfg.tiles_height_count; ++i) {
        tiles_row_height[i] = ((i+1) * encoder->in.height_in_lcu) / encoder->cfg.tiles_height_count -
                                   i * encoder->in.height_in_lcu / encoder->cfg.tiles_height_count;
      }
    } else {
      int32_t last_pos_in_px = 0;
      tiles_row_height[encoder->cfg.tiles_height_count - 1] = encoder->in.height_in_lcu;
      for (int i = 0; i < encoder->cfg.tiles_height_count - 1; ++i) {
        int32_t row_height_in_lcu = (cfg->tiles_height_split[i] - last_pos_in_px) / LCU_WIDTH;
        last_pos_in_px = cfg->tiles_height_split[i];
        tiles_row_height[i] = row_height_in_lcu;
        tiles_row_height[encoder->cfg.tiles_height_count - 1] -= row_height_in_lcu;
      }
      encoder->tiles_uniform_spacing_flag = 0;
    }

    //(6-5) in ITU-T Rec. H.265 (04/2013)
    tiles_col_bd[0] = 0;
    for (int i = 0; i < encoder->cfg.tiles_width_count; ++i) {
      tiles_col_bd[i+1] = tiles_col_bd[i] + tiles_col_width[i];
    }

    //(6-6) in ITU-T Rec. H.265 (04/2013)
    tiles_row_bd[0] = 0;
    for (int i = 0; i < encoder->cfg.tiles_height_count; ++i) {
      tiles_row_bd[i+1] = tiles_row_bd[i] + tiles_row_height[i];
    }

    //(6-7) in ITU-T Rec. H.265 (04/2013)
    //j == ctbAddrRs
    for (int j = 0; j < num_ctbs; ++j) {
      int tileX = 0, tileY = 0;
      int tbX = j % encoder->in.width_in_lcu;
      int tbY = j / encoder->in.width_in_lcu;

      for (int i = 0; i < encoder->cfg.tiles_width_count; ++i) {
        if (tbX >= tiles_col_bd[i]) tileX = i;
      }

      for (int i = 0; i < encoder->cfg.tiles_height_count; ++i) {
        if (tbY >= tiles_row_bd[i]) tileY = i;
      }

      tiles_ctb_addr_rs_to_ts[j] = 0;
      for (int i = 0; i < tileX; ++i) {
        tiles_ctb_addr_rs_to_ts[j] += tiles_row_height[tileY] * tiles_col_width[i];
      }
      for (int i = 0; i < tileY; ++i) {
        tiles_ctb_addr_rs_to_ts[j] += encoder->in.width_in_lcu * tiles_row_height[i];
      }
      tiles_ctb_addr_rs_to_ts[j] += (tbY - tiles_row_bd[tileY]) * tiles_col_width[tileX] +
                                     tbX - tiles_col_bd[tileX];
    }

    //(6-8) in ITU-T Rec. H.265 (04/2013)
    //Make reverse map from tile scan to raster scan
    for (int j = 0; j < num_ctbs; ++j) {
      tiles_ctb_addr_ts_to_rs[tiles_ctb_addr_rs_to_ts[j]] = j;
    }

    //(6-9) in ITU-T Rec. H.265 (04/2013)
    int tileIdx = 0;
    for (int j = 0; j < encoder->cfg.tiles_height_count; ++j) {
      for (int i = 0; i < encoder->cfg.tiles_width_count; ++i) {
        for (int y = tiles_row_bd[j]; y < tiles_row_bd[j+1]; ++y) {
          for (int x = tiles_col_bd[i]; x < tiles_col_bd[i+1]; ++x) {
            tiles_tile_id[tiles_ctb_addr_rs_to_ts[y * encoder->in.width_in_lcu + x]] = tileIdx;
          }
        }
        ++tileIdx;
      }
    }

    if (encoder->cfg.slices & KVZ_SLICES_WPP) {
      // Each WPP row will be put into a dependent slice.
      encoder->pps.dependent_slice_segments_enabled_flag = 1;
    }

    //Slices
    if (encoder->cfg.slices & KVZ_SLICES_TILES) {
      // Configure a single independent slice per tile.

      int *slice_addresses_in_ts;
      encoder->slice_count = encoder->cfg.tiles_width_count * encoder->cfg.tiles_height_count;
      encoder->slice_addresses_in_ts = slice_addresses_in_ts = MALLOC(int, encoder->slice_count);

      int slice_id = 0;
      for (int tile_row = 0; tile_row < encoder->cfg.tiles_height_count; ++tile_row) {
        for (int tile_col = 0; tile_col < encoder->cfg.tiles_width_count; ++tile_col) {
          int x = tiles_col_bd[tile_col];
          int y = tiles_row_bd[tile_row];
          int rs = y * encoder->in.width_in_lcu + x;
          int ts = tiles_ctb_addr_rs_to_ts[rs];
          slice_addresses_in_ts[slice_id] = ts;
          slice_id += 1;
        }
      }

    } else {
      int *slice_addresses_in_ts;
      encoder->slice_count = encoder->cfg.slice_count;
      if (encoder->slice_count == 0) {
        encoder->slice_count = 1;

        encoder->slice_addresses_in_ts = slice_addresses_in_ts =
          MALLOC(int, encoder->slice_count);
        if (!slice_addresses_in_ts) goto init_failed;

        slice_addresses_in_ts[0] = 0;

      } else {
        encoder->slice_addresses_in_ts = slice_addresses_in_ts =
          MALLOC(int, encoder->slice_count);
        if (!slice_addresses_in_ts) goto init_failed;

        if (!cfg->slice_addresses_in_ts) {
          slice_addresses_in_ts[0] = 0;
          for (int i = 1; i < encoder->slice_count; ++i) {
            slice_addresses_in_ts[i] = encoder->in.width_in_lcu * encoder->in.height_in_lcu * i / encoder->slice_count;
          }
        } else {
          for (int i = 0; i < encoder->slice_count; ++i) {
            slice_addresses_in_ts[i] = cfg->slice_addresses_in_ts[i];
          }
        }
      }
    }

#ifdef _DEBUG_PRINT_THREADING_INFO
    printf("Tiles columns width:");
    for (int i = 0; i < encoder->cfg.tiles_width_count; ++i) {
      printf(" %d", encoder->tiles_col_width[i]);
    }
    printf("\n");
    printf("Tiles row height:");
    for (int i = 0; i < encoder->cfg.tiles_height_count; ++i) {
      printf(" %d", encoder->tiles_row_height[i]);
    }
    printf("\n");
    //Print tile index map
    for (int y = 0; y < encoder->in.height_in_lcu; ++y) {
      for (int x = 0; x < encoder->in.width_in_lcu; ++x) {
        const int lcu_id_rs = y * encoder->in.width_in_lcu + x;
        const int lcu_id_ts = encoder->tiles_ctb_addr_rs_to_ts[lcu_id_rs];
        const char slice_start = kvz_lcu_at_slice_start(encoder, lcu_id_ts) ? '|' : ' ';
        const char slice_end = kvz_lcu_at_slice_end(encoder, lcu_id_ts)  ? '|' : ' ';

        printf("%c%03d%c", slice_start, encoder->tiles_tile_id[lcu_id_ts], slice_end);
      }
      printf("\n");
    }
    printf("\n");
    if (encoder->cfg.wpp) {
      printf("Wavefront Parallel Processing: enabled\n");
    } else {
      printf("Wavefront Parallel Processing: disabled\n");
    }
    printf("\n");
#endif //KVZ_DEBUG
  }

  for( size_t i = 0; i < KVZ_MAX_GOP_LAYERS; i++ )
  {
      if( encoder->cfg.pu_depth_inter.min[i] < 0 || cfg->pu_depth_inter.max[i] < 0 ) continue;
      assert( WITHIN( encoder->cfg.pu_depth_inter.min[i], PU_DEPTH_INTER_MIN, PU_DEPTH_INTER_MAX ) );
      assert( WITHIN( encoder->cfg.pu_depth_inter.max[i], PU_DEPTH_INTER_MIN, PU_DEPTH_INTER_MAX ) );

      if( encoder->cfg.pu_depth_intra.min[i] < 0 || cfg->pu_depth_intra.max[i] < 0 ) continue;
      assert( WITHIN( encoder->cfg.pu_depth_intra.min[i], PU_DEPTH_INTRA_MIN, PU_DEPTH_INTRA_MAX ) );
      assert( WITHIN( encoder->cfg.pu_depth_intra.max[i], PU_DEPTH_INTRA_MIN, PU_DEPTH_INTRA_MAX ) );
  }
  // Disable in-loop filters, sign hiding and transform skip when using
  // lossless coding.
  if (encoder->cfg.lossless) {
    encoder->cfg.deblock_enable  = false;
    encoder->cfg.sao_type        = false;
    encoder->cfg.signhide_enable = false;
    encoder->cfg.trskip_enable   = false;
  }

  // If fractional framerate is set, use that instead of the floating point framerate.
  if (encoder->cfg.framerate_num != 0) {
    encoder->vui.timing_info_present_flag = 1;
    encoder->vui.num_units_in_tick = encoder->cfg.framerate_denom;
    encoder->vui.time_scale = encoder->cfg.framerate_num;
    if (encoder->cfg.source_scan_type != KVZ_INTERLACING_NONE) {
      // when field_seq_flag=1, the time_scale and num_units_in_tick refer to
      // field rate rather than frame rate.
      encoder->vui.time_scale *= 2;
    }
  }

  if (encoder->cfg.vps_period >= 0) {
    encoder->cfg.vps_period = encoder->cfg.vps_period * encoder->cfg.intra_period;
  } else {
    encoder->cfg.vps_period = -1;
  }

  if(encoder->cfg.optional_key){
    encoder->cfg.optional_key = MALLOC(uint8_t,16);
    if (!encoder->cfg.optional_key) goto init_failed;
    memcpy(encoder->cfg.optional_key, cfg->optional_key, 16);
  }

  return encoder;

init_failed:
  kvz_encoder_control_free(encoder);
  return NULL;
}

/**
 * \brief Free an encoder control structure.
 */
void kvz_encoder_control_free(encoder_control_t *const encoder)
{
  if (!encoder) return;

  //Slices
  FREE_POINTER(encoder->slice_addresses_in_ts);

  //Tiles
  FREE_POINTER(encoder->tiles_col_width);
  FREE_POINTER(encoder->tiles_row_height);

  FREE_POINTER(encoder->tiles_col_bd);
  FREE_POINTER(encoder->tiles_row_bd);

  FREE_POINTER(encoder->tiles_ctb_addr_rs_to_ts);
  FREE_POINTER(encoder->tiles_ctb_addr_ts_to_rs);

  FREE_POINTER(encoder->tiles_tile_id);

  FREE_POINTER(encoder->cfg.roi.file_path);
  FREE_POINTER(encoder->cfg.optional_key);

  kvz_scalinglist_destroy(&encoder->scaling_list);

  kvz_threadqueue_free(encoder->threadqueue);
  encoder->threadqueue = NULL;

  kvz_close_rdcost_outfiles();

  if (encoder->roi_file) {
    fclose(encoder->roi_file);
  }

  free(encoder);
}

void kvz_encoder_control_input_init(encoder_control_t * const encoder,
                        const int32_t width, int32_t height)
{
  // Halve for interlaced content
  if (encoder->in.source_scan_type != 0) height /= 2;

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

  encoder->in.pixels_per_pic = encoder->in.width * encoder->in.height;


  #ifdef KVZ_DEBUG
  if (width != encoder->in.width || height != encoder->in.height) {
    printf("Picture buffer has been extended to be a multiple of the smallest block size:\r\n");
    printf("  Width = %d (%d), Height = %d (%d)\r\n", width, encoder->in.width, height,
           encoder->in.height);
  }
  #endif
}

/**
 * \brief Initialize GOP layer weights.
 * \return 1 on success, 0 on failure.
 *
 * Selects appropriate weights for layers according to the target bpp.
 * Only GOP structures with exactly four layers are supported with the.
 * exception of experimental GOP 16.
 */
static int encoder_control_init_gop_layer_weights(encoder_control_t * const encoder)
{

  kvz_gop_config const * const gop = encoder->cfg.gop;
  const int8_t gop_len = encoder->cfg.gop_len;

  int num_layers = 0;
  for (int i = 0; i < gop_len; ++i) {
    num_layers = MAX(gop[i].layer, num_layers);
  }

  switch (num_layers) {
    case 0:
    case 1:
      encoder->gop_layer_weights[0] = 1;
      break;

    // Use the first layers of the 4-layer weights.
    case 2:
    case 3:

    case 4:
      if (encoder->cfg.gop_lowdelay) {
        // These weights are based on http://doi.org/10.1109/TIP.2014.2336550
        // They are meant for lp-g4d3r4t1 gop, but work ok for others.
        if (encoder->target_avg_bpp <= 0.05) {
          encoder->gop_layer_weights[0] = 14;
          encoder->gop_layer_weights[1] = 3;
          encoder->gop_layer_weights[2] = 2;
          encoder->gop_layer_weights[3] = 1;
        } else if (encoder->target_avg_bpp <= 0.1) {
          encoder->gop_layer_weights[0] = 12;
          encoder->gop_layer_weights[1] = 3;
          encoder->gop_layer_weights[2] = 2;
          encoder->gop_layer_weights[3] = 1;
        } else if (encoder->target_avg_bpp <= 0.2) {
          encoder->gop_layer_weights[0] = 10;
          encoder->gop_layer_weights[1] = 3;
          encoder->gop_layer_weights[2] = 2;
          encoder->gop_layer_weights[3] = 1;
        } else {
          encoder->gop_layer_weights[0] = 6;
          encoder->gop_layer_weights[1] = 3;
          encoder->gop_layer_weights[2] = 2;
          encoder->gop_layer_weights[3] = 1;
        }
      } else {
        // These weights are from http://doi.org/10.1109/TIP.2014.2336550
        if (encoder->target_avg_bpp <= 0.05) {
          encoder->gop_layer_weights[0] = 30;
          encoder->gop_layer_weights[1] = 8;
          encoder->gop_layer_weights[2] = 4;
          encoder->gop_layer_weights[3] = 1;
        } else if (encoder->target_avg_bpp <= 0.1) {
          encoder->gop_layer_weights[0] = 25;
          encoder->gop_layer_weights[1] = 7;
          encoder->gop_layer_weights[2] = 4;
          encoder->gop_layer_weights[3] = 1;
        } else if (encoder->target_avg_bpp <= 0.2) {
          encoder->gop_layer_weights[0] = 20;
          encoder->gop_layer_weights[1] = 6;
          encoder->gop_layer_weights[2] = 4;
          encoder->gop_layer_weights[3] = 1;
        } else {
          encoder->gop_layer_weights[0] = 15;
          encoder->gop_layer_weights[1] = 5;
          encoder->gop_layer_weights[2] = 4;
          encoder->gop_layer_weights[3] = 1;
        }
      }
      break;
    case 5:
      if(!encoder->cfg.gop_lowdelay) {
        // These are obtained by running HM with RA GOP 16 collecting the ratio of bits spent for each
        // layer from the CTC sequences and then fitting power curve
        encoder->gop_layer_weights[0] = 13.0060187535 * pow(encoder->target_avg_bpp, -0.3727651453);
        encoder->gop_layer_weights[1] = 7.3654107392 * pow(encoder->target_avg_bpp, -0.0854329266);
        encoder->gop_layer_weights[2] = 3.6563990701 * pow(encoder->target_avg_bpp, -0.0576990493);
        encoder->gop_layer_weights[3] = 2.1486937288 * pow(encoder->target_avg_bpp, -0.0155389471);
        encoder->gop_layer_weights[4] = 1;        
      } 
      else {
        fprintf(stderr, "Unsupported amount of layers (%d) for lowdelay GOP\n", num_layers);
        return 0;
      }
      break;
    default:
      if (!encoder->cfg.gop_lowdelay && encoder->cfg.gop_len == 16) {
        fprintf(stdout, 
                "Rate control: Using experimental weights for GOP layers (%d)\n",
                num_layers);
        for (int i = 0; i < MAX_GOP_LAYERS; ++i) {
          encoder->gop_layer_weights[i] = (i == 0) ? 10 : 2;
        }
      } else {
        fprintf(stderr, "Unsupported number of GOP layers (%d)\n", num_layers);
        return 0;
      }
  }

  // Normalize weights so that the sum of weights in a GOP is one.
  double sum_weights = 0;
  for (int i = 0; i < gop_len; ++i) {
    sum_weights += encoder->gop_layer_weights[gop[i].layer - 1];
  }
  for (int i = 0; i < num_layers; ++i) {
    encoder->gop_layer_weights[i] /= sum_weights;
  }

  return 1;
}
