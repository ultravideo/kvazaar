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

#include "encoder.h"

// This define is required for M_PI on Windows.
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cfg.h"
#include "strategyselector.h"


/**
 * \brief Strength of QP adjustments when using adaptive QP for 360 video.
 *
 * Determined empirically.
 */
static const double ERP_AQP_STRENGTH = 3.0;


static int encoder_control_init_gop_layer_weights(encoder_control_t * const);

static int size_of_wpp_ends(int threads)
{
  // Based on the shape of the area where all threads can't yet run in parallel.
  return 4 * threads * threads - 2 * threads;
}

static int select_owf_auto(const kvz_config *const cfg)
{
  if (cfg->intra_period == 1) {
    if (cfg->wpp) {
      // If wpp is on, select owf such that less than 15% of the
      // frame is covered by the are threads can not work at the same time.
      const int lcu_width = CEILDIV(cfg->width, LCU_WIDTH);
      const int lcu_height = CEILDIV(cfg->height, LCU_WIDTH);

      // Find the largest number of threads per frame that satifies the
      // the condition: wpp start/stop inefficiency takes up  less than 15%
      // of frame area.
      int threads_per_frame = 1;
      const int wpp_treshold = lcu_width * lcu_height * 15 / 100;
      while ((threads_per_frame + 1) * 2 < lcu_width &&
        threads_per_frame + 1 < lcu_height &&
        size_of_wpp_ends(threads_per_frame + 1) < wpp_treshold) {
        ++threads_per_frame;
      }

      const int threads = MAX(cfg->threads, 1);
      const int frames = CEILDIV(threads, threads_per_frame);

      // Convert from number of parallel frames to number of additional frames.
      return CLIP(0, threads - 1, frames - 1);
    } else {
      // If wpp is not on, select owf such that there is enough
      // tiles for twice the number of threads.

      int tiles_per_frame = cfg->tiles_width_count * cfg->tiles_height_count;
      int threads = (cfg->threads > 1 ? cfg->threads : 1);
      int frames = CEILDIV(threads * 4, tiles_per_frame);

      // Limit number of frames to 1.25x the number of threads for the case
      // where there is only 1 tile per frame.
      frames = CLIP(1, threads * 4 / 3, frames);
      return frames - 1;
    }
  } else {
    // Try and estimate a good number of parallel frames for inter.
    const int lcu_width = CEILDIV(cfg->width, LCU_WIDTH);
    const int lcu_height = CEILDIV(cfg->height, LCU_WIDTH);
    int threads_per_frame = MIN(lcu_width / 2, lcu_height);
    int threads = cfg->threads;

    // If all threads fit into one frame, at least two parallel frames should
    // be used to reduce the effect of WPP spin-up and wind-down.
    int frames = 1;

    while (threads > 0 && threads_per_frame > 0) {
      frames += 1;
      threads -= threads_per_frame;
      threads_per_frame -= 2;
    }

    if (cfg->gop_len && cfg->gop_lowdelay && cfg->gop_lp_definition.t > 1) {
      // Temporal skipping makes every other frame very fast to encode so
      // more parallel frames should be used.
      frames *= 2;
    }
    return CLIP(0, cfg->threads * 2 - 1, frames - 1);
  }
}


static unsigned cfg_num_threads(void)
{
  unsigned cpus = kvz_g_hardware_flags.physical_cpu_count;
  unsigned fake_cpus = kvz_g_hardware_flags.logical_cpu_count - cpus;

  // Default to 4 if we don't know the number of CPUs.
  if (cpus == 0) return 4;

  // 1.5 times the number of physical cores seems to be a good compromise
  // when hyperthreading is available on Haswell.
  return cpus + fake_cpus / 2;
}

// ***********************************************
  // Modified for SHVC.
/**
 * \brief Return weight for 360 degree ERP video
 *
 * Returns the scaling factor of area from equirectangular projection to
 * spherical surface.
 *
 * \param y   y-coordinate of the pixel
 * \param h   height of the picture
 */
static double ws_weight(int y, int h)
{
  return cos((y - 0.5 * h + 0.5) * (M_PI / h));
}



/**
 * \brief Update ROI QPs for 360 video with equirectangular projection.
 *
 * Writes updated ROI parameters to encoder->cfg.roi.
 *
 * \param encoder       encoder control
 * \param orig_roi      original delta QPs or NULL
 * \param orig_width    width of orig_roi
 * \param orig_height   height of orig_roi
 */
static void init_erp_aqp_roi(encoder_control_t* encoder,
                             int8_t *orig_roi,
                             int32_t orig_width,
                             int32_t orig_height)
{
  // Update ROI with WS-PSNR delta QPs.
  int height = encoder->in.height_in_lcu;
  int width  = orig_roi ? orig_width : 1;

  int frame_height = encoder->in.real_height;

  encoder->cfg.roi.width  = width;
  encoder->cfg.roi.height = height;
  encoder->cfg.roi.dqps   = calloc(width * height, sizeof(orig_roi[0]));

  double total_weight = 0.0;
  for (int y = 0; y < frame_height; y++) {
    total_weight += ws_weight(y, frame_height);
  }

  for (int y_lcu = 0; y_lcu < height; y_lcu++) {
    int y_orig = LCU_WIDTH * y_lcu;
    int lcu_height = MIN(LCU_WIDTH, frame_height - y_orig);

    double lcu_weight = 0.0;
    for (int y = y_orig; y < y_orig + lcu_height; y++) {
      lcu_weight += ws_weight(y, frame_height);
    }
    // Normalize.
    lcu_weight = (lcu_weight * frame_height) / (total_weight * lcu_height);

    int8_t qp_delta = round(-ERP_AQP_STRENGTH * log2(lcu_weight));

    if (orig_roi) {
      // If a ROI array already exists, we copy the existing values to the
      // new array while adding qp_delta to each.
      int y_roi = y_lcu * orig_height / height;
      for (int x = 0; x < width; x++) {
        encoder->cfg.roi.dqps[x + y_lcu * width] =
          CLIP(-51, 51, orig_roi[x + y_roi * width] + qp_delta);
      }

    } else {
      // Otherwise, simply write qp_delta to the ROI array.
      encoder->cfg.roi.dqps[y_lcu] = qp_delta;
    }
  }
}


/**
 * \brief Allocate and initialize an encoder control structure.
 *
 * Initialize an encoder control structure for each layer in the cfg and
 * connect them through the next_enc_ctrl field. 
 *
 * \param cfg   encoder configuration
 * \return      initialized encoder control or NULL on failure
 */
encoder_control_t* kvz_encoder_control_init(const kvz_config *cfg)
{
  encoder_control_t *encoder = NULL;
  encoder_control_t *prev_enc = NULL;
  encoder_control_t *first_enc = NULL;

  if (!cfg) {
    fprintf(stderr, "Config object must not be null!\n");
    goto init_failed;
  }

  for (; cfg != NULL; cfg = cfg->next_cfg ) {

    // Make sure that the parameters make sense.
    if (!kvz_config_validate(cfg)) {
      goto init_failed;
    }

    encoder = calloc(1, sizeof(encoder_control_t));
    
    if (!encoder) {
      fprintf(stderr, "Failed to allocate encoder control.\n");
      goto init_failed;
    }

    if(prev_enc != NULL) {
      prev_enc->next_enc_ctrl = encoder;
    }
    else if(first_enc == NULL) {
      first_enc = encoder;
    }

    // Take a copy of the config.
    memcpy(&encoder->cfg, cfg, sizeof(encoder->cfg));
    // Set fields that are not copied to NULL.
    encoder->cfg.cqmfile = NULL;
    encoder->cfg.tiles_width_split = NULL;
    encoder->cfg.tiles_height_split = NULL;
    encoder->cfg.slice_addresses_in_ts = NULL;


    encoder->cfg.max_layers = NULL;
    encoder->cfg.max_input_layers = NULL;
    encoder->cfg.input_widths = NULL;
    encoder->cfg.input_heights = NULL;
    encoder->cfg.next_cfg = NULL;


    if (encoder->cfg.threads == -1) {
      encoder->cfg.threads = cfg_num_threads();
    }

    if (encoder->cfg.gop_len > 0) {
      if (encoder->cfg.gop_lowdelay) {
        kvz_config_process_lp_gop(&encoder->cfg);
      }
    }

    // Need to set owf before initializing threadqueue.
    if (encoder->cfg.owf < 0) {
      encoder->cfg.owf = select_owf_auto(&encoder->cfg);
      fprintf(stderr, "--owf=auto value set to %d.\n", encoder->cfg.owf);
    }
    if (encoder->cfg.source_scan_type != KVZ_INTERLACING_NONE) {
      // If using interlaced coding with OWF, the OWF has to be an even number
      // to ensure that the pair of fields will be output for the same picture.
      if (encoder->cfg.owf % 2 == 1) {
        encoder->cfg.owf += 1;
      }
    }

    encoder->threadqueue = MALLOC(threadqueue_queue_t, 1);
    if (!encoder->threadqueue ||
      !kvz_threadqueue_init(encoder->threadqueue,
      encoder->cfg.threads,
      encoder->cfg.owf > 0)) {
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
    if (cfg->cqmfile) {
      FILE* cqmfile = fopen(cfg->cqmfile, "rb");
      if (cqmfile) {
        kvz_scalinglist_parse(&encoder->scaling_list, cqmfile);
        fclose(cqmfile);
      }
      else {
        fprintf(stderr, "Could not open CQM file.\n");
        goto init_failed;
      }
    }
    kvz_scalinglist_process(&encoder->scaling_list, encoder->bitdepth);
    
    kvz_encoder_control_input_init(encoder, encoder->cfg.width, encoder->cfg.height);

    if (encoder->cfg.framerate_num != 0) {
      double framerate = encoder->cfg.framerate_num / (double)encoder->cfg.framerate_denom;
      encoder->target_avg_bppic = encoder->cfg.target_bitrate / framerate;
    }
    else {
      encoder->target_avg_bppic = encoder->cfg.target_bitrate / encoder->cfg.framerate;
    }
    encoder->target_avg_bpp = encoder->target_avg_bppic / encoder->in.pixels_per_pic;

    if (!encoder_control_init_gop_layer_weights(encoder)) {
      goto init_failed;
    }
    
  if (cfg->erp_aqp) {
    init_erp_aqp_roi(encoder,
                     cfg->roi.dqps,
                     cfg->roi.width,
                     cfg->roi.height);

  } else if (cfg->roi.dqps) {
    // Copy delta QP array for ROI coding.
    const size_t roi_size = encoder->cfg.roi.width * encoder->cfg.roi.height;
    encoder->cfg.roi.dqps = calloc(roi_size, sizeof(cfg->roi.dqps[0]));
    memcpy(encoder->cfg.roi.dqps,
           cfg->roi.dqps,
           roi_size * sizeof(*cfg->roi.dqps));

  }

  encoder->lcu_dqp_enabled = cfg->target_bitrate > 0 || encoder->cfg.roi.dqps;

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

      }
      else if (encoder->cfg.tiles_height_count > encoder->in.height_in_lcu) {
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
            (i + 1) * encoder->in.width_in_lcu / encoder->cfg.tiles_width_count -
            i    * encoder->in.width_in_lcu / encoder->cfg.tiles_width_count;
        }
      }
      else {
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
          tiles_row_height[i] = ((i + 1) * encoder->in.height_in_lcu) / encoder->cfg.tiles_height_count -
            i * encoder->in.height_in_lcu / encoder->cfg.tiles_height_count;
        }
      }
      else {
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
        tiles_col_bd[i + 1] = tiles_col_bd[i] + tiles_col_width[i];
      }

      //(6-6) in ITU-T Rec. H.265 (04/2013)
      tiles_row_bd[0] = 0;
      for (int i = 0; i < encoder->cfg.tiles_height_count; ++i) {
        tiles_row_bd[i + 1] = tiles_row_bd[i] + tiles_row_height[i];
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
          for (int y = tiles_row_bd[j]; y < tiles_row_bd[j + 1]; ++y) {
            for (int x = tiles_col_bd[i]; x < tiles_col_bd[i + 1]; ++x) {
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

      }
      else {
        int *slice_addresses_in_ts;
        encoder->slice_count = encoder->cfg.slice_count;
        if (encoder->slice_count == 0) {
          encoder->slice_count = 1;

          encoder->slice_addresses_in_ts = slice_addresses_in_ts =
            MALLOC(int, encoder->slice_count);
          if (!slice_addresses_in_ts) goto init_failed;

          slice_addresses_in_ts[0] = 0;

        }
        else {
          encoder->slice_addresses_in_ts = slice_addresses_in_ts =
            MALLOC(int, encoder->slice_count);
          if (!slice_addresses_in_ts) goto init_failed;

          if (!cfg->slice_addresses_in_ts) {
            slice_addresses_in_ts[0] = 0;
            for (int i = 1; i < encoder->slice_count; ++i) {
              slice_addresses_in_ts[i] = encoder->in.width_in_lcu * encoder->in.height_in_lcu * i / encoder->slice_count;
            }
          }
          else {
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
            }
      else {
        printf("Wavefront Parallel Processing: disabled\n");
      }
      printf("\n");
#endif //KVZ_DEBUG
          }

    assert(WITHIN(encoder->cfg.pu_depth_inter.min, PU_DEPTH_INTER_MIN, PU_DEPTH_INTER_MAX));
    assert(WITHIN(encoder->cfg.pu_depth_inter.max, PU_DEPTH_INTER_MIN, PU_DEPTH_INTER_MAX));
    assert(WITHIN(encoder->cfg.pu_depth_intra.min, PU_DEPTH_INTRA_MIN, PU_DEPTH_INTRA_MAX));
    assert(WITHIN(encoder->cfg.pu_depth_intra.max, PU_DEPTH_INTRA_MIN, PU_DEPTH_INTRA_MAX));

    // Disable in-loop filters, sign hiding and transform skip when using
    // lossless coding.
    if (encoder->cfg.lossless) {
      encoder->cfg.deblock_enable = false;
      encoder->cfg.sao_enable = false;
      encoder->cfg.signhide_enable = false;
      encoder->cfg.trskip_enable = false;
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
    }
    else {
      encoder->cfg.vps_period = -1;
    }

    //*********************************************
    //For scalable extension. TODO: Check that stuff from the cfg is copied properly since it is only copied now
    encoder->layer.layer_id = cfg->layer;
    encoder->layer.input_layer = cfg->input_layer;
    encoder->layer.max_layers = *cfg->max_layers;

    encoder->layer.num_layer_sets = encoder->layer.num_output_layer_sets = encoder->layer.max_layers;
    encoder->layer.list_modification_present_flag = (encoder->layer.layer_id > 0) && (cfg->ref_frames > 1) && (cfg->intra_period != 1) ? 1 : 0;
    encoder->layer.sps_ext_or_max_sub_layers_minus1 = 7;
    encoder->layer.multi_layer_ext_sps_flag = encoder->layer.layer_id != 0 && 
                                              encoder->layer.sps_ext_or_max_sub_layers_minus1 == 7;
    
    //encoder->layer.upscaling = NULL; //This is set later when the parameters have been set

    encoder->layer.input_width = (*cfg->input_widths)[encoder->layer.input_layer];
    encoder->layer.input_height = (*cfg->input_heights)[encoder->layer.input_layer];

    //Set scaling parameters
    //Prepare scaling parameters so that up/downscaling gives the correct parameters for up/downscaling from prev_layer/orig to current layer
    chroma_format_t csp = (chroma_format_t)KVZ_FORMAT2CSP(cfg->input_format);
    encoder->layer.downscaling = kvz_newScalingParameters(encoder->layer.input_width,
                                                      encoder->layer.input_height,
                                                      encoder->in.real_width,
                                                      encoder->in.real_height,
                                                      csp);
    if( prev_enc != NULL ){
      encoder->layer.upscaling = kvz_newScalingParameters(prev_enc->layer.upscaling.trgt_width,
                                                      prev_enc->layer.upscaling.trgt_height,
                                                      encoder->in.real_width,
                                                      encoder->in.real_height,
                                                      csp);
    }
    else {
      encoder->layer.upscaling = kvz_newScalingParameters(encoder->in.real_width,
                                                      encoder->in.real_height,
                                                      encoder->in.real_width,
                                                      encoder->in.real_height,
                                                      csp);
    }
    //Need to set the source (target?) to the padded size (because reasons) to conform with SHM. TODO: Trgt needs to be padded as well?
    //Scaling parameters need to be calculated for the true sizes.
    encoder->layer.upscaling.src_padding_x = (CU_MIN_SIZE_PIXELS - encoder->layer.upscaling.src_width % CU_MIN_SIZE_PIXELS) % CU_MIN_SIZE_PIXELS;
    encoder->layer.upscaling.src_padding_y = (CU_MIN_SIZE_PIXELS - encoder->layer.upscaling.src_height % CU_MIN_SIZE_PIXELS) % CU_MIN_SIZE_PIXELS;


    //*********************************************

    prev_enc = encoder;
  }

  return first_enc;

init_failed:
  if(encoder != NULL) kvz_encoder_control_free(encoder);
  if(prev_enc != NULL) prev_enc->next_enc_ctrl = NULL;
  for(encoder_control_t *enc = first_enc; enc != NULL; enc = (encoder_control_t*)enc->next_enc_ctrl)
    kvz_encoder_control_free(enc);
  return NULL;
}
// ***********************************************

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

  FREE_POINTER(encoder->cfg.roi.dqps);

  kvz_scalinglist_destroy(&encoder->scaling_list);

  if (encoder->threadqueue) {
    kvz_threadqueue_finalize(encoder->threadqueue);
  }
  FREE_POINTER(encoder->threadqueue);

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
 * Only GOP structures with exactly four layers are supported.
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

    default:
      fprintf(stderr, "Unsupported number of GOP layers (%d)\n", num_layers);
      return 0;
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
