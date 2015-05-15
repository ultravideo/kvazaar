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

/*
 * \file
 */

#include "encoder.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "tables.h"
#include "config.h"
#include "cabac.h"
#include "image.h"
#include "nal.h"
#include "context.h"
#include "transform.h"
#include "intra.h"
#include "inter.h"
#include "filter.h"
#include "search.h"
#include "sao.h"
#include "rdo.h"

static int encoder_control_init_gop_layer_weights(encoder_control_t * const);

int encoder_control_init(encoder_control_t * const encoder, const config_t * const cfg) {
  if (!cfg) {
    fprintf(stderr, "Config object must not be null!\n");
    return 0;
  }
  
  encoder->threadqueue = MALLOC(threadqueue_queue_t, 1);
  
  encoder->owf = cfg->owf;
    
  //Init threadqueue
  if (!encoder->threadqueue || !threadqueue_init(encoder->threadqueue, cfg->threads, encoder->owf > 0)) {
    fprintf(stderr, "Could not initialize threadqueue");
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
  encoder->full_intra_search = 0;
  
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
  
  encoder_control_input_init(encoder, cfg->width, cfg->height);

  encoder->target_avg_bppic = cfg->target_bitrate / cfg->framerate;
  encoder->target_avg_bpp = encoder->target_avg_bppic / encoder->in.pixels_per_pic;

  if (!encoder_control_init_gop_layer_weights(encoder)) {
    return 0;
  }
  
  //Tiles
  encoder->tiles_enable = (encoder->cfg->tiles_width_count > 0 || encoder->cfg->tiles_height_count > 0);
  {
    int i, j; //iteration variables
    const int num_ctbs = encoder->in.width_in_lcu * encoder->in.height_in_lcu;
    int tileIdx, x, y; //iterations variable for 6-9

    //Temporary pointers to allow encoder fields to be const
    int32_t *tiles_col_width, *tiles_row_height, *tiles_ctb_addr_rs_to_ts, *tiles_ctb_addr_ts_to_rs, *tiles_tile_id, *tiles_col_bd, *tiles_row_bd;
  
    if (encoder->cfg->tiles_width_count >= encoder->in.width_in_lcu) {
      fprintf(stderr, "Too many tiles (width)!\n");
      return 0;
    } else if (encoder->cfg->tiles_height_count >= encoder->in.height_in_lcu) {
      fprintf(stderr, "Too many tiles (height)!\n");
      return 0;
    }
    
    //Will be (perhaps) changed later
    encoder->tiles_uniform_spacing_flag = 1;
    
    //tilesn[x,y] contains the number of _separation_ between tiles, whereas the encoder needs the number of tiles.
    encoder->tiles_num_tile_columns = encoder->cfg->tiles_width_count + 1;
    encoder->tiles_num_tile_rows = encoder->cfg->tiles_height_count + 1;
    
    tiles_col_width = MALLOC(int32_t, encoder->tiles_num_tile_columns);
    tiles_row_height = MALLOC(int32_t, encoder->tiles_num_tile_rows);
    
    tiles_col_bd = MALLOC(int32_t, encoder->tiles_num_tile_columns + 1);
    tiles_row_bd = MALLOC(int32_t, encoder->tiles_num_tile_rows + 1);
    
    tiles_ctb_addr_rs_to_ts = MALLOC(int32_t, num_ctbs);
    tiles_ctb_addr_ts_to_rs = MALLOC(int32_t, num_ctbs);
    
    tiles_tile_id = MALLOC(int32_t, num_ctbs);
    
    //(6-3) and (6-4) in ITU-T Rec. H.265 (04/2013)
    if (!encoder->cfg->tiles_width_split) {
      for (i=0; i < encoder->tiles_num_tile_columns; ++i) {
        tiles_col_width[i] = ((i+1) * encoder->in.width_in_lcu) / encoder->tiles_num_tile_columns -
                                  i * encoder->in.width_in_lcu / encoder->tiles_num_tile_columns;
      }
    } else {
      int32_t last_pos_in_px = 0;
      tiles_col_width[encoder->tiles_num_tile_columns-1] = encoder->in.width_in_lcu;
      for (i=0; i < encoder->tiles_num_tile_columns - 1; ++i) {
        int32_t column_width_in_lcu = (cfg->tiles_width_split[i] - last_pos_in_px) / LCU_WIDTH;
        last_pos_in_px = cfg->tiles_width_split[i];
        tiles_col_width[i] = column_width_in_lcu;
        tiles_col_width[encoder->tiles_num_tile_columns - 1] -= column_width_in_lcu;
      }
      encoder->tiles_uniform_spacing_flag = 0;
    }
    
    if (!encoder->cfg->tiles_height_split) {
      for (i=0; i < encoder->tiles_num_tile_rows; ++i) {
        tiles_row_height[i] = ((i+1) * encoder->in.height_in_lcu) / encoder->tiles_num_tile_rows -
                                   i * encoder->in.height_in_lcu / encoder->tiles_num_tile_rows;
      }
    } else {
      int32_t last_pos_in_px = 0;
      tiles_row_height[encoder->tiles_num_tile_rows-1] = encoder->in.height_in_lcu;
      for (i=0; i < encoder->tiles_num_tile_rows - 1; ++i) {
        int32_t row_height_in_lcu = (cfg->tiles_height_split[i] - last_pos_in_px) / LCU_WIDTH;
        last_pos_in_px = cfg->tiles_height_split[i];
        tiles_row_height[i] = row_height_in_lcu;
        tiles_row_height[encoder->tiles_num_tile_rows - 1] -= row_height_in_lcu;
      }
      encoder->tiles_uniform_spacing_flag = 0;
    }
    
    //(6-5) in ITU-T Rec. H.265 (04/2013)
    tiles_col_bd[0] = 0;
    for (i = 0; i < encoder->tiles_num_tile_columns; ++i) {
      tiles_col_bd[i+1] = tiles_col_bd[i] + tiles_col_width[i];
    }
    
    //(6-6) in ITU-T Rec. H.265 (04/2013)
    tiles_row_bd[0] = 0;
    for (i = 0; i < encoder->tiles_num_tile_rows; ++i) {
      tiles_row_bd[i+1] = tiles_row_bd[i] + tiles_row_height[i];
    }

    //(6-7) in ITU-T Rec. H.265 (04/2013)
    //j == ctbAddrRs
    for (j = 0; j < num_ctbs; ++j) {
      int tileX = 0, tileY = 0;
      int tbX = j % encoder->in.width_in_lcu;
      int tbY = j / encoder->in.width_in_lcu;
      
      for (i = 0; i < encoder->tiles_num_tile_columns; ++i) {
        if (tbX >= tiles_col_bd[i]) tileX = i;
      }
      
      for (i = 0; i < encoder->tiles_num_tile_rows; ++i) {
        if (tbY >= tiles_row_bd[i]) tileY = i;
      }
      
      tiles_ctb_addr_rs_to_ts[j] = 0;
      for (i = 0; i < tileX; ++i) {
        tiles_ctb_addr_rs_to_ts[j] += tiles_row_height[tileY] * tiles_col_width[i];
      }
      for (i = 0; i < tileY; ++i) {
        tiles_ctb_addr_rs_to_ts[j] += encoder->in.width_in_lcu * tiles_row_height[i];
      }
      tiles_ctb_addr_rs_to_ts[j] += (tbY - tiles_row_bd[tileY]) * tiles_col_width[tileX] + 
                                     tbX - tiles_col_bd[tileX];
    }
    
    //(6-8) in ITU-T Rec. H.265 (04/2013)
    //Make reverse map from tile scan to raster scan
    for (j = 0; j < num_ctbs; ++j) {
      tiles_ctb_addr_ts_to_rs[tiles_ctb_addr_rs_to_ts[j]] = j;
    }
    
    //(6-9) in ITU-T Rec. H.265 (04/2013)
    tileIdx = 0;
    for (j=0; j < encoder->tiles_num_tile_rows; ++j) {
      for (i=0; i < encoder->tiles_num_tile_columns; ++i) {
        for (y = tiles_row_bd[j]; y < tiles_row_bd[j+1]; ++y) {
          for (x = tiles_col_bd[i]; x < tiles_col_bd[i+1]; ++x) {
            tiles_tile_id[tiles_ctb_addr_rs_to_ts[y * encoder->in.width_in_lcu + x]] = tileIdx;
          }
        }
        ++tileIdx;
      }
    }

    encoder->tiles_col_width = tiles_col_width;
    encoder->tiles_row_height = tiles_row_height;
    
    encoder->tiles_row_bd = tiles_row_bd;
    encoder->tiles_col_bd = tiles_col_bd;
    
    encoder->tiles_ctb_addr_rs_to_ts = tiles_ctb_addr_rs_to_ts;
    encoder->tiles_ctb_addr_ts_to_rs = tiles_ctb_addr_ts_to_rs;
    
    encoder->tiles_tile_id = tiles_tile_id;
    
    //Slices
    {
      int *slice_addresses_in_ts;
      encoder->slice_count = encoder->cfg->slice_count;
      if (encoder->slice_count == 0) {
        encoder->slice_count = 1;
        slice_addresses_in_ts = MALLOC(int, encoder->slice_count);
        slice_addresses_in_ts[0] = 0;
      } else {
        int i;
        slice_addresses_in_ts = MALLOC(int, encoder->slice_count);
        if (!encoder->cfg->slice_addresses_in_ts) {
          slice_addresses_in_ts[0] = 0;
          for (i=1; i < encoder->slice_count; ++i) {
            slice_addresses_in_ts[i] = encoder->in.width_in_lcu * encoder->in.height_in_lcu * i / encoder->slice_count;
          }
        } else {
          for (i=0; i < encoder->slice_count; ++i) {
            slice_addresses_in_ts[i] = encoder->cfg->slice_addresses_in_ts[i];
          }
        }
      }
      
      encoder->slice_addresses_in_ts = slice_addresses_in_ts;
    }
    
    encoder->wpp = encoder->cfg->wpp;

#ifdef _DEBUG
    printf("Tiles columns width:");
    for (i=0; i < encoder->tiles_num_tile_columns; ++i) {
      printf(" %d", encoder->tiles_col_width[i]);
    }
    printf("\n");
    printf("Tiles row height:");
    for (i=0; i < encoder->tiles_num_tile_rows; ++i) {
      printf(" %d", encoder->tiles_row_height[i]);
    }
    printf("\n");
    //Print tile index map
    for (y = 0; y < encoder->in.height_in_lcu; ++y) {
      for (x = 0; x < encoder->in.width_in_lcu; ++x) {
        const int lcu_id_rs = y * encoder->in.width_in_lcu + x;
        const int lcu_id_ts = encoder->tiles_ctb_addr_rs_to_ts[lcu_id_rs];
        const char slice_start = lcu_at_slice_start(encoder, lcu_id_ts) ? '|' : ' ';
        const char slice_end = lcu_at_slice_end(encoder, lcu_id_ts)  ? '|' : ' ';
        
        printf("%c%03d%c", slice_start, encoder->tiles_tile_id[lcu_id_ts], slice_end);
      }
      printf("\n");
    }
    printf("\n");
    if (encoder->wpp) {
      printf("Wavefront Parallel Processing: enabled\n");
    } else {
      printf("Wavefront Parallel Processing: disabled\n");
    }
    printf("\n");
#endif //_DEBUG
  }

  assert(WITHIN(cfg->pu_depth_inter.min, PU_DEPTH_INTER_MIN, PU_DEPTH_INTER_MAX));
  assert(WITHIN(cfg->pu_depth_inter.max, PU_DEPTH_INTER_MIN, PU_DEPTH_INTER_MAX));
  assert(WITHIN(cfg->pu_depth_intra.min, PU_DEPTH_INTRA_MIN, PU_DEPTH_INTRA_MAX));
  assert(WITHIN(cfg->pu_depth_intra.max, PU_DEPTH_INTRA_MIN, PU_DEPTH_INTRA_MAX));
  encoder->pu_depth_inter.min = cfg->pu_depth_inter.min;
  encoder->pu_depth_inter.max = cfg->pu_depth_inter.max;
  encoder->pu_depth_intra.min = cfg->pu_depth_intra.min;
  encoder->pu_depth_intra.max = cfg->pu_depth_intra.max;

  // input init (TODO: read from commandline / config)
  encoder->bitdepth = 8;
  encoder->in.video_format = FORMAT_420;

  // deblocking filter
  encoder->deblock_enable = (int8_t)encoder->cfg->deblock_enable;
  encoder->beta_offset_div2 = (int8_t)encoder->cfg->deblock_beta;
  encoder->tc_offset_div2 = (int8_t)encoder->cfg->deblock_tc;
  // SAO
  encoder->sao_enable = (int8_t)encoder->cfg->sao_enable;
  // RDO
  encoder->rdoq_enable = (int8_t)encoder->cfg->rdoq_enable;
  encoder->rdo = (int8_t)encoder->cfg->rdo;
  encoder->sign_hiding = encoder->cfg->signhide_enable;
  encoder->full_intra_search = (int8_t)encoder->cfg->full_intra_search;
  // TR SKIP
  encoder->trskip_enable = (int8_t)encoder->cfg->trskip_enable;
  encoder->tr_depth_intra = (int8_t)encoder->cfg->tr_depth_intra;
  // MOTION ESTIMATION
  encoder->fme_level = (int8_t)encoder->cfg->fme_level;
  // VUI
  encoder->vui.sar_width = (int16_t)encoder->cfg->vui.sar_width;
  encoder->vui.sar_height = (int16_t)encoder->cfg->vui.sar_height;
  encoder->vui.overscan = encoder->cfg->vui.overscan;
  encoder->vui.videoformat = encoder->cfg->vui.videoformat;
  encoder->vui.fullrange = encoder->cfg->vui.fullrange;
  encoder->vui.colorprim = encoder->cfg->vui.colorprim;
  encoder->vui.transfer = encoder->cfg->vui.transfer;
  encoder->vui.colormatrix = encoder->cfg->vui.colormatrix;
  encoder->vui.chroma_loc = (int8_t)encoder->cfg->vui.chroma_loc;
  // AUD
  encoder->aud_enable = (int8_t)encoder->cfg->aud_enable;

  encoder->vps_period = encoder->cfg->vps_period * encoder->cfg->intra_period;

  return 1;
}

int encoder_control_finalize(encoder_control_t * const encoder) {
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
  scalinglist_destroy(&encoder->scaling_list);
  
  if (!threadqueue_finalize(encoder->threadqueue)) {
    fprintf(stderr, "Could not initialize threadqueue");
    return 0;
  }
  
  FREE_POINTER(encoder->threadqueue);

  
  return 1;
}

void encoder_control_input_init(encoder_control_t * const encoder,
                        const int32_t width, const int32_t height)
{
  encoder->in.width = width;
  encoder->in.height = height;
  encoder->in.real_width = width;
  encoder->in.real_height = height;
  encoder->in.bitdepth = encoder->bitdepth;

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


  #ifdef _DEBUG
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

  gop_config_t const * const gop = encoder->cfg->gop;
  const int8_t gop_len = encoder->cfg->gop_len;

  int num_layers = 0;
  for (int i = 0; i < gop_len; ++i) {
    num_layers = MAX(gop[i].layer, num_layers);
  }

  switch (num_layers) {
    case 0:
      break;

    case 4:
      // These weights were copied from http://doi.org/10.1109/TIP.2014.2336550
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
