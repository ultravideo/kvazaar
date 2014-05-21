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
                       unsigned x_lcu, uint16_t y_lcu,
                       sao_info *sao_luma, sao_info *sao_chroma);

static void encoder_state_write_bitstream_leaf(encoder_state * const encoder_state);
static void worker_encoder_state_write_bitstream_leaf(void * opaque);

/*!
  \brief Initializes lambda-value for current QP

  Implementation closer to HM (Used HM12 as reference)
   - Still missing functionality when GOP and B-pictures are used
 */
void encoder_state_init_lambda(encoder_state * const encoder_state)
{
  double qp = encoder_state->global->QP;
  double lambda_scale = 1.0;
  double qp_temp      = qp - 12;
  double lambda;

  // Default QP-factor from HM config
  double qp_factor = 0.4624;

  if (encoder_state->global->slicetype == SLICE_I) {
    qp_factor=0.57*lambda_scale;
  }

  lambda = qp_factor*pow( 2.0, qp_temp/3.0 );

  if (encoder_state->global->slicetype != SLICE_I ) {
    lambda *= 0.95;
  }

  encoder_state->global->cur_lambda_cost = lambda;
}

static int lcu_at_slice_start(const encoder_control * const encoder, int lcu_addr_in_ts) {
  int i;
  assert(lcu_addr_in_ts >= 0 && lcu_addr_in_ts < encoder->in.height_in_lcu * encoder->in.width_in_lcu);
  if (lcu_addr_in_ts == 0) return 1;
  for (i = 0; i < encoder->slice_count; ++i) {
    if (encoder->slice_addresses_in_ts[i] == lcu_addr_in_ts) return 1;
  }
  return 0;
}

static int lcu_at_slice_end(const encoder_control * const encoder, int lcu_addr_in_ts) {
  int i;
  assert(lcu_addr_in_ts >= 0 && lcu_addr_in_ts < encoder->in.height_in_lcu * encoder->in.width_in_lcu);
  if (lcu_addr_in_ts == encoder->in.height_in_lcu * encoder->in.width_in_lcu - 1) return 1;
  for (i = 0; i < encoder->slice_count; ++i) {
    if (encoder->slice_addresses_in_ts[i] == lcu_addr_in_ts + 1) return 1;
  }
  return 0;
}

static int lcu_at_tile_start(const encoder_control * const encoder, int lcu_addr_in_ts) {
  assert(lcu_addr_in_ts >= 0 && lcu_addr_in_ts < encoder->in.height_in_lcu * encoder->in.width_in_lcu);
  if (lcu_addr_in_ts == 0) return 1;
  if (encoder->tiles_tile_id[lcu_addr_in_ts - 1] != encoder->tiles_tile_id[lcu_addr_in_ts]) {
    return 1;
  }
  return 0;
}

static int lcu_at_tile_end(const encoder_control * const encoder, int lcu_addr_in_ts) {
  assert(lcu_addr_in_ts >= 0 && lcu_addr_in_ts < encoder->in.height_in_lcu * encoder->in.width_in_lcu);
  if (lcu_addr_in_ts == encoder->in.height_in_lcu * encoder->in.width_in_lcu - 1) return 1;
  if (encoder->tiles_tile_id[lcu_addr_in_ts + 1] != encoder->tiles_tile_id[lcu_addr_in_ts]) {
    return 1;
  }
  return 0;
}

//Return 1 if the LCU is at the first row of a structure (tile or slice)
static int lcu_in_first_row(const encoder_state * const encoder_state, int lcu_addr_in_ts) {
  const int lcu_addr_in_rs = encoder_state->encoder_control->tiles_ctb_addr_ts_to_rs[lcu_addr_in_ts];
  
  if (lcu_addr_in_rs / encoder_state->encoder_control->in.width_in_lcu == encoder_state->tile->lcu_offset_y) {
    return 1;
  }
  
  if (lcu_addr_in_rs / encoder_state->encoder_control->in.width_in_lcu == encoder_state->slice->start_in_rs / encoder_state->encoder_control->in.width_in_lcu) {
    return 1;
  }
  
  //One row above is before the start of the slice => it's also a boundary
  if (lcu_addr_in_rs - encoder_state->encoder_control->in.width_in_lcu < encoder_state->slice->start_in_rs) {
    return 1;
  }
  
  return 0;
}

//Return 1 if the LCU is at the first row of a structure (tile or slice)
static int lcu_in_last_row(const encoder_state * const encoder_state, int lcu_addr_in_ts) {
  const int lcu_addr_in_rs = encoder_state->encoder_control->tiles_ctb_addr_ts_to_rs[lcu_addr_in_ts];
  
  if (lcu_addr_in_rs / encoder_state->encoder_control->in.width_in_lcu == encoder_state->tile->lcu_offset_y + encoder_state->tile->cur_pic->height_in_lcu - 1) {
    return 1;
  }
  
  if (lcu_addr_in_rs / encoder_state->encoder_control->in.width_in_lcu == encoder_state->slice->end_in_rs / encoder_state->encoder_control->in.width_in_lcu) {
    return 1;
  }
  
  //One row below is before the end of the slice => it's also a boundary
  if (lcu_addr_in_rs + encoder_state->encoder_control->in.width_in_lcu > encoder_state->slice->end_in_rs) {
    return 1;
  }
  
  return 0;
}


//Return 1 if the LCU is at the first column of a structure (tile or slice)
static int lcu_in_first_column(const encoder_state * const encoder_state, int lcu_addr_in_ts) {
  const int lcu_addr_in_rs = encoder_state->encoder_control->tiles_ctb_addr_ts_to_rs[lcu_addr_in_ts];
  
  //First column of tile?
  if (lcu_addr_in_rs % encoder_state->encoder_control->in.width_in_lcu == encoder_state->tile->lcu_offset_x) {
    return 1;
  }
  
  //Slice start may not be aligned with the tile, so we need to allow this
  if (lcu_addr_in_rs == encoder_state->slice->start_in_rs) {
    return 1;
  }
  
  return 0;
}

//Return 1 if the LCU is at the last column of a structure (tile or slice)
static int lcu_in_last_column(const encoder_state * const encoder_state, int lcu_addr_in_ts) {
  const int lcu_addr_in_rs = encoder_state->encoder_control->tiles_ctb_addr_ts_to_rs[lcu_addr_in_ts];
  
  //First column of tile?
  if (lcu_addr_in_rs % encoder_state->encoder_control->in.width_in_lcu == encoder_state->tile->lcu_offset_x + encoder_state->tile->cur_pic->width_in_lcu - 1) {
    return 1;
  }
  
  //Slice start may not be aligned with the tile, so we need to allow this
  if (lcu_addr_in_rs == encoder_state->slice->end_in_rs) {
    return 1;
  }
  
  return 0;
}

int encoder_control_init(encoder_control * const encoder, const config * const cfg) {
  if (!cfg) {
    fprintf(stderr, "Config object must not be null!\n");
    return 0;
  }
  
  encoder->threadqueue = MALLOC(threadqueue_queue, 1);
    
  //Init threadqueue
  if (!encoder->threadqueue || !threadqueue_init(encoder->threadqueue, cfg->threads)) {
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
  
  return 1;
}

int encoder_control_finalize(encoder_control * const encoder) {
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

void encoder_control_input_init(encoder_control * const encoder,
                        const int32_t width, const int32_t height)
{
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
  if (width != encoder->in.width || height != encoder->in.height) {
    printf("Picture buffer has been extended to be a multiple of the smallest block size:\r\n");
    printf("  Width = %d (%d), Height = %d (%d)\r\n", width, encoder->in.width, height,
           encoder->in.height);
  }
  #endif
}

static int encoder_state_config_global_init(encoder_state * const encoder_state) {
  encoder_state->global->ref = picture_list_init(MAX_REF_PIC_COUNT);
  if(!encoder_state->global->ref) {
    fprintf(stderr, "Failed to allocate the picture list!\n");
    return 0;
  }
  encoder_state->global->ref_list = REF_PIC_LIST_0;
  encoder_state->global->frame = 0;
  encoder_state->global->poc = 0;
  return 1;
}

static void encoder_state_config_global_finalize(encoder_state * const encoder_state) {
  picture_list_destroy(encoder_state->global->ref);
}



static int encoder_state_config_tile_init(encoder_state * const encoder_state, 
                                          const int lcu_offset_x, const int lcu_offset_y,
                                          const int width, const int height, const int width_in_lcu, const int height_in_lcu) {
  
  const encoder_control * const encoder = encoder_state->encoder_control;
  encoder_state->tile->cur_pic = picture_alloc(width, height, width_in_lcu, height_in_lcu);

  if (!encoder_state->tile->cur_pic) {
    printf("Error allocating picture!\r\n");
    return 0;
  }
  
  // Init coeff data table
  //FIXME: move them
  encoder_state->tile->cur_pic->coeff_y = MALLOC(coefficient, width * height);
  encoder_state->tile->cur_pic->coeff_u = MALLOC(coefficient, (width * height) >> 2);
  encoder_state->tile->cur_pic->coeff_v = MALLOC(coefficient, (width * height) >> 2);
  
  encoder_state->tile->lcu_offset_x = lcu_offset_x;
  encoder_state->tile->lcu_offset_y = lcu_offset_y;
  
  encoder_state->tile->lcu_offset_in_ts = encoder->tiles_ctb_addr_rs_to_ts[lcu_offset_x + lcu_offset_y * encoder->in.width_in_lcu];
  
  //Allocate buffers
  //order by row of (LCU_WIDTH * cur_pic->width_in_lcu) pixels
  encoder_state->tile->hor_buf_search = yuv_t_alloc(LCU_WIDTH * encoder_state->tile->cur_pic->width_in_lcu * encoder_state->tile->cur_pic->height_in_lcu);
  //order by column of (LCU_WIDTH * encoder_state->height_in_lcu) pixels (there is no more extra pixel, since we can use a negative index)
  encoder_state->tile->ver_buf_search = yuv_t_alloc(LCU_WIDTH * encoder_state->tile->cur_pic->height_in_lcu * encoder_state->tile->cur_pic->width_in_lcu);
  
  if (encoder->sao_enable) {
    encoder_state->tile->hor_buf_before_sao = yuv_t_alloc(LCU_WIDTH * encoder_state->tile->cur_pic->width_in_lcu * encoder_state->tile->cur_pic->height_in_lcu);
  } else {
    encoder_state->tile->hor_buf_before_sao = NULL;
  }
  
  if (encoder->wpp) {
    encoder_state->tile->wf_jobs = MALLOC(threadqueue_job*, encoder_state->tile->cur_pic->width_in_lcu * encoder_state->tile->cur_pic->height_in_lcu);
    if (!encoder_state->tile->wf_jobs) {
      printf("Error allocating wf_jobs array!\n");
      return 0;
    }
  } else {
    encoder_state->tile->wf_jobs = NULL;
  }
  
  encoder_state->tile->id = encoder->tiles_tile_id[encoder_state->tile->lcu_offset_in_ts];
  return 1;
}

static void encoder_state_config_tile_finalize(encoder_state * const encoder_state) {
  if (encoder_state->tile->hor_buf_before_sao) yuv_t_free(encoder_state->tile->hor_buf_before_sao);
  
  yuv_t_free(encoder_state->tile->hor_buf_search);
  yuv_t_free(encoder_state->tile->ver_buf_search);
  
  picture_free(encoder_state->tile->cur_pic);
  encoder_state->tile->cur_pic = NULL;
  
  FREE_POINTER(encoder_state->tile->wf_jobs);
}

static int encoder_state_config_slice_init(encoder_state * const encoder_state, 
                                          const int start_address_in_ts, const int end_address_in_ts) {
  int i = 0, slice_found=0;
  for (i = 0; i < encoder_state->encoder_control->slice_count; ++i) {
    if (encoder_state->encoder_control->slice_addresses_in_ts[i] == start_address_in_ts) {
      encoder_state->slice->id = i;
      slice_found = 1;
      break;
    }
  }
  assert(slice_found);
  encoder_state->slice->start_in_ts = start_address_in_ts;
  encoder_state->slice->end_in_ts = end_address_in_ts;
  
  encoder_state->slice->start_in_rs = encoder_state->encoder_control->tiles_ctb_addr_ts_to_rs[start_address_in_ts];
  encoder_state->slice->end_in_rs = encoder_state->encoder_control->tiles_ctb_addr_ts_to_rs[end_address_in_ts];
  return 1;
}

static void encoder_state_config_slice_finalize(encoder_state * const encoder_state) {
  //Nothing to do (yet?)
}

static int encoder_state_config_wfrow_init(encoder_state * const encoder_state, 
                                          const int lcu_offset_y) {
  
  encoder_state->wfrow->lcu_offset_y = lcu_offset_y;
  return 1;
}

static void encoder_state_config_wfrow_finalize(encoder_state * const encoder_state) {
  //Nothing to do (yet?)
}

#ifdef _DEBUG
static void encoder_state_dump_graphviz(const encoder_state * const encoder_state) {
  int i;
  
  if (!encoder_state->parent) {
    const encoder_control * const encoder = encoder_state->encoder_control;
    int y,x;
    //Empty lines (easier to copy-paste)
    printf("\n\n\n\n\n");
    //Some styling...
    printf("digraph EncoderStates {\n");
    printf(" fontname = \"Bitstream Vera Sans\"\n");
    printf(" fontsize = 8\n\n");
    printf(" node [\n");
    printf("  fontname = \"Bitstream Vera Sans\"\n");
    printf("  fontsize = 8\n");
    printf("  shape = \"record\"\n");
    printf(" ]\n\n");
    printf(" edge [\n");
    printf("  arrowtail = \"empty\"\n");
    printf(" ]\n\n");
    
    printf(" \"Map\" [\n");
    printf("  shape=plaintext\n");
    printf("  label = <<table cellborder=\"1\" cellspacing=\"0\" border=\"0\">");
    printf("<tr><td colspan=\"%d\" height=\"20\" valign=\"bottom\"><b>RS Map</b></td></tr>", encoder->in.width_in_lcu);
    for (y = 0; y < encoder->in.height_in_lcu; ++y) {
      printf("<tr>");
      for (x = 0; x < encoder->in.width_in_lcu; ++x) {
        const int lcu_id_rs = y * encoder->in.width_in_lcu + x;
        
        printf("<td>%d</td>", lcu_id_rs);
      }
      printf("</tr>");
    }
    printf("<tr><td colspan=\"%d\" height=\"20\" valign=\"bottom\"><b>TS Map</b></td></tr>", encoder->in.width_in_lcu);
    for (y = 0; y < encoder->in.height_in_lcu; ++y) {
      printf("<tr>");
      for (x = 0; x < encoder->in.width_in_lcu; ++x) {
        const int lcu_id_rs = y * encoder->in.width_in_lcu + x;
        const int lcu_id_ts = encoder->tiles_ctb_addr_rs_to_ts[lcu_id_rs];
        
        printf("<td>%d</td>", lcu_id_ts);
      }
      printf("</tr>");
    }
    printf("<tr><td colspan=\"%d\" height=\"20\" valign=\"bottom\"><b>Tile map</b></td></tr>", encoder->in.width_in_lcu);
    for (y = 0; y < encoder->in.height_in_lcu; ++y) {
      printf("<tr>");
      for (x = 0; x < encoder->in.width_in_lcu; ++x) {
        const int lcu_id_rs = y * encoder->in.width_in_lcu + x;
        const int lcu_id_ts = encoder->tiles_ctb_addr_rs_to_ts[lcu_id_rs];
        
        printf("<td>%d</td>", encoder->tiles_tile_id[lcu_id_ts]);
      }
      printf("</tr>");
    }
    printf("<tr><td colspan=\"%d\" height=\"20\" valign=\"bottom\"><b>Slice map</b></td></tr>", encoder->in.width_in_lcu);
    for (y = 0; y < encoder->in.height_in_lcu; ++y) {
      printf("<tr>");
      for (x = 0; x < encoder->in.width_in_lcu; ++x) {
        const int lcu_id_rs = y * encoder->in.width_in_lcu + x;
        const int lcu_id_ts = encoder->tiles_ctb_addr_rs_to_ts[lcu_id_rs];
        int slice_id = 0;
        
        //Not efficient, but who cares
        for (i=0; i < encoder->slice_count; ++i) {
          if (encoder->slice_addresses_in_ts[i] <= lcu_id_ts) {
            slice_id = i;
          }
        }
        
        printf("<td>%d</td>", slice_id);
      }
      printf("</tr>");
    }
    printf("</table>>\n ]\n");
  }
  
  printf(" \"%p\" [\n", encoder_state);
  printf("  label = \"{encoder_state|");
  printf("+ type=%c\\l", encoder_state->type);
  if (!encoder_state->parent || encoder_state->global != encoder_state->parent->global) {
    printf("|+ global\\l");
  }
  if (!encoder_state->parent || encoder_state->tile != encoder_state->parent->tile) {
    printf("|+ tile\\l");
    printf(" - id = %d\\l", encoder_state->tile->id);
    printf(" - lcu_offset_x = %d\\l", encoder_state->tile->lcu_offset_x);
    printf(" - lcu_offset_y = %d\\l", encoder_state->tile->lcu_offset_y);
    printf(" - lcu_offset_in_ts = %d\\l", encoder_state->tile->lcu_offset_in_ts);
  }
  if (!encoder_state->parent || encoder_state->slice != encoder_state->parent->slice) {
    printf("|+ slice\\l");
    printf(" - id = %d\\l", encoder_state->slice->id);
    printf(" - start_in_ts = %d\\l", encoder_state->slice->start_in_ts);
    printf(" - end_in_ts = %d\\l", encoder_state->slice->end_in_ts);
    printf(" - start_in_rs = %d\\l", encoder_state->slice->start_in_rs);
    printf(" - end_in_rs = %d\\l", encoder_state->slice->end_in_rs);
  }
  if (!encoder_state->parent || encoder_state->wfrow != encoder_state->parent->wfrow) {
    printf("|+ wfrow\\l");
    printf(" - lcu_offset_y = %d\\l", encoder_state->wfrow->lcu_offset_y);
  }
  printf("}\"\n");
  printf(" ]\n");
  
  if (encoder_state->parent) {
    printf(" \"%p\" -> \"%p\"\n", encoder_state->parent, encoder_state);
  }
  
  for (i = 0; encoder_state->children[i].encoder_control; ++i) {
    encoder_state_dump_graphviz(&encoder_state->children[i]);
  }
  
  if (!encoder_state->parent) {
    printf("}\n");
    //Empty lines (easier to copy-paste)
    printf("\n\n\n\n\n");
  }
}
#endif //_DEBUG

int encoder_state_init(encoder_state * const child_state, encoder_state * const parent_state) {
  //We require that, if parent_state is NULL:
  //child_state->encoder_control is set
  //
  //If parent_state is not NULL, the following variable should either be set to NULL,
  //in order to inherit from parent, or should point to a valid structure:
  //child_state->global
  //child_state->tile
  //child_state->slice
  //child_state->wfrow
  
  child_state->parent = parent_state;
  child_state->children = MALLOC(encoder_state, 1);
  child_state->children[0].encoder_control = NULL;
  
  if (!parent_state) {
    const encoder_control * const encoder = child_state->encoder_control;
    child_state->type = ENCODER_STATE_TYPE_MAIN;
    assert(child_state->encoder_control);
    child_state->global = MALLOC(encoder_state_config_global, 1);
    if (!child_state->global || !encoder_state_config_global_init(child_state)) {
      fprintf(stderr, "Could not initialize encoder_state->global!\n");
      return 0;
    }
    child_state->tile = MALLOC(encoder_state_config_tile, 1);
    if (!child_state->tile || !encoder_state_config_tile_init(child_state, 0, 0, encoder->in.width, encoder->in.height, encoder->in.width_in_lcu, encoder->in.height_in_lcu)) {
      fprintf(stderr, "Could not initialize encoder_state->tile!\n");
      return 0;
    }
    child_state->slice = MALLOC(encoder_state_config_slice, 1);
    if (!child_state->slice || !encoder_state_config_slice_init(child_state, 0, encoder->in.width_in_lcu * encoder->in.height_in_lcu - 1)) {
      fprintf(stderr, "Could not initialize encoder_state->slice!\n");
      return 0;
    }
    child_state->wfrow = MALLOC(encoder_state_config_wfrow, 1);
    if (!child_state->wfrow || !encoder_state_config_wfrow_init(child_state, 0)) {
      fprintf(stderr, "Could not initialize encoder_state->wfrow!\n");
      return 0;
    }
  } else {
    child_state->encoder_control = parent_state->encoder_control;
    if (!child_state->global) child_state->global = parent_state->global;
    if (!child_state->tile) child_state->tile = parent_state->tile;
    if (!child_state->slice) child_state->slice = parent_state->slice;
    if (!child_state->wfrow) child_state->wfrow = parent_state->wfrow;
  }
  
  //Allocate bitstream
  if (child_state->type == ENCODER_STATE_TYPE_MAIN) {
    //Main encoder outputs to file
    if (!bitstream_init(&child_state->stream, BITSTREAM_TYPE_FILE)) {
      fprintf(stderr, "Could not initialize stream!\n");
      return 0;
    }
    child_state->stream.file.output = child_state->encoder_control->out.file;
  } else {
    //Other encoders use a memory bitstream
    if (!bitstream_init(&child_state->stream, BITSTREAM_TYPE_MEMORY)) {
      fprintf(stderr, "Could not initialize stream!\n");
      return 0;
    }
  }
  
  // Set CABAC output bitstream
  child_state->cabac.stream = &child_state->stream;
  
  //Create sub-encoders
  {
    const encoder_control * const encoder = child_state->encoder_control;
    int child_count = 0;
    //We first check the type of this element.
    //If it's a MAIN, it can allow both slices or tiles as child
    //If it's a TILE, it can allow slices as child, if its parent is not a slice, or wavefront rows if there is no other children
    //If it's a SLICE, it can allow tiles as child, if its parent is not a tile, or wavefront rows if there is no other children
    //If it's a WAVEFRONT_ROW, it doesn't allow any children
    int children_allow_wavefront_row = 0;
    int children_allow_slice = 0;
    int children_allow_tile = 0;
    int range_start;
    
    int start_in_ts, end_in_ts;
    
    switch(child_state->type) {
      case ENCODER_STATE_TYPE_MAIN:
        children_allow_slice = 1;
        children_allow_tile = 1;
        start_in_ts = 0;
        end_in_ts = child_state->tile->cur_pic->width_in_lcu * child_state->tile->cur_pic->height_in_lcu;
        break;
      case ENCODER_STATE_TYPE_SLICE:
        assert(child_state->parent);
        if (child_state->parent->type != ENCODER_STATE_TYPE_TILE) children_allow_tile = 1;
        children_allow_wavefront_row = encoder->wpp;
        start_in_ts = child_state->slice->start_in_ts;
        end_in_ts = child_state->slice->end_in_ts;
        break;
      case ENCODER_STATE_TYPE_TILE:
        assert(child_state->parent);
        if (child_state->parent->type != ENCODER_STATE_TYPE_SLICE) children_allow_slice = 1;
        children_allow_wavefront_row = encoder->wpp;
        start_in_ts = child_state->tile->lcu_offset_in_ts;
        end_in_ts = child_state->tile->lcu_offset_in_ts + child_state->tile->cur_pic->width_in_lcu * child_state->tile->cur_pic->height_in_lcu;
        break;
      case ENCODER_STATE_TYPE_WAVEFRONT_ROW:
        //GCC tries to be too clever...
        start_in_ts = -1;
        end_in_ts = -1;
        break;
      default:
        fprintf(stderr, "Invalid encoder_state->type %d!\n", child_state->type);
        assert(0);
        return 0;
    }
    
    range_start = start_in_ts;
    //printf("%c-%p: start_in_ts=%d, end_in_ts=%d\n",child_state->type, child_state, start_in_ts, end_in_ts);
    while (range_start < end_in_ts && (children_allow_slice || children_allow_tile)) {
      encoder_state *new_child = NULL;
      int range_end_slice = range_start; //Will be incremented to get the range of the "thing"
      int range_end_tile = range_start; //Will be incremented to get the range of the "thing"
      
      int tile_allowed = lcu_at_tile_start(encoder, range_start) && children_allow_tile;
      int slice_allowed = lcu_at_slice_start(encoder, range_start) && children_allow_slice;
      
      //Find the smallest structure following the cursor
      if (slice_allowed) {
        while(!lcu_at_slice_end(encoder, range_end_slice)) {
          ++range_end_slice;
        }
      }
      
      if (tile_allowed) {
        while(!lcu_at_tile_end(encoder, range_end_tile)) {
          ++range_end_tile;
        }
      }
      
      //printf("range_start=%d, range_end_slice=%d, range_end_tile=%d, tile_allowed=%d, slice_allowed=%d end_in_ts=%d\n",range_start,range_end_slice,range_end_tile,tile_allowed,slice_allowed,end_in_ts);
      
      if ((!tile_allowed || (range_end_slice >= range_end_tile)) && !new_child && slice_allowed) {
        //Create a slice
        new_child = &child_state->children[child_count];
        new_child->encoder_control = encoder;
        new_child->type = ENCODER_STATE_TYPE_SLICE;
        new_child->global = child_state->global;
        new_child->tile = child_state->tile;
        new_child->wfrow = child_state->wfrow;
        new_child->slice = MALLOC(encoder_state_config_slice, 1);
        if (!new_child->slice || !encoder_state_config_slice_init(new_child, range_start, range_end_slice)) {
          fprintf(stderr, "Could not initialize encoder_state->slice!\n");
          return 0;
        }
      }
      
      if ((!slice_allowed || (range_end_slice < range_end_tile)) && !new_child && tile_allowed) {
        //Create a tile
        int tile_id = encoder->tiles_tile_id[range_start];
        int tile_x = tile_id % encoder->tiles_num_tile_columns;
        int tile_y = tile_id / encoder->tiles_num_tile_columns;
        
        int lcu_offset_x = encoder->tiles_col_bd[tile_x];
        int lcu_offset_y = encoder->tiles_row_bd[tile_y];
        int width_in_lcu = encoder->tiles_col_bd[tile_x+1]-encoder->tiles_col_bd[tile_x];
        int height_in_lcu = encoder->tiles_row_bd[tile_y+1]-encoder->tiles_row_bd[tile_y];
        int width = MIN(width_in_lcu * LCU_WIDTH, encoder->in.width - lcu_offset_x * LCU_WIDTH);
        int height = MIN(height_in_lcu * LCU_WIDTH, encoder->in.height - lcu_offset_y * LCU_WIDTH);
        
        new_child = &child_state->children[child_count];
        new_child->encoder_control = encoder;
        new_child->type = ENCODER_STATE_TYPE_TILE;
        new_child->global = child_state->global;
        new_child->tile = MALLOC(encoder_state_config_tile, 1);
        new_child->slice = child_state->slice;
        new_child->wfrow = child_state->wfrow;
        
        if (!new_child->tile || !encoder_state_config_tile_init(new_child, lcu_offset_x, lcu_offset_y, width, height, width_in_lcu, height_in_lcu)) {
          fprintf(stderr, "Could not initialize encoder_state->tile!\n");
          return 0;
        }
      }
      
      if (new_child) {
        child_state->children = realloc(child_state->children, sizeof(encoder_state) * (2+child_count));
        child_state->children[1+child_count].encoder_control = NULL;
        if (!child_state->children) {
          fprintf(stderr, "Failed to allocate memory for children...\n");
          return 0;
        }

        //Fix children parent (since we changed the address), except for the last one which is not ready yet
        {
          int i, j;
          for (i = 0; child_state->children[i].encoder_control && i < child_count; ++i) {
            for (j = 0; child_state->children[i].children[j].encoder_control; ++j) {
              child_state->children[i].children[j].parent = &child_state->children[i];
            }
            for (j = 0; j < child_state->children[i].lcu_order_count; ++j) {
              child_state->children[i].lcu_order[j].encoder_state = &child_state->children[i];
            }
            child_state->children[i].cabac.stream = &child_state->children[i].stream;
          }
        }
        
        if (!encoder_state_init(&child_state->children[child_count], child_state)) {
          fprintf(stderr, "Unable to init child...\n");
          return 0;
        }
        child_count += 1;
      }
      
      range_start = MAX(range_end_slice, range_end_tile) + 1;
    }
    
    //We create wavefronts only if we have no children
    if (children_allow_wavefront_row && child_count == 0) {
      int first_row = encoder->tiles_ctb_addr_ts_to_rs[start_in_ts] / encoder->in.width_in_lcu;
      int last_row = encoder->tiles_ctb_addr_ts_to_rs[start_in_ts] / encoder->in.width_in_lcu;
      int num_rows;
      int i;
      
      assert(!(children_allow_slice || children_allow_tile));
      assert(child_count == 0);
      
      for (i=start_in_ts; i<end_in_ts; ++i) {
        const int row = encoder->tiles_ctb_addr_ts_to_rs[i] / encoder->in.width_in_lcu;
        if (row < first_row) first_row = row;
        if (row > last_row) last_row = row;
      }
      
      num_rows = last_row - first_row + 1;
      
      //When entropy_coding_sync_enabled_flag is equal to 1 and the first coding tree block in a slice is not the first coding
      //tree block of a row of coding tree blocks in a tile, it is a requirement of bitstream conformance that the last coding tree
      //block in the slice shall belong to the same row of coding tree blocks as the first coding tree block in the slice.
      
      if (encoder->tiles_ctb_addr_ts_to_rs[start_in_ts] % encoder->in.width_in_lcu != child_state->tile->lcu_offset_x) {
        if (num_rows > 1) {
          fprintf(stderr, "Invalid: first CTB in slice %d is not at the tile %d edge, and the slice spans on more than one row.\n", child_state->slice->id, child_state->tile->id);
          return 0;
        }
      }
      
      //FIXME Do the same kind of check if we implement slice segments
    
      child_count = num_rows;
      child_state->children = realloc(child_state->children, sizeof(encoder_state) * (num_rows + 1));
      child_state->children[num_rows].encoder_control = NULL;
      
      for (i=0; i < num_rows; ++i) {
        encoder_state *new_child = &child_state->children[i];
        
        new_child->encoder_control = encoder;
        new_child->type = ENCODER_STATE_TYPE_WAVEFRONT_ROW;
        new_child->global = child_state->global;
        new_child->tile = child_state->tile;
        new_child->slice = child_state->slice;
        new_child->wfrow = MALLOC(encoder_state_config_wfrow, 1);
        
        if (!new_child->wfrow || !encoder_state_config_wfrow_init(new_child, i)) {
          fprintf(stderr, "Could not initialize encoder_state->wfrow!\n");
          return 0;
        }
        
        if (!encoder_state_init(new_child, child_state)) {
          fprintf(stderr, "Unable to init child...\n");
          return 0;
        }
      }
    }
    
    child_state->is_leaf = (child_count == 0);
    //This node is a leaf, compute LCU-order
    if (child_state->is_leaf) {
      //All LCU computations are relative to the tile
      //Remark: this could be optimized, but since it's run only once, it's better to do it in a understandable way.
      
      //By default, the full tile
      int i;
      int lcu_id;
      int lcu_start = 0;
      //End is the element AFTER the end (iterate < lcu_end)
      int lcu_end = child_state->tile->cur_pic->width_in_lcu * child_state->tile->cur_pic->height_in_lcu;
      
      //Restrict to the current slice if needed
      lcu_start = MAX(lcu_start, child_state->slice->start_in_ts - child_state->tile->lcu_offset_in_ts);
      lcu_end = MIN(lcu_end, child_state->slice->end_in_ts - child_state->tile->lcu_offset_in_ts + 1);
      
      //Restrict to the current wavefront row if needed
      if (child_state->type == ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
        lcu_start = MAX(lcu_start, (child_state->wfrow->lcu_offset_y) * child_state->tile->cur_pic->width_in_lcu);
        lcu_end = MIN(lcu_end, (child_state->wfrow->lcu_offset_y + 1) * child_state->tile->cur_pic->width_in_lcu);
      }
      
      child_state->lcu_order_count = lcu_end - lcu_start;
      child_state->lcu_order = MALLOC(lcu_order_element, child_state->lcu_order_count);
      assert(child_state->lcu_order);
      
      for (i = 0; i < child_state->lcu_order_count; ++i) {
        lcu_id = lcu_start + i;
        child_state->lcu_order[i].encoder_state = child_state;
        child_state->lcu_order[i].id = lcu_id;
        child_state->lcu_order[i].index = i;
        child_state->lcu_order[i].position.x = lcu_id % child_state->tile->cur_pic->width_in_lcu;
        child_state->lcu_order[i].position.y = lcu_id / child_state->tile->cur_pic->width_in_lcu;
        child_state->lcu_order[i].position_px.x = child_state->lcu_order[i].position.x * LCU_WIDTH;
        child_state->lcu_order[i].position_px.y = child_state->lcu_order[i].position.y * LCU_WIDTH;
        child_state->lcu_order[i].size.x = MIN(LCU_WIDTH, encoder->in.width - (child_state->tile->lcu_offset_x * LCU_WIDTH + child_state->lcu_order[i].position_px.x));
        child_state->lcu_order[i].size.y = MIN(LCU_WIDTH, encoder->in.height - (child_state->tile->lcu_offset_y * LCU_WIDTH + child_state->lcu_order[i].position_px.y));
        child_state->lcu_order[i].first_row = lcu_in_first_row(child_state, child_state->tile->lcu_offset_in_ts + lcu_id);
        child_state->lcu_order[i].last_row = lcu_in_last_row(child_state, child_state->tile->lcu_offset_in_ts + lcu_id);
        child_state->lcu_order[i].first_column = lcu_in_first_column(child_state, child_state->tile->lcu_offset_in_ts + lcu_id);
        child_state->lcu_order[i].last_column = lcu_in_last_column(child_state, child_state->tile->lcu_offset_in_ts + lcu_id);
        
        child_state->lcu_order[i].above = NULL;
        child_state->lcu_order[i].below = NULL;
        child_state->lcu_order[i].left = NULL;
        child_state->lcu_order[i].right = NULL;
        
        if (!child_state->lcu_order[i].first_row) {
          //Find LCU above
          if (child_state->type == ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
            int j;
            for (j=0; child_state->parent->children[j].encoder_control; ++j) {
              if (child_state->parent->children[j].wfrow->lcu_offset_y == child_state->wfrow->lcu_offset_y - 1) {
                int k;
                for (k=0; k < child_state->parent->children[j].lcu_order_count; ++k) {
                  if (child_state->parent->children[j].lcu_order[k].position.x == child_state->lcu_order[i].position.x) {
                    assert(child_state->parent->children[j].lcu_order[k].position.y == child_state->lcu_order[i].position.y - 1);
                    child_state->lcu_order[i].above = &child_state->parent->children[j].lcu_order[k];
                  }
                }
              }
            }
          } else {
            child_state->lcu_order[i].above = &child_state->lcu_order[i-child_state->tile->cur_pic->width_in_lcu];
          }
          assert(child_state->lcu_order[i].above);
          child_state->lcu_order[i].above->below = &child_state->lcu_order[i];
        }
        if (!child_state->lcu_order[i].first_column) {
          child_state->lcu_order[i].left = &child_state->lcu_order[i-1];
          assert(child_state->lcu_order[i].left->position.x == child_state->lcu_order[i].position.x - 1);
          child_state->lcu_order[i].left->right = &child_state->lcu_order[i];
        }
      }
    } else {
      child_state->lcu_order_count = 0;
      child_state->lcu_order = NULL;
    }
  }
  
  //Validate the structure
  if (child_state->type == ENCODER_STATE_TYPE_TILE) {
    if (child_state->tile->lcu_offset_in_ts < child_state->slice->start_in_ts) {
      fprintf(stderr, "Tile %d starts before slice %d, in which it should be included!\n", child_state->tile->id, child_state->slice->id);
      return 0;
    }
    if (child_state->tile->lcu_offset_in_ts + child_state->tile->cur_pic->width_in_lcu * child_state->tile->cur_pic->height_in_lcu - 1 > child_state->slice->end_in_ts) {
      fprintf(stderr, "Tile %d ends after slice %d, in which it should be included!\n", child_state->tile->id, child_state->slice->id);
      return 0;
    }
  }
  
  if (child_state->type == ENCODER_STATE_TYPE_SLICE) {
    if (child_state->slice->start_in_ts < child_state->tile->lcu_offset_in_ts) {
      fprintf(stderr, "Slice %d starts before tile %d, in which it should be included!\n", child_state->slice->id, child_state->tile->id);
      return 0;
    }
    if (child_state->slice->end_in_ts > child_state->tile->lcu_offset_in_ts + child_state->tile->cur_pic->width_in_lcu * child_state->tile->cur_pic->height_in_lcu - 1) {
      fprintf(stderr, "Slice %d ends after tile %d, in which it should be included!\n", child_state->slice->id, child_state->tile->id);
      return 0;
    }
  }
  
  
#ifdef _DEBUG
  if (!parent_state) encoder_state_dump_graphviz(child_state);
#endif //_DEBUG
  return 1;
}

void encoder_state_finalize(encoder_state * const encoder_state) {
  if (encoder_state->children) {
    int i=0;
    for (i = 0; encoder_state->children[i].encoder_control; ++i) {
      encoder_state_finalize(&encoder_state->children[i]);
    }
    
    FREE_POINTER(encoder_state->children);
  }
  
  FREE_POINTER(encoder_state->lcu_order);
  encoder_state->lcu_order_count = 0;
  
  if (!encoder_state->parent || (encoder_state->parent->wfrow != encoder_state->wfrow)) {
    encoder_state_config_wfrow_finalize(encoder_state);
    FREE_POINTER(encoder_state->wfrow);
  }
  
  if (!encoder_state->parent || (encoder_state->parent->slice != encoder_state->slice)) {
    encoder_state_config_slice_finalize(encoder_state);
    FREE_POINTER(encoder_state->slice);
  }
  
  if (!encoder_state->parent || (encoder_state->parent->tile != encoder_state->tile)) {
    encoder_state_config_tile_finalize(encoder_state);
    FREE_POINTER(encoder_state->tile);
  }
  
  if (!encoder_state->parent || (encoder_state->parent->global != encoder_state->global)) {
    encoder_state_config_global_finalize(encoder_state);
    FREE_POINTER(encoder_state->global);
  }
  
  bitstream_finalize(&encoder_state->stream);
}


static void encoder_state_clear_refs(encoder_state *encoder_state) {
  //FIXME: Do we need to handle children? At present they all share the same global
  while (encoder_state->global->ref->used_size) {
    picture_list_rem(encoder_state->global->ref, encoder_state->global->ref->used_size - 1);
  }

  encoder_state->global->poc = 0;
}

static void encoder_state_blit_pixels(const encoder_state * const target_enc, pixel * const target, const encoder_state * const source_enc, const pixel * const source, const int is_y_channel) {
  const int source_offset_x = source_enc->tile->lcu_offset_x * LCU_WIDTH;
  const int source_offset_y = source_enc->tile->lcu_offset_y * LCU_WIDTH;
  
  const int target_offset_x = target_enc->tile->lcu_offset_x * LCU_WIDTH;
  const int target_offset_y = target_enc->tile->lcu_offset_y * LCU_WIDTH;
  
  int source_stride = source_enc->tile->cur_pic->width;
  int target_stride = target_enc->tile->cur_pic->width;
  
  int width;
  int height;
  
  int source_offset;
  int target_offset;
  
  //Do nothing if the source and the destination is the same!
  if (source_enc->tile == target_enc->tile) return;

  if (is_y_channel) {
    target_offset = source_offset_x + source_offset_y * target_enc->tile->cur_pic->width;
    source_offset = target_offset_x + target_offset_y * source_enc->tile->cur_pic->width;
  } else {
    target_offset = source_offset_x/2 + source_offset_y/2 * target_enc->tile->cur_pic->width/2;
    source_offset = target_offset_x/2 + target_offset_y/2 * source_enc->tile->cur_pic->width/2;
  }
  
  if (target_enc->children) {
    //Use information from the source
    width = MIN(source_enc->tile->cur_pic->width_in_lcu * LCU_WIDTH, target_enc->tile->cur_pic->width - source_offset_x);
    height = MIN(source_enc->tile->cur_pic->height_in_lcu * LCU_WIDTH, target_enc->tile->cur_pic->height - source_offset_y);
  } else {
    //Use information from the target
    width = MIN(target_enc->tile->cur_pic->width_in_lcu * LCU_WIDTH, source_enc->tile->cur_pic->width - target_offset_x);
    height = MIN(target_enc->tile->cur_pic->height_in_lcu * LCU_WIDTH, source_enc->tile->cur_pic->height - target_offset_y);
  }
  
  if (!is_y_channel) {
    width /= 2;
    height /= 2;
    
    source_stride /= 2;
    target_stride /= 2;
  }
  
  //picture_blit_pixels(source + source_offset, target + target_offset, width, height, source_enc->cur_pic->width, target_enc->cur_pic->width);
  picture_blit_pixels(source + source_offset, target + target_offset, width, height, source_stride, target_stride);
}



static void write_aud(encoder_state * const encoder_state)
{
  bitstream * const stream = &encoder_state->stream;
  encode_access_unit_delimiter(encoder_state);
  nal_write(stream, AUD_NUT, 0, 1);
  bitstream_align(stream);
}

static void encoder_state_recdata_to_bufs(encoder_state * const encoder_state, const lcu_order_element * const lcu, yuv_t * const hor_buf, yuv_t * const ver_buf) {
  picture* const cur_pic = encoder_state->tile->cur_pic;
  
  if (hor_buf) {
    const int rdpx = lcu->position_px.x;
    const int rdpy = lcu->position_px.y + lcu->size.y - 1;
    const int by = lcu->position.y;
    
    //Copy the bottom row of this LCU to the horizontal buffer
    picture_blit_pixels(&cur_pic->y_recdata[rdpy * cur_pic->width + rdpx],
                        &hor_buf->y[lcu->position_px.x + by * cur_pic->width],
                        lcu->size.x, 1, cur_pic->width, cur_pic->width);
    picture_blit_pixels(&cur_pic->u_recdata[(rdpy/2) * cur_pic->width/2 + (rdpx/2)],
                        &hor_buf->u[lcu->position_px.x / 2 + by * cur_pic->width / 2],
                        lcu->size.x / 2, 1, cur_pic->width / 2, cur_pic->width / 2);
    picture_blit_pixels(&cur_pic->v_recdata[(rdpy/2) * cur_pic->width/2 + (rdpx/2)],
                        &hor_buf->v[lcu->position_px.x / 2 + by * cur_pic->width / 2],
                        lcu->size.x / 2, 1, cur_pic->width / 2, cur_pic->width / 2);
  }
  
  if (ver_buf) {
    const int rdpx = lcu->position_px.x + lcu->size.x - 1;
    const int rdpy = lcu->position_px.y;
    const int bx = lcu->position.x;
    
    
    //Copy the right row of this LCU to the vertical buffer.
    picture_blit_pixels(&cur_pic->y_recdata[rdpy * cur_pic->width + rdpx],
                        &ver_buf->y[lcu->position_px.y + bx * cur_pic->height],
                        1, lcu->size.y, cur_pic->width, 1);
    picture_blit_pixels(&cur_pic->u_recdata[(rdpy/2) * cur_pic->width/2 + (rdpx/2)],
                        &ver_buf->u[lcu->position_px.y / 2 + bx * cur_pic->height / 2],
                        1, lcu->size.y / 2, cur_pic->width / 2, 1);
    picture_blit_pixels(&cur_pic->v_recdata[(rdpy/2) * cur_pic->width/2 + (rdpx/2)],
                        &ver_buf->v[lcu->position_px.y / 2 + bx * cur_pic->height / 2],
                        1, lcu->size.y / 2, cur_pic->width / 2, 1);
  }
  
}


static void worker_encoder_state_encode_lcu(void * opaque) {
  const lcu_order_element * const lcu = opaque;
  encoder_state *encoder_state = lcu->encoder_state;
  const encoder_control * const encoder = encoder_state->encoder_control;
  picture* const cur_pic = encoder_state->tile->cur_pic;
  
  //This part doesn't write to bitstream, it's only search, deblock and sao
  
  search_lcu(encoder_state, lcu->position_px.x, lcu->position_px.y, encoder_state->tile->hor_buf_search, encoder_state->tile->ver_buf_search);
    
  encoder_state_recdata_to_bufs(encoder_state, lcu, encoder_state->tile->hor_buf_search, encoder_state->tile->ver_buf_search);

  if (encoder->deblock_enable) {
    filter_deblock_lcu(encoder_state, lcu->position_px.x, lcu->position_px.y);
  }

  if (encoder->sao_enable) {
    const int stride = cur_pic->width_in_lcu;
    sao_info *sao_luma = &cur_pic->sao_luma[lcu->position.y * stride + lcu->position.x];
    sao_info *sao_chroma = &cur_pic->sao_chroma[lcu->position.y * stride + lcu->position.x];
    init_sao_info(sao_luma);
    init_sao_info(sao_chroma);

    {
      sao_info *sao_top =  lcu->position.y != 0 ? &cur_pic->sao_luma[(lcu->position.y - 1) * stride + lcu->position.x] : NULL;
      sao_info *sao_left = lcu->position.x != 0 ? &cur_pic->sao_luma[lcu->position.y * stride + lcu->position.x -1] : NULL;
      sao_search_luma(encoder_state, cur_pic, lcu->position.x, lcu->position.y, sao_luma, sao_top, sao_left);
    }

    {
      sao_info *sao_top =  lcu->position.y != 0 ? &cur_pic->sao_chroma[(lcu->position.y - 1) * stride + lcu->position.x] : NULL;
      sao_info *sao_left = lcu->position.x != 0 ? &cur_pic->sao_chroma[lcu->position.y * stride + lcu->position.x - 1] : NULL;
      sao_search_chroma(encoder_state, cur_pic, lcu->position.x, lcu->position.y, sao_chroma, sao_top, sao_left);
    }

    // Merge only if both luma and chroma can be merged
    sao_luma->merge_left_flag = sao_luma->merge_left_flag & sao_chroma->merge_left_flag;
    sao_luma->merge_up_flag = sao_luma->merge_up_flag & sao_chroma->merge_up_flag;
  }
  
  
  //Now write data to bitstream (required to have a correct CABAC state)
  
  //First LCU, and we are in a slice. We need a slice header
  if (encoder_state->type == ENCODER_STATE_TYPE_SLICE && lcu->index == 0) {
    encode_slice_header(encoder_state);
    bitstream_align(&encoder_state->stream); 
  }
  
  //Encode SAO
  if (encoder->sao_enable) {
    encode_sao(encoder_state, lcu->position.x, lcu->position.y, &cur_pic->sao_luma[lcu->position.y * cur_pic->width_in_lcu + lcu->position.x], &cur_pic->sao_chroma[lcu->position.y * cur_pic->width_in_lcu + lcu->position.x]);
  }
  
  //Encode coding tree
  encode_coding_tree(encoder_state, lcu->position.x << MAX_DEPTH, lcu->position.y << MAX_DEPTH, 0);

  //Terminator
  if (lcu->index < encoder_state->lcu_order_count - 1) {
    //Since we don't handle slice segments, end of slice segment == end of slice
    //Always 0 since otherwise it would be split
    cabac_encode_bin_trm(&encoder_state->cabac, 0);  // end_of_slice_segment_flag
  }
  
  //Wavefronts need the context to be copied to the next row
  if (encoder_state->type == ENCODER_STATE_TYPE_WAVEFRONT_ROW && lcu->index == 1) {
    int j;
    //Find next encoder (next row)
    for (j=0; encoder_state->parent->children[j].encoder_control; ++j) {
      if (encoder_state->parent->children[j].wfrow->lcu_offset_y == encoder_state->wfrow->lcu_offset_y + 1) {
        //And copy context
        context_copy(&encoder_state->parent->children[j], encoder_state);
      }
    }
  }
  
  if (encoder->sao_enable && lcu->above) {
    //If we're not the first in the row
    if (lcu->above->left) {
      encoder_state_recdata_to_bufs(encoder_state, lcu->above->left, encoder_state->tile->hor_buf_before_sao, NULL);
    }
    //Latest LCU in the row, copy the data from the one above also
    if (!lcu->right) {
      encoder_state_recdata_to_bufs(encoder_state, lcu->above, encoder_state->tile->hor_buf_before_sao, NULL);
    }
  }
}

static void encoder_state_encode_leaf(encoder_state * const encoder_state) {
  const encoder_control * const encoder = encoder_state->encoder_control;
#ifndef NDEBUG
//  const unsigned long long int debug_bitstream_position = bitstream_tell(&(encoder_state->stream));
#endif
  
  int i = 0;
  
  assert(encoder_state->is_leaf);
  assert(encoder_state->lcu_order_count > 0);
  
  //If we're not using wavefronts, or we have a WAVEFRONT_ROW which is the single child of its parent, than we should not use parallelism
  if (encoder_state->type != ENCODER_STATE_TYPE_WAVEFRONT_ROW || (encoder_state->type == ENCODER_STATE_TYPE_WAVEFRONT_ROW && !encoder_state->parent->children[1].encoder_control)) {
    for (i = 0; i < encoder_state->lcu_order_count; ++i) {
      PERFORMANCE_MEASURE_START();

      worker_encoder_state_encode_lcu(&encoder_state->lcu_order[i]);

#ifdef _DEBUG
      {
        const lcu_order_element * const lcu = &encoder_state->lcu_order[i];
        PERFORMANCE_MEASURE_END(encoder_state->encoder_control->threadqueue, "type=search_lcu,frame=%d,tile=%d,slice=%d,position_x=%d,position_y=%d", encoder_state->global->frame, encoder_state->tile->id, encoder_state->slice->id, lcu->position.x + encoder_state->tile->lcu_offset_x, lcu->position.y + encoder_state->tile->lcu_offset_y);
      }
#endif //_DEBUG
    }
    
    if (encoder->sao_enable) {
      PERFORMANCE_MEASURE_START();
      sao_reconstruct_frame(encoder_state);
      PERFORMANCE_MEASURE_END(encoder_state->encoder_control->threadqueue, "type=sao_reconstruct_frame,frame=%d,tile=%d,slice=%d,row=%d-%d", encoder_state->global->frame, encoder_state->tile->id, encoder_state->slice->id, encoder_state->lcu_order[0].position.y + encoder_state->tile->lcu_offset_y, encoder_state->lcu_order[encoder_state->lcu_order_count-1].position.y + encoder_state->tile->lcu_offset_y);
    }
  } else {
    for (i = 0; i < encoder_state->lcu_order_count; ++i) {
      const lcu_order_element * const lcu = &encoder_state->lcu_order[i];
#ifdef _DEBUG
      char job_description[256];
      sprintf(job_description, "type=search_lcu,frame=%d,tile=%d,slice=%d,row=%d,position_x=%d,position_y=%d", encoder_state->global->frame, encoder_state->tile->id, encoder_state->slice->id, encoder_state->wfrow->lcu_offset_y, lcu->position.x + encoder_state->tile->lcu_offset_x, lcu->position.y + encoder_state->tile->lcu_offset_y);
#else
      char* job_description = NULL;
#endif
      encoder_state->tile->wf_jobs[lcu->id] = threadqueue_submit(encoder_state->encoder_control->threadqueue, worker_encoder_state_encode_lcu, (void*)lcu, 1, job_description);
      if (encoder_state->tile->wf_jobs[lcu->id]) {
        if (lcu->position.x > 0) {
          // Wait for the LCU on the left.
          threadqueue_job_dep_add(encoder_state->tile->wf_jobs[lcu->id], encoder_state->tile->wf_jobs[lcu->id - 1]);
        }
        if (lcu->position.y > 0) {
          if (lcu->position.x < encoder_state->tile->cur_pic->width_in_lcu - 1) {
            // Wait for the LCU to the top-right of this one.
            threadqueue_job_dep_add(encoder_state->tile->wf_jobs[lcu->id], encoder_state->tile->wf_jobs[lcu->id - encoder_state->tile->cur_pic->width_in_lcu + 1]);
          } else {
            // If there is no top-right LCU, wait for the one above.
            threadqueue_job_dep_add(encoder_state->tile->wf_jobs[lcu->id], encoder_state->tile->wf_jobs[lcu->id - encoder_state->tile->cur_pic->width_in_lcu]);
          }
        }
        threadqueue_job_unwait_job(encoder_state->encoder_control->threadqueue, encoder_state->tile->wf_jobs[lcu->id]);
      }
    }
  }


  
  //We should not have written to bitstream!
//  assert(debug_bitstream_position == bitstream_tell(&(encoder_state->stream)));
}

static void encoder_state_encode(encoder_state * const main_state);

static void worker_encoder_state_encode_children(void * opaque) {
  encoder_state *sub_state = opaque;
  encoder_state_encode(sub_state);
  if (sub_state->is_leaf) {
    if (sub_state->type != ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
      PERFORMANCE_MEASURE_START();
      encoder_state_write_bitstream_leaf(sub_state);
      PERFORMANCE_MEASURE_END(sub_state->encoder_control->threadqueue, "type=encoder_state_write_bitstream_leaf,frame=%d,tile=%d,slice=%d,row=%d-%d", sub_state->global->frame, sub_state->tile->id, sub_state->slice->id, sub_state->lcu_order[0].position.y + sub_state->tile->lcu_offset_y, sub_state->lcu_order[sub_state->lcu_order_count-1].position.y + sub_state->tile->lcu_offset_y);
    } else {
      threadqueue_job *job;
#ifdef _DEBUG
      char job_description[256];
      sprintf(job_description, "type=encoder_state_write_bitstream_leaf,frame=%d,tile=%d,slice=%d,row=%d", sub_state->global->frame, sub_state->tile->id, sub_state->slice->id, sub_state->wfrow->lcu_offset_y);
#else
      char* job_description = NULL;
#endif
      job = threadqueue_submit(sub_state->encoder_control->threadqueue, worker_encoder_state_write_bitstream_leaf, sub_state, 1, job_description);
      threadqueue_job_dep_add(job, sub_state->tile->wf_jobs[sub_state->wfrow->lcu_offset_y * sub_state->tile->cur_pic->width_in_lcu + sub_state->lcu_order_count - 1]);
      threadqueue_job_unwait_job(sub_state->encoder_control->threadqueue, job);
      return;
    }
  }
}

typedef struct {
  int y;
  const encoder_state * encoder_state;
} worker_sao_reconstruct_lcu_data;

// ./kvazaar -i /scratch/h265-encode/pedestrian_area_1080p25.yuv  --input-res 1920x1080 -o /tmp/out.h265 --qp 23 -p 60  --frames 10
// Processed 10 frames,    5063552 bits AVG PSNR: 42.9771 46.0609 48.0985
// Total time: 19.440 s.
void worker_sao_reconstruct_lcu(void *opaque) {
  worker_sao_reconstruct_lcu_data *data = opaque;
  picture * const cur_pic = data->encoder_state->tile->cur_pic;
  unsigned stride = cur_pic->width_in_lcu;
  int x;
  
  //TODO: copy only needed data
  pixel *new_y_data = MALLOC(pixel, cur_pic->width * cur_pic->height);
  pixel *new_u_data = MALLOC(pixel, (cur_pic->width * cur_pic->height) >> 2);
  pixel *new_v_data = MALLOC(pixel, (cur_pic->width * cur_pic->height) >> 2);
  
  const int offset = cur_pic->width * (data->y*LCU_WIDTH);
  const int offset_c = cur_pic->width/2 * (data->y*LCU_WIDTH_C);
  int num_pixels = cur_pic->width * (LCU_WIDTH + 2);
  
  if (num_pixels + offset > cur_pic->width * cur_pic->height) {
    num_pixels = cur_pic->width * cur_pic->height - offset;
  }
  
  memcpy(&new_y_data[offset], &cur_pic->y_recdata[offset], sizeof(pixel) * num_pixels);
  memcpy(&new_u_data[offset_c], &cur_pic->u_recdata[offset_c], sizeof(pixel) * num_pixels >> 2);
  memcpy(&new_v_data[offset_c], &cur_pic->v_recdata[offset_c], sizeof(pixel) * num_pixels >> 2);
  
  if (data->y>0) {
    //copy first row from buffer
    memcpy(&new_y_data[cur_pic->width * (data->y*LCU_WIDTH-1)], &data->encoder_state->tile->hor_buf_before_sao->y[cur_pic->width * (data->y-1)], cur_pic->width * sizeof(pixel));
    memcpy(&new_u_data[cur_pic->width/2 * (data->y*LCU_WIDTH_C-1)], &data->encoder_state->tile->hor_buf_before_sao->u[cur_pic->width/2 * (data->y-1)], cur_pic->width/2 * sizeof(pixel));
    memcpy(&new_v_data[cur_pic->width/2 * (data->y*LCU_WIDTH_C-1)], &data->encoder_state->tile->hor_buf_before_sao->v[cur_pic->width/2 * (data->y-1)], cur_pic->width/2 * sizeof(pixel));
  }
  //assertions to be sure everything's ok for the next line (don't bother with last one)
  /*  These assertions may not be true if the row are not processed in order. To avoid having an artificial dependency between rows, it's better to remove them.
  assert((data->y >= cur_pic->height_in_lcu - 1) || memcmp(&data->encoder_state->tile->hor_buf_before_sao->y[cur_pic->width * (data->y)], &cur_pic->y_recdata[cur_pic->width * ((data->y + 1)*LCU_WIDTH-1)], cur_pic->width * sizeof(pixel))==0);
  assert((data->y >= cur_pic->height_in_lcu - 1) || memcmp(&data->encoder_state->tile->hor_buf_before_sao->u[cur_pic->width/2 * (data->y)], &cur_pic->u_recdata[cur_pic->width/2 * ((data->y + 1)*LCU_WIDTH_C-1)], cur_pic->width/2 * sizeof(pixel))==0);
  assert((data->y >= cur_pic->height_in_lcu - 1) || memcmp(&data->encoder_state->tile->hor_buf_before_sao->v[cur_pic->width/2 * (data->y)], &cur_pic->v_recdata[cur_pic->width/2 * ((data->y + 1)*LCU_WIDTH_C-1)], cur_pic->width/2 * sizeof(pixel))==0);*/

  for (x = 0; x < cur_pic->width_in_lcu; x++) {
  // sao_do_rdo(encoder, lcu.x, lcu.y, sao_luma, sao_chroma);
    sao_info *sao_luma = &cur_pic->sao_luma[data->y * stride + x];
    sao_info *sao_chroma = &cur_pic->sao_chroma[data->y * stride + x];
    sao_reconstruct(data->encoder_state->encoder_control, cur_pic, new_y_data, x, data->y, sao_luma, COLOR_Y);
    sao_reconstruct(data->encoder_state->encoder_control, cur_pic, new_u_data, x, data->y, sao_chroma, COLOR_U);
    sao_reconstruct(data->encoder_state->encoder_control, cur_pic, new_v_data, x, data->y, sao_chroma, COLOR_V);
  }
  
  free(new_y_data);
  free(new_u_data);
  free(new_v_data);

  free(opaque);
}


static int tree_is_a_chain(const encoder_state * const encoder_state) {
  if (!encoder_state->children[0].encoder_control) return 1;
  if (encoder_state->children[1].encoder_control) return 0;
  return tree_is_a_chain(&encoder_state->children[0]);
}

static void encoder_state_encode(encoder_state * const main_state) {
  //If we have children, encode at child level
  if (main_state->children[0].encoder_control) {
    int i=0;
    //If we have only one child, than it cannot be the last split in tree
    int node_is_the_last_split_in_tree = (main_state->children[1].encoder_control != 0);
    
    for (i=0; main_state->children[i].encoder_control; ++i) {
      encoder_state *sub_state = &(main_state->children[i]);
      
      if (sub_state->tile != main_state->tile) {
        encoder_state_blit_pixels(sub_state, sub_state->tile->cur_pic->y_data, main_state, main_state->tile->cur_pic->y_data, 1);
        encoder_state_blit_pixels(sub_state, sub_state->tile->cur_pic->u_data, main_state, main_state->tile->cur_pic->u_data, 0);
        encoder_state_blit_pixels(sub_state, sub_state->tile->cur_pic->v_data, main_state, main_state->tile->cur_pic->v_data, 0);
      }
      
      //To be the last split, we require that every child is a chain
      node_is_the_last_split_in_tree = node_is_the_last_split_in_tree && tree_is_a_chain(&main_state->children[i]);
    }
    //If it's the latest split point
    if (node_is_the_last_split_in_tree) {
      for (i=0; main_state->children[i].encoder_control; ++i) {
        //If we don't have wavefronts, parallelize encoding of children.
        if (main_state->children[i].type != ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
#ifdef _DEBUG
          char job_description[256];
          switch (main_state->children[i].type) {
            case ENCODER_STATE_TYPE_TILE: 
              sprintf(job_description, "frame=%d,tile=%d,row=%d-%d,position_x=%d,position_y=%d", main_state->children[i].global->frame, main_state->children[i].tile->id, main_state->children[i].lcu_order[0].position.y + main_state->children[i].tile->lcu_offset_y, main_state->children[i].lcu_order[main_state->children[i].lcu_order_count-1].position.y + main_state->children[i].tile->lcu_offset_y, main_state->children[i].tile->lcu_offset_x, main_state->children[i].tile->lcu_offset_y);
              break;
            case ENCODER_STATE_TYPE_SLICE:
              sprintf(job_description, "frame=%d,slice=%d,start_in_ts=%d", main_state->children[i].global->frame, main_state->children[i].slice->id, main_state->children[i].slice->start_in_ts);
              break;
            default:
              sprintf(job_description, "frame=%d,invalid", main_state->children[i].global->frame);
              break;
          }
#else
          char* job_description = NULL;
#endif
          threadqueue_submit(main_state->encoder_control->threadqueue, worker_encoder_state_encode_children, &(main_state->children[i]), 0, job_description);
        } else {
          //Wavefront rows have parallelism at LCU level, so we should not launch multiple threads here!
          //FIXME: add an assert: we can only have wavefront children
          worker_encoder_state_encode_children(&(main_state->children[i]));
        }
      }
      
      //If children are wavefront, we need to reconstruct SAO
      if (main_state->encoder_control->sao_enable && main_state->children[0].type == ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
        int y;
        picture * const cur_pic = main_state->tile->cur_pic;
        threadqueue_job *previous_job = NULL;
        
        for (y = 0; y < cur_pic->height_in_lcu; ++y) {
          worker_sao_reconstruct_lcu_data *data = MALLOC(worker_sao_reconstruct_lcu_data, 1);
          threadqueue_job *job;
#ifdef _DEBUG
          char job_description[256];
          sprintf(job_description, "frame=%d,tile=%d,position_y=%d", main_state->global->frame, main_state->tile->id, y + main_state->tile->lcu_offset_y);
#else
          char* job_description = NULL;
#endif
          data->y = y;
          data->encoder_state = main_state;
          
          job = threadqueue_submit(main_state->encoder_control->threadqueue, worker_sao_reconstruct_lcu, data, 1, job_description);
          
          if (previous_job) {
            threadqueue_job_dep_add(job, previous_job);
          }
          previous_job = job;
          
          if (y < cur_pic->height_in_lcu - 1) {
            //Not last row: depend on the last LCU of the row below
            threadqueue_job_dep_add(job, main_state->tile->wf_jobs[(y + 1) * cur_pic->width_in_lcu + cur_pic->width_in_lcu - 1]);
          } else {
            //Last row: depend on the last LCU of the row
            threadqueue_job_dep_add(job, main_state->tile->wf_jobs[(y + 0) * cur_pic->width_in_lcu + cur_pic->width_in_lcu - 1]);
          }
          threadqueue_job_unwait_job(main_state->encoder_control->threadqueue, job);
          
        }
      }
      threadqueue_flush(main_state->encoder_control->threadqueue);
    } else {
      for (i=0; main_state->children[i].encoder_control; ++i) {
        worker_encoder_state_encode_children(&(main_state->children[i]));
      }
    }
    
    for (i=0; main_state->children[i].encoder_control; ++i) {
      encoder_state *sub_state = &(main_state->children[i]);
      if (sub_state->tile != main_state->tile) {
        encoder_state_blit_pixels(main_state, main_state->tile->cur_pic->y_recdata, sub_state, sub_state->tile->cur_pic->y_recdata, 1);
        encoder_state_blit_pixels(main_state, main_state->tile->cur_pic->u_recdata, sub_state, sub_state->tile->cur_pic->u_recdata, 0);
        encoder_state_blit_pixels(main_state, main_state->tile->cur_pic->v_recdata, sub_state, sub_state->tile->cur_pic->v_recdata, 0);
      }
    }
  } else {
    switch (main_state->type) {
      case ENCODER_STATE_TYPE_TILE:
      case ENCODER_STATE_TYPE_SLICE:
      case ENCODER_STATE_TYPE_WAVEFRONT_ROW:
        encoder_state_encode_leaf(main_state);
        break;
      default:
        fprintf(stderr, "Unsupported leaf type %c!\n", main_state->type);
        assert(0);
    }
  }
}

static void encoder_state_new_frame(encoder_state * const main_state) {
  int i;
  //FIXME Move this somewhere else!
  if (main_state->type == ENCODER_STATE_TYPE_MAIN) {
    const encoder_control * const encoder = main_state->encoder_control;
    
    const int is_first_frame = (main_state->global->frame == 0);
    const int is_i_radl = (encoder->cfg->intra_period == 1 && main_state->global->frame % 2 == 0);
    const int is_p_radl = (encoder->cfg->intra_period > 1 && (main_state->global->frame % encoder->cfg->intra_period) == 0);
    main_state->global->is_radl_frame = is_first_frame || is_i_radl || is_p_radl;
    
    if (main_state->global->is_radl_frame) {
      // Clear the reference list
      encoder_state_clear_refs(main_state);

      main_state->global->slicetype = SLICE_I;
      main_state->global->pictype = NAL_IDR_W_RADL;
    } else {
      main_state->global->slicetype = encoder->cfg->intra_period==1 ? SLICE_I : SLICE_P;
      main_state->global->pictype = NAL_TRAIL_R;
    }
  } else {
    //Clear the bitstream if it's not the main encoder
    bitstream_clear(&main_state->stream);
  }
  
  if (main_state->is_leaf) {
    //Leaf states have cabac and context
    cabac_start(&main_state->cabac);
    init_contexts(main_state, main_state->global->QP, main_state->global->slicetype);

    // Initialize lambda value(s) to use in search
    encoder_state_init_lambda(main_state);
  }
  for (i = 0; main_state->children[i].encoder_control; ++i) {
    encoder_state_new_frame(&main_state->children[i]);
  }
  

}

static void encoder_state_write_bitstream_main(encoder_state * const main_state) {
  const encoder_control * const encoder = main_state->encoder_control;
  bitstream * const stream = &main_state->stream;

  int i;


  if (main_state->global->is_radl_frame) {
    // Access Unit Delimiter (AUD)
    if (encoder->aud_enable)
      write_aud(main_state);

    // Video Parameter Set (VPS)
    nal_write(stream, NAL_VPS_NUT, 0, 1);
    encode_vid_parameter_set(main_state);
    bitstream_align(stream);

    // Sequence Parameter Set (SPS)
    nal_write(stream, NAL_SPS_NUT, 0, 1);
    encode_seq_parameter_set(main_state);
    bitstream_align(stream);

    // Picture Parameter Set (PPS)
    nal_write(stream, NAL_PPS_NUT, 0, 1);
    encode_pic_parameter_set(main_state);
    bitstream_align(stream);

    if (main_state->global->frame == 0) {
      // Prefix SEI
      nal_write(stream, PREFIX_SEI_NUT, 0, 0);
      encode_prefix_sei_version(main_state);
      bitstream_align(stream);
    }
  } else {
    // Access Unit Delimiter (AUD)
    if (encoder->aud_enable)
      write_aud(main_state);
  }

  {
    // Not quite sure if this is correct, but it seems to have worked so far
    // so I tried to not change it's behavior.
    int long_start_code = main_state->global->is_radl_frame || encoder->aud_enable ? 0 : 1;

    nal_write(stream,
              main_state->global->is_radl_frame ? NAL_IDR_W_RADL : NAL_TRAIL_R, 0, long_start_code);
  }
  
  for (i = 0; main_state->children[i].encoder_control; ++i) {
    //Append bitstream to main stream
    bitstream_append(&main_state->stream, &main_state->children[i].stream);
    //FIXME: Move this...
    bitstream_clear(&main_state->children[i].stream);
  }
  
  // Calculate checksum
  add_checksum(main_state);

  //FIXME: Why is this needed?
  main_state->tile->cur_pic->poc = main_state->global->poc;
}

static void worker_encoder_state_write_bitstream_leaf(void * opaque) {
  encoder_state_write_bitstream_leaf((encoder_state *) opaque);
}

static void encoder_state_write_bitstream_leaf(encoder_state * const encoder_state) {
  const encoder_control * const encoder = encoder_state->encoder_control;
  //Write terminator of the leaf
  assert(encoder_state->is_leaf);
  
  //Last LCU
  {
    const lcu_order_element * const lcu = &encoder_state->lcu_order[encoder_state->lcu_order_count - 1];
    const int lcu_addr_in_ts = lcu->id + encoder_state->tile->lcu_offset_in_ts;
    const int end_of_slice_segment_flag = lcu_at_slice_end(encoder, lcu_addr_in_ts);
  
    cabac_encode_bin_trm(&encoder_state->cabac, end_of_slice_segment_flag);  // end_of_slice_segment_flag
  
    if (!end_of_slice_segment_flag) {
      assert(lcu_at_tile_end(encoder, lcu_addr_in_ts) || lcu->position.x == (encoder_state->tile->cur_pic->width_in_lcu - 1));
      cabac_encode_bin_trm(&encoder_state->cabac, 1); // end_of_sub_stream_one_bit == 1
      cabac_flush(&encoder_state->cabac);
    } else {
      cabac_flush(&encoder_state->cabac);
      bitstream_align(&encoder_state->stream);
    }
  }
}

static void encoder_state_write_bitstream_tile(encoder_state * const main_state) {
  //If it's not a leaf, a tile is "nothing". We only have to write sub elements
  int i;
  for (i = 0; main_state->children[i].encoder_control; ++i) {
    //Append bitstream to main stream
    bitstream_append(&main_state->stream, &main_state->children[i].stream);
  }
}

static void encoder_state_write_bitstream_slice(encoder_state * const main_state) {
  int i;
  encode_slice_header(main_state);
  bitstream_align(&main_state->stream); 
  
  for (i = 0; main_state->children[i].encoder_control; ++i) {
    //Append bitstream to main stream
    bitstream_append(&main_state->stream, &main_state->children[i].stream);
  }
}


static void encoder_state_write_bitstream(encoder_state * const main_state) {
  int i;
  if (!main_state->is_leaf) {
    for (i=0; main_state->children[i].encoder_control; ++i) {
      encoder_state *sub_state = &(main_state->children[i]);
      encoder_state_write_bitstream(sub_state);
    }
    
    switch (main_state->type) {
      case ENCODER_STATE_TYPE_MAIN:
        encoder_state_write_bitstream_main(main_state);
        break;
      case ENCODER_STATE_TYPE_TILE:
        encoder_state_write_bitstream_tile(main_state);
        break;
      case ENCODER_STATE_TYPE_SLICE:
        encoder_state_write_bitstream_slice(main_state);
        break;
      default:
        fprintf(stderr, "Unsupported node type %c!\n", main_state->type);
        assert(0);
    }
  }
}

void encode_one_frame(encoder_state * const main_state)
{
  encoder_state_new_frame(main_state);
  encoder_state_encode(main_state);
  encoder_state_write_bitstream(main_state);
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
  unsigned array_width = encoder_state->tile->cur_pic->width;
  unsigned array_height = encoder_state->tile->cur_pic->height;

  if (width != array_width) {
    // In the case of frames not being aligned on 8 bit borders, bits need to be copied to fill them in.
    if (!read_and_fill_frame_data(file, width, height, array_width,
                                  encoder_state->tile->cur_pic->y_data) ||
        !read_and_fill_frame_data(file, width >> 1, height >> 1, array_width >> 1,
                                  encoder_state->tile->cur_pic->u_data) ||
        !read_and_fill_frame_data(file, width >> 1, height >> 1, array_width >> 1,
                                  encoder_state->tile->cur_pic->v_data))
      return 0;
  } else {
    // Otherwise the data can be read directly to the array.
    unsigned y_size = width * height;
    unsigned uv_size = (width >> 1) * (height >> 1);
    if (y_size  != fread(encoder_state->tile->cur_pic->y_data, sizeof(unsigned char),
                         y_size, file) ||
        uv_size != fread(encoder_state->tile->cur_pic->u_data, sizeof(unsigned char),
                         uv_size, file) ||
        uv_size != fread(encoder_state->tile->cur_pic->v_data, sizeof(unsigned char),
                         uv_size, file))
      return 0;
  }

  if (height != array_height) {
    fill_after_frame(height, array_width, array_height,
                     encoder_state->tile->cur_pic->y_data);
    fill_after_frame(height >> 1, array_width >> 1, array_height >> 1,
                     encoder_state->tile->cur_pic->u_data);
    fill_after_frame(height >> 1, array_width >> 1, array_height >> 1,
                     encoder_state->tile->cur_pic->v_data);
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
  const picture * const cur_pic = encoder_state->tile->cur_pic;
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
  uint8_t pic_type = encoder_state->global->slicetype == SLICE_I ? 0
                   : encoder_state->global->slicetype == SLICE_P ? 1
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
  const encoder_control * const encoder = encoder_state->encoder_control;
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
  WRITE_SE(stream, ((int8_t)encoder_state->global->QP)-26, "pic_init_qp_minus26");
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
  //FIXME: use encoder_control instead of cur_pic
  const picture * const cur_pic = encoder_state->tile->cur_pic;

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

void encoder_next_frame(encoder_state *encoder_state) {
  const encoder_control * const encoder = encoder_state->encoder_control;
  picture *old_pic;
  
  // Remove the ref pic (if present)
  if (encoder_state->global->ref->used_size == (uint32_t)encoder->cfg->ref_frames) {
    picture_list_rem(encoder_state->global->ref, encoder_state->global->ref->used_size-1);
  }
  // Add current picture as reference
  picture_list_add(encoder_state->global->ref, encoder_state->tile->cur_pic);
  // Allocate new memory to current picture
  old_pic = encoder_state->tile->cur_pic;
  // TODO: reuse memory from old reference
  encoder_state->tile->cur_pic = picture_alloc(encoder_state->tile->cur_pic->width, encoder_state->tile->cur_pic->height, encoder_state->tile->cur_pic->width_in_lcu, encoder_state->tile->cur_pic->height_in_lcu);

  //FIXME: does the coeff_* really belongs to cur_pic?
  // Copy pointer from the last cur_pic because we don't want to reallocate it
  MOVE_POINTER(encoder_state->tile->cur_pic->coeff_y,old_pic->coeff_y);
  MOVE_POINTER(encoder_state->tile->cur_pic->coeff_u,old_pic->coeff_u);
  MOVE_POINTER(encoder_state->tile->cur_pic->coeff_v,old_pic->coeff_v);
  
  picture_free(old_pic);

  encoder_state->global->frame++;
  encoder_state->global->poc++;
}

static void encoder_state_entry_points_explore(const encoder_state * const encoder_state, int * const r_count, int * const r_max_length) {
  int i;
  for (i = 0; encoder_state->children[i].encoder_control; ++i) {
    if (encoder_state->children[i].is_leaf) {
      const int my_length = bitstream_tell(&encoder_state->children[i].stream)/8;
      ++(*r_count);
      if (my_length > *r_max_length) {
        *r_max_length = my_length;
      }
    } else {
      encoder_state_entry_points_explore(&encoder_state->children[i], r_count, r_max_length);
    }
  }
}

static void encoder_state_entry_points_write(bitstream * const stream, const encoder_state * const encoder_state, const int num_entry_points, const int write_length, int * const r_count) {
  int i;
  for (i = 0; encoder_state->children[i].encoder_control; ++i) {
    if (encoder_state->children[i].is_leaf) {
      const int my_length = bitstream_tell(&encoder_state->children[i].stream)/8;
      ++(*r_count);
      //Don't write the last one
      if (*r_count < num_entry_points) {
        WRITE_U(stream, my_length - 1, write_length, "entry_point_offset-minus1")
      }
    } else {
      encoder_state_entry_points_write(stream, &encoder_state->children[i], num_entry_points, write_length, r_count);
    }
  }
}

static int num_bitcount(unsigned int n) {
  int pos = 0;
  if (n >= 1<<16) { n >>= 16; pos += 16; }
  if (n >= 1<< 8) { n >>=  8; pos +=  8; }
  if (n >= 1<< 4) { n >>=  4; pos +=  4; }
  if (n >= 1<< 2) { n >>=  2; pos +=  2; }
  if (n >= 1<< 1) {           pos +=  1; }
  return ((n == 0) ? (-1) : pos);
}

void encode_slice_header(encoder_state * const encoder_state)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  bitstream * const stream = &encoder_state->stream;

#ifdef _DEBUG
  printf("=========== Slice ===========\n");
#endif
  WRITE_U(stream, (encoder_state->slice->start_in_rs == 0), 1, "first_slice_segment_in_pic_flag");

  if (encoder_state->global->pictype >= NAL_BLA_W_LP
      && encoder_state->global->pictype <= NAL_RSV_IRAP_VCL23) {
    WRITE_U(stream, 1, 1, "no_output_of_prior_pics_flag");
  }

  WRITE_UE(stream, 0, "slice_pic_parameter_set_id");
  if (encoder_state->slice->start_in_rs > 0) {
    //For now, we don't support dependent slice segments
    //WRITE_U(stream, 0, 1, "dependent_slice_segment_flag");
    WRITE_UE(stream, encoder_state->slice->start_in_rs, "slice_segment_address");
  }

  WRITE_UE(stream, encoder_state->global->slicetype, "slice_type");

  // if !entropy_slice_flag

    //if output_flag_present_flag
      //WRITE_U(stream, 1, 1, "pic_output_flag");
    //end if
    //if( IdrPicFlag ) <- nal_unit_type == 5
  if (encoder_state->global->pictype != NAL_IDR_W_RADL
      && encoder_state->global->pictype != NAL_IDR_N_LP) {
      int j;
      int ref_negative = encoder_state->global->ref->used_size;
      int ref_positive = 0;
      WRITE_U(stream, encoder_state->global->poc&0xf, 4, "pic_order_cnt_lsb");
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
    WRITE_U(stream, 1, 1, "slice_sao_luma_flag");
    WRITE_U(stream, 1, 1, "slice_sao_chroma_flag");
  }

  if (encoder_state->global->slicetype != SLICE_I) {
      WRITE_U(stream, 1, 1, "num_ref_idx_active_override_flag");
        WRITE_UE(stream, encoder_state->global->ref->used_size-1, "num_ref_idx_l0_active_minus1");
      WRITE_UE(stream, 5-MRG_MAX_NUM_CANDS, "five_minus_max_num_merge_cand");
  }

  if (encoder_state->global->slicetype == SLICE_B) {
      WRITE_U(stream, 0, 1, "mvd_l1_zero_flag");
  }

  // Skip flags that are not present
  // if !entropy_slice_flag
    WRITE_SE(stream, 0, "slice_qp_delta");
    //WRITE_U(stream, 1, 1, "alignment");
   
  if (encoder->tiles_enable || encoder->wpp) {
    int num_entry_points = 0;
    int max_length_seen = 0;
    
    encoder_state_entry_points_explore(encoder_state, &num_entry_points, &max_length_seen);
    
    WRITE_UE(stream, num_entry_points - 1, "num_entry_point_offsets");
    if (num_entry_points > 0) {
      int entry_points_written = 0;
      int offset_len = num_bitcount(max_length_seen) + 1;
      WRITE_UE(stream, offset_len - 1, "offset_len_minus1");
      encoder_state_entry_points_write(stream, encoder_state, num_entry_points, offset_len, &entry_points_written); 
    }
  }
}


static void encode_sao_color(encoder_state * const encoder_state, sao_info *sao,
                             color_index color_i)
{
  cabac_data * const cabac = &encoder_state->cabac;
  sao_eo_cat i;

  // Skip colors with no SAO.
  //FIXME: for now, we always have SAO for all channels
  if (color_i == COLOR_Y && 0) return;
  if (color_i != COLOR_Y && 0) return;

  /// sao_type_idx_luma:   TR, cMax = 2, cRiceParam = 0, bins = {0, bypass}
  /// sao_type_idx_chroma: TR, cMax = 2, cRiceParam = 0, bins = {0, bypass}
  // Encode sao_type_idx for Y and U+V.
  if (color_i != COLOR_V) {
    cabac->ctx = &(cabac->ctx_sao_type_idx_model);;
    CABAC_BIN(cabac, sao->type != SAO_TYPE_NONE, "sao_type_idx");
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

static void encode_sao_merge_flags(encoder_state * const encoder_state, sao_info *sao, unsigned x_ctb, unsigned y_ctb)
{
  cabac_data * const cabac = &encoder_state->cabac;
  // SAO merge flags are not present for the first row and column.
  if (x_ctb > 0) {
    cabac->ctx = &(cabac->ctx_sao_merge_flag_model);
    CABAC_BIN(cabac, sao->merge_left_flag, "sao_merge_left_flag");
  }
  if (y_ctb > 0 && !sao->merge_left_flag) {
    cabac->ctx = &(cabac->ctx_sao_merge_flag_model);
    CABAC_BIN(cabac, sao->merge_up_flag, "sao_merge_up_flag");
  }
}

/**
 * \brief Encode SAO information.
 */
static void encode_sao(encoder_state * const encoder_state,
                       unsigned x_lcu, uint16_t y_lcu,
                       sao_info *sao_luma, sao_info *sao_chroma)
{
  // TODO: transmit merge flags outside sao_info
  encode_sao_merge_flags(encoder_state, sao_luma, x_lcu, y_lcu);

  // If SAO is merged, nothing else needs to be coded.
  if (!sao_luma->merge_left_flag && !sao_luma->merge_up_flag) {
    encode_sao_color(encoder_state, sao_luma, COLOR_Y);
    encode_sao_color(encoder_state, sao_chroma, COLOR_U);
    encode_sao_color(encoder_state, sao_chroma, COLOR_V);
  }
}


void encode_coding_tree(encoder_state * const encoder_state,
                        uint16_t x_ctb, uint16_t y_ctb, uint8_t depth)
{
  cabac_data * const cabac = &encoder_state->cabac;
  const picture * const cur_pic = encoder_state->tile->cur_pic;
  cu_info *cur_cu = &cur_pic->cu_array[x_ctb + y_ctb * (cur_pic->width_in_lcu << MAX_DEPTH)];
  uint8_t split_flag = GET_SPLITDATA(cur_cu, depth);
  uint8_t split_model = 0;
  
  //Absolute ctb
  uint16_t abs_x_ctb = x_ctb + (encoder_state->tile->lcu_offset_x * LCU_WIDTH) / (LCU_WIDTH >> MAX_DEPTH);
  uint16_t abs_y_ctb = y_ctb + (encoder_state->tile->lcu_offset_y * LCU_WIDTH) / (LCU_WIDTH >> MAX_DEPTH);

  // Check for slice border FIXME
  uint8_t border_x = ((encoder_state->encoder_control->in.width) < (abs_x_ctb * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth))) ? 1 : 0;
  uint8_t border_y = ((encoder_state->encoder_control->in.height) < (abs_y_ctb * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> depth))) ? 1 : 0;
  uint8_t border_split_x = ((encoder_state->encoder_control->in.width)  < ((abs_x_ctb + 1) * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> (depth + 1)))) ? 0 : 1;
  uint8_t border_split_y = ((encoder_state->encoder_control->in.height) < ((abs_y_ctb + 1) * (LCU_WIDTH >> MAX_DEPTH) + (LCU_WIDTH >> (depth + 1)))) ? 0 : 1;
  uint8_t border = border_x | border_y; /*!< are we in any border CU */

  // When not in MAX_DEPTH, insert split flag and split the blocks if needed
  if (depth != MAX_DEPTH) {
    // Implisit split flag when on border
    if (!border) {
      // Get left and top block split_flags and if they are present and true, increase model number
      if (x_ctb > 0 && GET_SPLITDATA(&(cur_pic->cu_array[x_ctb - 1 + y_ctb * (cur_pic->width_in_lcu << MAX_DEPTH)]), depth) == 1) {
        split_model++;
      }

      if (y_ctb > 0 && GET_SPLITDATA(&(cur_pic->cu_array[x_ctb + (y_ctb - 1) * (cur_pic->width_in_lcu << MAX_DEPTH)]), depth) == 1) {
        split_model++;
      }

      cabac->ctx = &(cabac->ctx_split_flag_model[split_model]);
      CABAC_BIN(cabac, split_flag, "SplitFlag");
    }

    if (split_flag || border) {
      // Split blocks and remember to change x and y block positions
      uint8_t change = 1<<(MAX_DEPTH-1-depth);
      encode_coding_tree(encoder_state, x_ctb, y_ctb, depth + 1); // x,y

      // TODO: fix when other half of the block would not be completely over the border
      if (!border_x || border_split_x) {
        encode_coding_tree(encoder_state, x_ctb + change, y_ctb, depth + 1);
      }
      if (!border_y || border_split_y) {
        encode_coding_tree(encoder_state, x_ctb, y_ctb + change, depth + 1);
      }
      if (!border || (border_split_x && border_split_y)) {
        encode_coding_tree(encoder_state, x_ctb + change, y_ctb + change, depth + 1);
      }
      return;
    }
  }



    // Encode skip flag
  if (encoder_state->global->slicetype != SLICE_I) {
    int8_t ctx_skip = 0; // uiCtxSkip = aboveskipped + leftskipped;
    int ui;
    int16_t num_cand = MRG_MAX_NUM_CANDS;
    // Get left and top skipped flags and if they are present and true, increase context number
    if (x_ctb > 0 && (&cur_pic->cu_array[x_ctb - 1 + y_ctb * (cur_pic->width_in_lcu << MAX_DEPTH)])->skipped) {
      ctx_skip++;
    }

    if (y_ctb > 0 && (&cur_pic->cu_array[x_ctb + (y_ctb - 1) * (cur_pic->width_in_lcu << MAX_DEPTH)])->skipped) {
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
  if (encoder_state->global->slicetype != SLICE_I) {
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
            if (encoder_state->global->ref->used_size != 1) { //encoder_state->ref_idx_num[uiRefListIdx] != 1)//NumRefIdx != 1)
              // parseRefFrmIdx
              int32_t ref_frame = cur_cu->inter.mv_ref;

              cabac->ctx = &(cabac->ctx_cu_ref_pic_model[0]);
              CABAC_BIN(cabac, (ref_frame != 0), "ref_frame_flag");

              if (ref_frame > 0) {
                int32_t i;
                int32_t ref_num = encoder_state->global->ref->used_size - 2;

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

            if (!(/*pcCU->getSlice()->getMvdL1ZeroFlag() &&*/ encoder_state->global->ref_list == REF_PIC_LIST_1 && cur_cu->inter.mv_dir == 3)) {
              const int32_t mvd_hor = cur_cu->inter.mvd[0];
              const int32_t mvd_ver = cur_cu->inter.mvd[1];
              const int8_t hor_abs_gr0 = mvd_hor != 0;
              const int8_t ver_abs_gr0 = mvd_ver != 0;
              const uint32_t mvd_hor_abs = abs(mvd_hor);
              const uint32_t mvd_ver_abs = abs(mvd_ver);

              cabac->ctx = &(cabac->ctx_cu_mvd_model[0]);
              CABAC_BIN(cabac, (mvd_hor != 0), "abs_mvd_greater0_flag_hor");
              CABAC_BIN(cabac, (mvd_ver != 0), "abs_mvd_greater0_flag_ver");

              cabac->ctx = &(cabac->ctx_cu_mvd_model[1]);

              if (hor_abs_gr0) {
                CABAC_BIN(cabac, (mvd_hor_abs>1), "abs_mvd_greater1_flag_hor");
              }

              if (ver_abs_gr0) {
                CABAC_BIN(cabac, (mvd_ver_abs>1), "abs_mvd_greater1_flag_ver");
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

    {
      int cbf = (cbf_is_set(cur_cu->cbf.y, depth) ||
                 cbf_is_set(cur_cu->cbf.u, depth) ||
                 cbf_is_set(cur_cu->cbf.v, depth));

      // Only need to signal coded block flag if not skipped or merged
      // skip = no coded residual, merge = coded residual
      if (!cur_cu->merged) {
        cabac->ctx = &(cabac->ctx_cu_qt_root_cbf_model);
        CABAC_BIN(cabac, cbf, "rqt_root_cbf");
      }
      // Code (possible) coeffs to bitstream

      if (cbf) {
        encode_transform_coeff(encoder_state, x_ctb * 2, y_ctb * 2, depth, 0, 0, 0);
      }
    }

    // END for each part
  } else if (cur_cu->type == CU_INTRA) {
    uint8_t intra_pred_mode[4] = {
      cur_cu->intra[0].mode, cur_cu->intra[1].mode,
      cur_cu->intra[2].mode, cur_cu->intra[3].mode };
      uint8_t intra_pred_mode_chroma = cur_cu->intra[0].mode_chroma;
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
        left_cu = &cur_pic->cu_array[x_ctb - 1 + y_ctb * (cur_pic->width_in_lcu << MAX_DEPTH)];
      }
      // Don't take the above CU across the LCU boundary.
      if (y_ctb > 0 && (y_ctb & 7) != 0) {
        above_cu = &cur_pic->cu_array[x_ctb + (y_ctb - 1) * (cur_pic->width_in_lcu << MAX_DEPTH)];
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

      if (intra_pred_mode_chroma == intra_pred_mode[0]) {
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

      // pred_mode == 5 mean intra_pred_mode_chroma is something that can't
      // be coded.
      assert(pred_mode != 5);

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

    encode_transform_coeff(encoder_state, x_ctb * 2, y_ctb * 2, depth, 0, 0, 0);
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
    assert(0);
    exit(1);
  }

   /* end prediction unit */
  /* end coding_unit */
}


coeff_scan_order_t get_scan_order(int8_t cu_type, int intra_mode, int depth)
{
  // Scan mode is diagonal, except for 4x4+8x8 luma and 4x4 chroma, where:
  // - angular 6-14 = vertical
  // - angular 22-30 = horizontal
  if (cu_type == CU_INTRA && depth >= 3) {
    if (intra_mode >= 6 && intra_mode <= 14) {
      return SCAN_VER;
    } else if (intra_mode >= 22 && intra_mode <= 30) {
      return SCAN_HOR;
    }
  }

  return SCAN_DIAG;
}


static void encode_transform_unit(encoder_state * const encoder_state,
                                  int x_pu, int y_pu, int depth)
{
  const picture * const cur_pic = encoder_state->tile->cur_pic;
  uint8_t width = LCU_WIDTH >> depth;
  uint8_t width_c = (depth == MAX_PU_DEPTH ? width : width / 2);

  int x_cu = x_pu / 2;
  int y_cu = y_pu / 2;
  cu_info *cur_cu = &cur_pic->cu_array[x_cu + y_cu * (cur_pic->width_in_lcu << MAX_DEPTH)];

  coefficient coeff_y[LCU_WIDTH*LCU_WIDTH+1];
  coefficient coeff_u[LCU_WIDTH*LCU_WIDTH>>2];
  coefficient coeff_v[LCU_WIDTH*LCU_WIDTH>>2];
  int32_t coeff_stride = cur_pic->width;

  int8_t scan_idx = get_scan_order(cur_cu->type, cur_cu->intra[PU_INDEX(x_pu, y_pu)].mode, depth);

  int cbf_y = cbf_is_set(cur_cu->cbf.y, depth + PU_INDEX(x_pu, y_pu));

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

  // CoeffNxN
  // Residual Coding
  if (cbf_y) {
    encode_coeff_nxn(encoder_state, coeff_y, width, 0, scan_idx, cur_cu->intra[PU_INDEX(x_pu, y_pu)].tr_skip);
  }

  if (depth == MAX_DEPTH + 1 && !(x_pu % 2 && y_pu % 2)) {
    // For size 4x4 luma transform the corresponding chroma transforms are
    // also of size 4x4 covering 8x8 luma pixels. The residual is coded
    // in the last transform unit so for the other ones, don't do anything.
    return;
  }

  if (cbf_is_set(cur_cu->cbf.u, depth) || cbf_is_set(cur_cu->cbf.v, depth)) {
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

    scan_idx = get_scan_order(cur_cu->type, cur_cu->intra[0].mode_chroma, depth);

    if (cbf_is_set(cur_cu->cbf.u, depth)) {
      encode_coeff_nxn(encoder_state, coeff_u, width_c, 2, scan_idx, 0);
    }

    if (cbf_is_set(cur_cu->cbf.v, depth)) {
      encode_coeff_nxn(encoder_state, coeff_v, width_c, 2, scan_idx, 0);
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
void encode_transform_coeff(encoder_state * const encoder_state, int32_t x_pu,int32_t y_pu,
                            int8_t depth, int8_t tr_depth, uint8_t parent_coeff_u, uint8_t parent_coeff_v)
{
  cabac_data * const cabac = &encoder_state->cabac;
  int32_t x_cu = x_pu / 2;
  int32_t y_cu = y_pu / 2;
  const picture * const cur_pic = encoder_state->tile->cur_pic;
  cu_info *cur_cu = &cur_pic->cu_array[x_cu + y_cu * (cur_pic->width_in_lcu << MAX_DEPTH)];

  // NxN signifies implicit transform split at the first transform level.
  // There is a similar implicit split for inter, but it is only used when
  // transform hierarchy is not in use.
  int intra_split_flag = (cur_cu->type == CU_INTRA && cur_cu->part_size == SIZE_NxN);

  // The implicit split by intra NxN is not counted towards max_tr_depth.
  int max_tr_depth = (cur_cu->type == CU_INTRA ? TR_DEPTH_INTRA + intra_split_flag : TR_DEPTH_INTER);

  int8_t split = (cur_cu->tr_depth > depth);

  const int cb_flag_y = cbf_is_set(cur_cu->cbf.y, depth + PU_INDEX(x_pu, y_pu));
  const int cb_flag_u = cbf_is_set(cur_cu->cbf.u, depth);
  const int cb_flag_v = cbf_is_set(cur_cu->cbf.v, depth);

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
    encode_transform_coeff(encoder_state, x_pu, y_pu, depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
    encode_transform_coeff(encoder_state, x_pu + pu_offset, y_pu,  depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
    encode_transform_coeff(encoder_state, x_pu, y_pu + pu_offset,  depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
    encode_transform_coeff(encoder_state, x_pu + pu_offset, y_pu + pu_offset,  depth + 1, tr_depth + 1, cb_flag_u, cb_flag_v);
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
    encode_transform_unit(encoder_state, x_pu, y_pu, depth);
  }
}

void encode_coeff_nxn(encoder_state * const encoder_state, coefficient *coeff, uint8_t width,
                      uint8_t type, int8_t scan_mode, int8_t tr_skip)
{
  const encoder_control * const encoder = encoder_state->encoder_control;
  cabac_data * const cabac = &encoder_state->cabac;
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
  encode_last_significant_xy(encoder_state, last_coeff_x, last_coeff_y, width, width,
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
      CABAC_BIN(cabac, sig_coeff_group, "coded_sub_block_flag");
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
          CABAC_BIN(cabac, sig, "sig_coeff_flag");
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
        CABAC_BIN(cabac, symbol, "coeff_abs_level_greater1_flag");

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
          CABAC_BIN(cabac, symbol, "coeff_abs_level_greater2_flag");
        }
      }

      if (be_valid && sign_hidden) {
        CABAC_BINS_EP(cabac, (coeff_signs >> 1), (num_non_zero - 1), "coeff_sign_flag");
      } else {
        CABAC_BINS_EP(cabac, coeff_signs, num_non_zero, "coeff_sign_flag");
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
void encode_last_significant_xy(encoder_state * const encoder_state,
                                uint8_t lastpos_x, uint8_t lastpos_y,
                                uint8_t width, uint8_t height,
                                uint8_t type, uint8_t scan)
{
  cabac_data * const cabac = &encoder_state->cabac;
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
    CABAC_BIN(cabac,1,"last_sig_coeff_x_prefix");
  }

  if (group_idx_x < g_group_idx[width - 1]) {
    cabac->ctx = &base_ctx_x[offset_x + (last_x >> shift_x)];
    CABAC_BIN(cabac,0,"last_sig_coeff_x_prefix");
  }

  // Last Y binarization
  for (last_y = 0; last_y < group_idx_y ; last_y++) {
    cabac->ctx = &base_ctx_y[offset_y + (last_y >> shift_y)];
    CABAC_BIN(cabac,1,"last_sig_coeff_y_prefix");
  }

  if (group_idx_y < g_group_idx[height - 1]) {
    cabac->ctx = &base_ctx_y[offset_y + (last_y >> shift_y)];
    CABAC_BIN(cabac,0,"last_sig_coeff_y_prefix");
  }

  // Last X
  if (group_idx_x > 3) {
    lastpos_x -= g_min_in_group[group_idx_x];

    for (i = ((group_idx_x - 2) >> 1) - 1; i >= 0; i--) {
      CABAC_BIN_EP(cabac,(lastpos_x>>i) & 1,"last_sig_coeff_x_suffix");
    }
  }

  // Last Y
  if (group_idx_y > 3) {
    lastpos_y -= g_min_in_group[group_idx_y];

    for (i = ((group_idx_y - 2) >> 1) - 1; i >= 0; i--) {
      CABAC_BIN_EP(cabac,(lastpos_y>>i) & 1,"last_sig_coeff_y_suffix");
    }
  }

  // end LastSignificantXY
}
