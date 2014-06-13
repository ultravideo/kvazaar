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

#include "encoder_state-ctors_dtors.h"

#include <stdlib.h>

#include "encoderstate.h"


static int encoder_state_config_global_init(encoder_state * const encoder_state) {
  encoder_state->global->ref = image_list_alloc(MAX_REF_PIC_COUNT);
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
  image_list_destroy(encoder_state->global->ref);
}

static int encoder_state_config_tile_init(encoder_state * const encoder_state, 
                                          const int lcu_offset_x, const int lcu_offset_y,
                                          const int width, const int height, const int width_in_lcu, const int height_in_lcu) {
  
  const encoder_control * const encoder = encoder_state->encoder_control;
  encoder_state->tile->frame = videoframe_alloc(width, height, 0);
  
  if (encoder_state->type == ENCODER_STATE_TYPE_MAIN) {
    //If not a parent, then we can avoid keeping a copy of the image
    encoder_state->tile->frame->source = image_alloc(encoder_state->tile->frame->width, encoder_state->tile->frame->height, 0);
    encoder_state->tile->frame->rec = image_alloc(encoder_state->tile->frame->width, encoder_state->tile->frame->height, 0);
  } else {
    encoder_state->tile->frame->source = NULL;
    encoder_state->tile->frame->rec = NULL;
  }

  if (!encoder_state->tile->frame) {
    printf("Error allocating videoframe!\r\n");
    return 0;
  }
  
  // Init coeff data table
  //FIXME: move them
  encoder_state->tile->frame->coeff_y = MALLOC(coefficient, width * height);
  encoder_state->tile->frame->coeff_u = MALLOC(coefficient, (width * height) >> 2);
  encoder_state->tile->frame->coeff_v = MALLOC(coefficient, (width * height) >> 2);
  
  encoder_state->tile->lcu_offset_x = lcu_offset_x;
  encoder_state->tile->lcu_offset_y = lcu_offset_y;
  
  encoder_state->tile->lcu_offset_in_ts = encoder->tiles_ctb_addr_rs_to_ts[lcu_offset_x + lcu_offset_y * encoder->in.width_in_lcu];
  
  //Allocate buffers
  //order by row of (LCU_WIDTH * frame->width_in_lcu) pixels
  encoder_state->tile->hor_buf_search = yuv_t_alloc(LCU_WIDTH * encoder_state->tile->frame->width_in_lcu * encoder_state->tile->frame->height_in_lcu);
  //order by column of (LCU_WIDTH * encoder_state->height_in_lcu) pixels (there is no more extra pixel, since we can use a negative index)
  encoder_state->tile->ver_buf_search = yuv_t_alloc(LCU_WIDTH * encoder_state->tile->frame->height_in_lcu * encoder_state->tile->frame->width_in_lcu);
  
  if (encoder->sao_enable) {
    encoder_state->tile->hor_buf_before_sao = yuv_t_alloc(LCU_WIDTH * encoder_state->tile->frame->width_in_lcu * encoder_state->tile->frame->height_in_lcu);
  } else {
    encoder_state->tile->hor_buf_before_sao = NULL;
  }
  
  if (encoder->wpp) {
    encoder_state->tile->wf_jobs = MALLOC(threadqueue_job*, encoder_state->tile->frame->width_in_lcu * encoder_state->tile->frame->height_in_lcu);
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
  
  if (encoder_state->tile->frame->source) image_free(encoder_state->tile->frame->source);
  if (encoder_state->tile->frame->rec) image_free(encoder_state->tile->frame->rec);

  videoframe_free(encoder_state->tile->frame);
  encoder_state->tile->frame = NULL;
  
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
  child_state->tqj_bitstream_written = NULL;
  child_state->tqj_recon_done = NULL;
  
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
        end_in_ts = child_state->tile->frame->width_in_lcu * child_state->tile->frame->height_in_lcu;
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
        end_in_ts = child_state->tile->lcu_offset_in_ts + child_state->tile->frame->width_in_lcu * child_state->tile->frame->height_in_lcu;
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
      int lcu_end = child_state->tile->frame->width_in_lcu * child_state->tile->frame->height_in_lcu;
      
      //Restrict to the current slice if needed
      lcu_start = MAX(lcu_start, child_state->slice->start_in_ts - child_state->tile->lcu_offset_in_ts);
      lcu_end = MIN(lcu_end, child_state->slice->end_in_ts - child_state->tile->lcu_offset_in_ts + 1);
      
      //Restrict to the current wavefront row if needed
      if (child_state->type == ENCODER_STATE_TYPE_WAVEFRONT_ROW) {
        lcu_start = MAX(lcu_start, (child_state->wfrow->lcu_offset_y) * child_state->tile->frame->width_in_lcu);
        lcu_end = MIN(lcu_end, (child_state->wfrow->lcu_offset_y + 1) * child_state->tile->frame->width_in_lcu);
      }
      
      child_state->lcu_order_count = lcu_end - lcu_start;
      child_state->lcu_order = MALLOC(lcu_order_element, child_state->lcu_order_count);
      assert(child_state->lcu_order);
      
      for (i = 0; i < child_state->lcu_order_count; ++i) {
        lcu_id = lcu_start + i;
        child_state->lcu_order[i].encoder_state = child_state;
        child_state->lcu_order[i].id = lcu_id;
        child_state->lcu_order[i].index = i;
        child_state->lcu_order[i].position.x = lcu_id % child_state->tile->frame->width_in_lcu;
        child_state->lcu_order[i].position.y = lcu_id / child_state->tile->frame->width_in_lcu;
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
            //For all previous wavefront rows
            for (j=0; &child_state->parent->children[j] != child_state && child_state->parent->children[j].encoder_control; ++j) {
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
            child_state->lcu_order[i].above = &child_state->lcu_order[i-child_state->tile->frame->width_in_lcu];
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
    if (child_state->tile->lcu_offset_in_ts + child_state->tile->frame->width_in_lcu * child_state->tile->frame->height_in_lcu - 1 > child_state->slice->end_in_ts) {
      fprintf(stderr, "Tile %d ends after slice %d, in which it should be included!\n", child_state->tile->id, child_state->slice->id);
      return 0;
    }
  }
  
  if (child_state->type == ENCODER_STATE_TYPE_SLICE) {
    if (child_state->slice->start_in_ts < child_state->tile->lcu_offset_in_ts) {
      fprintf(stderr, "Slice %d starts before tile %d, in which it should be included!\n", child_state->slice->id, child_state->tile->id);
      return 0;
    }
    if (child_state->slice->end_in_ts > child_state->tile->lcu_offset_in_ts + child_state->tile->frame->width_in_lcu * child_state->tile->frame->height_in_lcu - 1) {
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