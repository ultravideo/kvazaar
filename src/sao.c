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

#include "sao.h"

#include <string.h>
#include <stdlib.h>

#include "picture.h"

// Offsets of a and b in relation to c.
// dir_offset[dir][a or b]
// |       |   a   | a     |     a |
// | a c b |   c   |   c   |   c   |
// |       |   b   |     b | b     |
static const vector2d g_sao_edge_offsets[SAO_NUM_EO][2] = {
  { { -1, 0 }, { 1, 0 } },
  { { 0, -1 }, { 0, 1 } },
  { { -1, -1 }, { 1, 1 } },
  { { 1, -1 }, { -1, 1 } }
};



void init_sao_info(sao_info *sao) {
  sao->type = SAO_TYPE_NONE;
  sao->merge_left_flag = 0;
  sao->merge_up_flag = 0;
}

// Mapping of edge_idx values to eo-classes.
static const unsigned g_sao_eo_idx_to_eo_category[] = { 1, 2, 0, 3, 4 };
// Mapping relationships between a, b and c to eo_idx.
#define EO_IDX(a, b, c) (2 + SIGN3((c) - (a)) + SIGN3((c) - (b)))

/**
 * \brief Check merge conditions
 */
static int sao_check_merge(const sao_info *sao_candidate, int type,
                           int offsets[NUM_SAO_EDGE_CATEGORIES],
                           int band_position, int eo_class)
{
  if (sao_candidate->type == type) {
    if (offsets[1] == sao_candidate->offsets[1] &&
        offsets[2] == sao_candidate->offsets[2] &&
        offsets[3] == sao_candidate->offsets[3] &&
        offsets[4] == sao_candidate->offsets[4]) {
      // Type must be BAND or EDGE
      if ((type == SAO_TYPE_BAND && band_position == sao_candidate->band_position) ||
        (type == SAO_TYPE_EDGE && eo_class == sao_candidate->eo_class)) {
          return 1;
      }
    }
  }
  return 0;
}

/**
 * \brief calculate an array of intensity correlations for each intensity value
 */
static void calc_sao_offset_array(const sao_info *sao, int *offset)
{
  int val;
  int values = (1<<g_bitdepth);
  int shift = g_bitdepth-5;

  // Loop through all intensity values and construct an offset array
  for (val = 0; val < values; val++) {
    int cur_band = val>>shift;
    if(cur_band >= sao->band_position && cur_band < sao->band_position+4) {
      offset[val] = CLIP(0, values-1,val+sao->offsets[cur_band-sao->band_position+1]);
    } else {
      offset[val] = val;
    }
  }
}


/**
 * \param orig_data  Original pixel data. 64x64 for luma, 32x32 for chroma.
 * \param rec_data  Reconstructed pixel data. 64x64 for luma, 32x32 for chroma.
 * \param sao_bands an array of bands for original and reconstructed block
 */
static int calc_sao_band_offsets(int sao_bands[2][32], int offsets[4],
                                 int *band_position, int *rate)
{
  int band;
  int offset;
  int best_dist;
  int temp_dist;
  int dist[32];
  int temp_offsets[32];
  int temp_rate[32];
  int best_dist_pos = 0;

  memset(dist, 0, 32*sizeof(int));
  memset(temp_rate, 0, 32*sizeof(int));

  // Calculate distortion for each band using N*h^2 - 2*h*E
  for (band = 0; band < 32; band++) {
    best_dist = INT_MAX;
    offset = 0;
    if (sao_bands[1][band] != 0) {
      offset = (sao_bands[0][band] + (sao_bands[1][band] >> 1)) / sao_bands[1][band];
      offset = CLIP(-SAO_ABS_OFFSET_MAX, SAO_ABS_OFFSET_MAX, offset);
    }
    dist[band] = offset==0?0:INT_MAX;
    temp_offsets[band] = 0;
    while(offset != 0) {
      temp_dist = sao_bands[1][band]*offset*offset - 2*offset*sao_bands[0][band];

      // Calculate used rate to code this
      temp_rate[band] = abs(offset+2);
      // When max offset, no need to signal final bit
      if (abs(offset)==SAO_ABS_OFFSET_MAX) {
        temp_rate[band] --;
      }
      // Store best distortion and offset
      if(temp_dist < best_dist) {
        dist[band] = temp_dist;
        temp_offsets[band] = offset;
      }
      offset += (offset > 0) ? -1:1;
    }
  }

  best_dist = INT_MAX;
  //Find starting pos for best 4 band distortions
  for (band = 0; band < 28; band++) {
    temp_dist = dist[band] + dist[band+1] + dist[band+2] + dist[band+3];
    if(temp_dist < best_dist) {
      best_dist = temp_dist;
      best_dist_pos = band;
    }
  }
  *rate = temp_rate[best_dist_pos]+temp_rate[best_dist_pos+1]+
          temp_rate[best_dist_pos+2]+temp_rate[best_dist_pos+3];
  // Copy best offsets to output
  memcpy(offsets, &temp_offsets[best_dist_pos], 4*sizeof(int));

  *band_position = best_dist_pos;

  return best_dist;
}

/**
 * \param orig_data  Original pixel data. 64x64 for luma, 32x32 for chroma.
 * \param rec_data  Reconstructed pixel data. 64x64 for luma, 32x32 for chroma.
 * \param sao_bands an array of bands for original and reconstructed block
 */
static void calc_sao_bands(const pixel *orig_data, const pixel *rec_data,
                           int block_width, int block_height,
                           int sao_bands[2][32])
{
  int y, x;
  int shift = g_bitdepth-5;

  //Loop pixels and take top 5 bits to classify different bands
  for (y = 0; y < block_height; ++y) {
    for (x = 0; x < block_width; ++x) {
      sao_bands[0][rec_data[y * block_width + x]>>shift] += orig_data[y * block_width + x] - rec_data[y * block_width + x];
      sao_bands[1][rec_data[y * block_width + x]>>shift]++;
    }
  }
}


/**
 * \param orig_data  Original pixel data. 64x64 for luma, 32x32 for chroma.
 * \param rec_data  Reconstructed pixel data. 64x64 for luma, 32x32 for chroma.
 * \param dir_offsets
 * \param is_chroma  0 for luma, 1 for chroma. Indicates
 */
static void calc_sao_edge_dir(const pixel *orig_data, const pixel *rec_data,
                              int eo_class, int block_width, int block_height,
                              int cat_sum_cnt[2][NUM_SAO_EDGE_CATEGORIES])
{
  int y, x;
  vector2d a_ofs = g_sao_edge_offsets[eo_class][0];
  vector2d b_ofs = g_sao_edge_offsets[eo_class][1];
  // Arrays orig_data and rec_data are quarter size for chroma.

  // Don't sample the edge pixels because this function doesn't have access to
  // their neighbours.
  for (y = 1; y < block_height - 1; ++y) {
    for (x = 1; x < block_width - 1; ++x) {
      const pixel *c_data = &rec_data[y * block_width + x];
      pixel a = c_data[a_ofs.y * block_width + a_ofs.x];
      pixel c = c_data[0];
      pixel b = c_data[b_ofs.y * block_width + b_ofs.x];

      int eo_idx = EO_IDX(a, b, c);
      int eo_cat = g_sao_eo_idx_to_eo_category[eo_idx];

      cat_sum_cnt[0][eo_cat] += orig_data[y * block_width + x] - c;
      cat_sum_cnt[1][eo_cat] += 1;
    }
  }
}

static void sao_reconstruct_color(const pixel *rec_data, pixel *new_rec_data,
                                  const sao_info *sao,
                                  int stride, int new_stride,
                                  int block_width, int block_height)
{
  int y, x;
  vector2d a_ofs = g_sao_edge_offsets[sao->eo_class][0];
  vector2d b_ofs = g_sao_edge_offsets[sao->eo_class][1];
  // Arrays orig_data and rec_data are quarter size for chroma.


  if(sao->type == SAO_TYPE_BAND) {
    int offsets[1<<BIT_DEPTH];
    calc_sao_offset_array(sao, offsets);
    for (y = 0; y < block_height; ++y) {
      for (x = 0; x < block_width; ++x) {
        new_rec_data[y * new_stride + x] = offsets[rec_data[y * stride + x]];
      }
    }
  } else {
    // Don't sample the edge pixels because this function doesn't have access to
    // their neighbours.
    for (y = 0; y < block_height; ++y) {
      for (x = 0; x < block_width; ++x) {
        const pixel *c_data = &rec_data[y * stride + x];
        pixel *new_data = &new_rec_data[y * new_stride + x];
        pixel a = c_data[a_ofs.y * stride + a_ofs.x];
        pixel c = c_data[0];
        pixel b = c_data[b_ofs.y * stride + b_ofs.x];

        int eo_idx = EO_IDX(a, b, c);
        int eo_cat = g_sao_eo_idx_to_eo_category[eo_idx];

        new_data[0] = (pixel)CLIP(0, (1 << BIT_DEPTH) - 1, c_data[0] + sao->offsets[eo_cat]);
      }
    }
  }
}

/**
 * \brief Calculate dimensions of the buffer used by sao reconstruction.

 * \param pic  Picture.
 * \param sao  Sao parameters.
 * \param rec  Top-left corner of the LCU
 */
static void sao_calc_band_block_dims(const picture *pic, color_index color_i,
                                     vector2d *rec, vector2d *block)
{
  const int is_chroma = (color_i != COLOR_Y ? 1 : 0);
  int width = pic->width >> is_chroma;
  int height = pic->height >> is_chroma;
  int block_width = LCU_WIDTH >> is_chroma;


  // Handle right and bottom, taking care of non-LCU sized CUs.
  if (rec->y + block_width >= height) {
    if (rec->y + block_width >= height) {
      block->y = height - rec->y;
    }
  }
  if (rec->x + block_width >= width) {
    if (rec->x + block_width > width) {
      block->x = width - rec->x;
    }
  }

  rec->x = 0; rec->y = 0;
}

/**
 * \brief Calculate dimensions of the buffer used by sao reconstruction.
 *
 * This function calculates 4 vectors that can be used to make the temporary
 * buffers required by sao_reconstruct_color.
 *
 * Vector block is the area affected by sao. Vectors tr and br are top-left
 * margin and bottom-right margin, which contain pixels that are not modified
 * by the reconstruction of this LCU but are needed by the reconstruction.
 * Vector rec is the offset from the CU to the required pixel area.
 *
 * The margins are always either 0 or 1, depending on the direction of the
 * edge offset class.
 *
 * This also takes into account borders of the picture and non-LCU sized
 * CU's at the bottom and right of the picture.
 *
 * \ CU + rec
 *  +------+
 *  |\ tl  |
 *  | +--+ |
 *  | |\ block
 *  | | \| |
 *  | +--+ |
 *  |     \ br
 *  +------+
 *
 * \param pic  Picture.
 * \param sao  Sao parameters.
 * \param rec  Top-left corner of the LCU, modified to be top-left corner of
 */
static void sao_calc_edge_block_dims(const picture *pic, color_index color_i,
                                     const sao_info *sao, vector2d *rec,
                                     vector2d *tl, vector2d *br,
                                     vector2d *block)
{
  vector2d a_ofs = g_sao_edge_offsets[sao->eo_class][0];
  vector2d b_ofs = g_sao_edge_offsets[sao->eo_class][1];
  const int is_chroma = (color_i != COLOR_Y ? 1 : 0);
  int width = pic->width >> is_chroma;
  int height = pic->height >> is_chroma;
  int block_width = LCU_WIDTH >> is_chroma;

  // Handle top and left.
  if (rec->y == 0) {
    tl->y = 0;
    if (a_ofs.y == -1 || b_ofs.y == -1) {
      block->y -= 1;
      tl->y += 1;
    }
  }
  if (rec->x == 0) {
    tl->x = 0;
    if (a_ofs.x == -1 || b_ofs.x == -1) {
      block->x -= 1;
      tl->x += 1;
    }
  }

  // Handle right and bottom, taking care of non-LCU sized CUs.
  if (rec->y + block_width >= height) {
    br->y = 0;
    if (rec->y + block_width >= height) {
      block->y = height - rec->y;
    }
    if (a_ofs.y == 1 || b_ofs.y == 1) {
      block->y -= 1;
      br->y += 1;
    }
  }
  if (rec->x + block_width >= width) {
    br->x = 0;
    if (rec->x + block_width > width) {
      block->x = width - rec->x;
    }
    if (a_ofs.x == 1 || b_ofs.x == 1) {
      block->x -= 1;
      br->x += 1;
    }
  }

  rec->y = (rec->y == 0 ? 0 : -1);
  rec->x = (rec->x == 0 ? 0 : -1);
}

void sao_reconstruct(picture *pic, const pixel *old_rec,
                     unsigned x_ctb, unsigned y_ctb,
                     const sao_info *sao, color_index color_i)
{
  const int is_chroma = (color_i != COLOR_Y ? 1 : 0);
  const int pic_stride = pic->width >> is_chroma;
  const int lcu_stride = LCU_WIDTH >> is_chroma;
  const int buf_stride = lcu_stride + 2;

  pixel *recdata = pic->recdata[color_i];
  pixel buf_rec[(LCU_WIDTH + 2) * (LCU_WIDTH + 2)];
  pixel new_rec[LCU_WIDTH * LCU_WIDTH];
  // Calling CU_TO_PIXEL with depth 1 is the same as using block size of 32.
  pixel *lcu_rec = &recdata[CU_TO_PIXEL(x_ctb, y_ctb, is_chroma, pic_stride)];
  const pixel *old_lcu_rec = &old_rec[CU_TO_PIXEL(x_ctb, y_ctb, is_chroma, pic_stride)];

  vector2d ofs;
  vector2d tl = { 1, 1 };
  vector2d br = { 1, 1 };
  vector2d block = { LCU_WIDTH, LCU_WIDTH };

  if (sao->type == SAO_TYPE_NONE) {
    return;
  }

  ofs.x = x_ctb * lcu_stride;
  ofs.y = y_ctb * lcu_stride;
  block.x = lcu_stride;
  block.y = lcu_stride;
  if (sao->type == SAO_TYPE_BAND) {
    tl.x = 0; tl.y = 0;
    br.x = 0; br.y = 0;
    sao_calc_band_block_dims(pic, color_i, &ofs, &block);
  }
  else {
    sao_calc_edge_block_dims(pic, color_i, sao, &ofs, &tl, &br, &block);
  }

  // Data to tmp buffer.
  picture_blit_pixels(&old_lcu_rec[ofs.y * pic_stride + ofs.x],
                      buf_rec,
                      tl.x + block.x + br.x,
                      tl.y + block.y + br.y,
                      pic_stride, buf_stride);

  sao_reconstruct_color(&buf_rec[tl.y * buf_stride + tl.x],
                        &new_rec[(ofs.y + tl.y) * lcu_stride + ofs.x + tl.x],
                        sao,
                        buf_stride, lcu_stride,
                        block.x, block.y);

  // Copy reconstructed block from tmp buffer to rec image.
  picture_blit_pixels(&new_rec[(tl.y + ofs.y) * lcu_stride + (tl.x + ofs.x)],
                      &lcu_rec[(tl.y + ofs.y) * pic_stride + (tl.x + ofs.x)],
                      block.x, block.y, lcu_stride, pic_stride);
}


/**
 * \param data     Array of pointers to reference pixels.
 * \param recdata  Array of pointers to reconstructed pixels.
 * \param block_width   Width of the area to be examined.
 * \param block_height  Height of the area to be examined.
 * \param buf_cnt  Number of pointers data and recdata have.
 * \param sao_out  Output parameter for the best sao parameters.
 */
static void sao_search_best_mode(const pixel * data[], const pixel * recdata[],
                                 int block_width, int block_height,
                                 unsigned buf_cnt,
                                 sao_info *sao_out, sao_info *sao_top,
                                 sao_info *sao_left)
{
  sao_eo_class edge_class;
  // This array is used to calculate the mean offset used to minimize distortion.
  int cat_sum_cnt[2][NUM_SAO_EDGE_CATEGORIES];
  unsigned i = 0;
  int temp_rate = 0;
  memset(cat_sum_cnt, 0, sizeof(int) * 2 * NUM_SAO_EDGE_CATEGORIES);

  sao_out->ddistortion = INT_MAX;

  for (edge_class = SAO_EO0; edge_class <= SAO_EO3; ++edge_class) {
    int edge_offset[NUM_SAO_EDGE_CATEGORIES];
    int sum_ddistortion = 0;
    sao_eo_cat edge_cat;

    // Call calc_sao_edge_dir once for luma and twice for chroma.
    for (i = 0; i < buf_cnt; ++i) {
      calc_sao_edge_dir(data[i], recdata[i], edge_class,
                        block_width, block_height, cat_sum_cnt);
    }

    for (edge_cat = SAO_EO_CAT1; edge_cat <= SAO_EO_CAT4; ++edge_cat) {
      int cat_sum = cat_sum_cnt[0][edge_cat];
      int cat_cnt = cat_sum_cnt[1][edge_cat];

      // The optimum offset can be calculated by getting the minima of the
      // fast ddistortion estimation formula. The minima is the mean error
      // and we round that to the nearest integer.
      int offset = 0;
      if (cat_cnt != 0) {
        offset = (cat_sum + (cat_cnt >> 1)) / cat_cnt;
        offset = CLIP(-SAO_ABS_OFFSET_MAX, SAO_ABS_OFFSET_MAX, offset);
      }

      // Sharpening edge offsets can't be encoded, so set them to 0 here.
      if (edge_cat >= SAO_EO_CAT1 && edge_cat <= SAO_EO_CAT2 && offset < 0) {
        offset = 0;
      }
      if (edge_cat >= SAO_EO_CAT3 && edge_cat <= SAO_EO_CAT4 && offset > 0) {
        offset = 0;
      }

      edge_offset[edge_cat] = offset;
      // The ddistortion is amount by which the SSE of data changes. It should
      // be negative for all categories, if offset was chosen correctly.
      // ddistortion = N * h^2 - 2 * h * E, where N is the number of samples
      // and E is the sum of errors.
      // It basically says that all pixels that are not improved by offset
      // increase increase SSE by h^2 and all pixels that are improved by
      // offset decrease SSE by h*E.
      sum_ddistortion += cat_cnt * offset * offset - 2 * offset * cat_sum;

      // Calculate used rate to code this
      temp_rate += abs(offset+2);
      // When max offset, no need to signal final bit
      if (abs(offset)==SAO_ABS_OFFSET_MAX) {
        temp_rate --;
      }
    }
    // If can merge, use lower cost
    if((sao_top != NULL && sao_check_merge(sao_top, SAO_TYPE_EDGE, edge_offset,0, edge_class)) ||
       (sao_left != NULL && sao_check_merge(sao_left, SAO_TYPE_EDGE, edge_offset,0, edge_class))) {
          temp_rate = 1;
    }
    sum_ddistortion += (int)((double)temp_rate*(g_cur_lambda_cost+0.5));
    // SAO is not applied for category 0.
    edge_offset[SAO_EO_CAT0] = 0;

    // Choose the offset class that offers the least error after offset.
    if (sum_ddistortion < sao_out->ddistortion) {
      sao_out->eo_class = edge_class;
      sao_out->ddistortion = sum_ddistortion;
      memcpy(sao_out->offsets, edge_offset, sizeof(int) * NUM_SAO_EDGE_CATEGORIES);
    }
  }

  // Band offset
  {
    int sao_bands[2][32];
    int temp_offsets[4];
    int ddistortion;
    int temp_rate = 0;

    memset(sao_bands, 0, 2 * 32 * sizeof(int));
    for (i = 0; i < buf_cnt; ++i) {
      calc_sao_bands(data[i], recdata[i],block_width,
                     block_height,sao_bands);
    }

    ddistortion = calc_sao_band_offsets(sao_bands, temp_offsets, &sao_out->band_position, &temp_rate);

    if((sao_top != NULL && sao_check_merge(sao_top, SAO_TYPE_BAND, temp_offsets,sao_out->band_position, 0)) ||
       (sao_left != NULL && sao_check_merge(sao_left, SAO_TYPE_BAND, temp_offsets,sao_out->band_position, 0))) {
          temp_rate = 1;
    }
    ddistortion += (int)((double)temp_rate*(g_cur_lambda_cost+0.5));

    // Select band sao over edge sao when distortion is lower
    if (ddistortion < sao_out->ddistortion) {
      sao_out->type = SAO_TYPE_BAND;
      sao_out->ddistortion = ddistortion;
      memcpy(&sao_out->offsets[1], temp_offsets, sizeof(int) * 4);
    }
  }

  // When BD-rate would increase because of SAO, disable it
  if(sao_out->ddistortion >= 0) {
    sao_out->type = SAO_TYPE_NONE;
    return;
  }

  // Check for merge modes
  if (sao_top != NULL) {
    sao_out->merge_up_flag = sao_check_merge(sao_top, sao_out->type, sao_out->offsets,
                                             sao_out->band_position, sao_out->eo_class);
  }

  // Check for merge modes
  if (sao_left != NULL) {
    sao_out->merge_left_flag = sao_check_merge(sao_left, sao_out->type, sao_out->offsets,
                                               sao_out->band_position, sao_out->eo_class);
  }
}

 void sao_search_chroma(const picture *pic, unsigned x_ctb, unsigned y_ctb, sao_info *sao, sao_info *sao_top, sao_info *sao_left)
{
  int block_width  = (LCU_WIDTH / 2);
  int block_height = (LCU_WIDTH / 2);
  const pixel *orig_list[2];
  const pixel *rec_list[2];
  pixel orig[2][LCU_CHROMA_SIZE];
  pixel rec[2][LCU_CHROMA_SIZE];
  color_index color_i;

  // Check for right and bottom boundaries.
  if (x_ctb * (LCU_WIDTH / 2) + (LCU_WIDTH / 2) >= (unsigned)pic->width / 2) {
    block_width = (pic->width - x_ctb * LCU_WIDTH) / 2;
  }
  if (y_ctb * (LCU_WIDTH / 2) + (LCU_WIDTH / 2) >= (unsigned)pic->height / 2) {
    block_height = (pic->height - y_ctb * LCU_WIDTH) / 2;
  }

  sao->type = SAO_TYPE_EDGE;

  // Copy data to temporary buffers and init orig and rec lists to point to those buffers.
  for (color_i = COLOR_U; color_i <= COLOR_V; ++color_i) {
    pixel *data = &pic->data[color_i][CU_TO_PIXEL(x_ctb, y_ctb, 1, pic->width / 2)];
    pixel *recdata = &pic->recdata[color_i][CU_TO_PIXEL(x_ctb, y_ctb, 1, pic->width / 2)];
    picture_blit_pixels(data, orig[color_i - 1], block_width, block_height,
                        pic->width / 2, LCU_WIDTH / 2);
    picture_blit_pixels(recdata, rec[color_i - 1], block_width, block_height,
                        pic->width / 2, LCU_WIDTH / 2);
    orig_list[color_i - 1] = &orig[color_i - 1][0];
    rec_list[color_i - 1] = &rec[color_i - 1][0];
  }

  // Calculate
  sao_search_best_mode(orig_list, rec_list, block_width, block_height, 2, sao, sao_top, sao_left);
}

void sao_search_luma(const picture *pic, unsigned x_ctb, unsigned y_ctb, sao_info *sao, sao_info *sao_top, sao_info *sao_left)
{
  pixel orig[LCU_LUMA_SIZE];
  pixel rec[LCU_LUMA_SIZE];
  const pixel * orig_list[1] = { NULL };
  const pixel * rec_list[1] = { NULL };
  pixel *data = &pic->y_data[CU_TO_PIXEL(x_ctb, y_ctb, 0, pic->width)];
  pixel *recdata = &pic->y_recdata[CU_TO_PIXEL(x_ctb, y_ctb, 0, pic->width)];
  int block_width = LCU_WIDTH;
  int block_height = LCU_WIDTH;

  // Check for right and bottom boundaries.
  if (x_ctb * LCU_WIDTH + LCU_WIDTH >= (unsigned)pic->width) {
    block_width = pic->width - x_ctb * LCU_WIDTH;
  }
  if (y_ctb * LCU_WIDTH + LCU_WIDTH >= (unsigned)pic->height) {
    block_height = pic->height - y_ctb * LCU_WIDTH;
  }

  sao->type = SAO_TYPE_EDGE;

  // Fill temporary buffers with picture data.
  picture_blit_pixels(data, orig, block_width, block_height, pic->width, LCU_WIDTH);
  picture_blit_pixels(recdata, rec, block_width, block_height, pic->width, LCU_WIDTH);

  orig_list[0] = orig;
  rec_list[0] = rec;
  sao_search_best_mode(orig_list, rec_list, block_width, block_height, 1, sao, sao_top, sao_left);
}

void sao_reconstruct_frame(const encoder_control * const encoder)
{
  vector2d lcu;
  picture *pic = encoder->in.cur_pic;

  // These are needed because SAO needs the pre-SAO pixels form left and
  // top LCUs. Single pixel wide buffers, like what search_lcu takes, would
  // be enough though.
  pixel *new_y_data = MALLOC(pixel, pic->width * pic->height);
  pixel *new_u_data = MALLOC(pixel, (pic->width * pic->height) >> 2);
  pixel *new_v_data = MALLOC(pixel, (pic->width * pic->height) >> 2);
  memcpy(new_y_data, pic->y_recdata, sizeof(pixel) * pic->width * pic->height);
  memcpy(new_u_data, pic->u_recdata, sizeof(pixel) * (pic->width * pic->height) >> 2);
  memcpy(new_v_data, pic->v_recdata, sizeof(pixel) * (pic->width * pic->height) >> 2);

  for (lcu.y = 0; lcu.y < encoder->in.height_in_lcu; lcu.y++) {
    for (lcu.x = 0; lcu.x < encoder->in.width_in_lcu; lcu.x++) {
      unsigned stride = encoder->in.width_in_lcu;
      sao_info *sao_luma = &pic->sao_luma[lcu.y * stride + lcu.x];
      sao_info *sao_chroma = &pic->sao_chroma[lcu.y * stride + lcu.x];

      // sao_do_rdo(encoder, lcu.x, lcu.y, sao_luma, sao_chroma);
      sao_reconstruct(pic, new_y_data, lcu.x, lcu.y, sao_luma, COLOR_Y);
      sao_reconstruct(pic, new_u_data, lcu.x, lcu.y, sao_chroma, COLOR_U);
      sao_reconstruct(pic, new_v_data, lcu.x, lcu.y, sao_chroma, COLOR_V);
    }
  }

  free(new_y_data);
  free(new_u_data);
  free(new_v_data);
}
