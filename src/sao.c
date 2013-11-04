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

#include "sao.h"

#include <string.h>

#include "picture.h"



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
 * \param orig_data  Original pixel data. 64x64 for luma, 32x32 for chroma.
 * \param rec_data  Reconstructed pixel data. 64x64 for luma, 32x32 for chroma.
 * \param dir_offsets
 * \param is_chroma  0 for luma, 1 for chroma. Indicates 
 */
void calc_sao_edge_dir(const pixel *orig_data, const pixel *rec_data,
                       int eo_class, int block_width,
                       int cat_sum_cnt[2][NUM_SAO_EDGE_CATEGORIES])
{
  int y, x;
  vector2d a_ofs = g_sao_edge_offsets[eo_class][0];
  vector2d b_ofs = g_sao_edge_offsets[eo_class][1];
  // Arrays orig_data and rec_data are quarter size for chroma.

  // Don't sample the edge pixels because this function doesn't have access to
  // their neighbours.
  for (y = 1; y < block_width - 1; ++y) {
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

void sao_reconstruct_color(const pixel *rec_data, pixel *new_rec_data, const sao_info *sao, 
                           int stride, int new_stride, int block_width, int block_height)
{
  int y, x;
  vector2d a_ofs = g_sao_edge_offsets[sao->eo_class][0];
  vector2d b_ofs = g_sao_edge_offsets[sao->eo_class][1];
  // Arrays orig_data and rec_data are quarter size for chroma.

  // Don't sample the edge pixels because this function doesn't have access to
  // their neighbours.
  for (y = 1; y < block_height - 1; ++y) {
    for (x = 1; x < block_width - 1; ++x) {
      const pixel *c_data = &rec_data[y * stride + x];
      pixel *new_data = &new_rec_data[y * new_stride + x];
      pixel a = c_data[a_ofs.y * stride + a_ofs.x];
      pixel c = c_data[0];
      pixel b = c_data[b_ofs.y * stride + b_ofs.x];
      
      int eo_idx = EO_IDX(a, b, c);
      int eo_cat = g_sao_eo_idx_to_eo_category[eo_idx];

      new_data[0] = CLIP(0, (1 << BIT_DEPTH) - 1, c_data[0] + sao->offsets[eo_cat]);
    }
  }
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
 * Vector rec is the coordinate of the area required by sao reconstruction.
 *
 * The margins are always either 0 or 1, depending on the direction of the
 * edge offset class.
 *
 * This also takes into account borders of the picture and non-LCU sized
 * CU's at the bottom and right of the picture.
 * 
 * \ rec
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
void sao_calc_block_dims(const picture *pic, const sao_info *sao, vector2d *rec, 
                         vector2d *tl, vector2d *br, vector2d *block)
{
  vector2d a_ofs = g_sao_edge_offsets[sao->eo_class][0];
  vector2d b_ofs = g_sao_edge_offsets[sao->eo_class][1];

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
  if (rec->y + LCU_WIDTH >= pic->height) {
    br->y = 0;
    if (rec->y + LCU_WIDTH >= pic->height) {
      block->y = pic->height - rec->y;
    }
    if (a_ofs.y == 1 || b_ofs.y == 1) {
      block->y -= 1;
      br->y += 1;
    }
  }
  if (rec->x + LCU_WIDTH >= pic->width) {
    br->x = 0;
    if (rec->x + LCU_WIDTH > pic->width) {
      block->x = pic->width - rec->x;
    }
    if (a_ofs.x == 1 || b_ofs.y == 1) {
      block->x -= 1;
      br->x += 1;
    }
  }

  if (rec->y != 0) {
    rec->y -= 1;
  }
  if (rec->x != 0) {
    rec->x -= 1;
  }
}

void sao_reconstruct(picture *pic, pixel *new_y_data, unsigned x_ctb, unsigned y_ctb, 
                     const sao_info *sao_luma, const sao_info *sao_chroma)
{
  pixel rec_y[(LCU_WIDTH + 2) * (LCU_WIDTH + 2)];
  pixel new_rec_y[LCU_LUMA_SIZE];
  pixel *y_recdata = &pic->y_recdata[CU_TO_PIXEL(x_ctb, y_ctb, 0, pic->width)];
  pixel *new_y_recdata = &new_y_data[CU_TO_PIXEL(x_ctb, y_ctb, 0, pic->width)];

  int x = x_ctb * LCU_WIDTH, y = y_ctb * LCU_WIDTH;
  
  vector2d rec;
  vector2d tl = { 1, 1 };
  vector2d br = { 1, 1 };
  vector2d block = { LCU_WIDTH, LCU_WIDTH };

  rec.x = x_ctb * LCU_WIDTH;
  rec.y = y_ctb * LCU_WIDTH;

  sao_calc_block_dims(pic, sao_luma, &rec, &tl, &br, &block);

  // Data to tmp buffer.
  picture_blit_pixels(&pic->y_recdata[rec.y * pic->width + rec.x], rec_y,
                      tl.x + block.x + br.x,
                      tl.y + block.y + br.y,
                      pic->width, LCU_WIDTH + 2);

  picture_blit_pixels(y_recdata, new_rec_y, LCU_WIDTH, LCU_WIDTH, pic->width, LCU_WIDTH);

  sao_reconstruct_color(&rec_y[tl.y * (tl.x + block.x + br.x) + tl.x], 
                        new_rec_y, sao_luma, 
                        LCU_WIDTH + 2, LCU_WIDTH,
                        block.x, block.y);
  //sao_reconstruct_color(rec_u, sao_chroma, COLOR_U);
  //sao_reconstruct_color(rec_v, sao_chroma, COLOR_V);
  
  // Copy reconstructed block from tmp buffer to rec image.
  picture_blit_pixels(new_rec_y, new_y_recdata, LCU_WIDTH, LCU_WIDTH, LCU_WIDTH, pic->width);
}



void sao_search_best_mode(const pixel *data, const pixel *recdata, 
                          unsigned block_width, unsigned buf_size, unsigned buf_cnt,
                          sao_info *sao_out)
{
  sao_eo_class edge_class;
  // This array is used to calculate the mean offset used to minimize distortion.
  int cat_sum_cnt[2][NUM_SAO_EDGE_CATEGORIES];
  memset(cat_sum_cnt, 0, sizeof(int) * 2 * NUM_SAO_EDGE_CATEGORIES);

  sao_out->ddistortion = INT_MAX;

  for (edge_class = SAO_EO0; edge_class <= SAO_EO3; ++edge_class) {
    int edge_offset[NUM_SAO_EDGE_CATEGORIES];
    int sum_ddistortion = 0;
    sao_eo_cat edge_cat;
    unsigned i = 0;

    // Call calc_sao_edge_dir once for luma and twice for chroma.
    for (i = 0; i < buf_cnt; ++i) {
      calc_sao_edge_dir(data + i * buf_size, recdata + i * buf_size, edge_class, block_width, cat_sum_cnt);
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
      edge_offset[edge_cat] = offset;
      // The ddistortion is amount by which the SSE of data changes. It should
      // be negative for all categories, if offset was chosen correctly.
      // ddistortion = N * h^2 - 2 * h * E, where N is the number of samples 
      // and E is the sum of errors.
      // It basically says that all pixels that are not improved by offset
      // increase increase SSE by h^2 and all pixels that are improved by
      // offset decrease SSE by h*E.
      sum_ddistortion += cat_cnt * offset * offset - 2 * offset * cat_sum;
    }
    // SAO is not applied for category 0.
    edge_offset[SAO_EO_CAT0] = 0;

    // Choose the offset class that offers the least error after offset.
    if (sum_ddistortion < sao_out->ddistortion) {
      sao_out->eo_class = edge_class;
      sao_out->ddistortion = sum_ddistortion;
      memcpy(sao_out->offsets, edge_offset, sizeof(int) * NUM_SAO_EDGE_CATEGORIES);
    }
  }
}

 void sao_search_chroma(const picture *pic, unsigned x_ctb, unsigned y_ctb, sao_info *sao)
{
  
}

void sao_search_luma(const picture *pic, unsigned x_ctb, unsigned y_ctb, sao_info *sao)
{
  // These buffers are needed only until we switch to a LCU based data
  // structure for pixels. Then we can give pointers directly to that structure
  // without making copies.
  // It's 2-dimensional because sao_search_best_mode takes arguments as arrays.
  pixel orig_y[LCU_LUMA_SIZE];
  pixel rec_y[LCU_LUMA_SIZE];
  pixel *y_data = &pic->y_data[CU_TO_PIXEL(x_ctb, y_ctb, 0, pic->width)];
  pixel *y_recdata = &pic->y_recdata[CU_TO_PIXEL(x_ctb, y_ctb, 0, pic->width)];
  
  sao->offsets[SAO_EO_CAT0] = 0;
  sao->offsets[SAO_EO_CAT1] = 7;
  sao->offsets[SAO_EO_CAT2] = 7;
  sao->offsets[SAO_EO_CAT3] = -7;
  sao->offsets[SAO_EO_CAT4] = -7;
  sao->eo_class = SAO_EO0;
  sao->type = SAO_TYPE_EDGE;
  return;

  // Fill temporary buffers with picture data.
  picture_blit_pixels(y_data, orig_y, LCU_WIDTH, LCU_WIDTH, pic->width, LCU_WIDTH);
  picture_blit_pixels(y_recdata, rec_y, LCU_WIDTH, LCU_WIDTH, pic->width, LCU_WIDTH);

  sao_search_best_mode(orig_y, rec_y, LCU_WIDTH, LCU_LUMA_SIZE, 1, sao);
}
