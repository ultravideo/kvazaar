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

void sao_reconstruct_color(const pixel *rec_data, pixel *new_rec_data, const sao_info *sao, 
                           int stride, int new_stride, int block_width, int block_height)
{
  int y, x;
  vector2d a_ofs = g_sao_edge_offsets[sao->eo_class][0];
  vector2d b_ofs = g_sao_edge_offsets[sao->eo_class][1];
  // Arrays orig_data and rec_data are quarter size for chroma.

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
void sao_calc_block_dims(const picture *pic, color_index color_i, 
                         const sao_info *sao, vector2d *rec, 
                         vector2d *tl, vector2d *br, vector2d *block)
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

  pixel *recdata = (color_i == COLOR_Y ? pic->y_recdata : 
                    (color_i == COLOR_U ? pic->u_recdata : pic->v_recdata));
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
  sao_calc_block_dims(pic, color_i, sao, &ofs, &tl, &br, &block);

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



void sao_search_best_mode(const pixel *data[], const pixel *recdata[], 
                          int block_width, int block_height,
                          unsigned buf_cnt,
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
  pixel orig_u[LCU_CHROMA_SIZE];
  pixel rec_u[LCU_CHROMA_SIZE];
  pixel orig_v[LCU_CHROMA_SIZE];
  pixel rec_v[LCU_CHROMA_SIZE];
  pixel *orig[2] = { orig_u, orig_v };
  pixel *rec[2] = { rec_u, rec_v };
  pixel *u_data = &pic->u_data[CU_TO_PIXEL(x_ctb, y_ctb, 1, pic->width / 2)];
  pixel *u_recdata = &pic->u_recdata[CU_TO_PIXEL(x_ctb, y_ctb, 1, pic->width / 2)];
  pixel *v_data = &pic->v_data[CU_TO_PIXEL(x_ctb, y_ctb, 1, pic->width / 2)];
  pixel *v_recdata = &pic->v_recdata[CU_TO_PIXEL(x_ctb, y_ctb, 1, pic->width / 2)];
  int block_width  = (LCU_WIDTH / 2);
  int block_height = (LCU_WIDTH / 2);

  if (x_ctb * (LCU_WIDTH / 2) + (LCU_WIDTH / 2) >= (unsigned)pic->width / 2) {
    block_width = (pic->width - x_ctb * LCU_WIDTH) / 2;
  }
  if (y_ctb * (LCU_WIDTH / 2) + (LCU_WIDTH / 2) >= (unsigned)pic->height / 2) {
    block_height = (pic->height - y_ctb * LCU_WIDTH) / 2;
  }

  sao->type = SAO_TYPE_EDGE;

  // Fill temporary buffers with picture data.
  // These buffers are needed only until we switch to a LCU based data
  // structure for pixels. Then we can give pointers directly to that structure
  // without making copies.
  picture_blit_pixels(u_data, orig_u, block_width, block_height,
                      pic->width / 2, LCU_WIDTH / 2);
  picture_blit_pixels(v_data, orig_v, block_width, block_height, 
                      pic->width / 2, LCU_WIDTH / 2);
  picture_blit_pixels(u_recdata, rec_u, block_width, block_height,
                      pic->width / 2, LCU_WIDTH / 2);
  picture_blit_pixels(v_recdata, rec_v, block_width, block_height,
                      pic->width / 2, LCU_WIDTH / 2);

  sao_search_best_mode(orig, rec, block_width, block_height, 2, sao);
}

void sao_search_luma(const picture *pic, unsigned x_ctb, unsigned y_ctb, sao_info *sao)
{
  pixel orig_y[LCU_LUMA_SIZE];
  pixel rec_y[LCU_LUMA_SIZE];
  pixel *orig[1] = { orig_y };
  pixel *rec[1] = { rec_y };
  pixel *y_data = &pic->y_data[CU_TO_PIXEL(x_ctb, y_ctb, 0, pic->width)];
  pixel *y_recdata = &pic->y_recdata[CU_TO_PIXEL(x_ctb, y_ctb, 0, pic->width)];
  int block_width = LCU_WIDTH;
  int block_height = LCU_WIDTH;

  if (x_ctb * LCU_WIDTH + LCU_WIDTH >= (unsigned)pic->width) {
    block_width = pic->width - x_ctb * LCU_WIDTH;
  }
  if (y_ctb * LCU_WIDTH + LCU_WIDTH >= (unsigned)pic->height) {
    block_height = pic->height - y_ctb * LCU_WIDTH;
  }

  sao->type = SAO_TYPE_EDGE;

  // Fill temporary buffers with picture data.
  // These buffers are needed only until we switch to a LCU based data
  // structure for pixels. Then we can give pointers directly to that structure
  // without making copies.
  picture_blit_pixels(y_data, orig_y, block_width, block_height, pic->width, LCU_WIDTH);
  picture_blit_pixels(y_recdata, rec_y, block_width, block_height, pic->width, LCU_WIDTH);

  sao_search_best_mode(orig, rec, block_width, block_height, 1, sao);
}
