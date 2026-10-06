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

#include "lookahead.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "encoder.h"
#include "image.h"
#include "strategies/strategies-picture.h"

/** \brief Low resolution analysis block size in low resolution pixels. */
#define LA_BLOCK 8
/** \brief log2 of the number of low-res blocks in a CTU dimension. */
#define LA_CTU_SHIFT 2
/** \brief Padding around the low resolution picture. */
#define LA_PAD 48
/** \brief Maximum motion vector component in low resolution pixels. */
#define LA_SEARCH (LA_PAD - LA_BLOCK - 1)
/** \brief Size of the ring buffer for per frame complexity history. */
#define LA_HIST 1024
/** \brief Frames used in each direction for CRF complexity blurring. */
#define LA_BLUR 8
/** \brief Number of past frames used for the ABR complexity normalization. */
#define LA_ABR_PAST 96
/** \brief Lagrange multiplier for the motion vector cost in the low-res search. */
#define LA_MV_LAMBDA 4

typedef struct la_frame {
  kvz_picture *pic;
  int64_t index;

  kvz_pixel *buf;     //!< padded low resolution luma
  kvz_pixel *lowres;  //!< top-left pixel of the low resolution luma

  uint32_t *intra;    //!< intra cost of each block
  uint32_t *inter;    //!< min(inter cost, intra cost) of each block
  int16_t *mv;        //!< motion vector of each block (x, y) in low-res pixels
  float *prop;        //!< CU-tree propagate cost of each block

  uint64_t intra_sum;
  uint64_t inter_sum;

  bool has_prev;      //!< inter costs were computed against the previous frame
  bool scenecut;
} la_frame;

struct kvz_lookahead {
  const encoder_control_t *encoder;
  int depth;

  int width, height;  //!< full resolution size
  int lw, lh;         //!< low resolution size
  int bw, bh;         //!< size in blocks
  int nblocks;
  int stride;         //!< low resolution stride
  int buf_height;

  la_frame **queue;
  int count;
  int capacity;

  /** \brief The most recently released frame, needed for the next search. */
  la_frame *retired;

  int64_t next_index;

  // Scene cut state
  bool scenecut_enabled;
  double scenecut_thresh;
  int64_t keyint;
  int64_t keyint_min;
  int64_t last_key;
  int64_t segment_start;

  // CU-tree
  bool cutree;
  double cutree_strength;

  // CRF
  bool crf;
  double crf_value;
  double qcomp;
  double crf_base_cost;

  // Per frame history (ring buffer indexed by display index)
  double hist_cost[LA_HIST];
  bool hist_cut[LA_HIST];

  kvz_la_result *results;
};


static INLINE int hist_slot(int64_t index)
{
  return (int)(index % LA_HIST);
}


static void free_frame(la_frame *f)
{
  if (!f) return;
  kvz_image_free(f->pic);
  free(f->buf);
  free(f->intra);
  free(f->inter);
  free(f->mv);
  free(f->prop);
  free(f);
}


static la_frame *alloc_frame(const kvz_lookahead *la)
{
  la_frame *f = calloc(1, sizeof(la_frame));
  if (!f) return NULL;
  f->buf = MALLOC(kvz_pixel, la->stride * la->buf_height);
  f->intra = MALLOC(uint32_t, la->nblocks);
  f->inter = MALLOC(uint32_t, la->nblocks);
  f->mv = calloc(la->nblocks * 2, sizeof(int16_t));
  f->prop = MALLOC(float, la->nblocks);
  if (!f->buf || !f->intra || !f->inter || !f->mv || !f->prop) {
    free_frame(f);
    return NULL;
  }
  f->lowres = f->buf + LA_PAD * la->stride + LA_PAD;
  return f;
}


kvz_lookahead *kvz_lookahead_alloc(const encoder_control_t *encoder)
{
  const kvz_config *cfg = &encoder->cfg;
  kvz_lookahead *la = calloc(1, sizeof(kvz_lookahead));
  if (!la) return NULL;

  la->encoder = encoder;
  la->depth = MAX(1, cfg->lookahead);
  la->width = encoder->in.width;
  la->height = encoder->in.height;
  la->lw = (la->width + 1) / 2;
  la->lh = (la->height + 1) / 2;
  la->bw = CEILDIV(la->lw, LA_BLOCK);
  la->bh = CEILDIV(la->lh, LA_BLOCK);
  la->nblocks = la->bw * la->bh;
  la->stride = la->bw * LA_BLOCK + 2 * LA_PAD;
  la->buf_height = la->bh * LA_BLOCK + 2 * LA_PAD;

  la->capacity = la->depth + 4 * KVZ_MAX_GOP_LENGTH + 4;
  la->queue = calloc(la->capacity, sizeof(la_frame *));
  if (!la->queue) {
    free(la);
    return NULL;
  }

  // Scene cut detection
  la->scenecut_enabled = cfg->scenecut > 0 && cfg->intra_period != 1;
  la->scenecut_thresh = cfg->scenecut / 100.0;
  la->keyint = 0;
  if (cfg->intra_period > 0) {
    la->keyint = cfg->intra_period;
    if (cfg->gop_len > 0 && !cfg->gop_lowdelay && !cfg->open_gop) {
      // Closed GOP inserts an extra IDR after each intra period.
      la->keyint += 1;
    }
  }
  {
    const int64_t keyint_max = la->keyint > 0 ? la->keyint : 250;
    la->keyint_min = MAX(4, MIN(keyint_max / 10, 25));
  }

  la->cutree = cfg->cutree > 0;
  la->cutree_strength = cfg->cutree_strength;

  la->crf = cfg->crf >= 0;
  la->crf_value = cfg->crf;
  la->qcomp = cfg->qcomp;
  // Reference complexity per low-res block, the cost at which the frame QP
  // equals the CRF value. Empirically tuned.
  la->crf_base_cost = (cfg->gop_len > 0 && !cfg->gop_lowdelay) ? 60.0 : 40.0;

  return la;
}


void kvz_lookahead_free(kvz_lookahead *la)
{
  if (!la) return;
  for (int i = 0; i < la->count; i++) {
    free_frame(la->queue[i]);
  }
  free_frame(la->retired);
  free(la->queue);
  while (la->results) {
    kvz_la_result *next = la->results->next;
    kvz_lookahead_free_result(la->results);
    la->results = next;
  }
  free(la);
}


int kvz_lookahead_count(const kvz_lookahead *la)
{
  return la->count;
}


int kvz_lookahead_depth(const kvz_lookahead *la)
{
  return la->depth;
}


/**
 * \brief Downscale the luma of the picture by two and pad the result.
 */
static void downscale(const kvz_lookahead *la, const kvz_picture *pic, la_frame *f)
{
  const int w = MIN(pic->width, la->width);
  const int h = MIN(pic->height, la->height);
  const int stride = pic->stride;
  const kvz_pixel *src = pic->y;
  const int lw = la->lw;
  const int lh = la->lh;

  for (int y = 0; y < lh; y++) {
    const int y0 = MIN(2 * y, h - 1);
    const int y1 = MIN(2 * y + 1, h - 1);
    const kvz_pixel *r0 = &src[y0 * stride];
    const kvz_pixel *r1 = &src[y1 * stride];
    kvz_pixel *dst = &f->lowres[y * la->stride];
    for (int x = 0; x < lw; x++) {
      const int x0 = MIN(2 * x, w - 1);
      const int x1 = MIN(2 * x + 1, w - 1);
      dst[x] = (kvz_pixel)((r0[x0] + r0[x1] + r1[x0] + r1[x1] + 2) >> 2);
    }
  }

  // Pad horizontally.
  for (int y = 0; y < lh; y++) {
    kvz_pixel *row = &f->lowres[y * la->stride];
    const kvz_pixel left = row[0];
    const kvz_pixel right = row[lw - 1];
    for (int x = -LA_PAD; x < 0; x++) row[x] = left;
    for (int x = lw; x < la->stride - LA_PAD; x++) row[x] = right;
  }
  // Pad vertically.
  const kvz_pixel *first = &f->lowres[-LA_PAD];
  const kvz_pixel *last = &f->lowres[(lh - 1) * la->stride - LA_PAD];
  for (int y = -LA_PAD; y < 0; y++) {
    memcpy(&f->lowres[y * la->stride - LA_PAD], first, la->stride * sizeof(kvz_pixel));
  }
  for (int y = lh; y < la->buf_height - LA_PAD; y++) {
    memcpy(&f->lowres[y * la->stride - LA_PAD], last, la->stride * sizeof(kvz_pixel));
  }
}


/**
 * \brief Cheap intra cost of a block using DC, horizontal, vertical and
 * planar predictions from the neighbouring source pixels.
 */
static uint32_t intra_block_cost(const kvz_lookahead *la, const la_frame *f, int bx, int by)
{
  const int stride = la->stride;
  const kvz_pixel *orig = &f->lowres[by * LA_BLOCK * stride + bx * LA_BLOCK];
  const bool has_top = by > 0;
  const bool has_left = bx > 0;

  kvz_pixel top[LA_BLOCK + 1];
  kvz_pixel left[LA_BLOCK + 1];
  const kvz_pixel mid = 1 << (KVZ_BIT_DEPTH - 1);

  for (int i = 0; i < LA_BLOCK; i++) {
    top[i] = has_top ? orig[-stride + i] : (has_left ? orig[-1] : mid);
    left[i] = has_left ? orig[i * stride - 1] : (has_top ? orig[-stride] : mid);
  }
  top[LA_BLOCK] = has_top ? orig[-stride + LA_BLOCK - 1 + (bx + 1 < la->bw ? 1 : 0)] : top[LA_BLOCK - 1];
  left[LA_BLOCK] = has_left ? orig[LA_BLOCK * stride - 1 - (by + 1 < la->bh ? 0 : stride)] : left[LA_BLOCK - 1];

  kvz_pixel pred[LA_BLOCK * LA_BLOCK];
  uint32_t best = UINT32_MAX;

  // DC
  {
    int sum = 0;
    int n = 0;
    if (has_top)  { for (int i = 0; i < LA_BLOCK; i++) sum += top[i];  n += LA_BLOCK; }
    if (has_left) { for (int i = 0; i < LA_BLOCK; i++) sum += left[i]; n += LA_BLOCK; }
    const kvz_pixel dc = n ? (kvz_pixel)((sum + n / 2) / n) : mid;
    for (int i = 0; i < LA_BLOCK * LA_BLOCK; i++) pred[i] = dc;
    best = MIN(best, kvz_satd_any_size(LA_BLOCK, LA_BLOCK, pred, LA_BLOCK, orig, stride));
  }
  // Vertical
  if (has_top) {
    for (int y = 0; y < LA_BLOCK; y++) {
      for (int x = 0; x < LA_BLOCK; x++) pred[y * LA_BLOCK + x] = top[x];
    }
    best = MIN(best, kvz_satd_any_size(LA_BLOCK, LA_BLOCK, pred, LA_BLOCK, orig, stride));
  }
  // Horizontal
  if (has_left) {
    for (int y = 0; y < LA_BLOCK; y++) {
      for (int x = 0; x < LA_BLOCK; x++) pred[y * LA_BLOCK + x] = left[y];
    }
    best = MIN(best, kvz_satd_any_size(LA_BLOCK, LA_BLOCK, pred, LA_BLOCK, orig, stride));
  }
  // Planar
  if (has_top && has_left) {
    const int tr = top[LA_BLOCK];
    const int bl = left[LA_BLOCK];
    for (int y = 0; y < LA_BLOCK; y++) {
      for (int x = 0; x < LA_BLOCK; x++) {
        const int hor = (LA_BLOCK - 1 - x) * left[y] + (x + 1) * tr;
        const int ver = (LA_BLOCK - 1 - y) * top[x] + (y + 1) * bl;
        pred[y * LA_BLOCK + x] = (kvz_pixel)((hor + ver + LA_BLOCK) >> 4);
      }
    }
    best = MIN(best, kvz_satd_any_size(LA_BLOCK, LA_BLOCK, pred, LA_BLOCK, orig, stride));
  }

  // Constant for the mode signaling.
  return best + 4;
}


static INLINE int mv_bits(int v)
{
  v = abs(v);
  int bits = 1;
  while (v) {
    bits += 2;
    v >>= 1;
  }
  return bits;
}


static INLINE int median3(int a, int b, int c)
{
  return MAX(MIN(a, b), MIN(MAX(a, b), c));
}


static INLINE void clamp_mv(const kvz_lookahead *la, int bx, int by, int *mvx, int *mvy)
{
  // Keep the block inside the padded area.
  const int x = bx * LA_BLOCK;
  const int y = by * LA_BLOCK;
  const int min_x = -MIN(LA_SEARCH, x + LA_PAD - 1);
  const int min_y = -MIN(LA_SEARCH, y + LA_PAD - 1);
  const int max_x = MIN(LA_SEARCH, la->bw * LA_BLOCK + LA_PAD - LA_BLOCK - x - 1);
  const int max_y = MIN(LA_SEARCH, la->bh * LA_BLOCK + LA_PAD - LA_BLOCK - y - 1);
  *mvx = CLIP(min_x, max_x, *mvx);
  *mvy = CLIP(min_y, max_y, *mvy);
}


static INLINE uint32_t block_sad(const kvz_lookahead *la, const kvz_pixel *orig,
                                 const kvz_pixel *ref, int mvx, int mvy)
{
  return kvz_reg_sad(orig, ref + mvy * la->stride + mvx, LA_BLOCK, LA_BLOCK, la->stride, la->stride);
}


/**
 * \brief Motion search of one block against the previous frame.
 * \return inter cost (SATD + motion vector cost)
 */
static uint32_t inter_block_cost(const kvz_lookahead *la, la_frame *f, const la_frame *ref,
                                 int bx, int by)
{
  const int stride = la->stride;
  const int b = bx + by * la->bw;
  const kvz_pixel *orig = &f->lowres[by * LA_BLOCK * stride + bx * LA_BLOCK];
  const kvz_pixel *ref_blk = &ref->lowres[by * LA_BLOCK * stride + bx * LA_BLOCK];

  // Spatial predictors.
  int cand[6][2];
  int num_cand = 0;
  cand[num_cand][0] = 0; cand[num_cand][1] = 0; num_cand++;
  int pa[2] = { 0, 0 }, pb[2] = { 0, 0 }, pc[2] = { 0, 0 };
  if (bx > 0) {
    pa[0] = f->mv[2 * (b - 1)]; pa[1] = f->mv[2 * (b - 1) + 1];
    cand[num_cand][0] = pa[0]; cand[num_cand][1] = pa[1]; num_cand++;
  }
  if (by > 0) {
    pb[0] = f->mv[2 * (b - la->bw)]; pb[1] = f->mv[2 * (b - la->bw) + 1];
    cand[num_cand][0] = pb[0]; cand[num_cand][1] = pb[1]; num_cand++;
    if (bx + 1 < la->bw) {
      pc[0] = f->mv[2 * (b - la->bw + 1)]; pc[1] = f->mv[2 * (b - la->bw + 1) + 1];
      cand[num_cand][0] = pc[0]; cand[num_cand][1] = pc[1]; num_cand++;
    }
  }
  if (ref->has_prev) {
    // Temporal predictor, motion of the colocated block in the previous frame.
    cand[num_cand][0] = ref->mv[2 * b]; cand[num_cand][1] = ref->mv[2 * b + 1]; num_cand++;
  }
  const int pmv_x = median3(pa[0], pb[0], pc[0]);
  const int pmv_y = median3(pa[1], pb[1], pc[1]);

#define MV_COST(x, y) (LA_MV_LAMBDA * (mv_bits((x) - pmv_x) + mv_bits((y) - pmv_y)))

  int best_x = 0, best_y = 0;
  uint32_t best = UINT32_MAX;
  for (int i = 0; i < num_cand; i++) {
    int x = cand[i][0], y = cand[i][1];
    clamp_mv(la, bx, by, &x, &y);
    const uint32_t cost = block_sad(la, orig, ref_blk, x, y) + MV_COST(x, y);
    if (cost < best) {
      best = cost;
      best_x = x;
      best_y = y;
    }
  }

  // Hexagon search followed by a small diamond refinement.
  static const int hex[6][2] = { { -2, 0 }, { -1, 2 }, { 1, 2 }, { 2, 0 }, { 1, -2 }, { -1, -2 } };
  static const int dia[4][2] = { { -1, 0 }, { 0, 1 }, { 1, 0 }, { 0, -1 } };
  for (int iter = 0; iter < 16; iter++) {
    int cx = best_x, cy = best_y;
    for (int i = 0; i < 6; i++) {
      int x = cx + hex[i][0], y = cy + hex[i][1];
      int x_c = x, y_c = y;
      clamp_mv(la, bx, by, &x_c, &y_c);
      if (x_c != x || y_c != y) continue;
      const uint32_t cost = block_sad(la, orig, ref_blk, x, y) + MV_COST(x, y);
      if (cost < best) {
        best = cost;
        best_x = x;
        best_y = y;
      }
    }
    if (best_x == cx && best_y == cy) break;
  }
  {
    int cx = best_x, cy = best_y;
    for (int i = 0; i < 4; i++) {
      int x = cx + dia[i][0], y = cy + dia[i][1];
      int x_c = x, y_c = y;
      clamp_mv(la, bx, by, &x_c, &y_c);
      if (x_c != x || y_c != y) continue;
      const uint32_t cost = block_sad(la, orig, ref_blk, x, y) + MV_COST(x, y);
      if (cost < best) {
        best = cost;
        best_x = x;
        best_y = y;
      }
    }
  }

  f->mv[2 * b] = (int16_t)best_x;
  f->mv[2 * b + 1] = (int16_t)best_y;

  const uint32_t satd = kvz_satd_any_size(LA_BLOCK, LA_BLOCK,
                                          orig, stride,
                                          ref_blk + best_y * stride + best_x, stride);
  return satd + MV_COST(best_x, best_y);
#undef MV_COST
}


static void analyse_frame(const kvz_lookahead *la, la_frame *f, const la_frame *prev)
{
  f->intra_sum = 0;
  f->inter_sum = 0;
  f->has_prev = prev != NULL;

  for (int by = 0; by < la->bh; by++) {
    for (int bx = 0; bx < la->bw; bx++) {
      const int b = bx + by * la->bw;
      const uint32_t intra = intra_block_cost(la, f, bx, by);
      uint32_t inter = intra;
      if (prev) {
        inter = MIN(intra, inter_block_cost(la, f, prev, bx, by));
      } else {
        f->mv[2 * b] = 0;
        f->mv[2 * b + 1] = 0;
      }
      f->intra[b] = intra;
      f->inter[b] = inter;
      f->intra_sum += intra;
      f->inter_sum += inter;
    }
  }
}


/**
 * \brief Decide whether the frame starts a new scene.
 */
static void decide_scenecut(kvz_lookahead *la, la_frame *f)
{
  f->scenecut = false;

  if (f->index == 0) {
    la->last_key = 0;
    la->segment_start = 0;
    return;
  }

  const int64_t dist = f->index - la->last_key;

  if (la->scenecut_enabled && f->has_prev && dist >= 2) {
    const double thresh_max = la->scenecut_thresh;
    const double thresh_min = thresh_max * 0.25;
    const int64_t keyint_max = la->keyint > 0 ? la->keyint : 250;
    double bias;
    if (dist <= la->keyint_min / 4) {
      bias = thresh_min / 4;
    } else if (dist <= la->keyint_min) {
      bias = thresh_min * dist / la->keyint_min;
    } else {
      bias = thresh_min + (thresh_max - thresh_min) *
             (double)(dist - la->keyint_min) / MAX(1, keyint_max - la->keyint_min);
      bias = MIN(bias, thresh_max);
    }
    if ((double)f->inter_sum >= (1.0 - bias) * (double)f->intra_sum) {
      f->scenecut = true;
    }
  }

  if (f->scenecut) {
    la->last_key = f->index;
    la->segment_start = f->index;
  } else if (la->keyint > 0 && (f->index - la->segment_start) % la->keyint == 0) {
    // Periodic intra picture.
    la->last_key = f->index;
  }
}


void kvz_lookahead_push(kvz_lookahead *la, kvz_picture *pic)
{
  la_frame *f = alloc_frame(la);
  if (!f) {
    fprintf(stderr, "Lookahead: out of memory.\n");
    assert(0);
    return;
  }
  f->pic = kvz_image_copy_ref(pic);
  f->index = la->next_index++;

  downscale(la, pic, f);

  la_frame *prev = NULL;
  if (la->count > 0) {
    prev = la->queue[la->count - 1];
  } else if (la->retired && la->retired->index == f->index - 1) {
    prev = la->retired;
  }
  analyse_frame(la, f, prev);
  decide_scenecut(la, f);

  const int slot = hist_slot(f->index);
  la->hist_cost[slot] = (double)f->inter_sum / la->nblocks;
  la->hist_cut[slot] = f->scenecut || f->index == 0;
  if (f->index > 0 && !f->scenecut) {
    const int prev_slot = hist_slot(f->index - 1);
    if (la->hist_cut[prev_slot]) {
      // The P-cost of the first frame of a scene is meaningless, use the
      // temporal complexity of the following frame instead.
      la->hist_cost[prev_slot] = la->hist_cost[slot];
    }
  }

  if (la->count == la->capacity) {
    la->capacity *= 2;
    la->queue = realloc(la->queue, la->capacity * sizeof(la_frame *));
  }
  la->queue[la->count++] = f;
}


bool kvz_lookahead_next_is_scenecut(kvz_lookahead *la)
{
  return la->count > 0 && la->queue[0]->scenecut;
}


/**
 * \brief Propagate the block costs backwards through the queued frames.
 *
 * Simplified version of the x264 macroblock tree. Each frame is assumed to
 * be predicted from the previous frame in display order. The amount of
 * information a block inherits from its reference is propagated to the
 * reference blocks pointed to by the motion vector.
 */
static void cutree_propagate(kvz_lookahead *la)
{
  const int bw = la->bw;
  const int bh = la->bh;

  for (int i = 0; i < la->count; i++) {
    float *prop = la->queue[i]->prop;
    for (int b = 0; b < la->nblocks; b++) prop[b] = 0.0f;
  }

  for (int i = la->count - 1; i >= 1; i--) {
    const la_frame *f = la->queue[i];
    la_frame *ref = la->queue[i - 1];
    if (!f->has_prev || f->scenecut) continue;

    for (int by = 0; by < bh; by++) {
      for (int bx = 0; bx < bw; bx++) {
        const int b = bx + by * bw;
        const uint32_t intra = f->intra[b];
        const uint32_t inter = f->inter[b];
        if (inter >= intra) continue;

        const float amount = (intra + f->prop[b]) * (float)(intra - inter) / intra;

        const int x = bx * LA_BLOCK + f->mv[2 * b];
        const int y = by * LA_BLOCK + f->mv[2 * b + 1];
        const int x0 = x >> 3;  // floor division, LA_BLOCK == 8
        const int y0 = y >> 3;
        const int fx = x & (LA_BLOCK - 1);
        const int fy = y & (LA_BLOCK - 1);
        const float w00 = (float)((LA_BLOCK - fx) * (LA_BLOCK - fy)) / (LA_BLOCK * LA_BLOCK);
        const float w10 = (float)(fx * (LA_BLOCK - fy)) / (LA_BLOCK * LA_BLOCK);
        const float w01 = (float)((LA_BLOCK - fx) * fy) / (LA_BLOCK * LA_BLOCK);
        const float w11 = (float)(fx * fy) / (LA_BLOCK * LA_BLOCK);

        if (y0 >= 0 && y0 < bh) {
          if (x0 >= 0 && x0 < bw) ref->prop[x0 + y0 * bw] += amount * w00;
          if (x0 + 1 >= 0 && x0 + 1 < bw) ref->prop[x0 + 1 + y0 * bw] += amount * w10;
        }
        if (y0 + 1 >= 0 && y0 + 1 < bh) {
          if (x0 >= 0 && x0 < bw) ref->prop[x0 + (y0 + 1) * bw] += amount * w01;
          if (x0 + 1 >= 0 && x0 + 1 < bw) ref->prop[x0 + 1 + (y0 + 1) * bw] += amount * w11;
        }
      }
    }
  }
}


static double *cutree_ctu_offsets(const kvz_lookahead *la, const la_frame *f)
{
  const encoder_control_t *encoder = la->encoder;
  const int w_lcu = encoder->in.width_in_lcu;
  const int h_lcu = encoder->in.height_in_lcu;
  double *offsets = calloc(w_lcu * h_lcu, sizeof(double));
  if (!offsets) return NULL;

  for (int cy = 0; cy < h_lcu; cy++) {
    for (int cx = 0; cx < w_lcu; cx++) {
      double sum = 0.0;
      int n = 0;
      for (int by = cy << LA_CTU_SHIFT; by < MIN(la->bh, (cy + 1) << LA_CTU_SHIFT); by++) {
        for (int bx = cx << LA_CTU_SHIFT; bx < MIN(la->bw, (cx + 1) << LA_CTU_SHIFT); bx++) {
          const int b = bx + by * la->bw;
          const double intra = f->intra[b];
          sum += log2((intra + f->prop[b]) / intra);
          n++;
        }
      }
      offsets[cx + cy * w_lcu] = n ? -la->cutree_strength * sum / n : 0.0;
    }
  }
  return offsets;
}


/**
 * \brief Complexity of the frame blurred over the neighbouring frames of
 * the same scene.
 */
static double blurred_cost(const kvz_lookahead *la, const la_frame *f)
{
  const int64_t idx = f->index;
  double sum = la->hist_cost[hist_slot(idx)];
  double wsum = 1.0;

  // Future frames from the queue (queue[0] is f).
  double w = 1.0;
  for (int k = 1; k <= LA_BLUR && k < la->count; k++) {
    if (la->queue[k]->scenecut) break;
    w *= 0.5;
    sum += w * la->hist_cost[hist_slot(idx + k)];
    wsum += w;
  }
  // Past frames.
  w = 1.0;
  for (int k = 1; k <= LA_BLUR && idx - k >= 0 && k < LA_HIST - la->capacity; k++) {
    if (la->hist_cut[hist_slot(idx - k + 1)]) break;
    w *= 0.5;
    sum += w * la->hist_cost[hist_slot(idx - k)];
    wsum += w;
  }
  return sum / wsum;
}


kvz_picture *kvz_lookahead_pop(kvz_lookahead *la)
{
  if (la->count == 0) return NULL;

  la_frame *f = la->queue[0];

  kvz_la_result *res = calloc(1, sizeof(kvz_la_result));
  if (!res) return NULL;
  res->pic = f->pic;
  res->index = f->index;
  res->scenecut = f->scenecut;
  res->cost = la->hist_cost[hist_slot(f->index)];
  res->intra_cost = (double)f->intra_sum / la->nblocks;

  if (la->cutree) {
    cutree_propagate(la);
    res->ctu_offsets = cutree_ctu_offsets(la, f);
  }

  if (la->crf) {
    const double cost = MAX(1.0, blurred_cost(la, f));
    // qscale ~ complexity^(1 - qcomp), QP = 6 * log2(qscale) + const
    res->crf_qp = la->crf_value + 6.0 * (1.0 - la->qcomp) * log2(cost / la->crf_base_cost);
  }

  res->next = la->results;
  la->results = res;

  // Remove from the queue.
  memmove(&la->queue[0], &la->queue[1], (la->count - 1) * sizeof(la_frame *));
  la->count--;

  // Ownership of the picture moves to the caller. Keep the low resolution
  // data as it may be needed for the motion search of the next frame.
  kvz_picture *pic = f->pic;
  f->pic = NULL;
  free_frame(la->retired);
  la->retired = f;

  return pic;
}


kvz_la_result *kvz_lookahead_take_result(kvz_lookahead *la, const kvz_picture *pic)
{
  if (!la) return NULL;
  kvz_la_result **p = &la->results;
  while (*p) {
    if ((*p)->pic == pic) {
      kvz_la_result *res = *p;
      *p = res->next;
      res->next = NULL;
      return res;
    }
    p = &(*p)->next;
  }
  return NULL;
}


void kvz_lookahead_free_result(kvz_la_result *res)
{
  if (!res) return;
  free(res->ctu_offsets);
  free(res);
}


double kvz_lookahead_bit_weight(kvz_lookahead *la, int64_t first, int64_t last)
{
  if (!la || la->next_index == 0) return 1.0;

  const int64_t newest = la->next_index - 1;
  const int64_t oldest = MAX(0, newest - (LA_HIST - 1));
  first = CLIP(oldest, newest, first);
  last = CLIP(first, newest, last);

  double num = 0.0;
  for (int64_t i = first; i <= last; i++) {
    num += pow(MAX(1.0, la->hist_cost[hist_slot(i)]), la->qcomp);
  }
  num /= (double)(last - first + 1);

  const int64_t den_first = MAX(oldest, first - LA_ABR_PAST);
  double den = 0.0;
  for (int64_t i = den_first; i <= newest; i++) {
    den += pow(MAX(1.0, la->hist_cost[hist_slot(i)]), la->qcomp);
  }
  den /= (double)(newest - den_first + 1);

  if (den <= 0.0) return 1.0;
  return CLIP(0.5, 2.0, num / den);
}
