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
#include "threads.h"
#include "strategies/strategies-picture.h"
#include "temporal_filter.h"

/** \brief Low resolution analysis block size in low resolution pixels. */
#define LA_BLOCK 8
/** \brief log2 of the number of low-res blocks in a CTU dimension. */
#define LA_CTU_SHIFT 2
/** \brief Padding around the low resolution picture. */
#define LA_PAD 96
/** \brief Maximum motion vector component in low resolution pixels in the
 *         search between consecutive pictures. */
#define LA_SEARCH (48 - LA_BLOCK - 1)
/** \brief Maximum motion vector component in low resolution pixels in the
 *         searches over longer distances of the adaptive GOP. */
#define LA_SEARCH_FAR (LA_PAD - LA_BLOCK - 1)
/** \brief Size of the ring buffer for per frame complexity history. */
#define LA_HIST 1024
/** \brief Frames used in each direction for CRF complexity blurring. */
#define LA_BLUR 8
/** \brief Adaptive GOP: a mini-GOP is split in two when its estimated cost
 *         is this much higher than the cost of the two halves. */
#define LA_AGOP_SPLIT_MARGIN 1.25
/** \brief Adaptive GOP: maximum length of a flash in pictures. */
#define LA_FLASH_LEN 3
/** \brief Adaptive GOP: inter / intra cost ratio below which a picture
 *         after a scene cut is predicted well from the picture before it. */
#define LA_FLASH_RATIO 0.7
/** \brief Lagrange multiplier for the motion vector cost in the low-res search. */
#define LA_MV_LAMBDA 4
/** \brief Distance of the neighbours used by the temporal filter. */
#define LA_MCTF_RANGE 2

/** \brief Number of motion fields passed to the motion search in each
 *         direction (the maximum followed reference distance). */
#define LA_MV_DIST 16

/** \brief Maximum absolute CU-tree QP offset of a CTU. */
#define LA_CUTREE_MAX_OFFSET 4.0

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
  bool segstart;      //!< starts a new coded video sequence (IDR)
  bool flash_checked; //!< adaptive GOP: checked for a flash
  bool disturbed;     //!< adaptive GOP: flash picture, coded as non-reference
  int64_t segment_start; //!< display index of the first picture of the scene

  bool mctf_done;          //!< temporal filtering has been considered
  kvz_picture *filtered;   //!< temporally filtered picture or NULL
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

  optimized_sad_func_ptr_t sad8;  //!< 8 pixel wide SAD, NULL if not available

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
  double hist_intra[LA_HIST];
  bool hist_cut[LA_HIST];
  bool hist_key[LA_HIST];

  kvz_la_result *results;

  /** \brief Motion fields (full resolution) of the latest released
   *         pictures indexed by display index % LA_MV_DIST. */
  int16_t *mv_hist[LA_MV_DIST];
  int64_t mv_hist_index[LA_MV_DIST];

  // Temporal filter
  bool mctf;
  int mctf_qp;
  /** \brief Recently analysed original pictures, [0] is the newest. Used
   *         only by the analysis (or by kvz_lookahead_pop once everything
   *         has been analysed). */
  kvz_picture *ring_pics[2 * LA_MCTF_RANGE + 1];
  int64_t ring_index[2 * LA_MCTF_RANGE + 1];
  int64_t ring_segment[2 * LA_MCTF_RANGE + 1];

  // Analysis thread. The pictures are analysed in display order by a
  // dedicated thread, or in kvz_lookahead_push when threads are disabled.
  bool threaded;
  pthread_t thread;
  pthread_mutex_t lock;
  pthread_cond_t cond;
  bool stop;
  int64_t analysed;     //!< number of analysed pictures (in display order)

  /** \brief Display index of the last picture of the window used at the
   *         latest release. Only pictures up to this one are read by the
   *         rate control so that the results do not depend on timing. */
  int64_t window_last;

  // Adaptive GOP
  bool agop;              //!< --adaptive-gop
  int agop_max;           //!< maximum mini-GOP length
  int agop_table;         //!< length of the GOP table (--gop)
  bool agop_closed;       //!< periodic key pictures start a new segment (IDR)
  int64_t agop_seg;       //!< display index of the current IDR picture
  int64_t agop_last_key;  //!< display index of the latest key picture
  /** \brief Decision for the mini-GOP starting at queue[0]. */
  struct {
    bool valid;
    bool idr;             //!< queue[0] is an IDR picture alone
    bool key;             //!< the anchor is a periodic key picture (CRA)
    int n;                //!< number of pictures
    int n_struct;         //!< size of the GOP structure (n < n_struct if truncated)
  } plan;
  /** \brief Pictures of the current mini-GOP in coding order. */
  kvz_picture *rel_pic[KVZ_MAX_GOP_LENGTH];
  int rel_count, rel_pos;
  int64_t rel_coded;      //!< pictures released in the segment
  int64_t last_dts;
  bool have_dts;
  int64_t pts_hist[LA_HIST];
  // Scratch buffers of the mini-GOP decision
  int16_t *s_mv[4];
  int16_t *s_seed;
  uint32_t *s_cost[2];
};


int kvz_lookahead_gop_template(const kvz_config *cfg, int size, int pos, int span)
{
  if (size > cfg->gop_len) {
    if (pos == size) return cfg->gop_len + 1;
    if (span >= cfg->gop_len) return cfg->gop_len;
  }
  if (pos == size) {
    for (int i = 0; i < cfg->gop_len; i++) {
      if (cfg->gop[i].poc_offset == cfg->gop_len) return i;
    }
    return 0;
  }
  if (size == cfg->gop_len) {
    for (int i = 0; i < cfg->gop_len; i++) {
      if (cfg->gop[i].poc_offset == pos) return i;
    }
  }
  for (int i = 0; i < cfg->gop_len; i++) {
    const int p = cfg->gop[i].poc_offset;
    if ((p & -p) == span) return i;
  }
  return 0;
}


static INLINE int hist_slot(int64_t index)
{
  return (int)(index % LA_HIST);
}


/**
 * \brief Temporal complexity of an analysed picture.
 *
 * The P-cost of the first picture of a scene is meaningless, the complexity
 * of the following picture is used instead when it is within the window.
 */
static double frame_cost(const kvz_lookahead *la, int64_t index, int64_t last)
{
  const int slot = hist_slot(index);
  if (la->hist_cut[slot] && index + 1 <= last && !la->hist_cut[hist_slot(index + 1)]) {
    return la->hist_cost[hist_slot(index + 1)];
  }
  return la->hist_cost[slot];
}


static void free_frame(la_frame *f)
{
  if (!f) return;
  kvz_image_free(f->pic);
  kvz_image_free(f->filtered);
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


static void *analysis_thread(void *arg);
static void agop_plan(kvz_lookahead *la);


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
  la->sad8 = kvz_get_optimized_sad(LA_BLOCK);

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
    if (cfg->gop_len > 0 && !cfg->gop_lowdelay && !cfg->open_gop && !cfg->adaptive_gop) {
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
  la->crf_base_cost = (cfg->gop_len > 0 && !cfg->gop_lowdelay) ? 260.0 : 185.0;

  la->window_last = -1;

  la->agop = cfg->adaptive_gop > 0;
  if (la->agop) {
    la->agop_max = encoder->agop_max_len;
    la->agop_table = cfg->gop_len;
    la->agop_closed = !cfg->open_gop;
    for (int i = 0; i < 4; i++) la->s_mv[i] = MALLOC(int16_t, la->nblocks * 2);
    la->s_seed = MALLOC(int16_t, la->nblocks * 2);
    la->s_cost[0] = MALLOC(uint32_t, la->nblocks);
    la->s_cost[1] = MALLOC(uint32_t, la->nblocks);
  }

  la->mctf = cfg->mctf > 0;
  if (cfg->crf >= 0) {
    la->mctf_qp = (int)(cfg->crf + 0.5);
  } else if (cfg->target_bitrate > 0) {
    la->mctf_qp = 32;
  } else {
    la->mctf_qp = cfg->qp;
  }
  for (int i = 0; i < 2 * LA_MCTF_RANGE + 1; i++) la->ring_index[i] = -1;
  pthread_mutex_init(&la->lock, NULL);
  pthread_cond_init(&la->cond, NULL);
  for (int i = 0; i < LA_MV_DIST; i++) {
    la->mv_hist_index[i] = -1;
    la->mv_hist[i] = MALLOC(int16_t, la->nblocks * 2);
    if (!la->mv_hist[i]) {
      kvz_lookahead_free(la);
      return NULL;
    }
  }
  if (cfg->threads > 0) {
    if (pthread_create(&la->thread, NULL, analysis_thread, la) != 0) {
      fprintf(stderr, "Lookahead: failed to create the analysis thread.\n");
      kvz_lookahead_free(la);
      return NULL;
    }
    la->threaded = true;
  }

  return la;
}


void kvz_lookahead_free(kvz_lookahead *la)
{
  if (!la) return;
  if (la->threaded) {
    pthread_mutex_lock(&la->lock);
    la->stop = true;
    pthread_cond_broadcast(&la->cond);
    pthread_mutex_unlock(&la->lock);
    pthread_join(la->thread, NULL);
  }
  pthread_mutex_destroy(&la->lock);
  pthread_cond_destroy(&la->cond);
  for (int i = 0; i < la->count; i++) {
    free_frame(la->queue[i]);
  }
  for (int i = 0; i < 2 * LA_MCTF_RANGE + 1; i++) {
    kvz_image_free(la->ring_pics[i]);
  }
  free_frame(la->retired);
  free(la->queue);
  for (int i = la->rel_pos; i < la->rel_count; i++) kvz_image_free(la->rel_pic[i]);
  for (int i = 0; i < 4; i++) free(la->s_mv[i]);
  free(la->s_seed);
  free(la->s_cost[0]);
  free(la->s_cost[1]);
  for (int i = 0; i < LA_MV_DIST; i++) free(la->mv_hist[i]);
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


bool kvz_lookahead_ready(const kvz_lookahead *la, bool flush)
{
  if (!la->agop) {
    return la->count > la->depth || (flush && la->count > 0);
  }
  if (la->rel_pos < la->rel_count) return true;
  // The next mini-GOP is decided when its pictures and the window of
  // the last one are available.
  return la->count > la->depth + la->agop_max || (flush && la->count > 0);
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

  // Predictions: DC, vertical, horizontal and planar.
  kvz_pixel pred[4][LA_BLOCK * LA_BLOCK];
  int8_t valid[4] = { 1, has_top, has_left, has_top && has_left };

  {
    int sum = 0;
    int n = 0;
    if (has_top)  { for (int i = 0; i < LA_BLOCK; i++) sum += top[i];  n += LA_BLOCK; }
    if (has_left) { for (int i = 0; i < LA_BLOCK; i++) sum += left[i]; n += LA_BLOCK; }
    const kvz_pixel dc = n ? (kvz_pixel)((sum + n / 2) / n) : mid;
    for (int i = 0; i < LA_BLOCK * LA_BLOCK; i++) pred[0][i] = dc;
  }
  {
    const int tr = top[LA_BLOCK];
    const int bl = left[LA_BLOCK];
    for (int y = 0; y < LA_BLOCK; y++) {
      for (int x = 0; x < LA_BLOCK; x++) {
        pred[1][y * LA_BLOCK + x] = top[x];
        pred[2][y * LA_BLOCK + x] = left[y];
        const int hor = (LA_BLOCK - 1 - x) * left[y] + (x + 1) * tr;
        const int ver = (LA_BLOCK - 1 - y) * top[x] + (y + 1) * bl;
        pred[3][y * LA_BLOCK + x] = (kvz_pixel)((hor + ver + LA_BLOCK) >> 4);
      }
    }
  }
  const kvz_pixel *preds[4] = { pred[0], pred[1], pred[2], pred[3] };
  unsigned costs[4];
  kvz_satd_any_size_quad(LA_BLOCK, LA_BLOCK, preds, LA_BLOCK, orig, stride, 4, costs, valid);
  uint32_t best = UINT32_MAX;
  for (int i = 0; i < 4; i++) {
    if (valid[i]) best = MIN(best, costs[i]);
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


static INLINE void clamp_mv(const kvz_lookahead *la, int range, int bx, int by, int *mvx, int *mvy)
{
  // Keep the block inside the padded area.
  const int x = bx * LA_BLOCK;
  const int y = by * LA_BLOCK;
  const int min_x = -MIN(range, x + LA_PAD - 1);
  const int min_y = -MIN(range, y + LA_PAD - 1);
  const int max_x = MIN(range, la->bw * LA_BLOCK + LA_PAD - LA_BLOCK - x - 1);
  const int max_y = MIN(range, la->bh * LA_BLOCK + LA_PAD - LA_BLOCK - y - 1);
  *mvx = CLIP(min_x, max_x, *mvx);
  *mvy = CLIP(min_y, max_y, *mvy);
}


static INLINE uint32_t block_sad(const kvz_lookahead *la, const kvz_pixel *orig,
                                 const kvz_pixel *ref, int mvx, int mvy)
{
  if (la->sad8) {
    return la->sad8(orig, ref + mvy * la->stride + mvx, LA_BLOCK, la->stride, la->stride);
  }
  return kvz_reg_sad(orig, ref + mvy * la->stride + mvx, LA_BLOCK, LA_BLOCK, la->stride, la->stride);
}


/**
 * \brief Motion search of one low resolution block.
 *
 * Checks the candidates, followed by a hexagon search and a small diamond
 * refinement around the best one.
 *
 * \return SATD + motion vector cost of the best motion vector
 */
static uint32_t search_block(const kvz_lookahead *la, const kvz_pixel *orig,
                             const kvz_pixel *ref_blk, int bx, int by,
                             int cand[][2], int num_cand, int pmv_x, int pmv_y,
                             int range, int *mv_x, int *mv_y)
{
  const int stride = la->stride;

#define MV_COST(x, y) (LA_MV_LAMBDA * (mv_bits((x) - pmv_x) + mv_bits((y) - pmv_y)))

  int best_x = 0, best_y = 0;
  uint32_t best = UINT32_MAX;
  for (int i = 0; i < num_cand; i++) {
    int x = cand[i][0], y = cand[i][1];
    clamp_mv(la, range, bx, by, &x, &y);
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
      clamp_mv(la, range, bx, by, &x_c, &y_c);
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
      clamp_mv(la, range, bx, by, &x_c, &y_c);
      if (x_c != x || y_c != y) continue;
      const uint32_t cost = block_sad(la, orig, ref_blk, x, y) + MV_COST(x, y);
      if (cost < best) {
        best = cost;
        best_x = x;
        best_y = y;
      }
    }
  }

  *mv_x = best_x;
  *mv_y = best_y;

  const uint32_t satd = kvz_satd_any_size(LA_BLOCK, LA_BLOCK,
                                          orig, stride,
                                          ref_blk + best_y * stride + best_x, stride);
  return satd + MV_COST(best_x, best_y);
#undef MV_COST
}


/**
 * \brief Spatial motion vector candidates and predictor of a block from the
 *        already searched blocks of a motion field.
 * \return number of candidates
 */
static int spatial_candidates(const kvz_lookahead *la, const int16_t *mvf, int bx, int by,
                              int cand[][2], int *pmv_x, int *pmv_y)
{
  const int b = bx + by * la->bw;
  int num_cand = 0;
  int pa[2] = { 0, 0 }, pb[2] = { 0, 0 }, pc[2] = { 0, 0 };
  if (bx > 0) {
    pa[0] = mvf[2 * (b - 1)]; pa[1] = mvf[2 * (b - 1) + 1];
    cand[num_cand][0] = pa[0]; cand[num_cand][1] = pa[1]; num_cand++;
  }
  if (by > 0) {
    pb[0] = mvf[2 * (b - la->bw)]; pb[1] = mvf[2 * (b - la->bw) + 1];
    cand[num_cand][0] = pb[0]; cand[num_cand][1] = pb[1]; num_cand++;
    if (bx + 1 < la->bw) {
      pc[0] = mvf[2 * (b - la->bw + 1)]; pc[1] = mvf[2 * (b - la->bw + 1) + 1];
      cand[num_cand][0] = pc[0]; cand[num_cand][1] = pc[1]; num_cand++;
    }
  }
  *pmv_x = median3(pa[0], pb[0], pc[0]);
  *pmv_y = median3(pa[1], pb[1], pc[1]);
  return num_cand;
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

  // Zero, spatial and temporal candidates.
  int cand[6][2];
  int pmv_x, pmv_y;
  cand[0][0] = 0; cand[0][1] = 0;
  int num_cand = 1 + spatial_candidates(la, f->mv, bx, by, cand + 1, &pmv_x, &pmv_y);
  if (ref->has_prev) {
    // Temporal predictor, motion of the colocated block in the previous frame.
    cand[num_cand][0] = ref->mv[2 * b]; cand[num_cand][1] = ref->mv[2 * b + 1]; num_cand++;
  }

  int mv_x, mv_y;
  const uint32_t cost = search_block(la, orig, ref_blk, bx, by, cand, num_cand, pmv_x, pmv_y,
                                     LA_SEARCH, &mv_x, &mv_y);
  f->mv[2 * b] = (int16_t)mv_x;
  f->mv[2 * b + 1] = (int16_t)mv_y;
  return cost;
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
  f->segstart = false;

  if (f->index == 0) {
    la->last_key = 0;
    la->segment_start = 0;
    f->segstart = true;
    if (la->agop) la->hist_key[hist_slot(f->index)] = true;
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
    f->segstart = true;
  } else if (la->keyint > 0 && (f->index - la->segment_start) % la->keyint == 0) {
    // Periodic intra picture. With the adaptive GOP this is only the
    // nominal position used for the scene cut threshold, the picture is
    // chosen when the mini-GOPs are decided.
    la->last_key = f->index;
  }
  // With the adaptive GOP the periodic key pictures are marked when they
  // are chosen.
  la->hist_key[hist_slot(f->index)] = la->agop ? f->scenecut : la->last_key == f->index;
  f->segment_start = la->segment_start;
}


/**
 * \brief Analyse a queued picture. Pictures are analysed in display order.
 */
static INLINE bool is_mctf_target(const la_frame *f)
{
  return (f->index - f->segment_start) % 8 == 0;
}


/**
 * \brief Temporally filter a key picture with its neighbours in the ring of
 *        recently analysed pictures.
 */
static void mctf_filter(kvz_lookahead *la, la_frame *f)
{
  f->mctf_done = true;
  const kvz_picture *refs[2 * LA_MCTF_RANGE];
  int offsets[2 * LA_MCTF_RANGE];
  int n = 0;
  for (int d = -LA_MCTF_RANGE; d <= LA_MCTF_RANGE; d++) {
    if (d == 0) continue;
    for (int i = 0; i < 2 * LA_MCTF_RANGE + 1; i++) {
      if (la->ring_pics[i] && la->ring_index[i] == f->index + d &&
          la->ring_segment[i] == f->segment_start) {
        refs[n] = la->ring_pics[i];
        offsets[n++] = d;
        break;
      }
    }
  }
  const double strength = (f->index - f->segment_start) % 16 == 0 ? 1.5 : 0.95;
  f->filtered = kvz_temporal_filter(f->pic, refs, offsets, n, LA_MCTF_RANGE, strength, la->mctf_qp,
                                    MAX(1, la->encoder->cfg.threads));
}


/**
 * \brief Temporal filtering after picture k has been analysed.
 *
 * A key picture t is filtered as soon as its future neighbours are known:
 * when t + range has been analysed or when picture k starts a new scene.
 */
static void mctf_after_analysis(kvz_lookahead *la, la_frame *k)
{
  kvz_image_free(la->ring_pics[2 * LA_MCTF_RANGE]);
  for (int i = 2 * LA_MCTF_RANGE; i > 0; i--) {
    la->ring_pics[i] = la->ring_pics[i - 1];
    la->ring_index[i] = la->ring_index[i - 1];
    la->ring_segment[i] = la->ring_segment[i - 1];
  }
  la->ring_pics[0] = kvz_image_copy_ref(k->pic);
  la->ring_index[0] = k->index;
  la->ring_segment[0] = k->segment_start;

  for (int64_t t = k->index - LA_MCTF_RANGE; t <= k->index; t++) {
    if (t < 0) continue;
    la_frame *f = NULL;
    pthread_mutex_lock(&la->lock);
    if (la->count > 0 && t >= la->queue[0]->index && t < la->queue[0]->index + la->count) {
      f = la->queue[t - la->queue[0]->index];
    }
    pthread_mutex_unlock(&la->lock);
    if (!f || f->mctf_done || !is_mctf_target(f)) continue;
    if (t == k->index - LA_MCTF_RANGE || (k->scenecut && t < k->index)) {
      mctf_filter(la, f);
    }
  }
}


static void analyse_picture(kvz_lookahead *la, la_frame *f)
{
  // The previous picture is either queued or the latest released one. It
  // has been analysed and can not be freed before this picture is analysed.
  pthread_mutex_lock(&la->lock);
  la_frame *prev = NULL;
  if (la->count > 0 && la->queue[0]->index < f->index) {
    prev = la->queue[f->index - 1 - la->queue[0]->index];
  } else if (la->retired && la->retired->index == f->index - 1) {
    prev = la->retired;
  }
  pthread_mutex_unlock(&la->lock);

  downscale(la, f->pic, f);
  analyse_frame(la, f, prev);
  decide_scenecut(la, f);

  const int slot = hist_slot(f->index);
  la->hist_cost[slot] = (double)f->inter_sum / la->nblocks;
  la->hist_intra[slot] = (double)f->intra_sum / la->nblocks;
  la->hist_cut[slot] = f->scenecut || f->index == 0;

  if (la->mctf) {
    mctf_after_analysis(la, f);
  }

  pthread_mutex_lock(&la->lock);
  la->analysed = f->index + 1;
  pthread_cond_broadcast(&la->cond);
  pthread_mutex_unlock(&la->lock);
}


static void *analysis_thread(void *arg)
{
  kvz_lookahead *la = arg;
  pthread_mutex_lock(&la->lock);
  for (;;) {
    la_frame *f = NULL;
    while (!la->stop) {
      if (la->count > 0) {
        const int64_t pos = la->analysed - la->queue[0]->index;
        if (pos >= 0 && pos < la->count) {
          f = la->queue[pos];
          break;
        }
      }
      pthread_cond_wait(&la->cond, &la->lock);
    }
    if (la->stop) break;
    pthread_mutex_unlock(&la->lock);
    analyse_picture(la, f);
    pthread_mutex_lock(&la->lock);
  }
  pthread_mutex_unlock(&la->lock);
  return NULL;
}


/** \brief Wait until the picture with the given display index is analysed. */
static void wait_analysed(kvz_lookahead *la, int64_t index)
{
  pthread_mutex_lock(&la->lock);
  while (la->analysed <= index) {
    pthread_cond_wait(&la->cond, &la->lock);
  }
  pthread_mutex_unlock(&la->lock);
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

  pthread_mutex_lock(&la->lock);
  if (la->count == la->capacity) {
    la->capacity *= 2;
    la->queue = realloc(la->queue, la->capacity * sizeof(la_frame *));
  }
  la->queue[la->count++] = f;
  pthread_cond_broadcast(&la->cond);
  pthread_mutex_unlock(&la->lock);

  if (!la->threaded) {
    analyse_picture(la, f);
  }
}


bool kvz_lookahead_next_is_scenecut(kvz_lookahead *la)
{
  if (la->agop) {
    // A new segment starts only at the start of a mini-GOP.
    if (la->rel_pos < la->rel_count || la->count == 0) return false;
    agop_plan(la);
    return la->plan.idr && la->queue[0]->index > 0;
  }
  if (la->count == 0) return false;
  wait_analysed(la, la->queue[0]->index);
  return la->queue[0]->scenecut;
}


/**
 * \brief Propagate the block costs backwards through the queued frames.
 *
 * Simplified version of the x264 macroblock tree. Each frame is assumed to
 * be predicted from the previous frame in display order. The amount of
 * information a block inherits from its reference is propagated to the
 * reference blocks pointed to by the motion vector.
 */
static void cutree_propagate(kvz_lookahead *la, int window)
{
  const int bw = la->bw;
  const int bh = la->bh;

  for (int i = 0; i < window; i++) {
    float *prop = la->queue[i]->prop;
    for (int b = 0; b < la->nblocks; b++) prop[b] = 0.0f;
  }

  for (int i = window - 1; i >= 1; i--) {
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


/**
 * \brief CU-tree QP offsets of the frame for blocks of 64 >> shift pixels.
 */
static double *cutree_offsets(const kvz_lookahead *la, const la_frame *f, int shift)
{
  const encoder_control_t *encoder = la->encoder;
  const int size = LCU_WIDTH >> shift;
  const int w_lcu = CEILDIV(encoder->in.width, size);
  const int h_lcu = CEILDIV(encoder->in.height, size);
  const int blk_shift = LA_CTU_SHIFT - shift;
  double *offsets = calloc(w_lcu * h_lcu, sizeof(double));
  if (!offsets) return NULL;

  for (int cy = 0; cy < h_lcu; cy++) {
    for (int cx = 0; cx < w_lcu; cx++) {
      double sum = 0.0;
      int n = 0;
      for (int by = cy << blk_shift; by < MIN(la->bh, (cy + 1) << blk_shift); by++) {
        for (int bx = cx << blk_shift; bx < MIN(la->bw, (cx + 1) << blk_shift); bx++) {
          const int b = bx + by * la->bw;
          const double intra = f->intra[b];
          sum += log2((intra + f->prop[b]) / intra);
          n++;
        }
      }
      offsets[cx + cy * w_lcu] = n ? -la->cutree_strength * sum / n : 0.0;
    }
  }

  // Only the distribution of the bits inside the frame is changed. The frame
  // level QP is already adapted to the temporal structure by the GOP QP
  // offsets, and applying the frame level part of the CU-tree on top of them
  // was found to hurt the coding efficiency.
  const int num_lcus = w_lcu * h_lcu;
  double mean = 0.0;
  for (int i = 0; i < num_lcus; i++) mean += offsets[i];
  mean /= num_lcus;
  for (int i = 0; i < num_lcus; i++) {
    offsets[i] = CLIP(-LA_CUTREE_MAX_OFFSET, LA_CUTREE_MAX_OFFSET, offsets[i] - mean);
  }
  return offsets;
}


/**
 * \brief Complexity of the frame blurred over the neighbouring frames of
 * the same scene.
 */
static double blurred_cost(const kvz_lookahead *la, const la_frame *f, int window)
{
  const int64_t idx = f->index;
  const int64_t last = la->window_last;
  double sum = frame_cost(la, idx, last);
  double wsum = 1.0;

  // Future frames from the window (queue[0] is f).
  double w = 1.0;
  for (int k = 1; k <= LA_BLUR && k < window; k++) {
    if (la->queue[k]->scenecut) break;
    w *= 0.5;
    sum += w * frame_cost(la, idx + k, last);
    wsum += w;
  }
  // Past frames.
  w = 1.0;
  for (int k = 1; k <= LA_BLUR && idx - k >= 0 && k < LA_HIST - la->capacity; k++) {
    if (la->hist_cut[hist_slot(idx - k + 1)]) break;
    w *= 0.5;
    sum += w * frame_cost(la, idx - k, last);
    wsum += w;
  }
  return sum / wsum;
}


/** \brief Copy the motion field of a frame in full resolution pixels. */
static void copy_motion(const kvz_lookahead *la, const la_frame *f, int16_t *dst)
{
  for (int i = 0; i < la->nblocks * 2; i++) {
    dst[i] = (int16_t)(f->mv[i] * 2);
  }
}


/**
 * \brief Collect the motion fields of the previous and following pictures
 *        of the scene for the motion search of the encoder.
 */
static void collect_motion(kvz_lookahead *la, const la_frame *f, int window, kvz_la_result *res)
{
  const int n = la->nblocks * 2;
  res->mv_bw = la->bw;
  res->mv_bh = la->bh;

  // Keep the field of this picture for the following ones.
  const int slot = (int)(f->index % LA_MV_DIST);
  if (f->has_prev) {
    copy_motion(la, f, la->mv_hist[slot]);
    la->mv_hist_index[slot] = f->index;
  } else {
    la->mv_hist_index[slot] = -1;
  }

  int num_past = 0;
  while (num_past < LA_MV_DIST &&
         la->mv_hist_index[(f->index - num_past) % LA_MV_DIST] == f->index - num_past &&
         f->index - num_past >= 0) {
    num_past++;
  }
  int num_future = 0;
  while (num_future < LA_MV_DIST && num_future + 1 < window &&
         la->queue[num_future + 1]->has_prev && !la->queue[num_future + 1]->scenecut) {
    num_future++;
  }

  if (num_past > 0) {
    res->past_mvs = MALLOC(int16_t, (size_t)n * num_past);
    if (res->past_mvs) {
      res->num_past_mvs = num_past;
      for (int k = 0; k < num_past; k++) {
        memcpy(res->past_mvs + (size_t)k * n, la->mv_hist[(f->index - k) % LA_MV_DIST],
               n * sizeof(int16_t));
      }
    }
  }
  if (num_future > 0) {
    res->future_mvs = MALLOC(int16_t, (size_t)n * num_future);
    if (res->future_mvs) {
      res->num_future_mvs = num_future;
      for (int k = 0; k < num_future; k++) {
        copy_motion(la, la->queue[k + 1], res->future_mvs + (size_t)k * n);
      }
    }
  }
}


/**
 * \brief Run the final analysis of the oldest queued picture and remove it
 *        from the queue.
 *
 * \param use_filter  code the temporally filtered picture if there is one
 * \param res_out     the analysis result of the picture
 */
static kvz_picture *pop_one(kvz_lookahead *la, bool use_filter, kvz_la_result **res_out)
{
  if (la->count == 0) return NULL;

  // The window is the oldest depth pictures. The newest picture may still be
  // under analysis while the encoder continues.
  const int window = MIN(la->count, la->depth);
  wait_analysed(la, la->queue[window - 1]->index);
  la->window_last = la->queue[window - 1]->index;

  la_frame *f = la->queue[0];
  la->pts_hist[hist_slot(f->index)] = f->pic->pts;

  kvz_la_result *res = calloc(1, sizeof(kvz_la_result));
  if (!res) return NULL;
  if (res_out) *res_out = res;
  res->pic = f->pic;
  res->index = f->index;
  res->scenecut = f->scenecut;
  res->cost = frame_cost(la, f->index, la->window_last);
  res->intra_cost = (double)f->intra_sum / la->nblocks;

  if (la->cutree) {
    cutree_propagate(la, window);
    res->ctu_offsets = cutree_offsets(la, f, 0);
    if (la->encoder->cfg.cutree_qg32) {
      res->qg_offsets = cutree_offsets(la, f, 1);
    }
  }

  res->blurred_cost = MAX(1.0, blurred_cost(la, f, window));

  // Motion fields for the motion search of the encoder.
  collect_motion(la, f, window, res);
  if (la->crf) {
    // qscale ~ complexity^(1 - qcomp), QP = 6 * log2(qscale) + const
    res->crf_qp = la->crf_value + 6.0 * (1.0 - la->qcomp) * log2(res->blurred_cost / la->crf_base_cost);
  }

  if (la->mctf && !f->mctf_done && is_mctf_target(f)) {
    // The end of the sequence was reached before the future neighbours.
    // Everything has been analysed, so the ring is not used by the thread.
    mctf_filter(la, f);
  }
  kvz_picture *filtered = f->filtered;
  f->filtered = NULL;
  if (filtered && !use_filter) {
    kvz_image_free(filtered);
    filtered = NULL;
  }

  if (filtered) {
    // The encoder codes the filtered picture, the original one is kept for
    // the output (e.g. PSNR computation).
    res->orig = f->pic;
    res->pic = filtered;
    f->pic = filtered;
  }

  res->next = la->results;
  la->results = res;

  // Remove from the queue.
  pthread_mutex_lock(&la->lock);
  memmove(&la->queue[0], &la->queue[1], (la->count - 1) * sizeof(la_frame *));
  la->count--;

  // Ownership of the picture moves to the caller. Keep the low resolution
  // data as it may be needed for the motion search of the next frame.
  kvz_picture *pic = f->pic;
  f->pic = NULL;
  la_frame *old = la->retired;
  la->retired = f;
  pthread_mutex_unlock(&la->lock);
  free_frame(old);

  return pic;
}


/** \brief Analysed frame with the given display index (queued or retired). */
static la_frame *frame_at(const kvz_lookahead *la, int64_t index)
{
  if (la->count > 0 && index >= la->queue[0]->index && index < la->queue[0]->index + la->count) {
    return la->queue[index - la->queue[0]->index];
  }
  if (la->retired && la->retired->index == index) return la->retired;
  return NULL;
}


/**
 * \brief Motion of each block of picture `from` to picture `to` (to < from)
 *        by following the motion between consecutive pictures.
 */
static void trajectory(const kvz_lookahead *la, int64_t from, int64_t to, int16_t *out)
{
  const int bw = la->bw, bh = la->bh;
  for (int by = 0; by < bh; by++) {
    for (int bx = 0; bx < bw; bx++) {
      int vx = 0, vy = 0;
      for (int64_t j = from; j > to; j--) {
        const la_frame *g = frame_at(la, j);
        if (!g || !g->has_prev) break;
        const int cx = CLIP(0, bw - 1, (bx * LA_BLOCK + vx + LA_BLOCK / 2) >> 3);
        const int cy = CLIP(0, bh - 1, (by * LA_BLOCK + vy + LA_BLOCK / 2) >> 3);
        vx += g->mv[2 * (cx + cy * bw)];
        vy += g->mv[2 * (cx + cy * bw) + 1];
      }
      const int b = bx + by * bw;
      out[2 * b] = (int16_t)CLIP(-LA_SEARCH_FAR, LA_SEARCH_FAR, vx);
      out[2 * b + 1] = (int16_t)CLIP(-LA_SEARCH_FAR, LA_SEARCH_FAR, vy);
    }
  }
}


/**
 * \brief Cost of predicting picture cur from picture ref.
 *
 * \param seed   motion vector candidate of each block
 * \param negate use the negated seed
 * \param mvf    output motion field
 * \param costs  output cost of each block (min of inter and intra), or NULL
 * \return sum of the block costs
 */
static double pair_cost(const kvz_lookahead *la, const la_frame *cur, const la_frame *ref,
                        const int16_t *seed, bool negate, int16_t *mvf, uint32_t *costs)
{
  const int stride = la->stride;
  double sum = 0.0;
  for (int by = 0; by < la->bh; by++) {
    for (int bx = 0; bx < la->bw; bx++) {
      const int b = bx + by * la->bw;
      const kvz_pixel *orig = &cur->lowres[by * LA_BLOCK * stride + bx * LA_BLOCK];
      const kvz_pixel *ref_blk = &ref->lowres[by * LA_BLOCK * stride + bx * LA_BLOCK];
      int cand[6][2];
      int pmv_x, pmv_y;
      cand[0][0] = 0; cand[0][1] = 0;
      cand[1][0] = negate ? -seed[2 * b] : seed[2 * b];
      cand[1][1] = negate ? -seed[2 * b + 1] : seed[2 * b + 1];
      const int num_cand = 2 + spatial_candidates(la, mvf, bx, by, cand + 2, &pmv_x, &pmv_y);
      int mv_x, mv_y;
      const uint32_t inter = search_block(la, orig, ref_blk, bx, by, cand, num_cand,
                                          pmv_x, pmv_y, LA_SEARCH_FAR, &mv_x, &mv_y);
      mvf[2 * b] = (int16_t)mv_x;
      mvf[2 * b + 1] = (int16_t)mv_y;
      const uint32_t c = MIN(inter, cur->intra[b]);
      if (costs) costs[b] = c;
      sum += c;
    }
  }
  return sum;
}


/**
 * \brief Cost of a bi-predicted picture from the forward and backward
 *        searches. Each block uses the cheapest of intra, forward, backward
 *        and the average of the two predictions.
 */
static double bi_cost(const kvz_lookahead *la, const la_frame *cur,
                      const la_frame *ref0, const int16_t *mv0, const uint32_t *cost0,
                      const la_frame *ref1, const int16_t *mv1, const uint32_t *cost1)
{
  const int stride = la->stride;
  double sum = 0.0;
  for (int by = 0; by < la->bh; by++) {
    for (int bx = 0; bx < la->bw; bx++) {
      const int b = bx + by * la->bw;
      const int off = by * LA_BLOCK * stride + bx * LA_BLOCK;
      const kvz_pixel *p0 = &ref0->lowres[off + mv0[2 * b + 1] * stride + mv0[2 * b]];
      const kvz_pixel *p1 = &ref1->lowres[off + mv1[2 * b + 1] * stride + mv1[2 * b]];
      kvz_pixel avg[LA_BLOCK * LA_BLOCK];
      for (int y = 0; y < LA_BLOCK; y++) {
        for (int x = 0; x < LA_BLOCK; x++) {
          avg[y * LA_BLOCK + x] = (kvz_pixel)((p0[y * stride + x] + p1[y * stride + x] + 1) >> 1);
        }
      }
      // Both motion vectors are signalled, approximated by a constant.
      const uint32_t bi = kvz_satd_any_size(LA_BLOCK, LA_BLOCK, &cur->lowres[off], stride,
                                            avg, LA_BLOCK) + 2 * LA_MV_LAMBDA * 4;
      sum += MIN(MIN(cost0[b], cost1[b]), bi);
    }
  }
  return sum;
}


/**
 * \brief Whether the mini-GOP from picture s (the previous anchor) to s + n
 *        should be split into two mini-GOPs of n / 2 pictures.
 *
 * Without the split, the anchor s + n is predicted from s and the middle
 * picture is bi-predicted from s and s + n. With the split, both s + n / 2
 * and s + n are anchors predicted over n / 2 pictures. The rest of the
 * pictures have the same structure in both cases.
 */
static bool agop_split(kvz_lookahead *la, int64_t s, int n)
{
  const int64_t m = s + n / 2;
  const int64_t e = s + n;
  const la_frame *fs = frame_at(la, s);
  const la_frame *fm = frame_at(la, m);
  const la_frame *fe = frame_at(la, e);
  if (!fs || !fm || !fe || !fm->has_prev || !fe->has_prev) return false;

  int16_t *mv_full = la->s_mv[0];
  int16_t *mv_fwd = la->s_mv[1];
  int16_t *mv_half = la->s_mv[2];
  int16_t *mv_bwd = la->s_mv[3];

  trajectory(la, e, s, la->s_seed);
  const double full = pair_cost(la, fe, fs, la->s_seed, false, mv_full, NULL);

  trajectory(la, m, s, la->s_seed);
  const double first = pair_cost(la, fm, fs, la->s_seed, false, mv_fwd, la->s_cost[0]);

  trajectory(la, e, m, la->s_seed);
  const double second = pair_cost(la, fe, fm, la->s_seed, false, mv_half, NULL);
  pair_cost(la, fm, fe, la->s_seed, true, mv_bwd, la->s_cost[1]);
  const double mid = bi_cost(la, fm, fs, mv_fwd, la->s_cost[0], fe, mv_bwd, la->s_cost[1]);

  // Longer mini-GOPs were better unless the prediction over the longer
  // distance is clearly worse, so the split needs a clear margin.
  return full + mid > LA_AGOP_SPLIT_MARGIN * (first + second);
}


/** \brief Coding order of a mini-GOP: hierarchical bisection. */
static void agop_order(kvz_la_result **res, int avail, int lo, int hi, int8_t *chain,
                       int chain_len, int *order, int *num)
{
  if (hi - lo < 2) return;
  const int mid = (lo + hi) / 2;
  if (mid <= avail) {
    kvz_la_result *r = res[mid - 1];
    r->gop_span = (hi - lo) / 2;
    // The pictures at the ends of the enclosing intervals are kept as
    // references, as in the GOP tables.
    r->gop_num_refs = 0;
    for (int i = 0; i < chain_len && r->gop_num_refs < KVZ_LA_MAX_GOP_REFS; i++) {
      r->gop_refs[r->gop_num_refs++] = (int8_t)(chain[i] - mid);
    }
    order[(*num)++] = mid;
  }
  // A non-reference picture is not a reference of the pictures below it.
  const bool ref = mid > avail || !res[mid - 1]->gop_nonref;
  if (ref) chain[chain_len] = (int8_t)mid;
  const int len = chain_len + (ref ? 1 : 0);
  agop_order(res, avail, lo, mid, chain, len, order, num);
  agop_order(res, avail, mid, hi, chain, len, order, num);
}


/**
 * \brief Adaptive GOP: whether the scene cut at picture c is a flash.
 *
 * The scene cut is a flash (or another short disturbance) if one of the
 * following pictures is predicted well from the picture before the cut.
 */
static int agop_is_flash(kvz_lookahead *la, int64_t c, int64_t last)
{
  const la_frame *before = frame_at(la, c - 1);
  if (!before) return 0;
  memset(la->s_seed, 0, la->nblocks * 2 * sizeof(int16_t));
  for (int k = 1; k <= LA_FLASH_LEN && c + k <= last; k++) {
    const la_frame *g = frame_at(la, c + k);
    if (!g) break;
    const double cost = pair_cost(la, g, before, la->s_seed, false, la->s_mv[0], NULL);
    if (cost < LA_FLASH_RATIO * (double)g->intra_sum) return k;
  }
  return 0;
}


/**
 * \brief Adaptive GOP: remove the scene cuts that are flashes.
 *
 * The pictures of a flash are marked as disturbed. They are coded as
 * non-reference pictures (except an unavoidable anchor) with the QP of the
 * highest GOP layer, and the mini-GOP anchors avoid them.
 */
static void agop_remove_flashes(kvz_lookahead *la)
{
  const int64_t last = la->queue[la->count - 1]->index;
  for (int i = 0; i < MIN(la->agop_max + 1, la->count); i++) {
    la_frame *g = la->queue[i];
    if (!g->scenecut || g->flash_checked) continue;
    g->flash_checked = true;
    const int64_t until = MIN(g->index + LA_FLASH_LEN, last);
    wait_analysed(la, until);
    const int len = agop_is_flash(la, g->index, until);
    if (len == 0) continue;
    // The pictures of the flash are coded in the scene.
    for (int k = 0; k <= LA_FLASH_LEN && i + k < la->count; k++) {
      la_frame *h = la->queue[i + k];
      if (k < len) h->disturbed = true;
      if (!h->scenecut) continue;
      h->scenecut = false;
      h->segstart = false;
      h->flash_checked = true;
      la->hist_cut[hist_slot(h->index)] = false;
      la->hist_key[hist_slot(h->index)] = false;
    }
  }
}


/**
 * \brief Decide the mini-GOP that starts with the oldest queued picture.
 *
 * The mini-GOP ends at the latest before the next segment start (scene
 * cut). Its length is halved as long as the cost model prefers two shorter
 * mini-GOPs. A periodic key picture is placed at the first anchor within
 * the last --gop / 2 pictures of the intra period so that it does not break
 * the mini-GOP structure. With a closed GOP the key picture is an IDR picture
 * that starts a new segment.
 */
static void agop_plan(kvz_lookahead *la)
{
  if (la->plan.valid) return;
  la_frame *f0 = la->queue[0];
  wait_analysed(la, MIN(f0->index + la->agop_max, la->queue[la->count - 1]->index));

  agop_remove_flashes(la);

  la->plan.valid = true;
  la->plan.n = 1;
  la->plan.n_struct = 1;
  la->plan.key = false;
  la->plan.idr = f0->segstart;
  if (la->plan.idr) return;

  const int64_t s = f0->index - 1;
  int seg_n = 0;
  for (int i = 0; i < MIN(la->agop_max, la->count); i++) {
    if (i > 0 && la->queue[i]->segstart) break;
    seg_n = i + 1;
  }
  const int64_t due = la->keyint > 0 ? la->agop_last_key + la->keyint : INT64_MAX;

  if (due - s > seg_n && seg_n < la->agop_table) {
    // The segment ends within a mini-GOP of the GOP table length. The
    // pictures are coded with the structure of a full mini-GOP from which
    // the missing pictures are left out, as with the GOP tables. Nothing
    // refers to them after the segment, so no anchor is needed.
    la->plan.n_struct = la->agop_table;
    la->plan.n = seg_n;
    return;
  }

  const int max_n = (int)MIN(seg_n, due - s);
  int n = 1;
  while (n * 2 <= max_n) n *= 2;
  while (n > 1 && agop_split(la, s, n)) n /= 2;
  if (la->queue[n - 1]->disturbed) {
    // The anchor is a reference of the following mini-GOPs, so a flash
    // picture is avoided: the nearest length with another anchor.
    int best = n;
    for (int m = n / 2; m >= 1 && best == n; m /= 2) {
      if (!la->queue[m - 1]->disturbed) best = m;
    }
    for (int m = n * 2; m <= max_n && best == n; m *= 2) {
      if (!la->queue[m - 1]->disturbed) best = m;
    }
    n = best;
  }
  la->plan.n = n;
  la->plan.n_struct = n;

  const int64_t e = s + n;
  if (la->keyint > 0 && e > due - la->agop_table / 2) {
    if (la->agop_closed) {
      // The key picture starts a new segment, the pictures before it form
      // a truncated mini-GOP.
      la->queue[n - 1]->segstart = true;
      if (n == 1) {
        la->plan.idr = true;
      } else {
        la->plan.n = n - 1;
        la->plan.n_struct = n;
      }
    } else {
      la->plan.key = true;
    }
  }
}


/**
 * \brief Analyse the pictures of the planned mini-GOP and put them in the
 *        release list in coding order.
 */
static void agop_next(kvz_lookahead *la)
{
  agop_plan(la);
  const int n = la->plan.n;
  const int n_struct = la->plan.n_struct;
  const bool key = la->plan.key;
  const bool idr = la->plan.idr;
  la->plan.valid = false;

  kvz_picture *pics[KVZ_MAX_GOP_LENGTH];
  kvz_la_result *res[KVZ_MAX_GOP_LENGTH];
  bool nonref[KVZ_MAX_GOP_LENGTH];
  for (int i = 0; i < n; i++) {
    const int pos = i + 1;
    // Flash pictures other than the anchor are not used as references.
    nonref[i] = la->queue[0]->disturbed && !(pos == n && n == n_struct);
    // Only the pictures high in the hierarchy are temporally filtered.
    const bool use_filter = !nonref[i] && ((pos & -pos) >= 4 || (pos == n && n == n_struct));
    res[i] = NULL;
    pics[i] = pop_one(la, use_filter, &res[i]);
    assert(pics[i] && res[i]);
  }

  if (idr) {
    la->agop_seg = res[0]->index;
    la->agop_last_key = res[0]->index;
    la->rel_coded = 0;
  }
  if (key) {
    la->agop_last_key = res[n - 1]->index;
  }
  if (idr || key) {
    la->hist_key[hist_slot(res[n - 1]->index)] = true;
  }
  for (int i = 0; i < n; i++) {
    res[i]->gop_size = n_struct;
    res[i]->gop_pos = i + 1;
    res[i]->gop_span = n_struct;
    res[i]->gop_key = false;
    res[i]->gop_poc = (int32_t)(res[i]->index - la->agop_seg);
    res[i]->gop_nonref = nonref[i];
  }
  res[n - 1]->gop_key = key;

  int order[KVZ_MAX_GOP_LENGTH];
  int num = 0;
  if (n == n_struct) order[num++] = n;
  int8_t chain[8] = { 0, (int8_t)n_struct };
  agop_order(res, n, 0, n_struct, chain, 2, order, &num);
  assert(num == n);

  // Decoding timestamps: the decoding order is delayed by the maximum
  // reordering of the stream.
  const int reorder = la->agop_max - 1;
  int64_t dur = 1;
  if (res[n - 1]->index > la->agop_seg) {
    dur = MAX(1, la->pts_hist[hist_slot(la->agop_seg + 1)] - la->pts_hist[hist_slot(la->agop_seg)]);
  }
  for (int k = 0; k < n; k++) {
    kvz_picture *pic = pics[order[k] - 1];
    const int64_t c = la->rel_coded++;
    int64_t dts;
    if (c >= reorder) {
      dts = la->pts_hist[hist_slot(la->agop_seg + c - reorder)] - 1;
    } else {
      dts = la->pts_hist[hist_slot(la->agop_seg)] - (reorder - c) * dur;
    }
    if (la->have_dts && dts <= la->last_dts) dts = la->last_dts + 1;
    la->last_dts = dts;
    la->have_dts = true;
    pic->dts = dts;
    la->rel_pic[k] = pic;
  }
  la->rel_count = n;
  la->rel_pos = 0;
}


kvz_picture *kvz_lookahead_pop(kvz_lookahead *la)
{
  if (!la->agop) {
    return pop_one(la, true, NULL);
  }
  if (la->rel_pos == la->rel_count) {
    if (la->count == 0) return NULL;
    agop_next(la);
  }
  kvz_picture *pic = la->rel_pic[la->rel_pos];
  la->rel_pic[la->rel_pos++] = NULL;
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
  free(res->past_mvs);
  free(res->future_mvs);
  kvz_image_free(res->orig);
  free(res->qg_offsets);
  free(res->ctu_offsets);
  free(res);
}


int64_t kvz_lookahead_newest(const kvz_lookahead *la)
{
  return la ? la->window_last : -1;
}


bool kvz_lookahead_frame_info(const kvz_lookahead *la, int64_t index,
                              double *cost, double *intra, bool *key)
{
  if (!la || index < 0 || index > la->window_last || index <= la->window_last - LA_HIST + la->capacity) {
    return false;
  }
  const int slot = hist_slot(index);
  *cost = MAX(1.0, frame_cost(la, index, la->window_last));
  *intra = MAX(1.0, la->hist_intra[slot]);
  *key = la->hist_key[slot];
  return true;
}

