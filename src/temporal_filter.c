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

#include "temporal_filter.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "image.h"
#include "threads.h"
#include "strategies/strategies-picture.h"

/** \brief Block size of the final motion vectors. */
#define TF_BLOCK 8
/** \brief Padding of the motion search planes in pixels. */
#define TF_PAD 128
/** \brief Search range at quarter resolution. */
#define TF_RANGE_QUARTER 8

#define TF_MAX_REFS 8
/** \brief Largest integer search range done with kvz_bilinear_ssd_multi. */
#define TF_MAX_SEARCH_RANGE 2
/** \brief Neighbours with a larger block error are not used for the block. */
#define TF_MAX_BLOCK_ERROR 100
/** \brief Maximum number of threads used for filtering one picture. */
#define TF_MAX_JOBS 8

typedef struct {
  kvz_pixel *buf;
  kvz_pixel *p;   //!< top-left pixel
  int w, h, stride;
} tf_plane;

typedef struct {
  int x, y;       //!< motion vector in units of 1 / factor pixels of the level
  uint64_t err;
} tf_mv;


static bool plane_alloc(tf_plane *pl, int w, int h)
{
  pl->w = w;
  pl->h = h;
  pl->stride = w + 2 * TF_PAD;
  pl->buf = MALLOC(kvz_pixel, pl->stride * (h + 2 * TF_PAD));
  pl->p = pl->buf ? pl->buf + TF_PAD * pl->stride + TF_PAD : NULL;
  return pl->buf != NULL;
}


static void plane_pad(tf_plane *pl)
{
  for (int y = 0; y < pl->h; y++) {
    kvz_pixel *row = &pl->p[y * pl->stride];
    for (int x = -TF_PAD; x < 0; x++) row[x] = row[0];
    for (int x = pl->w; x < pl->w + TF_PAD; x++) row[x] = row[pl->w - 1];
  }
  for (int y = -TF_PAD; y < 0; y++) {
    memcpy(&pl->p[y * pl->stride - TF_PAD], &pl->p[-TF_PAD], pl->stride * sizeof(kvz_pixel));
  }
  for (int y = pl->h; y < pl->h + TF_PAD; y++) {
    memcpy(&pl->p[y * pl->stride - TF_PAD], &pl->p[(pl->h - 1) * pl->stride - TF_PAD],
           pl->stride * sizeof(kvz_pixel));
  }
}


static bool plane_from_picture(tf_plane *pl, const kvz_pixel *src, int w, int h, int stride)
{
  if (!plane_alloc(pl, w, h)) return false;
  for (int y = 0; y < h; y++) {
    memcpy(&pl->p[y * pl->stride], &src[y * stride], w * sizeof(kvz_pixel));
  }
  plane_pad(pl);
  return true;
}


static bool plane_downsample(tf_plane *dst, const tf_plane *src)
{
  if (!plane_alloc(dst, src->w / 2, src->h / 2)) return false;
  for (int y = 0; y < dst->h; y++) {
    const kvz_pixel *r0 = &src->p[2 * y * src->stride];
    const kvz_pixel *r1 = r0 + src->stride;
    for (int x = 0; x < dst->w; x++) {
      dst->p[y * dst->stride + x] =
        (kvz_pixel)((r0[2 * x] + r0[2 * x + 1] + r1[2 * x] + r1[2 * x + 1] + 2) >> 2);
    }
  }
  plane_pad(dst);
  return true;
}


static void plane_free(tf_plane *pl)
{
  free(pl->buf);
  pl->buf = pl->p = NULL;
}


/**
 * \brief Error (SSE) of a block for a motion vector in units of 1 / f pixels.
 *
 * The motion vector must keep the block and one extra pixel inside the
 * padded area (see clamp_mv).
 */
static uint64_t block_error(const tf_plane *org, const tf_plane *ref, int bx, int by, int bs,
                            int mvx, int mvy, int f, uint64_t best, bool sad)
{
  // Floor division, f is 4 in all calls.
  const int ix = f == 4 ? mvx >> 2 : (mvx >= 0 ? mvx / f : -((-mvx + f - 1) / f));
  const int iy = f == 4 ? mvy >> 2 : (mvy >= 0 ? mvy / f : -((-mvy + f - 1) / f));
  const int fx = mvx - ix * f;
  const int fy = mvy - iy * f;
  const int w = MIN(bs, org->w - bx);
  const int h = MIN(bs, org->h - by);
  const int os = org->stride, rs = ref->stride;
  const kvz_pixel *o = &org->p[by * os + bx];
  const kvz_pixel *r = &ref->p[(by + iy) * rs + bx + ix];
  uint64_t err = 0;
  if (!fx && !fy && sad) {
    // Integer positions of the coarse search use the SIMD SAD.
    return kvz_reg_sad(o, r, w, h, os, rs);
  }
  if (w == bs && h == bs && (bs == 8 || bs == 16)) {
    // Full blocks use the SIMD SSD, quarter-pel positions the SIMD
    // bilinear interpolation + SSD.
    if (!fx && !fy) {
      return kvz_pixels_calc_ssd(o, r, os, rs, bs);
    }
    if (f == 4) {
      return kvz_bilinear_ssd(o, os, r, rs, bs, fx, fy);
    }
  }
  if (!fx && !fy) {
    for (int y = 0; y < h; y++) {
      for (int x = 0; x < w; x++) {
        const int diff = o[x] - r[x];
        err += diff * diff;
      }
      if (err > best) return err;
      o += os;
      r += rs;
    }
  } else {
    const int w00 = (f - fx) * (f - fy), w10 = fx * (f - fy), w01 = (f - fx) * fy, w11 = fx * fy;
    const int den = f * f;
    for (int y = 0; y < h; y++) {
      for (int x = 0; x < w; x++) {
        const int pred = (w00 * r[x] + w10 * r[x + 1] + w01 * r[x + rs] + w11 * r[x + rs + 1] + den / 2) / den;
        const int diff = o[x] - pred;
        err += diff * diff;
      }
      if (err > best) return err;
      o += os;
      r += rs;
    }
  }
  return err;
}


/** \brief Keep the motion vector (quarter pixels) inside the padded area. */
static INLINE void clamp_mv(int bs, int *x, int *y)
{
  const int lim = 4 * (TF_PAD - bs - 2);
  *x = CLIP(-lim, lim, *x);
  *y = CLIP(-lim, lim, *y);
}


/**
 * \brief Motion estimation of one resolution level.
 *
 * \param prev       motion field of the previous (coarser) level or NULL
 * \param prev_bs    block size of the previous level in its own pixels
 * \param bs         block size of this level
 * \param range      integer search range around the best candidate
 * \param subpel     search quarter pixel positions
 */
static void motion_level(const tf_plane *org, const tf_plane *ref, tf_mv *mvs, int bs,
                         const tf_mv *prev, int prev_bs, int prev_w_blocks, int prev_h_blocks,
                         int prev_factor,
                         int range, bool subpel)
{
  const int wb = CEILDIV(org->w, bs);
  const int hb = CEILDIV(org->h, bs);
  const int f = 4;  // all vectors in quarter pixels of this level

  const int full_range = range;
  for (int by = 0; by < hb; by++) {
    for (int bx = 0; bx < wb; bx++) {
      range = full_range;
      const int px = bx * bs;
      const int py = by * bs;
      tf_mv best = { 0, 0, UINT64_MAX };

      // Candidates
      tf_mv cand[8];
      int nc = 0;
      cand[nc].x = 0; cand[nc].y = 0; nc++;
      if (prev) {
        // Block of the coarser level covering this block and its neighbours.
        const int cx = (px / prev_factor) / prev_bs;
        const int cy = (py / prev_factor) / prev_bs;
        for (int dy = -1; dy <= 1; dy++) {
          for (int dx = -1; dx <= 1; dx++) {
            if (abs(dx) + abs(dy) > 1) continue;
            const int x = cx + dx, y = cy + dy;
            if (x < 0 || y < 0 || x >= prev_w_blocks || y >= prev_h_blocks) continue;
            const tf_mv *m = &prev[y * prev_w_blocks + x];
            if (m->err == UINT64_MAX) continue;
            cand[nc].x = m->x * prev_factor;
            cand[nc].y = m->y * prev_factor;
            nc++;
          }
        }
      }
      if (bx > 0) { cand[nc] = mvs[by * wb + bx - 1]; nc++; }
      if (by > 0 && nc < 8) { cand[nc] = mvs[(by - 1) * wb + bx]; nc++; }

      for (int i = 0; i < nc; i++) {
        // Round to integer pixels.
        int x = ((cand[i].x + (cand[i].x >= 0 ? 2 : -2)) / 4) * 4;
        int y = ((cand[i].y + (cand[i].y >= 0 ? 2 : -2)) / 4) * 4;
        clamp_mv(bs, &x, &y);
        const uint64_t e = block_error(org, ref, px, py, bs, x, y, f, best.err, !subpel);
        if (e < best.err) { best.x = x; best.y = y; best.err = e; }
      }

      if (!prev && bx > 0 && by > 0) {
        // Coarsest level: predictive search. The neighbours (incl. the
        // top-right one) are candidates, followed by a hexagon search with
        // a final +-1 refinement instead of the exhaustive search.
        if (bx + 1 < wb) {
          const tf_mv *tr = &mvs[(by - 1) * wb + bx + 1];
          int x = tr->x, y = tr->y;
          clamp_mv(bs, &x, &y);
          const uint64_t e = block_error(org, ref, px, py, bs, x, y, f, best.err, true);
          if (e < best.err) { best.x = x; best.y = y; best.err = e; }
        }
        static const int hex[6][2] = { { -2, 0 }, { -1, 2 }, { 1, 2 }, { 2, 0 }, { 1, -2 }, { -1, -2 } };
        for (int iter = 0; iter < 4 * range; iter++) {
          const int sx = best.x, sy = best.y;
          for (int i = 0; i < 6; i++) {
            const int x = sx + 4 * hex[i][0], y = sy + 4 * hex[i][1];
            if (abs(x) > 4 * (TF_PAD - bs - 2) || abs(y) > 4 * (TF_PAD - bs - 2)) continue;
            const uint64_t e = block_error(org, ref, px, py, bs, x, y, f, best.err, true);
            if (e < best.err) { best.x = x; best.y = y; best.err = e; }
          }
          if (best.x == sx && best.y == sy) break;
        }
        range = 1;
      }

      // Integer search around the best candidate. On the finest level the
      // errors of the positions around the best one are computed together.
      const bool multi = subpel && bs == 8 && px + bs <= org->w && py + bs <= org->h;
      const kvz_pixel *org_blk = &org->p[py * org->stride + px];
      const kvz_pixel *ref_blk = &ref->p[py * ref->stride + px];
      int pos_x[(2 * TF_MAX_SEARCH_RANGE + 1) * (2 * TF_MAX_SEARCH_RANGE + 1)];
      int pos_y[(2 * TF_MAX_SEARCH_RANGE + 1) * (2 * TF_MAX_SEARCH_RANGE + 1)];
      unsigned pos_err[(2 * TF_MAX_SEARCH_RANGE + 1) * (2 * TF_MAX_SEARCH_RANGE + 1)];
      const int cx = best.x, cy = best.y;
      if (multi && range <= TF_MAX_SEARCH_RANGE) {
        int np = 0;
        for (int dy = -range; dy <= range; dy++) {
          for (int dx = -range; dx <= range; dx++) {
            const int x = cx + 4 * dx, y = cy + 4 * dy;
            if (x == cx && y == cy) continue;
            if (abs(x) > 4 * (TF_PAD - bs - 2) || abs(y) > 4 * (TF_PAD - bs - 2)) continue;
            pos_x[np] = x;
            pos_y[np++] = y;
          }
        }
        kvz_bilinear_ssd_multi(org_blk, org->stride, ref_blk, ref->stride, bs, pos_x, pos_y, np, pos_err);
        for (int i = 0; i < np; i++) {
          if (pos_err[i] < best.err) { best.x = pos_x[i]; best.y = pos_y[i]; best.err = pos_err[i]; }
        }
      } else {
        for (int dy = -range; dy <= range; dy++) {
          for (int dx = -range; dx <= range; dx++) {
            const int x = cx + 4 * dx, y = cy + 4 * dy;
            if (x == cx && y == cy) continue;
            if (abs(x) > 4 * (TF_PAD - bs - 2) || abs(y) > 4 * (TF_PAD - bs - 2)) continue;
            const uint64_t e = block_error(org, ref, px, py, bs, x, y, f, best.err, !subpel);
            if (e < best.err) { best.x = x; best.y = y; best.err = e; }
          }
        }
      }

      if (subpel) {
        for (int step = 2; step >= 1; step >>= 1) {
          const int sx = best.x, sy = best.y;
          int np = 0;
          for (int dy = -step; dy <= step; dy += step) {
            for (int dx = -step; dx <= step; dx += step) {
              if (!dx && !dy) continue;
              if (abs(sx + dx) > 4 * (TF_PAD - bs - 2) || abs(sy + dy) > 4 * (TF_PAD - bs - 2)) continue;
              if (multi) {
                pos_x[np] = sx + dx;
                pos_y[np++] = sy + dy;
                continue;
              }
              const uint64_t e = block_error(org, ref, px, py, bs, sx + dx, sy + dy, f, best.err, false);
              if (e < best.err) { best.x = sx + dx; best.y = sy + dy; best.err = e; }
            }
          }
          if (multi) {
            kvz_bilinear_ssd_multi(org_blk, org->stride, ref_blk, ref->stride, bs, pos_x, pos_y, np, pos_err);
            for (int i = 0; i < np; i++) {
              if (pos_err[i] < best.err) { best.x = pos_x[i]; best.y = pos_y[i]; best.err = pos_err[i]; }
            }
          }
        }
      }
      mvs[by * wb + bx] = best;
    }
  }
}


/**
 * \brief Hierarchical motion estimation, returns 8x8 block vectors in quarter
 *        pixels of the full resolution.
 */
static tf_mv *motion_estimation(const tf_plane org[3], const tf_plane ref[3])
{
  const int w = org[0].w, h = org[0].h;
  // Level 2: quarter resolution, 4x4 blocks (16x16 at full resolution)
  const int bs2 = 4, bs1 = 8, bs0 = 16;
  const int wb2 = CEILDIV(org[2].w, bs2), hb2 = CEILDIV(org[2].h, bs2);
  const int wb1 = CEILDIV(org[1].w, bs1), hb1 = CEILDIV(org[1].h, bs1);
  const int wb0 = CEILDIV(w, bs0), hb0 = CEILDIV(h, bs0);
  const int wbf = CEILDIV(w, TF_BLOCK), hbf = CEILDIV(h, TF_BLOCK);

  tf_mv *m2 = MALLOC(tf_mv, wb2 * hb2);
  tf_mv *m1 = MALLOC(tf_mv, wb1 * hb1);
  tf_mv *m0 = MALLOC(tf_mv, wb0 * hb0);
  tf_mv *mf = MALLOC(tf_mv, wbf * hbf);
  if (!m2 || !m1 || !m0 || !mf) {
    free(m2); free(m1); free(m0); free(mf);
    return NULL;
  }

  motion_level(&org[2], &ref[2], m2, bs2, NULL, 0, 0, 0, 1, TF_RANGE_QUARTER, false);
  motion_level(&org[1], &ref[1], m1, bs1, m2, bs2, wb2, hb2, 2, 2, false);
  motion_level(&org[0], &ref[0], m0, bs0, m1, bs1, wb1, hb1, 2, 2, false);
  motion_level(&org[0], &ref[0], mf, TF_BLOCK, m0, bs0, wb0, hb0, 1, 1, true);

  free(m2); free(m1); free(m0);
  return mf;
}


/** \brief Shared data of the filtering of one picture. */
typedef struct {
  const kvz_picture *orig;
  const kvz_picture *const *refs;
  int w, h, num_planes, shift_w, shift_h, wbf, hbf;
  const tf_plane *org_l;
  kvz_pixel *pred[TF_MAX_REFS][3];  //!< motion compensated neighbours
  int *blk_error[TF_MAX_REFS];      //!< per 8x8 block prediction error
  int *blk_noise[TF_MAX_REFS];      //!< per 8x8 block noise estimate
  bool ok[TF_MAX_REFS];

  // Filtering
  int num_refs;
  kvz_picture *out;
  float rw_strength[TF_MAX_REFS];
  const float *lut[3];              //!< weight LUTs of each plane
  int lut_size, max_val;
} tf_ctx;


/**
 * \brief Motion search, motion compensation and block statistics of one
 *        neighbour.
 */
static void tf_process_ref(tf_ctx *t, int r)
{
  const kvz_picture *orig = t->orig;
  const kvz_picture *ref = t->refs[r];
  const int w = t->w, h = t->h, wbf = t->wbf, hbf = t->hbf;
  const int shift_w = t->shift_w, shift_h = t->shift_h;
  t->ok[r] = false;

  tf_plane ref_l[3];
  memset(ref_l, 0, sizeof(ref_l));
  bool ok = plane_from_picture(&ref_l[0], ref->y, w, h, ref->stride) &&
            plane_downsample(&ref_l[1], &ref_l[0]) &&
            plane_downsample(&ref_l[2], &ref_l[1]);
  tf_mv *mvs = ok ? motion_estimation(t->org_l, ref_l) : NULL;
  ok = ok && mvs;

  for (int c = 0; ok && c < t->num_planes; c++) {
    const int pw = c ? w >> shift_w : w;
    const int ph = c ? h >> shift_h : h;
    const int stride = c ? ref->stride_c : ref->stride;
    tf_plane rp;
    if (c == 0) {
      rp = ref_l[0];
    } else if (!plane_from_picture(&rp, ref->data[c], pw, ph, stride)) {
      ok = false;
      break;
    }
    kvz_pixel *pred = MALLOC(kvz_pixel, pw * ph);
    t->pred[r][c] = pred;
    if (!pred) {
      ok = false;
      if (c) plane_free(&rp);
      break;
    }
    // Motion vector in units of 1 / (4 << shift) pixels of this plane.
    const int fx_den = 4 << (c ? shift_w : 0);
    const int fy_den = 4 << (c ? shift_h : 0);
    // Each 8x8 luma block covers (8 >> shift) chroma pixels.
    const int bw_c = TF_BLOCK >> (c ? shift_w : 0);
    const int bh_c = TF_BLOCK >> (c ? shift_h : 0);
    const int sh_x = c ? 2 + shift_w : 2;  // log2(fx_den)
    const int sh_y = c ? 2 + shift_h : 2;
    for (int by = 0; by < ph; by += bh_c) {
      for (int bx = 0; bx < pw; bx += bw_c) {
        const int lx = c ? bx << shift_w : bx;
        const int ly = c ? by << shift_h : by;
        const tf_mv *mv = &mvs[(ly / TF_BLOCK) * wbf + lx / TF_BLOCK];
        const int ix = mv->x >> sh_x;  // floor
        const int iy = mv->y >> sh_y;
        const int fx = mv->x - ix * fx_den;
        const int fy = mv->y - iy * fy_den;
        const int w00 = (fx_den - fx) * (fy_den - fy), w10 = fx * (fy_den - fy);
        const int w01 = (fx_den - fx) * fy, w11 = fx * fy;
        const int shift = sh_x + sh_y;
        const int round = 1 << (shift - 1);
        const int x_end = MIN(pw, bx + bw_c);
        const int y_end = MIN(ph, by + bh_c);
        if (bx + ix > -TF_PAD + 1 && x_end + ix < rp.w + TF_PAD - 2 &&
            by + iy >= -TF_PAD + 1 && y_end - 1 + iy <= rp.h + TF_PAD - 2) {
          // The block and the extra column and row are inside the padding.
          kvz_tf_bilinear_pred(&pred[by * pw + bx], pw, &rp.p[(by + iy) * rp.stride + bx + ix], rp.stride,
                               x_end - bx, y_end - by, w00, w10, w01, w11, shift);
          continue;
        }
        for (int y = by; y < y_end; y++) {
          const int yy = CLIP(-TF_PAD + 1, rp.h + TF_PAD - 2, y + iy);
          const kvz_pixel *r0 = &rp.p[yy * rp.stride];
          const kvz_pixel *r1 = r0 + rp.stride;
          kvz_pixel *dst = &pred[y * pw];
          if (bx + ix > -TF_PAD + 1 && x_end + ix < rp.w + TF_PAD - 2) {
            for (int x = bx; x < x_end; x++) {
              const int xx = x + ix;
              dst[x] = (kvz_pixel)((w00 * r0[xx] + w10 * r0[xx + 1] + w01 * r1[xx] + w11 * r1[xx + 1] + round) >> shift);
            }
          } else {
            for (int x = bx; x < x_end; x++) {
              const int xx = CLIP(-TF_PAD + 1, rp.w + TF_PAD - 2, x + ix);
              dst[x] = (kvz_pixel)((w00 * r0[xx] + w10 * r0[xx + 1] + w01 * r1[xx] + w11 * r1[xx + 1] + round) >> shift);
            }
          }
        }
      }
    }
    if (c) plane_free(&rp);
  }
  free(mvs);
  for (int l = 0; l < 3; l++) plane_free(&ref_l[l]);
  if (!ok) return;

  // Block error and noise as in the VTM temporal filter.
  int *blk_error = MALLOC(int, wbf * hbf);
  int *blk_noise = MALLOC(int, wbf * hbf);
  t->blk_error[r] = blk_error;
  t->blk_noise[r] = blk_noise;
  if (!blk_error || !blk_noise) return;
  const kvz_pixel *pr = t->pred[r][0];
  for (int by = 0; by < hbf; by++) {
    for (int bx = 0; bx < wbf; bx++) {
      const int x0 = bx * TF_BLOCK, y0 = by * TF_BLOCK;
      const int bw = MIN(TF_BLOCK, w - x0), bh = MIN(TF_BLOCK, h - y0);
      int64_t st[5];
      kvz_tf_block_stats(&orig->y[y0 * orig->stride + x0], orig->stride, &pr[y0 * w + x0], w, bw, bh, st);
      const int64_t sum = st[0], sum2 = st[1], sse = st[2], var_d = st[3], diffsum = st[4];
      const int n = bw * bh;
      const double var_o = (double)sum2 - (double)sum * (double)sum / n;
      const double scale = (double)(1 << (2 * (KVZ_BIT_DEPTH - 8)));
      blk_noise[by * wbf + bx] = (int)((15.0 * 2.0 * (double)var_d / scale + 5.0) / ((double)diffsum / scale + 5.0) + 0.5);
      blk_error[by * wbf + bx] = (int)(20.0 * (((double)sse / scale + 5.0) / (var_o / scale + 5.0)) +
                                       ((double)sse / scale / n) / 50.0 + 0.5);
    }
  }
  t->ok[r] = true;
}


/**
 * \brief Filter the rows of 8x8 luma blocks [by_start, by_end) of all planes.
 */
static void tf_filter_rows(tf_ctx *t, int by_start, int by_end)
{
  const kvz_picture *orig = t->orig;
  const int wbf = t->wbf, num_refs = t->num_refs, max_val = t->max_val;
  for (int c = 0; c < t->num_planes; c++) {
    const int sx = c ? t->shift_w : 0, sy = c ? t->shift_h : 0;
    const int pw = t->w >> sx;
    const int ph = t->h >> sy;
    const kvz_pixel *src = orig->data[c];
    const int src_stride = c ? orig->stride_c : orig->stride;
    kvz_pixel *dst = t->out->data[c];
    const int dst_stride = c ? t->out->stride_c : t->out->stride;
    const float *lut = t->lut[c];
    const int lut_size = t->lut_size;

    // The weights of a reference depend on the 8x8 luma block only, so
    // they are computed once per block.
    const int bw_c = TF_BLOCK >> sx;
    const int bh_c = TF_BLOCK >> sy;
    for (int by = by_start; by < by_end; by++) {
      const int y0 = by * bh_c;
      if (y0 >= ph) break;
      for (int x0 = 0; x0 < pw; x0 += bw_c) {
        const int b = by * wbf + (x0 << sx) / TF_BLOCK;
        int min_error = INT32_MAX;
        for (int r = 0; r < num_refs; r++) {
          min_error = MIN(min_error, t->blk_error[r][b]);
        }
        float rw[TF_MAX_REFS];
        const float *rlut[TF_MAX_REFS];
        for (int r = 0; r < num_refs; r++) {
          const int noise = t->blk_noise[r][b];
          const int error = t->blk_error[r][b];
          double ww = 1.0;
          ww *= noise < 25 ? 1.0 : 1.2;
          ww *= error < 50 ? 1.2 : (error > 100 ? 0.8 : 1.0);
          ww *= (min_error + 1) / (double)(error + 1);
          // sw = (noise < 25 ? 1.3 : 0.8) * (error < 50 ? 1.3 : 1.0)
          const int sw_idx = noise < 25 ? (error < 50 ? 3 : 2) : (error < 50 ? 1 : 0);
          // A neighbour that predicts the block badly (e.g. over a fade or
          // a transition) is not used for it at all. With only the relative
          // weighting above, the block would be averaged with neighbours
          // that are all poor predictions.
          if (error > TF_MAX_BLOCK_ERROR) ww = 0.0;
          rw[r] = (float)(t->rw_strength[r] * ww);
          rlut[r] = &lut[sw_idx * lut_size + max_val];
        }
        const int y_end = MIN(ph, y0 + bh_c);
        const int x_end = MIN(pw, x0 + bw_c);
        const kvz_pixel *preds[TF_MAX_REFS];
        for (int r = 0; r < num_refs; r++) preds[r] = &t->pred[r][c][y0 * pw + x0];
        kvz_tf_filter_block(&dst[y0 * dst_stride + x0], dst_stride, &src[y0 * src_stride + x0], src_stride,
                            preds, pw, num_refs, rw, rlut, x_end - x0, y_end - y0, max_val);
      }
    }
  }
}


typedef struct {
  tf_ctx *t;
  int kind;    //!< 0: process references, 1: filter rows
  int first, last, step;
} tf_job;


static void *tf_job_run(void *arg)
{
  tf_job *j = arg;
  if (j->kind == 0) {
    for (int r = j->first; r < j->last; r += j->step) tf_process_ref(j->t, r);
  } else {
    tf_filter_rows(j->t, j->first, j->last);
  }
  return NULL;
}


/**
 * \brief Run the jobs, using up to num_jobs - 1 extra threads. The last job
 *        runs on the calling thread. Falls back to running the jobs serially
 *        if a thread can not be created.
 */
static void tf_run_jobs(tf_job *jobs, int num_jobs)
{
  pthread_t th[TF_MAX_JOBS];
  bool started[TF_MAX_JOBS];
  for (int i = 0; i < num_jobs - 1; i++) {
    started[i] = pthread_create(&th[i], NULL, tf_job_run, &jobs[i]) == 0;
    if (!started[i]) tf_job_run(&jobs[i]);
  }
  tf_job_run(&jobs[num_jobs - 1]);
  for (int i = 0; i < num_jobs - 1; i++) {
    if (started[i]) pthread_join(th[i], NULL);
  }
}


kvz_picture *kvz_temporal_filter(const kvz_picture *orig,
                                 const kvz_picture *const *refs,
                                 const int *offsets,
                                 int num_refs,
                                 int range,
                                 double strength,
                                 int qp,
                                 int threads)
{
  static const double ref_strengths[3][4] = {
    // abs(offset) 1, 2, 3, 4
    { 0.85, 0.57, 0.41, 0.33 },  // all 2 * range neighbours available
    { 1.13, 0.97, 0.81, 0.57 },  // range neighbours available
    { 0.30, 0.30, 0.30, 0.30 }   // otherwise
  };
  const double chroma_factor = 0.55;
  const double sigma_multiplier = 9.0;
  const double sigma_zero_point = 10.0;

  num_refs = MIN(num_refs, TF_MAX_REFS);
  if (num_refs <= 0 || qp <= sigma_zero_point) return NULL;
  threads = CLIP(1, TF_MAX_JOBS, threads);

  tf_ctx t;
  memset(&t, 0, sizeof(t));
  t.orig = orig;
  t.refs = refs;
  t.num_refs = num_refs;
  t.w = orig->width;
  t.h = orig->height;
  t.num_planes = orig->chroma_format == KVZ_CSP_400 ? 1 : 3;
  t.shift_w = (orig->chroma_format == KVZ_CSP_420 || KVZ_IS_422(orig->chroma_format)) ? 1 : 0;
  t.shift_h = orig->chroma_format == KVZ_CSP_420 ? 1 : 0;
  t.wbf = CEILDIV(t.w, TF_BLOCK);
  t.hbf = CEILDIV(t.h, TF_BLOCK);

  kvz_picture *out = kvz_image_alloc(orig->chroma_format, t.w, t.h);
  if (!out) return NULL;
  out->pts = orig->pts;
  out->dts = orig->dts;
  out->interlacing = orig->interlacing;
  t.out = out;

  // Planes for the motion search: full, half and quarter resolution luma.
  tf_plane org_l[3];
  memset(org_l, 0, sizeof(org_l));
  bool ok = plane_from_picture(&org_l[0], orig->y, t.w, t.h, orig->stride) &&
            plane_downsample(&org_l[1], &org_l[0]) &&
            plane_downsample(&org_l[2], &org_l[1]);
  t.org_l = org_l;

  // The neighbours are processed in parallel. The results do not depend on
  // the number of threads.
  if (ok) {
    tf_job jobs[TF_MAX_JOBS];
    const int nj = MIN(threads, num_refs);
    for (int i = 0; i < nj; i++) {
      jobs[i] = (tf_job){ &t, 0, i, num_refs, nj };
    }
    tf_run_jobs(jobs, nj);
    for (int r = 0; r < num_refs; r++) ok = ok && t.ok[r];
  }

  float *lut = NULL;
  if (ok) {
    // Strength of each neighbour.
    const int strength_idx = num_refs == 2 * range ? 0 : (num_refs == range ? 1 : 2);
    for (int r = 0; r < num_refs; r++) {
      t.rw_strength[r] = (float)ref_strengths[strength_idx][MIN(3, abs(offsets[r]) - 1)];
    }

    t.max_val = (1 << KVZ_BIT_DEPTH) - 1;
    const double bd_scale = (double)(1 << (KVZ_BIT_DEPTH - 8));
    t.lut_size = 2 * t.max_val + 1;
    lut = MALLOC(float, 5 * t.lut_size * t.num_planes);
    ok = lut != NULL;

    for (int c = 0; c < t.num_planes && lut; c++) {
      const double sigma_sq = c == 0 ?
        (qp - sigma_zero_point) * (qp - sigma_zero_point) * sigma_multiplier : 30.0 * 30.0;
      const double weight_scaling = strength * (c == 0 ? 0.4 : chroma_factor);
      // LUTs for the sigma scalings 0.8, 1.04, 1.3 and 1.69 (and 1.0 when
      // the block adaptation is not used).
      static const double sw_values[5] = { 0.8, 1.04, 1.3, 1.69, 1.0 };
      float *l = &lut[5 * t.lut_size * c];
      for (int i = 0; i < 5; i++) {
        for (int d = -t.max_val; d <= t.max_val; d++) {
          const double diff = d / bd_scale;
          l[i * t.lut_size + d + t.max_val] = (float)(weight_scaling * exp(-diff * diff / (2.0 * sw_values[i] * sigma_sq)));
        }
      }
      t.lut[c] = l;
    }
  }
  if (ok) {
    // Filter bands of block rows in parallel.
    tf_job jobs[TF_MAX_JOBS];
    const int nj = MIN(threads, t.hbf);
    for (int i = 0; i < nj; i++) {
      jobs[i] = (tf_job){ &t, 1, t.hbf * i / nj, t.hbf * (i + 1) / nj, 1 };
    }
    tf_run_jobs(jobs, nj);
  }
  free(lut);

  for (int r = 0; r < TF_MAX_REFS; r++) {
    for (int c = 0; c < 3; c++) free(t.pred[r][c]);
    free(t.blk_error[r]);
    free(t.blk_noise[r]);
  }
  for (int l = 0; l < 3; l++) plane_free(&org_l[l]);

  if (!ok) {
    kvz_image_free(out);
    return NULL;
  }
  return out;
}
