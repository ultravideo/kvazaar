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
#include "strategies/strategies-picture.h"

/** \brief Block size of the final motion vectors. */
#define TF_BLOCK 8
/** \brief Padding of the motion search planes in pixels. */
#define TF_PAD 128
/** \brief Search range at quarter resolution. */
#define TF_RANGE_QUARTER 8

#define TF_MAX_REFS 8

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
  const int ix = mvx >= 0 ? mvx / f : -((-mvx + f - 1) / f);
  const int iy = mvy >= 0 ? mvy / f : -((-mvy + f - 1) / f);
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

  for (int by = 0; by < hb; by++) {
    for (int bx = 0; bx < wb; bx++) {
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

      // Integer search around the best candidate.
      const int cx = best.x, cy = best.y;
      for (int dy = -range; dy <= range; dy++) {
        for (int dx = -range; dx <= range; dx++) {
          const int x = cx + 4 * dx, y = cy + 4 * dy;
          if (x == cx && y == cy) continue;
          if (abs(x) > 4 * (TF_PAD - bs - 2) || abs(y) > 4 * (TF_PAD - bs - 2)) continue;
          const uint64_t e = block_error(org, ref, px, py, bs, x, y, f, best.err, !subpel);
          if (e < best.err) { best.x = x; best.y = y; best.err = e; }
        }
      }

      if (subpel) {
        for (int step = 2; step >= 1; step >>= 1) {
          const int sx = best.x, sy = best.y;
          for (int dy = -step; dy <= step; dy += step) {
            for (int dx = -step; dx <= step; dx += step) {
              if (!dx && !dy) continue;
              if (abs(sx + dx) > 4 * (TF_PAD - bs - 2) || abs(sy + dy) > 4 * (TF_PAD - bs - 2)) continue;
              const uint64_t e = block_error(org, ref, px, py, bs, sx + dx, sy + dy, f, best.err, false);
              if (e < best.err) { best.x = sx + dx; best.y = sy + dy; best.err = e; }
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


kvz_picture *kvz_temporal_filter(const kvz_picture *orig,
                                 const kvz_picture *const *refs,
                                 const int *offsets,
                                 int num_refs,
                                 int range,
                                 double strength,
                                 int qp)
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

  const int w = orig->width, h = orig->height;
  const int num_planes = orig->chroma_format == KVZ_CSP_400 ? 1 : 3;
  const int shift_w = (orig->chroma_format == KVZ_CSP_420 || KVZ_IS_422(orig->chroma_format)) ? 1 : 0;
  const int shift_h = orig->chroma_format == KVZ_CSP_420 ? 1 : 0;

  kvz_picture *out = kvz_image_alloc(orig->chroma_format, w, h);
  if (!out) return NULL;
  out->pts = orig->pts;
  out->dts = orig->dts;
  out->interlacing = orig->interlacing;

  // Planes for the motion search: full, half and quarter resolution luma.
  tf_plane org_l[3];
  memset(org_l, 0, sizeof(org_l));
  bool ok = plane_from_picture(&org_l[0], orig->y, w, h, orig->stride) &&
            plane_downsample(&org_l[1], &org_l[0]) &&
            plane_downsample(&org_l[2], &org_l[1]);

  // Motion compensated predictions of each neighbour for each plane.
  kvz_pixel *pred[TF_MAX_REFS][3];
  memset(pred, 0, sizeof(pred));
  // Per 8x8 block prediction error and noise estimate of each neighbour.
  const int wbf = CEILDIV(w, TF_BLOCK), hbf = CEILDIV(h, TF_BLOCK);
  int *blk_error[TF_MAX_REFS];
  int *blk_noise[TF_MAX_REFS];
  memset(blk_error, 0, sizeof(blk_error));
  memset(blk_noise, 0, sizeof(blk_noise));

  for (int r = 0; ok && r < num_refs; r++) {
    const kvz_picture *ref = refs[r];
    tf_plane ref_l[3];
    memset(ref_l, 0, sizeof(ref_l));
    ok = plane_from_picture(&ref_l[0], ref->y, w, h, ref->stride) &&
         plane_downsample(&ref_l[1], &ref_l[0]) &&
         plane_downsample(&ref_l[2], &ref_l[1]);
    tf_mv *mvs = ok ? motion_estimation(org_l, ref_l) : NULL;
    ok = ok && mvs;
    const int wbf = CEILDIV(w, TF_BLOCK);

    for (int c = 0; ok && c < num_planes; c++) {
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
      pred[r][c] = MALLOC(kvz_pixel, pw * ph);
      if (!pred[r][c]) {
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
          for (int y = by; y < y_end; y++) {
            const int yy = CLIP(-TF_PAD + 1, rp.h + TF_PAD - 2, y + iy);
            const kvz_pixel *r0 = &rp.p[yy * rp.stride];
            const kvz_pixel *r1 = r0 + rp.stride;
            kvz_pixel *dst = &pred[r][c][y * pw];
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

    if (ok) {
      // Block error and noise as in the VTM temporal filter.
      blk_error[r] = MALLOC(int, wbf * hbf);
      blk_noise[r] = MALLOC(int, wbf * hbf);
      if (!blk_error[r] || !blk_noise[r]) { ok = false; break; }
      const kvz_pixel *pr = pred[r][0];
      for (int by = 0; by < hbf; by++) {
        for (int bx = 0; bx < wbf; bx++) {
          const int x0 = bx * TF_BLOCK, y0 = by * TF_BLOCK;
          const int bw = MIN(TF_BLOCK, w - x0), bh = MIN(TF_BLOCK, h - y0);
          double sum = 0, sum2 = 0, sse = 0, diffsum = 0, var_d = 0;
          for (int y = y0; y < y0 + bh; y++) {
            for (int x = x0; x < x0 + bw; x++) {
              const int o = orig->y[y * orig->stride + x];
              const int d = o - pr[y * w + x];
              sum += o;
              sum2 += o * o;
              sse += d * d;
              if (x + 1 < x0 + bw && y + 1 < y0 + bh) {
                const int dr = orig->y[y * orig->stride + x + 1] - pr[y * w + x + 1];
                const int dd = orig->y[(y + 1) * orig->stride + x] - pr[(y + 1) * w + x];
                var_d += d * d;
                diffsum += (dr - d) * (dr - d) + (dd - d) * (dd - d);
              }
            }
          }
          const int n = bw * bh;
          const double var_o = sum2 - sum * sum / n;
          const double scale = (double)(1 << (2 * (KVZ_BIT_DEPTH - 8)));
          blk_noise[r][by * wbf + bx] = (int)((15.0 * 2.0 * var_d / scale + 5.0) / (diffsum / scale + 5.0) + 0.5);
          blk_error[r][by * wbf + bx] = (int)(20.0 * ((sse / scale + 5.0) / (var_o / scale + 5.0)) +
                                              (sse / scale / n) / 50.0 + 0.5);
        }
      }
    }
  }

  if (ok) {
    // Strength of each neighbour.
    const int strength_idx = num_refs == 2 * range ? 0 : (num_refs == range ? 1 : 2);
    double ref_strength[TF_MAX_REFS];
    for (int r = 0; r < num_refs; r++) {
      ref_strength[r] = ref_strengths[strength_idx][MIN(3, abs(offsets[r]) - 1)];
    }

    const int max_val = (1 << KVZ_BIT_DEPTH) - 1;
    const double bd_scale = (double)(1 << (KVZ_BIT_DEPTH - 8));
    const int lut_size = 2 * max_val + 1;
    double *lut = MALLOC(double, 5 * lut_size);

    for (int c = 0; c < num_planes && lut; c++) {
      const int pw = c ? w >> shift_w : w;
      const int ph = c ? h >> shift_h : h;
      const kvz_pixel *src = orig->data[c];
      const int src_stride = c ? orig->stride_c : orig->stride;
      kvz_pixel *dst = out->data[c];
      const int dst_stride = c ? out->stride_c : out->stride;

      const double sigma_sq = c == 0 ?
        (qp - sigma_zero_point) * (qp - sigma_zero_point) * sigma_multiplier : 30.0 * 30.0;
      const double weight_scaling = strength * (c == 0 ? 0.4 : chroma_factor);
      // LUTs for the sigma scalings 0.8, 1.04, 1.3 and 1.69 (and 1.0 when
      // the block adaptation is not used).
      static const double sw_values[5] = { 0.8, 1.04, 1.3, 1.69, 1.0 };
      for (int i = 0; i < 5; i++) {
        for (int d = -max_val; d <= max_val; d++) {
          const double diff = d / bd_scale;
          lut[i * lut_size + d + max_val] = weight_scaling * exp(-diff * diff / (2.0 * sw_values[i] * sigma_sq));
        }
      }

      for (int y = 0; y < ph; y++) {
        for (int x = 0; x < pw; x++) {
          const int org_val = src[y * src_stride + x];
          double sum = org_val;
          double wsum = 1.0;
          const int b = ((c ? y << shift_h : y) / TF_BLOCK) * wbf + (c ? x << shift_w : x) / TF_BLOCK;
          int min_error = INT32_MAX;
          for (int r = 0; r < num_refs; r++) {
            if (blk_error[r]) min_error = MIN(min_error, blk_error[r][b]);
          }
          for (int r = 0; r < num_refs; r++) {
            const int ref_val = pred[r][c][y * pw + x];
            double ww = 1.0;
            int sw_idx = 4;
            if (blk_error[r]) {
              const int noise = blk_noise[r][b];
              const int error = blk_error[r][b];
              ww *= noise < 25 ? 1.0 : 1.2;
              ww *= error < 50 ? 1.2 : (error > 100 ? 0.8 : 1.0);
              ww *= (min_error + 1) / (double)(error + 1);
              // sw = (noise < 25 ? 1.3 : 0.8) * (error < 50 ? 1.3 : 1.0)
              sw_idx = noise < 25 ? (error < 50 ? 3 : 2) : (error < 50 ? 1 : 0);
            }
            const double wgt = ref_strength[r] * ww * lut[sw_idx * lut_size + ref_val - org_val + max_val];
            sum += wgt * ref_val;
            wsum += wgt;
          }
          dst[y * dst_stride + x] = (kvz_pixel)CLIP(0, max_val, (int)(sum / wsum + 0.5));
        }
      }
    }
    if (!lut) ok = false;
    free(lut);
  }

  for (int r = 0; r < TF_MAX_REFS; r++) {
    for (int c = 0; c < 3; c++) free(pred[r][c]);
    free(blk_error[r]);
    free(blk_noise[r]);
  }
  for (int l = 0; l < 3; l++) plane_free(&org_l[l]);

  if (!ok) {
    kvz_image_free(out);
    return NULL;
  }
  return out;
}
