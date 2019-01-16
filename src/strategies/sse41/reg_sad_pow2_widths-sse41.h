#ifndef REG_SAD_POW2_WIDTHS_SSE41_H_
#define REG_SAD_POW2_WIDTHS_SSE41_H_

#include <immintrin.h>
#include "kvazaar.h"

static INLINE uint32_t reg_sad_w4(const kvz_pixel * const data1, const kvz_pixel * const data2,
                           const int32_t height, const uint32_t stride1,
                           const uint32_t stride2)
{
  __m128i sse_inc = _mm_setzero_si128();
  int32_t y;

  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  for (y = 0; y < height_fourline_groups; y += 4) {
    __m128i a = _mm_cvtsi32_si128(*(uint32_t *)(data1 + y * stride1));
    __m128i b = _mm_cvtsi32_si128(*(uint32_t *)(data2 + y * stride2));

    a = _mm_insert_epi32(a, *(const uint32_t *)(data1 + (y + 1) * stride1), 1);
    b = _mm_insert_epi32(b, *(const uint32_t *)(data2 + (y + 1) * stride2), 1);
    a = _mm_insert_epi32(a, *(const uint32_t *)(data1 + (y + 2) * stride1), 2);
    b = _mm_insert_epi32(b, *(const uint32_t *)(data2 + (y + 2) * stride2), 2);
    a = _mm_insert_epi32(a, *(const uint32_t *)(data1 + (y + 3) * stride1), 3);
    b = _mm_insert_epi32(b, *(const uint32_t *)(data2 + (y + 3) * stride2), 3);

    __m128i curr_sads = _mm_sad_epu8(a, b);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads);
  }
  if (height_residual_lines) {
    for (; y < height; y++) {
      __m128i a = _mm_cvtsi32_si128(*(const uint32_t *)(data1 + y * stride1));
      __m128i b = _mm_cvtsi32_si128(*(const uint32_t *)(data2 + y * stride2));

      __m128i curr_sads = _mm_sad_epu8(a, b);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads);
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);

  return _mm_cvtsi128_si32(sad);
}

static INLINE uint32_t reg_sad_w8(const kvz_pixel * const data1, const kvz_pixel * const data2,
                           const int32_t height, const uint32_t stride1,
                           const uint32_t stride2)
{
  __m128i sse_inc = _mm_setzero_si128();
  uint64_t result = 0;
  int32_t y;

  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  for (y = 0; y < height_fourline_groups; y += 4) {
    __m128d a_d = _mm_setzero_pd();
    __m128d b_d = _mm_setzero_pd();
    __m128d c_d = _mm_setzero_pd();
    __m128d d_d = _mm_setzero_pd();

    a_d = _mm_loadl_pd(a_d, (const double *)(data1 + (y + 0) * stride1));
    b_d = _mm_loadl_pd(b_d, (const double *)(data2 + (y + 0) * stride2));
    a_d = _mm_loadh_pd(a_d, (const double *)(data1 + (y + 1) * stride1));
    b_d = _mm_loadh_pd(b_d, (const double *)(data2 + (y + 1) * stride2));

    c_d = _mm_loadl_pd(c_d, (const double *)(data1 + (y + 2) * stride1));
    d_d = _mm_loadl_pd(d_d, (const double *)(data2 + (y + 2) * stride2));
    c_d = _mm_loadh_pd(c_d, (const double *)(data1 + (y + 3) * stride1));
    d_d = _mm_loadh_pd(d_d, (const double *)(data2 + (y + 3) * stride2));

    __m128i a = _mm_castpd_si128(a_d);
    __m128i b = _mm_castpd_si128(b_d);
    __m128i c = _mm_castpd_si128(c_d);
    __m128i d = _mm_castpd_si128(d_d);

    __m128i curr_sads_ab = _mm_sad_epu8(a, b);
    __m128i curr_sads_cd = _mm_sad_epu8(c, d);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
  }
  if (height_residual_lines) {
    for (; y < height; y++) {
      __m64 a = *(__m64 *)(data1 + y * stride1);
      __m64 b = *(__m64 *)(data2 + y * stride2);
      __m64 sads = _mm_sad_pu8(a, b);
      result += (uint64_t)sads;
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);

  result += _mm_cvtsi128_si32(sad);
  return result;
}

static INLINE uint32_t reg_sad_w12(const kvz_pixel * const data1, const kvz_pixel * const data2,
                            const int32_t height, const uint32_t stride1,
                            const uint32_t stride2)
{
  __m128i sse_inc = _mm_setzero_si128();
  int32_t y;
  for (y = 0; y < height; y++) {
    __m128i a = _mm_loadu_si128((const __m128i *)(data1 + y * stride1));
    __m128i b = _mm_loadu_si128((const __m128i *)(data2 + y * stride2));

    __m128i b_masked  = _mm_blend_epi16(a, b, 0x3f);
    __m128i curr_sads = _mm_sad_epu8   (a, b_masked);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads);
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);
  return _mm_cvtsi128_si32(sad);
}

static INLINE uint32_t reg_sad_w16(const kvz_pixel * const data1, const kvz_pixel * const data2,
                            const int32_t height, const uint32_t stride1,
                            const uint32_t stride2)
{
  __m128i sse_inc = _mm_setzero_si128();
  int32_t y;

  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  for (y = 0; y < height_fourline_groups; y += 4) {
    __m128i a = _mm_loadu_si128((const __m128i *)(data1 + (y + 0) * stride1));
    __m128i b = _mm_loadu_si128((const __m128i *)(data2 + (y + 0) * stride2));
    __m128i c = _mm_loadu_si128((const __m128i *)(data1 + (y + 1) * stride1));
    __m128i d = _mm_loadu_si128((const __m128i *)(data2 + (y + 1) * stride2));
    __m128i e = _mm_loadu_si128((const __m128i *)(data1 + (y + 2) * stride1));
    __m128i f = _mm_loadu_si128((const __m128i *)(data2 + (y + 2) * stride2));
    __m128i g = _mm_loadu_si128((const __m128i *)(data1 + (y + 3) * stride1));
    __m128i h = _mm_loadu_si128((const __m128i *)(data2 + (y + 3) * stride2));

    __m128i curr_sads_ab = _mm_sad_epu8(a, b);
    __m128i curr_sads_cd = _mm_sad_epu8(c, d);
    __m128i curr_sads_ef = _mm_sad_epu8(e, f);
    __m128i curr_sads_gh = _mm_sad_epu8(g, h);

    sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_ef);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_gh);
  }
  if (height_residual_lines) {
    for (; y < height; y++) {
      __m128i a = _mm_loadu_si128((const __m128i *)(data1 + (y + 0) * stride1));
      __m128i b = _mm_loadu_si128((const __m128i *)(data2 + (y + 0) * stride2));

      __m128i curr_sads = _mm_sad_epu8(a, b);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads);
    }
  }

  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);
  return _mm_cvtsi128_si32(sad);
}

static INLINE uint32_t reg_sad_w24(const kvz_pixel * const data1, const kvz_pixel * const data2,
                            const int32_t height, const uint32_t stride1,
                            const uint32_t stride2)
{
  __m128i sse_inc = _mm_setzero_si128();
  int32_t y;

  const int32_t height_doublelines = height & ~1;
  const int32_t height_parity      = height &  1;

  for (y = 0; y < height_doublelines; y += 2) {
    __m128i a = _mm_loadu_si128((const __m128i *)(data1 + (y + 0) * stride1));
    __m128i b = _mm_loadu_si128((const __m128i *)(data2 + (y + 0) * stride2));
    __m128i c = _mm_loadu_si128((const __m128i *)(data1 + (y + 1) * stride1));
    __m128i d = _mm_loadu_si128((const __m128i *)(data2 + (y + 1) * stride2));

    __m128d e_d = _mm_setzero_pd();
    __m128d f_d = _mm_setzero_pd();

    e_d = _mm_loadl_pd(e_d, (const double *)(data1 + (y + 0) * stride1 + 16));
    f_d = _mm_loadl_pd(f_d, (const double *)(data2 + (y + 0) * stride2 + 16));
    e_d = _mm_loadh_pd(e_d, (const double *)(data1 + (y + 1) * stride1 + 16));
    f_d = _mm_loadh_pd(f_d, (const double *)(data2 + (y + 1) * stride2 + 16));

    __m128i e = _mm_castpd_si128(e_d);
    __m128i f = _mm_castpd_si128(f_d);

    __m128i curr_sads_1 = _mm_sad_epu8(a, b);
    __m128i curr_sads_2 = _mm_sad_epu8(c, d);
    __m128i curr_sads_3 = _mm_sad_epu8(e, f);

    sse_inc = _mm_add_epi64(sse_inc, curr_sads_1);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_2);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_3);
  }
  if (height_parity) {
    __m128i a = _mm_loadu_si128   ((const __m128i *) (data1 + y * stride1));
    __m128i b = _mm_loadu_si128   ((const __m128i *) (data2 + y * stride2));
    __m128i c = _mm_cvtsi64_si128(*(const uint64_t *)(data1 + y * stride1 + 8));
    __m128i d = _mm_cvtsi64_si128(*(const uint64_t *)(data2 + y * stride2 + 8));

    __m128i curr_sads_1 = _mm_sad_epu8(a, b);
    __m128i curr_sads_2 = _mm_sad_epu8(c, d);

    sse_inc = _mm_add_epi64(sse_inc, curr_sads_1);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_2);
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);
  return _mm_cvtsi128_si32(sad);
}

static INLINE uint32_t reg_sad_arbitrary(const kvz_pixel * const data1, const kvz_pixel * const data2,
                                  const int32_t width, const int32_t height, const uint32_t stride1,
                                  const uint32_t stride2)
{
  int32_t y, x;
  __m128i sse_inc = _mm_setzero_si128();
  
  // Bytes in block in 128-bit blocks per each scanline, and remainder
  const int32_t width_xmms             = width  & ~15;
  const int32_t width_residual_pixels  = width  &  15;

  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  const __m128i rds    = _mm_set1_epi8 (width_residual_pixels);
  const __m128i ns     = _mm_setr_epi8 (0,  1,  2,  3,  4,  5,  6,  7,
                                        8,  9,  10, 11, 12, 13, 14, 15);
  const __m128i rdmask = _mm_cmpgt_epi8(rds, ns);

  for (x = 0; x < width_xmms; x += 16) {
    for (y = 0; y < height_fourline_groups; y += 4) {
      __m128i a = _mm_loadu_si128((const __m128i *)(data1 + (y + 0) * stride1 + x));
      __m128i b = _mm_loadu_si128((const __m128i *)(data2 + (y + 0) * stride2 + x));
      __m128i c = _mm_loadu_si128((const __m128i *)(data1 + (y + 1) * stride1 + x));
      __m128i d = _mm_loadu_si128((const __m128i *)(data2 + (y + 1) * stride2 + x));
      __m128i e = _mm_loadu_si128((const __m128i *)(data1 + (y + 2) * stride1 + x));
      __m128i f = _mm_loadu_si128((const __m128i *)(data2 + (y + 2) * stride2 + x));
      __m128i g = _mm_loadu_si128((const __m128i *)(data1 + (y + 3) * stride1 + x));
      __m128i h = _mm_loadu_si128((const __m128i *)(data2 + (y + 3) * stride2 + x));

      __m128i curr_sads_ab = _mm_sad_epu8(a, b);
      __m128i curr_sads_cd = _mm_sad_epu8(c, d);
      __m128i curr_sads_ef = _mm_sad_epu8(e, f);
      __m128i curr_sads_gh = _mm_sad_epu8(g, h);

      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ef);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_gh);
    }
    if (height_residual_lines) {
      for (; y < height; y++) {
        __m128i a = _mm_loadu_si128((const __m128i *)(data1 + y * stride1 + x));
        __m128i b = _mm_loadu_si128((const __m128i *)(data2 + y * stride2 + x));

        __m128i curr_sads = _mm_sad_epu8(a, b);

        sse_inc = _mm_add_epi64(sse_inc, curr_sads);
      }
    }
  }

  if (width_residual_pixels) {
    for (y = 0; y < height_fourline_groups; y += 4) {
      __m128i a = _mm_loadu_si128((const __m128i *)(data1 + (y + 0) * stride1 + x));
      __m128i b = _mm_loadu_si128((const __m128i *)(data2 + (y + 0) * stride2 + x));
      __m128i c = _mm_loadu_si128((const __m128i *)(data1 + (y + 1) * stride1 + x));
      __m128i d = _mm_loadu_si128((const __m128i *)(data2 + (y + 1) * stride2 + x));
      __m128i e = _mm_loadu_si128((const __m128i *)(data1 + (y + 2) * stride1 + x));
      __m128i f = _mm_loadu_si128((const __m128i *)(data2 + (y + 2) * stride2 + x));
      __m128i g = _mm_loadu_si128((const __m128i *)(data1 + (y + 3) * stride1 + x));
      __m128i h = _mm_loadu_si128((const __m128i *)(data2 + (y + 3) * stride2 + x));

      __m128i b_masked     = _mm_blendv_epi8(a, b, rdmask);
      __m128i d_masked     = _mm_blendv_epi8(c, d, rdmask);
      __m128i f_masked     = _mm_blendv_epi8(e, f, rdmask);
      __m128i h_masked     = _mm_blendv_epi8(g, h, rdmask);

      __m128i curr_sads_ab = _mm_sad_epu8   (a, b_masked);
      __m128i curr_sads_cd = _mm_sad_epu8   (c, d_masked);
      __m128i curr_sads_ef = _mm_sad_epu8   (e, f_masked);
      __m128i curr_sads_gh = _mm_sad_epu8   (g, h_masked);

      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ef);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_gh);
    }
    if (height_residual_lines) {
      for (; y < height; y++) {
        __m128i a = _mm_loadu_si128((const __m128i *)(data1 + y * stride1 + x));
        __m128i b = _mm_loadu_si128((const __m128i *)(data2 + y * stride2 + x));

        __m128i b_masked  = _mm_blendv_epi8(a, b, rdmask);
        __m128i curr_sads = _mm_sad_epu8   (a, b_masked);

        sse_inc = _mm_add_epi64(sse_inc, curr_sads);
      }
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);

  return _mm_cvtsi128_si32(sad);
}

#endif
