#ifndef REG_SAD_POW2_WIDTHS_SSE41_H_
#define REG_SAD_POW2_WIDTHS_SSE41_H_

#include <immintrin.h>
#include "kvazaar.h"

static INLINE uint32_t reg_sad_w0(const kvz_pixel * const data1, const kvz_pixel * const data2,
                           const int32_t height, const uint32_t stride1,
                           const uint32_t stride2)
{
  return 0;
}

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
    __m128i c = _mm_cvtsi64_si128(*(const uint64_t *)(data1 + y * stride1 + 16));
    __m128i d = _mm_cvtsi64_si128(*(const uint64_t *)(data2 + y * stride2 + 16));

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

static uint32_t ver_sad_w4(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                           int32_t height, uint32_t stride)
{
  __m128i ref_row = _mm_set1_epi32(*(const uint32_t *)ref_data);
  __m128i sse_inc = _mm_setzero_si128();
  int32_t y;

  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  for (y = 0; y < height_fourline_groups; y += 4) {
    __m128i a = _mm_cvtsi32_si128(*(uint32_t *)(pic_data + y * stride));

    a = _mm_insert_epi32(a, *(const uint32_t *)(pic_data + (y + 1) * stride), 1);
    a = _mm_insert_epi32(a, *(const uint32_t *)(pic_data + (y + 2) * stride), 2);
    a = _mm_insert_epi32(a, *(const uint32_t *)(pic_data + (y + 3) * stride), 3);

    __m128i curr_sads = _mm_sad_epu8(a, ref_row);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads);
  }
  if (height_residual_lines) {
    // Only pick the last dword, because we're comparing single dwords (lines)
    ref_row = _mm_bsrli_si128(ref_row, 12);

    for (; y < height; y++) {
      __m128i a = _mm_cvtsi32_si128(*(const uint32_t *)(pic_data + y * stride));

      __m128i curr_sads = _mm_sad_epu8(a, ref_row);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads);
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);

  return _mm_cvtsi128_si32(sad);
}

static uint32_t ver_sad_w8(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                           int32_t height, uint32_t stride)
{
  const __m128i ref_row = _mm_set1_epi64x(*(const uint64_t *)ref_data);
  __m128i sse_inc = _mm_setzero_si128();
  uint64_t result = 0;
  int32_t y;

  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  for (y = 0; y < height_fourline_groups; y += 4) {
    __m128d a_d = _mm_setzero_pd();
    __m128d c_d = _mm_setzero_pd();

    a_d = _mm_loadl_pd(a_d, (const double *)(pic_data + (y + 0) * stride));
    a_d = _mm_loadh_pd(a_d, (const double *)(pic_data + (y + 1) * stride));

    c_d = _mm_loadl_pd(c_d, (const double *)(pic_data + (y + 2) * stride));
    c_d = _mm_loadh_pd(c_d, (const double *)(pic_data + (y + 3) * stride));

    __m128i a = _mm_castpd_si128(a_d);
    __m128i c = _mm_castpd_si128(c_d);

    __m128i curr_sads_ab = _mm_sad_epu8(a, ref_row);
    __m128i curr_sads_cd = _mm_sad_epu8(c, ref_row);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
  }
  if (height_residual_lines) {
    __m64 b = (__m64)_mm_cvtsi128_si64(ref_row);

    for (; y < height; y++) {
      __m64 a = *(__m64 *)(pic_data + y * stride);
      __m64 sads = _mm_sad_pu8(a, b);
      result += (uint64_t)sads;
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);

  result += _mm_cvtsi128_si32(sad);
  return result;
}

static uint32_t ver_sad_w12(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                            int32_t height, uint32_t stride)
{
  const __m128i ref_row = _mm_loadu_si128((__m128i *)ref_data);
  __m128i sse_inc = _mm_setzero_si128();
  int32_t y;

  for (y = 0; y < height; y++) {
    __m128i a = _mm_loadu_si128((const __m128i *)(pic_data + y * stride));

    __m128i a_masked  = _mm_blend_epi16(ref_row, a, 0x3f);
    __m128i curr_sads = _mm_sad_epu8   (ref_row, a_masked);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads);
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);
  return _mm_cvtsi128_si32(sad);
}

static uint32_t ver_sad_w16(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                            int32_t height, uint32_t stride)
{
  const __m128i ref_row = _mm_loadu_si128((__m128i *)ref_data);
  __m128i sse_inc       = _mm_setzero_si128();
  int32_t y;

  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  for (y = 0; y < height_fourline_groups; y += 4) {
    __m128i pic_row_1   = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * stride));
    __m128i pic_row_2   = _mm_loadu_si128((__m128i *)(pic_data + (y + 1) * stride));
    __m128i pic_row_3   = _mm_loadu_si128((__m128i *)(pic_data + (y + 2) * stride));
    __m128i pic_row_4   = _mm_loadu_si128((__m128i *)(pic_data + (y + 3) * stride));

    __m128i curr_sads_1 = _mm_sad_epu8   (pic_row_1, ref_row);
    __m128i curr_sads_2 = _mm_sad_epu8   (pic_row_2, ref_row);
    __m128i curr_sads_3 = _mm_sad_epu8   (pic_row_3, ref_row);
    __m128i curr_sads_4 = _mm_sad_epu8   (pic_row_4, ref_row);

    sse_inc = _mm_add_epi64(sse_inc, curr_sads_1);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_2);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_3);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_4);
  }
  if (height_residual_lines) {
    for (; y < height; y++) {
      __m128i pic_row   = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * stride));
      __m128i curr_sads = _mm_sad_epu8   (pic_row, ref_row);

      sse_inc = _mm_add_epi64(sse_inc, curr_sads);
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);

  return _mm_cvtsi128_si32(sad);
}

static uint32_t ver_sad_arbitrary(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                  int32_t width, int32_t height, uint32_t stride)
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
    const __m128i ref_row = _mm_loadu_si128((__m128i *)(ref_data + x));
    for (y = 0; y < height_fourline_groups; y += 4) {
      __m128i a = _mm_loadu_si128((const __m128i *)(pic_data + (y + 0) * stride + x));
      __m128i c = _mm_loadu_si128((const __m128i *)(pic_data + (y + 1) * stride + x));
      __m128i e = _mm_loadu_si128((const __m128i *)(pic_data + (y + 2) * stride + x));
      __m128i g = _mm_loadu_si128((const __m128i *)(pic_data + (y + 3) * stride + x));

      __m128i curr_sads_ab = _mm_sad_epu8(ref_row, a);
      __m128i curr_sads_cd = _mm_sad_epu8(ref_row, c);
      __m128i curr_sads_ef = _mm_sad_epu8(ref_row, e);
      __m128i curr_sads_gh = _mm_sad_epu8(ref_row, g);

      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ef);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_gh);
    }
    if (height_residual_lines) {
      for (; y < height; y++) {
        __m128i a = _mm_loadu_si128((const __m128i *)(pic_data + y * stride + x));

        __m128i curr_sads = _mm_sad_epu8(a, ref_row);

        sse_inc = _mm_add_epi64(sse_inc, curr_sads);
      }
    }
  }

  if (width_residual_pixels) {
    const __m128i ref_row = _mm_loadu_si128((__m128i *)(ref_data + x));
    for (y = 0; y < height_fourline_groups; y += 4) {
      __m128i a = _mm_loadu_si128((const __m128i *)(pic_data + (y + 0) * stride + x));
      __m128i c = _mm_loadu_si128((const __m128i *)(pic_data + (y + 1) * stride + x));
      __m128i e = _mm_loadu_si128((const __m128i *)(pic_data + (y + 2) * stride + x));
      __m128i g = _mm_loadu_si128((const __m128i *)(pic_data + (y + 3) * stride + x));

      __m128i a_masked     = _mm_blendv_epi8(ref_row, a, rdmask);
      __m128i c_masked     = _mm_blendv_epi8(ref_row, c, rdmask);
      __m128i e_masked     = _mm_blendv_epi8(ref_row, e, rdmask);
      __m128i g_masked     = _mm_blendv_epi8(ref_row, g, rdmask);

      __m128i curr_sads_ab = _mm_sad_epu8   (ref_row, a_masked);
      __m128i curr_sads_cd = _mm_sad_epu8   (ref_row, c_masked);
      __m128i curr_sads_ef = _mm_sad_epu8   (ref_row, e_masked);
      __m128i curr_sads_gh = _mm_sad_epu8   (ref_row, g_masked);

      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ef);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_gh);
    }
    if (height_residual_lines) {
      for (; y < height; y++) {
        __m128i a = _mm_loadu_si128((const __m128i *)(pic_data + y * stride + x));

        __m128i a_masked  = _mm_blendv_epi8(ref_row, a, rdmask);
        __m128i curr_sads = _mm_sad_epu8   (ref_row, a_masked);

        sse_inc = _mm_add_epi64(sse_inc, curr_sads);
      }
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);

  return _mm_cvtsi128_si32(sad);
}

static uint32_t hor_sad_sse41_w4(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                 int32_t height, uint32_t pic_stride, uint32_t ref_stride,
                                 uint32_t left, uint32_t right)
{
  int32_t leftoff = left;
  int8_t border_idx;
  if (left)
    border_idx = left;
  else
    border_idx = 3 - right;

  // Dualword (ie. line) base indexes, ie. the edges the lines read will be
  // clamped towards
  const __m128i dwbaseids   = _mm_setr_epi8(0, 0, 0, 0, 4, 4, 4, 4,
                                            8, 8, 8, 8, 12, 12, 12, 12);

  const __m128i border_idxs = _mm_set1_epi8(border_idx);
  const __m128i ns          = _mm_setr_epi8(0,  1,  2,  3,  4,  5,  6,  7,
                                            8,  9,  10, 11, 12, 13, 14, 15);
  __m128i epol_mask;
  if (left) {
    __m128i mask1 = _mm_sub_epi8(ns,    border_idxs);
    epol_mask     = _mm_max_epi8(mask1, dwbaseids);
  } else {
    if (right != 4) {
      __m128i border_idxs_linewise = _mm_add_epi8(border_idxs, dwbaseids);
      epol_mask = _mm_min_epi8(ns, border_idxs_linewise);
    } else {
      epol_mask = dwbaseids;
      leftoff = -1;
    }
  }
  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  __m128i sse_inc = _mm_setzero_si128();
  int32_t y;
  for (y = 0; y < height_fourline_groups; y += 4) {
    __m128i a = _mm_cvtsi32_si128(*(const uint32_t *)(pic_data + y * pic_stride));
    __m128i b = _mm_cvtsi32_si128(*(const uint32_t *)(ref_data + y * ref_stride + leftoff));

    a = _mm_insert_epi32(a, *(const uint32_t *)(pic_data + (y + 1) * pic_stride),           1);
    b = _mm_insert_epi32(b, *(const uint32_t *)(ref_data + (y + 1) * ref_stride + leftoff), 1);
    a = _mm_insert_epi32(a, *(const uint32_t *)(pic_data + (y + 2) * pic_stride),           2);
    b = _mm_insert_epi32(b, *(const uint32_t *)(ref_data + (y + 2) * ref_stride + leftoff), 2);
    a = _mm_insert_epi32(a, *(const uint32_t *)(pic_data + (y + 3) * pic_stride),           3);
    b = _mm_insert_epi32(b, *(const uint32_t *)(ref_data + (y + 3) * ref_stride + leftoff), 3);

    __m128i b_epol    = _mm_shuffle_epi8(b,       epol_mask);
    __m128i curr_sads = _mm_sad_epu8    (a,       b_epol);
            sse_inc   = _mm_add_epi64   (sse_inc, curr_sads);
  }
  if (height_residual_lines) {
    for (; y < height; y++) {
      __m128i a = _mm_cvtsi32_si128(*(const uint32_t *)(pic_data + y * pic_stride));
      __m128i b = _mm_cvtsi32_si128(*(const uint32_t *)(ref_data + y * ref_stride + leftoff));

      __m128i b_epol = _mm_shuffle_epi8(b, epol_mask);
      __m128i curr_sads = _mm_sad_epu8 (a, b_epol);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads);
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);

  return _mm_cvtsi128_si32(sad);
}

static uint32_t hor_sad_sse41_w8(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                 int32_t height, uint32_t pic_stride, uint32_t ref_stride,
                                 uint32_t left, uint32_t right)
{
  int32_t leftoff = left;
  int8_t border_idx;
  if (left)
    border_idx = left;
  else
    border_idx = 7 - right;

  // Quadword (ie. line) base indexes, ie. the edges the lines read will be
  // clamped towards; higher qword (lower line) bytes tend towards 8 and lower
  // qword (higher line) bytes towards 0
  const __m128i qwbaseids   = _mm_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0,
                                            8, 8, 8, 8, 8, 8, 8, 8);

  const __m128i border_idxs = _mm_set1_epi8(border_idx);
  const __m128i ns          = _mm_setr_epi8(0,  1,  2,  3,  4,  5,  6,  7,
                                            8,  9,  10, 11, 12, 13, 14, 15);
  __m128i epol_mask;
  if (left) {
    __m128i mask1     = _mm_sub_epi8(ns,    border_idxs);
    epol_mask         = _mm_max_epi8(mask1, qwbaseids);
  } else {
    if (right != 8) {
      __m128i border_idxs_linewise = _mm_add_epi8(border_idxs, qwbaseids);
      epol_mask = _mm_min_epi8(ns, border_idxs_linewise);
    } else {
      epol_mask = qwbaseids;
      leftoff = -1;
    }
  }
  const __m64 epol_mask_64 = (__m64)_mm_cvtsi128_si64(epol_mask);

  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  uint64_t result = 0;
  __m128i sse_inc = _mm_setzero_si128();
  int32_t y;
  for (y = 0; y < height_fourline_groups; y += 4) {
    __m128d a_d = _mm_setzero_pd();
    __m128d b_d = _mm_setzero_pd();
    __m128d c_d = _mm_setzero_pd();
    __m128d d_d = _mm_setzero_pd();

    a_d = _mm_loadl_pd(a_d, (const double *)(pic_data + (y + 0) * pic_stride));
    b_d = _mm_loadl_pd(b_d, (const double *)(ref_data + (y + 0) * ref_stride + leftoff));
    a_d = _mm_loadh_pd(a_d, (const double *)(pic_data + (y + 1) * pic_stride));
    b_d = _mm_loadh_pd(b_d, (const double *)(ref_data + (y + 1) * ref_stride + leftoff));

    c_d = _mm_loadl_pd(c_d, (const double *)(pic_data + (y + 2) * pic_stride));
    d_d = _mm_loadl_pd(d_d, (const double *)(ref_data + (y + 2) * ref_stride + leftoff));
    c_d = _mm_loadh_pd(c_d, (const double *)(pic_data + (y + 3) * pic_stride));
    d_d = _mm_loadh_pd(d_d, (const double *)(ref_data + (y + 3) * ref_stride + leftoff));

    __m128i a = _mm_castpd_si128(a_d);
    __m128i b = _mm_castpd_si128(b_d);
    __m128i c = _mm_castpd_si128(c_d);
    __m128i d = _mm_castpd_si128(d_d);

    __m128i b_epol = _mm_shuffle_epi8(b, epol_mask);
    __m128i d_epol = _mm_shuffle_epi8(d, epol_mask);

    __m128i curr_sads_ab = _mm_sad_epu8(a, b_epol);
    __m128i curr_sads_cd = _mm_sad_epu8(c, d_epol);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
  }
  if (height_residual_lines) {
    for (; y < height; y++) {
      __m64 a = *(__m64 *)(pic_data + y * pic_stride);
      __m64 b = *(__m64 *)(ref_data + y * ref_stride + leftoff);

      __m64 b_epol = _mm_shuffle_pi8(b, epol_mask_64);
      __m64 sads = _mm_sad_pu8(a, b_epol);
      result += (uint64_t)sads;
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);
  result += _mm_cvtsi128_si32(sad);
  return result;
}

/*
 * left and right measure how many pixels of one horizontal scanline will be
 * outside either the left or the right screen border. For blocks straddling
 * the left border, read the scanlines starting from the left border instead,
 * and use the extrapolation mask to essentially move the pixels right while
 * copying the left border pixel to the vector positions that logically point
 * outside of the buffer.
 *
 * For blocks straddling the right border, just read over the right border,
 * and extrapolate all pixels beyond the border idx to copy the value of the
 * border pixel. An exception is right == width (leftmost reference pixel is
 * one place right from the right border, it's ugly because the pixel to
 * extrapolate from is located at relative X offset -1), abuse the left border
 * aligning functionality instead to actually read starting from the valid
 * border pixel, and use a suitable mask to fill all the other pixels with
 * that value.
 */
static uint32_t hor_sad_sse41_w16(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                  int32_t height, uint32_t pic_stride, uint32_t ref_stride,
                                  const uint32_t left, const uint32_t right)
{
  // right is the number of overhanging pixels in the vector, so it has to be
  // handled this way to produce the index of last valid (border) pixel
  const int32_t right_border_idx = 15 - right;
  const int32_t border_idx       = left ? left : right_border_idx;

  const __m128i ns               = _mm_setr_epi8(0,  1,  2,  3,  4,  5,  6,  7,
                                                 8,  9,  10, 11, 12, 13, 14, 15);
  const __m128i zero             = _mm_setzero_si128();

  // Dirty hack alert! If right == block_width (ie. the entire vector is
  // outside the frame), move the block offset one pixel to the left (so
  // that the leftmost pixel in vector is actually the valid border pixel
  // from which we want to extrapolate), and use an epol mask that will
  // simply stretch the pixel all over the vector.
  //
  // To avoid a branch here:
  // The mask will be -1 (0xffffffff) for border_idx -1 and 0 for >= 0
  const int32_t border_idx_negative = border_idx >> 31;
  const int32_t leftoff             = border_idx_negative | left;

  __m128i right_border_idxs = _mm_set1_epi8((int8_t)right_border_idx);
  __m128i left_128          = _mm_set1_epi8((int8_t)left);

  // If we're straddling the left border, right_border_idx is 15 and the first
  // operation does nothing. If right border, left is 0 and the second
  // operation does nothing.
  __m128i mask_right        = _mm_min_epi8 (ns,         right_border_idxs);
  __m128i mask1             = _mm_sub_epi8 (mask_right, left_128);

  // If right == 16 (we're completely outside the frame), right_border_idx is
  // -1 and so is mask1. Clamp negative values to zero and as discussed
  // earlier, adjust the load offset instead to load the "-1'st" pixel and
  // using an all-zero shuffle mask, broadcast it all over the vector.
  const __m128i epol_mask = _mm_max_epi8(mask1, zero);

  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  __m128i sse_inc = _mm_setzero_si128();
  int32_t y;
  for (y = 0; y < height_fourline_groups; y += 4) {
    __m128i a = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * pic_stride));
    __m128i b = _mm_loadu_si128((__m128i *)(ref_data + (y + 0) * ref_stride + leftoff));
    __m128i c = _mm_loadu_si128((__m128i *)(pic_data + (y + 1) * pic_stride));
    __m128i d = _mm_loadu_si128((__m128i *)(ref_data + (y + 1) * ref_stride + leftoff));
    __m128i e = _mm_loadu_si128((__m128i *)(pic_data + (y + 2) * pic_stride));
    __m128i f = _mm_loadu_si128((__m128i *)(ref_data + (y + 2) * ref_stride + leftoff));
    __m128i g = _mm_loadu_si128((__m128i *)(pic_data + (y + 3) * pic_stride));
    __m128i h = _mm_loadu_si128((__m128i *)(ref_data + (y + 3) * ref_stride + leftoff));

    __m128i b_epol = _mm_shuffle_epi8(b, epol_mask);
    __m128i d_epol = _mm_shuffle_epi8(d, epol_mask);
    __m128i f_epol = _mm_shuffle_epi8(f, epol_mask);
    __m128i h_epol = _mm_shuffle_epi8(h, epol_mask);

    __m128i curr_sads_ab = _mm_sad_epu8(a, b_epol);
    __m128i curr_sads_cd = _mm_sad_epu8(c, d_epol);
    __m128i curr_sads_ef = _mm_sad_epu8(e, f_epol);
    __m128i curr_sads_gh = _mm_sad_epu8(g, h_epol);

    sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_ef);
    sse_inc = _mm_add_epi64(sse_inc, curr_sads_gh);
  }
  if (height_residual_lines) {
    for (; y < height; y++) {
      __m128i a = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * pic_stride));
      __m128i b = _mm_loadu_si128((__m128i *)(ref_data + (y + 0) * ref_stride + leftoff));
      __m128i b_epol = _mm_shuffle_epi8(b, epol_mask);
      __m128i curr_sads = _mm_sad_epu8(a, b_epol);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads);
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);
  return _mm_cvtsi128_si32(sad);
}

static uint32_t hor_sad_sse41_arbitrary(const kvz_pixel *pic_data, const kvz_pixel *ref_data,
                                        int32_t width, int32_t height, uint32_t pic_stride,
                                        uint32_t ref_stride, uint32_t left, uint32_t right)
{
  const size_t xmm_width = 16;
  const __m128i xmm_widths = _mm_set1_epi8(xmm_width);

  // Bytes in block in 128-bit blocks per each scanline, and remainder
  const int32_t width_xmms             = width  & ~(xmm_width - 1);
  const int32_t width_residual_pixels  = width  &  (xmm_width - 1);

  const int32_t height_fourline_groups = height & ~3;
  const int32_t height_residual_lines  = height &  3;

  __m128i ns = _mm_setr_epi8(0,  1,  2,  3,  4,  5,  6,  7,
                             8,  9,  10, 11, 12, 13, 14, 15);

  const __m128i rds    = _mm_set1_epi8 (width_residual_pixels);
  const __m128i rdmask = _mm_cmpgt_epi8(rds, ns);

  int32_t border_idx;
  __m128i is_right_border = _mm_setzero_si128();
  if (left) {
    border_idx = left;
  } else {
    border_idx = width - (right + 1);
    is_right_border = _mm_cmpeq_epi8(is_right_border, is_right_border);
  }
  const __m128i epol_src_idx = _mm_set1_epi8(border_idx);

  int32_t x, y;
  __m128i sse_inc = _mm_setzero_si128();
  __m128i epol_mask;
  for (x = 0; x < width_xmms; x += xmm_width) {

    // This is a dirty hack, but it saves us an easily predicted branch! It
    // also marks the first or last valid pixel (the border one) for
    // extrapolating, but that makes no difference since the pixels marked
    // for extrapolation will always be written over with that exact pixel's
    // value.
    epol_mask = _mm_cmpgt_epi8(epol_src_idx, ns);
    epol_mask = _mm_xor_si128 (epol_mask,    is_right_border);

    for (y = 0; y < height_fourline_groups; y += 4) {
      __m128i a = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * pic_stride + x));
      __m128i b = _mm_loadu_si128((__m128i *)(ref_data + (y + 0) * ref_stride + x));
      __m128i c = _mm_loadu_si128((__m128i *)(pic_data + (y + 1) * pic_stride + x));
      __m128i d = _mm_loadu_si128((__m128i *)(ref_data + (y + 1) * ref_stride + x));
      __m128i e = _mm_loadu_si128((__m128i *)(pic_data + (y + 2) * pic_stride + x));
      __m128i f = _mm_loadu_si128((__m128i *)(ref_data + (y + 2) * ref_stride + x));
      __m128i g = _mm_loadu_si128((__m128i *)(pic_data + (y + 3) * pic_stride + x));
      __m128i h = _mm_loadu_si128((__m128i *)(ref_data + (y + 3) * ref_stride + x));

      __m128i border_px_b  = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 0) * ref_stride + border_idx));
      __m128i border_px_d  = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 1) * ref_stride + border_idx));
      __m128i border_px_f  = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 2) * ref_stride + border_idx));
      __m128i border_px_h  = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 3) * ref_stride + border_idx));
      __m128i b_epol       = _mm_blendv_epi8(b, border_px_b, epol_mask);
      __m128i d_epol       = _mm_blendv_epi8(d, border_px_d, epol_mask);
      __m128i f_epol       = _mm_blendv_epi8(f, border_px_f, epol_mask);
      __m128i h_epol       = _mm_blendv_epi8(h, border_px_h, epol_mask);

      __m128i curr_sads_ab = _mm_sad_epu8(a, b_epol);
      __m128i curr_sads_cd = _mm_sad_epu8(c, d_epol);
      __m128i curr_sads_ef = _mm_sad_epu8(e, f_epol);
      __m128i curr_sads_gh = _mm_sad_epu8(g, h_epol);

      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ef);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_gh);
    }
    if (height_residual_lines) {
      for (; y < height; y++) {
        __m128i a = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * pic_stride + x));
        __m128i b = _mm_loadu_si128((__m128i *)(ref_data + (y + 0) * ref_stride + x));

        __m128i border_px_b  = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 0) * ref_stride + border_idx));
        __m128i b_epol       = _mm_blendv_epi8(b, border_px_b, epol_mask);

        __m128i curr_sads_ab = _mm_sad_epu8(a, b_epol);

        sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
      }
    }
    ns = _mm_add_epi8(ns, xmm_widths);
  }
  if (width_residual_pixels) {
    epol_mask = _mm_cmpgt_epi8(epol_src_idx, ns);
    epol_mask = _mm_xor_si128 (epol_mask,    is_right_border);

    for (y = 0; y < height_fourline_groups; y += 4) {
      __m128i a = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * pic_stride + x));
      __m128i b = _mm_loadu_si128((__m128i *)(ref_data + (y + 0) * ref_stride + x));
      __m128i c = _mm_loadu_si128((__m128i *)(pic_data + (y + 1) * pic_stride + x));
      __m128i d = _mm_loadu_si128((__m128i *)(ref_data + (y + 1) * ref_stride + x));
      __m128i e = _mm_loadu_si128((__m128i *)(pic_data + (y + 2) * pic_stride + x));
      __m128i f = _mm_loadu_si128((__m128i *)(ref_data + (y + 2) * ref_stride + x));
      __m128i g = _mm_loadu_si128((__m128i *)(pic_data + (y + 3) * pic_stride + x));
      __m128i h = _mm_loadu_si128((__m128i *)(ref_data + (y + 3) * ref_stride + x));

      __m128i border_px_b  = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 0) * ref_stride + border_idx));
      __m128i border_px_d  = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 1) * ref_stride + border_idx));
      __m128i border_px_f  = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 2) * ref_stride + border_idx));
      __m128i border_px_h  = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 3) * ref_stride + border_idx));

      __m128i b_epol_1     = _mm_blendv_epi8(b, border_px_b, epol_mask);
      __m128i d_epol_1     = _mm_blendv_epi8(d, border_px_d, epol_mask);
      __m128i f_epol_1     = _mm_blendv_epi8(f, border_px_f, epol_mask);
      __m128i h_epol_1     = _mm_blendv_epi8(h, border_px_h, epol_mask);

      __m128i b_epol_2     = _mm_blendv_epi8(a, b_epol_1,    rdmask);
      __m128i d_epol_2     = _mm_blendv_epi8(c, d_epol_1,    rdmask);
      __m128i f_epol_2     = _mm_blendv_epi8(e, f_epol_1,    rdmask);
      __m128i h_epol_2     = _mm_blendv_epi8(g, h_epol_1,    rdmask);

      __m128i curr_sads_ab = _mm_sad_epu8(a, b_epol_2);
      __m128i curr_sads_cd = _mm_sad_epu8(c, d_epol_2);
      __m128i curr_sads_ef = _mm_sad_epu8(e, f_epol_2);
      __m128i curr_sads_gh = _mm_sad_epu8(g, h_epol_2);

      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_cd);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_ef);
      sse_inc = _mm_add_epi64(sse_inc, curr_sads_gh);
    }
    if (height_residual_lines) {
      for (; y < height; y++) {
        __m128i a = _mm_loadu_si128((__m128i *)(pic_data + (y + 0) * pic_stride + x));
        __m128i b = _mm_loadu_si128((__m128i *)(ref_data + (y + 0) * ref_stride + x));

        __m128i border_px_b  = _mm_set1_epi8  (*(uint8_t *)(ref_data + (y + 0) * ref_stride + border_idx));
        __m128i b_epol_1     = _mm_blendv_epi8(b, border_px_b, epol_mask);
        __m128i b_epol_2     = _mm_blendv_epi8(a, b_epol_1,    rdmask);

        __m128i curr_sads_ab = _mm_sad_epu8(a, b_epol_2);

        sse_inc = _mm_add_epi64(sse_inc, curr_sads_ab);
      }
    }
  }
  __m128i sse_inc_2 = _mm_shuffle_epi32(sse_inc, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad       = _mm_add_epi64    (sse_inc, sse_inc_2);
  return _mm_cvtsi128_si32(sad);
}

#endif
