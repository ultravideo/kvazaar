#ifndef REG_SAD_POW2_WIDTHS_AVX2_H_
#define REG_SAD_POW2_WIDTHS_AVX2_H_

#include "strategies/sse41/reg_sad_pow2_widths-sse41.h"
#include "kvazaar.h"

static uint32_t reg_sad_w16_avx2(const kvz_pixel * const data1, const kvz_pixel * const data2,
                                 const int32_t height, const uint32_t stride1,
                                 const uint32_t stride2)
{
  __m256i avx_inc = _mm256_setzero_si256();
  int32_t y;

  const int32_t height_ymm_bytes = height & ~1;
  const int32_t height_parity    = height &  1;

  for (y = 0; y < height_ymm_bytes; y += 2) {
    __m128i a_up = _mm_loadu_si128((const __m128i *)(data1 + (y + 0) * stride1));
    __m128i b_up = _mm_loadu_si128((const __m128i *)(data2 + (y + 0) * stride2));
    __m128i a_dn = _mm_loadu_si128((const __m128i *)(data1 + (y + 1) * stride1));
    __m128i b_dn = _mm_loadu_si128((const __m128i *)(data2 + (y + 1) * stride2));

    __m256i a = _mm256_inserti128_si256(_mm256_castsi128_si256(a_up), a_dn, 1);
    __m256i b = _mm256_inserti128_si256(_mm256_castsi128_si256(b_up), b_dn, 1);

    __m256i curr_sads = _mm256_sad_epu8(a, b);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads);
  }
  __m128i inchi = _mm256_extracti128_si256(avx_inc, 1);
  __m128i inclo = _mm256_castsi256_si128  (avx_inc);

  if (height_parity) {
    __m128i a = _mm_loadu_si128 ((__m128i *)(data1 + y * stride1));
    __m128i b = _mm_loadu_si128 ((__m128i *)(data2 + y * stride2));
    __m128i sads = _mm_sad_epu8 (a, b);
    inclo        = _mm_add_epi64(inclo, sads);
  }
  __m128i sum_1 = _mm_add_epi64    (inclo, inchi);
  __m128i sum_2 = _mm_shuffle_epi32(sum_1, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad   = _mm_add_epi64    (sum_1, sum_2);

  return _mm_cvtsi128_si32(sad);
}

static uint32_t reg_sad_w32(const kvz_pixel * const data1, const kvz_pixel * const data2,
                            const int32_t height, const uint32_t stride1,
                            const uint32_t stride2)
{
  __m256i avx_inc = _mm256_setzero_si256();
  int32_t y;

  for (y = 0; y < height; y++) {
    __m256i a = _mm256_loadu_si256((const __m256i *)(data1 + y * stride1));
    __m256i b = _mm256_loadu_si256((const __m256i *)(data2 + y * stride2));

    __m256i curr_sads = _mm256_sad_epu8(a, b);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads);
  }
  __m128i inchi = _mm256_extracti128_si256(avx_inc, 1);
  __m128i inclo = _mm256_castsi256_si128  (avx_inc);

  __m128i sum_1 = _mm_add_epi64    (inclo, inchi);
  __m128i sum_2 = _mm_shuffle_epi32(sum_1, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad   = _mm_add_epi64    (sum_1, sum_2);

  return _mm_cvtsi128_si32(sad);
}

static uint32_t reg_sad_w64(const kvz_pixel * const data1, const kvz_pixel * const data2,
                            const int32_t height, const uint32_t stride1,
                            const uint32_t stride2)
{
  __m256i avx_inc = _mm256_setzero_si256();
  int32_t y;

  for (y = 0; y < height; y++) {
    __m256i a = _mm256_loadu_si256((const __m256i *)(data1 + y * stride1));
    __m256i b = _mm256_loadu_si256((const __m256i *)(data2 + y * stride2));
    __m256i c = _mm256_loadu_si256((const __m256i *)(data1 + y * stride1 + 32));
    __m256i d = _mm256_loadu_si256((const __m256i *)(data2 + y * stride2 + 32));

    __m256i curr_sads_1 = _mm256_sad_epu8(a, b);
    __m256i curr_sads_2 = _mm256_sad_epu8(c, d);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_1);
    avx_inc = _mm256_add_epi64(avx_inc, curr_sads_2);
  }
  __m128i inchi = _mm256_extracti128_si256(avx_inc, 1);
  __m128i inclo = _mm256_castsi256_si128  (avx_inc);

  __m128i sum_1 = _mm_add_epi64    (inclo, inchi);
  __m128i sum_2 = _mm_shuffle_epi32(sum_1, _MM_SHUFFLE(1, 0, 3, 2));
  __m128i sad   = _mm_add_epi64    (sum_1, sum_2);

  return _mm_cvtsi128_si32(sad);
}

#endif
