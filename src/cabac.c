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

#include "cabac.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


const uint8_t g_auc_next_state_mps[128] =
{
    2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17,
   18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,
   34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,
   50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,
   66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,
   82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,
   98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113,
  114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 124, 125, 126, 127
};

const uint8_t g_auc_next_state_lps[128] =
{
   1,  0,  0,  1,  2,  3,  4,  5,  4,  5,  8,  9,  8,  9,  10,  11,
  12, 13, 14, 15, 16, 17, 18, 19, 18, 19, 22, 23, 22, 23,  24,  25,
  26, 27, 26, 27, 30, 31, 30, 31, 32, 33, 32, 33, 36, 37,  36,  37,
  38, 39, 38, 39, 42, 43, 42, 43, 44, 45, 44, 45, 46, 47,  48,  49,
  48, 49, 50, 51, 52, 53, 52, 53, 54, 55, 54, 55, 56, 57,  58,  59,
  58, 59, 60, 61, 60, 61, 60, 61, 62, 63, 64, 65, 64, 65,  66,  67,
  66, 67, 66, 67, 68, 69, 68, 69, 70, 71, 70, 71, 70, 71,  72,  73,
  72, 73, 72, 73, 74, 75, 74, 75, 74, 75, 76, 77, 76, 77, 126, 127
};

const uint8_t g_auc_lpst_table[64][4] =
{
  {128, 176, 208, 240},  {128, 167, 197, 227},  {128, 158, 187, 216},  {123, 150, 178, 205},  {116, 142, 169, 195},
  {111, 135, 160, 185},  {105, 128, 152, 175},  {100, 122, 144, 166},  { 95, 116, 137, 158},  { 90, 110, 130, 150},
  { 85, 104, 123, 142},  { 81,  99, 117, 135},  { 77,  94, 111, 128},  { 73,  89, 105, 122},  { 69,  85, 100, 116},
  { 66,  80,  95, 110},  { 62,  76,  90, 104},  { 59,  72,  86,  99},  { 56,  69,  81,  94},  { 53,  65,  77,  89},
  { 51,  62,  73,  85},  { 48,  59,  69,  80},  { 46,  56,  66,  76},  { 43,  53,  63,  72},  { 41,  50,  59,  69},
  { 39,  48,  56,  65},  { 37,  45,  54,  62},  { 35,  43,  51,  59},  { 33,  41,  48,  56},  { 32,  39,  46,  53},
  { 30,  37,  43,  50},  { 29,  35,  41,  48},  { 27,  33,  39,  45},  { 26,  31,  37,  43},  { 24,  30,  35,  41},
  { 23,  28,  33,  39},  { 22,  27,  32,  37},  { 21,  26,  30,  35},  { 20,  24,  29,  33},  { 19,  23,  27,  31},
  { 18,  22,  26,  30},  { 17,  21,  25,  28},  { 16,  20,  23,  27},  { 15,  19,  22,  25},  { 14,  18,  21,  24},
  { 14,  17,  20,  23},  { 13,  16,  19,  22},  { 12,  15,  18,  21},  { 12,  14,  17,  20},  { 11,  14,  16,  19},
  { 11,  13,  15,  18},  { 10,  12,  15,  17},  { 10,  12,  14,  16},  {  9,  11,  13,  15},  {  9,  11,  12,  14},
  {  8,  10,  12,  14},  {  8,   9,  11,  13},  {  7,   9,  11,  12},  {  7,   9,  10,  12},  {  7,   8,  10,  11},
  {  6,   8,   9,  11},  {  6,   7,   9,  10},  {  6,   7,   8,   9},  {  2,   2,   2,   2}
};

const uint8_t g_auc_renorm_table[32] =
{
  6, 5, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

cabac_data cabac;

/**
 * \brief Initialize struct cabac_ctx.
 */
void ctx_init(cabac_ctx *ctx, uint32_t qp, uint32_t init_value)
{
  int slope = (init_value >> 4) * 5 - 45;
  int offset = ((init_value & 15) << 3) - 16;
  int init_state = MIN(MAX(1, ((slope * (int)qp) >> 4) + offset), 126);
  uint8_t mp_state = (init_state >= 64) ? 1 : 0;

  if (mp_state) {
    ctx->uc_state = (init_state - 64) << 1 + mp_state;
  } else {
    ctx->uc_state = (63 - init_state) << 1;
  }
  ctx->bins_coded = 0;
}

/**
 * \brief Initialize struct cabac_data.
 */
void cabac_init(cabac_data* data)
{
  data->frac_bits = 0;
  data->bin_count_increment = 0;
  data->bins_coded = 0;
}

/**
 * \brief Initialize struct cabac_data.
 */
void cabac_start(cabac_data *data)
{
  data->low = 0;
  data->range = 510;
  data->bits_left = 23;
  data->num_buffered_bytes = 0;
  data->buffered_byte = 0xff;
}

/**
 * \brief
 */
void cabac_encode_bin(cabac_data *data, uint32_t bin_value)
{
  uint32_t lps;
  
  data->bins_coded += data->bin_count_increment;
  data->ctx->bins_coded = 1;
  
  lps = g_auc_lpst_table[CTX_STATE(data->ctx)][(data->range >> 6) & 3];
  data->range -= lps;
  
  // Not the Most Probable Symbol?
  if (bin_value != CTX_MPS(data->ctx)) {
    int num_bits = g_auc_renorm_table[lps >> 3];
    data->low = (data->low + data->range) << num_bits;
    data->range = lps << num_bits;
    
    CTX_UPDATE_LPS(data->ctx);
    
    data->bits_left -= num_bits;
  } else {
    CTX_UPDATE_MPS(data->ctx);
    if (data->range >= 256) return;
    
    data->low <<= 1;
    data->range <<= 1;
    data->bits_left--;
  }
  
  if (data->bits_left < 12) {
    cabac_write(data);
  }
}

/**
 * \brief
 */
void cabac_write(cabac_data *data)
{
  uint32_t lead_byte = data->low >> (24 - data->bits_left);
  data->bits_left += 8;
  data->low &= 0xffffffffu >> data->bits_left;
  
  if (lead_byte == 0xff) {
    data->num_buffered_bytes++;
  } else {
    if (data->num_buffered_bytes > 0) {
      uint32_t carry = lead_byte >> 8;
      uint32_t byte = data->buffered_byte + carry;
      data->buffered_byte = lead_byte & 0xff;
      bitstream_put(data->stream, byte, 8);

      byte = (0xff + carry) & 0xff;
      while (data->num_buffered_bytes > 1) {
        bitstream_put(data->stream, byte, 8);
        data->num_buffered_bytes--;
      }
    } else {
      data->num_buffered_bytes = 1;
      data->buffered_byte = lead_byte;
    }
  }
}

/**
 * \brief
 */
void cabac_finish(cabac_data *data)
{
  if (data->low >> (32 - data->bits_left)) {
    bitstream_put(data->stream,data->buffered_byte + 1, 8);
    while (data->num_buffered_bytes > 1) {
      bitstream_put(data->stream, 0, 8);
      data->num_buffered_bytes--;
    }
    data->low -= 1 << (32 - data->bits_left);
  } else {
    if (data->num_buffered_bytes > 0) {
      bitstream_put(data->stream,data->buffered_byte, 8);
    }
    while (data->num_buffered_bytes > 1) {
      bitstream_put(data->stream, 0xff, 8);
      data->num_buffered_bytes--;
    }
  }
  bitstream_put(data->stream, data->low >> 8, 24 - data->bits_left);
}

/*!
  \brief Encode terminating bin
  \param binValue bin value
*/
void cabac_encode_bin_trm(cabac_data *data, uint8_t bin_value)
{
  data->bins_coded += data->bin_count_increment;
  data->range -= 2;
  if(bin_value) {
    data->low += data->range;
    data->low <<= 7;
    data->range = 2 << 7;
    data->bits_left -= 7;
  } else if (data->range >= 256) {
    return;
  } else {
    data->low <<= 1;
    data->range <<= 1;
    data->bits_left--;
  }
  
  if (data->bits_left < 12) {
    cabac_write(data);
  }
}

/**
 * \brief
 */
void cabac_flush(cabac_data *data)
{
  cabac_encode_bin_trm(data, 1);
  cabac_finish(data);
  bitstream_put(data->stream, 1, 1);
  bitstream_align_zero(data->stream);
  cabac_start(data);
}

/**
 * \brief
 */
void cabac_encode_bin_ep(cabac_data *data, uint32_t bin_value)
{
  data->bins_coded += data->bin_count_increment;
  data->low <<= 1;
  if (bin_value) {
    data->low += data->range;
  }
  data->bits_left--;

  if (data->bits_left < 12) {
    cabac_write(data);
  }
}

/**
 * \brief
 */
void cabac_encode_bins_ep(cabac_data *data, uint32_t bin_values, int num_bins)
{
  uint32_t pattern;
  data->bins_coded += num_bins & -data->bin_count_increment;

  while (num_bins > 8) {
    num_bins -= 8;
    pattern = bin_values >> num_bins;
    data->low <<= 8;
    data->low += data->range * pattern;
    bin_values -= pattern << num_bins;
    data->bits_left -= 8;
    
    if(data->bits_left < 12) {
      cabac_write(data);
    }
  }
  
  data->low <<= num_bins;
  data->low += data->range * bin_values;
  data->bits_left -= num_bins;
  
  if (data->bits_left < 12) {
    cabac_write(data);
  }
}

/**
 * \brief Coding of coeff_abs_level_minus3.
 * \param symbol Value of coeff_abs_level_minus3.
 * \param r_param Reference to Rice parameter.
 */
void cabac_write_coeff_remain(cabac_data *cabac, uint32_t symbol, uint32_t r_param)
{
  int32_t code_number = symbol;
  uint32_t length;
  if (code_number < (3 << r_param)) {
    length = code_number >> r_param;
    cabac_encode_bins_ep(cabac, (1 << (length + 1)) - 2 , length + 1);
    cabac_encode_bins_ep(cabac, (code_number % (1 << r_param)), r_param);
  } else {
    length = r_param;
    code_number = code_number - (3 << r_param);
    while (code_number >= (1 << length)) {
      code_number -= 1 << length;
      ++length;
    }
    cabac_encode_bins_ep(cabac, (1 << (3 + length + 1 - r_param)) - 2, 3 + length + 1 - r_param);
    cabac_encode_bins_ep(cabac, code_number, length);
  }
}

/**
 * \brief
 */
void cabac_write_unary_max_symbol(cabac_data *data, cabac_ctx *ctx, uint32_t symbol, int32_t offset, uint32_t max_symbol)
{
  int8_t code_last = max_symbol > symbol;

  if (!max_symbol) return;
  
  data->ctx = &ctx[0];
  cabac_encode_bin(data, symbol ? 1 : 0);
  
  if (!symbol) return;
  
  while (--symbol) {
    data->ctx = &ctx[offset];
    cabac_encode_bin(data, 1);
  }
  if (code_last) {
    data->ctx = &ctx[offset];
    cabac_encode_bin(data, 0);
  }
}

/**
 * \brief
 */
void cabac_write_ep_ex_golomb(cabac_data *data, uint32_t symbol, uint32_t count)
{
  uint32_t bins = 0;
  int32_t num_bins = 0;
  
  while (symbol >= (uint32_t)(1 << count)) {
    bins = 2 * bins + 1;
    ++num_bins;
    symbol -= 1 << count;
    ++count;
  }
  bins = 2 * bins;
  ++num_bins;
  
  bins = (bins << count) | symbol;
  num_bins += count;
  
  cabac_encode_bins_ep(data, bins, num_bins);
}
