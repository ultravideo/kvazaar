#ifndef CABAC_H_
#define CABAC_H_
/**
 * \file
 * \brief The Content Adaptive Binary Arithmetic Coder (CABAC).
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "global.h"

#include "bitstream.h"
#include "context.h"


// Types
typedef struct
{
  cabac_ctx *ctx;
  uint32_t   low;
  uint32_t   range;
  uint32_t   buffered_byte;
  int32_t    num_buffered_bytes;
  int32_t    bits_left;
  bitstream *stream;
} cabac_data;


// Globals
extern const uint8_t g_auc_next_state_mps[128];
extern const uint8_t g_auc_next_state_lps[128];
extern const uint8_t g_auc_lpst_table[64][4];
extern const uint8_t g_auc_renorm_table[32];
extern cabac_data cabac;


// Functions
void cabac_start(cabac_data *data);
void cabac_encode_bin(cabac_data *data, uint32_t bin_value);
void cabac_encode_bin_ep(cabac_data *data, uint32_t bin_value);
void cabac_encode_bins_ep(cabac_data *data, uint32_t bin_values, int num_bins);
void cabac_encode_bin_trm(cabac_data *data, uint8_t bin_value);
void cabac_write(cabac_data *data);
void cabac_finish(cabac_data *data);
void cabac_flush(cabac_data *data);
void cabac_write_coeff_remain(cabac_data *cabac, uint32_t symbol, 
                              uint32_t r_param);
void cabac_write_ep_ex_golomb(cabac_data *data, uint32_t symbol, 
                              uint32_t count);
void cabac_write_unary_max_symbol(cabac_data *data, cabac_ctx *ctx, 
                                  uint32_t symbol, int32_t offset, 
                                  uint32_t max_symbol);


// Macros
#define CTX_STATE(ctx) (ctx->uc_state >> 1)
#define CTX_MPS(ctx) (ctx->uc_state & 1)
#define CTX_UPDATE_LPS(ctx) { (ctx)->uc_state = g_auc_next_state_lps[ (ctx)->uc_state ]; }
#define CTX_UPDATE_MPS(ctx) { (ctx)->uc_state = g_auc_next_state_mps[ (ctx)->uc_state ]; }
#ifdef VERBOSE
  #define CABAC_BIN(data, value, name) { \
    uint32_t prev_state = (data)->ctx->uc_state; \
    cabac_encode_bin(data, value); \
    printf("%s = %d prev_state=%d state=%d\n", \
           name, value, prev_state, (data)->ctx->uc_state); }

  #define CABAC_BINS_EP(data, value, bins, name) { \
    uint32_t prev_state = (data)->ctx->uc_state; \
    cabac_encode_bins_ep(data, value, bins); \
    printf("%s = %d prev_state=%d state=%d\n", \
           name, value, prev_state, (data)->ctx->uc_state); }

  #define CABAC_BIN_EP(data, value, name) { \
    uint32_t prev_state = (data)->ctx->uc_state; \
    cabac_encode_bin_ep(data, value); \
    printf("%s = %d prev_state=%d state=%d\n", \
           name, value, prev_state, (data)->ctx->uc_state); }
#else
  #define CABAC_BIN(data, value, name) \
    cabac_encode_bin(data, value);
  #define CABAC_BINS_EP(data, value, bins, name) \
    cabac_encode_bins_ep(data, value, bins);
  #define CABAC_BIN_EP(data, value, name) \
    cabac_encode_bin_ep(data, value);
#endif

#endif
