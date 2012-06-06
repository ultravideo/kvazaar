/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file cabac.c
    \brief CABAC
    \author Marko Viitanen
    \date 2012-06
    
    Content-adaptive binary arithmetic coder
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "bitstream.h"
#include "cabac.h"


const uint8_t g_aucNextStateMPS[ 128 ] =
{
  2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
  18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
  34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
  50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65,
  66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81,
  82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97,
  98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113,
  114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 124, 125, 126, 127
};

const uint8_t g_aucNextStateLPS[ 128 ] =
{
  1, 0, 0, 1, 2, 3, 4, 5, 4, 5, 8, 9, 8, 9, 10, 11,
  12, 13, 14, 15, 16, 17, 18, 19, 18, 19, 22, 23, 22, 23, 24, 25,
  26, 27, 26, 27, 30, 31, 30, 31, 32, 33, 32, 33, 36, 37, 36, 37,
  38, 39, 38, 39, 42, 43, 42, 43, 44, 45, 44, 45, 46, 47, 48, 49,
  48, 49, 50, 51, 52, 53, 52, 53, 54, 55, 54, 55, 56, 57, 58, 59,
  58, 59, 60, 61, 60, 61, 60, 61, 62, 63, 64, 65, 64, 65, 66, 67,
  66, 67, 66, 67, 68, 69, 68, 69, 70, 71, 70, 71, 70, 71, 72, 73,
  72, 73, 72, 73, 74, 75, 74, 75, 74, 75, 76, 77, 76, 77, 126, 127
};

const uint32_t g_entropyBits[128] =
{
  0x07b23, 0x085f9, 0x074a0, 0x08cbc, 0x06ee4, 0x09354, 0x067f4, 0x09c1b, 0x060b0, 0x0a62a, 0x05a9c, 0x0af5b, 0x0548d, 0x0b955, 0x04f56, 0x0c2a9,
  0x04a87, 0x0cbf7, 0x045d6, 0x0d5c3, 0x04144, 0x0e01b, 0x03d88, 0x0e937, 0x039e0, 0x0f2cd, 0x03663, 0x0fc9e, 0x03347, 0x10600, 0x03050, 0x10f95,
  0x02d4d, 0x11a02, 0x02ad3, 0x12333, 0x0286e, 0x12cad, 0x02604, 0x136df, 0x02425, 0x13f48, 0x021f4, 0x149c4, 0x0203e, 0x1527b, 0x01e4d, 0x15d00,
  0x01c99, 0x166de, 0x01b18, 0x17017, 0x019a5, 0x17988, 0x01841, 0x18327, 0x016df, 0x18d50, 0x015d9, 0x19547, 0x0147c, 0x1a083, 0x0138e, 0x1a8a3,
  0x01251, 0x1b418, 0x01166, 0x1bd27, 0x01068, 0x1c77b, 0x00f7f, 0x1d18e, 0x00eda, 0x1d91a, 0x00e19, 0x1e254, 0x00d4f, 0x1ec9a, 0x00c90, 0x1f6e0,
  0x00c01, 0x1fef8, 0x00b5f, 0x208b1, 0x00ab6, 0x21362, 0x00a15, 0x21e46, 0x00988, 0x2285d, 0x00934, 0x22ea8, 0x008a8, 0x239b2, 0x0081d, 0x24577,
  0x007c9, 0x24ce6, 0x00763, 0x25663, 0x00710, 0x25e8f, 0x006a0, 0x26a26, 0x00672, 0x26f23, 0x005e8, 0x27ef8, 0x005ba, 0x284b5, 0x0055e, 0x29057,
  0x0050c, 0x29bab, 0x004c1, 0x2a674, 0x004a7, 0x2aa5e, 0x0046f, 0x2b32f, 0x0041f, 0x2c0ad, 0x003e7, 0x2ca8d, 0x003ba, 0x2d323, 0x0010c, 0x3bfbb
};

const uint8_t g_aucLPSTable[64][4] =
{
  { 128, 176, 208, 240},  { 128, 167, 197, 227},  { 128, 158, 187, 216},  { 123, 150, 178, 205},  { 116, 142, 169, 195},
  { 111, 135, 160, 185},  { 105, 128, 152, 175},  { 100, 122, 144, 166},  {  95, 116, 137, 158},  {  90, 110, 130, 150},
  {  85, 104, 123, 142},  {  81,  99, 117, 135},  {  77,  94, 111, 128},  {  73,  89, 105, 122},  {  69,  85, 100, 116},
  {  66,  80,  95, 110},  {  62,  76,  90, 104},  {  59,  72,  86,  99},  {  56,  69,  81,  94},  {  53,  65,  77,  89},
  {  51,  62,  73,  85},  {  48,  59,  69,  80},  {  46,  56,  66,  76},  {  43,  53,  63,  72},  {  41,  50,  59,  69},
  {  39,  48,  56,  65},  {  37,  45,  54,  62},  {  35,  43,  51,  59},  {  33,  41,  48,  56},  {  32,  39,  46,  53},
  {  30,  37,  43,  50},  {  29,  35,  41,  48},  {  27,  33,  39,  45},  {  26,  31,  37,  43},  {  24,  30,  35,  41},
  {  23,  28,  33,  39},  {  22,  27,  32,  37},  {  21,  26,  30,  35},  {  20,  24,  29,  33},  {  19,  23,  27,  31},
  {  18,  22,  26,  30},  {  17,  21,  25,  28},  {  16,  20,  23,  27},  {  15,  19,  22,  25},  {  14,  18,  21,  24},
  {  14,  17,  20,  23},  {  13,  16,  19,  22},  {  12,  15,  18,  21},  {  12,  14,  17,  20},  {  11,  14,  16,  19},
  {  11,  13,  15,  18},  {  10,  12,  15,  17},  {  10,  12,  14,  16},  {   9,  11,  13,  15},  {   9,  11,  12,  14},
  {   8,  10,  12,  14},  {   8,   9,  11,  13},  {   7,   9,  11,  12},  {   7,   9,  10,  12},  {   7,   8,  10,  11},
  {   6,   8,   9,  11},  {   6,   7,   9,  10},  {   6,   7,   8,   9},  {   2,   2,   2,   2}
};

const uint8_t g_aucRenormTable[32] = {  6,  5,  4,  4,  3,  3,  3,  3,  2,  2,  2,  2,  2,  2,  2,  2,
                                        1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1};


uint8_t g_nextState[128][2];

cabac_data cabac;


void cxt_init(cabac_ctx* ctx, uint32_t qp, uint32_t initValue )
{
  int  slope       = (initValue>>4)*5 - 45;
  int  offset      = ((initValue&15)<<3)-16;
  int  initState   =  MIN( MAX( 1, ( ( ( slope * qp ) >> 4 ) + offset ) ), 126 );
  uint8_t mpState  = (initState >= 64 )?1:0;
  ctx->ucState     = ( (mpState? (initState - 64):(63 - initState)) <<1) + mpState;
  ctx->binsCoded   = 0;
}


void cxt_buildNextStateTable()
{
  int i,j;
  for (i = 0; i < 128; i++)
  {
    for (j = 0; j < 2; j++)
    {
      g_nextState[i][j] = ((i&1) == j) ? g_aucNextStateMPS[i] : g_aucNextStateLPS[i];
    }
  }
}

void ctx_update(cabac_ctx* ctx, int val ) { ctx->ucState = g_nextState[ctx->ucState][val]; }
void ctx_update_LPS(cabac_ctx* ctx) { ctx->ucState = g_aucNextStateLPS[ ctx->ucState ]; }
void ctx_update_MPS(cabac_ctx* ctx) { ctx->ucState = g_aucNextStateMPS[ ctx->ucState ]; }

void cabac_init(cabac_data* data)
{
  data->fracBits = 0;

  cxt_buildNextStateTable();
}

void cabac_start(cabac_data* data)
{
  data->uiLow            = 0;
  data->uiRange          = 510;
  data->bitsLeft         = 23;
  data->numBufferedBytes = 0;
  data->bufferedByte     = 0xff;
}



void cabac_encodeBin(cabac_data* data, uint32_t binValue )
{
  uint32_t uiLPS;
  data->uiBinsCoded += data->binCountIncrement;
  data->ctx.binsCoded = 1;
  
  uiLPS   = g_aucLPSTable[ CTX_STATE(data->ctx) ][ ( data->uiRange >> 6 ) & 3 ];
  data->uiRange    -= uiLPS;
  
  if( binValue != CTX_MPS(data->ctx) )
  {
    int numBits = g_aucRenormTable[ uiLPS >> 3 ];
    data->uiLow     = ( data->uiLow + data->uiRange ) << numBits;
    data->uiRange   = uiLPS << numBits;
    
    ctx_update_LPS(&data->ctx);
    
    data->bitsLeft -= numBits;
  }
  else
  {
    ctx_update_MPS(&data->ctx);
    if (  data->uiRange >= 256 )
    {
      return;
    }
    
    data->uiLow <<= 1;
    data->uiRange <<= 1;
    data->bitsLeft--;
  }
  
  if(data->bitsLeft < 12)
  {
    cabac_write(data);
  }
}

void cabac_write(cabac_data* data)
{
  uint32_t leadByte = data->uiLow >> (24 - data->bitsLeft);
  data->bitsLeft += 8;
  data->uiLow &= 0xffffffffu >> data->bitsLeft;
  
  if ( leadByte == 0xff )
  {
    data->numBufferedBytes++;
  }
  else
  {
    if ( data->numBufferedBytes > 0 )
    {
      uint32_t carry = leadByte >> 8;
      uint32_t byte = data->bufferedByte + carry;
      data->bufferedByte = leadByte & 0xff;
      bitstream_put(data->stream,byte,8);

      byte = ( 0xff + carry ) & 0xff;
      while ( data->numBufferedBytes > 1 )
      {
        bitstream_put(data->stream,byte,8);
        data->numBufferedBytes--;
      }
    }
    else
    {
      data->numBufferedBytes = 1;
      data->bufferedByte = leadByte;
    }
  }
}

void cabac_encodeFlush(cabac_data* data, uint8_t end )
{
  data->uiRange = 2;

  data->uiLow  += 2;
  data->uiLow <<= 7;
  data->uiRange = 2 << 7;
  data->bitsLeft -= 7;
  if(data->bitsLeft < 12)
  {
    cabac_write(data);
  }
  cabac_finish(data);

  if(!end)
  {
    bitstream_put(data->stream,1,1);
  }
}

void cabac_finish(cabac_data* data)
{
  if ( data->uiLow >> ( 32 - data->bitsLeft ) )
  {
    bitstream_put(data->stream,data->bufferedByte + 1, 8 );
    while ( data->numBufferedBytes > 1 )
    {
      bitstream_put(data->stream,0, 8 );
      data->numBufferedBytes--;
    }
    data->uiLow -= 1 << ( 32 - data->bitsLeft );
  }
  else
  {
    if ( data->numBufferedBytes > 0 )
    {
      bitstream_put(data->stream,data->bufferedByte, 8 );
    }
    while ( data->numBufferedBytes > 1 )
    {
      bitstream_put(data->stream, 0xff, 8 );
      data->numBufferedBytes--;
    }
  }
  bitstream_put(data->stream, data->uiLow >> 8, 24 - data->bitsLeft );
}

/**
 * \brief Encode terminating bin
 *
 * \param binValue bin value
 */
void cabac_encodeBinTrm(cabac_data* data, uint32_t binValue )
{
  data->uiBinsCoded += data->binCountIncrement;
  data->uiRange -= 2;
  if( binValue )
  {
    data->uiLow  += data->uiRange;
    data->uiLow <<= 7;
    data->uiRange = 2 << 7;
    data->bitsLeft -= 7;
  }
  else if ( data->uiRange >= 256 )
  {
    return;
  }
  else
  {
    data->uiLow   <<= 1;
    data->uiRange <<= 1;
    data->bitsLeft--;
  }
  
  if(data->bitsLeft < 12)
  {
    cabac_write(data);
  }
}

void cabac_flush(cabac_data* data)
{
  cabac_encodeBinTrm(data,1);
  cabac_finish(data);
  bitstream_put(data->stream,1,1);

  cabac_start(data);
}

void cabac_encodeBinEP(cabac_data* data, uint32_t binValue )
{
  data->uiBinsCoded += data->binCountIncrement;
  data->uiLow <<= 1;
  if( binValue )
  {
    data->uiLow += data->uiRange;
  }
  data->bitsLeft--;

  if(data->bitsLeft < 12)
  {
    cabac_write(data);
  }
}

void cabac_encodeBinsEP(cabac_data* data, uint32_t binValues, int numBins )
{
  uint32_t pattern;
  data->uiBinsCoded += numBins & -data->binCountIncrement;

  while ( numBins > 8 )
  {
    numBins -= 8;
    pattern = binValues >> numBins;
    data->uiLow <<= 8;
    data->uiLow += data->uiRange * pattern;
    binValues -= pattern << numBins;
    data->bitsLeft -= 8;
    
    if(data->bitsLeft < 12)
    {
      cabac_write(data);
    }
  }
  
  data->uiLow <<= numBins;
  data->uiLow += data->uiRange * binValues;
  data->bitsLeft -= numBins;
  
  if(data->bitsLeft < 12)
  {
    cabac_write(data);
  }
}