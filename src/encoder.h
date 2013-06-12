/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file encoder.h
    \brief Encoding related functions
    \author Marko Viitanen
    \date 2012-06
    
    Structures for encoding
*/
#ifndef __ENCODER_H
#define __ENCODER_H

#include "picture.h"
#include "bitstream.h"

typedef struct encoder_control;

/* ToDo: add ME data */
typedef struct
{
  void (*IME)();
  void (*FME)();
  int range;
 
} encoder_me;

enum { FORMAT_400 = 0, FORMAT_420, FORMAT_422, FORMAT_444 };

/* Input info struct */
typedef struct
{
  FILE* file;
  int32_t width;  /*!< \brief input picture width */
  int32_t height; /*!< \brief input picture height */
  int32_t height_in_LCU; /*!< \brief input picture width in LCU*/
  int32_t width_in_LCU;  /*!< \brief input picture height in LCU */
  picture cur_pic;
  int8_t video_format;
  int8_t bitdepth;  /*!< \brief input bit depth (8,10) */
} encoder_input;

/* Encoder control options, the main struct */
typedef struct
{
  int32_t frame;
  config *cfg;
  encoder_input in;
  encoder_me me;
  bitstream* stream;
  FILE *output;
  picture_list *ref;
  int8_t QP;
  int8_t bitdepth;

  /* Filtering */
  int8_t deblock_enable; /*!< \brief Flag to enable deblocking filter */
  int8_t sao_enable;     /*!< \brief Flag to enable sample adaptive offset filter */
  int8_t betaOffsetdiv2; /*!< \brief (deblocking) beta offset (div 2), range -6...6 */
  int8_t tcOffsetdiv2;   /*!< \brief (deblocking)tc offset (div 2), range -6...6 */
} encoder_control;

typedef struct
{
  int8_t idx;
  uint8_t *base;
  uint8_t *baseU;
  uint8_t *baseV;
  
  uint8_t *recbase;
  uint8_t *recbaseU;
  uint8_t *recbaseV;
  
  int16_t *pred;
  int16_t *predU;
  int16_t *predV;

  int32_t base_stride;
  int32_t recbase_stride;
  int32_t pred_stride;
  
  /* ToDo: unify luma+chroma arrays */
  int16_t *coeff[3];
  int8_t cb_top[3];
  int8_t cb[4];
  int8_t intraPredMode;
  int8_t intraPredModeChroma;
  int32_t split[4];

  int32_t xCtb,yCtb;

} transform_info;

void init_tables(void);
void init_encoder_control(encoder_control* control,bitstream* output);
void init_encoder_input(encoder_input* input,FILE* inputfile, int32_t width, int32_t height);
void encode_one_frame(encoder_control* encoder);


void encode_seq_parameter_set(encoder_control* encoder);
void encode_pic_parameter_set(encoder_control* encoder);
void encode_vid_parameter_set(encoder_control* encoder);
void encode_slice_data(encoder_control* encoder);
void encode_slice_header(encoder_control* encoder);
void encode_coding_tree(encoder_control* encoder,uint16_t xCtb,uint16_t yCtb, uint8_t depth);
void encode_lastSignificantXY(encoder_control* encoder,uint8_t lastpos_x, uint8_t lastpos_y, uint8_t width, uint8_t height, uint8_t type, uint8_t scan);
void encode_CoeffNxN(encoder_control* encoder,int16_t* coeff, uint8_t width, uint8_t type, int8_t scanMode);
void encode_transform_tree(encoder_control* encoder,transform_info* ti,uint8_t depth);
void encode_transform_coeff(encoder_control* encoder,transform_info* ti,int8_t depth, int8_t trDepth);

extern int16_t g_lambda_cost[55];
extern uint32_t* g_auiSigLastScan[3][7];
int8_t g_aucConvertToBit[LCU_WIDTH+1];
static int8_t g_bitDepth     = 8;
static int8_t g_uiBitIncrement = 0;

#define MAX_NUM_SPU_W ((1<<(MAX_DEPTH))/4)
static uint32_t g_auiZscanToRaster [ MAX_NUM_SPU_W*MAX_NUM_SPU_W ] = { 0, };
static uint32_t g_auiRasterToZscan [ MAX_NUM_SPU_W*MAX_NUM_SPU_W ] = { 0, };
static const uint8_t g_uiGroupIdx[ 32 ]    = {0,1,2,3,4,4,5,5,6,6,6,6,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9};
static const uint8_t g_uiMinInGroup[ 10 ]  = {0,1,2,3,4,6,8,12,16,24};
static uint32_t g_sigLastScanCG32x32[ 64 ] = 
{  0, 8, 1,16, 9, 2,24,17,
  10, 3,32,25,18,11, 4,40,
  33,26,19,12, 5,48,41,34,
  27,20,13, 6,56,49,42,35,
  28,21,14, 7,57,50,43,36,
  29,22,15,58,51,44,37,30,
  23,59,52,45,38,31,60,53,
  46,39,61,54,47,62,55,63 };

static const uint32_t g_sigLastScan8x8[ 3 ][ 4 ] =
{ {0, 2, 1, 3},
  {0, 1, 2, 3},
  {0, 2, 1, 3}
};

// 
//4 8 16 32 64 128
//0 1  2  3  4   5
static const uint8_t g_toBits[129] =
{  
  0,
  0,0,0,0,
  0,0,0,1,
  0,0,0,0,0,0,0,2,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5
};
#define TOBITS(len) g_toBits[len]


#define C1FLAG_NUMBER               8 /*!< maximum number of largerThan1 flag coded in one chunk */
#define C2FLAG_NUMBER               1 /*!< maximum number of largerThan2 flag coded in one chunk */

enum COEFF_SCAN_TYPE
{
  SCAN_DIAG = 0,         /*!< up-right diagonal scan */
  SCAN_HOR,              /*!< horizontal first scan  */
  SCAN_VER               /*!< vertical first scan    */
};


#endif