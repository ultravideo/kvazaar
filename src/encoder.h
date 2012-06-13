/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
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

//ToDo: add ME data
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
  uint32_t width;
  uint32_t height;
  uint32_t height_in_LCU;
  uint32_t width_in_LCU;
  picture cur_pic;
  uint8_t video_format;
} encoder_input;

typedef struct
{
  uint32_t frame;
  config *cfg;
  encoder_input in;
  encoder_me me;
  bitstream* stream;
  FILE *output;
  picture_list *ref;
  uint8_t QP;
} encoder_control;

void init_encoder_control(encoder_control* control,bitstream* output);
void init_encoder_input(encoder_input* input,FILE* inputfile, uint32_t width, uint32_t height);
void encode_one_frame(encoder_control* encoder);


void encode_seq_parameter_set(encoder_control* encoder);
void encode_pic_parameter_set(encoder_control* encoder);
void encode_slice_data(encoder_control* encoder);
void encode_slice_header(encoder_control* encoder);
void encode_coding_tree(encoder_control* encoder,uint16_t xCtb,uint16_t yCtb, uint8_t depth);


static const uint8_t  INIT_SPLIT_FLAG[3][3] =  
                       { { 107,  139,  126 }, { 107,  139,  126 },  { 139,  141,  157 } };

static const uint8_t INIT_INTRA_PRED_MODE[3] = { 183,154,184 };

static const uint8_t INIT_CHROMA_PRED_MODE[3][2] = { { 152,  139 }, { 152,  139 }, {  63,  139 } };

#define CNU 154
static const uint8_t INIT_TRANS_SUBDIV_FLAG[3][4] = 
{
  { CNU,  153,  138,  138 }, 
  { CNU,  124,  138,   94 }, 
  { CNU,  224,  167,  122 }
};

static const uint8_t INIT_QT_CBF[3][6] =  
{
  { 153,  111,  CNU,  149,   92,  167 },
  { 153,  111,  CNU,  149,  107,  167 },
  { 111,  141,  CNU,   94,  138,  182 }
};

static const uint8_t INIT_SIG_CG_FLAG[3][4] =  
 {  { 121,  140,  61,  154  },  { 121,  140,  61,  154 }, {  91,  171,  134,  141  } };

static const uint8_t INIT_SIG_FLAG[3][45] = 
{{170,154,139,153,139,123,123, 63,124,153,153,152,152,152,137,152,137,137,166,183,140,136,153,154,170,153,138,138,122,121,122,121,167,153,167,136,121,122,136,121,122,91,151,183,140,},
 {155,154,139,153,139,123,123, 63,153,153,153,152,152,152,137,152,137,122,166,183,140,136,153,154,170,153,123,123,107,121,107,121,167,153,167,136,149,107,136,121,122,91,151,183,140,},
 {111,111,125,110,110, 94,124,108,124,139,139,139,168,124,138,124,138,107,107,125,141,179,153,125,140,139,182,182,152,136,152,136,153,182,137,149,192,152,224,136,31,136,136,139,111,} };

static const uint8_t INIT_LAST[3][30] =  
{
  { 125,  110,  124,  110,   95,   94,  125,  111,  111,   79,  125,  126,  111,  111,   79,
    108,  123,   93,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU  }, 
  { 125,  110,   94,  110,   95,   79,  125,  111,  110,   78,  110,  111,  111,   95,   94,
    108,  123,  108,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU  }, 
  { 110,  110,  124,  125,  140,  153,  125,  127,  140,  109,  111,  143,  127,  111,   79,
    108,  123,   63,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU,  CNU  }
};

static const uint8_t g_uiGroupIdx[ 32 ]   = {0,1,2,3,4,4,5,5,6,6,6,6,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9};
static const uint8_t g_uiMinInGroup[ 10 ] = {0,1,2,3,4,6,8,12,16,24};

#endif