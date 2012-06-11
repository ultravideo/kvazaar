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

/* Input info struct */
typedef struct
{
  FILE* file;
  uint32_t width;
  uint32_t height;
  uint32_t height_in_LCU;
  uint32_t width_in_LCU;
  picture cur_pic;
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

static const uint8_t INIT_QT_CBF[3][8] =  
{
  { 153,  111,  CNU,  CNU,  CNU,  149,   92,  167 }, 
  { 153,  111,  CNU,  CNU,  CNU,  149,  107,  167 }, 
  { 111,  141,  CNU,  CNU,  CNU,   94,  138,  182 }
};



#endif