/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file encoder.c
    \brief Encoding related functions
    \author Marko Viitanen
    \date 2012-06
    
    Encoder main level
*/
/* Suppress some windows warnings */
#ifdef WIN32
  #define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "config.h"
#include "encoder.h"
#include "cabac.h"
#include "picture.h"

void init_encoder_control(encoder_control* control,bitstream* output) {control->stream = output;};
void init_encoder_input(encoder_input* input,FILE* inputfile, uint32_t width, uint32_t height)
{
  input->file = inputfile;
  input->width = width;
  input->height = height;

  input->cur_pic.width = width;
  input->cur_pic.height = height;
  input->cur_pic.referenced = 0;
  /* Allocate buffers */
  input->cur_pic.yData = (uint8_t *)malloc(width*height);
  input->cur_pic.uData = (uint8_t *)malloc((width*height)>>2);
  input->cur_pic.vData = (uint8_t *)malloc((width*height)>>2);
};


void encode_one_frame(encoder_control* encoder)
{


}