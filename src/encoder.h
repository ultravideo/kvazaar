/*! \file encoder.h
    \brief Encoding related functions
    \author Marko Viitanen
    \date 2012-06
    
    Structures for encoding
*/
#ifndef _ENCODER_H
#define _ENCODER_H

//ToDo: add ME data
typedef struct
{
  void (*IME)(encoder_control* encoder);
 
} encoder_me;


typedef struct
{
  FILE* file;
  uint32_t width;
  uint32_t height;
  uint32_t height_in_LCU;
  uint32_t width_in_LCU; 
} encoder_input;

typedef struct
{
  encoder_input in;
  encoder_me me;
  FILE* output; 
} encoder_control;

init_encoder_control(encoder_control* control,FILE* output) {control->output = output;};
init_encoder_input(encoder_input* input,FILE* inputfile, uint32_t width, uint32_t height) {input->file = inputfile; input->width = width; input->height = height;};

#endif