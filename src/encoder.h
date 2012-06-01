/*! \file encoder.h
    \brief Encoding related functions
    \author Marko Viitanen
    \date 2012-06
    
    Structures for encoding
*/
#ifndef _ENCODER_H
#define _ENCODER_H


typedef struct
{
   void (*IME)(encoder_control* encoder);
 
} encoder_me;


typedef struct
{
    FILE *file;

 
} encoder_input;

typedef struct
{
  encoder_input in;
  encoder_me me;
  FILE *output; 
} encoder_control;



#endif