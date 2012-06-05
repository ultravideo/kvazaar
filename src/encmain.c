/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Computer Systems.
 */

/*! \file decmain.c
    \brief main file for the Decoder
    \author Marko Viitanen
    \date 2012-05
    
    This file contains main() function
*/

/*! \mainpage HEVC Encoder
 *
 * \section Coding style
 *
 * Coding style is explained in it's own document.
 *
 * \section usage_sec Usage
 *
 * \subsection encode_subsec Basic Decoding:
 * Use encmain.exe -i input.yuv -o output.hevc
 * 
 * \subsection options_subsec All program options:
 *            - -i <filename>: input
 *            - -o <filename>: output
 *            - -w <width>: frame width
 *            - -h <height>: frame height
 *            - -n <n>: encode only n frames
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
 
 
 /*!
     \brief Program main function.
     \param argc Argument count from commandline
     \param argv Argument list
     \return Program exit state
 */
  int main(int argc, char* argv[])
  {
    uint32_t curFrame = 0;
    config *cfg  = NULL;       /* Global configuration */
    FILE *input  = NULL;
    FILE *output = NULL;
    encoder_control* encoder = (encoder_control*)malloc(sizeof(encoder_control));;
 
    /* Handle configuration */
    cfg = config_alloc();
    
    /* If problem with configuration, shutdown */
    if(!config_init(cfg) || !config_read(cfg,argc,argv))
    {
      fprintf(stderr, "/***********************************************/\r\n");
      fprintf(stderr, " *           HEVC Encoder v. " VERSION_STRING "*\r\n");
      fprintf(stderr, " *     Tampere University of Technology  2012  *\r\n");
      fprintf(stderr, "/***********************************************/\r\n\r\n");
      
      fprintf(stderr, "Usage:\r\n");
      fprintf(stderr, "encmain -i <input> -w <width> -h <height> -o <output>\r\n");
      fprintf(stderr, "Optional parameters:\r\n");
      fprintf(stderr, "      -n <frames> : number of frames to decode\r\n");
      fprintf(stderr, "      -s <frames> : number of frames to skip from the beginning\r\n");

      config_destroy(cfg);
      return EXIT_FAILURE;
    }

	  printf("Input: %s, output: %s\n", cfg->input, cfg->output);
    printf("  Video size: %dx%d\n", cfg->width, cfg->height);

    /* Open input file and check that it was opened correctly */
    input = fopen(cfg->input, "rb");
    if(input == NULL)
    {
      fprintf(stderr, "Could not open input file, shutting down!\n");
      config_destroy(cfg);
      return EXIT_FAILURE;
    }
    
    /* Open output file and check that it was opened correctly */
    output = fopen(cfg->output, "wb");
    if(output == NULL)
    {
      fprintf(stderr, "Could not open output file, shutting down!\n");
      config_destroy(cfg);
      return EXIT_FAILURE;
    }

    /* Initialization */
    init_exp_golomb(4096*8);
    cabac_init(&cabac);
    init_encoder_control(encoder, (bitstream*)malloc(sizeof(bitstream)));

    /* Init bitstream */
    bitstream_init(encoder->stream);
    encoder->stream->buffer_pos = 0;
    encoder->stream->output = 0;
    bitstream_alloc(encoder->stream, 1024*1024);

    /* Config pointer to encoder struct */
    encoder->cfg = cfg;
    /* Set output file */
    encoder->output = output;
    /* Set CABAC output bitstream */
    cabac.stream = encoder->stream;

    /* input init */
    encoder->frame = 0;
    init_encoder_input(&encoder->in, input, cfg->width, cfg->height);

    /* Start coding cycle */
    while(!feof(input) && (!cfg->frames || curFrame < cfg->frames))
    {      
      /* Read one frame from the input */
      fread(encoder->in.cur_pic.yData, cfg->width*cfg->height,1,input);
      fread(encoder->in.cur_pic.uData, cfg->width*cfg->height/4,1,input);
      fread(encoder->in.cur_pic.vData, cfg->width*cfg->height/4,1,input);
      encode_one_frame(encoder);
      encoder->frame++;
    }
    /* Coding finished */

    printf(" Processed %d frames\n", encoder->frame-1);

    fclose(input);
    fclose(output);

    config_destroy(cfg);

    return EXIT_SUCCESS;
  }