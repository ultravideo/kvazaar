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
 
 
 /*!
     \brief Program main function.
     \param argc Argument count from commandline
     \param argv Argument list
     \return Program exit state
 */
  int main(int argc, char* argv[])
  {

    config *cfg  = NULL;       /* Global configuration */
    FILE *input  = NULL;
    FILE *output = NULL;
    encoder_control* encoder = (encoder_control*)malloc(sizeof(encoder_control));;
 
    /* Handle configuration */
    cfg = config_alloc();
    
    /* If problem with configuration, shutdown */
    if(!config_init(cfg) || !config_read(cfg,argc,argv))
    {
      fprintf(stderr, "/////////////////////////////////////////////////\r\n");
      fprintf(stderr, "//           HEVC Encoder v. " VERSION_STRING "//\r\n");
      fprintf(stderr, "//      Tampere University of Technology  2012 //\r\n");
      fprintf(stderr, "/////////////////////////////////////////////////\r\n\r\n");
      
      fprintf(stderr, "Usage:\r\n");
      fprintf(stderr, "encmain -i <input> -w <width> -h <height> -o <output>\r\n");
      fprintf(stderr, "Optional parameters:\r\n");
      fprintf(stderr, "      -n <frames> : number of frames to decode\r\n");
      fprintf(stderr, "      -s <frames> : number of frames to skip from the beginning\r\n");

      config_destroy(cfg);
      return EXIT_FAILURE;
    }

	  printf("Input: %s, output: %s\r\n", cfg->input, cfg->output);
    printf("  Video size: %dx%d\r\n", cfg->width, cfg->height);

    /* Open input file and check that it was opened correctly */
    input = fopen(cfg->input, "rb");
    if(input == NULL)
    {
      fprintf(stderr, "Couldn't open input file!\r\n");
      config_destroy(cfg);
      return EXIT_FAILURE;
    }
    
    /* Open output file and check that it was opened correctly */
    output = fopen(cfg->output, "wb");
    if(output == NULL)
    {
      fprintf(stderr, "Couldn't open output file!\r\n");
      config_destroy(cfg);
      return EXIT_FAILURE;
    }

    /* Initialization */
    cabac_init(&cabac);
    //ToDo: add bitstream
    //cabac.stream = 
    init_encoder_control(encoder, output);
    init_encoder_input(&encoder->in, input, 320, 240);

    fclose(input);
    fclose(output);

    config_destroy(cfg);

    return EXIT_SUCCESS;
  }