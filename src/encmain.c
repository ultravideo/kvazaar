/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
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
 #include "transform.h"
 
 /* Assembly optimizations */
#ifndef X64
 #include "x86/test.h"
#else
 #include "x64/test64.h"
#endif
 
 /*!
     \brief Program main function.
     \param argc Argument count from commandline
     \param argv Argument list
     \return Program exit state
 */
  int main(int argc, char* argv[])
  {    
    int ecx = 0,edx =0;
    /* CPU feature bits */
    enum { BIT_SSE3 = 0,BIT_SSSE3 = 9, BIT_SSE41 = 19, BIT_SSE42 = 20, BIT_MMX = 24, BIT_SSE = 25, BIT_SSE2 = 26, BIT_AVX = 28};
    uint32_t curFrame = 0;
    config *cfg  = NULL;       /* Global configuration */
    FILE *input  = NULL;
    FILE *output = NULL;
    double PSNR[3] = { 0.0, 0.0, 0.0 };
    fpos_t curpos = 0;
    fpos_t lastpos = 0;
    #ifdef _DEBUG
    FILE *recout = fopen("encrec.yuv","wb");
    #endif
    encoder_control* encoder = (encoder_control*)malloc(sizeof(encoder_control));


 
    /* Handle configuration */
    cfg = config_alloc();
    
    /* If problem with configuration, shutdown */
    if(!config_init(cfg) || !config_read(cfg,argc,argv))
    {
      fprintf(stderr, "/***********************************************/\r\n");
      fprintf(stderr, " *           HEVC Encoder v. " VERSION_STRING "*\r\n");
      fprintf(stderr, " *     Tampere University of Technology  2013  *\r\n");
      fprintf(stderr, "/***********************************************/\r\n\r\n");
      
      fprintf(stderr, "Usage:\r\n");
      fprintf(stderr, "encmain -i <input> -w <width> -h <height> -o <output>\r\n");
      fprintf(stderr, "Optional parameters:\r\n");
      fprintf(stderr, "      -n <frames> : number of frames to decode\r\n");
      fprintf(stderr, "      -s <frames> : number of frames to skip from the beginning\r\n");

      config_destroy(cfg);
      return EXIT_FAILURE;
    }
    /* CPU id */
    
    #ifndef X64
    //cpuId32(&ecx,&edx);
    #else
    cpuId64(&ecx,&edx);
    #endif
    printf("CPU features enabled: ");
    /* EDX */
    if(edx & (1<<BIT_MMX))  printf("MMX ");
    if(edx & (1<<BIT_SSE))  printf("SSE ");
    if(edx & (1<<BIT_SSE2)) printf("SSE2 ");
    /* ECX */
    if(ecx & (1<<BIT_SSE3))  printf("SSE3 ");
    if(ecx & (1<<BIT_SSSE3)) printf("SSSE3 ");
    if(ecx & (1<<BIT_SSE41)) printf("SSE4.1 ");
    if(ecx & (1<<BIT_SSE42)) printf("SSE4.2 ");
    if(ecx & (1<<BIT_AVX))   printf("AVX ");
    printf("\r\n");
    

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
    init_tables();
    init_exp_golomb(4096*8);
    cabac_init(&cabac);
    scalinglist_init();
    init_encoder_control(encoder, (bitstream*)malloc(sizeof(bitstream))); 
    encoder->ref = picture_list_init(MAX_REF_PIC_COUNT);
    
    /* Init bitstream */
    bitstream_init(encoder->stream);
    encoder->stream->buffer_pos = 0;
    encoder->stream->output = 0;
    /* Alloc 1MB */
    bitstream_alloc(encoder->stream, 1024*2*cfg->width);

    /* Config pointer to encoder struct */
    encoder->cfg = cfg;
    /* Set output file */
    encoder->output = output;
    /* Set CABAC output bitstream */
    cabac.stream = encoder->stream;

    /* input init (TODO: read from commandline / config) */
    encoder->bitdepth = 8;
    encoder->frame    = 0;
    encoder->QP       = 32;
    encoder->in.video_format = FORMAT_420;
    /* deblocking */
    encoder->deblock_enable  = 1;
    encoder->betaOffsetdiv2  = 0;
    encoder->tcOffsetdiv2    = 0;
    /* SAO */
    encoder->sao_enable = 0;

    init_encoder_input(&encoder->in, input, cfg->width, cfg->height);

    /* Start coding cycle */
    while(!feof(input) && (!cfg->frames || encoder->frame < cfg->frames))
    {
      /* Read one frame from the input */
      read_one_frame(input, encoder);

      /* Clear reconstruction buffers (not needed, for debugging) */
      /*
      memset(encoder->in.cur_pic->yRecData, 0, cfg->width*cfg->height);
      memset(encoder->in.cur_pic->uRecData, 128, cfg->width*cfg->height>>2);
      memset(encoder->in.cur_pic->vRecData, 128, cfg->width*cfg->height>>2);
      */
      /* /////////////THE ACTUAL CODING HAPPENDS HERE\\\\\\\\\\\\\\\\\\\ */
      encode_one_frame(encoder);
      /* ////////////CODING NOW DONE\\\\\\\\\\\\\\\\\ */

      #ifdef _DEBUG
      /* Write reconstructed frame out */     
      fwrite(encoder->in.cur_pic->yRecData,cfg->width*cfg->height,1,recout);
      fwrite(encoder->in.cur_pic->uRecData,cfg->width*cfg->height>>2,1,recout);
      fwrite(encoder->in.cur_pic->vRecData,cfg->width*cfg->height>>2,1,recout);
      #endif
      {
        int32_t diff;
        double temp_PSNR[3];
        fgetpos(output,&curpos);
        diff = (int32_t)(curpos-lastpos);
        lastpos = curpos;

        temp_PSNR[0] = imagePSNR(encoder->in.cur_pic->yData,encoder->in.cur_pic->yRecData,cfg->width,cfg->height);
        temp_PSNR[1] = imagePSNR(encoder->in.cur_pic->uData,encoder->in.cur_pic->uRecData,cfg->width>>1,cfg->height>>1);
        temp_PSNR[2] = imagePSNR(encoder->in.cur_pic->vData,encoder->in.cur_pic->vRecData,cfg->width>>1,cfg->height>>1);
        
        printf("POC %4d (%c-frame) %10d bits PSNR: %2.4f %2.4f %2.4f\n", encoder->frame, "BPI"[encoder->in.cur_pic->slicetype%3],diff<<3,
                                                        temp_PSNR[0],temp_PSNR[1],temp_PSNR[2]);
        PSNR[0]+=temp_PSNR[0];
        PSNR[1]+=temp_PSNR[1];
        PSNR[2]+=temp_PSNR[2];
      }
      encoder->frame++;
    }
    /* Coding finished */
    fgetpos(output,&curpos);

    printf(" Processed %d frames, %10d bits AVG PSNR: %2.4f %2.4f %2.4f\n", encoder->frame,((int32_t)curpos)<<3,PSNR[0]/encoder->frame,PSNR[1]/encoder->frame,PSNR[2]/encoder->frame);

    fclose(input);
    fclose(output);
    #ifdef _DEBUG
    fclose(recout);
    #endif

    /* Deallocating */
    config_destroy(cfg);
    scalinglist_destroy();

    return EXIT_SUCCESS;
  }