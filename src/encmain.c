/**
 * \file
 * \brief User interface for the encoder.
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 *
 * TODO: Check that these usage instructions are correct.
 * \subsection options_subsec All program options:
 *            - -i <filename>: input
 *            - -o <filename>: output
 *            - -w <width>: frame width
 *            - -h <height>: frame height
 *            - -n <n>: encode only n frames
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "config.h"
#include "encoder.h"
#include "cabac.h"
#include "picture.h"
#include "transform.h"
 
// Assembly optimization headers
#ifdef X86_64
  #include "x64/test64.h" 
#else
  #include "x86/test.h" 
#endif
 
/**
 * \brief Program main function.
 * \param argc Argument count from commandline
 * \param argv Argument list
 * \return Program exit state
 */
int main(int argc, char *argv[])
{    
  int ecx = 0,edx =0;
  /* CPU feature bits */
  enum { BIT_SSE3 = 0,BIT_SSSE3 = 9, BIT_SSE41 = 19, BIT_SSE42 = 20, BIT_MMX = 24, BIT_SSE = 25, BIT_SSE2 = 26, BIT_AVX = 28};
  uint32_t cur_frame = 0;
  config *cfg  = NULL; //!< Global configuration
  FILE *input  = NULL; //!< input file (YUV)
  FILE *output = NULL; //!< output file (HEVC NAL stream)
  double psnr[3] = { 0.0, 0.0, 0.0 };
  fpos_t curpos  = 0;
  fpos_t lastpos = 0;
  #ifdef _DEBUG
  FILE *recout = fopen("encrec.yuv","wb"); //!< reconstructed YUV output (only on debug mode)
  #endif
  encoder_control *encoder = (encoder_control*)malloc(sizeof(encoder_control));
    
  // Handle configuration
  cfg = config_alloc();
  
  // If problem with configuration, print banner and shutdown
  if (!config_init(cfg) || !config_read(cfg,argc,argv)) {
    fprintf(stderr, "/***********************************************/\r\n");
    fprintf(stderr, " *           HEVC Encoder v. " VERSION_STRING "*\r\n");
    fprintf(stderr, " *     Tampere University of Technology  2013  *\r\n");
    fprintf(stderr, "/***********************************************/\r\n\r\n");
      
    fprintf(stderr, "Usage:\r\n");
    fprintf(stderr, "hevc_encoder -i <input> -w <width> -h <height> -o <output>\r\n");
    fprintf(stderr, "Optional parameters:\r\n");
    fprintf(stderr, "      -n <frames> : number of frames to decode\r\n");
    fprintf(stderr, "      -s <frames> : number of frames to skip from the beginning\r\n");

    config_destroy(cfg);
    return EXIT_FAILURE;
  }

  // Dig CPU features with cpuid
  #ifdef X86_64
  cpuId64(&ecx,&edx);
  #else
  cpuId32(&ecx,&edx);
  #endif
  printf("CPU features enabled: ");
  // EDX
  if (edx & (1<<BIT_MMX))  printf("MMX ");
  if (edx & (1<<BIT_SSE))  printf("SSE ");
  if (edx & (1<<BIT_SSE2)) printf("SSE2 ");
  // ECX
  if (ecx & (1<<BIT_SSE3))  printf("SSE3 ");
  if (ecx & (1<<BIT_SSSE3)) printf("SSSE3 ");
  if (ecx & (1<<BIT_SSE41)) printf("SSE4.1 ");
  if (ecx & (1<<BIT_SSE42)) printf("SSE4.2 ");
  if (ecx & (1<<BIT_AVX))   printf("AVX ");
  printf("\r\n");
    

	printf("Input: %s, output: %s\n", cfg->input, cfg->output);
  printf("  Video size: %dx%d\n", cfg->width, cfg->height);

  // Open input file and check that it was opened correctly
  input = fopen(cfg->input, "rb");
  if (input == NULL) {
    fprintf(stderr, "Could not open input file, shutting down!\n");
    config_destroy(cfg);
    return EXIT_FAILURE;
  }
    
  // Open output file and check that it was opened correctly
  output = fopen(cfg->output, "wb");
  if (output == NULL) {
    fprintf(stderr, "Could not open output file, shutting down!\n");
    config_destroy(cfg);
    return EXIT_FAILURE;
  }

  // Initialization
  init_tables();
  init_exp_golomb(4096*8); //Allocate and init exp golomb table
  scalinglist_init();
  init_encoder_control(encoder, (bitstream*)malloc(sizeof(bitstream))); 
  encoder->ref = picture_list_init(MAX_REF_PIC_COUNT);
    
  // Init bitstream
  bitstream_init(encoder->stream);
  encoder->stream->buffer_pos = 0;
  encoder->stream->output = 0;

  // Alloc 2kB*width for bitstream buffer (for one coded frame)
  bitstream_alloc(encoder->stream, 1024*2*cfg->width);

  // Config pointer to encoder struct
  encoder->cfg = cfg;
  // Set output file
  encoder->output = output;
  // Set CABAC output bitstream
  cabac.stream = encoder->stream;

  // input init (TODO: read from commandline / config)
  encoder->bitdepth = 8;
  encoder->frame    = 0;
  encoder->QP       = 32;
  encoder->in.video_format = FORMAT_420;
  // deblocking filter
  encoder->deblock_enable  = 1;
  encoder->beta_offset_div2  = 0;
  encoder->tc_offset_div2    = 0;
  // SAO
  encoder->sao_enable = 0;

  init_encoder_input(&encoder->in, input, cfg->width, cfg->height);

  // Start coding cycle while data on input and not on the last frame
  while(!feof(input) && (!cfg->frames || encoder->frame < cfg->frames)) {
    int32_t diff;
    double temp_psnr[3];

    // Read one frame from the input
    read_one_frame(input, encoder);

    // The actual coding happens here, after this function we have a coded frame
    encode_one_frame(encoder);    

    #ifdef _DEBUG
    // Write reconstructed frame out (for debugging purposes)
    fwrite(encoder->in.cur_pic->y_recdata, cfg->width * cfg->height, 1, recout);
    fwrite(encoder->in.cur_pic->u_recdata, (cfg->width * cfg->height)>>2, 1, recout);
    fwrite(encoder->in.cur_pic->v_recdata, (cfg->width * cfg->height)>>2, 1, recout);
    #endif

    // Calculate the bytes pushed to output for this frame
    fgetpos(output,&curpos);
    diff = (int32_t)(curpos-lastpos);
    lastpos = curpos;

    // PSNR calculations
    temp_psnr[0] = image_psnr(encoder->in.cur_pic->y_data, encoder->in.cur_pic->y_recdata, cfg->width, cfg->height);
    temp_psnr[1] = image_psnr(encoder->in.cur_pic->u_data, encoder->in.cur_pic->u_recdata, cfg->width>>1, cfg->height>>1);
    temp_psnr[2] = image_psnr(encoder->in.cur_pic->v_data, encoder->in.cur_pic->v_recdata, cfg->width>>1, cfg->height>>1);

    printf("POC %4d (%c-frame) %10d bits PSNR: %2.4f %2.4f %2.4f\n", encoder->frame,
           "BPI"[encoder->in.cur_pic->slicetype%3], diff<<3,
           temp_psnr[0], temp_psnr[1], temp_psnr[2]);

    // Increment total PSNR
    psnr[0] += temp_psnr[0];
    psnr[1] += temp_psnr[1];
    psnr[2] += temp_psnr[2];


    // TODO: add more than one reference

    // Remove the ref pic (if present)
    picture_list_rem(encoder->ref, 0, 1);
    // Add current picture as reference
    picture_list_add(encoder->ref, encoder->in.cur_pic);
    // Allocate new memory to current picture
    // TODO: reuse memory from old reference
    encoder->in.cur_pic = picture_init(encoder->in.width, encoder->in.height, encoder->in.width_in_lcu, encoder->in.height_in_lcu);

    encoder->frame++;
  }
  // Coding finished
  fgetpos(output,&curpos);

  // Print statistics of the coding
  printf(" Processed %d frames, %10d bits AVG PSNR: %2.4f %2.4f %2.4f\n", encoder->frame, ((int32_t)curpos)<<3,
         psnr[0] / encoder->frame, psnr[1] / encoder->frame, psnr[2] / encoder->frame);

  fclose(input);
  fclose(output);
  #ifdef _DEBUG
  fclose(recout);
  #endif

  // Deallocating
  config_destroy(cfg);
  scalinglist_destroy();
  picture_list_destroy(encoder->ref);

  return EXIT_SUCCESS;
}
