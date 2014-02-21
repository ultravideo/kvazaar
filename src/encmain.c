/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2014 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation.
 *
 * Kvazaar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

/*
 * \file
 *
 */

#ifdef _WIN32
/* The following two defines must be located before the inclusion of any system header files. */
#define WINVER       0x0500
#define _WIN32_WINNT 0x0500
#include <io.h>       /* _setmode() */
#include <fcntl.h>    /* _O_BINARY */
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

// Assembly optimization headers
#include "x86/cpu.h"

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
  FILE *cqmfile = NULL; //!< HM-compatible CQM file
  encoder_control *encoder = NULL; //!< Encoder control struct
  double psnr[3] = { 0.0, 0.0, 0.0 };
  uint64_t curpos  = 0;
  uint64_t lastpos = 0;
  #ifdef _DEBUG
  FILE *recout = fopen("encrec_832x480_60.yuv","wb"); //!< reconstructed YUV output (only on debug mode)
  #endif

  // Stdin and stdout need to be binary for input and output to work.
  // Stderr needs to be text mode to convert \n to \r\n in Windows.
  #ifdef _WIN32
      _setmode( _fileno( stdin ),  _O_BINARY );
      _setmode( _fileno( stdout ), _O_BINARY );
      _setmode( _fileno( stderr ), _O_TEXT );
  #endif

  // Handle configuration
  cfg = config_alloc();

  // If problem with configuration, print banner and shutdown
  if (!cfg || !config_init(cfg) || !config_read(cfg,argc,argv)) {
    fprintf(stderr,
            "/***********************************************/\n"
            " *   Kvazaar HEVC Encoder v. " VERSION_STRING "             *\n"
            " *     Tampere University of Technology 2014   *\n"
            "/***********************************************/\n\n");

    fprintf(stderr,
            "Usage:\n"
            "kvazaar -i <input> -w <width> -h <height> -o <output>\n"
            "\n"
            "Optional parameters:\n"
            "      -n, --frames <integer>     : number of frames to code [all]\n"
            "      -q, --qp <integer>         : Quantization Parameter [32]\n"
            "      -p, --period <integer>     : Period of intra pictures [0]\n"
            "                                     0: only first picture is intra\n"
            "                                     1: all pictures are intra\n"
            "                                     2-N: every Nth picture is intra\n"
            "      -r, --ref <integer>        : Reference frames, range 1..15 [3]\n"
            "          --no-deblock           : Disable deblocking filter\n"
            "          --deblock <beta:tc>    : Deblocking filter parameters\n"
            "                                   beta and tc range is -6..6 [0:0]\n"
            "          --no-sao               : Disable sample adaptive offset\n"
            "          --aud                  : Use access unit delimiters\n"
            "          --cqmfile <string>     : Custom Quantization Matrices from a file\n"
            "\n"
            "  Video Usability Information:\n"
            "          --sar <width:height>   : Specify Sample Aspect Ratio\n"
            "          --overscan <string>    : Specify crop overscan setting [\"undef\"]\n"
            "                                     - undef, show, crop\n"
            "          --videoformat <string> : Specify video format [\"undef\"]\n"
            "                                     - component, pal, ntsc, secam, mac, undef\n"
            "          --range <string>       : Specify color range [\"tv\"]\n"
            "                                     - tv, pc\n"
            "          --colorprim <string>   : Specify color primaries [\"undef\"]\n"
            "                                     - undef, bt709, bt470m, bt470bg,\n"
            "                                       smpte170m, smpte240m, film, bt2020\n"
            "          --transfer <string>    : Specify transfer characteristics [\"undef\"]\n"
            "                                     - undef, bt709, bt470m, bt470bg,\n"
            "                                       smpte170m, smpte240m, linear, log100,\n"
            "                                       log316, iec61966-2-4, bt1361e,\n"
            "                                       iec61966-2-1, bt2020-10, bt2020-12\n"
            "          --colormatrix <string> : Specify color matrix setting [\"undef\"]\n"
            "                                     - undef, bt709, fcc, bt470bg, smpte170m,\n"
            "                                       smpte240m, GBR, YCgCo, bt2020nc, bt2020c\n"
            "          --chromaloc <integer>  : Specify chroma sample location (0 to 5) [0]\n");

    if (cfg)
      config_destroy(cfg);

    return EXIT_FAILURE;
  }

  // Dig CPU features with cpuid
  kvz_cpu_cpuid(&ecx,&edx);
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

  // Check if the input file name is a dash, this means stdin
  if (!strcmp(cfg->input, "-")) {
    input = stdin;
  } else {
    // Otherwise we try to open the input file
    input = fopen(cfg->input, "rb");
  }

  // Check that input was opened correctly
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

  encoder = init_encoder_control(cfg);
  if (!encoder)
    return EXIT_FAILURE;

  // Set output file
  encoder->output = output;

  // input init (TODO: read from commandline / config)
  encoder->bitdepth = 8;
  encoder->frame    = 0;
  encoder->QP       = encoder->cfg->qp;
  encoder->in.video_format = FORMAT_420;
  // deblocking filter
  encoder->deblock_enable   = encoder->cfg->deblock_enable;
  encoder->beta_offset_div2 = encoder->cfg->deblock_beta;
  encoder->tc_offset_div2   = encoder->cfg->deblock_tc;
  // SAO
  encoder->sao_enable = encoder->cfg->sao_enable;
  // VUI
  encoder->vui.sar_width   = encoder->cfg->vui.sar_width;
  encoder->vui.sar_height  = encoder->cfg->vui.sar_height;
  encoder->vui.overscan    = encoder->cfg->vui.overscan;
  encoder->vui.videoformat = encoder->cfg->vui.videoformat;
  encoder->vui.fullrange   = encoder->cfg->vui.fullrange;
  encoder->vui.colorprim   = encoder->cfg->vui.colorprim;
  encoder->vui.transfer    = encoder->cfg->vui.transfer;
  encoder->vui.colormatrix = encoder->cfg->vui.colormatrix;
  encoder->vui.chroma_loc  = encoder->cfg->vui.chroma_loc;
  // AUD
  encoder->aud_enable = encoder->cfg->aud_enable;
  // CQM
  cqmfile = cfg->cqmfile ? fopen(cfg->cqmfile, "rb") : NULL;
  encoder->cqmfile = cqmfile;

  init_encoder_input(&encoder->in, input, cfg->width, cfg->height);

  printf("Input: %s, output: %s\n", cfg->input, cfg->output);
  printf("  Video size: %dx%d\n (%dx%d)",
         encoder->in.width, encoder->in.height,
         encoder->in.real_width, encoder->in.real_height);

  // Only the code that handles conformance window coding needs to know
  // the real dimensions. As a quick fix for broken non-multiple of 8 videos,
  // change the input values here to be the real values. For a real fix
  // encoder.in probably needs to be merged into cfg.
  cfg->width = encoder->in.width;
  cfg->height = encoder->in.height;

  // Init coeff data table
  encoder->in.cur_pic->coeff_y = MALLOC(coefficient, cfg->width * cfg->height);
  encoder->in.cur_pic->coeff_u = MALLOC(coefficient, (cfg->width * cfg->height) >> 2);
  encoder->in.cur_pic->coeff_v = MALLOC(coefficient, (cfg->width * cfg->height) >> 2);

  // Init predicted data table
  encoder->in.cur_pic->pred_y = MALLOC(pixel, cfg->width * cfg->height);
  encoder->in.cur_pic->pred_u = MALLOC(pixel, (cfg->width * cfg->height) >> 2);
  encoder->in.cur_pic->pred_v = MALLOC(pixel, (cfg->width * cfg->height) >> 2);

  // Start coding cycle while data on input and not on the last frame
  while(!cfg->frames || encoder->frame < cfg->frames) {
    int32_t diff;
    double temp_psnr[3];

    // Read one frame from the input
    if (!read_one_frame(input, encoder)) {
      if (!feof(input))
        printf("Failed to read a frame %d\n", encoder->frame);
      break;
    }

    // The actual coding happens here, after this function we have a coded frame
    encode_one_frame(encoder);

    #ifdef _DEBUG
    // Write reconstructed frame out (for debugging purposes)
    fwrite(encoder->in.cur_pic->y_recdata, cfg->width * cfg->height, 1, recout);
    fwrite(encoder->in.cur_pic->u_recdata, (cfg->width * cfg->height)>>2, 1, recout);
    fwrite(encoder->in.cur_pic->v_recdata, (cfg->width * cfg->height)>>2, 1, recout);
    #endif

    // Calculate the bytes pushed to output for this frame
    fgetpos(output,(fpos_t*)&curpos);
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
    if (encoder->ref->used_size == encoder->cfg->ref_frames) {
	    picture_list_rem(encoder->ref, encoder->ref->used_size-1, 1);
    }
    // Add current picture as reference
    picture_list_add(encoder->ref, encoder->in.cur_pic);
    // Allocate new memory to current picture
    // TODO: reuse memory from old reference
    encoder->in.cur_pic = picture_init(encoder->in.width, encoder->in.height, encoder->in.width_in_lcu, encoder->in.height_in_lcu);

    // Copy pointer from the last cur_pic because we don't want to reallocate it
    MOVE_POINTER(encoder->in.cur_pic->coeff_y,encoder->ref->pics[0]->coeff_y);
    MOVE_POINTER(encoder->in.cur_pic->coeff_u,encoder->ref->pics[0]->coeff_u);
    MOVE_POINTER(encoder->in.cur_pic->coeff_v,encoder->ref->pics[0]->coeff_v);

    MOVE_POINTER(encoder->in.cur_pic->pred_y,encoder->ref->pics[0]->pred_y);
    MOVE_POINTER(encoder->in.cur_pic->pred_u,encoder->ref->pics[0]->pred_u);
    MOVE_POINTER(encoder->in.cur_pic->pred_v,encoder->ref->pics[0]->pred_v);

    encoder->frame++;
    encoder->poc++;
  }
  // Coding finished
  fgetpos(output,(fpos_t*)&curpos);

  // Print statistics of the coding
  printf(" Processed %d frames, %10d bits AVG PSNR: %2.4f %2.4f %2.4f\n", encoder->frame, ((int32_t)curpos)<<3,
         psnr[0] / encoder->frame, psnr[1] / encoder->frame, psnr[2] / encoder->frame);

  fclose(input);
  fclose(output);
  if(cqmfile != NULL) fclose(cqmfile);
  #ifdef _DEBUG
  fclose(recout);
  #endif

  // Deallocating
  config_destroy(cfg);
  scalinglist_destroy();
  picture_list_destroy(encoder->ref);
  picture_destroy(encoder->in.cur_pic);
  FREE_POINTER(encoder->in.cur_pic);
  bitstream_free(encoder->stream);
  FREE_POINTER(encoder->stream);
  free(encoder);
  free_tables();
  FREE_POINTER(g_exp_table);

  return EXIT_SUCCESS;
}
