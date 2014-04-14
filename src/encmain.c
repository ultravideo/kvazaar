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
#include <time.h>

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
  config *cfg  = NULL; //!< Global configuration
  FILE *input  = NULL; //!< input file (YUV)
  FILE *output = NULL; //!< output file (HEVC NAL stream)
  FILE *cqmfile = NULL; //!< HM-compatible CQM file
  encoder_control *encoder = NULL; //!< Encoder control struct
  double psnr[3] = { 0.0, 0.0, 0.0 };
  uint64_t curpos  = 0;
  uint64_t lastpos = 0;
  FILE *recout = NULL; //!< reconstructed YUV output, --debug
  clock_t start_time = clock();

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
            "kvazaar -i <input> --input-res <width>x<height> -o <output>\n"
            "\n"
            "Optional parameters:\n"
            "      -n, --frames <integer>     : Number of frames to code [all]\n"
            "      --seek <integer>           : First frame to code [0]\n"
            "      --input-res <int>x<int>    : Input resolution (width x height)\n"
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
            "          --no-rdoq              : Disable RDO quantization\n"
            "          --rd <integer>         : Rate-Distortion Optimization level [1]\n"
            "                                     0: no RDO\n"
            "                                     1: estimated RDO\n"
            "                                     2: full RDO\n"
            "          --no-transform-skip    : Disable transform skip\n"
            "          --aud                  : Use access unit delimiters\n"
            "          --cqmfile <string>     : Custom Quantization Matrices from a file\n"
            "          --debug <string>       : Output encoders reconstruction.\n"
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
            "          --chromaloc <integer>  : Specify chroma sample location (0 to 5) [0]\n"
            "\n"
            "  Deprecated parameters: (might be removed at some point)\n"
            "     Use --input-res:\n"
            "       -w, --width               : Width of input in pixels\n"
            "       -h, --height              : Height of input in pixels\n");

    if (cfg)
      config_destroy(cfg);

    return EXIT_FAILURE;
  }

  // Add dimensions to the reconstructions file name.
  if (cfg->debug != NULL) {
    char dim_str[50]; // log10(2^64) < 20, so this should suffice. I hate C.
    size_t left_len, right_len;
    sprintf(dim_str, "_%dx%d.yuv", cfg->width, cfg->height);
    left_len = strlen(cfg->debug);
    right_len = strlen(dim_str);
    cfg->debug = realloc(cfg->debug, left_len + right_len + 1);
    if (!cfg->debug) {
      fprintf(stderr, "realloc failed!\n");
      return EXIT_FAILURE;
    }
    strcpy(cfg->debug + left_len, dim_str);
  }

  // Do more validation to make sure the parameters we have make sense.
  if (!config_validate(cfg)) {
    config_destroy(cfg);
    return EXIT_FAILURE;
  }

  // Dig CPU features with cpuid
  kvz_cpu_cpuid(&ecx,&edx);
  fprintf(stderr, "CPU features enabled: ");
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
  fprintf(stderr, "\n");

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

  if (cfg->debug != NULL) {
    recout = fopen(cfg->debug, "wb");
    if (recout == NULL) {
      fprintf(stderr, "Could not open reconstruction file (%s), shutting down!\n", cfg->debug);
      config_destroy(cfg);
      return EXIT_FAILURE;
    }
  }

  encoder = init_encoder_control(cfg);
  if (!encoder)
    return EXIT_FAILURE;

  // Set output file

  encoder->output = output;
  ((bitstream_file*) encoder->stream)->output = output;
  // input init (TODO: read from commandline / config)
  encoder->bitdepth = 8;
  encoder->frame    = 0;
  encoder->QP       = (int8_t)encoder->cfg->qp;
  encoder->in.video_format = FORMAT_420;
  // deblocking filter
  encoder->deblock_enable   = (int8_t)encoder->cfg->deblock_enable;
  encoder->beta_offset_div2 = (int8_t)encoder->cfg->deblock_beta;
  encoder->tc_offset_div2   = (int8_t)encoder->cfg->deblock_tc;
  // SAO
  encoder->sao_enable = (int8_t)encoder->cfg->sao_enable;
  // RDO
  encoder->rdoq_enable = (int8_t)encoder->cfg->rdoq_enable;
  encoder->rdo         = (int8_t)encoder->cfg->rdo;
  // TR SKIP
  encoder->trskip_enable = (int8_t)encoder->cfg->trskip_enable;
  // VUI
  encoder->vui.sar_width   = (int16_t)encoder->cfg->vui.sar_width;
  encoder->vui.sar_height  = (int16_t)encoder->cfg->vui.sar_height;
  encoder->vui.overscan    = encoder->cfg->vui.overscan;
  encoder->vui.videoformat = encoder->cfg->vui.videoformat;
  encoder->vui.fullrange   = encoder->cfg->vui.fullrange;
  encoder->vui.colorprim   = encoder->cfg->vui.colorprim;
  encoder->vui.transfer    = encoder->cfg->vui.transfer;
  encoder->vui.colormatrix = encoder->cfg->vui.colormatrix;
  encoder->vui.chroma_loc  = (int8_t)encoder->cfg->vui.chroma_loc;
  // AUD
  encoder->aud_enable = (int8_t)encoder->cfg->aud_enable;
  // CQM
  cqmfile = cfg->cqmfile ? fopen(cfg->cqmfile, "rb") : NULL;
  encoder->cqmfile = cqmfile;

  init_encoder_input(&encoder->in, input, cfg->width, cfg->height);

  fprintf(stderr, "Input: %s, output: %s\n", cfg->input, cfg->output);
  fprintf(stderr, "  Video size: %dx%d (input=%dx%d)\n",
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

    // Skip '--seek' frames before input.
    // This block can be moved outside this while loop when there is a
    // mechanism to skip the while loop on error.
    if (encoder->frame == 0 && cfg->seek > 0) {
      int frame_bytes = cfg->width * cfg->height * 3 / 2;
      int error = 0;

      if (!strcmp(cfg->input, "-")) {
        // Input is stdin.
        int i;
        for (i = 0; !error && i < cfg->seek; ++i) {
          error = !read_one_frame(input, encoder);
        }
      } else {
        // input is a file. We hope. Proper detection is OS dependent.
        error = fseek(input, cfg->seek * frame_bytes, SEEK_CUR);
      }
      if (error && !feof(input)) {
        fprintf(stderr, "Failed to seek %d frames.\n", cfg->seek);
        break;
      }
    }

    // Read one frame from the input
    if (!read_one_frame(input, encoder)) {
      if (!feof(input))
        fprintf(stderr, "Failed to read a frame %d\n", encoder->frame);
      break;
    }

    // The actual coding happens here, after this function we have a coded frame
    encode_one_frame(encoder);

    if (cfg->debug != NULL) {
      // Write reconstructed frame out.
      // Use conformance-window dimensions instead of internal ones.
      const int width = encoder->in.width;
      const int out_width = encoder->in.real_width;
      const int out_height = encoder->in.real_height;
      int y;
      const pixel *y_rec = encoder->in.cur_pic->y_recdata;
      const pixel *u_rec = encoder->in.cur_pic->u_recdata;
      const pixel *v_rec = encoder->in.cur_pic->v_recdata;

      for (y = 0; y < out_height; ++y) {
        fwrite(&y_rec[y * width], sizeof(*y_rec), out_width, recout);
      }
      for (y = 0; y < out_height / 2; ++y) {
        fwrite(&u_rec[y * width / 2], sizeof(*u_rec), out_width / 2, recout);
      }
      for (y = 0; y < out_height / 2; ++y) {
        fwrite(&v_rec[y * width / 2], sizeof(*v_rec), out_width / 2, recout);
      }
    }

    // Calculate the bytes pushed to output for this frame
    fgetpos(output,(fpos_t*)&curpos);
    diff = (int32_t)(curpos-lastpos);
    lastpos = curpos;

    // PSNR calculations
    temp_psnr[0] = image_psnr(encoder->in.cur_pic->y_data, encoder->in.cur_pic->y_recdata, cfg->width, cfg->height);
    temp_psnr[1] = image_psnr(encoder->in.cur_pic->u_data, encoder->in.cur_pic->u_recdata, cfg->width>>1, cfg->height>>1);
    temp_psnr[2] = image_psnr(encoder->in.cur_pic->v_data, encoder->in.cur_pic->v_recdata, cfg->width>>1, cfg->height>>1);

    fprintf(stderr, "POC %4d (%c-frame) %10d bits PSNR: %2.4f %2.4f %2.4f\n", encoder->frame,
           "BPI"[encoder->in.cur_pic->slicetype%3], diff<<3,
           temp_psnr[0], temp_psnr[1], temp_psnr[2]);

    // Increment total PSNR
    psnr[0] += temp_psnr[0];
    psnr[1] += temp_psnr[1];
    psnr[2] += temp_psnr[2];


    // TODO: add more than one reference

    // Remove the ref pic (if present)
    if (encoder->ref->used_size == (uint32_t)encoder->cfg->ref_frames) {
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
  fprintf(stderr, " Processed %d frames, %10llu bits AVG PSNR: %2.4f %2.4f %2.4f\n", encoder->frame, (long long unsigned int) curpos<<3,
         psnr[0] / encoder->frame, psnr[1] / encoder->frame, psnr[2] / encoder->frame);
  fprintf(stderr, " Total time: %.3f s.\n", ((float)(clock() - start_time)) / CLOCKS_PER_SEC);

  fclose(input);
  fclose(output);
  if(cqmfile != NULL) fclose(cqmfile);
  if(recout != NULL) fclose(recout);

  // Deallocating
  config_destroy(cfg);
  scalinglist_destroy();
  picture_list_destroy(encoder->ref);
  picture_destroy(encoder->in.cur_pic);
  FREE_POINTER(encoder->in.cur_pic);
  free_bitstream(encoder->stream);
  free(encoder);
  free_tables();
  free_exp_golomb();

  return EXIT_SUCCESS;
}
