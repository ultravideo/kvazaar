/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2015 Tampere University of Technology and others (see
 * COPYING file).
 *
 * Kvazaar is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * Kvazaar is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
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

#include "kvazaar_internal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "checkpoint.h"
#include "global.h"
#include "config.h"
#include "threadqueue.h"
#include "encoder.h"
#include "encoderstate.h"
#include "cli.h"
#include "yuv_io.h"

/**
 * \brief Open a file for reading.
 *
 * If the file is "-", stdin is used.
 *
 * \param filename  name of the file to open or "-"
 * \return          the opened file or NULL if opening fails
 */
static FILE* open_input_file(const char* filename)
{
  if (!strcmp(filename, "-")) return stdin;
  return fopen(filename, "rb");
}

/**
 * \brief Open a file for writing.
 *
 * If the file is "-", stdout is used.
 *
 * \param filename  name of the file to open or "-"
 * \return          the opened file or NULL if opening fails
 */
static FILE* open_output_file(const char* filename)
{
  if (!strcmp(filename, "-")) return stdout;
  return fopen(filename, "wb");
}

/**
 * \brief Program main function.
 * \param argc Argument count from commandline
 * \param argv Argument list
 * \return Program exit state
 */
int main(int argc, char *argv[])
{
  int retval = EXIT_SUCCESS;

  cmdline_opts_t *opts = NULL; //!< Command line options
  kvz_encoder* enc = NULL;
  FILE *input  = NULL; //!< input file (YUV)
  FILE *output = NULL; //!< output file (HEVC NAL stream)
  FILE *recout = NULL; //!< reconstructed YUV output, --debug
  clock_t start_time = clock();
  clock_t encoding_start_cpu_time;
  CLOCK_T encoding_start_real_time;
  
  clock_t encoding_end_cpu_time;
  CLOCK_T encoding_end_real_time;

  // Stdin and stdout need to be binary for input and output to work.
  // Stderr needs to be text mode to convert \n to \r\n in Windows.
  #ifdef _WIN32
      _setmode( _fileno( stdin ),  _O_BINARY );
      _setmode( _fileno( stdout ), _O_BINARY );
      _setmode( _fileno( stderr ), _O_TEXT );
  #endif
      
  CHECKPOINTS_INIT();

  const kvz_api * const api = kvz_api_get(8);

  opts = cmdline_opts_parse(api, argc, argv);
  // If problem with command line options, print banner and shutdown.
  if (!opts) {
    print_version();
    print_help();

    goto exit_failure;
  }

  input = open_input_file(opts->input);
  if (input == NULL) {
    fprintf(stderr, "Could not open input file, shutting down!\n");
    goto exit_failure;
  }

  output = open_output_file(opts->output);
  if (output == NULL) {
    fprintf(stderr, "Could not open output file, shutting down!\n");
    goto exit_failure;
  }

  if (opts->debug != NULL) {
    recout = open_output_file(opts->debug);
    if (recout == NULL) {
      fprintf(stderr, "Could not open reconstruction file (%s), shutting down!\n", opts->debug);
      goto exit_failure;
    }
  }

  enc = api->encoder_open(opts->config);
  if (!enc) {
    fprintf(stderr, "Failed to open encoder.\n");
    goto exit_failure;
  }

  encoder_control_t *encoder = enc->control;
  
  fprintf(stderr, "Input: %s, output: %s\n", opts->input, opts->output);
  fprintf(stderr, "  Video size: %dx%d (input=%dx%d)\n",
         encoder->in.width, encoder->in.height,
         encoder->in.real_width, encoder->in.real_height);

  if (opts->seek > 0 && !yuv_io_seek(input, opts->seek, opts->config->width, opts->config->height)) {
    fprintf(stderr, "Failed to seek %d frames.\n", opts->seek);
    goto exit_failure;
  }
  encoder->vui.field_seq_flag = encoder->cfg->source_scan_type != 0;
  encoder->vui.frame_field_info_present_flag = encoder->cfg->source_scan_type != 0;

  //Now, do the real stuff
  {

    GET_TIME(&encoding_start_real_time);
    encoding_start_cpu_time = clock();

    uint64_t bitstream_length = 0;
    uint32_t frames_read = 0;
    uint32_t frames_done = 0;
    double psnr_sum[3] = { 0.0, 0.0, 0.0 };

    int8_t field_parity = 0;
    kvz_picture *frame_in = NULL;

    for (;;) {

      kvz_picture *img_in = NULL;

      if (!feof(input) && (opts->frames == 0 || frames_read < opts->frames || field_parity == 1) ) {
        // Try to read an input frame.
        if(field_parity == 0) frame_in = api->picture_alloc(opts->config->width, opts->config->height);
        if (!frame_in) {
          fprintf(stderr, "Failed to allocate image.\n");
          goto exit_failure;
        }

        if (field_parity == 0){
          if (yuv_io_read(input, opts->config->width, opts->config->height, frame_in)) {
            frames_read += 1;
            img_in = frame_in;
          } else {
            // EOF or some error
            api->picture_free(img_in);
            img_in = NULL;
            if (!feof(input)) {
              fprintf(stderr, "Failed to read a frame %d\n", frames_read);
              goto exit_failure;
            }
          }
        }

        if (encoder->cfg->source_scan_type != 0){
          img_in = api->picture_alloc(encoder->in.width, encoder->in.height);
          yuv_io_extract_field(frame_in, encoder->cfg->source_scan_type, field_parity, img_in);
          if (field_parity == 1) api->picture_free(frame_in);
          field_parity ^= 1; //0->1 or 1->0
        }
      }

      kvz_data_chunk* chunks_out = NULL;
      kvz_picture *img_rec = NULL;
      uint32_t len_out = 0;
      if (!api->encoder_encode(enc, img_in, &chunks_out, &len_out, &img_rec)) {
        fprintf(stderr, "Failed to encode image.\n");
        api->picture_free(img_in);
        goto exit_failure;
      }

      if (chunks_out == NULL && img_in == NULL) {
        // We are done since there is no more input and output left.
        break;
      }

      if (chunks_out != NULL) {
        uint64_t written = 0;
        // Write data into the output file.
        for (kvz_data_chunk *chunk = chunks_out;
             chunk != NULL;
             chunk = chunk->next) {
          assert(written + chunk->len <= len_out);
          if (fwrite(chunk->data, sizeof(uint8_t), chunk->len, output) != chunk->len) {
            fprintf(stderr, "Failed to write data to file.\n");
            api->picture_free(img_in);
            api->chunk_free(chunks_out);
            goto exit_failure;
          }
          written += chunk->len;
        }
        fflush(output);

        bitstream_length += len_out;

        // Compute and print stats.

        // Number of the state that was finished is one less than
        // enc->out_state_num.
        encoder_state_t *state = &enc->states[
          (enc->out_state_num + enc->num_encoder_states - 1) %
          enc->num_encoder_states
        ];
        double frame_psnr[3] = { 0.0, 0.0, 0.0 };
        encoder_compute_stats(state, frame_psnr);

        if (recout) {
          // Since chunks_out was not NULL, img_rec should have been set.
          assert(img_rec);
          if (!yuv_io_write(recout,
                            img_rec,
                            opts->config->width,
                            opts->config->height)) {
            fprintf(stderr, "Failed to write reconstructed picture!\n");
          }
        }

        frames_done += 1;
        psnr_sum[0] += frame_psnr[0];
        psnr_sum[1] += frame_psnr[1];
        psnr_sum[2] += frame_psnr[2];

        print_frame_info(state, frame_psnr);
      }

      api->picture_free(img_in);
      api->chunk_free(chunks_out);
      api->picture_free(img_rec);
    }

    GET_TIME(&encoding_end_real_time);
    encoding_end_cpu_time = clock();

    threadqueue_flush(encoder->threadqueue);
    // Coding finished

    // Print statistics of the coding
    fprintf(stderr, " Processed %d frames, %10llu bits",
            frames_done,
            (long long unsigned int)bitstream_length * 8);
    if (frames_done > 0) {
      fprintf(stderr, " AVG PSNR: %2.4f %2.4f %2.4f",
              psnr_sum[0] / frames_done,
              psnr_sum[1] / frames_done,
              psnr_sum[2] / frames_done);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, " Total CPU time: %.3f s.\n", ((float)(clock() - start_time)) / CLOCKS_PER_SEC);

    {
      double encoding_time = ( (double)(encoding_end_cpu_time - encoding_start_cpu_time) ) / (double) CLOCKS_PER_SEC;
      double wall_time = CLOCK_T_AS_DOUBLE(encoding_end_real_time) - CLOCK_T_AS_DOUBLE(encoding_start_real_time);
      fprintf(stderr, " Encoding time: %.3f s.\n", encoding_time);
      fprintf(stderr, " Encoding wall time: %.3f s.\n", wall_time);
      fprintf(stderr, " Encoding CPU usage: %.2f%%\n", encoding_time/wall_time*100.f);
      fprintf(stderr, " FPS: %.2f\n", ((double)frames_done)/wall_time);
    }
  }

  goto done;

exit_failure:
  retval = EXIT_FAILURE;

done:
  // deallocate structures
  if (enc) api->encoder_close(enc);
  if (opts) cmdline_opts_free(api, opts);

  // close files
  if (input)  fclose(input);
  if (output) fclose(output);
  if (recout) fclose(recout);

  CHECKPOINTS_FINALIZE();

  return retval;
}
