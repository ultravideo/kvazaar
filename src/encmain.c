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
#include <fcntl.h>    /* _O_BINARY */
#include <io.h>       /* _setmode() */
#endif

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // IWYU pragma: keep for CLOCKS_PER_SEC

#include <errno.h>

#include "checkpoint.h"
#include "cli.h"
#include "encoder.h"
#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"
#include "kvazaar_internal.h"
#include "threads.h"
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

static unsigned get_padding(unsigned width_or_height){
  if (width_or_height % CU_MIN_SIZE_PIXELS){
    return CU_MIN_SIZE_PIXELS - (width_or_height % CU_MIN_SIZE_PIXELS);
  }else{
    return 0;
  }
}

/**
 * \brief Value that is printed instead of PSNR when SSE is zero.
 */
static const double MAX_PSNR = 999.99;
static const double MAX_SQUARED_ERROR = (double)PIXEL_MAX * (double)PIXEL_MAX;

/**
 * \brief Calculates image PSNR value
 *
 * \param src   source picture
 * \param rec   reconstructed picture
 * \prama psnr  returns the PSNR
 */
static void compute_psnr(const kvz_picture *const src,
                         const kvz_picture *const rec,
                         double psnr[3])
{
  assert(src->width  == rec->width);
  assert(src->height == rec->height);

  int32_t pixels = src->width * src->height;
  int colors = rec->chroma_format == KVZ_CSP_400 ? 1 : 3;
  double sse[3] = { 0.0 };

  for (int32_t c = 0; c < colors; ++c) {
    int32_t num_pixels = pixels;
    if (c != COLOR_Y) {
      num_pixels >>= 2;
    }
    for (int32_t i = 0; i < num_pixels; ++i) {
      const int32_t error = src->data[c][i] - rec->data[c][i];
      sse[c] += error * error;
    }

    // Avoid division by zero
    if (sse[c] == 0.0) {
      psnr[c] = MAX_PSNR;
    } else {
      psnr[c] = 10.0 * log10(num_pixels * MAX_SQUARED_ERROR / sse[c]);
    }
  }
}

typedef struct {
  // Semaphores for synchronization.
  kvz_sem_t* available_input_slots;
  kvz_sem_t* filled_input_slots;

  // Parameters passed from main thread to input thread.
  FILE* input;
  const kvz_api *api;
  const cmdline_opts_t *opts;
  const encoder_control_t *encoder;
  uint8_t padding_x;
  uint8_t padding_y;

  // Picture and thread status passed from input thread to main thread.
  kvz_picture *img_in;
  int retval;

   // ***********************************************
   // Modified for SHVC
   const int input_layer;
   // ***********************************************

} input_handler_args;

#define RETVAL_RUNNING 0
#define RETVAL_FAILURE 1
#define RETVAL_EOF 2

/**
* \brief Handles input reading in a thread
*
* \param in_args  pointer to argument struct
*/
static void* input_read_thread(void* in_args)
{

  // Reading a frame works as follows:
  // - read full frame
  // if progressive: set read frame as output
  // if interlaced:
  // - allocate two fields and fill them according to field order
  // - deallocate the initial full frame

  input_handler_args* args = (input_handler_args*)in_args;
  kvz_picture *frame_in = NULL;
  int retval = RETVAL_RUNNING;
  int frames_read = 0;

  for (;;) {
    // Each iteration of this loop puts either a single frame or a field into
    // args->img_in for main thread to process.

    bool input_empty = !(args->opts->frames == 0 // number of frames to read is unknown
                         || frames_read < args->opts->frames); // not all frames have been read
    if (feof(args->input) || input_empty) {
      retval = RETVAL_EOF;
      goto done;
    }

    // ***********************************************
  // Modified for SHVC
    enum kvz_chroma_format csp = KVZ_FORMAT2CSP(args->opts->config->input_format);
    frame_in = args->api->picture_alloc_csp(csp,
                                            args->opts->config->shared->input_widths[args->input_layer]  + args->padding_x,
                                            args->opts->config->shared->input_heights[args->input_layer] + args->padding_y);

    if (!frame_in) {
      fprintf(stderr, "Failed to allocate image.\n");
      retval = RETVAL_FAILURE;
      goto done;
    }

    // Set PTS to make sure we pass it on correctly.
    frame_in->pts = frames_read;

    bool read_success = yuv_io_read(args->input,
                                    args->opts->config->shared->input_widths[args->input_layer],
                                    args->opts->config->shared->input_heights[args->input_layer],
                                    args->encoder->cfg.input_bitdepth,
                                    args->encoder->bitdepth,
                                    frame_in);
    if (!read_success) {
      // reading failed
      if (feof(args->input)) {
        // When looping input, re-open the file and re-read data.
        if (args->opts->loop_input && args->input != stdin) {
          fclose(args->input);
          args->input = fopen(args->opts->input[args->input_layer], "rb");
          if (args->input == NULL) {
            fprintf(stderr, "Could not re-open input file, shutting down!\n");
            retval = RETVAL_FAILURE;
            goto done;
          }
          bool read_success = yuv_io_read(args->input,
                                          args->opts->config->shared->input_widths[args->input_layer],
                                          args->opts->config->shared->input_heights[args->input_layer],
                                          args->encoder->cfg.input_bitdepth,
                                          args->encoder->bitdepth,
                                          frame_in);
          if (!read_success) {
            fprintf(stderr, "Could not re-open input file, shutting down!\n");
            retval = RETVAL_FAILURE;
            goto done;
          }
        } else {
          retval = RETVAL_EOF;
          goto done;
        }
      } else {
        fprintf(stderr, "Failed to read a frame %d\n", frames_read);
        retval = RETVAL_FAILURE;
        goto done;
      }
    }
    // **********************************************

    frames_read++;

    if (args->encoder->cfg.source_scan_type != 0) {
      // Set source scan type for frame, so that it will be turned into fields.
      frame_in->interlacing = args->encoder->cfg.source_scan_type;
    }

    // Wait until main thread is ready to receive the next frame.
    kvz_sem_wait(args->available_input_slots);
    args->img_in = frame_in;
    args->retval = retval;
    // Unlock main_thread_mutex to notify main thread that the new img_in
    // and retval have been placed to args.
    kvz_sem_post(args->filled_input_slots);

    frame_in = NULL;
  }

done:
  // Wait until main thread is ready to receive the next frame.
  kvz_sem_wait(args->available_input_slots);
  args->img_in = NULL;
  args->retval = retval;
  // Unlock main_thread_mutex to notify main thread that the new img_in
  // and retval have been placed to args.
  kvz_sem_post(args->filled_input_slots);

  // Do some cleaning up.
  args->api->picture_free(frame_in);

  pthread_exit(NULL);
  return NULL;
}
  
// ***********************************************
// Modified for SHVC
//Helper function for deallocating chained kvz_pictures
static void free_chained_img(const kvz_api *const api, kvz_picture *img )
{
  if( img == NULL ) return;
  while( img != img->base_image ) {
    kvz_picture *next_img = img->base_image;
    img->base_image = img;
    api->picture_free(img);
    img = next_img;
  }
  img->base_image = img;
  api->picture_free(img);
}
// ***********************************************

void output_recon_pictures(const kvz_api *const api,
                           FILE *recout,
                           kvz_picture *buffer[KVZ_MAX_GOP_LENGTH],
                           int *buffer_size,
                           uint64_t *next_pts,
                           unsigned width,
                           unsigned height)
{
  bool picture_written;
  do {
    picture_written = false;
    for (int i = 0; i < *buffer_size; i++) {

      kvz_picture *pic = buffer[i];
      if (pic->pts == *next_pts) {
        // Output the picture and remove it.
        if (!yuv_io_write(recout, pic, width, height)) {
          fprintf(stderr, "Failed to write reconstructed picture!\n");
        }
        api->picture_free(pic);
        picture_written = true;
        (*next_pts)++;

        // Move rest of the pictures one position backward.
        for (i++; i < *buffer_size; i++) {
          buffer[i - 1] = buffer[i];
          buffer[i] = NULL;
        }
        (*buffer_size)--;
      }
    }
  } while (picture_written);
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
  FILE **input  = NULL; //!< input files (YUV)
  FILE *output = NULL; //!< output file (HEVC NAL stream)
  FILE **recout = NULL; //!< reconstructed YUV outputs, --debug
  clock_t start_time = clock();
  clock_t encoding_start_cpu_time;
  KVZ_CLOCK_T encoding_start_real_time;

  clock_t encoding_end_cpu_time;
  KVZ_CLOCK_T encoding_end_real_time;

  // ***********************************************
  // Modified for SHVC
  
  // PTS of the reconstructed picture that should be output next.
  // Only used with --debug.
  uint64_t *next_recon_pts = NULL;
  // Buffer for storing reconstructed pictures that are not to be output
  // yet (i.e. in wrong order because GOP is used).
  // Only used with --debug.
  kvz_picture *(*recon_buffer)[KVZ_MAX_GOP_LENGTH] = NULL; //Feels so wrong but should mean recon_buffer is a pointer to kvz_picture* [KVZ_MAX_GOP_LENGTH] -arrays
  int *recon_buffer_size = NULL;

  //Need to declare these here because goto may jump over the initialization
  kvz_frame_info* info_out = NULL;
  uint32_t *len_out = NULL;
  pthread_t *input_threads = NULL;
 
  input_handler_args *in_args = NULL;

  uint64_t *substream_lengths = NULL; //Store total layer bitstream lengths
  double (*layer_psnr_sum)[3] = NULL; //Per layer psnr sums
  //*************************************************

  // Semaphores for synchronizing the input reader thread and the main
  // thread.
  //
  // available_input_slots tells whether the main thread is currently using
  // input_handler_args.img_in. (0 = in use, 1 = not in use)
  //
  // filled_input_slots tells whether there is a new input picture (or NULL
  // if the input has ended) in input_handler_args.img_in placed by the
  // input reader thread. (0 = no new image, 1 = one new image)
  //
  kvz_sem_t **available_input_slots = NULL;
  kvz_sem_t **filled_input_slots = NULL;

#ifdef _WIN32
  // Stderr needs to be text mode to convert \n to \r\n in Windows.
  setmode( _fileno( stderr ), _O_TEXT );
#endif

  CHECKPOINTS_INIT();

  const kvz_api * const api = kvz_api_get(8);

  opts = cmdline_opts_parse(api, argc, argv);
  // If problem with command line options, print banner and shutdown.
  if (!opts) {
    print_usage();

    goto exit_failure;
  }
  if (opts->version) {
    print_version();
    goto done;
  }
  if (opts->help) {
    print_help();
    goto done;
  }

  // ***********************************************
  // Modified for SHVC
  input = calloc(opts->num_inputs, sizeof(FILE*));
  for (int8_t i = 0; i < opts->num_inputs; i++) {
    input[i] = open_input_file(opts->input[i]);
    if (input[i] == NULL) {
      fprintf(stderr, "Could not open input file, shutting down!\n");
      goto exit_failure;
    }
  }

  output = open_output_file(opts->output);
  if (output == NULL) {
    fprintf(stderr, "Could not open output file, shutting down!\n");
    goto exit_failure;
  }

#ifdef _WIN32
  // Set stdin and stdout to binary for pipes.
  if (*input == stdin) {
    _setmode(_fileno(stdin), _O_BINARY);
  }
  if (output == stdout) {
    _setmode(_fileno(stdout), _O_BINARY);
  }
#endif

  recout = calloc(opts->num_debugs, sizeof(FILE*));
  for (int8_t i = 0; i < opts->num_debugs; i++) {
    if (opts->debug[i] != NULL) {
      recout[i] = open_output_file(opts->debug[i]);
      if (recout[i] == NULL) {
        fprintf(stderr, "Could not open reconstruction file (%s), shutting down!\n", opts->debug[i]);
        goto exit_failure;
      }
    }
  }
  
  enc = api->encoder_open(opts->config);
  if (!enc) {
    fprintf(stderr, "Failed to open encoder.\n");
    goto exit_failure;
  }

  const encoder_control_t *encoder = enc->control;

  for ( int i = 0; i < opts->num_inputs; i++) {
    fprintf(stderr, "Input layer %d:\n", i);
    fprintf(stderr, "  Input: %s, output: %s\n", opts->input[i], opts->output);
    fprintf(stderr, "    Video input size: %dx%d\n",
      opts->config->shared->input_widths[i], opts->config->shared->input_heights[i]);

    if (opts->seek > 0 && !yuv_io_seek(input[i], opts->seek, opts->config->shared->input_widths[i], opts->config->shared->input_heights[i])) {
      fprintf(stderr, "Failed to seek %d frames.\n", opts->seek);
      goto exit_failure;
    }
  }
  for( const encoder_control_t* ctrl = encoder; ctrl != NULL ; ctrl = ctrl->next_enc_ctrl) {
    fprintf(stderr, "Layer %d:\n", ctrl->layer.layer_id);
    fprintf(stderr, "  Input layer: %d\n", ctrl->layer.input_layer);
    fprintf(stderr, "    Video size: %dx%d (input=%dx%d)\n",
      ctrl->in.width, ctrl->in.height,
      ctrl->in.real_width, ctrl->in.real_height);
    if( ctrl->layer.layer_id < opts->num_debugs ){
      fprintf(stderr, "  Debug output: %s\n", opts->debug[ctrl->layer.layer_id]);
    }
  }
//******************************************
  
    
  // ***********************************************
  // Modified for SHVC
  //Allocate space for some stuff
  info_out = malloc(sizeof(kvz_frame_info) * opts->config->shared->max_layers);
  len_out = calloc( opts->config->shared->max_layers, sizeof(uint32_t)); // Each layer has their own len_out

  input_threads = malloc(sizeof(pthread_t)*opts->num_inputs);

  //Allocate semaphores
  available_input_slots = calloc(opts->num_inputs, sizeof(kvz_sem_t*));
  filled_input_slots = calloc(opts->num_inputs, sizeof(kvz_sem_t*));

  in_args = calloc(opts->num_inputs,sizeof(input_handler_args));

  //Allocate stuff for the debug output buffers etc.
  next_recon_pts = calloc(opts->num_debugs,sizeof(uint64_t));
  recon_buffer = calloc(opts->num_debugs,sizeof(kvz_picture*[KVZ_MAX_GOP_LENGTH])); //Feels so wrong but should mean recon_buffer is a pointer to kvz_picture* [KVZ_MAX_GOP_LENGTH] -arrays
  recon_buffer_size = calloc(opts->num_debugs,sizeof(int));

  substream_lengths = calloc( opts->config->shared->max_layers,sizeof(uint64_t));
  layer_psnr_sum = calloc( opts->config->shared->max_layers,sizeof(double[3]));

  //Now, do the real stuff
  {

    KVZ_GET_TIME(&encoding_start_real_time);
    encoding_start_cpu_time = clock();

    uint64_t bitstream_length = 0;
    uint32_t frames_done = 0;
    double psnr_sum[3] = { 0.0, 0.0, 0.0 };

    for (int i = 0; i < opts->num_inputs; i++) {

      uint8_t padding_x = get_padding( opts->config->shared->input_widths[i] );
      uint8_t padding_y = get_padding( opts->config->shared->input_heights[i] );

      available_input_slots[i] = calloc(1, sizeof(kvz_sem_t));
      filled_input_slots[i]    = calloc(1, sizeof(kvz_sem_t));
      kvz_sem_init(available_input_slots[i], 0);
      kvz_sem_init(filled_input_slots[i],    0);

      in_args[i].available_input_slots = available_input_slots[i];
      in_args[i].filled_input_slots = filled_input_slots[i];

      in_args[i].input = input[i];
      in_args[i].api = api;
      in_args[i].opts = opts;
      in_args[i].encoder = encoder;
      in_args[i].padding_x = padding_x;
      in_args[i].padding_y = padding_y;
      
      in_args[i].img_in = NULL;
      in_args[i].retval = RETVAL_RUNNING;
      
      
      if (pthread_create(&input_threads[i], NULL, input_read_thread, (void*)&in_args[i]) != 0) {
        fprintf(stderr, "pthread_create failed!\n");
        assert(0);
        return 0;
      }
    }
    // ***********************************************
    kvz_picture *cur_in_img = NULL;
    for (;;) {

      // ***********************************************
      // Modified for SHVC
      //TODO: Move relevant de/allocation to done tag.
      kvz_picture *cur_img = NULL;
      for (int i = 0; i < opts->num_inputs; i++) {
        if (in_args[i].retval == RETVAL_RUNNING) {
          // Increase available_input_slots so that the input thread can
          // write the new img_in and retval to in_args.
          kvz_sem_post(available_input_slots[i]);
          // Wait until the input thread has updated in_args and then
          // decrease filled_input_slots.
          kvz_sem_wait(filled_input_slots[i]);

          if( cur_in_img == NULL ){
            cur_in_img = in_args[i].img_in;
            cur_img = cur_in_img;
          } else {
            cur_img->base_image = in_args[i].img_in;
            cur_img = cur_img->base_image;
          }
          in_args[i].img_in = NULL;
        } else if(in_args[i].retval == RETVAL_EOF) {
          cur_in_img = NULL;
          break;
        } else if (in_args[i].retval == EXIT_FAILURE) {
          goto exit_failure;
        }
      }

      kvz_data_chunk* chunks_out = NULL;
      kvz_picture *img_rec = NULL;
      kvz_picture *img_src = NULL;
      
      if (!api->encoder_encode(enc,
                               cur_in_img,
                               &chunks_out,
                               len_out,
                               &img_rec,
                               &img_src,
                               info_out)) {
        fprintf(stderr, "Failed to encode image.\n");
        free_chained_img(api, cur_in_img);
        goto exit_failure;
      }
      // ***********************************************

      if (chunks_out == NULL && cur_in_img == NULL) {
        // We are done since there is no more input and output left.
        break;
      }

      // ***********************************************
      // Modified for SHVC
      uint32_t tot_len_out = 0;
      for (int i = 0; i < opts->config->shared->max_layers; ++i) {
        tot_len_out += len_out[i];
        substream_lengths[i] += len_out[i];
      }

      if (chunks_out != NULL) {
        uint64_t written = 0;
        // Write data into the output file.
        for (kvz_data_chunk *chunk = chunks_out;
             chunk != NULL;
             chunk = chunk->next) {
          assert(written + chunk->len <= tot_len_out);
          if (fwrite(chunk->data, sizeof(uint8_t), chunk->len, output) != chunk->len) {
            fprintf(stderr, "Failed to write data to file.\n");
            free_chained_img(api, cur_in_img);
            api->chunk_free(chunks_out);
            goto exit_failure;
          }
          written += chunk->len;
        }
        fflush(output);

        bitstream_length += tot_len_out;

        // ***********************************************
        // Modified for SHVC

        // Compute and print stats.

        double frame_psnr[3] = { 0.0, 0.0, 0.0 };
        kvz_config *cfg = opts->config;
        kvz_picture *cur_src = img_src;
        kvz_picture *cur_rec = img_rec; img_rec = NULL; //Set img_rec later to the first rec that is not output
        int layer_id = 0;

        for (layer_id = 0; cfg != NULL; layer_id++) { 
          if (cfg->calc_psnr && cfg->source_scan_type == KVZ_INTERLACING_NONE) {
            // Do not compute PSNR for interlaced frames, because img_rec does not contain
            // the deinterlaced frame yet.
            compute_psnr(cur_src, cur_rec, frame_psnr);
          }

          //Write recout only for as many layers as rec outs were given
          if (layer_id < opts->num_debugs) {
            // Since chunks_out was not NULL, img_rec should have been set.
            assert(cur_rec);

            //Check that recon file has been set
            if(!recout[layer_id]) {
              printf("Error: Recout file not set for layer %d.\n",layer_id);
              free_chained_img(api, cur_in_img);
              api->chunk_free(chunks_out);
              free_chained_img(api, cur_rec);
              free_chained_img(api, img_src);
              goto exit_failure;
            }
            
            // Move img_rec to the recon buffer.
            assert(recon_buffer_size[layer_id] < KVZ_MAX_GOP_LENGTH);
            recon_buffer[layer_id][recon_buffer_size[layer_id]++] = cur_rec;
            //cur_rec = NULL;
            //Cur rec might be freed in output_recon_pictures so we need to move to the next picture here
            kvz_picture *tmp = cur_rec;
            cur_rec = cur_rec->base_image;
            tmp->base_image = tmp;

            // Try to output some reconstructed pictures.
            output_recon_pictures(api,
                                  recout[layer_id],
                                  recon_buffer[layer_id],
                                  &recon_buffer_size[layer_id],
                                  &next_recon_pts[layer_id],
                                  cfg->width,
                                  cfg->height);
          }
          else {
            if( img_rec == NULL) {
              //Should be the first time in the loop the else branch has been taken and cur_rec is the first rec not output
              img_rec = cur_rec;
            }
            cur_rec = cur_rec->base_image;
          }

          frames_done += 1;
          psnr_sum[0] += frame_psnr[0];
          psnr_sum[1] += frame_psnr[1];
          psnr_sum[2] += frame_psnr[2];

          //Set per layer psnr_sum
          layer_psnr_sum[layer_id][0] += frame_psnr[0];
          layer_psnr_sum[layer_id][1] += frame_psnr[1];
          layer_psnr_sum[layer_id][2] += frame_psnr[2];

          print_frame_info(&(info_out[layer_id]), frame_psnr, len_out[layer_id]); 

          //Update stuff. The images of different layers should be chained together using base_image;
          cfg = cfg->next_cfg;
          if(cur_src) cur_src = cur_src->base_image;
        }
      }

      free_chained_img(api, cur_in_img); cur_in_img = NULL;
      api->chunk_free(chunks_out);
      free_chained_img(api, img_rec); //img_rec should contain image not output. Others should be freed by output_recon_pictures.
      img_rec = NULL;  
      free_chained_img(api, img_src); img_src = NULL;
    }
    
    KVZ_GET_TIME(&encoding_end_real_time);
    encoding_end_cpu_time = clock();
    // Coding finished

    // All reconstructed pictures should have been output.
    for (int i = 0; i < opts->num_debugs; i++) {
      assert(recon_buffer_size[i] == 0);
    }
    // Print statistics of the coding
    fprintf(stderr, " Processed %d frames over %d layer(s), %10llu bits",
            frames_done / opts->config->shared->max_layers,
            opts->config->shared->max_layers,
            (long long unsigned int)bitstream_length * 8);
    if (frames_done > 0) {
      // ***********************************************
      // Modified for SHVC.
      fprintf(stderr, " AVG PSNR: %2.4f %2.4f %2.4f",
              psnr_sum[0] / frames_done,
              psnr_sum[1] / frames_done,
              psnr_sum[2] / frames_done);
      // ***********************************************
    }
    fprintf(stderr, "\n");
    //Print layer stats if more than two layers
    if ( opts->config->shared->max_layers > 1) {
      for (int i = 0; i < opts->config->shared->max_layers; ++i) {
        fprintf(stderr, "  Layer %d: %10llu bits,",
          i,
          (long long unsigned int)substream_lengths[i] * 8);
        if (frames_done > 0) {
          fprintf(stderr, " AVG PSNR: %2.4f %2.4f %2.4f",
            layer_psnr_sum[i][0] / frames_done * opts->config->shared->max_layers,
            layer_psnr_sum[i][1] / frames_done * opts->config->shared->max_layers,
            layer_psnr_sum[i][2] / frames_done * opts->config->shared->max_layers );
        }
        fprintf(stderr, "\n");
      }
    }

    fprintf(stderr, " Total CPU time: %.3f s.\n", ((float)(clock() - start_time)) / CLOCKS_PER_SEC);

    {
      double encoding_time = ( (double)(encoding_end_cpu_time - encoding_start_cpu_time) ) / (double) CLOCKS_PER_SEC;
      double wall_time = KVZ_CLOCK_T_AS_DOUBLE(encoding_end_real_time) - KVZ_CLOCK_T_AS_DOUBLE(encoding_start_real_time);
      fprintf(stderr, " Encoding time: %.3f s.\n", encoding_time);
      fprintf(stderr, " Encoding wall time: %.3f s.\n", wall_time);
      fprintf(stderr, " Encoding CPU usage: %.2f%%\n", encoding_time/wall_time*100.f);
      fprintf(stderr, " FPS: %.2f\n", ((double)frames_done)/wall_time);
    }
    for (int i = 0; i < opts->num_inputs; i++) {
      pthread_join(input_threads[i], NULL);
    }

    //********************************
  }

  goto done;

exit_failure:
  retval = EXIT_FAILURE;

done:
  // ***********************************************
  // Modified for SHVC

  // deallocate structures
  if (enc) api->encoder_close(enc);
  
  // close files
  if (opts != NULL && input != NULL) {
    for (int8_t i = 0; i < opts->num_inputs; i++) {
      if (input[i])  fclose(input[i]);
    }
    for (int8_t i = 0; i < opts->num_debugs; i++) {
      if (recout[i]) fclose(recout[i]);
    }
  }
  FREE_POINTER(input);
  if (output) fclose(output);  
  FREE_POINTER(recout);

  //Free some stuff
  FREE_POINTER(info_out);
  FREE_POINTER(len_out);
  FREE_POINTER(in_args);
  FREE_POINTER(input_threads);

  // destroy semaphores
  if (opts != NULL && available_input_slots != NULL && filled_input_slots != NULL) {
    for (int i = 0; i < opts->num_inputs; i++) {
      if (available_input_slots[i]) kvz_sem_destroy(available_input_slots[i]);
      if (filled_input_slots[i])    kvz_sem_destroy(filled_input_slots[i]);
      FREE_POINTER(available_input_slots[i]);
      FREE_POINTER(filled_input_slots[i]);
    }
  }
  FREE_POINTER(available_input_slots);
  FREE_POINTER(filled_input_slots);

  FREE_POINTER(next_recon_pts);
  FREE_POINTER(recon_buffer);
  FREE_POINTER(recon_buffer_size);

  FREE_POINTER(layer_psnr_sum);
  FREE_POINTER(substream_lengths);
  
  if (opts) cmdline_opts_free(api, opts);
  // ***********************************************

  CHECKPOINTS_FINALIZE();

  return retval;
}
