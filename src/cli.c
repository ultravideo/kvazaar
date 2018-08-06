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

#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>

static const char short_options[] = "i:o:d:w:h:n:q:p:r:";
static const struct option long_options[] = {
  { "input",              required_argument, NULL, 'i' },
  { "output",             required_argument, NULL, 'o' },
  { "debug",              required_argument, NULL, 'd' },
  { "width",              required_argument, NULL, 'w' }, // deprecated
  { "height",             required_argument, NULL, 'h' }, // deprecated
  { "frames",             required_argument, NULL, 'n' },
  { "qp",                 required_argument, NULL, 'q' },
  { "period",             required_argument, NULL, 'p' },
  { "ref",                required_argument, NULL, 'r' },
  { "vps-period",         required_argument, NULL, 0 },
  { "input-res",          required_argument, NULL, 0 },
  { "input-fps",          required_argument, NULL, 0 },
  { "deblock",            required_argument, NULL, 0 },
  { "no-deblock",               no_argument, NULL, 0 },
  { "sao",                optional_argument, NULL, 0 },
  { "no-sao",                   no_argument, NULL, 0 },
  { "rdoq",                     no_argument, NULL, 0 },
  { "no-rdoq",                  no_argument, NULL, 0 },
  { "signhide",                 no_argument, NULL, 0 },
  { "no-signhide",              no_argument, NULL, 0 },
  { "smp",                      no_argument, NULL, 0 },
  { "no-smp",                   no_argument, NULL, 0 },
  { "amp",                      no_argument, NULL, 0 },
  { "no-amp",                   no_argument, NULL, 0 },
  { "rd",                 required_argument, NULL, 0 },
  { "full-intra-search",        no_argument, NULL, 0 },
  { "no-full-intra-search",     no_argument, NULL, 0 },
  { "transform-skip",           no_argument, NULL, 0 },
  { "no-transform-skip",        no_argument, NULL, 0 },
  { "tr-depth-intra",     required_argument, NULL, 0 },
  { "me",                 required_argument, NULL, 0 },
  { "subme",              required_argument, NULL, 0 },
  { "source-scan-type",   required_argument, NULL, 0 },
  { "sar",                required_argument, NULL, 0 },
  { "overscan",           required_argument, NULL, 0 },
  { "videoformat",        required_argument, NULL, 0 },
  { "range",              required_argument, NULL, 0 },
  { "colorprim",          required_argument, NULL, 0 },
  { "transfer",           required_argument, NULL, 0 },
  { "colormatrix",        required_argument, NULL, 0 },
  { "chromaloc",          required_argument, NULL, 0 },
  { "aud",                      no_argument, NULL, 0 },
  { "no-aud",                   no_argument, NULL, 0 },
  { "cqmfile",            required_argument, NULL, 0 },
  { "seek",               required_argument, NULL, 0 },
  { "tiles",              required_argument, NULL, 0 },
  { "tiles-width-split",  required_argument, NULL, 0 },
  { "tiles-height-split", required_argument, NULL, 0 },
  { "wpp",                      no_argument, NULL, 0 },
  { "no-wpp",                   no_argument, NULL, 0 },
  { "owf",                required_argument, NULL, 0 },
  { "slices",             required_argument, NULL, 0 },
  { "threads",            required_argument, NULL, 0 },
  { "cpuid",              optional_argument, NULL, 0 },
  { "no-cpuid",                 no_argument, NULL, 0 },
  { "pu-depth-inter",     required_argument, NULL, 0 },
  { "pu-depth-intra",     required_argument, NULL, 0 },
  { "info",                     no_argument, NULL, 0 },
  { "no-info",                  no_argument, NULL, 0 },
  { "gop",                required_argument, NULL, 0 },
  { "bipred",                   no_argument, NULL, 0 },
  { "no-bipred",                no_argument, NULL, 0 },
  { "bitrate",            required_argument, NULL, 0 },
  { "preset",             required_argument, NULL, 0 },
  { "mv-rdo",                   no_argument, NULL, 0 },
  { "no-mv-rdo",                no_argument, NULL, 0 },
  { "psnr",                     no_argument, NULL, 0 },
  { "no-psnr",                  no_argument, NULL, 0 },
  { "version",                  no_argument, NULL, 0 },
  { "help",                     no_argument, NULL, 0 },
  { "loop-input",               no_argument, NULL, 0 },
  { "mv-constraint",      required_argument, NULL, 0 },
  { "hash",               required_argument, NULL, 0 },
  {"cu-split-termination",required_argument, NULL, 0 },
  { "crypto",             required_argument, NULL, 0 },
  { "key",                required_argument, NULL, 0 },
  { "me-early-termination",required_argument, NULL, 0 },
  { "intra-rdo-et",             no_argument, NULL, 0 },
  { "no-intra-rdo-et",          no_argument, NULL, 0 },
  { "lossless",                 no_argument, NULL, 0 },
  { "no-lossless",              no_argument, NULL, 0 },
  { "tmvp",                     no_argument, NULL, 0 },
  { "no-tmvp",                  no_argument, NULL, 0 },
  { "rdoq-skip",                no_argument, NULL, 0 },
  { "no-rdoq-skip",             no_argument, NULL, 0 },
  { "input-bitdepth",     required_argument, NULL, 0 },
  { "input-format",       required_argument, NULL, 0 },
  { "implicit-rdpcm",           no_argument, NULL, 0 },
  { "no-implicit-rdpcm",        no_argument, NULL, 0 },
  { "roi",                required_argument, NULL, 0 },
  { "erp-aqp",                  no_argument, NULL, 0 },
  { "no-erp-aqp",               no_argument, NULL, 0 },
  { "level",              required_argument, NULL, 0 },
  { "force-level",        required_argument, NULL, 0 },
  { "high-tier",                no_argument, NULL, 0 },
  { "me-steps",           required_argument, NULL, 0 },
  { "fast-residual-cost", required_argument, NULL, 0 },
  //*********************************************
  //For scalable extension.
  { "layer",                    no_argument, NULL, 0 }, //New layer
  { "layer-res",          required_argument, NULL, 0 }, //Set resolution of layer
  { "input-layer",        required_argument, NULL, 0 }, //Manualy set the input layer the current layer uses
  { "ilr",                required_argument, NULL, 0 }, //Number of ILR
  { "print-es-hierarchy",       no_argument, NULL, 0 }, //Print encoder state hierarchy
  //*********************************************
  {0, 0, 0, 0}
};

/**
* \brief Try to detect resolution from file name automatically
*
* \param file_name    file name to get dimensions from
* \param out_width    detected width
* \param out_height   detected height
* \return      1 if the resolution is set, 0 on fail
*/
static int select_input_res_auto(const char *file_name, int32_t *out_width, int32_t *out_height)
{
  if (!file_name) return 0;

  // Find the last delimiter char ( / or \ ). Hope the other kind is not used in the name.
  // If delim is not found, set pointer to the beginning.
  unsigned char* sub_str = (unsigned char*)MAX(strrchr(file_name, '/'), strrchr(file_name, '\\'));
  if (!sub_str) sub_str = (unsigned char*)file_name;

  int success = 0;
  // Try if the substring starts with "<int>x<int>" without either of them being 0
  do {
    success = (sscanf((char*)sub_str, "%dx%d%*s", out_width, out_height) == 2);
    success &= (*out_width > 0 && *out_height > 0);
    // Move to the next char until a digit is found or the string ends
    do{
      ++sub_str;
    } while (*sub_str != 0 && !isdigit(*sub_str));
    // Continue until "<int>x<int>" is found or the string ends
  } while (*sub_str != 0 && !success);

  return success;
}

/**
 * \brief Parse command line arguments.
 * \param argc  Number of arguments
 * \param argv  Argument list
 * \return      Pointer to the parsed options, or NULL on failure.
 */
cmdline_opts_t* cmdline_opts_parse(const kvz_api *const api, int argc, char *argv[])
{
  int ok = 1;
  cmdline_opts_t *opts = calloc(1, sizeof(cmdline_opts_t));
  if (!opts) {
    fprintf(stderr, "Memory error: Could not allocate or init opts.\n");
    ok = 0;
    goto done;
  }
  //*********************************************
  //For scalable extension. TODO: Add error checking
  //Initialize input and debug arrays
  opts->input = calloc(1, sizeof(char*));
  opts->num_inputs = 1;
  opts->debug = NULL;//calloc(1, sizeof(char*));
  opts->num_debugs = 0;
  //*********************************************

  opts->config = api->config_alloc();
  if (!opts->config || !api->config_init(opts->config)) {
    fprintf(stderr, "Memory error: Could not allocate or init config.\n");
    ok = 0;
    goto done;
  }

  // Parse command line options
  for (optind = 0;;) {
    int long_options_index = -1;

    int c = getopt_long(argc, argv, short_options, long_options, &long_options_index);
    if (c == -1)
      break;

    if (long_options_index < 0) {
      int i;
      for (i = 0; long_options[i].name; i++)
        if (long_options[i].val == c) {
            long_options_index = i;
            break;
        }
      if (long_options_index < 0) {
        // getopt_long already printed an error message
        ok = 0;
        goto done;
      }
    }

    //*********************************************
    //For scalable extension. TODO: Add error cheching
    const char* name = long_options[long_options_index].name;
    if (!strcmp(name, "input")) {
      if (*opts->input != NULL) {
       /* fprintf(stderr, "Input error: More than one input file given.\n");
        ok = 0;
        goto done;*/
        //Increase input array size
        char **tmp = realloc(opts->input, ++opts->num_inputs * sizeof(char*));
        if (tmp == NULL) {
          fprintf(stderr, "Memory error: Could not realloc input array.\n");
          ok = 0;
          goto done;
        }
        opts->input = tmp;

      }
      opts->input[opts->num_inputs-1] = strdup(optarg);
      api->config_parse(opts->config, name, optarg); //Send input for parsing so that input layers can be set correctly.
    } else if (!strcmp(name, "output")) {
      if (opts->output) {
        fprintf(stderr, "Input error: More than one output file given.\n");
        ok = 0;
        goto done;
      }
      opts->output = strdup(optarg);
    } else if (!strcmp(name, "debug")) {
      if (opts->debug != NULL) {
        /*fprintf(stderr, "Input error: More than one debug output file given.\n");
        ok = 0;
        goto done;*/
        char **tmp = realloc(opts->debug, ++opts->num_debugs * sizeof(char*));
        if (tmp == NULL) {
          fprintf(stderr, "Memory error: Could not realloc input array.\n");
          ok = 0;
          goto done;
        }
        opts->debug = tmp;
      }
      else {
        opts->debug = calloc(1, sizeof(char*));
        if (opts->debug == NULL) {
          fprintf(stderr, "Memory error: Could not alloc debug array.\n");
          ok = 0;
          goto done;
        }
        opts->num_debugs = 1;
      }
      opts->debug[opts->num_debugs-1] = strdup(optarg);
    } else if (!strcmp(name, "seek")) {
      opts->seek = atoi(optarg);
    } else if (!strcmp(name, "frames")) {
      opts->frames = atoi(optarg);
    } else if (!strcmp(name, "version")) {
      opts->version = true;
      goto done;
    } else if (!strcmp(name, "help")) {
      opts->help = true;
      goto done;
    } else if (!strcmp(name, "loop-input")) {
      opts->loop_input = true;
    } else if (!api->config_parse(opts->config, name, optarg)) {
      fprintf(stderr, "invalid argument: %s=%s\n", name, optarg);
      ok = 0;
      goto done;
    }
  }
  //*********************************************
  
  // Check for extra arguments.
  if (argc - optind > 0) {
    fprintf(stderr, "Input error: Extra argument found: \"%s\"\n", argv[optind]);
    ok = 0;
    goto done;
  }

  // Check that the required files were defined
  if (opts->input == NULL || opts->output == NULL) {
    fprintf(stderr, "Input error: Need to specify input and output files.\n");
    ok = 0;
    goto done;
  }

  if (opts->config->vps_period < 0) {
    // Disabling parameter sets is only possible when using Kvazaar as
    // a library.
    fprintf(stderr, "Input error: vps_period must be non-negative\n");
    ok = 0;
    goto done;
  }
  //*********************************************
  //For scalable extension. TODO: move to cfg parsing?
  // Set resolution automatically if necessary
  for (int i = 0; i < opts->num_inputs; i++) {
    //Automatically set layer size to match the respective input layer size
    kvz_config *cfg = opts->config;
    if (cfg->shared != NULL){
      if (cfg->shared->input_widths[i] == 0 && cfg->shared->input_heights[i] == 0) {
        if (!select_input_res_auto(opts->input[i], &cfg->shared->input_widths[i], &cfg->shared->input_heights[i])) {
          fprintf(stderr, "Input error: No size found for input layer %d\n", i);
          ok = 0;
        }
      }
    } else if (cfg->width == 0 && cfg->height == 0) {
      if(!select_input_res_auto(opts->input[i], &opts->config->width, &opts->config->height)){
        fprintf(stderr, "Input error: No size found for input\n");
        ok = 0;
      }
    }

    while( ok && cfg != NULL ) {
      if (cfg->shared != NULL) {
        //Default input_layer to highest input layer
        if (cfg->input_layer == -1) {
          cfg->input_layer = cfg->shared->max_input_layers - 1;
        }
        // If layer width/height not set use selected input layers size
        if (cfg->input_layer == i && cfg->width == 0 && cfg->height == 0) {
          cfg->width = cfg->shared->input_widths[i];
          cfg->height = cfg->shared->input_heights[i];
        }
      }
      cfg = cfg->next_cfg;
    }

  }
  //*********************************************

done:

  if (!ok) {
    cmdline_opts_free(api, opts);
    opts = NULL;
  }

  return opts;
}


/**
 * \brief Deallocate a cmdline_opts_t structure.
 */
void cmdline_opts_free(const kvz_api *const api, cmdline_opts_t *opts)
{
  if (opts) {
    //*********************************************
    //For scalable extension
    for (size_t i = 0; i < opts->num_inputs; i++) {
      FREE_POINTER(opts->input[i]);
    }
    FREE_POINTER(opts->input);
    FREE_POINTER(opts->output);
    for (size_t i = 0; i < opts->num_debugs; i++) {
      FREE_POINTER(opts->debug[i]);
    }
    FREE_POINTER(opts->debug);
    //*********************************************
    api->config_destroy(opts->config);
    opts->config = NULL;
  }
  FREE_POINTER(opts);
}


void print_usage(void)
{
  fprintf(stdout,
    "Kvazaar usage: -i and --input-res to set input, -o to set output\n"
    "               --help for more information\n");
}


void print_version(void)
{
  fprintf(stdout,
    "Kvazaar " VERSION_STRING "\n"
    "Kvazaar license: LGPL version 2\n");
}


void print_help(void)
{
  fprintf(stdout,
    "Usage:\n"
    "kvazaar -i <input> --input-res <width>x<height> -o <output>\n"
    "\n"
    /* Word wrap to this width to stay under 80 characters (including ") *************/
    "Required:\n"
    "  -i, --input <filename>     : Input file\n"
    "      --input-res <res>      : Input resolution [auto]\n"
    "                                   - auto: Detect from file name.\n"
    "                                   - <int>x<int>: width times height\n"
    "  -o, --output <filename>    : Output file\n"
    "\n"
    /* Word wrap to this width to stay under 80 characters (including ") *************/
    "Presets:\n"
    "      --preset <preset>      : Set options to a preset [medium]\n"
    "                                   - ultrafast, superfast, veryfast, faster,\n"
    "                                     fast, medium, slow, slower, veryslow\n"
    "                                     placebo\n"
    "\n"
    /* Word wrap to this width to stay under 80 characters (including ") *************/
    "Input:\n"
    "  -n, --frames <integer>     : Number of frames to code [all]\n"
    "      --seek <integer>       : First frame to code [0]\n"
    "      --input-fps <num>[/<denom>] : Frame rate of the input video [25]\n"
    "      --source-scan-type <string> : Source scan type [progressive]\n"
    "                                   - progressive: Progressive scan\n"
    "                                   - tff: Top field first\n"
    "                                   - bff: Bottom field first\n"
    "      --input-format <string> : P420 or P400 [P420]\n"
    "      --input-bitdepth <int> : 8-16 [8]\n"
    "      --loop-input           : Re-read input file forever.\n"
    "\n"
    /* Word wrap to this width to stay under 80 characters (including ") *************/
    "Options:\n"
    "      --help                 : Print this help message and exit.\n"
    "      --version              : Print version information and exit.\n"
    "      --(no-)aud             : Use access unit delimiters. [disabled]\n"
    "      --debug <filename>     : Output internal reconstruction.\n"
    "      --(no-)cpuid           : Enable runtime CPU optimizations. [enabled]\n"
    "      --hash <string>        : Decoded picture hash [checksum]\n"
    "                                   - none: 0 bytes\n"
    "                                   - checksum: 18 bytes\n"
    "                                   - md5: 56 bytes\n"
    "      --(no-)psnr            : Calculate PSNR for frames. [enabled]\n"
    "      --(no-)info            : Add encoder info SEI. [enabled]\n"
    "      --crypto <string>      : Selective encryption. Crypto support must be\n"
    "                               enabled at compile-time. Can be 'on' or 'off' or\n"
    "                               a list of features separated with a '+'. [off]\n"
    "                                   - on: Enable all encryption features.\n"
    "                                   - off: Disable selective encryption.\n"
    "                                   - mvs: Motion vector magnitudes.\n"
    "                                   - mv_signs: Motion vector signs.\n"
    "                                   - trans_coeffs: Coefficient magnitudes.\n"
    "                                   - trans_coeff_signs: Coefficient signs.\n"
    "                                   - intra_pred_modes: Intra prediction modes.\n"
    "      --key <string>         : Encryption key [16,213,27,56,255,127,242,112,\n"
    "                                               97,126,197,204,25,59,38,30]\n"
    "\n"
    /* Word wrap to this width to stay under 80 characters (including ") *************/
    "Video structure:\n"
    "  -q, --qp <integer>         : Quantization parameter [22]\n"
    "  -p, --period <integer>     : Period of intra pictures [64]\n"
    "                                   - 0: Only first picture is intra.\n"
    "                                   - 1: All pictures are intra.\n"
    "                                   - N: Every Nth picture is intra.\n"
    "      --vps-period <integer> : How often the video parameter set is re-sent [0]\n"
    "                                   - 0: Only send VPS with the first frame.\n"
    "                                   - N: Send VPS with every Nth intra frame.\n"
    "  -r, --ref <integer>        : Number of reference frames, in range 1..15 [4]\n"
    "      --gop <string>         : GOP structure [8]\n"
    "                                   - 0: Disabled\n"
    "                                   - 8: B-frame pyramid of length 8\n"
    "                                   - lp-<string>: Low-delay P-frame GOP\n"
    "                                     (e.g. lp-g8d4t2, see README)\n"
    "      --cqmfile <filename>   : Read custom quantization matrices from a file.\n"
    "      --bitrate <integer>    : Target bitrate [0]\n"
    "                                   - 0: Disable rate control.\n"
    "                                   - N: Target N bits per second.\n"
    "      --(no-)lossless        : Use lossless coding. [disabled]\n"
    "      --mv-constraint <string> : Constrain movement vectors. [none]\n"
    "                                   - none: No constraint\n"
    "                                   - frametile: Constrain within the tile.\n"
    "                                   - frametilemargin: Constrain even more.\n"
    "      --roi <filename>       : Use a delta QP map for region of interest.\n"
    "                               Reads an array of delta QP values from a text\n"
    "                               file. The file format is: width and height of\n"
    "                               the QP delta map followed by width*height delta\n"
    "                               QP values in raster order. The map can be of any\n"
    "                               size and will be scaled to the video size.\n"
    "      --(no-)erp-aqp         : Use adaptive QP for 360 degree video with\n"
    "                               equirectangular projection. [disabled]\n"
    "      --level <number>       : Use the given HEVC level in the output and give\n"
    "                               an error if level limits are exceeded. [6.2]\n"
    "                                   - 1, 2, 2.1, 3, 3.1, 4, 4.1, 5, 5.1, 5.2, 6,\n"
    "                                     6.1, 6.2\n"
    "      --force-level <number> : Same as --level but warnings instead of errors.\n"
    "      --high-tier            : Used with --level. Use high tier bitrate limits\n"
    "                               instead of the main tier limits during encoding.\n"
    "                               High tier requires level 4 or higher.\n"
    "\n"
    /* Word wrap to this width to stay under 80 characters (including ") *************/
    "Compression tools:\n"
    "      --(no-)deblock <beta:tc> : Deblocking filter. [0:0]\n"
    "                                   - beta: Between -6 and 6\n"
    "                                   - tc: Between -6 and 6\n"
    "      --sao <string>         : Sample Adaptive Offset [full]\n"
    "                                   - off: SAO disabled\n"
    "                                   - band: Band offset only\n"
    "                                   - edge: Edge offset only\n"
    "                                   - full: Full SAO\n"
    "      --(no-)rdoq            : Rate-distortion optimized quantization [enabled]\n"
    "      --(no-)rdoq-skip       : Skip RDOQ for 4x4 blocks. [disabled]\n"
    "      --(no-)signhide        : Sign hiding [disabled]\n"
    "      --(no-)smp             : Symmetric motion partition [disabled]\n"
    "      --(no-)amp             : Asymmetric motion partition [disabled]\n"
    "      --rd <integer>         : Intra mode search complexity [0]\n"
    "                                   - 0: Skip intra if inter is good enough.\n"
    "                                   - 1: Rough intra mode search with SATD.\n"
    "                                   - 2: Refine intra mode search with SSE.\n"
    "                                   - 3: Try all intra modes and enable intra\n"
    "                                        chroma mode search.\n"
    "      --(no-)mv-rdo          : Rate-distortion optimized motion vector costs\n"
    "                               [disabled]\n"
    "      --(no-)full-intra-search : Try all intra modes during rough search.\n"
    "                               [disabled]\n"
    "      --(no-)transform-skip  : Try transform skip [disabled]\n"
    "      --me <string>          : Integer motion estimation algorithm [hexbs]\n"
    "                                   - hexbs: Hexagon Based Search\n"
    "                                   - tz:    Test Zone Search\n"
    "                                   - full:  Full Search\n"
    "                                   - full8, full16, full32, full64\n"
    "                                   - dia:   Diamond Search\n"
    "      --me-steps <integer>   : Motion estimation search step limit. Only\n"
    "                               affects 'hexbs' and 'dia'. [-1]\n"
    "      --subme <integer>      : Fractional pixel motion estimation level [4]\n"
    "                                   - 0: Integer motion estimation only\n"
    "                                   - 1: + 1/2-pixel horizontal and vertical\n"
    "                                   - 2: + 1/2-pixel diagonal\n"
    "                                   - 3: + 1/4-pixel horizontal and vertical\n"
    "                                   - 4: + 1/4-pixel diagonal\n"
    "      --pu-depth-inter <int>-<int> : Inter prediction units sizes [0-3]\n"
    "                                   - 0, 1, 2, 3: from 64x64 to 8x8\n"
    "      --pu-depth-intra <int>-<int> : Intra prediction units sizes [1-4]\n"
    "                                   - 0, 1, 2, 3, 4: from 64x64 to 4x4\n"
    "      --tr-depth-intra <int> : Transform split depth for intra blocks [0]\n"
    "      --(no-)bipred          : Bi-prediction [disabled]\n"
    "      --cu-split-termination <string> : CU split search termination [zero]\n"
    "                                   - off: Don't terminate early.\n"
    "                                   - zero: Terminate when residual is zero.\n"
    "      --me-early-termination <string> : Motion estimation termination [on]\n"
    "                                   - off: Don't terminate early.\n"
    "                                   - on: Terminate early.\n"
    "                                   - sensitive: Terminate even earlier.\n"
    "      --fast-residual-cost <int> : Skip CABAC cost for residual coefficients\n"
    "                                   when QP is below the limit. [0]\n"
    "      --(no-)intra-rdo-et    : Check intra modes in rdo stage only until\n"
    "                               a zero coefficient CU is found. [disabled]\n"
    "      --(no-)implicit-rdpcm  : Implicit residual DPCM. Currently only supported\n"
    "                               with lossless coding. [disabled]\n"
    "      --(no-)tmvp            : Temporal motion vector prediction [enabled]\n"
    "\n"
    /* Word wrap to this width to stay under 80 characters (including ") *************/
    "Parallel processing:\n"
    "      --threads <integer>    : Number of threads to use [auto]\n"
    "                                   - 0: Process everything with main thread.\n"
    "                                   - N: Use N threads for encoding.\n"
    "                                   - auto: Select automatically.\n"
    "      --owf <integer>        : Frame-level parallelism [auto]\n"
    "                                   - N: Process N+1 frames at a time.\n"
    "                                   - auto: Select automatically.\n"
    "      --(no-)wpp             : Wavefront parallel processing. [enabled]\n"
    "                               Enabling tiles automatically disables WPP.\n"
    "                               To enable WPP with tiles, re-enable it after\n"
    "                               enabling tiles. Enabling wpp with tiles is,\n"
    "                               however, an experimental feature since it is\n"
    "                               not supported in any HEVC profile.\n"
    "      --tiles <int>x<int>    : Split picture into width x height uniform tiles.\n"
    "      --tiles-width-split <string>|u<int> :\n"
    "                                   - <string>: A comma-separated list of tile\n"
    "                                               column pixel coordinates.\n"
    "                                   - u<int>: Number of tile columns of uniform\n"
    "                                             width.\n"
    "      --tiles-height-split <string>|u<int> :\n"
    "                                   - <string>: A comma-separated list of tile row\n"
    "                                               column pixel coordinates.\n"
    "                                   - u<int>: Number of tile rows of uniform\n"
    "                                             height.\n"
    "      --slices <string>      : Control how slices are used.\n"
    "                                   - tiles: Put tiles in independent slices.\n"
    "                                   - wpp: Put rows in dependent slices.\n"
    "                                   - tiles+wpp: Do both.\n"
    "\n"
    /* Word wrap to this width to stay under 80 characters (including ") *************/
    "Video Usability Information:\n"
    "      --sar <width:height>   : Specify sample aspect ratio\n"
    "      --overscan <string>    : Specify crop overscan setting [undef]\n"
    "                                   - undef, show, crop\n"
    "      --videoformat <string> : Specify video format [undef]\n"
    "                                   - undef, component, pal, ntsc, secam, mac\n"
    "      --range <string>       : Specify color range [tv]\n"
    "                                   - tv, pc\n"
    "      --colorprim <string>   : Specify color primaries [undef]\n"
    "                                   - undef, bt709, bt470m, bt470bg,\n"
    "                                     smpte170m, smpte240m, film, bt2020\n"
    "      --transfer <string>    : Specify transfer characteristics [undef]\n"
    "                                   - undef, bt709, bt470m, bt470bg,\n"
    "                                     smpte170m, smpte240m, linear, log100,\n"
    "                                     log316, iec61966-2-4, bt1361e,\n"
    "                                     iec61966-2-1, bt2020-10, bt2020-12\n"
    "      --colormatrix <string> : Specify color matrix setting [undef]\n"
    "                                   - undef, bt709, fcc, bt470bg, smpte170m,\n"
    "                                     smpte240m, GBR, YCgCo, bt2020nc, bt2020c\n"
    "      --chromaloc <integer>  : Specify chroma sample location (0 to 5) [0]\n"
    "\n"
    /* Word wrap to this width to stay under 80 characters (including ") *************/
    "Scalability extension (experimental):\n"
    "      --layer                : Add a new layer (EL). Parameters for this layer\n"
    "                               should be given after this tag. The base layer\n"
    "                               (BL) is created implicitly and parameters for\n"
    "                               it should be given before the first --layer tag.\n"
    "      --layer-res <res>      : Specify the size for the current layer. If no\n"
    "                               size is given, the size is assumed to be the\n"
    "                               same as the input size. When the layer size and\n"
    "                               input size differ, the input is scaled to match\n"
    "                               the layer size.\n"
    "      --input-layer <integer>: Explicitly specify which input layer the current\n"
    "                               layer should use. If not given tries to infer\n"
    "                               it. If no input is given in the current layer,\n"
    "                               use the latest input layer. Declaring a new\n"
    "                               input will associate the it with the current\n"
    "                               layer if not set before.\n"
    "  -i, --input <string>       : Declare a new input layer and specify the\n"
    "                               filename that should be read. See Required\n"
    "                               section for more info.\n"
    "      --input-res <res>      : Specify the size for the latest input.\n"
    "                               See Required section for more info.\n"
    "      --ilr <integer>        : Set the number of inter layer reference frames\n"
    "                               for the current layer. Currently only the value\n"
    "                               of 1 is supported. [0]\n"
    "      --debug <string>       : Specify the filename for reconstruction\n"
    "                               output. Each layer will have a debug file\n"
    "                               assosiated with it in order. See Options section\n"
    "                               for more info.\n"
    "  note: The scalability functionality is still WIP, so compatibility with other\n"
    "        parameters is limited. At least GOP and Parallel tools are not supported\n"
    "        with Scalability features. Also, only a two layer (BL+EL) configuration\n"
    "        have been tested, with 2x spatial and PSNR scalability.\n"
    "\n"
    /* Word wrap to this width to stay under 80 characters (including ") ************/
    "Deprecated parameters: (might be removed at some point)\n"
    "  -w, --width <integer>       : Use --input-res.\n"
    "  -h, --height <integer>      : Use --input-res.\n");
}

// ***********************************************
// Modified for SHVC
void print_frame_info(const kvz_frame_info *const info,
                      const double frame_psnr[3],
                      const uint32_t bytes,
                      const bool print_psnr)
{
  fprintf(stderr, "POC %4d LId %d TId %d QP %2d (%c-frame) %10d bits",
          info->poc,
          info->lid,
          info->tid,
          info->qp,
          "BPI"[info->slice_type % 3],
          bytes << 3);
  if (print_psnr) {
    fprintf(stderr, " PSNR Y %2.4f U %2.4f V %2.4f",
            frame_psnr[0], frame_psnr[1], frame_psnr[2]);
  }

  if (info->slice_type != KVZ_SLICE_I) {
    // Print reference picture lists
    fprintf(stderr, " [L0 ");
    for (int j = 0; j < info->ref_list_len[0]; j++) {
      fprintf(stderr, "%d ", info->ref_list[0][j]);
    }
    fprintf(stderr, "] [L1 ");
    for (int j = 0; j < info->ref_list_len[1]; j++) {
      fprintf(stderr, "%d ", info->ref_list[1][j]);
    }
    fprintf(stderr, "]");
  }

  fprintf(stderr, "\n");
}
// ***********************************************