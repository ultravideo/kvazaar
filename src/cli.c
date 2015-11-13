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
  { "width",              required_argument, NULL, 'w' },
  { "height",             required_argument, NULL, 'h' }, // deprecated
  { "frames",             required_argument, NULL, 'n' }, // deprecated
  { "qp",                 required_argument, NULL, 'q' },
  { "period",             required_argument, NULL, 'p' },
  { "ref",                required_argument, NULL, 'r' },
  { "vps-period",         required_argument, NULL, 0 },
  { "input-res",          required_argument, NULL, 0 },
  { "input-fps",          required_argument, NULL, 0 },
  { "deblock",                  no_argument, NULL, 0 },
  { "no-deblock",               no_argument, NULL, 0 },
  { "sao",                      no_argument, NULL, 0 },
  { "no-sao",                   no_argument, NULL, 0 },
  { "rdoq",                     no_argument, NULL, 0 },
  { "no-rdoq",                  no_argument, NULL, 0 },
  { "signhide",                 no_argument, NULL, 0 },
  { "no-signhide",              no_argument, NULL, 0 },
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
  { "tiles-width-split",  required_argument, NULL, 0 },
  { "tiles-height-split", required_argument, NULL, 0 },
  { "wpp",                      no_argument, NULL, 0 },
  { "no-wpp",                   no_argument, NULL, 0 },
  { "owf",                required_argument, NULL, 0 },
  { "slice-addresses",    required_argument, NULL, 0 },
  { "threads",            required_argument, NULL, 0 },
  { "cpuid",              required_argument, NULL, 0 },
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
    ok = 0;
    goto done;
  }

  opts->config = api->config_alloc();
  if (!opts->config || !api->config_init(opts->config)) {
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

    const char* name = long_options[long_options_index].name;
    if (!strcmp(name, "input")) {
      if (opts->input) {
        fprintf(stderr, "Input error: More than one input file given.\n");
        ok = 0;
        goto done;
      }
      opts->input = strdup(optarg);
    } else if (!strcmp(name, "output")) {
      if (opts->output) {
        fprintf(stderr, "Input error: More than one output file given.\n");
        ok = 0;
        goto done;
      }
      opts->output = strdup(optarg);
    } else if (!strcmp(name, "debug")) {
      if (opts->debug) {
        fprintf(stderr, "Input error: More than one debug output file given.\n");
        ok = 0;
        goto done;
      }
      opts->debug = strdup(optarg);
    } else if (!strcmp(name, "seek")) {
      opts->seek = atoi(optarg);
    } else if (!strcmp(name, "frames")) {
      opts->frames = atoi(optarg);
    } else if (!api->config_parse(opts->config, name, optarg)) {
      fprintf(stderr, "invalid argument: %s=%s\n", name, optarg);
      ok = 0;
      goto done;
    }
  }

  // Check for extra arguments.
  if (argc - optind > 0) {
    fprintf(stderr, "Input error: Extra argument found: \"%s\"\n", argv[optind]);
    ok = 0;
    goto done;
  }

  // Check that the required files were defined
  if (opts->input == NULL || opts->output == NULL) {
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

  // Set resolution automatically if necessary
  if (opts->config->width == 0 && opts->config->width == 0){
    ok = select_input_res_auto(opts->input, &opts->config->width, &opts->config->height);
    goto done;
  }

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
    FREE_POINTER(opts->input);
    FREE_POINTER(opts->output);
    FREE_POINTER(opts->debug);
    api->config_destroy(opts->config);
    opts->config = NULL;
  }
  FREE_POINTER(opts);
}


void print_version(void)
{
  fprintf(stderr,
    "/***********************************************/\n"
    " *   Kvazaar HEVC Encoder v. " VERSION_STRING "             *\n"
    " *     Tampere University of Technology 2015   *\n"
    "/***********************************************/\n\n");
}


void print_help(void)
{
  fprintf(stderr,
    "Usage:\n"
    "kvazaar -i <input> --input-res <width>x<height> -o <output>\n"
    "\n"
    "Optional parameters:\n"
    "      -n, --frames <integer>     : Number of frames to code [all]\n"
    "      --seek <integer>           : First frame to code [0]\n"
    "      --input-res <int>x<int>    : Input resolution (width x height) or\n"
    "                  auto           : try to detect from file name [auto]\n"
    "      --input-fps <number>       : Framerate of the input video [25.0]\n"
    "      -q, --qp <integer>         : Quantization Parameter [32]\n"
    "      -p, --period <integer>     : Period of intra pictures [0]\n"
    "                                     0: only first picture is intra\n"
    "                                     1: all pictures are intra\n"
    "                                     2-N: every Nth picture is intra\n"
    "          --vps-period <integer> : Specify how often the video parameter set is\n"
    "                                   re-sent. [0]\n"
    "                                     0: only send VPS with the first frame\n"
    "                                     1: send VPS with every intra frame\n"
    "                                     N: send VPS with every Nth intra frame\n"
    "      -r, --ref <integer>        : Reference frames, range 1..15 [3]\n"
    "          --no-deblock           : Disable deblocking filter\n"
    "          --deblock <beta:tc>    : Deblocking filter parameters\n"
    "                                   beta and tc range is -6..6 [0:0]\n"
    "          --no-sao               : Disable sample adaptive offset\n"
    "          --no-rdoq              : Disable RDO quantization\n"
    "          --no-signhide          : Disable sign hiding in quantization\n"
    "          --rd <integer>         : Rate-Distortion Optimization level [1]\n"
    "                                     0: no RDO\n"
    "                                     1: estimated RDO\n"
    "                                     2: full RDO\n"
    "          --mv-rdo               : Enable Rate-Distortion Optimized motion vector costs\n"
    "          --full-intra-search    : Try all intra modes.\n"
    "          --me <string>          : Set integer motion estimation algorithm [\"hexbs\"]\n"
    "                                     \"hexbs\": Hexagon Based Search (faster)\n"
    "                                     \"tz\":    Test Zone Search (better quality)\n"
    "          --no-transform-skip    : Disable transform skip\n"
    "          --aud                  : Use access unit delimiters\n"
    "          --cqmfile <string>     : Custom Quantization Matrices from a file\n"
    "          --debug <string>       : Output encoders reconstruction.\n"
    "          --cpuid <integer>      : Disable runtime cpu optimizations with value 0.\n"
    "          --subme <integer>      : Set fractional pixel motion estimation level [1].\n"
    "                                     0: only integer motion estimation\n"
    "                                     1: fractional pixel motion estimation enabled\n"
    "          --source-scan-type <string> : Set source scan type [\"progressive\"].\n"
    "                                     \"progressive\": progressive scan\n"
    "                                     \"tff\": top field first\n"
    "                                     \"bff\": bottom field first\n"
    "          --pu-depth-inter <int>-<int> : Range for sizes of inter prediction units to try.\n"
    "                                     0: 64x64, 1: 32x32, 2: 16x16, 3: 8x8\n"
    "          --pu-depth-intra <int>-<int> : Range for sizes of intra prediction units to try.\n"
    "                                     0: 64x64, 1: 32x32, 2: 16x16, 3: 8x8, 4: 4x4\n"
    "          --no-info              : Don't add information about the encoder to settings.\n"
    "          --gop <int>           : Length of Group of Pictures, must be 8 or 0 [0]\n"
    "          --bipred               : Enable bi-prediction search\n"
    "          --bitrate <integer>    : Target bitrate. [0]\n"
    "                                     0: disable rate-control\n"
    "                                     N: target N bits per second\n"
    "          --preset <string>      : Use preset\n"
    "                                     ultrafast, superfast,veryfast, faster,\n"
    "                                     fast, medium, slow, slower, veryslow, placebo\n"
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
    "  Parallel processing:\n"
    "          --threads <integer>    : Maximum number of threads to use.\n"
    "                                   Disable threads if set to 0.\n"
    "\n"
    "  Tiles:\n"
    "          --tiles-width-split <string>|u<int> : \n"
    "                                   Specifies a comma separated list of pixel\n"
    "                                   positions of tiles columns separation coordinates.\n"
    "                                   Can also be u followed by and a single int n,\n"
    "                                   in which case it produces columns of uniform width.\n"
    "          --tiles-height-split <string>|u<int> : \n"
    "                                   Specifies a comma separated list of pixel\n"
    "                                   positions of tiles rows separation coordinates.\n"
    "                                   Can also be u followed by and a single int n,\n"
    "                                   in which case it produces rows of uniform height.\n"
    "\n"
    "  Wpp:\n"
    "          --wpp                  : Enable wavefront parallel processing\n"
    "          --owf <integer>|auto   : Number of parallel frames to process. 0 to disable.\n"
    "\n"
    "  Slices:\n"
    "          --slice-addresses <string>|u<int>: \n"
    "                                   Specifies a comma separated list of LCU\n"
    "                                   positions in tile scan order of tile separations.\n"
    "                                   Can also be u followed by and a single int n,\n"
    "                                   in which case it produces uniform slice length.\n"
    "\n"
    "  Deprecated parameters: (might be removed at some point)\n"
    "     Use --input-res:\n"
    "       -w, --width               : Width of input in pixels\n"
    "       -h, --height              : Height of input in pixels\n");
}


void print_frame_info(const kvz_frame_info *const info,
                      const double frame_psnr[3],
                      const uint32_t bytes)
{
  fprintf(stderr, "POC %4d QP %2d (%c-frame) %10d bits PSNR: %2.4f %2.4f %2.4f",
          info->poc,
          info->qp,
          "BPI"[info->slice_type % 3],
          bytes << 3,
          frame_psnr[0], frame_psnr[1], frame_psnr[2]);

  if (info->slice_type != KVZ_SLICE_I) {
    // Print reference picture lists
    fprintf(stderr, " [L0 ");
    for (int j = info->ref_list_len[0] - 1; j >= 0; j--) {
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
