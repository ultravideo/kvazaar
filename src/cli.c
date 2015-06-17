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

#include "config.h"

#include <stdio.h>

#include "encoderstate.h"


void print_version(void)
{
  fprintf(stderr,
    "/***********************************************/\n"
    " *   Kvazaar HEVC Encoder v. " VERSION_STRING "             *\n"
    " *     Tampere University of Technology 2014   *\n"
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
    "      --input-res <int>x<int>    : Input resolution (width x height)\n"
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


void print_frame_info(encoder_state_t *state, double frame_psnr[3])
{
  fprintf(stderr, "POC %4d QP %2d (%c-frame) %10d bits PSNR: %2.4f %2.4f %2.4f",
          state->global->poc,
          state->global->QP,
          "BPI"[state->global->slicetype % 3], state->stats_bitstream_length << 3,
          frame_psnr[0], frame_psnr[1], frame_psnr[2]);

  // Print reference picture lists
  if (state->global->slicetype != SLICE_I) {
    int j, ref_list[2] = { 0, 0 }, ref_list_poc[2][16];
    // List all pocs of lists
    for (j = 0; j < state->global->ref->used_size; j++) {
      if (state->global->ref->pocs[j] < state->global->poc) {
        ref_list_poc[0][ref_list[0]] = state->global->ref->pocs[j];
        ref_list[0]++;
      } else {
        ref_list_poc[1][ref_list[1]] = state->global->ref->pocs[j];
        ref_list[1]++;
      }
    }
    encoder_ref_insertion_sort(ref_list_poc[0], ref_list[0]);
    encoder_ref_insertion_sort(ref_list_poc[1], ref_list[1]);

    fprintf(stderr, " [L0 ");
    for (j = ref_list[0] - 1; j >= 0; j--) {
      fprintf(stderr, "%d ", ref_list_poc[0][j]);
    }
    fprintf(stderr, "] [L1 ");
    for (j = 0; j < ref_list[1]; j++) {
      fprintf(stderr, "%d ", ref_list_poc[1][j]);
    }
    fprintf(stderr, "]");
  }

  fprintf(stderr, "\n");
}
