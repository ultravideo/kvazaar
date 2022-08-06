#ifndef CLI_H_
#define CLI_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

/**
 * \file
 * Command line interface
 */

#include "global.h" // IWYU pragma: keep
#include "kvazaar.h"

typedef struct cmdline_opts_t {
  /** \brief Input filename */
  char *input;
  /** \brief Output filename */
  char *output;
  /** \brief Debug output */
  char *debug;
  /** \brief Number of input frames to skip */
  int32_t seek;
  /** \brief Number of frames to encode */
  int32_t frames;
  /** \brief Encoder configuration */
  kvz_config *config;
  /** \brief Show help message and exit */
  bool help;
  /** \brief Show version information and exit */
  bool version;
  /** \brief Whether to loop input */
  bool loop_input;
} cmdline_opts_t;

cmdline_opts_t* cmdline_opts_parse(const kvz_api *api, int argc, char *argv[]);
void cmdline_opts_free(const kvz_api *api, cmdline_opts_t *opts);

void print_usage(void);
void print_version(void);
void print_help(void);
void print_frame_info(const kvz_frame_info *const info,
                      const double frame_psnr[3],
                      const uint32_t bytes,
                      const bool print_psnr,
                      const double avg_qp);

#endif
