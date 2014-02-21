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
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "extras/getopt.h"

/**
 * \brief Allocate memory for config object
 * \return pointer to allocated memory
 */
config *config_alloc(void)
{
  config *cfg = (config *)malloc(sizeof(config));
  if (!cfg) {
    fprintf(stderr, "Failed to allocate a config object!\n");
    return cfg;
  }

  memset(cfg, 0, sizeof(config));
  return cfg;
}

/**
 * \brief Initialize config structure
 * \param cfg config object
 * \return 1 on success, 0 on failure
 */
int config_init(config *cfg)
{
  cfg->input           = NULL;
  cfg->output          = NULL;
  cfg->debug           = NULL;
  cfg->frames          = 0;
  cfg->width           = 320;
  cfg->height          = 240;
  cfg->qp              = 32;
  cfg->intra_period    = 0;
  cfg->deblock_enable  = 1;
  cfg->deblock_beta    = 0;
  cfg->deblock_tc      = 0;
  cfg->sao_enable      = 1;
  cfg->vui.sar_width   = 0;
  cfg->vui.sar_height  = 0;
  cfg->vui.overscan    = 0; /* undef */
  cfg->vui.videoformat = 5; /* undef */
  cfg->vui.fullrange   = 0; /* limited range */
  cfg->vui.colorprim   = 2; /* undef */
  cfg->vui.transfer    = 2; /* undef */
  cfg->vui.colormatrix = 2; /* undef */
  cfg->vui.chroma_loc  = 0; /* left center */
  cfg->aud_enable      = 0;
  cfg->cqmfile         = NULL;
  cfg->ref_frames      = DEFAULT_REF_PIC_COUNT;

  return 1;
}

/**
 * \brief Free memory allocated to the config
 * \param cfg config object
 * \return 1 on success, 0 on failure
 */
int config_destroy(config *cfg)
{
  FREE_POINTER(cfg->input);
  FREE_POINTER(cfg->output);
  FREE_POINTER(cfg->cqmfile);
  free(cfg);

  return 1;
}

/**
 * \brief Allocates memory space for a string, and copies it
 * \param char * string to copy
 * \return a pointer to the copied string on success, null on failure
 */
static char *copy_string(const char *string)
{
  // Allocate +1 for \0
  char *allocated_string = (char *)malloc(strlen(string) + 1);
  if (!allocated_string) {
    fprintf(stderr, "Failed to allocate a string!\n");
    return allocated_string;
  }

  // Copy the string to the new buffer
  memcpy(allocated_string, string, strlen(string) + 1);

  return allocated_string;
}

static int atobool(const char *str)
{
  if (!strcmp(str, "1")    ||
      !strcmp(str, "true") ||
      !strcmp(str, "yes"))
    return 1;
  if (!strcmp(str, "0")     ||
      !strcmp(str, "false") ||
      !strcmp(str, "no"))
    return 0;
  return 0;
}

static int parse_enum(const char *arg, const char * const *names, int8_t *dst)
{
  int i;
  for (i = 0; names[i]; i++)
    if (!strcmp(arg, names[i])) {
      *dst = i;
      return 1;
    }

  return 0;
}

static int config_parse(config *cfg, const char *name, const char *value)
{
  static const char * const overscan_names[]    = { "undef", "show", "crop", NULL };
  static const char * const videoformat_names[] = { "component", "pal", "ntsc", "secam", "mac", "undef", NULL };
  static const char * const range_names[]       = { "tv", "pc", NULL };
  static const char * const colorprim_names[]   = { "", "bt709", "undef", "", "bt470m", "bt470bg", "smpte170m",
                                                    "smpte240m", "film", "bt2020", NULL };
  static const char * const transfer_names[]    = { "", "bt709", "undef", "", "bt470m", "bt470bg", "smpte170m",
                                                    "smpte240m", "linear", "log100", "log316", "iec61966-2-4",
                                                    "bt1361e", "iec61966-2-1", "bt2020-10", "bt2020-12", NULL };
  static const char * const colormatrix_names[] = { "GBR", "bt709", "undef", "", "fcc", "bt470bg", "smpte170m",
                                                    "smpte240m", "YCgCo", "bt2020nc", "bt2020c", NULL };

  int i;
  int error = 0;

  if (!name)
    return 0;
  if (!value)
    value = "true";

  if ((!strncmp(name, "no-", 3) && (i = 3))) {
    name += i;
    value = atobool(value) ? "false" : "true";
  }

#define OPT(STR) else if (!strcmp(name, STR))
  if (0);
  OPT("input")
    cfg->input = copy_string(value);
  OPT("output")
    cfg->output = copy_string(value);
  OPT("debug")
    cfg->debug = copy_string(value);
  OPT("width")
    cfg->width = atoi(value);
  OPT("height")
    cfg->height = atoi(value);
  OPT("frames")
    cfg->frames = atoi(value);
  OPT("qp")
    cfg->qp = atoi(value);
  OPT("period")
    cfg->intra_period = atoi(value);
  OPT("ref") {
    cfg->ref_frames = atoi(value);
    if (cfg->ref_frames  < 1 || cfg->ref_frames >= MAX_REF_PIC_COUNT) {
      fprintf(stderr, "--ref out of range [1..15], set to 3\n");
      cfg->ref_frames = 3;
    }
  }
  OPT("deblock") {
    int beta, tc;
    if (2 == sscanf(value, "%d:%d", &beta, &tc)) {
      cfg->deblock_enable = 1;
      cfg->deblock_beta   = beta;
      cfg->deblock_tc     = tc;
    } else if (sscanf(value, "%d", &beta)) {
      cfg->deblock_enable = 1;
      cfg->deblock_beta   = beta;
      cfg->deblock_tc     = cfg->deblock_beta;
    } else
      cfg->deblock_enable = atobool(value);

    if (cfg->deblock_beta  < -6 || cfg->deblock_beta  > 6) {
      fprintf(stderr, "--deblock beta parameter out of range [-6..6], set to 0\n");
      cfg->deblock_beta = 0;
    }
    if (cfg->deblock_tc < -6 || cfg->deblock_tc > 6) {
      fprintf(stderr, "--deblock tc parameter out of range [-6..6], set to 0\n");
      cfg->deblock_tc = 0;
    }
  }
  OPT("sao")
    cfg->sao_enable = atobool(value);
  OPT("sar") {
      int sar_width, sar_height;
      if (2 == sscanf(value, "%d:%d", &sar_width, &sar_height)) {
        cfg->vui.sar_width  = sar_width;
        cfg->vui.sar_height = sar_height;
      } else
        error = 1;
  }
  OPT("overscan")
    error = !parse_enum(value, overscan_names, &cfg->vui.overscan);
  OPT("videoformat")
    error = !parse_enum(value, videoformat_names, &cfg->vui.videoformat);
  OPT("range")
    error = !parse_enum(value, range_names, &cfg->vui.fullrange);
  OPT("colorprim")
    error = !parse_enum(value, colorprim_names, &cfg->vui.colorprim);
  OPT("transfer")
    error = !parse_enum(value, transfer_names, &cfg->vui.transfer);
  OPT("colormatrix")
    error = !parse_enum(value, colormatrix_names, &cfg->vui.colormatrix);
  OPT("chromaloc") {
      cfg->vui.chroma_loc = atoi(value);
      if (cfg->vui.chroma_loc < 0 || cfg->vui.chroma_loc > 5) {
        fprintf(stderr, "--chromaloc parameter out of range [0..5], set to 0\n");
        cfg->vui.chroma_loc = 0;
      }
  }
  OPT("aud")
    cfg->aud_enable = atobool(value);
  OPT("cqmfile")
    cfg->cqmfile = copy_string(value);
  else
    return 0;
#undef OPT

  return error ? 0 : 1;
}

/**
 * \brief Read configuration options from argv to the config struct
 * \param cfg config object
 * \param argc argument count
 * \param argv argument list
 * \return 1 on success, 0 on failure
 */
int config_read(config *cfg,int argc, char *argv[])
{
  static char short_options[] = "i:o:d:w:h:n:q:p:r:";
  static struct option long_options[] =
  {
    { "input",              required_argument, NULL, 'i' },
    { "output",             required_argument, NULL, 'o' },
    { "debug",              required_argument, NULL, 'd' },
    { "width",              required_argument, NULL, 'w' },
    { "height",             required_argument, NULL, 'h' },
    { "frames",             required_argument, NULL, 'n' },
    { "qp",                 required_argument, NULL, 'q' },
    { "period",             required_argument, NULL, 'p' },
    { "ref",                required_argument, NULL, 'r' },
    { "no-deblock",               no_argument, NULL, 0 },
    { "deblock",            required_argument, NULL, 0 },
    { "no-sao",                   no_argument, NULL, 0 },
    { "sar",                required_argument, NULL, 0 },
    { "overscan",           required_argument, NULL, 0 },
    { "videoformat",        required_argument, NULL, 0 },
    { "range",              required_argument, NULL, 0 },
    { "colorprim",          required_argument, NULL, 0 },
    { "transfer",           required_argument, NULL, 0 },
    { "colormatrix",        required_argument, NULL, 0 },
    { "chromaloc",          required_argument, NULL, 0 },
    { "aud",                      no_argument, NULL, 0 },
    { "cqmfile",            required_argument, NULL, 0 },
    {0, 0, 0, 0}
  };

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
        return 0;
      }
    }

    if (!config_parse(cfg, long_options[long_options_index].name, optarg)) {
      const char *name = long_options_index > 0 ? long_options[long_options_index].name : argv[optind-2];
      fprintf(stderr, "invalid argument: %s = %s\r\n", name, optarg );
      return 0;
    }
  }

  // Check that the required files were defined
  if(cfg->input == NULL || cfg->output == NULL) return 0;

  return 1;
}
