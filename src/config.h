#ifndef CONFIG_H_
#define CONFIG_H_
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
 * \brief Handles parsing and storing of configuration of the encoder.
 */

#include "kvazaar.h"
#include "global.h"

/* Function definitions */
kvz_config *config_alloc(void);
int config_init(kvz_config *cfg);
int config_destroy(kvz_config *cfg);
int config_read(kvz_config *cfg,int argc, char *argv[]);
int config_parse(kvz_config *cfg, const char *name, const char *value);
int config_validate(const kvz_config *cfg);
int config_set_owf_auto(kvz_config *cfg);

#endif
