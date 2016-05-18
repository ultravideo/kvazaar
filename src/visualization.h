#ifndef VISUALIZAITON_H_
#define VISUALIZAITON_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (C) 2013-2016 Tampere University of Technology and others (see
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

#include "global.h"

#if KVZ_VISUALIZATION

#include <SDL.h>
#include <SDL_ttf.h>
#include <math.h>

#include "encoderstate.h"
#include "threads.h"
#include "threadqueue.h"
#include "cu.h"

void *kvz_visualization_eventloop(void *temp);

void kvz_visualization_init(int width, int height);
void kvz_visualization_free();

void kvz_visualization_frame_init(encoder_control_t *encoder, kvz_picture *img_in);

bool kvz_visualization_draw_block(const encoder_state_t *state, lcu_t *lcu, cu_info_t *cur_cu, int x, int y, int depth);
bool kvz_visualization_draw_block_with_delay(const encoder_state_t *state, lcu_t *lcu, cu_info_t *cur_cu, int x, int y, int depth);

void kvz_visualization_mv_draw_lcu(encoder_state_t * const state, int x, int y, lcu_t *lcu);
void kvz_visualization_mv_clear_lcu(encoder_state_t * const state, int x, int y);

#endif // KVZ_VISUALIZATION

#endif // VISUALIZAITON_H_
