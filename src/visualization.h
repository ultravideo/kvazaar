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

#if KVZ_VISUALIZATION == 1

#include <SDL.h>
#include <SDL_ttf.h>
#include <math.h>

#include "encoderstate.h"
#include "threadqueue.h"
#include "cu.h"

static int INFO_WIDTH = 480;
static int INFO_HEIGHT = 240;

extern SDL_Renderer *renderer, *info_renderer;
extern SDL_Window *window, *info_window;
extern SDL_Surface *screen, *pic;
extern SDL_Texture *overlay, *overlay_blocks, *overlay_intra, *overlay_inter[2], *overlay_hilight;
extern int screen_w, screen_h;
extern int sdl_draw_blocks;
extern int sdl_draw_intra;
extern int sdl_block_info;
extern pthread_mutex_t sdl_mutex;
extern kvz_pixel *sdl_pixels_hilight;
extern kvz_pixel *sdl_pixels_RGB;
extern kvz_pixel *sdl_pixels_RGB_intra_dir;
extern kvz_pixel *sdl_pixels_RGB_inter[2];
extern kvz_pixel *sdl_pixels;
extern kvz_pixel *sdl_pixels_u;
extern kvz_pixel *sdl_pixels_v;
extern int32_t sdl_delay;
extern SDL_Surface *textSurface;
extern SDL_Texture *text;

extern cu_info_t *sdl_cu_array;
extern TTF_Font *font;


void *eventloop_main(void *temp);

void kvz_visualization_init(int width, int height);
void kvz_visualization_free();

void kvz_visualization_frame_init(encoder_control_t *encoder, kvz_picture *img_in);

void kvz_visualization_draw_block(const encoder_state_t *state, lcu_t *lcu, cu_info_t *cur_cu, int x, int y, int depth);

#define PUTPIXEL_Y(pixel_x, pixel_y, color_y) sdl_pixels_RGB[luma_index + (pixel_x) + (pixel_y)*pic_width] = color_y;
#define PUTPIXEL_U(pixel_x, pixel_y, color_u) sdl_pixels_u[chroma_index + (pixel_x>>1) + (pixel_y>>1)*(pic_width>>1)] = color_u;
#define PUTPIXEL_V(pixel_x, pixel_y, color_v) sdl_pixels_v[chroma_index + (pixel_x>>1) + (pixel_y>>1)*(pic_width>>1)] = color_v;
#define PUTPIXEL(pixel_x, pixel_y, color_r, color_g, color_b, color_alpha) sdl_pixels_RGB[index_RGB + (pixel_x<<2) + (pixel_y)*(pic_width<<2)+3] = color_alpha; \
  sdl_pixels_RGB[index_RGB + (pixel_x<<2) + (pixel_y)*(pic_width<<2) +2] = color_r; \
  sdl_pixels_RGB[index_RGB + (pixel_x<<2) + (pixel_y)*(pic_width<<2) +1] = color_g; \
  sdl_pixels_RGB[index_RGB + (pixel_x<<2) + (pixel_y)*(pic_width<<2) +0] = color_b;

#define PUTPIXEL_intra(pixel_x, pixel_y, color_r, color_g, color_b, color_alpha) sdl_pixels_RGB_intra_dir[index_RGB + (pixel_x<<2) + (pixel_y)*(pic_width<<2)+3] = color_alpha; \
  sdl_pixels_RGB_intra_dir[index_RGB + (pixel_x<<2) + (pixel_y)*(pic_width<<2) +2] = color_r; \
  sdl_pixels_RGB_intra_dir[index_RGB + (pixel_x<<2) + (pixel_y)*(pic_width<<2) +1] = color_g; \
  sdl_pixels_RGB_intra_dir[index_RGB + (pixel_x<<2) + (pixel_y)*(pic_width<<2) +0] = color_b;
#define PUTPIXEL_YUV(pixel_x, pixel_y, color_y, color_u, color_v) PUTPIXEL_Y(pixel_x,pixel_y, color_y); PUTPIXEL_U(pixel_x,pixel_y, color_u); PUTPIXEL_V(pixel_x,pixel_y, color_v);

static INLINE void kvz_putpixel(kvz_pixel *buffer, int pic_width, int x, int y, kvz_pixel color_r, kvz_pixel color_g, kvz_pixel color_b, kvz_pixel color_a)
{
  int index = (x + y * pic_width) * 4;
  buffer[index + 0] = color_b;
  buffer[index + 1] = color_g;
  buffer[index + 2] = color_r;
  buffer[index + 3] = color_a;
}

static const uint32_t frame_r[8] = { 0, 128, 100, 128, 255, 128, 255, 128 };
static const uint32_t frame_g[8] = { 255, 128, 100, 128, 255, 128, 0, 128 };
static const uint32_t frame_b[8] = { 0, 128, 255, 128, 0, 128, 100, 128 };

static void draw_line(int pic_width, int index_RGB, int x1, int y1, int x2, int y2, int color_r, int color_g, int color_b)
{
  float temp_x = x1; float temp_y = y1;
  int len = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
  if (len > 0) {
    float x_off = ((float)(x2 - x1)) / (float)len;
    float y_off = (y2 - y1) / (float)len;
    for (int i = 0; i < len; i++) {
      int xx1 = temp_x;
      int yy1 = temp_y;
      PUTPIXEL(xx1, yy1, color_r, color_g, color_b, 255);
      temp_x += x_off; temp_y += y_off;
    }
  }
}

static void draw_mv(kvz_pixel *buffer, int pic_width, int pic_height, int x1, int y1, int x2, int y2, int color_r, int color_g, int color_b)
{
  int frac_x = x1 << 4;
  int frac_y = y1 << 4;
  int num_samples = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
  num_samples = num_samples == 0 ? 1 : num_samples;

  const int x_off = ((x2 - x1) << 4) / num_samples;
  const int y_off = ((y2 - y1) << 4) / num_samples;
  for (int i = 0; i <= num_samples; ++i, frac_x += x_off, frac_y += y_off) {
    int x = frac_x >> 4;
    int y = frac_y >> 4;
    if (x < 0 || x >= pic_width || y < 0 || y >= pic_height) break;
    kvz_putpixel(buffer, pic_width, x, y, color_r, color_g, color_b, 255);
  }

  // Mark the origin of the motion vector.
  kvz_putpixel(buffer, pic_width, x1, y1, 255 - color_r, 255 - color_g, 255 - color_b, 255);
}

static void kvz_visualization_delay()
{
  volatile int64_t i = 0;
  if (sdl_delay) {
    //SDL_Delay(sdl_delay);
    int64_t wait_cycles = pow(2, sdl_delay) * 1000;
    while (i < wait_cycles) {
      ++i;
    }
  }
}

#endif

#endif // VISUALIZAITON_H_
