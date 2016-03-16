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


#if KVZ_VISUALIZATION == 1

#include "threadqueue.h"
#include <SDL.h>
extern SDL_Renderer *renderer;
extern SDL_Surface *screen, *pic;
extern SDL_Texture *overlay, *overlay_blocks, *overlay_intra;
extern int screen_w, screen_h;
extern int sdl_draw_blocks;
extern pthread_mutex_t sdl_mutex;
extern kvz_pixel *sdl_pixels_RGB;
extern kvz_pixel *sdl_pixels_RGB_intra_dir;
extern kvz_pixel *sdl_pixels;
extern kvz_pixel *sdl_pixels_u;
extern kvz_pixel *sdl_pixels_v;
extern int32_t sdl_delay;
extern cu_info_t *sdl_cu_array;

#define PTHREAD_LOCK(l) if (pthread_mutex_lock((l)) != 0) { fprintf(stderr, "pthread_mutex_lock(%s) failed!\n", #l); assert(0); return 0; }
#define PTHREAD_UNLOCK(l) if (pthread_mutex_unlock((l)) != 0) { fprintf(stderr, "pthread_mutex_unlock(%s) failed!\n", #l); assert(0); return 0; }

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
#endif

#endif // VISUALIZAITON_H_
