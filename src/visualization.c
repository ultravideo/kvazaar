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

/*
 * \file Realtime SDL visualizer code.
 */

#include "visualization.h"

#if KVZ_VISUALIZATION

#include <SDL.h>
#include <SDL_ttf.h>
#include <math.h>

#include "cu.h"
#include "encoderstate.h"
#include "threads.h"
#include "threadqueue.h"


static int INFO_WIDTH = 480;
static int INFO_HEIGHT = 240;
static SDL_Renderer *renderer, *info_renderer;
static SDL_Window *window, *info_window = NULL;
static SDL_Texture *overlay, *overlay_blocks, *overlay_intra, *overlay_inter[2], *overlay_hilight;
static int screen_w, screen_h;
static int sdl_draw_blocks = 1;
static int sdl_draw_intra = 1;
static int sdl_block_info = 0;
static pthread_mutex_t sdl_mutex;
static kvz_pixel *sdl_pixels_hilight;
static kvz_pixel *sdl_pixels_RGB;
static kvz_pixel *sdl_pixels_RGB_intra_dir;
static kvz_pixel *sdl_pixels_RGB_inter[2];
static kvz_pixel *sdl_pixels;
static kvz_pixel *sdl_pixels_u;
static kvz_pixel *sdl_pixels_v;
static int32_t sdl_delay;
static SDL_Surface* textSurface;
static SDL_Texture* text;
static cu_info_t *sdl_cu_array;
static TTF_Font* font;

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
    // Busy loop, because normal sleep is not fine grained enough.
    int64_t wait_cycles = pow(2, sdl_delay) * 1000;
    while (i < wait_cycles) {
      ++i;
    }
  }
}

static void sdl_force_redraw(int locked)
{
  if (locked) {
    SDL_RenderClear(renderer);
    SDL_RenderCopy(renderer, overlay, NULL, NULL);
    if (sdl_draw_blocks)
      SDL_RenderCopy(renderer, overlay_blocks, NULL, NULL);
    if (sdl_draw_intra) {
      SDL_RenderCopy(renderer, overlay_intra, NULL, NULL);
      SDL_RenderCopy(renderer, overlay_inter[0], NULL, NULL);
      SDL_RenderCopy(renderer, overlay_inter[1], NULL, NULL);
    }


    if (sdl_block_info) {
      SDL_RenderCopy(renderer, overlay_hilight, NULL, NULL);
    }
  }
}

static void sdl_render_multiline_text(char* text, int x, int y)
{
  SDL_Color White = { 255, 255, 255 };
  SDL_Surface* temp_surface = TTF_RenderText_Solid(font, text, White);
  SDL_Rect src, dst;
  src.x = 0; src.y = 0; src.w = temp_surface->w; src.h = temp_surface->h;
  dst.x = x; dst.y = y; dst.w = temp_surface->w; dst.h = temp_surface->h;
  SDL_BlitSurface(temp_surface, &src, textSurface, &dst);
}

static void *kvz_visualization_eventloop(void* temp)
{

  int sdl_fader = 0;
  int sdl_faded_overlay = 0;

  int mouse_x = 0, mouse_y = 0;

  /* Initialize the display */

  window = SDL_CreateWindow(
    "Kvazaar",                  // window title
    SDL_WINDOWPOS_UNDEFINED,           // initial x position
    SDL_WINDOWPOS_UNDEFINED,           // initial y position
    screen_w, screen_h,
    SDL_WINDOW_RESIZABLE
    );

  if (window == NULL) {
    fprintf(stderr, "Couldn't set %dx%dx%d video mode: %s\n",
      screen_w, screen_h, 0, SDL_GetError());
    SDL_Quit();
    exit(1);
  }

  sdl_delay = 0;

  int height_in_lcu = screen_h / LCU_WIDTH;
  int width_in_lcu = screen_w / LCU_WIDTH;

  // Add one extra LCU when image not divisible by LCU_WIDTH
  if (height_in_lcu * LCU_WIDTH < screen_h) {
    height_in_lcu++;
  }
  if (width_in_lcu * LCU_WIDTH < screen_w) {
    width_in_lcu++;
  }

  // LCUs are split into 4x4 regions. Note that smallest "CU" is here 4x4.
  sdl_cu_array = MALLOC(cu_info_t, (height_in_lcu * LCU_CU_WIDTH)*(width_in_lcu * LCU_CU_WIDTH));
  FILL_ARRAY(sdl_cu_array, 0, (height_in_lcu * LCU_CU_WIDTH)*(width_in_lcu * LCU_CU_WIDTH));

  // Set the window manager title bar
  renderer = SDL_CreateRenderer(window, -1, 0);
  // Create overlays for reconstruction and reconstruction with block borders
  overlay = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_IYUV, SDL_TEXTUREACCESS_STREAMING, screen_w, screen_h);
  overlay_blocks = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, screen_w, screen_h);
  overlay_intra = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, screen_w, screen_h);
  overlay_inter[0] = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, screen_w, screen_h);
  overlay_inter[1] = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, screen_w, screen_h);
  overlay_hilight = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, screen_w, screen_h);


  SDL_SetTextureBlendMode(overlay_hilight, SDL_BLENDMODE_BLEND);
  SDL_SetTextureBlendMode(overlay_intra, SDL_BLENDMODE_BLEND);
  SDL_SetTextureBlendMode(overlay_inter[0], SDL_BLENDMODE_BLEND);
  SDL_SetTextureBlendMode(overlay_inter[1], SDL_BLENDMODE_BLEND);
  SDL_SetTextureBlendMode(overlay_blocks, SDL_BLENDMODE_BLEND);
  SDL_SetTextureBlendMode(overlay, SDL_BLENDMODE_BLEND);
  sdl_pixels_RGB = (kvz_pixel*)malloc(screen_w*screen_h * 4);
  memset(sdl_pixels_RGB, 0, (screen_w*screen_h * 4));

  sdl_pixels_hilight = (kvz_pixel*)malloc(screen_w*screen_h * 4);
  memset(sdl_pixels_hilight, 0, (screen_w*screen_h * 4));

  sdl_pixels_RGB_intra_dir = (kvz_pixel*)malloc(screen_w*screen_h * 4);
  memset(sdl_pixels_RGB_intra_dir, 0, (screen_w*screen_h * 4));

  sdl_pixels_RGB_inter[0] = (kvz_pixel*)malloc(screen_w*screen_h * 4);
  memset(sdl_pixels_RGB_inter[0], 0, (screen_w*screen_h * 4));
  sdl_pixels_RGB_inter[1] = (kvz_pixel*)malloc(screen_w*screen_h * 4);
  memset(sdl_pixels_RGB_inter[1], 0, (screen_w*screen_h * 4));

  sdl_pixels = (kvz_pixel*)malloc(screen_w*screen_h * 2 * sizeof(kvz_pixel));
  sdl_pixels_u = sdl_pixels + screen_w*screen_h;
  sdl_pixels_v = sdl_pixels_u + (screen_w*screen_h >> 2);

  if (overlay == NULL || overlay_blocks == NULL) {
    fprintf(stderr, "Couldn't create overlay: %s\n", SDL_GetError());
    SDL_Quit();
    exit(1);
  }

  if (TTF_Init() == -1) { printf("SDL_ttf could not initialize! SDL_ttf Error: %s\n", TTF_GetError()); }

  font = TTF_OpenFont("arial.ttf", 24);
  if (!font) {
    printf("TTF_OpenFont: %s\n", TTF_GetError());
    // handle error
  }

  Uint32 rmask, gmask, bmask, amask;
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
  rmask = 0xff000000;
  gmask = 0x00ff0000;
  bmask = 0x0000ff00;
  amask = 0x000000ff;
#else
  rmask = 0x000000ff;
  gmask = 0x0000ff00;
  bmask = 0x00ff0000;
  amask = 0xff000000;
#endif

  textSurface = SDL_CreateRGBSurface(0, INFO_WIDTH, INFO_HEIGHT, 32, rmask, gmask, bmask, amask);
  SDL_SetSurfaceBlendMode(textSurface, SDL_BLENDMODE_BLEND);
  cu_info_t* selected_cu = NULL;


  kvz_mutex_unlock(&sdl_mutex);
  int locked = 0;
  for (;;) {
    SDL_Event event;
    for (;;) {

      while (SDL_PollEvent(&event)) {
        if (event.type == SDL_MOUSEMOTION) {
          // Update mouse coordinates
          SDL_GetMouseState(&mouse_x, &mouse_y);
          // If mouse inside the perimeters of the frame
          if (mouse_x > 0 && mouse_x < screen_w && mouse_y > 0 && mouse_y < screen_h) {
            cu_info_t *over_cu = &sdl_cu_array[(mouse_x / SCU_WIDTH) + (mouse_y / SCU_WIDTH)*(width_in_lcu * LCU_CU_WIDTH)];
            char* cu_types[5] = { "64x64", "32x32", "16x16", "8x8", "4x4" };

            // If block has changed
            if (sdl_block_info && selected_cu != over_cu) {
              char temp[128];
              selected_cu = over_cu;
              sprintf(temp, "Block type: Intra\nIntra mode: %d", over_cu->intra.mode);

              // Clear the hilight buffer
              memset(sdl_pixels_hilight, 0, (screen_w*screen_h * 4));
              int depth = over_cu->part_size == SIZE_2Nx2N ? over_cu->depth : 3;
              int tr_depth = over_cu->tr_depth;
              int block_border_x = (mouse_x / (LCU_WIDTH >> (depth))) *(LCU_WIDTH >> (depth));
              int block_border_y = (mouse_y / (LCU_WIDTH >> (depth))) *(LCU_WIDTH >> (depth));

              int block_border_tr_x = (mouse_x / (LCU_WIDTH >> (tr_depth))) *(LCU_WIDTH >> (tr_depth));
              int block_border_tr_y = (mouse_y / (LCU_WIDTH >> (tr_depth))) *(LCU_WIDTH >> (tr_depth));

              for (int y = block_border_y; y < block_border_y + (LCU_WIDTH >> over_cu->depth); y++) {
                for (int x = block_border_x; x < block_border_x + (LCU_WIDTH >> over_cu->depth); x++) {
                  kvz_putpixel(sdl_pixels_hilight, screen_w, x, y, 255, 255, 255, 128);
                }
              }
              if (over_cu->depth != over_cu->tr_depth) {
                for (int y = block_border_tr_y; y < block_border_tr_y + (LCU_WIDTH >> over_cu->tr_depth); y++) {
                  for (int x = block_border_tr_x; x < block_border_tr_x + (LCU_WIDTH >> over_cu->tr_depth); x++) {
                    kvz_putpixel(sdl_pixels_hilight, screen_w, x, y, 100, 100, 255, 128);
                  }
                }
              }
              SDL_FillRect(textSurface, NULL, SDL_MapRGBA(textSurface->format, 0, 0, 0, 0));

              if (over_cu->type == CU_INTRA) {
                sprintf(temp, "Type: Intra");
                sdl_render_multiline_text(temp, 0, 0);
                sprintf(temp, "Size: %s", cu_types[over_cu->depth]);
                sdl_render_multiline_text(temp, 0, 20);
                sprintf(temp, "Tr-size: %s", cu_types[over_cu->tr_depth]);
                sdl_render_multiline_text(temp, 0, 40);
                if (over_cu->part_size == SIZE_2Nx2N) {
                  sprintf(temp, "Intra mode: %d", over_cu->intra.mode);
                  sdl_render_multiline_text(temp, 0, 60);
                } else {
                  assert(over_cu->depth == MAX_DEPTH); // Might need an update for VVC
                  uint8_t modes[4];
                  uint16_t cu_width = CU_WIDTH_FROM_DEPTH(over_cu->depth);
                  uint16_t cu_x = mouse_x >> over_cu->depth << over_cu->depth;
                  uint16_t cu_y = mouse_y >> over_cu->depth << over_cu->depth;
                  for (int pu_i = 0; pu_i < 4; ++pu_i) {
                    uint16_t pu_off_x = PU_GET_X(SIZE_NxN, cu_width, 0, pu_i);
                    uint16_t pu_off_y = PU_GET_Y(SIZE_NxN, cu_width, 0, pu_i);
                    uint16_t pu_x = cu_x + pu_off_x;
                    uint16_t pu_y = cu_y + pu_off_y;
                    cu_info_t *pu = &sdl_cu_array[(pu_x / SCU_WIDTH) + (pu_y / SCU_WIDTH)*(width_in_lcu * LCU_CU_WIDTH)];
                    modes[pu_i] = pu->intra.mode;
                  }

                  sprintf(temp, "Intra mode: %d, %d, %d, %d", modes[0], modes[1]
                    , modes[2], modes[3]);
                  sdl_render_multiline_text(temp, 0, 60);
                }
              }
              if (over_cu->type == CU_INTER) {
                sprintf(temp, "Type: Inter %s", over_cu->skipped ? "SKIP" : over_cu->merged ? "MERGE" : "");
                sdl_render_multiline_text(temp, 0, 0);
                sprintf(temp, "Size: %s", cu_types[over_cu->depth]);
                sdl_render_multiline_text(temp, 0, 20);
                sprintf(temp, "Tr-size: %s", cu_types[over_cu->tr_depth]);
                sdl_render_multiline_text(temp, 0, 40);
                sprintf(temp, "Dir: %d", over_cu->inter.mv_dir);
                sdl_render_multiline_text(temp, 0, 60);
                if (over_cu->inter.mv_dir == 3) {
                  sprintf(temp, "MV[L0]: %.3f, %.3f MV[L1]: %.3f, %.3f",
                    (float)over_cu->inter.mv[0][0] / 4.0, (float)over_cu->inter.mv[0][1] / 4.0,
                    (float)over_cu->inter.mv[1][0] / 4.0, (float)over_cu->inter.mv[1][1] / 4.0);
                } else if (over_cu->inter.mv_dir & 1) {
                  sprintf(temp, "MV[L0]: %.3f, %.3f",
                    (float)over_cu->inter.mv[0][0] / 4.0, (float)over_cu->inter.mv[0][1] / 4.0);
                } else if (over_cu->inter.mv_dir & 2) {
                  sprintf(temp, "MV[L1]: %.3f, %.3f",
                    (float)over_cu->inter.mv[1][0] / 4.0, (float)over_cu->inter.mv[1][1] / 4.0);
                }
                sdl_render_multiline_text(temp, 0, 80);
              }

              if (text) SDL_DestroyTexture(text);
              text = SDL_CreateTextureFromSurface(info_renderer, textSurface);
              SDL_SetTextureBlendMode(text, SDL_BLENDMODE_BLEND);



              SDL_Rect rect;
              rect.w = screen_w; rect.h = screen_h; rect.x = 0; rect.y = 0;
              SDL_UpdateTexture(overlay_hilight, &rect, sdl_pixels_hilight, screen_w * 4);

              sdl_force_redraw(locked);
            }
          }
        }
        if (event.type == SDL_KEYDOWN) {
          if (event.key.keysym.sym == SDLK_RETURN) {
            Uint32 flags = SDL_GetWindowFlags(window) ^ SDL_WINDOW_FULLSCREEN_DESKTOP;
            if (SDL_SetWindowFullscreen(window, flags) < 0) {
              fprintf(stderr, "Toggling fullscreen failed.\n");
            }
          }
          if (event.key.keysym.sym == SDLK_d) {
            sdl_draw_blocks = sdl_draw_blocks ? 0 : 1;
            sdl_force_redraw(locked);
          }
          if (event.key.keysym.sym == SDLK_b) {
            sdl_block_info = sdl_block_info ? 0 : 1;

            if (sdl_block_info) {
              info_window = SDL_CreateWindow(
                "Info",                  // window title
                SDL_WINDOWPOS_UNDEFINED,           // initial x position
                SDL_WINDOWPOS_UNDEFINED,           // initial y position
                INFO_WIDTH, INFO_HEIGHT,
                0
                );
              info_renderer = SDL_CreateRenderer(info_window, -1, 0);
            } else {
              SDL_DestroyRenderer(info_renderer);
              info_renderer = NULL;
              SDL_DestroyWindow(info_window);
              info_window = NULL;
            }
          }
          if (event.key.keysym.sym == SDLK_i) {
            sdl_draw_intra = sdl_draw_intra ? 0 : 1;
            sdl_force_redraw(locked);
          }
          if (event.key.keysym.sym == SDLK_f) sdl_fader = !sdl_fader;
          if (event.key.keysym.sym == SDLK_o) {
            sdl_faded_overlay = !sdl_faded_overlay;
            SDL_SetTextureAlphaMod(overlay_blocks, sdl_faded_overlay ? 150 : 255);
            sdl_force_redraw(locked);
          }

          if (event.key.keysym.sym == SDLK_KP_PLUS || event.key.keysym.sym == SDLK_COMMA) sdl_delay += 1;
          if (event.key.keysym.sym == SDLK_KP_MINUS || event.key.keysym.sym == SDLK_PERIOD) { sdl_delay -= 1; if (sdl_delay < 0) sdl_delay = 0; }
          if (event.key.keysym.sym == SDLK_p) {
            if (locked) {
              locked = 0;
              kvz_mutex_unlock(&sdl_mutex);
            } else {
              locked = 1;
              kvz_mutex_lock(&sdl_mutex);
            }
            sdl_force_redraw(locked);
          }
        }
        if (event.type == SDL_QUIT) {
          SDL_Quit();
          exit(1);
        }
      }

      if (!locked) {
        kvz_mutex_lock(&sdl_mutex);

        if (sdl_fader) {
          for (int i = 0; i < screen_w*screen_h * 4; i += 4) {
            if (sdl_pixels_RGB[i + 3]) {
              sdl_pixels_RGB[i + 3] = MAX(0, sdl_pixels_RGB[i + 3] - 2);
            }
          }
          SDL_Rect rect;
          rect.w = screen_w; rect.h = screen_h; rect.x = 0; rect.y = 0;
          SDL_UpdateTexture(overlay_blocks, &rect, sdl_pixels_RGB, screen_w * 4);
        }

        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, overlay, NULL, NULL);
        if (sdl_draw_blocks)
          SDL_RenderCopy(renderer, overlay_blocks, NULL, NULL);
        if (sdl_draw_intra) {
          SDL_RenderCopy(renderer, overlay_intra, NULL, NULL);
          SDL_RenderCopy(renderer, overlay_inter[0], NULL, NULL);
          SDL_RenderCopy(renderer, overlay_inter[1], NULL, NULL);
        }
        if (sdl_block_info) {
          SDL_RenderCopy(renderer, overlay_hilight, NULL, NULL);
        }
        SDL_RenderPresent(renderer);

        kvz_mutex_unlock(&sdl_mutex);
      } else {
        SDL_RenderPresent(renderer);
      }

      if (sdl_block_info) {
        SDL_Rect renderQuad;
        renderQuad.w = textSurface->w; renderQuad.h = textSurface->h; renderQuad.x = 0; renderQuad.y = 0;
        SDL_RenderClear(info_renderer);
        SDL_RenderCopy(info_renderer, text, NULL, &renderQuad);
        SDL_RenderPresent(info_renderer);
      }

      SDL_Delay(10); // Limit loop CPU usage
    }
  }
}

void kvz_visualization_init(int width, int height)
{
  screen_w = width;
  screen_h = height;
  pthread_t sdl_thread;

  pthread_mutex_init(&sdl_mutex, NULL);

  // Lock for eventloop thread to unlock
  kvz_mutex_lock(&sdl_mutex);

  if (SDL_Init(SDL_INIT_EVERYTHING) < 0) {
    fprintf(stderr, "Couldn't initialize SDL: %s\n", SDL_GetError());
    exit(EXIT_FAILURE);
  }

  if (pthread_create(&sdl_thread, NULL, kvz_visualization_eventloop, NULL) != 0) {
    fprintf(stderr, "pthread_create failed!\n");
    exit(EXIT_FAILURE);
  }

  // Wait for eventloop to handle opening the window etc
  kvz_mutex_lock(&sdl_mutex);
  kvz_mutex_unlock(&sdl_mutex);
}

void kvz_visualization_free()
{
  free(sdl_pixels);
  free(sdl_pixels_RGB);
  SDL_Quit();
}


static void render_image(encoder_control_t *encoder, kvz_picture *image)
{
  kvz_mutex_lock(&sdl_mutex);
  
  memcpy(sdl_pixels, image->y, (encoder->cfg.width * encoder->cfg.height));
  memcpy(sdl_pixels_u, image->u, (encoder->cfg.width * encoder->cfg.height) >> 2);
  memcpy(sdl_pixels_v, image->v, (encoder->cfg.width * encoder->cfg.height) >> 2);

  SDL_Rect rect;
  rect.w = screen_w; rect.h = screen_h; rect.x = 0; rect.y = 0;
  SDL_UpdateYUVTexture(overlay, &rect, sdl_pixels, encoder->cfg.width, sdl_pixels_u, encoder->cfg.width >> 1, sdl_pixels_v, encoder->cfg.width >> 1);
  SDL_RenderClear(renderer);
  SDL_RenderCopy(renderer, overlay, NULL, NULL);
  SDL_RenderPresent(renderer);

  kvz_mutex_unlock(&sdl_mutex);
}


void kvz_visualization_frame_init(encoder_control_t *encoder, kvz_picture *target_img)
{
  // This function is called for every frame from the main thread.
  static bool is_first_frame = true;
  if (is_first_frame) {
    render_image(encoder, target_img);
    is_first_frame = false;
  }
}

/**
 * \brief Draw block borders.
 * \return true if block was drawn, false otherwise
 */
bool kvz_visualization_draw_block(const encoder_state_t *state, lcu_t *lcu, cu_info_t *cur_cu, int x, int y, int depth)
{
  const int cu_width = LCU_WIDTH >> depth;

  if (!(x + cu_width <= state->tile->frame->source->width && y + cu_width <= state->tile->frame->source->height)) {
    return false;
  }
  
  SDL_Rect rect;

  const encoder_control_t *ctrl = state->encoder_control;
  kvz_picture * const pic = state->tile->frame->source;

  const vector2d_t frame_px = {
    x + state->tile->lcu_offset_x * LCU_WIDTH,
    y + state->tile->lcu_offset_y * LCU_WIDTH
  };
  const vector2d_t lcu_px = { x & (LCU_WIDTH - 1), y & (LCU_WIDTH - 1) };
  const vector2d_t lcu_px_c = { lcu_px.x / 2, lcu_px.y / 2 };

  const int pic_width = state->encoder_control->cfg.width;
  const int visible_width = MIN(frame_px.x + cu_width, pic->width) - frame_px.x;
  const int visible_height = MIN(frame_px.y + cu_width, pic->height) - frame_px.y;
  const int index_RGB = (x + y * pic_width +
    state->tile->lcu_offset_x*LCU_WIDTH +
    state->tile->lcu_offset_y *LCU_WIDTH * pic_width)<<2;
  const int luma_index = x + y * pic_width +
    state->tile->lcu_offset_x*LCU_WIDTH +
    state->tile->lcu_offset_y *LCU_WIDTH * pic_width;
  const int chroma_index = (x / 2) + (y / 2) * (pic_width / 2) +
    state->tile->lcu_offset_x*(LCU_WIDTH / 2) +
    state->tile->lcu_offset_y *(LCU_WIDTH / 2) * (pic_width / 2);

  if ((cur_cu->depth == 0) || cur_cu->depth == depth || !(depth < ctrl->cfg.pu_depth_intra.max || depth < ctrl->cfg.pu_depth_inter.max)) {
    for (int row = 0; row < visible_height; ++row) {
      memcpy(&sdl_pixels[frame_px.x + (frame_px.y + row) * pic_width],
             &lcu->rec.y[lcu_px.x + (lcu_px.y + row) * LCU_WIDTH],
             visible_width * sizeof(kvz_pixel));
    }
    for (int row = 0; row < visible_height / 2; ++row) {
      memcpy(&sdl_pixels_u[frame_px.x / 2 + (frame_px.y / 2 + row) * pic_width / 2],
             &lcu->rec.u[lcu_px.x / 2 + (lcu_px.y / 2 + row) * LCU_WIDTH_C],
             visible_width / 2 * sizeof(kvz_pixel));
    }
    for (int row = 0; row < visible_height / 2; ++row) {
      memcpy(&sdl_pixels_v[frame_px.x / 2 + (frame_px.y / 2 + row) * pic_width / 2],
             &lcu->rec.v[lcu_px.x / 2 + (lcu_px.y / 2 + row) * LCU_WIDTH_C],
             visible_width / 2 * sizeof(kvz_pixel));
    }

    // Clear RGB buffer area
    {
      int temp_y;
      for (temp_y = 0; temp_y < cu_width; temp_y++) {
        memset(&sdl_pixels_RGB[index_RGB + (temp_y*pic_width << 2)], 0, cu_width << 2);
        memset(&sdl_pixels_RGB_intra_dir[index_RGB + (temp_y*pic_width << 2)], 0, cu_width << 2);          
      }
    }


    {
      const int width_cu = cur_cu->part_size == SIZE_2Nx2N ? LCU_CU_WIDTH >> cur_cu->depth : 1;
      const int x_cu = (x / SCU_WIDTH);
      const int y_cu = (y / SCU_WIDTH);
      int temp_x, temp_y;
      // Set mode in every CU covered by part_mode in this depth.
      for (temp_y = y_cu; temp_y < y_cu + width_cu; ++temp_y) {
        for (temp_x = x_cu; temp_x < x_cu + width_cu; ++temp_x) {
          cu_info_t *cu = &sdl_cu_array[temp_x + temp_y * (state->tile->frame->width_in_lcu * LCU_CU_WIDTH)];
          memcpy(cu, cur_cu, sizeof(cu_info_t));
        }
      }
    }

    //if (cu_width > 4 || (!(x & 7) && !(y & 7))) 

    uint8_t framemod = state->frame->num % 8;
    
    {
      int temp_x;
      // Add block borders
      if ((y + cu_width) % 8 == 0) {
        for (temp_x = 0; temp_x < cu_width; temp_x++) {
          kvz_putpixel(sdl_pixels_RGB + index_RGB, pic_width, temp_x, cu_width - 1, frame_r[framemod], frame_g[framemod], frame_b[framemod], 255);
        }
      }
      if ((x + cu_width) % 8 == 0) {
        int y;
        for (y = 0; y < cu_width; y++) {
          kvz_putpixel(sdl_pixels_RGB + index_RGB, pic_width, cu_width - 1, y, frame_r[framemod], frame_g[framemod], frame_b[framemod], 255);
        }
      }
    }

    // Intra directions
    if (cur_cu->type == CU_INTRA) {
      int mode = cur_cu->intra.mode;
      // These define end points for lines originating from the center of the block.
      const int line_from_center_8x8_x[] = { 8, 8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -7, -6, -5, -4, -3, -2, -1,  0,  1,  2,  3,  4,  5,  6,  7,  8};
      const int line_from_center_8x8_y[] = {-8,-8,  8,  7,  6,  5,  4,  3,  2,  1,  0, -1, -2, -3, -4, -5, -6, -7, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8};
      int mode_line_x = line_from_center_8x8_x[mode];
      int mode_line_y = line_from_center_8x8_y[mode];
      int center = cu_width / 2 - 1;
      
      if (mode == 0) { // Planar
        int viz_width = cu_width == 4 ? cu_width / 4 : cu_width / 8;
        
        for (int i = -viz_width; i < viz_width + 1; i++) {
          int xx = center + i * mode_line_x / 8;
          int yy = center + i * mode_line_y / 8;
          kvz_putpixel(sdl_pixels_RGB + index_RGB, pic_width, xx, yy, 255, 255, 255, 255);
        }
        for (int i = -viz_width; i < 0; i++) {
          int xx = center + i * mode_line_x / 8;
          int yy = center - i * mode_line_y / 8;
          kvz_putpixel(sdl_pixels_RGB + index_RGB, pic_width, xx, yy, 255, 255, 255, 255);
        }
      } else if (mode == 1) { // DC
        int viz_width = cu_width == 4 ? cu_width / 4 : cu_width / 8;
        for (int i = -viz_width; i < viz_width + 1; i++) {
          int xx = center + i * mode_line_x / 8;
          int yy = center + i * mode_line_y / 8;
          kvz_putpixel(sdl_pixels_RGB + index_RGB, pic_width, xx, yy, 255, 255, 255, 255);
        }
        for (int i = -viz_width; i < viz_width + 1; i++) {
          int xx = center + i * mode_line_x / 8;
          int yy = center - i * mode_line_y / 8;
          kvz_putpixel(sdl_pixels_RGB + index_RGB, pic_width, xx, yy, 255, 255, 255, 255);
        }
      } else { // Angular
        for (int i = -cu_width / 4; i < cu_width / 4 + 1; i++) {
          int xx = center + i * mode_line_x / 8;
          int yy = center + i * mode_line_y / 8;
          kvz_putpixel(sdl_pixels_RGB_intra_dir + index_RGB, pic_width, xx, yy, 255, 255, 255, 255);
        }
      }
    }
  }

  rect.w = cu_width; rect.h = cu_width; rect.x = x + state->tile->lcu_offset_x*LCU_WIDTH; rect.y = y + state->tile->lcu_offset_y*LCU_WIDTH;
  SDL_UpdateYUVTexture(overlay, &rect, sdl_pixels + luma_index, pic_width, sdl_pixels_u + chroma_index, pic_width >> 1, sdl_pixels_v + chroma_index, pic_width >> 1);
  SDL_UpdateTexture(overlay_blocks, &rect, sdl_pixels_RGB+index_RGB, pic_width * 4);
  SDL_UpdateTexture(overlay_intra, &rect, sdl_pixels_RGB_intra_dir + index_RGB, pic_width * 4);

  return true;
}


bool kvz_visualization_draw_block_with_delay(const encoder_state_t *state, lcu_t *lcu, cu_info_t *cur_cu, int x, int y, int depth)
{
  const int cu_width = LCU_WIDTH >> depth;
  if (!(x + cu_width <= state->tile->frame->source->width &&
        y + cu_width <= state->tile->frame->source->height))
  {
    return false;
  }

  kvz_mutex_lock(&sdl_mutex);

  // If drawing and delay are not within a single mutex lock, it will look choppy.
  bool block_drawn = kvz_visualization_draw_block(state, lcu, cur_cu, x, y, depth);
  if (block_drawn) {
    kvz_visualization_delay();
  }

  kvz_mutex_unlock(&sdl_mutex);
  return block_drawn;
}


static void vis_mv(encoder_state_t * const state, int x, int y, lcu_t *lcu, int depth, int mv_dir)
{
  const int cu_width = LCU_WIDTH >> depth;
  const cu_info_t *cur_cu = LCU_GET_CU_AT_PX(lcu, SUB_SCU(x), SUB_SCU(y));
  const int poc = state->frame->poc;
  
  // The right and bottom borders are located within the block,
  // so the center of the block is (width / 2 - 1).
  const int x1 = x + cu_width / 2 - 1 + state->tile->lcu_offset_x * LCU_WIDTH;
  const int y1 = y + cu_width / 2 - 1 + state->tile->lcu_offset_y * LCU_WIDTH;
  const int x2 = CLIP(0, screen_w - 1, x1 + ((cur_cu->inter.mv[mv_dir][0] + 2) >> 2));
  const int y2 = CLIP(0, screen_h - 1, y1 + ((cur_cu->inter.mv[mv_dir][1] + 2) >> 2));

  const int ref_idx = MIN(2, cur_cu->inter.mv_ref[mv_dir]);
  const int ref_poc = state->frame->ref->pocs[ref_idx];
  const int ref_framemod = ref_poc % 8;

  if (x1 != x2 || y1 != y2) {
    draw_mv(sdl_pixels_RGB_inter[poc % 2], screen_w, screen_h, x1, y1, x2, y2, frame_r[ref_framemod], frame_g[ref_framemod], frame_b[ref_framemod]);
  }

  vector2d_t tl = { screen_w - 1, screen_h - 1 };
  vector2d_t br = { 0, 0 };
  tl.x = MIN(MIN(tl.x, x1), x2);
  tl.y = MIN(MIN(tl.y, y1), y2);
  br.x = MAX(MAX(br.x, x1), x2);
  br.y = MAX(MAX(br.y, y1), y2);

  tl.x = CLIP(0, screen_w - 1, tl.x);
  tl.y = CLIP(0, screen_h - 1, tl.y);
  br.x = CLIP(0, screen_w - 2, br.x);
  br.y = CLIP(0, screen_h - 2, br.y);

  SDL_Rect rect = {
    tl.x,
    tl.y,
    br.x - tl.x + 1,
    br.y - tl.y + 1,
  };
  
  if (rect.w > 1 || rect.h > 1) {
    SDL_UpdateTexture(overlay_inter[poc % 2], &rect, sdl_pixels_RGB_inter[poc % 2] + (tl.x + tl.y * screen_w) * 4, screen_h * 4);
  }
}

static void recur_vis_mv(encoder_state_t * const state, int x0, int y0, lcu_t *lcu, int depth)
{
  int cu_width = LCU_WIDTH >> depth;
  cu_info_t *cur_cu = LCU_GET_CU_AT_PX(lcu, SUB_SCU(x0), SUB_SCU(y0));

  if (cur_cu->type != CU_INTER) return;

  const int x1 = x0 + cu_width / 2;
  const int y1 = y0 + cu_width / 2;
  if (cur_cu->depth > depth) {
    recur_vis_mv(state, x0, y0, lcu, depth + 1);
    recur_vis_mv(state, x1, y0, lcu, depth + 1);
    recur_vis_mv(state, x0, y1, lcu, depth + 1);
    recur_vis_mv(state, x1, y1, lcu, depth + 1);
    return;
  } else {
    if (cur_cu->inter.mv_dir & 1) {
      vis_mv(state, x0, y0, lcu, depth, 0);
    }
    if (cur_cu->inter.mv_dir & 2) {
      vis_mv(state, x0, y0, lcu, depth, 1);
    }
  }
}

void kvz_visualization_mv_draw_lcu(encoder_state_t * const state, int x, int y, lcu_t *lcu)
{
  kvz_mutex_lock(&sdl_mutex);
  recur_vis_mv(state, x, y, lcu, 0);
  kvz_mutex_unlock(&sdl_mutex);
}

void kvz_visualization_mv_clear_lcu(encoder_state_t * const state, int x, int y)
{
  kvz_mutex_lock(&sdl_mutex);

  const int poc = state->frame->poc;
  const int lcu_width = (x + 64 >= screen_w ? screen_w - x : 64);
  const int lcu_height = (y + 64 >= screen_h ? screen_h - y : 64);

  kvz_pixel *buffer = sdl_pixels_RGB_inter[(poc + 1) % 2];
  for (int lcu_y = 0; lcu_y < lcu_height; ++lcu_y) {
    int index = (x + (y + lcu_y) * screen_w) * 4;
    memset(&buffer[index], 0, lcu_width * 4);
  }

  SDL_Rect lcu_rect = {
    x, y,
    lcu_width, lcu_height
  };
  SDL_UpdateTexture(overlay_inter[(poc + 1) % 2], &lcu_rect, buffer + (x + y * screen_w) * 4, screen_w * 4);

  kvz_mutex_unlock(&sdl_mutex);
}

#endif
