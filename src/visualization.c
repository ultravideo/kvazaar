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

#if KVZ_VISUALIZATION == 1

#include "encoderstate.h"

SDL_Renderer *renderer, *info_renderer;
SDL_Window *window, *info_window = NULL;
SDL_Surface *screen, *pic;
SDL_Texture *overlay, *overlay_blocks, *overlay_intra, *overlay_inter[2], *overlay_hilight;
int screen_w, screen_h;
int sdl_draw_blocks = 1;
int sdl_draw_intra = 1;
int sdl_block_info = 0;
pthread_mutex_t sdl_mutex;
kvz_pixel *sdl_pixels_hilight;
kvz_pixel *sdl_pixels_RGB;
kvz_pixel *sdl_pixels_RGB_intra_dir;
kvz_pixel *sdl_pixels_RGB_inter[2];
kvz_pixel *sdl_pixels;
kvz_pixel *sdl_pixels_u;
kvz_pixel *sdl_pixels_v;
int32_t sdl_delay;
SDL_Surface* textSurface;
SDL_Texture* text;

cu_info_t *sdl_cu_array;
TTF_Font* font;


kvz_visualization_init(int width, int height)
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

  if (pthread_create(&sdl_thread, NULL, eventloop_main, NULL) != 0) {
    fprintf(stderr, "pthread_create failed!\n");
    exit(EXIT_FAILURE);
  }
}

kvz_visualization_free()
{
  free(sdl_pixels);
  free(sdl_pixels_RGB);
  SDL_Quit();
}

void kvz_visualization_frame_init(encoder_control_t *encoder, kvz_picture *img_in)
{
  static bool allocated = false;
  if (allocated) {
    return;
  } else {
    allocated = true;
  }

  kvz_mutex_lock(&sdl_mutex);
  // Copy original frame with darkened colors
  for (int y = 0; y < encoder->cfg->height; y++) {
    for (int x = 0; x < encoder->cfg->width; x++) {
      int16_t pix_value = img_in->y[x + y*encoder->cfg->width] - 10;
      if (pix_value < 0) pix_value = 0;
      sdl_pixels[x + y*encoder->cfg->width] = sdl_pixels[x + y*encoder->cfg->width] = pix_value;
    }
  }

  // Copy chroma to both overlays
  memcpy(sdl_pixels_u, img_in->u, (encoder->cfg->width*encoder->cfg->height) >> 2);
  memcpy(sdl_pixels_v, img_in->v, (encoder->cfg->width*encoder->cfg->height) >> 2);

  // ToDo: block overlay
  //memcpy(overlay_blocks->pixels[1], img_in->u, (encoder->cfg->width*encoder->cfg->height) >> 2);
  //memcpy(overlay_blocks->pixels[2], img_in->v, (encoder->cfg->width*encoder->cfg->height) >> 2);

  SDL_Rect rect;
  rect.w = screen_w; rect.h = screen_h; rect.x = 0; rect.y = 0;
  SDL_UpdateYUVTexture(overlay, &rect, sdl_pixels, encoder->cfg->width, sdl_pixels_u, encoder->cfg->width >> 1, sdl_pixels_v, encoder->cfg->width >> 1);
  SDL_RenderClear(renderer);
  SDL_RenderCopy(renderer, overlay, NULL, NULL);
  SDL_RenderPresent(renderer);
  kvz_mutex_unlock(&sdl_mutex);
}


#define PUTPIXEL_hilight(pixel_x, pixel_y, color_r, color_g, color_b, color_alpha) sdl_pixels_hilight[(pixel_x<<2) + (pixel_y)*(screen_w<<2)+3] = color_alpha; \
  sdl_pixels_hilight[(pixel_x<<2) + (pixel_y)*(screen_w<<2) +2] = color_r; \
  sdl_pixels_hilight[(pixel_x<<2) + (pixel_y)*(screen_w<<2) +1] = color_g; \
  sdl_pixels_hilight[(pixel_x<<2) + (pixel_y)*(screen_w<<2) +0] = color_b;

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

void *eventloop_main(void* temp)
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

  sdl_cu_array = MALLOC(cu_info_t, (height_in_lcu << MAX_DEPTH)*(width_in_lcu << MAX_DEPTH));
  FILL_ARRAY(sdl_cu_array, 0, (height_in_lcu << MAX_DEPTH)*(width_in_lcu << MAX_DEPTH));

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
  SDL_Color White = { 255, 255, 255 };

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
    while (1) {

      while (SDL_PollEvent(&event)) {
        if (event.type == SDL_MOUSEMOTION) {
          // Update mouse coordinates
          SDL_GetMouseState(&mouse_x, &mouse_y);
          // If mouse inside the perimeters of the frame
          if (mouse_x > 0 && mouse_x < screen_w && mouse_y > 0 && mouse_y < screen_h) {

            cu_info_t *over_cu = &sdl_cu_array[(mouse_x / (LCU_WIDTH >> MAX_DEPTH)) + (mouse_y / (LCU_WIDTH >> MAX_DEPTH))*(width_in_lcu << MAX_DEPTH)];
            char* cu_types[5] = { "64x64", "32x32", "16x16", "8x8", "4x4" };

            // If block has changed
            if (sdl_block_info && selected_cu != over_cu) {
              char temp[128];
              selected_cu = over_cu;
              sprintf(temp, "Block type: Intra\nIntra mode: %d", over_cu->intra->mode);

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
                  PUTPIXEL_hilight(x, y, 255, 255, 255, 128);
                }
              }
              if (over_cu->depth != over_cu->tr_depth) {
                for (int y = block_border_tr_y; y < block_border_tr_y + (LCU_WIDTH >> over_cu->tr_depth); y++) {
                  for (int x = block_border_tr_x; x < block_border_tr_x + (LCU_WIDTH >> over_cu->tr_depth); x++) {
                    PUTPIXEL_hilight(x, y, 100, 100, 255, 128);
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
                  sprintf(temp, "Intra mode: %d", over_cu->intra->mode);
                  sdl_render_multiline_text(temp, 0, 60);
                } else {
                  sprintf(temp, "Intra mode: %d, %d, %d, %d", over_cu->intra[0].mode, over_cu->intra[1].mode
                    , over_cu->intra[2].mode, over_cu->intra[3].mode);
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

          if (event.key.keysym.sym == SDLK_KP_PLUS) sdl_delay += 1;
          if (event.key.keysym.sym == SDLK_KP_MINUS) { sdl_delay -= 1; if (sdl_delay < 0) sdl_delay = 0; }
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

void kvz_visualization_draw_block(const encoder_state_t *state, lcu_t *lcu, cu_info_t *cur_cu, int x, int y, int depth)
{
  const int cu_width = LCU_WIDTH >> depth;

  if (!(x + cu_width <= state->tile->frame->source->width && y + cu_width <= state->tile->frame->source->height)) return;
  
  SDL_Rect rect;

  encoder_control_t *ctrl = state->encoder_control;
  kvz_picture * const pic = state->tile->frame->source;

  
  const int pic_width = state->encoder_control->cfg->width;
  const int pic_height = state->encoder_control->cfg->height;
  const int x_max = MIN(x + cu_width, pic->width) - x;
  const int y_max = MIN(y + cu_width, pic->height) - y;
  const int index_RGB = (x + y * pic_width +
    state->tile->lcu_offset_x*LCU_WIDTH +
    state->tile->lcu_offset_y *LCU_WIDTH * pic_width)<<2;
  const int luma_index = x + y * pic_width +
    state->tile->lcu_offset_x*LCU_WIDTH +
    state->tile->lcu_offset_y *LCU_WIDTH * pic_width;
  const int chroma_index = (x / 2) + (y / 2) * (pic_width / 2) +
    state->tile->lcu_offset_x*(LCU_WIDTH / 2) +
    state->tile->lcu_offset_y *(LCU_WIDTH / 2) * (pic_width / 2);

  if ((cur_cu->depth == 0) || cur_cu->depth == depth || !(depth < ctrl->pu_depth_intra.max || depth < ctrl->pu_depth_inter.max)) {
    kvz_pixels_blit(&lcu->rec.y[(x & 63) + (y & 63)*LCU_WIDTH], &sdl_pixels[luma_index],
      x_max, y_max, LCU_WIDTH, pic_width);
    kvz_pixels_blit(&lcu->rec.u[(x & 63) / 2 + (y & 63)*LCU_WIDTH / 4], &sdl_pixels_u[chroma_index],
      x_max / 2, y_max / 2, LCU_WIDTH / 2, pic_width / 2);
    kvz_pixels_blit(&lcu->rec.v[(x & 63) / 2 + (y & 63)*LCU_WIDTH / 4], &sdl_pixels_v[chroma_index],
      x_max / 2, y_max / 2, LCU_WIDTH / 2, pic_width / 2);

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
      const int x_cu = (x / (LCU_WIDTH >> MAX_DEPTH));
      const int y_cu = (y / (LCU_WIDTH >> MAX_DEPTH));
      int temp_x, temp_y;
      // Set mode in every CU covered by part_mode in this depth.
      for (temp_y = y_cu; temp_y < y_cu + width_cu; ++temp_y) {
        for (temp_x = x_cu; temp_x < x_cu + width_cu; ++temp_x) {
          cu_info_t *cu = &sdl_cu_array[temp_x + temp_y *  (state->tile->frame->width_in_lcu << MAX_DEPTH)];
          memcpy(cu, cur_cu, sizeof(cu_info_t));
        }
      }
    }

    //if (cu_width > 4 || (!(x & 7) && !(y & 7))) 

    uint8_t framemod = state->global->frame % 8;
      
    {
      int temp_x;
      // Add block borders
      if ((y + cu_width) % 8 == 0) {
        for (temp_x = 0; temp_x < cu_width; temp_x++) {
          PUTPIXEL(temp_x, (cu_width - 1), frame_r[framemod], frame_g[framemod], frame_b[framemod], 255);
        }
      }
      if ((x + cu_width) % 8 == 0) {
        int y;
        for (y = 0; y < cu_width; y++) {
          PUTPIXEL((cu_width - 1), y, frame_r[framemod], frame_g[framemod], frame_b[framemod], 255);
        }
      }
    }

    // Intra directions
    if (cur_cu->type == CU_INTRA) {
      int i = 1;
      int mode = cur_cu->intra[PU_INDEX(x / 4, y / 4)].mode;
      const int x_off[] = { 8, 8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -7, -6, -5, -4, -3, -2, -1,  0,  1,  2,  3,  4,  5,  6,  7,  8};
      const int y_off[] = {-8,-8,  8,  7,  6,  5,  4,  3,  2,  1,  0, -1, -2, -3, -4, -5, -6, -7, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8};        
      if (mode == 0) { // Planar
        int viz_width = cu_width == 4 ? cu_width / 4 : cu_width / 8;
        for (i = -viz_width; i < viz_width + 1; i++) {
          PUTPIXEL(((cu_width >> 1) + ((i*x_off[mode]) >> 3) - 1), ((cu_width >> 1) + ((i*y_off[mode]) >> 3) - 1), 255, 255, 255, 255);
        }
        for (i = -viz_width; i < 0; i++) {
          PUTPIXEL(((cu_width >> 1) + ((i*x_off[mode]) >> 3) - 1), ((cu_width >> 1) + ((-i*y_off[mode]) >> 3) - 1), 255, 255, 255, 255);
        }
      } else if (mode == 1) { // DC
        int viz_width = cu_width == 4 ? cu_width / 4 : cu_width / 8;
        for (i = -viz_width; i < viz_width + 1; i++) {
          PUTPIXEL(((cu_width >> 1) + ((i*x_off[mode]) >> 3) - 1), ((cu_width >> 1) + ((i*y_off[mode]) >> 3) - 1), 255, 255, 255, 255);
        }
        for (i = -viz_width; i < viz_width + 1; i++) {
          PUTPIXEL(((cu_width >> 1) + ((i*x_off[mode]) >> 3) - 1), ((cu_width >> 1) + ((-i*y_off[mode]) >> 3) - 1), 255, 255, 255, 255);
        }
      } else { // Angular
        for (i = -cu_width / 4; i < cu_width / 4 + 1; i++) {
          PUTPIXEL_intra(((cu_width >> 1) + ((i*x_off[mode]) >> 3) - 1), ((cu_width >> 1) + ((i*y_off[mode]) >> 3) - 1), 255, 255, 255, 255);
        }
      }
    }
  }
  rect.w = cu_width; rect.h = cu_width; rect.x = x + state->tile->lcu_offset_x*LCU_WIDTH; rect.y = y + state->tile->lcu_offset_y*LCU_WIDTH;
  SDL_UpdateYUVTexture(overlay, &rect, sdl_pixels + luma_index, pic_width, sdl_pixels_u + chroma_index, pic_width >> 1, sdl_pixels_v + chroma_index, pic_width >> 1);
  SDL_UpdateTexture(overlay_blocks, &rect, sdl_pixels_RGB+index_RGB, pic_width * 4);
  SDL_UpdateTexture(overlay_intra, &rect, sdl_pixels_RGB_intra_dir + index_RGB, pic_width * 4);
}

#endif
