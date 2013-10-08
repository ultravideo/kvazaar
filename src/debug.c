/**
 * \file
 * 
 * \author Marko Viitanen ( fador@iki.fi ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 * \author Ari Koivula ( ari@koivu.la ), 
 *         Tampere University of Technology,
 *         Department of Pervasive Computing.
 */

#include "debug.h"

#include <stdio.h>


/**
 * Return a FILE* to a file that can be used with render_cu_file and close_cu_file.
 */
FILE * open_cu_file(char *filename) {
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "<?xml version='1.0' encoding='UTF-8' ?>\r\n"
          "<html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en'>\r\n"
          "<head><link rel='stylesheet' type='text/css' href='cu_style.css' /></head><body>");
  return fp;
}

/**
 * Close the FILE* returned by open_cu_file.
 */
void close_cu_file(FILE *fp) {
  fprintf(fp, "</body></html>");
  fclose(fp);
}

void yuv2rgb(unsigned char yuv[3], unsigned char rgb[3])
{
  int y = yuv[0];
  int u = yuv[1];
  int v = yuv[2];

  int r = 1.164 * y + 1.596 * (v - 128);
  int g = 1.165 * y - 0.392 * (u - 128) - 0.813 * (v - 128);
  int b = 1.164 * y + 2.017 * (u - 128);

  rgb[0] = CLIP(0, 255, r);
  rgb[1] = CLIP(0, 255, g);
  rgb[2] = CLIP(0, 255, b);
}

/**
 * Print information about the Coding Unit (CU) into the FILE* provided by open_cu_file.
 */
unsigned render_cu_file(encoder_control *encoder, picture *pic, 
                        unsigned depth, uint16_t xCtb, uint16_t yCtb, FILE *fp)
{
  cu_info *cu = &pic->cu_array[depth][xCtb + yCtb * (pic->width_in_lcu<<MAX_DEPTH)];
  cu_info *final_cu = &pic->cu_array[MAX_DEPTH][xCtb + yCtb * (pic->width_in_lcu<<MAX_DEPTH)];
  unsigned lambda_cost = (4 * g_lambda_cost[encoder->QP]) << 4;
  unsigned sum = 0;
  unsigned best_cost = -1;
  char type = cu->type == CU_INTRA ? 'I' : 'P';
  unsigned x = xCtb * CU_MIN_SIZE_PIXELS;
  unsigned y = yCtb * CU_MIN_SIZE_PIXELS;
  unsigned luma = y * pic->width + x;
  unsigned chroma = (y >> 1) * (pic->width >> 1) + (x >> 1);
  unsigned char yuv[3] = { 0, 0, 0 };
  unsigned char rgb[3] = { 0, 0, 0 };

  if (x >= pic->width || y >= pic->height) {
    // Don't output anything for CU's completely outside the botders.
    return 0;
  }

  if (encoder->ref->used_size > 0) {
    const picture *ref_pic = encoder->ref->pics[0];
    yuv[0] = ref_pic->y_recdata[luma];
    yuv[1] = ref_pic->u_recdata[chroma];
    yuv[2] = ref_pic->v_recdata[chroma];
    yuv2rgb(yuv, rgb);
  }

  // Enclose everything in a table with the assumption that this function is
  // called from left to right and from top to down.
  if (depth == 0) {
    if (yCtb == 0 && xCtb == 0) {
      fprintf(fp, "<table><tr><td>");
    } else if (xCtb == 0) {
      fprintf(fp, "</td></tr><tr><td>");
    } else if (xCtb == NO_SCU_IN_LCU(pic->width_in_lcu)
               && yCtb == NO_SCU_IN_LCU(pic->height_in_lcu)) {
      fprintf(fp, "</td></tr></table>");
    } else {
      fprintf(fp, "</td><td>");
    }
  }

  fprintf(fp, 
    "\n<table class='d%u' bgcolor='#%02x%02x%02x'><tr><td colspan='2'>"
    "%u (%u, %u), %c, "
    "c=%u, mv=(%d, %d)</td></tr>\n", 
    depth, rgb[0], rgb[1], rgb[2],
    depth, xCtb, yCtb, (cu->type == CU_INTRA ? 'I' : 'P'),
    cu->inter.cost, cu->inter.mv[0], cu->inter.mv[1]);


  if(depth != MAX_INTER_SEARCH_DEPTH)
  {
    /* Split blocks and remember to change x and y block positions */
    uint8_t change = 1<<(MAX_DEPTH-1-depth);

    fprintf(fp, "<tr><td>");
    sum += render_cu_file(encoder, pic, depth + 1, xCtb, yCtb, fp);
    fprintf(fp, "</td><td>");
    sum += render_cu_file(encoder, pic, depth + 1, xCtb + change, yCtb, fp);
    fprintf(fp, "</td></tr>");

    fprintf(fp, "<tr><td>");
    sum += render_cu_file(encoder, pic, depth + 1, xCtb, yCtb + change, fp);
    fprintf(fp, "</td><td>");
    sum += render_cu_file(encoder, pic, depth + 1, xCtb + change, yCtb + change, fp);
    fprintf(fp, "</td></tr>");

    fprintf(fp, "<tr><td colspan='2'>sum=%u, sum+lambda=%u</td></tr>",
      sum, sum + lambda_cost);
    if (sum + lambda_cost < cu->inter.cost) {
      best_cost = sum + lambda_cost;
    } else {
      best_cost = cu->inter.cost;
    }
  } else {
    best_cost = cu->inter.cost;
  }

  if (depth == 0) {
    fprintf(fp, 
      "<tr><td colspan='2'>"
      "best depth=%u, %c, "
      "c=%u, mv=(%d, %d)</td></tr>\n"
      "</td></tr>", 
      final_cu->depth, (final_cu->type == CU_INTRA ? 'I' : 'P'),
      final_cu->inter.cost, final_cu->inter.mv[0], final_cu->inter.mv[1]);
  }

  fprintf(fp, "</table>");
  return best_cost;
}
