//Main function for testing the scaler class

#include <stdlib.h>
#include <stdio.h>

#include "scaler.h"

#define in_y_width 8
#define in_y_height 4
#define in_cb_width 4
#define in_cb_height 2

#define out_y_width 4
#define out_y_height 2

void printPicBuffer(pic_buffer_t* buffer)
{
  for (int i = 0; i < buffer->height; i++) {
    for (int j = 0; j < buffer->width; j++) {
      printf("%i ", buffer->data[buffer->width*i + j]);
    }
    printf("\n");
  }
}

void printout( yuv_buffer_t* buffer )
{
  printf("Y:\n");
  printPicBuffer(buffer->y);
  printf("Cb:\n");
  printPicBuffer(buffer->u);
  printf("Cr:\n");
  printPicBuffer(buffer->v);
}

int main()
{
  //Create a simple "picture" to debug scaler

  uint8_t y_data[in_y_width*in_y_height] = {
       25, 64, 240, 40 , 16, 250, 0, 42,
       125, 164, 240, 140 , 16, 250, 3, 12,
       /*25, 164, 20, 40 , 16, 250, 0, 4,
       25, 14, 140, 50 , 16, 205, 234, 44,
       57, 82, 34, 90, 65, 90, 44, 22,
       89, 92, 0, 71, 61, 78, 109, 100,*/
       0, 124, 78, 56, 29, 0, 4, 8,
       4, 7, 56, 12, 49, 7, 2, 0
  };

  uint8_t cb_data[in_cb_width*in_cb_height] = {
       240, 40 , 16, 7,
       /*164, 240, 16, 7,
       16, 16, 16, 7,*/
       7, 35, 79, 5
  };

  uint8_t cr_data[in_cb_width*in_cb_height] = {
       40, 140 , 16, 7,
       /*135, 40 , 16, 6,
       16, 16, 16, 5,*/
       8, 54, 21, 4
  };


  int is_420 = in_y_width != in_cb_width ? 1 : 0;
  yuv_buffer_t* pic = newYuvBuffer_uint8(y_data, cb_data, cr_data, in_y_width, in_y_height, is_420);
  scaling_parameter_t param = newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height);

  yuv_buffer_t* scaled = yuvDownscaling(pic, &param, is_420);
  printout(scaled);

  //Free memory
  deallocateYuvBuffer(pic);
  deallocateYuvBuffer(scaled);

  system("Pause");
  return EXIT_SUCCESS;
}