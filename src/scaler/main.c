//Main function for testing the scaler class

#include <stdlib.h>
#include <stdio.h>

#include "scaler.h"
#include <math.h>

#define IN_Y_W 8
#define IN_Y_H 4
#define in_cb_width 4
#define in_cb_height 2

#define OUT_Y_W 4
#define OUT_Y_H 2

#define BUFF_SIZE 256

#ifdef _MSC_VER
#if _MSC_VER <= 1800
#define snprintf(...) _snprintf_s(__VA_ARGS__)
#endif
#endif

void printPicBuffer(pic_buffer_t* buffer)
{
  for (int i = 0; i < buffer->height; i++) {
    for (int j = 0; j < buffer->width; j++) {
      printf("%i ", buffer->data[buffer->width * i + j]);
    }
    printf("\n");
  }
}

void printout(yuv_buffer_t* buffer)
{
  printf("Y:\n");
  printPicBuffer(buffer->y);
  printf("Cb:\n");
  printPicBuffer(buffer->u);
  printf("Cr:\n");
  printPicBuffer(buffer->v);
}

//Copy data to a output array
void copyBack(pic_data_t* dst, uint8_t* src, int size)
{
  for (int i = 0; i < size; i++) {
    dst[i] = (pic_data_t)src[i];
  }
}

void copyFrom(uint8_t* dst, pic_data_t* src, int size)
{
  for (int i = 0; i < size; i++) {
    dst[i] = (uint8_t)src[i];
  }
}

/**
* \brief Read a single frame from a file.
*
* Read luma and chroma values from file. Extend pixels if the image buffer
* is larger than the input image.
*
* \param file          input file
* \param input_width   width of the input video in pixels
* \param input_height  height of the input video in pixels
* \param img_out       image buffer
*
* \return              1 on success, 0 on failure
*/
int yuv_io_read(FILE* file,
                unsigned input_width, unsigned input_height,
                yuv_buffer_t* img_out)
{
  const unsigned y_size = input_width * input_height;
  const unsigned uv_input_width = input_width / 2;
  const unsigned uv_input_height = input_height / 2;
  const unsigned uv_size = uv_input_width * uv_input_height;

  if (input_width == img_out->y->width) {
    // No need to extend pixels.
    const size_t pixel_size = sizeof(unsigned char);
    unsigned char* buffer = malloc(sizeof(unsigned char) * y_size);

    if (fread(buffer, pixel_size, y_size, file) != y_size) {
      free(buffer);
      return 0;
    }
    copyBack(img_out->y->data, buffer, y_size);
    if (fread(buffer, pixel_size, uv_size, file) != uv_size) {
      free(buffer);
      return 0;
    }
    copyBack(img_out->u->data, buffer, uv_size);
    if (fread(buffer, pixel_size, uv_size, file) != uv_size) {
      free(buffer);
      return 0;
    }
    copyBack(img_out->v->data, buffer, uv_size);
  }


  return 1;
}

/**
* \brief Write a single frame to a file.
*
* \param file           output file
* \param img            image to output
* \param output_width   width of the output in pixels
* \param output_height  height of the output in pixels
*
* \return              1 on success, 0 on failure
*/
int yuv_io_write(FILE* file,
                 const yuv_buffer_t* img,
                 unsigned output_width, unsigned output_height)
{
  const int width = img->y->width;
  const unsigned y_size = output_width * output_height;
  const unsigned uv_input_width = output_width / 2;
  const unsigned uv_input_height = output_height / 2;
  const unsigned uv_size = uv_input_width * uv_input_height;
  //printf("chroma size: %i", _msize(img->u->data)/sizeof(pic_data_t));
  unsigned char* buffer = malloc(sizeof(unsigned char) * y_size);
  copyFrom(buffer, img->y->data, y_size);

  for (unsigned int y = 0; y < output_height; ++y) {
    fwrite(&buffer[y * width], sizeof(unsigned char), output_width, file);
    // TODO: Check that fwrite succeeded.
  }
  copyFrom(buffer, img->u->data, uv_size);
  for (unsigned int y = 0; y < output_height / 2; ++y) {
    fwrite(&buffer[y * (width / 2)], sizeof(unsigned char), output_width / 2, file);
  }
  copyFrom(buffer, img->v->data, uv_size);
  for (unsigned int y = 0; y < output_height / 2; ++y) {
    fwrite(&buffer[y * (width / 2)], sizeof(unsigned char), output_width / 2, file);
  }

  return 1;
}

/**
* \brief Downscale given kvz_picture. Use sizes in kvz_picture for scaling
*/
/*void kvzDownscaling(yuv_buffer_t* in, yuv_buffer_t* out)
{
  //Create picture buffers based on given kvz_pictures
  int32_t in_y_width = in->y->width;
  int32_t in_y_height = in->y->height;
  int32_t out_y_width = out->y->width;
  int32_t out_y_height = out->y->height;

  //assumes 420
  int is_420 = in->y->width != in->u->width ? 1 : 0;
  scaling_parameter_t param = kvz_newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, CHROMA_420);

  yuv_buffer_t* scaled = scale(in, &param, is_420);

  //copy data to out
  pic_data_t* tmp = out->y->data;
  out->y->data = scaled->y->data;
  scaled->y->data = tmp;

  tmp = out->u->data;
  out->u->data = scaled->u->data;
  scaled->u->data = tmp;

  tmp = out->v->data;
  out->v->data = scaled->v->data;
  scaled->v->data = tmp;

  //Free memory
  deallocateYuvBuffer(scaled);
}*/

/**
* \brief Use sizes in kvz_picture for scaling. Can handle up and downscaling
*/
void kvzScaling(yuv_buffer_t* in, yuv_buffer_t** out)
{
  //Create picture buffers based on given kvz_pictures
  int32_t in_y_width = in->y->width;
  int32_t in_y_height = in->y->height;
  int32_t out_y_width = (*out)->y->width;
  int32_t out_y_height = (*out)->y->height;

  //assumes 420
  //int is_420 = in->y->width != in->u->width ? 1 : 0;
  scaling_parameter_t param = kvz_newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, CHROMA_420);

  *out = yuvScaling(in, &param, *out);
}

void _kvzScaling(yuv_buffer_t* in, yuv_buffer_t** out)
{
  //Create picture buffers based on given kvz_pictures
  int32_t in_y_width = in->y->width;
  int32_t in_y_height = in->y->height;
  int32_t out_y_width = (*out)->y->width;
  int32_t out_y_height = (*out)->y->height;

  //assumes 420
  //int is_420 = in->y->width != in->u->width ? 1 : 0;
  scaling_parameter_t param = kvz_newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, CHROMA_420);

  *out = kvz_yuvScaling_(in, &param, *out);
}


/**
 * \brief Calculates image PSNR value
 *
 * \param src   source picture
 * \param rec   reconstructed picture
 * \prama psnr  returns the PSNR
 */
static void compute_psnr(const yuv_buffer_t *const src,
                         const yuv_buffer_t *const rec,
                         double psnr[3])
{
  int32_t pixels = src->y->width * src->y->height;
  pic_data_t* src_data[3] = {src->y->data,src->u->data,src->v->data};
  pic_data_t* rec_data[3] = {rec->y->data,rec->u->data,rec->v->data};

  for (int32_t c = 0; c < 3; ++c) {
    int32_t num_pixels = pixels;
    if (c != 0) {
      num_pixels >>= 2;
    }
    psnr[c] = 0;
    for (int32_t i = 0; i < num_pixels; ++i) {
      const int32_t error = src_data[c][i] - rec_data[c][i];
      psnr[c] += error * error;
    }

    // Avoid division by zero
    if (psnr[c] == 0) psnr[c] = 99.0;
    psnr[c] = 10 * log10((num_pixels * 255.0*255.0) / ((double)psnr[c]));;
  }
}

pic_buffer_t* diffComponent(pic_buffer_t* pic1, pic_buffer_t* pic2)
{
  int width = pic1->width;
  int height = pic1->height;
  int row = 0;
  pic_buffer_t* res = newPictureBuffer(width, height, 0);

  for( int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      res->data[row + j] = pic1->data[row + j] - pic2->data[row + j];
    }
    row += width;
  }

  return res;
}

yuv_buffer_t* yuvDiff(yuv_buffer_t* yuv1, yuv_buffer_t* yuv2)
{
  yuv_buffer_t* res = kvz_newYuvBuffer(yuv1->y->width, yuv1->y->height, yuv1->format, 0);

  res->y = diffComponent(yuv1->y, yuv2->y);
  res->u = diffComponent(yuv1->u, yuv2->u);
  res->v = diffComponent(yuv1->v, yuv2->v);

  return res;
}

double meanColorComponent(pic_buffer_t* pic)
{
  int size = pic->width*pic->height;
  long int sum = 0;

  for(int i = 0; i < size; i++) {
    sum += abs(pic->data[i]);
  }
  return (double)sum/(double)size;
}

void printMeanColor(yuv_buffer_t* yuv)
{
  printf("Y avg color: %f\n",meanColorComponent(yuv->y));
  printf("U avg color: %f\n",meanColorComponent(yuv->u));
  printf("V avg color: %f\n",meanColorComponent(yuv->v));
}

double sumComponent(pic_buffer_t* pic)
{
  int size = pic->width*pic->height;
  long int sum = 0;

  for(int i = 0; i < size; i++) {
    sum += abs(pic->data[i]);
  }
  return (double)sum;
}

int isZero(yuv_buffer_t* yuv)
{
  return sumComponent(yuv->y) != 0 ? 0 : (sumComponent(yuv->u) != 0 ? 0 : (sumComponent(yuv->v) != 0 ? 0 : 1));
}

int isSame(yuv_buffer_t* yuv1, yuv_buffer_t* yuv2 )
{
  yuv_buffer_t* result = yuvDiff(yuv1, yuv2);
  
  int res = isZero(result);

  deallocateYuvBuffer(result);

  return res;
}

/*void test1()
{
  //Create a simple "picture" to debug scaler

  uint8_t y_data[IN_Y_W*IN_Y_H] = {
    25, 64, 240, 40, 16, 250, 0, 42,
    125, 164, 240, 140, 16, 250, 3, 12,
    //25, 164, 20, 40 , 16, 250, 0, 4,
    //25, 14, 140, 50 , 16, 205, 234, 44,
    //57, 82, 34, 90, 65, 90, 44, 22,
    //89, 92, 0, 71, 61, 78, 109, 100,
    0, 124, 78, 56, 29, 0, 4, 8,
    4, 7, 56, 12, 49, 7, 2, 0
  };

  uint8_t cb_data[in_cb_width*in_cb_height] = {
    240, 40, 16, 7,
    //164, 240, 16, 7,
    //16, 16, 16, 7,
    7, 35, 79, 5
  };

  uint8_t cr_data[in_cb_width*in_cb_height] = {
    40, 140, 16, 7,
    //135, 40 , 16, 6,
    //16, 16, 16, 5,
    8, 54, 21, 4
  };


  int is_420 = IN_Y_W != in_cb_width ? 1 : 0;
  yuv_buffer_t* pic = newYuvBuffer_uint8(y_data, cb_data, cr_data, IN_Y_W, IN_Y_H, is_420);
  scaling_parameter_t param = kvz_newScalingParameters(IN_Y_W, IN_Y_H, OUT_Y_W, OUT_Y_H, IN_Y_W != in_cb_width ? CHROMA_420 : CHROMA_444);

  yuv_buffer_t* scaled = yuvDownscaling(pic, &param, is_420);
  printout(scaled);

  //Free memory
  deallocateYuvBuffer(pic);
  deallocateYuvBuffer(scaled);
}*/

/*void test2()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width = 601;
  int32_t out_height = 602;
  int framerate = 24;
  
  const char* file_name_format = "Kimono1_%ix%i_%i.yuv";
  
  char in_file_name[BUFF_SIZE];
  sprintf(in_file_name, "Kimono1_%ix%i_%i.yuv", in_width, in_height, framerate);
  
  char out_file_name[BUFF_SIZE];
  sprintf(out_file_name, "Kimono1_%ix%i_%i_s.yuv", out_width, out_height, framerate);

  FILE* file = fopen(in_file_name,"rb");
  if (file == NULL) {
    perror("FIle open failed");
    printf("File name: %s", in_file_name);
  }
  
  yuv_buffer_t* data = newYuvBuffer_uint8(NULL, NULL, NULL, in_width, in_height, 1);
  yuv_buffer_t* out = newYuvBuffer_uint8(NULL, NULL, NULL, out_width, out_height, 1);

  if (yuv_io_read(file, in_width, in_height, data)) {
    
    kvzDownscaling(data, out);
    
    FILE* out_file = fopen(out_file_name, "wb");
    yuv_io_write(out_file, out, out->y->width, out->y->height);
    fclose(out_file);
  }
  else {
    perror("Read fail");
    if (feof(file)) perror("End of file");
    if (ferror(file)) perror("File error");
  }

  deallocateYuvBuffer(data);
  deallocateYuvBuffer(out);
  fclose(file);
}*/

void test3()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width = 100;
  int32_t out_height = 290;
  int framerate = 24;

  //const char* file_name_format = "Kimono1_%ix%i_%i.yuv";

  char in_file_name[BUFF_SIZE];
  sprintf(in_file_name, "Kimono1_%ix%i_%i.yuv", in_width, in_height, framerate);

  char out_file_name[BUFF_SIZE];
  sprintf(out_file_name, "Kimono1_%ix%i_%i_s.yuv", out_width, out_height, framerate);

  FILE* file = fopen(in_file_name, "rb");
  if (file == NULL) {
    perror("FIle open failed");
    printf("File name: %s", in_file_name);
  }

  yuv_buffer_t* data = kvz_newYuvBuffer(in_width, in_height, CHROMA_420, 0);
  yuv_buffer_t* out = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);

  if (yuv_io_read(file, in_width, in_height, data)) {

    _kvzScaling(data, &out);

    FILE* out_file = fopen(out_file_name, "wb");
    yuv_io_write(out_file, out, out->y->width, out->y->height);
    fclose(out_file);
  }
  else {
    perror("Read fail");
    if (feof(file)) perror("End of file");
    if (ferror(file)) perror("File error");
  }

  deallocateYuvBuffer(data);
  deallocateYuvBuffer(out);
  fclose(file);
}

//Compare "rounded" and "unrounded" scaling
void test4()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width = 500;
  int32_t out_height = 550;
  int framerate = 24;

  //const char* file_name_format = "Kimono1_%ix%i_%i.yuv";

  char in_file_name[BUFF_SIZE];
  sprintf(in_file_name, "Kimono1_%ix%i_%i.yuv", in_width, in_height, framerate);
  
  FILE* file = fopen(in_file_name, "rb");
  if (file == NULL) {
    perror("FIle open failed");
    printf("File name: %s", in_file_name);
  }

  yuv_buffer_t* data1 = kvz_newYuvBuffer(in_width, in_height, CHROMA_420, 0);
  yuv_buffer_t* data2 = NULL;
  yuv_buffer_t* out1 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  yuv_buffer_t* out2 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);

  if (yuv_io_read(file, in_width, in_height, data1)) {

    data2 = cloneYuvBuffer(data1);

    kvzScaling(data1, &out1);
    _kvzScaling(data2, &out2);

    yuv_buffer_t* result = yuvDiff(out1, out2);
    printout(result);
    printMeanColor(result);

    double psnr[3] = {0.0,0.0,0.0};
    compute_psnr(out1, out2, psnr);

    printf("PSNR:\n Y: %f\n U: %f\n V: %f\n", psnr[0], psnr[1], psnr[2]);

    deallocateYuvBuffer(result);

    //FILE* out_file = fopen(out_file_name, "wb");
    //yuv_io_write(out_file, out1, out1->y->width, out1->y->height);
    //fclose(out_file);
  }
  else {
    perror("Read fail");
    if (feof(file)) perror("End of file");
    if (ferror(file)) perror("File error");
  }

  deallocateYuvBuffer(data1);
  deallocateYuvBuffer(out1);
  deallocateYuvBuffer(data2);
  deallocateYuvBuffer(out2);
  fclose(file);
}

//Check different widths and heights to make sure the tested methods produce the same result
void test5()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width_beg = 100;
  int32_t out_height_beg = 100;
  int32_t max_width = 3000;
  int32_t max_height = 3000;

  const int w_step = 10;
  const int h_step = 10;
  int framerate = 24;

  //const char* file_name_format = "Kimono1_%ix%i_%i.yuv";

  char in_file_name[BUFF_SIZE];
  sprintf(in_file_name, "Kimono1_%ix%i_%i.yuv", in_width, in_height, framerate);

  FILE* file = fopen(in_file_name, "rb");
  if (file == NULL) {
    perror("FIle open failed");
    printf("File name: %s", in_file_name);
  }

  yuv_buffer_t* data = kvz_newYuvBuffer(in_width, in_height, CHROMA_420, 0);
  yuv_buffer_t* data1 = NULL;
  yuv_buffer_t* data2 = NULL;

  int all_same = 1;

  if (yuv_io_read(file, in_width, in_height, data)) {
    //Iterate over different widths and heights
    for (int w = out_width_beg; w < max_width; w += w_step) {
      for (int h = out_height_beg; h < max_height; h += h_step) {
        data1 = cloneYuvBuffer(data);
        data2 = cloneYuvBuffer(data);

        yuv_buffer_t* out1 = kvz_newYuvBuffer(w, h, CHROMA_420, 0);
        yuv_buffer_t* out2 = kvz_newYuvBuffer(w, h, CHROMA_420, 0);

        kvzScaling(data1, &out1);
        _kvzScaling(data2, &out2);

        if (isSame(out1,out2)) {
          fprintf(stderr, "\rCur size: %dx%d\n", w, h);
        }
        else {
          printf("Not same at %ix%i size.\n", w, h);
          all_same = 0;
        }
        deallocateYuvBuffer(data1);
        deallocateYuvBuffer(data2);
        deallocateYuvBuffer(out1);
        deallocateYuvBuffer(out2);
      }
    }
  }
  else {
    perror("Read fail");
    if (feof(file)) perror("End of file");
    if (ferror(file)) perror("File error");
  }

  if(all_same) printf("All were same.\n");

  deallocateYuvBuffer(data);
  
  fclose(file);
}

//Scale videos
void vscaling()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width = 960;
  int32_t out_height = 540;
  int framerate = 24;
  int frames = 300;

  //const char* file_name_format = "Kimono1_%ix%i_%i.yuv";

  char in_file_name[BUFF_SIZE];
  sprintf(in_file_name, "Kimono1_%ix%i_%i.yuv", in_width, in_height, framerate);

  char out_file_name[BUFF_SIZE];
  sprintf(out_file_name, "Kimono1_%ix%i_%i_s.yuv", out_width, out_height, framerate);
  
  FILE* out_file = fopen(out_file_name, "wb");
  FILE* file = fopen(in_file_name, "rb");
  if (file == NULL || out_file == NULL) {
    perror("FIle open failed");
    printf("File name: %s", in_file_name);
  }

  yuv_buffer_t* data = kvz_newYuvBuffer(in_width, in_height, CHROMA_420, 0);
  yuv_buffer_t* out = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  int i = 0;
  
  while (yuv_io_read(file, in_width, in_height, data) && frames > i) {

      kvzScaling(data, &out);


      yuv_io_write(out_file, out, out->y->width, out->y->height);

      printf("Frame number %i\n",++i);
  }

  deallocateYuvBuffer(data);
  deallocateYuvBuffer(out);
  fclose(file);
  fclose(out_file);
}

int main()
{
  vscaling();

  system("Pause");
  return EXIT_SUCCESS;
}
