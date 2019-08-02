//Main function for testing the scaler class

#include <stdlib.h>
#include <stdio.h>

#include "scaler.h"
#include "scaler-avx2.h"

#include <math.h>
#include <memory.h>

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

static void printPicBuffer(pic_buffer_t* buffer)
{
  for (int i = 0; i < buffer->height; i++) {
    for (int j = 0; j < buffer->width; j++) {
      printf("%i ", buffer->data[buffer->width * i + j]);
    }
    printf("\n");
  }
}

static void printout(yuv_buffer_t* buffer)
{
  printf("Y:\n");
  printPicBuffer(buffer->y);
  printf("Cb:\n");
  printPicBuffer(buffer->u);
  printf("Cr:\n");
  printPicBuffer(buffer->v);
}

//Copy data to a output array
static void copyBack(pic_data_t* dst, uint8_t* src, int size)
{
  for (int i = 0; i < size; i++) {
    dst[i] = (pic_data_t)src[i];
  }
}

static void copyFrom(uint8_t* dst, pic_data_t* src, int size)
{
  for (int i = 0; i < size; i++) {
    dst[i] = (uint8_t)src[i];
  }
}

#define COPY_CUSTOM(dst, src, size, dst_type, src_type) \
  for (int i = 0; i < (size); i++) {\
    ((dst_type *)(dst))[i] = (dst_type)(((src_type *)(src))[i]);\
  }\

static int opaque_piccmp(opaque_pic_buffer_t *b1, opaque_pic_buffer_t *b2)
{
  int ret_val = 1;
  for (int i = 0; i < b1->width * b1->height; i++) {
    long b1_val, b2_val;

    switch (b1->depth)
    {
      case sizeof(char):
        b1_val = (long)(((unsigned char *)(b1->data))[i]);
        break;

      case sizeof(short) :
        b1_val = (long)(((unsigned short *)(b1->data))[i]);
        break;

      case sizeof(int) :
        b1_val = (long)(((unsigned int *)(b1->data))[i]);
        break;

      default:
      break;
    }

    switch (b2->depth)
    {
      case sizeof(char) :
        b2_val = (long)(((unsigned char *)(b2->data))[i]);
        break;

      case sizeof(short) :
        b2_val = (long)(((unsigned short *)(b2->data))[i]);
        break;

      case sizeof(int) :
        b2_val = (long)(((unsigned int *)(b2->data))[i]);
        break;

      default:
        break;
    }

    ret_val &= b1_val == b2_val;
  }

  return ret_val;
}

static int opaque_yuvcmp(opaque_yuv_buffer_t *yuv1, opaque_yuv_buffer_t *yuv2)
{
  int ret_val = 1;
  ret_val &= opaque_piccmp(yuv1->y, yuv2->y);
  ret_val &= opaque_piccmp(yuv1->u, yuv2->u);
  ret_val &= opaque_piccmp(yuv1->v, yuv2->v);
  return ret_val;
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
static int yuv_io_read(FILE* file,
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
static int yuv_io_write(FILE* file,
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
static void kvzScaling(yuv_buffer_t* in, yuv_buffer_t** out)
{
  //Create picture buffers based on given kvz_pictures
  int32_t in_y_width = in->y->width;
  int32_t in_y_height = in->y->height;
  int32_t out_y_width = (*out)->y->width;
  int32_t out_y_height = (*out)->y->height;

  //assumes 420
  //int is_420 = in->y->width != in->u->width ? 1 : 0;
  scaling_parameter_t param = kvz_newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, CHROMA_420);

  *out = kvz_yuvScaling(in, &param, *out);
}

static void kvzScaling_ver(yuv_buffer_t* in, yuv_buffer_t** out, int ver)
{
  resample_func *func;
  switch (ver) {
  case 2:
    func = kvz_alt_resample_func;
    break;

  default:
  case 1:
    func = kvz_default_resample_func;
  }

  //Create picture buffers based on given kvz_pictures
  int32_t in_y_width = in->y->width;
  int32_t in_y_height = in->y->height;
  int32_t out_y_width = (*out)->y->width;
  int32_t out_y_height = (*out)->y->height;

  //assumes 420
  //int is_420 = in->y->width != in->u->width ? 1 : 0;
  scaling_parameter_t param = kvz_newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, CHROMA_420);

  *out = kvz_yuvScaling_adapter(in, &param, *out, func);
}

static void kvzBlockScaling(yuv_buffer_t* in, yuv_buffer_t** out)
{
  //Create picture buffers based on given kvz_pictures
  int32_t in_y_width = in->y->width;
  int32_t in_y_height = in->y->height;
  int32_t out_y_width = (*out)->y->width;
  int32_t out_y_height = (*out)->y->height;
  int part = 1;

  //assumes 420
  //int is_420 = in->y->width != in->u->width ? 1 : 0;
  scaling_parameter_t param = kvz_newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, CHROMA_420);


  int block_width = out_y_width >> part;
  int block_height = out_y_height >> part;

  for (size_t y = 0; y < (1<<part); y++) {
    for (size_t x = 0; x < (1<<part); x++) {
      int block_x = block_width * x;
      int block_y = block_height * y;

      kvz_yuvBlockScaling(in, &param, *out, block_x, block_y,
        block_width, block_height);
    }
  }
}

static void kvzBlockStepScaling(yuv_buffer_t* in, yuv_buffer_t** out)
{
  //Create picture buffers based on given kvz_pictures
  int32_t in_y_width = in->y->width;
  int32_t in_y_height = in->y->height;
  int32_t out_y_width = (*out)->y->width;
  int32_t out_y_height = (*out)->y->height;
  int part = 1;

  //assumes 420
  //int is_420 = in->y->width != in->u->width ? 1 : 0;
  scaling_parameter_t param = kvz_newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, CHROMA_420);


  int block_width = out_y_width >> part;
  int block_height = out_y_height >> part;

  yuv_buffer_t* tmp = kvz_newYuvBuffer(out_y_width, in_y_height, CHROMA_420, 0);

  int32_t in_parts = (in_y_height + block_height - 1) / block_height;

  //Horizontal
  for (size_t y = 0; y < in_parts; y++) {
    for (size_t x = 0; x < (1 << part); x++) {
      int block_x = block_width * x;
      int block_y = block_height * y;
      int bh = min(block_height, in_y_height-block_y);

      kvz_yuvBlockStepScaling(tmp, in, &param, block_x, block_y,
        block_width, bh, 0);
    }
  }

  //Vertical
  for (size_t y = 0; y < (1 << part); y++) {
    for (size_t x = 0; x < (1 << part); x++) {
      int block_x = block_width * x;
      int block_y = block_height * y;

      kvz_yuvBlockStepScaling(*out, tmp, &param, block_x, block_y,
        block_width, block_height, 1);
    }
  }

  kvz_deallocateYuvBuffer(tmp);
}

static void kvzOpaqueBlockStepScaling(opaque_yuv_buffer_t* in, opaque_yuv_buffer_t** out, int tmp_depth)
{
  //Create picture buffers based on given kvz_pictures
  int32_t in_y_width = in->y->width;
  int32_t in_y_height = in->y->height;
  int32_t out_y_width = (*out)->y->width;
  int32_t out_y_height = (*out)->y->height;
  int part = 1;

  //assumes 420
  //int is_420 = in->y->width != in->u->width ? 1 : 0;
  scaling_parameter_t param = kvz_newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, CHROMA_420);

  int block_width = out_y_width >> part;
  int block_height = out_y_height >> part;

  int luma_size = out_y_width * in_y_height;
  int chroma_size = luma_size >> 2;

  opaque_yuv_buffer_t* tmp;

  switch (tmp_depth)
  {
    case sizeof(pic_data_t) :
      tmp = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, out_y_width, in_y_height, out_y_width, CHROMA_420, sizeof(pic_data_t));
      break;

    case sizeof(short) :
      tmp = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, out_y_width, in_y_height, out_y_width, CHROMA_420, sizeof(short));
      break;

  default:
    return;
    break;
  }

  int32_t in_parts = (in_y_height + block_height - 1) / block_height;

  //Horizontal
  for (size_t y = 0; y < in_parts; y++) {
    for (size_t x = 0; x < (1 << part); x++) {
      int block_x = block_width * x;
      int block_y = block_height * y;
      int bh = min(block_height, in_y_height - block_y);

      kvz_opaqueYuvBlockStepScaling_adapter(tmp, in, &param, block_x, block_y,
        block_width, bh, 0, kvz_opaque_block_step_resample_func);
    }
  }

  //Vertical
  for (size_t y = 0; y < (1 << part); y++) {
    for (size_t x = 0; x < (1 << part); x++) {
      int block_x = block_width * x;
      int block_y = block_height * y;

      kvz_opaqueYuvBlockStepScaling_adapter(*out, tmp, &param, block_x, block_y,
        block_width, block_height, 1, kvz_opaque_block_step_resample_func);
    }
  }

  kvz_deallocateOpaqueYuvBuffer(tmp, 1);
  //kvz_deallocateYuvBuffer(tmp);
}

static void kvzOpaqueBlockStepScalingAvx2(opaque_yuv_buffer_t* in, opaque_yuv_buffer_t** out, int tmp_depth)
{
  //Create picture buffers based on given kvz_pictures
  int32_t in_y_width = in->y->width;
  int32_t in_y_height = in->y->height;
  int32_t out_y_width = (*out)->y->width;
  int32_t out_y_height = (*out)->y->height;
  int part = 1;

  //assumes 420
  //int is_420 = in->y->width != in->u->width ? 1 : 0;
  scaling_parameter_t param = kvz_newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, CHROMA_420);

  int block_width = out_y_width >> part;
  int block_height = out_y_height >> part;

  int luma_size = out_y_width * in_y_height;
  int chroma_size = luma_size >> 2;

  opaque_yuv_buffer_t* tmp;

  switch (tmp_depth)
  {
    case sizeof(pic_data_t) :
      tmp = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, out_y_width, in_y_height, out_y_width, CHROMA_420, sizeof(pic_data_t));
      break;

      case sizeof(short) :
        tmp = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, out_y_width, in_y_height, out_y_width, CHROMA_420, sizeof(short));
        break;

      default:
        return;
        break;
  }

  int32_t in_parts = (in_y_height + block_height - 1) / block_height;

  //Horizontal
  for (size_t y = 0; y < in_parts; y++) {
    for (size_t x = 0; x < (1 << part); x++) {
      int block_x = block_width * x;
      int block_y = block_height * y;
      int bh = min(block_height, in_y_height - block_y);

      kvz_opaqueYuvBlockStepScaling_adapter(tmp, in, &param, block_x, block_y,
        block_width, bh, 0, kvz_opaque_block_step_resample_func_avx2);
    }
  }

  //Vertical
  for (size_t y = 0; y < (1 << part); y++) {
    for (size_t x = 0; x < (1 << part); x++) {
      int block_x = block_width * x;
      int block_y = block_height * y;

      kvz_opaqueYuvBlockStepScaling_adapter(*out, tmp, &param, block_x, block_y,
        block_width, block_height, 1, kvz_opaque_block_step_resample_func_avx2);
    }
  }

  kvz_deallocateOpaqueYuvBuffer(tmp, 1);
  //kvz_deallocateYuvBuffer(tmp);
}

static void kvzScaling_avx2(yuv_buffer_t* in, yuv_buffer_t** out)
{
  //Create picture buffers based on given kvz_pictures
  int32_t in_y_width = in->y->width;
  int32_t in_y_height = in->y->height;
  int32_t out_y_width = (*out)->y->width;
  int32_t out_y_height = (*out)->y->height;

  //assumes 420
  //int is_420 = in->y->width != in->u->width ? 1 : 0;
  scaling_parameter_t param = kvz_newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, CHROMA_420);

  *out = kvz_yuvScaling_adapter(in, &param, *out, kvz_default_resample_func);
}

static void kvzBlockStepScaling_avx2(yuv_buffer_t* in, yuv_buffer_t** out, int ver)
{
  //Create picture buffers based on given kvz_pictures
  int32_t in_y_width = in->y->width;
  int32_t in_y_height = in->y->height;
  int32_t out_y_width = (*out)->y->width;
  int32_t out_y_height = (*out)->y->height;
  int part = 1;

  resample_block_step_func *func;
  switch (ver) {
  case 3:
    func = kvz_alt2_block_step_resample_func_avx2;
    break;

  case 2:
    func = kvz_alt1_block_step_resample_func_avx2;
    break;

  default:
  case 1:
    func = kvz_default_block_step_resample_func_avx2;
  }

  //assumes 420
  //int is_420 = in->y->width != in->u->width ? 1 : 0;
  scaling_parameter_t param = kvz_newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, CHROMA_420);


  int block_width = out_y_width >> part;
  int block_height = out_y_height >> part;

  yuv_buffer_t* tmp = kvz_newYuvBuffer(out_y_width, in_y_height, CHROMA_420, 0);

  int32_t in_parts = (in_y_height + block_height - 1) / block_height;

  //Horizontal
  for (size_t y = 0; y < in_parts; y++) {
    for (size_t x = 0; x < (1 << part); x++) {
      int block_x = block_width * x;
      int block_y = block_height * y;
      int bh = min(block_height, in_y_height - block_y);

      kvz_yuvBlockStepScaling_adapter(tmp, in, &param, block_x, block_y,
        block_width, bh, 0, func);
    }
  }

  //Vertical
  for (size_t y = 0; y < (1 << part); y++) {
    for (size_t x = 0; x < (1 << part); x++) {
      int block_x = block_width * x;
      int block_y = block_height * y;

      kvz_yuvBlockStepScaling_adapter(*out, tmp, &param, block_x, block_y,
        block_width, block_height, 1, func);
    }
  }

  kvz_deallocateYuvBuffer(tmp);
}

static void _kvzScaling(yuv_buffer_t* in, yuv_buffer_t** out)
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

static pic_buffer_t* diffComponent(pic_buffer_t* pic1, pic_buffer_t* pic2)
{
  int width = pic1->width;
  int height = pic1->height;
  int row = 0;
  pic_buffer_t* res = kvz_newPictureBuffer(width, height, 0);

  for( int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      res->data[row + j] = pic1->data[row + j] - pic2->data[row + j];
    }
    row += width;
  }

  return res;
}

static yuv_buffer_t* yuvDiff(yuv_buffer_t* yuv1, yuv_buffer_t* yuv2)
{
  yuv_buffer_t* res = kvz_newYuvBuffer(yuv1->y->width, yuv1->y->height, yuv1->format, 0);

  res->y = diffComponent(yuv1->y, yuv2->y);
  res->u = diffComponent(yuv1->u, yuv2->u);
  res->v = diffComponent(yuv1->v, yuv2->v);

  return res;
}

static double meanColorComponent(pic_buffer_t* pic)
{
  int size = pic->width*pic->height;
  long int sum = 0;

  for(int i = 0; i < size; i++) {
    sum += abs(pic->data[i]);
  }
  return (double)sum/(double)size;
}

static void printMeanColor(yuv_buffer_t* yuv)
{
  printf("Y avg color: %f\n",meanColorComponent(yuv->y));
  printf("U avg color: %f\n",meanColorComponent(yuv->u));
  printf("V avg color: %f\n",meanColorComponent(yuv->v));
}

static double sumComponent(pic_buffer_t* pic)
{
  int size = pic->width*pic->height;
  long int sum = 0;

  for(int i = 0; i < size; i++) {
    sum += abs(pic->data[i]);
  }
  return (double)sum;
}

static int isZero(yuv_buffer_t* yuv)
{
  return sumComponent(yuv->y) != 0 ? 0 : (sumComponent(yuv->u) != 0 ? 0 : (sumComponent(yuv->v) != 0 ? 0 : 1));
}

static int isSame(yuv_buffer_t* yuv1, yuv_buffer_t* yuv2 )
{
  yuv_buffer_t* result = yuvDiff(yuv1, yuv2);
  
  int res = isZero(result);

  kvz_deallocateYuvBuffer(result);

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

static void test3()
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

  kvz_deallocateYuvBuffer(data);
  kvz_deallocateYuvBuffer(out);
  fclose(file);
}

//Compare "rounded" and "unrounded" scaling
static void test4()
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

    data2 = kvz_cloneYuvBuffer(data1);

    kvzScaling(data1, &out1);
    _kvzScaling(data2, &out2);

    yuv_buffer_t* result = yuvDiff(out1, out2);
    printout(result);
    printMeanColor(result);

    double psnr[3] = {0.0,0.0,0.0};
    compute_psnr(out1, out2, psnr);

    printf("PSNR:\n Y: %f\n U: %f\n V: %f\n", psnr[0], psnr[1], psnr[2]);

    kvz_deallocateYuvBuffer(result);

    //FILE* out_file = fopen(out_file_name, "wb");
    //yuv_io_write(out_file, out1, out1->y->width, out1->y->height);
    //fclose(out_file);
  }
  else {
    perror("Read fail");
    if (feof(file)) perror("End of file");
    if (ferror(file)) perror("File error");
  }

  kvz_deallocateYuvBuffer(data1);
  kvz_deallocateYuvBuffer(out1);
  kvz_deallocateYuvBuffer(data2);
  kvz_deallocateYuvBuffer(out2);
  fclose(file);
}

//Check different widths and heights to make sure the tested methods produce the same result
static void test5()
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
        data1 = kvz_cloneYuvBuffer(data);
        data2 = kvz_cloneYuvBuffer(data);

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
        kvz_deallocateYuvBuffer(data1);
        kvz_deallocateYuvBuffer(data2);
        kvz_deallocateYuvBuffer(out1);
        kvz_deallocateYuvBuffer(out2);
      }
    }
  }
  else {
    perror("Read fail");
    if (feof(file)) perror("End of file");
    if (ferror(file)) perror("File error");
  }

  if(all_same) printf("All were same.\n");

  kvz_deallocateYuvBuffer(data);
  
  fclose(file);
}

//Scale videos
static void vscaling()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width = 264;
  int32_t out_height = 130;
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

  kvz_deallocateYuvBuffer(data);
  kvz_deallocateYuvBuffer(out);
  fclose(file);
  fclose(out_file);
}

static void validate_test()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width = in_width << 1;//264;
  int32_t out_height = in_height << 1;//130;
  int32_t out_chroma_width = out_width >> 1;
  int32_t out_chroma_height = out_height >> 1;
  int framerate = 24;
  int frames = 10;

  //const char* file_name_format = "Kimono1_%ix%i_%i.yuv";

  char in_file_name[BUFF_SIZE];
  sprintf(in_file_name, "Kimono1_%ix%i_%i.yuv", in_width, in_height, framerate);

  char out_file_name1[BUFF_SIZE];
  char out_file_name2[BUFF_SIZE];
  sprintf(out_file_name1, "Kimono1_%ix%i_%i_ref.yuv", out_width, out_height, framerate);
  sprintf(out_file_name2, "Kimono1_%ix%i_%i_block.yuv", out_width, out_height, framerate);

  FILE* out_file1 = fopen(out_file_name1, "wb");
  FILE* out_file2 = fopen(out_file_name2, "wb");
  FILE* file = fopen(in_file_name, "rb");
  if (file == NULL || out_file1 == NULL || out_file2 == NULL) {
    perror("File open failed");
    printf("File name: %s", in_file_name);
  }

  yuv_buffer_t* data = kvz_newYuvBuffer(in_width, in_height, CHROMA_420, 0);
  yuv_buffer_t* out1 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  yuv_buffer_t* out2 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  int i = 0;

  while (yuv_io_read(file, in_width, in_height, data) && frames > i) {

    kvzScaling(data, &out1);
    kvzBlockScaling(data, &out2);

    if( memcmp(out1->y->data, out2->y->data, sizeof(pic_data_t)*out_width*out_height) != 0 || memcmp(out1->u->data, out2->u->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0 || memcmp(out1->v->data, out2->v->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0 ){
      printf("Frame %i differs\n", i+1);
    }

    yuv_io_write(out_file1, out1, out1->y->width, out1->y->height);
    yuv_io_write(out_file2, out2, out2->y->width, out2->y->height);

    printf("Frame number %i\r", ++i);
  }
  printf("Wrote %i frames.\n", i);

  kvz_deallocateYuvBuffer(data);
  kvz_deallocateYuvBuffer(out1);
  kvz_deallocateYuvBuffer(out2);
  fclose(file);
  fclose(out_file1);
  fclose(out_file2);

}

static void validate_test2()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width = in_width << 1;//264;
  int32_t out_height = in_height << 1;//130;
  int32_t out_chroma_width = out_width >> 1;
  int32_t out_chroma_height = out_height >> 1;
  int framerate = 24;
  int frames = 10;

  //const char* file_name_format = "Kimono1_%ix%i_%i.yuv";

  char in_file_name[BUFF_SIZE];
  sprintf(in_file_name, "Kimono1_%ix%i_%i.yuv", in_width, in_height, framerate);

  char out_file_name1[BUFF_SIZE];
  char out_file_name2[BUFF_SIZE];
  sprintf(out_file_name1, "Kimono1_%ix%i_%i_ref.yuv", out_width, out_height, framerate);
  sprintf(out_file_name2, "Kimono1_%ix%i_%i_block.yuv", out_width, out_height, framerate);

  FILE* out_file1 = fopen(out_file_name1, "wb");
  FILE* out_file2 = fopen(out_file_name2, "wb");
  FILE* file = fopen(in_file_name, "rb");
  if (file == NULL || out_file1 == NULL || out_file2 == NULL) {
    perror("File open failed");
    printf("File name: %s", in_file_name);
  }

  yuv_buffer_t* data = kvz_newYuvBuffer(in_width, in_height, CHROMA_420, 0);
  yuv_buffer_t* out1 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  yuv_buffer_t* out2 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  int i = 0;

  while (yuv_io_read(file, in_width, in_height, data) && frames > i) {

    kvzScaling(data, &out1);
    kvzBlockStepScaling(data, &out2);

    if (memcmp(out1->y->data, out2->y->data, sizeof(pic_data_t)*out_width*out_height) != 0 || memcmp(out1->u->data, out2->u->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0 || memcmp(out1->v->data, out2->v->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0) {
      printf("Frame %i differs\n", i + 1);
    }

    yuv_io_write(out_file1, out1, out1->y->width, out1->y->height);
    yuv_io_write(out_file2, out2, out2->y->width, out2->y->height);

    printf("Frame number %i\r", ++i);
  }
  printf("Wrote %i frames.\n", i);

  kvz_deallocateYuvBuffer(data);
  kvz_deallocateYuvBuffer(out1);
  kvz_deallocateYuvBuffer(out2);
  fclose(file);
  fclose(out_file1);
  fclose(out_file2);

}

static void validate_test3()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width = in_width << 1;//264;
  int32_t out_height = in_height << 1;//130;
  int32_t out_chroma_width = out_width >> 1;
  int32_t out_chroma_height = out_height >> 1;
  int framerate = 24;
  int frames = 10;

  //const char* file_name_format = "Kimono1_%ix%i_%i.yuv";

  char in_file_name[BUFF_SIZE];
  sprintf(in_file_name, "Kimono1_%ix%i_%i.yuv", in_width, in_height, framerate);

  char out_file_name1[BUFF_SIZE];
  char out_file_name11[BUFF_SIZE];
  char out_file_name2[BUFF_SIZE];
  char out_file_name3[BUFF_SIZE];
  char out_file_name4[BUFF_SIZE];
  sprintf(out_file_name1, "Kimono1_%ix%i_%i_ref.yuv", out_width, out_height, framerate);
  sprintf(out_file_name11, "Kimono1_%ix%i_%i_avx2.yuv", out_width, out_height, framerate);
  sprintf(out_file_name2, "Kimono1_%ix%i_%i_block_avx2.yuv", out_width, out_height, framerate);
  sprintf(out_file_name3, "Kimono1_%ix%i_%i_block_avx2_2.yuv", out_width, out_height, framerate);
  sprintf(out_file_name4, "Kimono1_%ix%i_%i_block_avx2_3.yuv", out_width, out_height, framerate);

  FILE* out_file1 = fopen(out_file_name1, "wb");
  FILE* out_file11 = fopen(out_file_name11, "wb");
  FILE* out_file2 = fopen(out_file_name2, "wb");
  FILE* out_file3 = fopen(out_file_name3, "wb");
  FILE* out_file4 = fopen(out_file_name4, "wb");
  FILE* file = fopen(in_file_name, "rb");
  if (file == NULL || out_file1 == NULL || out_file2 == NULL || out_file3 == NULL || out_file11 == NULL || out_file4 == NULL) {
    perror("File open failed");
    printf("File name: %s", in_file_name);
  }

  yuv_buffer_t* data = kvz_newYuvBuffer(in_width, in_height, CHROMA_420, 0);
  yuv_buffer_t* out1 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  yuv_buffer_t* out11 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  yuv_buffer_t* out2 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  yuv_buffer_t* out3 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  yuv_buffer_t* out4 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  int i = 0;

  while (yuv_io_read(file, in_width, in_height, data) && frames > i) {

    kvzScaling(data, &out1);
    kvzScaling_avx2(data, &out11);
    kvzBlockStepScaling_avx2(data, &out2, 1);
    kvzBlockStepScaling_avx2(data, &out3, 2);
    kvzBlockStepScaling_avx2(data, &out4, 3);


    if (memcmp(out1->y->data, out11->y->data, sizeof(pic_data_t)*out_width*out_height) != 0 || memcmp(out1->u->data, out11->u->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0 || memcmp(out1->v->data, out11->v->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0) {
      printf("Frame %i differs in avx2\n", i + 1);
    }

    if (memcmp(out1->y->data, out2->y->data, sizeof(pic_data_t)*out_width*out_height) != 0 || memcmp(out1->u->data, out2->u->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0 || memcmp(out1->v->data, out2->v->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0) {
      printf("Frame %i differs in block avx2\n", i + 1);
    }

    if (memcmp(out1->y->data, out3->y->data, sizeof(pic_data_t)*out_width*out_height) != 0 || memcmp(out1->u->data, out3->u->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0 || memcmp(out1->v->data, out3->v->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0) {
      printf("Frame %i differs in block avx2_2\n", i + 1);
    }

    if (memcmp(out1->y->data, out4->y->data, sizeof(pic_data_t)*out_width*out_height) != 0 || memcmp(out1->u->data, out4->u->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0 || memcmp(out1->v->data, out4->v->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0) {
      printf("Frame %i differs in block avx2_3\n", i + 1);
    }

    yuv_io_write(out_file1, out1, out1->y->width, out1->y->height);
    yuv_io_write(out_file11, out11, out11->y->width, out11->y->height);
    yuv_io_write(out_file2, out2, out2->y->width, out2->y->height);
    yuv_io_write(out_file3, out3, out3->y->width, out3->y->height);
    yuv_io_write(out_file4, out4, out4->y->width, out4->y->height);

    printf("Frame number %i\r", ++i);
  }
  printf("Wrote %i frames.\n", i);

  kvz_deallocateYuvBuffer(data);
  kvz_deallocateYuvBuffer(out1);
  kvz_deallocateYuvBuffer(out11);
  kvz_deallocateYuvBuffer(out2);
  kvz_deallocateYuvBuffer(out3);
  kvz_deallocateYuvBuffer(out4);
  fclose(file);
  fclose(out_file1);
  fclose(out_file11);
  fclose(out_file2);
  fclose(out_file3);
  fclose(out_file4);

}

static void opaque_validate_test_avx2()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width = in_width << 1;//264;
  int32_t out_height = in_height << 1;//130;
  int32_t out_chroma_width = out_width >> 1;
  int32_t out_chroma_height = out_height >> 1;
  int framerate = 24;
  int frames = 10;

  int is_downscaling = out_width < in_width;
  //const char* file_name_format = "Kimono1_%ix%i_%i.yuv";

  char in_file_name[BUFF_SIZE];
  sprintf(in_file_name, "Kimono1_%ix%i_%i.yuv", in_width, in_height, framerate);

  char out_file_name1[BUFF_SIZE];
  char out_file_name2[BUFF_SIZE];
  char out_file_name3[BUFF_SIZE];
  char out_file_name4[BUFF_SIZE];
  sprintf(out_file_name1, "Kimono1_%ix%i_%i_ref.yuv", out_width, out_height, framerate);
  sprintf(out_file_name2, "Kimono1_%ix%i_%i_block_char_char.yuv", out_width, out_height, framerate);
  sprintf(out_file_name3, "Kimono1_%ix%i_%i_block_short_short.yuv", out_width, out_height, framerate);
  sprintf(out_file_name4, "Kimono1_%ix%i_%i_block_int_int.yuv", out_width, out_height, framerate);

  FILE* out_file1 = fopen(out_file_name1, "wb");
  FILE* out_file2 = fopen(out_file_name2, "wb");
  FILE* out_file3 = fopen(out_file_name3, "wb");
  FILE* out_file4 = fopen(out_file_name4, "wb");
  FILE* file = fopen(in_file_name, "rb");
  if (file == NULL || out_file1 == NULL || out_file2 == NULL || out_file3 == NULL || out_file4 == NULL) {
    perror("File open failed");
    printf("File name: %s", in_file_name);
  }

  yuv_buffer_t* data = kvz_newYuvBuffer(in_width, in_height, CHROMA_420, 0);
  opaque_yuv_buffer_t* data2 = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, in_width, in_height, in_width, CHROMA_420, sizeof(char));
  opaque_yuv_buffer_t* data3 = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, in_width, in_height, in_width, CHROMA_420, sizeof(short));
  opaque_yuv_buffer_t* data4 = kvz_newOpaqueYuvBuffer(data->y->data, data->u->data, data->v->data, in_width, in_height, in_width, CHROMA_420, sizeof(pic_data_t));
  yuv_buffer_t* out1 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  opaque_yuv_buffer_t* op_out1 = kvz_newOpaqueYuvBuffer(out1->y->data, out1->u->data, out1->v->data, out_width, out_height, out_width, CHROMA_420, sizeof(pic_data_t));
  opaque_yuv_buffer_t* out2 = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, out_width, out_height, out_width, CHROMA_420, sizeof(char));
  opaque_yuv_buffer_t* out3 = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, out_width, out_height, out_width, CHROMA_420, sizeof(short));
  opaque_yuv_buffer_t* out4 = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, out_width, out_height, out_width, CHROMA_420, sizeof(pic_data_t));
  int i = 0;

  opaque_yuv_buffer_t* out_array[3] = { out2, out3, out4 };

  while (yuv_io_read(file, in_width, in_height, data) && frames > i) {

    //COPY_CUSTOM(data2->y->data, data->y->data, in_width * in_height, unsigned char, pic_data_t);
    //COPY_CUSTOM(data2->u->data, data->u->data, (in_width * in_height) >> 2, unsigned char, pic_data_t);
    //COPY_CUSTOM(data2->v->data, data->v->data, (in_width * in_height) >> 2, unsigned char, pic_data_t);
    copyFrom((uint8_t *)data2->y->data, data->y->data, in_width * in_height);
    copyFrom((uint8_t *)data2->u->data, data->u->data, (in_width * in_height) >> 2);
    copyFrom((uint8_t *)data2->v->data, data->v->data, (in_width * in_height) >> 2);

    COPY_CUSTOM(data3->y->data, data->y->data, in_width * in_height, short, pic_data_t);
    COPY_CUSTOM(data3->u->data, data->u->data, (in_width * in_height) >> 2, short, pic_data_t);
    COPY_CUSTOM(data3->v->data, data->v->data, (in_width * in_height) >> 2, short, pic_data_t);

    kvzScaling(data, &out1);

    for (int k = 0; k < 3; k++)
    {
      for (int j = is_downscaling ? sizeof(pic_data_t) : sizeof(short); j <= sizeof(pic_data_t); j++)
      {
        kvzOpaqueBlockStepScalingAvx2(data2, &out_array[k], j);
        kvzOpaqueBlockStepScalingAvx2(data3, &out_array[(k + 1) % 3], j);
        kvzOpaqueBlockStepScalingAvx2(data4, &out_array[(k + 2) % 3], j);

        if (opaque_yuvcmp(op_out1, out2) == 0) {
          printf("Frame %i differs. In depth %i, tmp depth %i, out depth 1\n", i + 1, (k % 3) + 1, j);
        }

        if (opaque_yuvcmp(op_out1, out3) == 0) {
          printf("Frame %i differs. In depth %i, tmp depth %i, out depth 2\n", i + 1, (k + 1) % 3 + 1, j);
        }

        if (opaque_yuvcmp(op_out1, out4) == 0) {
          printf("Frame %i differs. In depth %i, tmp depth %i, out depth 3\n", i + 1, (k + 2) % 3 + 1, j);
        }
      }
    }
    //yuv_io_write(out_file1, out1, out1->y->width, out1->y->height);
    //yuv_io_write(out_file2, (yuv_buffer_t *)out2, out2->y->width, out2->y->height);
    //yuv_io_write(out_file3, (yuv_buffer_t *)out3, out3->y->width, out3->y->height);
    //yuv_io_write(out_file4, (yuv_buffer_t *)out4, out4->y->width, out4->y->height);

    printf("Frame number %i\r", ++i);
  }
  printf("Wrote %i frames.\n", i);

  kvz_deallocateYuvBuffer(data);
  kvz_deallocateOpaqueYuvBuffer(data2, 1);
  kvz_deallocateOpaqueYuvBuffer(data3, 1);
  kvz_deallocateOpaqueYuvBuffer(data4, 0);
  kvz_deallocateYuvBuffer(out1);
  kvz_deallocateOpaqueYuvBuffer(op_out1, 0);
  kvz_deallocateOpaqueYuvBuffer(out2, 1);
  kvz_deallocateOpaqueYuvBuffer(out3, 1);
  kvz_deallocateOpaqueYuvBuffer(out4, 1);
  fclose(file);
  fclose(out_file1);
  fclose(out_file2);
  fclose(out_file3);
  fclose(out_file4);
}

static void opaque_validate_test()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width = in_width << 1;//264;
  int32_t out_height = in_height << 1;//130;
  int32_t out_chroma_width = out_width >> 1;
  int32_t out_chroma_height = out_height >> 1;
  int framerate = 24;
  int frames = 10;

  int is_downscaling = out_width < in_width;
  //const char* file_name_format = "Kimono1_%ix%i_%i.yuv";

  char in_file_name[BUFF_SIZE];
  sprintf(in_file_name, "Kimono1_%ix%i_%i.yuv", in_width, in_height, framerate);

  char out_file_name1[BUFF_SIZE];
  char out_file_name2[BUFF_SIZE];
  char out_file_name3[BUFF_SIZE];
  char out_file_name4[BUFF_SIZE];
  sprintf(out_file_name1, "Kimono1_%ix%i_%i_ref.yuv", out_width, out_height, framerate);
  sprintf(out_file_name2, "Kimono1_%ix%i_%i_block_char_char.yuv", out_width, out_height, framerate);
  sprintf(out_file_name3, "Kimono1_%ix%i_%i_block_short_short.yuv", out_width, out_height, framerate);
  sprintf(out_file_name4, "Kimono1_%ix%i_%i_block_int_int.yuv", out_width, out_height, framerate);

  FILE* out_file1 = fopen(out_file_name1, "wb");
  FILE* out_file2 = fopen(out_file_name2, "wb");
  FILE* out_file3 = fopen(out_file_name3, "wb");
  FILE* out_file4 = fopen(out_file_name4, "wb");
  FILE* file = fopen(in_file_name, "rb");
  if (file == NULL || out_file1 == NULL || out_file2 == NULL || out_file3 == NULL || out_file4 == NULL) {
    perror("File open failed");
    printf("File name: %s", in_file_name);
  }

  yuv_buffer_t* data = kvz_newYuvBuffer(in_width, in_height, CHROMA_420, 0);
  opaque_yuv_buffer_t* data2 = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, in_width, in_height, in_width, CHROMA_420, sizeof(char));
  opaque_yuv_buffer_t* data3 = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, in_width, in_height, in_width, CHROMA_420, sizeof(short));
  opaque_yuv_buffer_t* data4 = kvz_newOpaqueYuvBuffer(data->y->data, data->u->data, data->v->data, in_width, in_height, in_width, CHROMA_420, sizeof(pic_data_t));
  yuv_buffer_t* out1 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  opaque_yuv_buffer_t* op_out1 = kvz_newOpaqueYuvBuffer(out1->y->data, out1->u->data, out1->v->data, out_width, out_height, out_width, CHROMA_420, sizeof(pic_data_t));
  opaque_yuv_buffer_t* out2 = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, out_width, out_height, out_width, CHROMA_420, sizeof(char));
  opaque_yuv_buffer_t* out3 = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, out_width, out_height, out_width, CHROMA_420, sizeof(short));
  opaque_yuv_buffer_t* out4 = kvz_newOpaqueYuvBuffer(NULL, NULL, NULL, out_width, out_height, out_width, CHROMA_420, sizeof(pic_data_t));
  int i = 0;

  opaque_yuv_buffer_t* out_array[3] = { out2, out3, out4 };

  while (yuv_io_read(file, in_width, in_height, data) && frames > i) {

    //COPY_CUSTOM(data2->y->data, data->y->data, in_width * in_height, unsigned char, pic_data_t);
    //COPY_CUSTOM(data2->u->data, data->u->data, (in_width * in_height) >> 2, unsigned char, pic_data_t);
    //COPY_CUSTOM(data2->v->data, data->v->data, (in_width * in_height) >> 2, unsigned char, pic_data_t);
    copyFrom((uint8_t *)data2->y->data, data->y->data, in_width * in_height);
    copyFrom((uint8_t *)data2->u->data, data->u->data, (in_width * in_height) >> 2);
    copyFrom((uint8_t *)data2->v->data, data->v->data, (in_width * in_height) >> 2);

    COPY_CUSTOM(data3->y->data, data->y->data, in_width * in_height, short, pic_data_t);
    COPY_CUSTOM(data3->u->data, data->u->data, (in_width * in_height) >> 2, short, pic_data_t);
    COPY_CUSTOM(data3->v->data, data->v->data, (in_width * in_height) >> 2, short, pic_data_t);

    kvzScaling(data, &out1);
    
    for (int k = 0; k < 3; k++)
    {
      for (int j = is_downscaling ? sizeof(pic_data_t) : sizeof(short); j <= sizeof(pic_data_t); j++)
      {
        kvzOpaqueBlockStepScaling(data2, &out_array[k], j);
        kvzOpaqueBlockStepScaling(data3, &out_array[(k+1)%3], j);
        kvzOpaqueBlockStepScaling(data4, &out_array[(k+2)%3], j);

        if (opaque_yuvcmp(op_out1, out2) == 0) {
          printf("Frame %i differs. In depth %i, tmp depth %i, out depth 1\n", i + 1, (k % 3) + 1, j);
        }

        if (opaque_yuvcmp(op_out1, out3) == 0) {
          printf("Frame %i differs. In depth %i, tmp depth %i, out depth 2\n", i + 1, (k + 1) % 3 + 1, j);
        }

        if (opaque_yuvcmp(op_out1, out4) == 0) {
          printf("Frame %i differs. In depth %i, tmp depth %i, out depth 3\n", i + 1, (k + 2) % 3 + 1, j);
        }
      }
    }
    //yuv_io_write(out_file1, out1, out1->y->width, out1->y->height);
    //yuv_io_write(out_file2, (yuv_buffer_t *)out2, out2->y->width, out2->y->height);
    //yuv_io_write(out_file3, (yuv_buffer_t *)out3, out3->y->width, out3->y->height);
    //yuv_io_write(out_file4, (yuv_buffer_t *)out4, out4->y->width, out4->y->height);

    printf("Frame number %i\r", ++i);
  }
  printf("Wrote %i frames.\n", i);

  kvz_deallocateYuvBuffer(data);
  kvz_deallocateOpaqueYuvBuffer(data2, 1);
  kvz_deallocateOpaqueYuvBuffer(data3, 1);
  kvz_deallocateOpaqueYuvBuffer(data4, 0);
  kvz_deallocateYuvBuffer(out1);
  kvz_deallocateOpaqueYuvBuffer(op_out1, 0);
  kvz_deallocateOpaqueYuvBuffer(out2, 1);
  kvz_deallocateOpaqueYuvBuffer(out3, 1);
  kvz_deallocateOpaqueYuvBuffer(out4, 1);
  fclose(file);
  fclose(out_file1);
  fclose(out_file2);
  fclose(out_file3);
  fclose(out_file4);
}

static void validate_test4()
{
  int32_t in_width = 1920;
  int32_t in_height = 1080;
  int32_t out_width = in_width << 1;//264;
  int32_t out_height = in_height << 1;//130;
  int32_t out_chroma_width = out_width >> 1;
  int32_t out_chroma_height = out_height >> 1;
  int framerate = 24;
  int frames = 10;

  //const char* file_name_format = "Kimono1_%ix%i_%i.yuv";

  char in_file_name[BUFF_SIZE];
  sprintf(in_file_name, "Kimono1_%ix%i_%i.yuv", in_width, in_height, framerate);

  char out_file_name1[BUFF_SIZE];
  char out_file_name2[BUFF_SIZE];
  sprintf(out_file_name1, "Kimono1_%ix%i_%i_ref.yuv", out_width, out_height, framerate);
  sprintf(out_file_name2, "Kimono1_%ix%i_%i_block.yuv", out_width, out_height, framerate);

  FILE* out_file1 = fopen(out_file_name1, "wb");
  FILE* out_file2 = fopen(out_file_name2, "wb");
  FILE* file = fopen(in_file_name, "rb");
  if (file == NULL || out_file1 == NULL || out_file2 == NULL) {
    perror("File open failed");
    printf("File name: %s", in_file_name);
  }

  yuv_buffer_t* data = kvz_newYuvBuffer(in_width, in_height, CHROMA_420, 0);
  yuv_buffer_t* out1 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  yuv_buffer_t* out2 = kvz_newYuvBuffer(out_width, out_height, CHROMA_420, 0);
  int i = 0;

  while (yuv_io_read(file, in_width, in_height, data) && frames > i) {

    kvzScaling_ver(data, &out1, 1);
    kvzScaling_ver(data, &out2, 2);

    if (memcmp(out1->y->data, out2->y->data, sizeof(pic_data_t)*out_width*out_height) != 0 || memcmp(out1->u->data, out2->u->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0 || memcmp(out1->v->data, out2->v->data, sizeof(pic_data_t)*out_chroma_height*out_chroma_width) != 0) {
      printf("Frame %i differs\n", i + 1);
    }

    yuv_io_write(out_file1, out1, out1->y->width, out1->y->height);
    yuv_io_write(out_file2, out2, out2->y->width, out2->y->height);

    printf("Frame number %i\r", ++i);
  }
  printf("Wrote %i frames.\n", i);

  kvz_deallocateYuvBuffer(data);
  kvz_deallocateYuvBuffer(out1);
  kvz_deallocateYuvBuffer(out2);
  fclose(file);
  fclose(out_file1);
  fclose(out_file2);

}

int main()
{
  //vscaling();
  //validate_test3();
  opaque_validate_test_avx2();
  //int r = test_avx();
  //printf("%d", r);

  system("Pause");
  return EXIT_SUCCESS;
}
