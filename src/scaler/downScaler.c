/*
 * downScaler.c - Matlab integration for testing downscaling
 * 
 * Take an input YCbCr data array (NxM) and scale it to the
 * specified size. Size should be of the form [width height]
 *
 * Calling syntax:
 *    [Y,Cb,Cr] = downScaler( inYdata, outYsize,
 *                            inCbData, outCbSize,
 *                            inCrData, outCrSize )
*/
#include "mex.h"
#include "scaler.h"
#include "scaler.c"
#include "stdlib.h"
#include "matrix.h"

#define STR(s) #s
#define ERROR_(type,msg) mexErrMsgIdAndTxt(STR(type),msg)
#define ERROR(type,msg) ERROR_(MyToolbox:arrayProduct:##type,msg)

//Copy data to a output array
void copyBack(uint8_t* dst, int* src, int size)
{
  for (int i = 0; i < size; i++) {
    dst[i] = (uint8_t)src[i];
  }
}

//Transform a matlab 2d array form to "normal" C form.
//Matlab has colums first in memory. Rows expected.
uint8_t* col2rowMO(uint8_t* input, int width, int height)
{
  uint8_t* output = malloc(sizeof(uint8_t)*width*height);
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      output[width*i + j] = input[height*j + i];
    }
  }

  return output;
}

//Switch from row major order to col major order
int* row2colMO(int* input, int width, int height)
{
  int* output = malloc(sizeof(int)*width*height);
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      output[height*i + j] = input[width*j + i];
    }
  }

  return output;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //Check correct input and output parameter number
  if (nlhs != 3) {
    ERROR(nlhs, "Three output required.");
  }

  if (nrhs != 6) {
    ERROR(nrhs, "Nine input required.");
  }

  //Test input data type
  if (!mxIsUint8(prhs[0]) || !mxIsUint8(prhs[2]) || !mxIsUint8(prhs[4])) {
    ERROR(notUint8, "Input data needs to be Uint8");
  }

  //Check correct input sizes
  if (mxGetN(prhs[1]) != 2 || mxGetM(prhs[1]) != 1 ||
      mxGetN(prhs[3]) != 2 || mxGetM(prhs[3]) != 1 ||
      mxGetN(prhs[5]) != 2 || mxGetM(prhs[5]) !=1 ) {
    ERROR(notSizeVector, "The input sizes need to be 1x2 vectors");
  }

  //Get Sizes
  int in_y_width = (int)mxGetN(prhs[0]);
  int in_y_height = (int)mxGetM(prhs[0]);
  int in_cb_width = (int)mxGetN(prhs[2]);
  int in_cb_height = (int)mxGetM(prhs[2]);
  int in_cr_width = (int)mxGetN(prhs[4]);
  int in_cr_height = (int)mxGetM(prhs[4]);

  int out_y_width = ((uint32_t*)mxGetData(prhs[1]))[1];
  int out_y_height = ((uint32_t*)mxGetData(prhs[1]))[0];
  int out_cb_width = ((uint32_t*)mxGetData(prhs[3]))[1];
  int out_cb_height = ((uint32_t*)mxGetData(prhs[3]))[0];
  int out_cr_width = ((uint32_t*)mxGetData(prhs[5]))[1];
  int out_cr_height = ((uint32_t*)mxGetData(prhs[5]))[0];

  uint8_t* y_data = (uint8_t*)mxGetData(prhs[0]);
  uint8_t* cb_data = (uint8_t*)mxGetData(prhs[2]);
  uint8_t* cr_data = (uint8_t*)mxGetData(prhs[4]);

  //Swithch major order
  y_data = col2rowMO(y_data, in_y_width, in_y_height);
  cb_data = col2rowMO(cb_data, in_cb_width, in_cb_height);
  cr_data = col2rowMO(cr_data, in_cr_width, in_cr_height);

  //Define output
  plhs[0] = mxCreateNumericMatrix(out_y_height, out_y_width, mxUINT8_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(out_cb_height, out_cb_width, mxUINT8_CLASS, mxREAL);
  plhs[2] = mxCreateNumericMatrix(out_cr_height, out_cr_width, mxUINT8_CLASS, mxREAL);

  //Do actual scaling stuff
  chroma_format_t is_420 = in_y_width != in_cb_width ? CHROMA_420 : CHROMA_444;
  yuv_buffer_t* pic = newYuvBuffer_uint8(y_data, cb_data, cr_data, in_y_width, in_y_height, is_420);
  scaling_parameter_t param = newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height, in_y_width != in_cb_width ? CHROMA_420 : CHROMA_444);

  yuv_buffer_t* scaled = yuvScaling(pic, &param, NULL);//scale(pic, &param, is_420);
  if (scaled == NULL) {
    ERROR(scalingError, "Scaling failed. Tried to up- and downscale at the same time or unable to allocate memory");
  }
  //ERROR(debug, "got past downscaling");
  free(y_data);
  free(cb_data);
  free(cr_data);

  //Go back to col major order
  y_data = row2colMO(scaled->y->data, out_y_width, out_y_height);
  cb_data = row2colMO(scaled->u->data, out_cb_width, out_cb_height);
  cr_data = row2colMO(scaled->v->data, out_cr_width, out_cr_height);

  //Copy results to output
  copyBack((uint8_t*)mxGetData(plhs[0]), y_data, out_y_height*out_y_width);
  copyBack((uint8_t*)mxGetData(plhs[1]), cb_data, out_cb_height*out_cb_width);
  copyBack((uint8_t*)mxGetData(plhs[2]), cr_data, out_cr_height*out_cr_width);

  //ERROR(debug, "Got to the end");

  //Free memory appears problematic. Need to use matlabs memory deallocation?
  free(y_data);
  free(cb_data);
  free(cr_data);
  deallocateYuvBuffer(pic);
  deallocateYuvBuffer(scaled);
}