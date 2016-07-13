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

  //Define output
  plhs[0] = mxCreateNumericMatrix(out_y_height, out_y_width, mxUINT8_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(out_cb_height, out_cb_width, mxUINT8_CLASS, mxREAL);
  plhs[2] = mxCreateNumericMatrix(out_cr_height, out_cr_width, mxUINT8_CLASS, mxREAL);

  //Do actual scaling stuff
  int is_420 = in_y_width != in_cb_width ? 1 : 0;
  yuv_buffer_t* pic = newYuvBuffer_uint8(y_data, cb_data, cr_data, in_y_width, in_y_height, is_420);
  scaling_parameter_t param = newScalingParameters(in_y_width, in_y_height, out_y_width, out_y_height);

  yuv_buffer_t* scaled = scale(pic, &param, is_420);
  
  //ERROR(debug, "got past downscaling");

  //Copy results to output
  copyBack((uint8_t*)mxGetData(plhs[0]), scaled->y->data, out_y_height*out_y_width);
  copyBack((uint8_t*)mxGetData(plhs[1]), scaled->u->data, out_cb_height*out_cb_width);
  copyBack((uint8_t*)mxGetData(plhs[2]), scaled->v->data, out_cr_height*out_cr_width);

  //ERROR(debug, "Got to the end");

  //Free memory appears problematic. Need to use matlabs memory deallocation?
  deallocateYuvBuffer(pic);
  deallocateYuvBuffer(scaled);
}