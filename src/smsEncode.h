#pragma once

#include <fftw3.h>
#include <omp.h>
#include <string.h>
#include <sys/unistd.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>
#include "ArrayTemplates.hpp"
#include "mri_data.h"

#ifndef PI
#define PI 3.14159265359
#endif

class smsEncode {
 public:
  // Kernel Types
  enum { TRIANGLE_KERNEL,
         KAISER_KERNEL,
         SINC_KERNEL,
         POLY_KERNEL };

  NDarray::Array<NDarray::Array<complex<float>, 3>, 2>
      k3d_grid; /*Actual Gridding Space*/

  // Controls for phase encode / 2D directions
  int fft_in_x;
  int fft_in_y;
  int grid_in_x;
  int grid_in_y;
  int grid_in_z;

  bool pruned_fft;

  // Grid Testing for Double
  int double_grid;

  /*Overgridding Factor*/
  float grid_x;
  float grid_y;

  /*Overgridding Crop values*/
  int og_sx;
  int og_sy;
  int og_ex;
  int og_ey;

  float grid_scale_x;
  float grid_scale_y;
  float grid_scale_z;

  float sms_slice_phase;
  int sms_factor;

  // arma::matrix Size
  int Nx;
  int Ny;
  int Nz;
  int Ne;
  int Nt;

  // Grid Size
  int Sx;
  int Sy;
  int Sz;
  int kernel_type;
  float overgrid;

  /*Kaiser Bessel Beta - Calculated*/
  float betaX;
  float betaY;
  float betaZ;

  // Discrete Gridding Kernel Variables
  NDarray::Array<float, 1> winx;
  NDarray::Array<float, 1> winy;
  float dwinX;
  float dwinY;
  float dwinZ;
  float sms_fwhm;
  float grid_modX;
  float grid_modY;
  float grid_modZ;
  NDarray::Array<float, 1> grid_filterX;
  NDarray::Array<float, 1> grid_filterY;
  NDarray::Array<float, 1> grid_filterZ;

  NDarray::Array<float, 1> z_position;

  // FFT
  fftwf_plan fft_plan;
  fftwf_plan ifft_plan;

  float k_rad;
  int time_grid;

  smsEncode();

  void alloc_grid();
  void read_commandline(int numarg, char **pstring);
  void precalc_gridding(int Nz, int Ny, int Nx, int Nt, int Ne, float sms_gap,
                        MRI_DATA &);
  void precalc_kernel(void);

  void do_fft(void);
  void do_ifft(void);

  // Main Calls with sensitivity map
  void forward(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &X,
               NDarray::Array<complex<float>, 3> &smap,
               NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &data,
               NDarray::Array<NDarray::Array<float, 3>, 2> &kx,
               NDarray::Array<NDarray::Array<float, 3>, 2> &ky,
               NDarray::Array<NDarray::Array<float, 3>, 2> &kz,
               NDarray::Array<NDarray::Array<float, 3>, 2> &kw,
               NDarray::Array<NDarray::Array<float, 3>, 3> &z);

  void backward(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &X,
                NDarray::Array<complex<float>, 3> &smap,
                NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &data,
                NDarray::Array<NDarray::Array<float, 3>, 2> &kx,
                NDarray::Array<NDarray::Array<float, 3>, 2> &ky,
                NDarray::Array<NDarray::Array<float, 3>, 2> &kz,
                NDarray::Array<NDarray::Array<float, 3>, 2> &kw,
                NDarray::Array<NDarray::Array<float, 3>, 3> &z);

  // For accumulating images
  void accumulate(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &X,
                  NDarray::Array<complex<float>, 3> &smap);

  // For setting images
  void set_image(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &X,
                 NDarray::Array<complex<float>, 3> &smap);

  void chop_grid_forward(
      NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &data,
      NDarray::Array<NDarray::Array<float, 3>, 2> &kx,
      NDarray::Array<NDarray::Array<float, 3>, 2> &ky,
      NDarray::Array<NDarray::Array<float, 3>, 2> &kz,
      NDarray::Array<NDarray::Array<float, 3>, 2> &kw,
      NDarray::Array<NDarray::Array<float, 3>, 3> &z);

  void chop_grid_backward(
      NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &data,
      NDarray::Array<NDarray::Array<float, 3>, 2> &kx,
      NDarray::Array<NDarray::Array<float, 3>, 2> &ky,
      NDarray::Array<NDarray::Array<float, 3>, 2> &kz,
      NDarray::Array<NDarray::Array<float, 3>, 2> &kw,
      NDarray::Array<NDarray::Array<float, 3>, 3> &z);

  static float bessi0(float);
  static void loadKernelTable(NDarray::Array<float, 1> &out);
  static void help_message(void);

 private:
};
