#ifndef hGRIDFFTCOILTHREADED
#define hGRIDFFTCOILTHREADED

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

typedef NDarray::Array<complex<float>, 3> Complex3D;
typedef NDarray::Array<NDarray::Array<complex<float>, 3>, 1> Complex4D_GRID;
typedef NDarray::Array<NDarray::Array<complex<float>, 3>, 1> Complex4D;
typedef NDarray::Array<NDarray::Array<complex<float>, 3>, 2> Complex5D;

class gridFFT_CoilThreaded {
 public:
  // Kernel Types
  enum { TRIANGLE_KERNEL,
         KAISER_KERNEL,
         SINC_KERNEL,
         POLY_KERNEL };

  Complex4D_GRID k3d_grid; /*Actual Gridding Space*/
  Complex4D_GRID image;    /*Complex Storage Space*/

  // Controls for phase encode / 2D directions
  int fft_in_x;
  int fft_in_y;
  int fft_in_z;
  int grid_in_x;
  int grid_in_y;
  int grid_in_z;

  bool pruned_fft;
  bool fast_fft_plan;

  // Grid Testing for Double
  int double_grid;

  /*Overgridding Factor*/
  float grid_x;
  float grid_y;
  float grid_z;

  /*Overgridding Crop values*/
  int og_sx;
  int og_sy;
  int og_sz;
  int og_ex;
  int og_ey;
  int og_ez;

  float grid_scale_x;
  float grid_scale_y;
  float grid_scale_z;

  bool sms_flag;
  float sms_factor;

  // arma::matrix Size
  int Nx;
  int Ny;
  int Nz;
  int Num_Coils;

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
  NDarray::Array<float, 1> winz;
  float dwinX;
  float dwinY;
  float dwinZ;
  float grid_modX;
  float grid_modY;
  float grid_modZ;
  NDarray::Array<float, 1> grid_filterX;
  NDarray::Array<float, 1> grid_filterY;
  NDarray::Array<float, 1> grid_filterZ;

  // FFT
  fftwf_plan fft_plan;
  fftwf_plan ifft_plan;

  float k_rad;
  int time_grid;

  gridFFT_CoilThreaded();

  void alloc_grid();
  void read_commandline(int numarg, const char **pstring);
  void precalc_gridding(int Nz, int Ny, int Nx, MRI_DATA &);
  void precalc_kernel(void);

  void do_fft(void);
  void do_ifft(void);
  void deapp(void);

  void loadKernelTable(NDarray::Array<float, 1> &out);

  // Main Calls with sensitivity map
  void forward(NDarray::Array<complex<float>, 3> &X, const Complex4D &smap,
               const Complex4D &data, const NDarray::Array<float, 3> &kx,
               const NDarray::Array<float, 3> &ky,
               const NDarray::Array<float, 3> &kz,
               const NDarray::Array<float, 3> &kw);

  void backward(const NDarray::Array<complex<float>, 3> &X,
                const Complex4D &smap, Complex4D &data,
                const NDarray::Array<float, 3> &kx,
                const NDarray::Array<float, 3> &ky,
                const NDarray::Array<float, 3> &kz,
                const NDarray::Array<float, 3> &kw);

  void backward_residual(const NDarray::Array<complex<float>, 3> &X,
                         const Complex4D &smap, Complex4D &data,
                         const NDarray::Array<float, 3> &kx,
                         const NDarray::Array<float, 3> &ky,
                         const NDarray::Array<float, 3> &kz,
                         const NDarray::Array<float, 3> &kw,
                         const Complex4D &dataSub);

  // Main Calls without sensitivity map
  void forward(NDarray::Array<complex<float>, 3> &X, const Complex4D &data,
               const NDarray::Array<float, 3> &kx,
               const NDarray::Array<float, 3> &ky,
               const NDarray::Array<float, 3> &kz,
               const NDarray::Array<float, 3> &kw);

  void backward(const NDarray::Array<complex<float>, 3> &X, Complex4D &data,
                const NDarray::Array<float, 3> &kx,
                const NDarray::Array<float, 3> &ky,
                const NDarray::Array<float, 3> &kz,
                const NDarray::Array<float, 3> &kw);

  // Special call for sum of squares
  void forward_sos(NDarray::Array<complex<float>, 3> &X, const Complex4D &data,
                   const NDarray::Array<float, 3> &kx,
                   const NDarray::Array<float, 3> &ky,
                   const NDarray::Array<float, 3> &kz,
                   const NDarray::Array<float, 3> &kw);

  // For accumulating images
  void accumulate(NDarray::Array<complex<float>, 3> &X, const Complex4D &smap);
  void accumulate(NDarray::Array<complex<float>, 3> &X);
  void accumulate_sos(NDarray::Array<complex<float>, 3> &X);

  // For setting images
  void set_image(const NDarray::Array<complex<float>, 3> &X,
                 const Complex4D &smap);
  void set_image(const NDarray::Array<complex<float>, 3> &X);
  void zero(void);

  void chop_grid_forward(const Complex4D &data,
                         const NDarray::Array<float, 3> &kx,
                         const NDarray::Array<float, 3> &ky,
                         const NDarray::Array<float, 3> &kz,
                         const NDarray::Array<float, 3> &kw);
  void chop_grid_backward(Complex4D &data, const NDarray::Array<float, 3> &kx,
                          const NDarray::Array<float, 3> &ky,
                          const NDarray::Array<float, 3> &kz,
                          const NDarray::Array<float, 3> &kw,
                          const Complex4D &diff_data, bool);
  static float bessi0(float);

  static void help_message(void);

 private:
};

#endif
