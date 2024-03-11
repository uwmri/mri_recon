#ifndef hDCFGRID
#define hDCFGRID

#include <fftw3.h>
#include <omp.h>
#include <string.h>
#include <sys/unistd.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "ArrayTemplates.hpp"
#include "mri_data.h"

#ifndef PI
#define PI 3.14159265359
#endif

class DCFgridFFT {
 public:
  // Kernel Types
  enum { POLY_KERNEL,
         POLY_KERNEL2 };

  // Controls for phase encode / 2D directions
  int grid_in_x;
  int grid_in_y;
  int grid_in_z;

  /*Overgridding Factor*/
  float grid_x;
  float grid_y;
  float grid_z;

  float grid_scale_x;
  float grid_scale_y;
  float grid_scale_z;

  // arma::matrix Size
  int Nx;
  int Ny;
  int Nz;

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
  float dwinX;
  float dwinY;
  float dwinZ;
  float grid_modX;
  float grid_modY;
  float grid_modZ;
  NDarray::Array<float, 1> grid_filterX;
  NDarray::Array<float, 1> grid_filterY;
  NDarray::Array<float, 1> grid_filterZ;

  int time_grid;

  float acc;

  DCFgridFFT();

  void precalc_kernel(void);
  void loadKernelTable(NDarray::Array<float, 1>& out);
  void loadKernelTable2(NDarray::Array<float, 1>& out);

  // For setting images
  void forward(NDarray::Array<float, 3>& image,
               const NDarray::Array<float, 3>& data,
               const NDarray::Array<float, 3>& kx,
               const NDarray::Array<float, 3>& ky,
               const NDarray::Array<float, 3>& kz);
  void backward(NDarray::Array<float, 3>& image, NDarray::Array<float, 3>& data,
                const NDarray::Array<float, 3>& kx,
                const NDarray::Array<float, 3>& ky,
                const NDarray::Array<float, 3>& kz);

  void scale_kw(NDarray::Array<float, 3>& data,
                const NDarray::Array<float, 3>& kx,
                const NDarray::Array<float, 3>& ky,
                const NDarray::Array<float, 3>& kz);

 private:
};

#endif