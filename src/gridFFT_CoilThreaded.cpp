/************************************************
Gridding and FFT Libraries for K-space to Image Domain Tranformation

Initial Author:
        Kevin M. Johnson

Usage Example:
        // Do Once
        gridFFT_CoilThreaded gridding;
        gridding.read_commandline(argc,argv);
        gridding.precalc_gridding(64,64,64,48,3,KAISER);

        // Use for transform
        gridding.forward()


*************************************************/

#include "gridFFT_CoilThreaded.h"
#include "io_templates.hpp"
#include "tictoc.hpp"
using namespace NDarray;

//----------------------------------------
// Constructor - Sets Default Vals
//----------------------------------------
gridFFT_CoilThreaded::gridFFT_CoilThreaded() {
  grid_x = -1;
  grid_y = -1;
  grid_z = -1;
  overgrid = 1.375;
  kernel_type = KAISER_KERNEL;
  betaX = 12;
  betaY = 12;
  betaZ = 12;
  dwinX = -1;
  dwinY = -1;
  dwinZ = -1;
  fft_plan = NULL;
  ifft_plan = NULL;
  fft_in_x = -1;
  fft_in_y = -1;
  fft_in_z = -1;
  grid_in_x = 1;
  grid_in_y = 1;
  grid_in_z = -1;
  grid_scale_x = 1.0;
  grid_scale_y = 1.0;
  grid_scale_z = 1.0;
  k_rad = 9999.0;
  time_grid = 0;
  double_grid = 0;
  Num_Coils = 1;
  pruned_fft = true;

  sms_flag = 0;
  sms_factor = 2;
  fast_fft_plan = false;
}

//----------------------------------------
// Allocate Memory for Gridding
//----------------------------------------

void gridFFT_CoilThreaded::alloc_grid() {
  // Alloc contain used to grid data to
  Array<Array<complex<float>, 3>, 1> temp =
      Alloc4DContainer<complex<float> >(Sx, Sy, Sz, Num_Coils);
  k3d_grid.reference(temp);

  // Alloc a reference to central matrix
  image.resize(Num_Coils);
  for (int coil = 0; coil < Num_Coils; coil++) {
    Array<complex<float>, 3> image2 =
        k3d_grid(coil)(Range(og_sx, og_ex - 1), Range(og_sy, og_ey - 1),
                       Range(og_sz, og_ez - 1));
    image(coil).reference(image2);
  }
}

// ----------------------
// Help Message
// ----------------------
void gridFFT_CoilThreaded::help_message(void) {
  cout << "----------------------------------------------" << endl;
  cout << "   Gridding Control " << endl;
  cout << "----------------------------------------------" << endl;

  cout << "Control" << endl;
  help_flag("-kaiser", "use kaiser bessel kernel");
  help_flag("-triangle", "use triangle kernel");
  help_flag("-sinc", "use sinc kernel");
  help_flag("-overgrid []", "overgrid by factor []");
  help_flag("-fast_grid", "no overgrid,traingle kernel");
  help_flag("-time_grid", "output times for gridding");

  cout << "Fine Control" << endl;
  help_flag("-grid_in_x []", "0=use nearest neighbor interpolation in x");
  help_flag("-grid_in_y []", "0=use nearest neighbor interpolation in y");
  help_flag("-grid_in_z []", "0=use nearest neighbor interpolation in z");
  help_flag("-fft_in_x []", "0=no fft in x");
  help_flag("-fft_in_y []", "0=no fft in y");
  help_flag("-fft_in_z []", "0=no fft in z");
  help_flag("-dwinX []", "Size of kernel in x");
  help_flag("-dwinY []", "Size of kernel in y");
  help_flag("-dwinZ []", "Size of kernel in z");
  help_flag("-double_grid", "Use complex<double> for forward gridding (slow)");
  help_flag("-pruned_fftf", "Used pruned FFT (experimental)");
}

//----------------------------------------
// Parse Command Line Args
//----------------------------------------

void gridFFT_CoilThreaded::read_commandline(int numarg, const char **pstring) {
#define trig_flag(num, name, val)             \
  }                                           \
  else if (strcmp(name, pstring[pos]) == 0) { \
    val = num;
#define float_flag(name, val)                 \
  }                                           \
  else if (strcmp(name, pstring[pos]) == 0) { \
    pos++;                                    \
    val = atof(pstring[pos]);
#define int_flag(name, val)                   \
  }                                           \
  else if (strcmp(name, pstring[pos]) == 0) { \
    pos++;                                    \
    val = atoi(pstring[pos]);

  for (int pos = 0; pos < numarg; pos++) {
    if (strcmp("-h", pstring[pos]) == 0) {
      float_flag("-overgrid", overgrid);
      float_flag("-dwinX", dwinX);
      float_flag("-dwinY", dwinY);
      float_flag("-dwinZ", dwinZ);
      float_flag("-grid_x", grid_x);
      float_flag("-grid_y", grid_y);
      float_flag("-grid_z", grid_z);

      float_flag("-grid_scale_x", grid_scale_x);
      float_flag("-grid_scale_y", grid_scale_y);
      float_flag("-grid_scale_z", grid_scale_z);

      int_flag("-grid_in_x", grid_in_x);
      int_flag("-grid_in_y", grid_in_y);
      int_flag("-grid_in_z", grid_in_z);
      int_flag("-fft_in_x", fft_in_x);
      int_flag("-fft_in_y", fft_in_y);
      int_flag("-fft_in_z", fft_in_z);

      trig_flag(true, "-sms_flag", sms_flag);
      int_flag("-sms_factor", sms_factor);

      trig_flag(KAISER_KERNEL, "-kaiser", kernel_type);
      trig_flag(TRIANGLE_KERNEL, "-triangle", kernel_type);
      trig_flag(SINC_KERNEL, "-sinc", kernel_type);
      trig_flag(POLY_KERNEL, "-poly_kernel", kernel_type);

      trig_flag(1, "-time_grid", time_grid);
      trig_flag(1, "-double_grid", double_grid);

      trig_flag(true, "-pruned_fft", pruned_fft);
      // Special Copies
    } else if (strcmp("-fast_grid", pstring[pos]) == 0) {
      overgrid = 1.0;
      kernel_type = TRIANGLE_KERNEL;
      dwinX = 1.0;
      dwinY = 1.0;
      dwinZ = 1.0;
    }
  }
}

/* The kernel's radius FOV product is the length,
 * in pixels, to the first truncation point.
 */

void gridFFT_CoilThreaded::loadKernelTable(Array<float, 1> &out) {
  int len = out.numElements();
  double c0 = 1.;
  double c1 = 0.04522831;
  double c2 = -3.36020304;
  double c3 = 1.12417012;
  double c4 = 2.82448025;
  double c5 = -1.63447764;

  for (int i = 0; i < len; i++) {
    double x = double(i) / double(len);
    double x2 = x * x;
    double x4 = x2 * x2;
    double x3 = x * x2;
    double x5 = x * x4;
    out(i) = (float)(c0 + c1 * x + c2 * x2 + c3 * x3 + c4 * x4 + c5 * x5);
  }

  return;
}

//----------------------------------------
//    Setup for Gridding
//----------------------------------------

void gridFFT_CoilThreaded::precalc_kernel(void) {
  // How many pts in kernel per delta k for kernel lookup table
  grid_modX = 600;
  grid_modY = 600;
  grid_modZ = 600;

  // ------------------------------------------------
  //    Kernel Calculations
  // ------------------------------------------------

  switch (kernel_type) {
    case (TRIANGLE_KERNEL): {
      // Kernel Half Size
      dwinX = (dwinX == -1) ? (1.0) : (dwinX);
      dwinY = (dwinY == -1) ? (1.0) : (dwinY);
      dwinZ = (dwinZ == -1) ? (1.0) : (dwinZ);

      // Grid Length for precomputed kernel
      int grid_lengthX = (int)((float)dwinX * (float)grid_modX);
      int grid_lengthY = (int)((float)dwinY * (float)grid_modY);
      int grid_lengthZ = (int)((float)dwinZ * (float)grid_modZ);

      // Alloc Lookup Table Structs for Gridding
      grid_filterX.resize(grid_lengthX + 10);
      grid_filterY.resize(grid_lengthY + 10);
      grid_filterZ.resize(grid_lengthZ + 10);
      grid_filterX = 0.0;
      grid_filterY = 0.0;
      grid_filterZ = 0.0;

      // Compute Seperable Kernel
      for (int i = 0; i < (grid_lengthX + 1); i++) {
        float grid_pos = (float)i / (float)grid_lengthX;
        grid_filterX(i) = 1.0 - grid_pos;
      }

      for (int i = 0; i < (grid_lengthY + 1); i++) {
        float grid_pos = (float)i / (float)grid_lengthY;
        grid_filterY(i) = 1.0 - grid_pos;
      }

      for (int i = 0; i < (grid_lengthZ + 1); i++) {
        float grid_pos = (float)i / (float)grid_lengthZ;
        grid_filterZ(i) = 1.0 - grid_pos;
      }
    } break;

    case (KAISER_KERNEL): {
      // Kernel Half Size
      dwinX = (dwinX == -1) ? (2.5) : (dwinX);
      dwinY = (dwinY == -1) ? (2.5) : (dwinY);
      dwinZ = (dwinZ == -1) ? (2.5) : (dwinZ);

      // Grid Length for precomputed kernel
      int grid_lengthX = (int)((float)dwinX * (float)grid_modX);
      int grid_lengthY = (int)((float)dwinY * (float)grid_modY);
      int grid_lengthZ = (int)((float)dwinZ * (float)grid_modZ);

      // Alloc Lookup Table Structs for Gridding
      grid_filterX.resize(grid_lengthX + 10);
      grid_filterY.resize(grid_lengthY + 10);
      grid_filterZ.resize(grid_lengthZ + 10);
      grid_filterX = 0.0;
      grid_filterY = 0.0;
      grid_filterZ = 0.0;

      // Get optimal Beta per Beatty et al
      float act_grid_x = grid_x * grid_scale_x;
      float act_grid_y = grid_y * grid_scale_y;
      float act_grid_z = grid_z * grid_scale_z;

      betaX = PI * sqrtf((dwinX * dwinX) / (act_grid_x * act_grid_x) *
                             (act_grid_x - 0.5) * (act_grid_x - 0.5) -
                         0.8);
      betaY = PI * sqrtf((dwinY * dwinY) / (act_grid_y * act_grid_y) *
                             (act_grid_y - 0.5) * (act_grid_y - 0.5) -
                         0.8);
      betaZ = PI * sqrtf((dwinZ * dwinZ) / (act_grid_z * act_grid_z) *
                             (act_grid_z - 0.5) * (act_grid_z - 0.5) -
                         0.8);

      float beta_minX = sqrt(pow(PI * 2 * dwinX / act_grid_x, 2.0) - 9.6752);
      float beta_minY = sqrt(pow(PI * 2 * dwinY / act_grid_y, 2.0) - 9.6752);
      float beta_minZ = sqrt(pow(PI * 2 * dwinZ / act_grid_z, 2.0) - 9.6752);
      betaX = (beta_minX > betaX) ? (beta_minX) : (betaX);
      betaY = (beta_minY > betaY) ? (beta_minY) : (betaY);
      betaZ = (beta_minZ > betaZ) ? (beta_minZ) : (betaZ);

      // Compute Seperable Kernels
      for (int i = 0; i < (grid_lengthX + 1); i++) {
        float grid_pos = ((float)i) / ((float)grid_lengthX);
        float grid_arg = sqrtf(1.0 - grid_pos * grid_pos);
        grid_filterX(i) = bessi0(betaX * grid_arg) / bessi0(betaX);
      }

      for (int i = 0; i < (grid_lengthY + 1); i++) {
        float grid_pos = ((float)i) / ((float)grid_lengthY);
        float grid_arg = sqrtf(1.0 - grid_pos * grid_pos);
        grid_filterY(i) = bessi0(betaY * grid_arg) / bessi0(betaY);
      }

      for (int i = 0; i < (grid_lengthZ + 1); i++) {
        float grid_pos = ((float)i) / ((float)grid_lengthZ);
        float grid_arg = sqrtf(1.0 - grid_pos * grid_pos);
        grid_filterZ(i) = bessi0(betaZ * grid_arg) / bessi0(betaZ);
      }

    } break;

    case (SINC_KERNEL): {
      // Kernel Half Size
      dwinX = (dwinX == -1) ? (2.5) : (dwinX);
      dwinY = (dwinY == -1) ? (2.5) : (dwinY);
      dwinZ = (dwinZ == -1) ? (2.5) : (dwinZ);

      // Grid Length for precomputed kernel
      int grid_lengthX = (int)((float)dwinX * (float)grid_modX);
      int grid_lengthY = (int)((float)dwinY * (float)grid_modY);
      int grid_lengthZ = (int)((float)dwinZ * (float)grid_modZ);

      // Alloc Lookup Table Structs for Gridding
      grid_filterX.resize(grid_lengthX + 10);
      grid_filterY.resize(grid_lengthY + 10);
      grid_filterZ.resize(grid_lengthZ + 10);
      grid_filterX = 0.0;
      grid_filterY = 0.0;
      grid_filterZ = 0.0;

      // Compute Seperable Kernels
      for (int i = 0; i < (grid_lengthX + 1); i++) {
        float grid_pos = ((float)i * (float)dwinX / (float)grid_scale_x) /
                         ((float)grid_lengthX);
        if (grid_pos < 0.01) {
          grid_filterX(i) = 1.0;
        } else {
          grid_filterX(i) = (float)(sin(grid_pos * PI) / (grid_pos * PI));
        }
        grid_filterX(i) *=
            (0.54 + 0.46 * cos(PI * (float)i / ((float)grid_lengthZ)));
      }

      for (int i = 0; i < (grid_lengthY + 1); i++) {
        float grid_pos = ((float)i * (float)dwinY / (float)grid_scale_y) /
                         ((float)grid_lengthY);
        if (grid_pos < 0.01) {
          grid_filterY(i) = 1.0;
        } else {
          grid_filterY(i) = (float)(sin(grid_pos * PI) / (grid_pos * PI));
        }
        grid_filterY(i) *=
            (0.54 + 0.46 * cos(PI * (float)i / ((float)grid_lengthY)));
      }

      for (int i = 0; i < (grid_lengthZ + 1); i++) {
        float grid_pos = ((float)i * (float)dwinZ / (float)grid_scale_z) /
                         ((float)grid_lengthZ);
        if (grid_pos < 0.01) {
          grid_filterZ(i) = 1.0;
        } else {
          grid_filterZ(i) = (float)(sin(grid_pos * PI) / (grid_pos * PI));
        }
        grid_filterZ(i) *=
            (0.54 + 0.46 * cos(PI * (float)i / ((float)grid_lengthZ)));
      }

    } break;

    case (POLY_KERNEL): {
      // Kernel Half Size
      dwinX = (dwinX == -1) ? (2) : (dwinX);
      dwinY = (dwinY == -1) ? (2) : (dwinY);
      dwinZ = (dwinZ == -1) ? (2) : (dwinZ);

      // Grid Length for precomputed kernel
      int grid_lengthX = (int)((float)dwinX * (float)grid_modX);
      int grid_lengthY = (int)((float)dwinY * (float)grid_modY);
      int grid_lengthZ = (int)((float)dwinZ * (float)grid_modZ);

      // Alloc Lookup Table Structs for Gridding
      grid_filterX.resize(grid_lengthX + 10);
      grid_filterY.resize(grid_lengthY + 10);
      grid_filterZ.resize(grid_lengthZ + 10);
      grid_filterX = 0.0;
      grid_filterY = 0.0;
      grid_filterZ = 0.0;

      loadKernelTable(grid_filterX);
      loadKernelTable(grid_filterY);
      loadKernelTable(grid_filterZ);

    } break;
  }

  // Normalize
  grid_filterX *= 0.5 * grid_modX / sum(grid_filterX);
  grid_filterY *= 0.5 * grid_modY / sum(grid_filterY);
  grid_filterZ *= 0.5 * grid_modZ / sum(grid_filterZ);
}

void gridFFT_CoilThreaded::precalc_gridding(int NzT, int NyT, int NxT,
                                            MRI_DATA &data) {
  this->Nx = NxT;
  this->Ny = NyT;
  this->Nz = NzT;
  this->Num_Coils = data.Num_Coils;

  // ---------------------------------------------
  // Determine what needs to be grid
  // ---------------------------------------------

  grid_in_x = (data.trajectory_type(0) == MRI_DATA::NONCARTESIAN) ? (1) : (0);
  grid_in_y = (data.trajectory_type(1) == MRI_DATA::NONCARTESIAN) ? (1) : (0);
  grid_in_z = (data.trajectory_type(2) == MRI_DATA::NONCARTESIAN) ? (1) : (0);

  if (fft_in_x == -1) {
    fft_in_x = data.dft_needed(0);
  }

  if (fft_in_y == -1) {
    fft_in_y = data.dft_needed(1);
  }

  if (fft_in_z == -1) {
    fft_in_z = data.dft_needed(2);
  }

  // Get rounded Gridding ratio*
  if (grid_in_x == 1) {
    if (grid_x == -1) {
      grid_x = 16.0 * ceil((overgrid * (float)Nx) / 16.0) / (float)Nx;
    }
  } else {
    grid_x = 1;
  }

  if (grid_in_y == 1) {
    if (grid_y == -1) {
      grid_y = 16.0 * ceil((overgrid * (float)Ny) / 16.0) / (float)Ny;
    }
  } else {
    grid_y = 1;
  }

  if (grid_in_z == 1) {
    if (grid_z == -1) {
      grid_z = 16.0 * ceil((overgrid * (float)Nz) / 16.0) / (float)Nz;
    }
  } else {
    grid_z = 1;
  }
  // Compute Grid Size
  Sz = (int)(grid_z * Nz);
  Sy = (int)(grid_y * Ny);
  Sx = (int)(grid_x * Nx);

  precalc_kernel();

  // ------------------------------------------------
  //    Image Domain Calcs (crop + deapp)
  // ------------------------------------------------

  // Calculations to Determine Crop Positions
  og_sx = (int)((float)Nx * (grid_x - 1) / 2.0);
  og_sy = (int)((float)Ny * (grid_y - 1) / 2.0);
  og_sz = (int)((float)Nz * (grid_z - 1) / 2.0);
  og_ex = og_sx + Nx;
  og_ey = og_sy + Ny;
  og_ez = og_sz + Nz;

  printf("\n\nCoil (coil = %d) Threaded Gridding Kernel Info\n",
         this->Num_Coils);
  printf("Dwin 		%f %f %f\n", dwinX, dwinY, dwinZ);
  printf("Mod 		%f %f %f\n", grid_modX, grid_modY, grid_modZ);
  printf("Og %d-%d x %d-%d x %d-%d\n", og_sx, og_ex, og_sy, og_ey, og_sz,
         og_ez);
  printf("Grid in x=%d, y=%d, z=%d\n", grid_in_x, grid_in_y, grid_in_z);
  printf("FFT in x=%d, y=%d, z=%d\n", fft_in_x, fft_in_y, fft_in_z);

  // Deapp Windows
  winx.resize(Sx);
  if ((fft_in_x == 1) && (grid_in_x == 1)) {
    winx = 0.0;
    for (int i = 0; i < Sx; i++) {
      float ipos = i - (float)Sx / 2.0;
      for (int grid_pos = 0; grid_pos < dwinX * grid_modX; grid_pos++) {
        // Fourier Transform of Kernel
        winx(i) +=
            2 * cos(2 * PI * ipos * grid_pos / (float)grid_modX / (float)Sx) *
            grid_filterX(grid_pos);
      }
      winx(i) = (float)grid_modX / winx(i);

      // Put chopping + zeroing in window to save time
      float fact = ((float)(2 * ((i) % 2) - 1));
      winx(i) *= fact / Sx;
      winx(i) *= (i < og_sx) ? (0.0) : (1.0);
      winx(i) *= (i > og_ex - 1) ? (0.0) : (1.0);
    }
    winx /= max(winx);
  } else if (fft_in_x) {
    for (int i = 0; i < Sx; i++) {
      float fact = ((float)(2 * ((i) % 2) - 1));
      winx(i) = fact / Sx;
    }
  } else {
    winx = 1.0;
  }

  winy.resize(Sy);
  if ((fft_in_y == 1) && (grid_in_y == 1)) {
    winy = 0.0;
    for (int i = 0; i < Sy; i++) {
      float ipos = i - (float)Sy / 2.0;
      for (int grid_pos = 0; grid_pos < dwinY * grid_modY; grid_pos++) {
        winy(i) +=
            2 * cos(2 * PI * ipos * grid_pos / (float)grid_modY / (float)Sy) *
            grid_filterY(grid_pos);
      }

      winy(i) = (float)grid_modY / winy(i);
      float fact = ((float)(2 * ((i) % 2) - 1));
      winy(i) *= fact / Sy;
      winy(i) *= (i < og_sy) ? (0.0) : (1.0);
      winy(i) *= (i > og_ey - 1) ? (0.0) : (1.0);
    }
    winy /= max(winy);
  } else if (fft_in_y) {
    for (int i = 0; i < Sy; i++) {
      float fact = ((float)(2 * ((i) % 2) - 1));
      winy(i) = fact / Sy;
    }
  } else {
    winy = 1.0;
  }

  winz.resize(Sz);
  if ((fft_in_z == 1) && (grid_in_z == 1)) {
    winz = 0.0;
    for (int i = 0; i < Sz; i++) {
      float ipos = i - (float)Sz / 2.0;
      for (int grid_pos = 0; grid_pos < dwinZ * grid_modZ; grid_pos++) {
        winz(i) +=
            2 * cos(2 * PI * ipos * grid_pos / (float)grid_modZ / (float)Sz) *
            grid_filterZ(grid_pos);
      }
      winz(i) = (float)grid_modZ / winz(i);
      float fact = ((float)(2 * ((i) % 2) - 1));
      winz(i) *= fact / Sz;
      winz(i) *= (i < og_sz) ? (0.0) : (1.0);
      winz(i) *= (i > og_ez - 1) ? (0.0) : (1.0);
    }
    winz /= max(winz);
  } else if (fft_in_z) {
    for (int i = 0; i < Sz; i++) {
      float fact = ((float)(2 * ((i) % 2) - 1));
      winz(i) = fact / Sz;
    }
  } else {
    winz = 1.0;
  }

  // Allocate Memory
  cout << "Alloc Grid" << endl;
  alloc_grid();
}

void gridFFT_CoilThreaded::do_fft(void) {
  // Smallest array
  if (fft_in_z == 1) {
    tictoc T;
    if (time_grid) T.tic();

    // Create Plan
    fftwf_plan plan;
#pragma omp critical
    {
      // Get a plan but never use
      fftwf_complex *ptr =
          (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * Sz);
      plan = fftwf_plan_dft_1d(Sz, ptr, ptr, FFTW_FORWARD, FFTW_MEASURE);
      fftwf_free(ptr);
    }

    // nested
    int DIMS[2];
    DIMS[0] = Ny;
    DIMS[1] = this->Num_Coils;
    long Nlines = DIMS[0] * DIMS[1];

#pragma omp parallel for
    for (long index = 0; index < Nlines; index++) {
      // Get the actual position
      int stride[2];
      nested_workaround(index, DIMS, stride, 2);
      int j = stride[0] + og_sy;
      int coil = stride[1];

      // Storage block
      Array<complex<float>, 2> TEMP(Sz, Nx, ColumnMajorArray<2>());

      // Gather (i is continous)
      for (int k = 0; k < Sz; k++) {
        for (int i = og_sx; i < og_ex; i++) {
          TEMP(k, i - og_sx) = k3d_grid(coil)(i, j, k);
        }
      }

      // FFT ( k is continous in this array but its small)
      for (int i = 0; i < Nx; i++) {
        fftwf_complex *data = reinterpret_cast<fftwf_complex *>(&TEMP(0, i));
        fftwf_execute_dft(plan, data, data);
      }

      // Scatter
      for (int k = 0; k < Sz; k++) {
        for (int i = og_sx; i < og_ex; i++) {
          k3d_grid(coil)(i, j, k) = TEMP(k, i - og_sx);
        }
      }
    }

    // Cleanup
    fftwf_destroy_plan(plan);
    if (time_grid) cout << "[FFTz " << T << "]";
  }

  // Bigger
  if (fft_in_y == 1) {
    tictoc T;
    if (time_grid) T.tic();

    // Create Plan
    fftwf_plan plan;
#pragma omp critical
    {
      // Get a plan but never use
      fftwf_complex *ptr =
          (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * Sy);
      plan = fftwf_plan_dft_1d(Sy, ptr, ptr, FFTW_FORWARD, FFTW_MEASURE);
      fftwf_free(ptr);
    }

    // nested
    int DIMS[2];
    DIMS[0] = Sz;
    DIMS[1] = this->Num_Coils;
    long Nlines = DIMS[0] * DIMS[1];

#pragma omp parallel for
    for (long index = 0; index < Nlines; index++) {
      // Get the actual position
      int stride[2];
      nested_workaround(index, DIMS, stride, 2);
      int k = stride[0];
      int coil = stride[1];

      // Copy data
      fftwf_complex *data_ptr =
          (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * Sy);
      complex<float> *data = reinterpret_cast<complex<float> *>(data_ptr);
      for (int i = og_sx; i < og_ex; i++) {
        for (int j = 0; j < Sy; j++) {
          data[j] = k3d_grid(coil)(i, j, k);
        }

        // FFT
        fftwf_execute_dft(plan, data_ptr, data_ptr);

        // Copy back
        for (int j = 0; j < Sy; j++) {
          k3d_grid(coil)(i, j, k) = data[j];
        }
      }
      fftwf_free(data_ptr);
    }

    // Cleanup
    fftwf_destroy_plan(plan);
    if (time_grid) cout << "[FFTy " << T << "]";
  }

  // Biggest
  if (fft_in_x == 1) {
    tictoc T;
    if (time_grid) T.tic();

    // Create Plan
    fftwf_plan plan;
#pragma omp critical
    {
      // Get a plan but never use
      fftwf_complex *ptr =
          (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * Sx);
      plan = fftwf_plan_dft_1d(Sx, ptr, ptr, FFTW_FORWARD, FFTW_MEASURE);
      fftwf_free(ptr);
    }

    // nested
    int DIMS[3];
    DIMS[0] = Sy;
    DIMS[1] = Sz;
    DIMS[2] = this->Num_Coils;
    long Nlines = DIMS[0] * DIMS[1] * DIMS[2];

#pragma omp parallel for
    for (long index = 0; index < Nlines; index++) {
      // Get the actual position
      int stride[3];
      nested_workaround(index, DIMS, stride, 3);
      int j = stride[0];
      int k = stride[1];
      int coil = stride[2];

      // i is contigous, just grab pointer and FFT
      fftwf_complex *data_ptr =
          reinterpret_cast<fftwf_complex *>(&k3d_grid(coil)(0, j, k));
      fftwf_execute_dft(plan, data_ptr, data_ptr);
    }

    // Cleanup
    fftwf_destroy_plan(plan);

    if (time_grid) cout << "[FFTx " << T << "]";
  }
}

void gridFFT_CoilThreaded::do_ifft(void) {
  // Biggest
  if (fft_in_x == 1) {
    tictoc T;
    if (time_grid) T.tic();

    int N = Sx;

    // Create Plan
    fftwf_plan plan;
#pragma omp critical
    {
      // Get a plan but never use
      fftwf_complex *ptr =
          (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * N);
      plan = fftwf_plan_dft_1d(N, ptr, ptr, FFTW_BACKWARD, FFTW_MEASURE);
      fftwf_free(ptr);
    }

    // nested
    int DIMS[3];
    DIMS[0] = Sy;
    DIMS[1] = Sz;
    DIMS[2] = this->Num_Coils;
    long Nlines = DIMS[0] * DIMS[1] * DIMS[2];

#pragma omp parallel for
    for (long index = 0; index < Nlines; index++) {
      // Get the actual position
      int stride[3];
      nested_workaround(index, DIMS, stride, 3);
      int j = stride[0];
      int k = stride[1];
      int coil = stride[2];

      // i is contigous, just grab pointer and FFT
      fftwf_complex *data_ptr =
          reinterpret_cast<fftwf_complex *>(&k3d_grid(coil)(0, j, k));
      fftwf_execute_dft(plan, data_ptr, data_ptr);
    }

    // Cleanup
    fftwf_destroy_plan(plan);

    if (time_grid) cout << "[FFTx " << T << "]";
  }

  // Bigger
  if (fft_in_y == 1) {
    tictoc T;
    if (time_grid) T.tic();

    int N = Sy;

    // Create Plan
    fftwf_plan plan;
#pragma omp critical
    {
      // Get a plan but never use
      fftwf_complex *ptr =
          (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * N);
      plan = fftwf_plan_dft_1d(N, ptr, ptr, FFTW_BACKWARD, FFTW_MEASURE);
      fftwf_free(ptr);
    }

    // nested
    int DIMS[2];
    DIMS[0] = Sz;
    DIMS[1] = this->Num_Coils;
    long Nlines = DIMS[0] * DIMS[1];

#pragma omp parallel for
    for (long index = 0; index < Nlines; index++) {
      // Get the actual position
      int stride[2];
      nested_workaround(index, DIMS, stride, 2);
      int kk = stride[0];
      int coil = stride[1];

      // Copy data
      fftwf_complex *data_ptr =
          (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * N);
      complex<float> *data = reinterpret_cast<complex<float> *>(data_ptr);
      for (int i = og_sx; i < og_ex; i++) {
        for (int j = 0; j < N; j++) {
          data[j] = k3d_grid(coil)(i, j, kk);
        }

        // FFT
        fftwf_execute_dft(plan, data_ptr, data_ptr);

        // Copy back
        for (int j = 0; j < N; j++) {
          k3d_grid(coil)(i, j, kk) = data[j];
        }
      }
      fftwf_free(data_ptr);
    }

    // Cleanup
    fftwf_destroy_plan(plan);
    if (time_grid) cout << "[FFTy " << T << "]";
  }

  // Smallest array
  if (fft_in_z == 1) {
    tictoc T;
    if (time_grid) T.tic();
    int N = Sz;

    // Create Plan
    fftwf_plan plan;
#pragma omp critical
    {
      // Get a plan but never use
      fftwf_complex *ptr =
          (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * N);
      plan = fftwf_plan_dft_1d(N, ptr, ptr, FFTW_BACKWARD, FFTW_MEASURE);
      fftwf_free(ptr);
    }

    // nested
    int DIMS[2];
    DIMS[0] = Ny;
    DIMS[1] = this->Num_Coils;
    long Nlines = DIMS[0] * DIMS[1];

#pragma omp parallel for
    for (long index = 0; index < Nlines; index++) {
      // Get the actual position
      int stride[2];
      nested_workaround(index, DIMS, stride, 2);
      int j = stride[0] + og_sy;
      int coil = stride[1];

      // Storage block
      Array<complex<float>, 2> TEMP(N, Nx, ColumnMajorArray<2>());

      // Gather (i is continous)
      for (int k = 0; k < N; k++) {
        for (int i = og_sx; i < og_ex; i++) {
          TEMP(k, i - og_sx) = k3d_grid(coil)(i, j, k);
        }
      }

      // FFT ( k is continous in this array but its small)
      for (int i = 0; i < Nx; i++) {
        fftwf_complex *data = reinterpret_cast<fftwf_complex *>(&TEMP(0, i));
        fftwf_execute_dft(plan, data, data);
      }

      // Scatter
      for (int k = 0; k < N; k++) {
        for (int i = og_sx; i < og_ex; i++) {
          k3d_grid(coil)(i, j, k) = TEMP(k, i - og_sx);
        }
      }
    }

    // Cleanup
    fftwf_destroy_plan(plan);
    if (time_grid) cout << "[FFTz " << T << "]";
  }
}

void gridFFT_CoilThreaded::zero(void) {
  // nested
  int *N = new int[2];
  N[0] = Sz;
  N[1] = Num_Coils;

  // nested
  long Ns = Sz * Num_Coils;
#pragma omp parallel for
  for (long index = 0; index < Ns; index++) {
    int stride[2];
    nested_workaround(index, N, stride, 2);
    int k = stride[0];
    int coil = stride[1];

    for (int j = 0; j < Sy; j++) {
      for (int i = 0; i < Sx; i++) {
        k3d_grid(coil)(i, j, k) = complex<float>(0.0, 0.0);
      }
    }
  }
}

/**
 * Forward gridding sos
 * @param X image to be accumulated
 * @param data raw k-space data
 * @param kx corrdinate in x
 * @param ky corrdinate in y
 * @param kz corrdinate in z
 * @param kw k-space weight
 */
void gridFFT_CoilThreaded::forward_sos(Array<complex<float>, 3> &X,
                                       const Complex4D &data,
                                       const Array<float, 3> &kx,
                                       const Array<float, 3> &ky,
                                       const Array<float, 3> &kz,
                                       const Array<float, 3> &kw) {
  tictoc T;
  if (time_grid) T.tic();
  zero();
  if (time_grid) cout << "Forward::zero:" << T;

  if (time_grid) T.tic();
  chop_grid_forward(data, kx, ky, kz, kw);  // Grid data to K-Space
  if (time_grid) cout << ",grid:" << T;

  if (time_grid) T.tic();
  do_ifft();
  if (time_grid) cout << ",fft:" << T;

  if (time_grid) T.tic();
  deapp();  // Deapp and multiply by sensitivity map
  if (time_grid) cout << ",deapp:" << T;

  if (time_grid) T.tic();
  accumulate_sos(X);  // Deapp,multiply by sensitivity map, and copy
  if (time_grid) cout << ",accumulate:" << T << endl
                      << flush;
}

/**
 * Forward gridding
 * @param X image to be accumulated
 * @param data raw k-space data
 * @param kx corrdinate in x
 * @param ky corrdinate in y
 * @param kz corrdinate in z
 * @param kw k-space weight
 */
void gridFFT_CoilThreaded::forward(Array<complex<float>, 3> &X,
                                   const Complex4D &data,
                                   const Array<float, 3> &kx,
                                   const Array<float, 3> &ky,
                                   const Array<float, 3> &kz,
                                   const Array<float, 3> &kw) {
  tictoc T;
  if (time_grid) T.tic();
  zero();
  if (time_grid) cout << "Forward::zero:" << T;

  if (time_grid) T.tic();
  chop_grid_forward(data, kx, ky, kz, kw);  // Grid data to K-Space
  if (time_grid) cout << ",grid:" << T;

  if (time_grid) T.tic();
  do_ifft();
  if (time_grid) cout << ",fft:" << T;

  if (time_grid) T.tic();
  accumulate(X);  // Deapp,multiply by sensitivity map, and copy
  if (time_grid) cout << ",accumulate:" << T << endl;
}

/**
 * Forward gridding
 * @param X image to be accumulated
 * @param sensitivity map
 * @param data raw k-space data
 * @param kx corrdinate in x
 * @param ky corrdinate in y
 * @param kz corrdinate in z
 * @param kw k-space weight
 */
void gridFFT_CoilThreaded::forward(Array<complex<float>, 3> &X,
                                   const Complex4D &smap, const Complex4D &data,
                                   const Array<float, 3> &kx,
                                   const Array<float, 3> &ky,
                                   const Array<float, 3> &kz,
                                   const Array<float, 3> &kw) {
  tictoc T;
  if (time_grid) T.tic();
  zero();
  if (time_grid) cout << "Forward::zero:" << T;

  if (time_grid) T.tic();
  chop_grid_forward(data, kx, ky, kz, kw);  // Grid data to K-Space
  if (time_grid) cout << ",grid:" << T;

  if (time_grid) T.tic();
  do_ifft();
  if (time_grid) cout << ",fft:" << T;

  if (time_grid) T.tic();
  accumulate(X, smap);  // Deapp,multiply by sensitivity map, and copy
  if (time_grid) cout << ",accumulate:" << T << endl;
}

/**
 * Backward gridding with sensitivity map
 * @param X image to be accumulated
 * @param sensitivity map
 * @param data raw k-space data
 * @param kx corrdinate in x
 * @param ky corrdinate in y
 * @param kz corrdinate in z
 * @param kw k-space weight
 */
void gridFFT_CoilThreaded::backward(const Array<complex<float>, 3> &X,
                                    const Complex4D &smap, Complex4D &data,
                                    const Array<float, 3> &kx,
                                    const Array<float, 3> &ky,
                                    const Array<float, 3> &kz,
                                    const Array<float, 3> &kw) {
  tictoc T;

  if (time_grid) T.tic();
  set_image(X, smap);  // Copy image to gridding
  if (time_grid) cout << ",copy:" << T << std::flush;

  if (time_grid) T.tic();
  do_fft();
  if (time_grid) cout << ",ifft:" << T << std::flush;

  if (time_grid) T.tic();
  Complex4D temp;
  chop_grid_backward(data, kx, ky, kz, kw, temp, false);  // Inverse gridding
  if (time_grid) cout << ",igrid:" << T << endl
                      << std::flush;
}

/**
 * Backward gridding with sensitivity map
 * @param X image to be accumulated
 * @param sensitivity map
 * @param data raw k-space data
 * @param kx corrdinate in x
 * @param ky corrdinate in y
 * @param kz corrdinate in z
 * @param kw k-space weight
 */
void gridFFT_CoilThreaded::backward_residual(
    const Array<complex<float>, 3> &X, const Complex4D &smap, Complex4D &data,
    const Array<float, 3> &kx, const Array<float, 3> &ky,
    const Array<float, 3> &kz, const Array<float, 3> &kw,
    const Complex4D &diff_data) {
  tictoc T;

  if (time_grid) T.tic();
  set_image(X, smap);  // Copy image to gridding
  if (time_grid) cout << ",copy:" << T << std::flush;

  if (time_grid) T.tic();
  do_fft();
  if (time_grid) cout << ",ifft:" << T << std::flush;

  if (time_grid) T.tic();
  chop_grid_backward(data, kx, ky, kz, kw, diff_data,
                     true);  // Inverse gridding
  if (time_grid) cout << ",igrid:" << T << endl
                      << std::flush;
}

/**
 * Backward gridding without sensitivity map
 * @param X image to be accumulated
 * @param data raw k-space data
 * @param kx corrdinate in x
 * @param ky corrdinate in y
 * @param kz corrdinate in z
 * @param kw k-space weight
 */
void gridFFT_CoilThreaded::backward(const Array<complex<float>, 3> &X,
                                    Complex4D &data, const Array<float, 3> &kx,
                                    const Array<float, 3> &ky,
                                    const Array<float, 3> &kz,
                                    const Array<float, 3> &kw) {
  tictoc T;

  if (time_grid) T.tic();
  set_image(X);  // Copy image to gridding
  if (time_grid) cout << ",copy:" << T;

  if (time_grid) T.tic();
  do_fft();
  if (time_grid) cout << ",ifft:" << T;

  if (time_grid) T.tic();
  Complex4D temp;
  chop_grid_backward(data, kx, ky, kz, kw, temp, false);  // Inverse gridding
  if (time_grid) cout << ",igrid:" << T << endl;
}

//----------------------------------------
//    Crop from Gridding Matrix to Image
//----------------------------------------

void gridFFT_CoilThreaded::deapp(void) {
  for (int coil = 0; coil < this->Num_Coils; coil++) {
#pragma omp parallel for
    for (int k = 0; k < Nz; k++) {
      float wtz = winz(k + og_sz);
      for (int j = 0; j < Ny; j++) {
        float wty = wtz * winy(j + og_sy);
        for (int i = 0; i < Nx; i++) {
          image(coil)(i, j, k) *=
              wty *
              winx(i +
                   og_sx); /* Don't just sum coils when no sense map is given*/
        }
      }
    }
  }
}

void gridFFT_CoilThreaded::accumulate(Array<complex<float>, 3> &X,
                                      const Complex4D &smap) {
  for (int coil = 0; coil < this->Num_Coils; coil++) {
#pragma omp parallel for
    for (int k = 0; k < Nz; k++) {
      float wtz = winz(k + og_sz);

      for (int j = 0; j < Ny; j++) {
        float wty = wtz * winy(j + og_sy);

        for (int i = 0; i < Nx; i++) {
          float wt = wty * winx(i + og_sx);

          X(i, j, k) += wt * image(coil)(i, j, k) * conj(smap(coil)(i, j, k));
        }
      }
    }
  }
}

void gridFFT_CoilThreaded::accumulate(Array<complex<float>, 3> &X) {
  for (int coil = 0; coil < this->Num_Coils; coil++) {
#pragma omp parallel for
    for (int k = 0; k < Nz; k++) {
      float wtz = winz(k + og_sz);

      for (int j = 0; j < Ny; j++) {
        float wty = wtz * winy(j + og_sy);

        for (int i = 0; i < Nx; i++) {
          float wt = wty * winx(i + og_sx);

          X(i, j, k) += wt * image(coil)(i, j, k);
        }
      }
    }
  }
}

void gridFFT_CoilThreaded::accumulate_sos(Array<complex<float>, 3> &X) {
  for (int coil = 0; coil < this->Num_Coils; coil++) {
#pragma omp parallel for
    for (int k = 0; k < Nz; k++) {
      for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
          // Acumulate in a thread safe manner
          X(i, j, k) += image(coil)(i, j, k) * conj(image(coil)(i, j, k));
          ;
        }
      }
    }
  }
}

void gridFFT_CoilThreaded::set_image(const Array<complex<float>, 3> &X,
                                     const Complex4D &smap) {
  for (int coil = 0; coil < this->Num_Coils; coil++) {
#pragma omp parallel for
    for (int k = 0; k < Sz; k++) {
      // Zero the Edges
      if ((k < og_sz) || (k > (og_ez - 1))) {
        for (int j = 0; j < Sy; j++) {
          for (int i = 0; i < Sx; i++) {
            k3d_grid(coil)(i, j, k) = complex<float>(0.0, 0.0);
          }
        }
      } else {
        // Loop y
        float wtz = winz(k);
        for (int j = 0; j < Sy; j++) {
          if ((j < og_sy) || (j > (og_ey - 1))) {
            for (int i = 0; i < Sx; i++) {
              k3d_grid(coil)(i, j, k) = complex<float>(0.0, 0.0);
            }
          } else {
            float wty = wtz * winy(j);
            for (int i = 0; i < Sx; i++) {
              if ((i < og_sx) || (i > (og_ex - 1))) {
                k3d_grid(coil)(i, j, k) = complex<float>(0.0, 0.0);
              } else {
                float wt = wty * winx(i);
                k3d_grid(coil)(i, j, k) =
                    wt * X(i - og_sx, j - og_sy, k - og_sz) *
                    smap(coil)(i - og_sx, j - og_sy, k - og_sz);
              }
            }
          }
        }
      }
    }  // Loop k
  }    // Loop coil
}

void gridFFT_CoilThreaded::set_image(const Array<complex<float>, 3> &X) {
  for (int coil = 0; coil < this->Num_Coils; coil++) {
#pragma omp parallel for
    for (int k = 0; k < Sz; k++) {
      // Zero the Edges
      if ((k < og_sz) || (k > (og_ex - 1))) {
        for (int j = 0; j < Sy; j++) {
          for (int i = 0; i < Sx; i++) {
            image(coil)(i, j, k) = complex<float>(0.0, 0.0);
          }
        }
      } else {
        // Loop y
        float wtz = winz(k);
        for (int j = 0; j < Sy; j++) {
          if ((j < og_sy) || (j > (og_ey - 1))) {
            for (int i = 0; i < Sx; i++) {
              image(coil)(i, j, k) = complex<float>(0.0, 0.0);
            }
          } else {
            float wty = wtz * winy(j);
            for (int i = 0; i < Sx; i++) {
              if ((i < og_sx) || (i > (og_ex - 1))) {
                image(coil)(i, j, k) = complex<float>(0.0, 0.0);
              } else {
                float wt = wty * winx(i);
                image(coil)(i, j, k) = wt * X(i - og_sx, j - og_sy, k - og_sz);
              }
            }
          }
        }
      }
    }  // Loop k
  }    // Loop coil
}

// -------------------------------------------------------
//  This is the main function for Gridding.  Assumes data
//  is already density compensated, etc.
// -------------------------------------------------------

void gridFFT_CoilThreaded::chop_grid_forward(const Complex4D &dataA,
                                             const Array<float, 3> &kxA,
                                             const Array<float, 3> &kyA,
                                             const Array<float, 3> &kzA,
                                             const Array<float, 3> &kwA) {
  float cx = Sx / 2;
  float cy = Sy / 2;
  float cz = Sz / 2;

  // nested
  int *N = new int[4];
  N[0] = dataA(0).length(firstDim);
  N[1] = dataA(0).length(secondDim);
  N[2] = dataA(0).length(thirdDim);
  N[3] = this->Num_Coils;
  long Npts = N[0];
  Npts *= N[1];
  Npts *= N[2];
  Npts *= N[3];

#pragma omp parallel for
  for (long index = 0; index < Npts; index++) {
    // Get the actual position
    int stride[4];
    nested_workaround(index, N, stride, 4);
    int ii = stride[0];
    int jj = stride[1];
    int kk = stride[2];
    int coil = stride[3];

    float kw = kwA(ii, jj, kk);
    if (kw == 0) {
      continue;
    }

    float kx = kxA(ii, jj, kk);
    float ky = kyA(ii, jj, kk);
    float kz = kzA(ii, jj, kk);

    // Calculate the exact kspace sample point in
    // dimension flag->grid* kspace that this point (i,j)
    // is contributing too.

    // Compute Coordinates + Check
    float dkx = kx * grid_x + cx;
    int sx;
    int ex;
    if (grid_in_x) {
      sx = max((int)ceil(dkx - dwinX), 0);
      ex = min((int)floor(dkx + dwinX), Sx - 1);
    } else {
      sx = (int)(dkx);
      ex = sx;
    }
    if (sx >= Sx) continue;
    if (ex < 0) continue;

    float dky = ky * grid_y + cy;
    int sy;
    int ey;
    if (grid_in_y) {
      sy = max((int)ceil(dky - dwinY), 0);
      ey = min((int)floor(dky + dwinY), Sy - 1);
    } else {
      sy = (int)(dky);
      ey = sy;
    }
    if (sy >= Sy) continue;
    if (ey < 0) continue;

    float dkz = kz * grid_z + cz;
    int sz;
    int ez;
    if (grid_in_z) {
      sz = max((int)ceil(dkz - dwinZ), 0);
      ez = min((int)floor(dkz + dwinZ), Sz - 1);
    } else {
      sz = (int)(dkz);
      ez = sz;
    }
    if (sz >= Sz) continue;
    if (ez < 0) continue;

    float kr = 0.0;
    if (fft_in_x) {
      kr += (kx * kx);
    }

    if (fft_in_y) {
      kr += (ky * ky);
    }

    if (fft_in_z) {
      kr += (kz * kz);
    }
    kw *= exp(-kr / (2.0 * k_rad * k_rad));

    /*This is the main loop - most time is spent here*/
    for (int lz = sz; lz <= ez; lz++) {
      float delz = fabs(grid_modZ * (dkz - (float)lz));
      float dz = delz - (float)((int)delz);
      float wtz = grid_filterZ((int)delz) * (1.0 - dz) +
                  grid_filterZ((int)delz + 1) * dz;
      if (!grid_in_z) {
        wtz = 1.0;
      }

      if (fft_in_z) {
        wtz *= ((float)(2 * (lz % 2) - 1));
      }

      for (int ly = sy; ly <= ey; ly++) {
        float dely = fabs(grid_modY * (dky - (float)ly));
        float dy = dely - (float)((int)dely);
        float wty = wtz * (grid_filterY((int)dely) * (1.0 - dy) +
                           grid_filterY((int)dely + 1) * dy);
        if (fft_in_y) {
          wty *= ((float)(2 * (ly % 2) - 1));
        }

        for (int lx = sx; lx <= ex; lx++) {
          float delx = fabs(grid_modX * (dkx - (float)lx));
          float dx = delx - (float)((int)delx);
          float wtx = wty * (grid_filterX((int)delx) * (1.0 - dx) +
                             grid_filterX((int)delx + 1) * dx);

          if (fft_in_x) {
            wtx *= ((float)(2 * (lx % 2) - 1));
          }

          // for(int coil=0; coil < this->Num_Coils; coil++){
          {
            complex<float> temp = dataA(coil)(ii, jj, kk);

            // Density Comp
            temp *= kw;

            // Do not grid zeros
            if (temp == complex<float>(0, 0)) continue;

            complex<float> temp2 = wtx * temp;
            float RD = real(temp2);
            float ID = imag(temp2);
            float *I = reinterpret_cast<float *>(&k3d_grid(coil)(lx, ly, lz));
            float *R = I++;

// Prevent Race conditions in multi-threaded
#pragma omp atomic
            *R += RD;

#pragma omp atomic
            *I += ID;
          }

          /*This Memory Access is the Bottleneck - Also not thread safe!*/
          // k3d_grid.vals[lz][ly][lx]+=temp2;
        } /* end lz loop */
      }   /* end ly */
    }     /* end lx */
  }       /* end data loop */

  return;
}

void gridFFT_CoilThreaded::chop_grid_backward(
    Complex4D &dataA, const Array<float, 3> &kxA, const Array<float, 3> &kyA,
    const Array<float, 3> &kzA, const Array<float, 3> &kwA,
    const Complex4D &diff_dataA, bool sub_data_flag) {
  float cx = Sx / 2;
  float cy = Sy / 2;
  float cz = Sz / 2;

  // nested
  int *N = new int[4];
  N[0] = dataA(0).length(firstDim);
  N[1] = dataA(0).length(secondDim);
  N[2] = dataA(0).length(thirdDim);
  N[3] = this->Num_Coils;
  long Npts = N[0];
  Npts *= N[1];
  Npts *= N[2];
  Npts *= N[3];

#pragma omp parallel for
  for (long index = 0; index < Npts; index++) {
    // Get the actual position
    int stride[4];
    nested_workaround(index, N, stride, 4);
    int ii = stride[0];
    int jj = stride[1];
    int kk = stride[2];
    int coil = stride[3];

    // Storage for value
    complex<float> temp(0.0, 0.0);

    // Do not grid zeros
    float kw = kwA(ii, jj, kk);
    if (kw == 0.0) {
      dataA(coil)(ii, jj, kk) = complex<float>(0.0, 0.0);
      continue;
    }

    float kx = kxA(ii, jj, kk);
    float ky = kyA(ii, jj, kk);
    float kz = kzA(ii, jj, kk);

    // Calculate the exact kspace sample point in
    // dimension flag->grid* kspace that this point (i,j)
    // is contributing too.

    // Compute Coordinates + Check
    float dkx = kx * grid_x + cx;
    int sx;
    int ex;
    if (grid_in_x) {
      sx = max((int)ceil(dkx - dwinX), 0);
      ex = min((int)floor(dkx + dwinX), Sx - 1);
    } else {
      sx = (int)(dkx);
      ex = sx;
    }

    if (sx >= Sx) {
      dataA(coil)(ii, jj, kk) = complex<float>(0.0, 0.0);
      continue;
    }

    if (ex < 0) {
      dataA(coil)(ii, jj, kk) = complex<float>(0.0, 0.0);
      continue;
    }

    float dky = ky * grid_y + cy;
    int sy;
    int ey;
    if (grid_in_y) {
      sy = max((int)ceil(dky - dwinY), 0);
      ey = min((int)floor(dky + dwinY), Sy - 1);
    } else {
      sy = (int)(dky);
      ey = sy;
    }

    if (sy >= Sy) {
      dataA(coil)(ii, jj, kk) = complex<float>(0.0, 0.0);
      continue;
    }

    if (ey < 0) {
      dataA(coil)(ii, jj, kk) = complex<float>(0.0, 0.0);
      continue;
    }

    float dkz = kz * grid_z + cz;
    int sz;
    int ez;
    if (grid_in_z) {
      sz = max((int)ceil(dkz - dwinZ), 0);
      ez = min((int)floor(dkz + dwinZ), Sz - 1);
    } else {
      sz = (int)(dkz);
      ez = sz;
    }
    if (sz >= Sz) {
      dataA(coil)(ii, jj, kk) = complex<float>(0.0, 0.0);
      continue;
    }

    if (ez < 0) {
      dataA(coil)(ii, jj, kk) = complex<float>(0.0, 0.0);
      continue;
    }

    /*This is the main loop - most time is spent here*/
    for (int lz = sz; lz <= ez; lz++) {
      float delz = fabs(grid_modZ * (dkz - (float)lz));
      float dz = delz - (float)((int)delz);
      float wtz = grid_filterZ((int)delz) * (1.0 - dz) +
                  grid_filterZ((int)delz + 1) * dz;
      if (!grid_in_z) {
        wtz = 1.0;
      }

      if (fft_in_z) {
        wtz *= ((float)(2 * (lz % 2) - 1));
      }

      for (int ly = sy; ly <= ey; ly++) {
        float dely = fabs(grid_modY * (dky - (float)ly));
        float dy = dely - (float)((int)dely);
        float wty = wtz * (grid_filterY((int)dely) * (1.0 - dy) +
                           grid_filterY((int)dely + 1) * dy);
        if (fft_in_y) {
          wty *= ((float)(2 * (ly % 2) - 1));
        }

        for (int lx = sx; lx <= ex; lx++) {
          float delx = fabs(grid_modX * (dkx - (float)lx));
          float dx = delx - (float)((int)delx);
          float wtx = wty * (grid_filterX((int)delx) * (1.0 - dx) +
                             grid_filterX((int)delx + 1) * dx);

          if (fft_in_x) {
            wtx *= ((float)(2 * (lx % 2) - 1));
          }

          /*This Memory Access is the Bottleneck*/
          temp += wtx * k3d_grid(coil)(lx, ly, lz);

        } /* end lz loop */
      }   /* end ly */
    }     /* end lx */

    if (sub_data_flag) {
      dataA(coil)(ii, jj, kk) = temp - diff_dataA(coil)(ii, jj, kk);
    } else {
      dataA(coil)(ii, jj, kk) = temp;
    }

  } /* end data loop */
  return;
}

// For Kaiser Bessel Window
float gridFFT_CoilThreaded::bessi0(float x) {
  float ax, ans;
  double y;

  if ((ax = fabs(x)) < 3.75) {
    y = x / 3.75;
    y *= y;
    ans =
        1.0 +
        y * (3.5156229 +
             y * (3.0899424 +
                  y * (1.2067492 +
                       y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
  } else {
    y = 3.75 / ax;
    ans = (exp((double)ax) / sqrt((double)ax)) *
          (0.39894228 +
           y * (0.1328592e-1 +
                y * (0.225319e-2 +
                     y * (-0.157565e-2 +
                          y * (0.916281e-2 +
                               y * (-0.2057706e-1 +
                                    y * (0.2635537e-1 +
                                         y * (-0.1647633e-1 +
                                              y * 0.392377e-2))))))));
  }
  return (ans);
}
