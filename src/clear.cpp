#include "clear.h"
#include "io_templates.hpp"

using arma::cx_mat;
using arma::fmat;
using arma::uvec;
using arma::vec;
using namespace std;
using namespace NDarray;

// ----------------------
// Help Message
// ----------------------
void LOWRANKCOIL::help_message(void) {
  cout << "----------------------------------------------" << endl;
  cout << "  Control for Low Rank Operators " << endl;
  cout << "----------------------------------------------" << endl;

  help_flag("-clear_nx []", "kernel size in x");
  help_flag("-clear_ny []", "kernel size in y");
  help_flag("-clear_nz []", "kernel size in z");
  help_flag("-clear_alpha_coil []", "regularization factor");
  help_flag("-clear_alpha_time []", "regularization factor time");
  help_flag("-clear_alpha_encode []", "regularization factor space");
  help_flag("-clear_amp []", "control jumps of -clear_alpha during iter");
  help_flag("-clear_t2 []", "control jumps of -clear_alpha during iter");
}

LOWRANKCOIL::LOWRANKCOIL() {
  // Kernel Size Defaults
  block_size_x = 4;
  block_size_y = 4;
  block_size_z = 4;
  block_iter = 4;

  clear_alpha_coil = 0.0;
  clear_alpha_time = 0.0;
  clear_alpha_encode = 0.0;
  debug = 0;
  smax = 1.0;
  clear_normalized = false;

  clear_amp = 4;
  clear_t2 = 6;
}

LOWRANKCOIL::LOWRANKCOIL(int numarg, const char **pstring) {
  // Kernel Size Defaults
  block_size_x = 4;
  block_size_y = 4;
  block_size_z = 4;
  block_iter = 4;

  clear_alpha_coil = 0.01;
  clear_alpha_time = 0.0;
  clear_alpha_encode = 0.0;
  debug = 0;
  smax = 1.0;
  clear_normalized = true;
  clear_amp = 4;
  clear_t2 = 6;

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
      int_flag("-clear_nx", block_size_x);
      int_flag("-clear_ny", block_size_y);
      int_flag("-clear_nz", block_size_z);
      int_flag("-clear_block_iter", block_iter);
      float_flag("-clear_alpha_coil", clear_alpha_coil);
      float_flag("-clear_alpha_time", clear_alpha_time);
      float_flag("-clear_alpha_encode", clear_alpha_encode);
      float_flag("-clear_amp", clear_amp);
      float_flag("-clear_t2", clear_t2);
    }
  }
}

//-----------------------------------------------------
// Clear Threshold Estimate
//
//  Just estimate the max singular value over the matrix
//-----------------------------------------------------

void LOWRANKCOIL::update_threshold(Array<Array<complex<float>, 3>, 3> &image, int dim, int iter) {
  // Shorthand
  int Nt = image.extent(firstDim);
  int Ne = image.extent(secondDim);
  int Ncoils = image.extent(thirdDim);
  int Nx = image(0, 0, 0).extent(firstDim);
  int Ny = image(0, 0, 0).extent(secondDim);
  int Nz = image(0, 0, 0).extent(thirdDim);

  // Blocks shouldn't be larger than the dimension
  block_size_x = (block_size_x > Nx) ? (Nx) : (block_size_x);
  block_size_y = (block_size_y > Ny) ? (Ny) : (block_size_y);
  block_size_z = (block_size_z > Nz) ? (Nz) : (block_size_z);

  int N = image.extent(dim);
  int Np = image.numElements() * block_size_x * block_size_y * block_size_z / N;

  cout << "Actual Block Size = " << block_size_x << " x " << block_size_y
       << " x " << block_size_z << endl;
  cout << "Getting Low Rank threshold (N=" << N << ")(Np = " << Np << ")"
       << endl;
  smax = 0.0;

  int block_Nx = (int)(Nx / block_size_x);
  int block_Ny = (int)(Ny / block_size_y);
  int block_Nz = (int)(Nz / block_size_z);

  int total_blocks = block_Nx * block_Ny * block_Nz;
  cout << "Total Block Size" << total_blocks << " ( " << block_Nx << ","
       << block_Ny << "," << block_Nz << ")" << endl;

  vec sblocks(total_blocks);

#pragma omp parallel for
  for (int block = 0; block < total_blocks; block++) {
    // Nested parallelism workaround
    int i = (int)(block % block_Nx);

    int temp = (int)((float)block / (float)block_Nx);
    int j = (int)(temp % block_Ny);
    int k = (int)((float)temp / (float)(block_Ny));

    k *= block_size_z;
    j *= block_size_y;
    i *= block_size_x;

    cx_mat A;
    A.zeros(Np, N);  // Pixels x N

    fmat S;
    S.zeros(Np, N);  // Pixels x N

    //-----------------------------------------------------
    //   Block section
    //-----------------------------------------------------
    switch (dim) {
      case (0): {
        for (int t = 0; t < Nt; t++) {
          int count = 0;
          for (int c = 0; c < Ncoils; c++) {
            for (int e = 0; e < Ne; e++) {
              for (int kk = k; kk < k + block_size_z; kk++) {
                for (int jj = j; jj < j + block_size_y; jj++) {
                  for (int ii = i; ii < i + block_size_x; ii++) {
                    A(count, t) = image(t, e, c)(ii, jj, kk);
                    count++;
                  }
                }
              }
            }
          }
        }
      } break;

      case (1): {
        for (int e = 0; e < Ne; e++) {
          int count = 0;
          for (int c = 0; c < Ncoils; c++) {
            for (int t = 0; t < Nt; t++) {
              for (int kk = k; kk < k + block_size_z; kk++) {
                for (int jj = j; jj < j + block_size_y; jj++) {
                  for (int ii = i; ii < i + block_size_x; ii++) {
                    A(count, t) = image(t, e, c)(ii, jj, kk);
                    count++;
                  }
                }
              }
            }
          }
        }
      } break;

      case (2): {
        for (int c = 0; c < Ncoils; c++) {
          int count = 0;
          for (int e = 0; e < Ne; e++) {
            for (int t = 0; t < Nt; t++) {
              for (int kk = k; kk < k + block_size_z; kk++) {
                for (int jj = j; jj < j + block_size_y; jj++) {
                  for (int ii = i; ii < i + block_size_x; ii++) {
                    A(count, t) = image(t, e, c)(ii, jj, kk);
                    count++;
                  }
                }
              }
            }
          }
        }
      } break;
    }  // Switch Dim

    // SVD
    vec s;
    arma::svd(s, A);

    sblocks(block) = s(0);

  }  // Block (threaded)

  float amp_scale = (1 + clear_amp * exp(-((float)iter) / clear_t2));

  vec nonzero_sblocks = sblocks.elem(find(sblocks));
  vec sorted_sblocks = sort(nonzero_sblocks);
  smax = median(sorted_sblocks) * amp_scale;
  cout << "Amp scale " << amp_scale << endl;
  // cout << "Max singular value is " << max(sblocks) << endl;
  cout << "Median singular value is" << smax << endl;
}

//-----------------------------------------------------
// Clear Threshold Estimate
//
//  Just estimate the max singular value over the matrix
//-----------------------------------------------------
void LOWRANKCOIL::update_threshold(Array<Array<complex<float>, 3>, 2> &image,
                                   int dim, int iter) {
  // Shorthand
  int Nt = image.extent(firstDim);
  int Ne = image.extent(secondDim);
  int Nx = image(0, 0).extent(firstDim);
  int Ny = image(0, 0).extent(secondDim);
  int Nz = image(0, 0).extent(thirdDim);

  // Blocks shouldn't be larger than the dimension
  block_size_x = (block_size_x > Nx) ? (Nx) : (block_size_x);
  block_size_y = (block_size_y > Ny) ? (Ny) : (block_size_y);
  block_size_z = (block_size_z > Nz) ? (Nz) : (block_size_z);

  int N = image.extent(dim);
  int Np =
      image.extent((dim + 1) % 2) * block_size_x * block_size_y * block_size_z;

  cout << "Actual Block Size = " << block_size_x << " x " << block_size_y
       << " x " << block_size_z << endl;
  cout << "Getting Low Rank threshold (N=" << N << ")(Np = " << Np << ")"
       << endl;
  smax = 0.0;

  int block_Nx = (int)(Nx / block_size_x);
  int block_Ny = (int)(Ny / block_size_y);
  int block_Nz = (int)(Nz / block_size_z);

  int total_blocks = block_Nx * block_Ny * block_Nz;
  cout << "Total Block Size" << total_blocks << " ( " << block_Nx << ","
       << block_Ny << "," << block_Nz << ")" << endl;

  vec sblocks(total_blocks);

  arma::Col<int> act_i(total_blocks);
  arma::Col<int> act_j(total_blocks);
  arma::Col<int> act_k(total_blocks);
  {
    int count = 0;
    for (int i = 0; i < block_Nx; i++) {
      for (int j = 0; j < block_Ny; j++) {
        for (int k = 0; k < block_Nz; k++) {
          act_i(count) = i;
          act_j(count) = j;
          act_k(count) = k;
          count++;
        }
      }
    }
  }

#pragma omp parallel for
  for (int block = 0; block < total_blocks; block++) {
    // Nested parallelism workaround
    int i = act_i(block) * block_size_x;
    int j = act_j(block) * block_size_y;
    int k = act_k(block) * block_size_z;

    // cout << "B = " << block << "(" << i << "," << j << "," << k << ")" <<
    // endl;

    // Set to one to prevent  downstream error
    cx_mat A;
    A.ones(Np, N);  // Pixels x N

    //-----------------------------------------------------
    //   Block section
    //-----------------------------------------------------
    switch (dim) {
      case (0): {
        for (int t = 0; t < Nt; t++) {
          int count = 0;
          for (int e = 0; e < Ne; e++) {
            for (int kk = k; kk < k + block_size_z; kk++) {
              for (int jj = j; jj < j + block_size_y; jj++) {
                for (int ii = i; ii < i + block_size_x; ii++) {
                  A(count, t) = image(t, e)(ii, jj, kk);
                  count++;
                }
              }
            }
          }
        }
      } break;

      case (1): {
        for (int e = 0; e < Ne; e++) {
          int count = 0;
          for (int t = 0; t < Nt; t++) {
            for (int kk = k; kk < k + block_size_z; kk++) {
              for (int jj = j; jj < j + block_size_y; jj++) {
                for (int ii = i; ii < i + block_size_x; ii++) {
                  A(count, e) = image(t, e)(ii, jj, kk);
                  count++;
                }
              }
            }
          }
        }
      } break;
    }  // Switch Dim

    // SVD
    vec s;
    arma::svd(s, A);

    sblocks(block) = s(0);

  }  // Block (threaded)

  float amp_scale = (1 + clear_amp * exp(-((float)iter) / clear_t2));

  vec nonzero_sblocks = sblocks.elem(find(sblocks));
  vec sorted_sblocks = sort(nonzero_sblocks);
  smax = median(sorted_sblocks) * amp_scale;
  cout << "Amp scale " << amp_scale << endl;
  // cout << "Max singular value is " << max(sblocks) << endl;
  cout << "Median singular value is" << smax << endl;
}

//-----------------------------------------------------
//  Clear Threshold Call
//-----------------------------------------------------
void LOWRANKCOIL::thresh(Array<Array<complex<float>, 3>, 2> &image, int dim) {
  // Shorthand
  int Nt = image.extent(firstDim);
  int Ne = image.extent(secondDim);
  int Nx = image(0, 0).extent(firstDim);
  int Ny = image(0, 0).extent(secondDim);
  int Nz = image(0, 0).extent(thirdDim);

  int block_Nx = (int)(Nx / block_size_x);
  int block_Ny = (int)(Ny / block_size_y);
  int block_Nz = (int)(Nz / block_size_z);
  int total_blocks = block_Nx * block_Ny * block_Nz;

  int N = image.extent(dim);
  int Np = image.extent((dim + 1) % 2) * block_size_x * block_size_y * block_size_z;

  cout << " Low Rank threshold (N=" << N << ")(Np = " << Np << ")" << endl;

  double clear_alpha = 0.0;
  switch (dim) {
    case (0): {
      clear_alpha = clear_alpha_time;
    } break;
    case (1): {
      clear_alpha = clear_alpha_encode;
    } break;
  }

  if (clear_alpha == 0.0) {
    return;
  }
  clear_alpha *= smax;
  cout << "Clear alpha " << clear_alpha << endl;

  int count = 0;
  arma::Col<int> act_i(total_blocks);
  arma::Col<int> act_j(total_blocks);
  arma::Col<int> act_k(total_blocks);
  for (int i = 0; i < block_Nx; i++) {
    for (int j = 0; j < block_Ny; j++) {
      for (int k = 0; k < block_Nz; k++) {
        act_i(count) = i;
        act_j(count) = j;
        act_k(count) = k;
        count++;
      }
    }
  }

  for (int iteration = 0; iteration < block_iter; iteration++) {
    cout << "Block iter =" << iteration << endl;

    int shiftx = rand() % block_size_x - block_size_x / 2;
    int shifty = rand() % block_size_y - block_size_y / 2;
    int shiftz = rand() % block_size_z - block_size_z / 2;
    if (block_iter == 1) {
      shiftx = 0;
      shifty = 0;
      shiftz = 0;
    }

#pragma omp parallel for
    for (int block = 0; block < total_blocks; block++) {
      // cout << "Block = " << block << endl << flush;

      // Nested parallelism workaround
      int i = act_i(block) * block_size_x;
      int j = act_j(block) * block_size_y;
      int k = act_k(block) * block_size_z;

      // cout << "Index = " << i << "," << j << "," << k << endl << flush;
      // cout << "Allocating for " << Np << " x " << N << endl << flush;

      arma::cx_mat A;
      A.zeros(Np, N);  // Pixels x Coils

      //-----------------------------------------------------
      //   Block Gather
      //-----------------------------------------------------

      // cout << "Gather (Block =" << i << "," << j << "," << k << ")" << endl<<
      // flush<< flush; cout << "	    (Size  =" << block_size_x << "," <<
      // block_size_y << "," << block_size_z << ")" << endl<< flush;

      switch (dim) {
        case (0): {
          for (int t = 0; t < Nt; t++) {
            int count = 0;
            for (int e = 0; e < Ne; e++) {
              for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz);
                   kk++) {
                for (int jj = (j + shifty); jj < (j + block_size_y + shifty);
                     jj++) {
                  for (int ii = (i + shiftx); ii < (i + block_size_x + shiftx);
                       ii++) {
                    int px = (ii + Nx) % Nx;
                    int py = (jj + Ny) % Ny;
                    int pz = (kk + Nz) % Nz;
                    A(count, t) = image(t, e)(px, py, pz);
                    count++;
                  }
                }
              }
            }
          }
        } break;

        case (1): {
          for (int e = 0; e < Ne; e++) {
            int count = 0;
            for (int t = 0; t < Nt; t++) {
              for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz);
                   kk++) {
                for (int jj = (j + shifty); jj < (j + block_size_y + shifty);
                     jj++) {
                  for (int ii = (i + shiftx); ii < (i + block_size_x + shiftx);
                       ii++) {
                    int px = (ii + Nx) % Nx;
                    int py = (jj + Ny) % Ny;
                    int pz = (kk + Nz) % Nz;
                    A(count, e) = image(t, e)(px, py, pz);
                    count++;
                  }
                }
              }
            }
          }
        } break;
      }  // Switch Dim

      // --------------------------------------------------------
      // SVD Threshold
      // --------------------------------------------------------

      /* Do SVD*/
      arma::cx_mat U;
      arma::cx_mat V;
      arma::vec s;
      arma::svd_econ(U, s, V, A, "both", "std");

      for (vec::iterator miter = s.begin(); miter != s.end(); miter++) {
        (*miter) = max((*miter) - clear_alpha / block_iter, 0.0);
      }

      // Reconstruct
      arma::cx_mat X = U * diagmat(s) * V.t();

      // --------------------------------------------------------
      //  Block Scatter
      // --------------------------------------------------------

      switch (dim) {
        case (0): {
          for (int t = 0; t < Nt; t++) {
            int count = 0;
            for (int e = 0; e < Ne; e++) {
              for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz);
                   kk++) {
                for (int jj = (j + shifty); jj < (j + block_size_y + shifty);
                     jj++) {
                  for (int ii = (i + shiftx); ii < (i + block_size_x + shiftx);
                       ii++) {
                    int px = (ii + Nx) % Nx;
                    int py = (jj + Ny) % Ny;
                    int pz = (kk + Nz) % Nz;
                    image(t, e)(px, py, pz) = X(count, t);
                    count++;
                  }
                }
              }
            }
          }
        } break;

        case (1): {
          for (int e = 0; e < Ne; e++) {
            int count = 0;
            for (int t = 0; t < Nt; t++) {
              for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz);
                   kk++) {
                for (int jj = (j + shifty); jj < (j + block_size_y + shifty);
                     jj++) {
                  for (int ii = (i + shiftx); ii < (i + block_size_x + shiftx);
                       ii++) {
                    int px = (ii + Nx) % Nx;
                    int py = (jj + Ny) % Ny;
                    int pz = (kk + Nz) % Nz;
                    image(t, e)(px, py, pz) = X(count, e);
                    count++;
                  }
                }
              }
            }
          }
        } break;
      }  // Switch Dim

    }  // Block Loop

  }  // iteration
}

//-----------------------------------------------------
//  Clear Threshold Call
//-----------------------------------------------------
void LOWRANKCOIL::thresh(Array<Array<complex<float>, 3>, 3> &image, int dim) {
  // Shorthand
  int Nt = image.extent(firstDim);
  int Ne = image.extent(secondDim);
  int Ncoils = image.extent(thirdDim);
  int Nx = image(0, 0).extent(firstDim);
  int Ny = image(0, 0).extent(secondDim);
  int Nz = image(0, 0).extent(thirdDim);

  int block_Nx = (int)(Nx / block_size_x);
  int block_Ny = (int)(Ny / block_size_y);
  int block_Nz = (int)(Nz / block_size_z);
  int total_blocks = block_Nx * block_Ny * block_Nz;

  int N = image.extent(dim);
  int Np = image.numElements() * block_size_x * block_size_y * block_size_z / N;

  cout << " Low Rank threshold (N=" << N << ")(Np = " << Np << ")" << endl;

  double clear_alpha = 0.0;
  switch (dim) {
    case (0): {
      clear_alpha = clear_alpha_time;
    } break;
    case (1): {
      clear_alpha = clear_alpha_encode;
    } break;
    case (2): {
      clear_alpha = clear_alpha_coil;
    } break;
  }

  if (clear_alpha == 0.0) {
    return;
  }
  clear_alpha *= smax;

  /*Pre allocate Matrices*/
  int Nthreads = omp_get_max_threads();
  cx_mat *AA = new cx_mat[Nthreads];
  cx_mat *UU = new cx_mat[Nthreads];
  cx_mat *SS = new cx_mat[Nthreads];
  cx_mat *VV = new cx_mat[Nthreads];

  for (int t = 0; t < Nthreads; t++) {
    AA[t].zeros(Np, N);   // Pixels x Coils
    UU[t].zeros(Np, Np);  // Pixels x Pixels
    SS[t].zeros(Np, N);   // Pixels x coils
    VV[t].zeros(N, N);    // Coils x Coils
  }

  for (int iteration = 0; iteration < block_iter; iteration++) {
    int shiftx = rand() % block_size_x - block_size_x / 2;
    int shifty = rand() % block_size_y - block_size_y / 2;
    int shiftz = rand() % block_size_z - block_size_z / 2;

#pragma omp parallel for
    for (int block = 0; block < total_blocks; block++) {
      int thread = omp_get_thread_num();

      // Nested parallelism workaround
      int i = (int)(block % block_Nx);
      int temp = (int)((float)block / (float)block_Nx);
      int j = (int)(temp % block_Ny);
      int k = (int)((float)temp / (float)(block_Ny));

      k *= block_size_z;
      j *= block_size_y;
      i *= block_size_x;

      //-----------------------------------------------------
      //   Block Gather
      //-----------------------------------------------------
      switch (dim) {
        case (0): {
          for (int t = 0; t < Nt; t++) {
            int count = 0;
            for (int c = 0; c < Ncoils; c++) {
              for (int e = 0; e < Ne; e++) {
                for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz);
                     kk++) {
                  for (int jj = (j + shifty); jj < (j + block_size_y + shifty);
                       jj++) {
                    for (int ii = (i + shiftx);
                         ii < (i + block_size_x + shiftx); ii++) {
                      int px = (ii + Nx) % Nx;
                      int py = (jj + Ny) % Ny;
                      int pz = (kk + Nz) % Nz;
                      AA[thread](count, t) = image(t, e, c)(px, py, pz);
                      count++;
                    }
                  }
                }
              }
            }
          }
        } break;

        case (1): {
          for (int e = 0; e < Ne; e++) {
            int count = 0;
            for (int c = 0; c < Ncoils; c++) {
              for (int t = 0; t < Nt; t++) {
                for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz);
                     kk++) {
                  for (int jj = (j + shifty); jj < (j + block_size_y + shifty);
                       jj++) {
                    for (int ii = (i + shiftx);
                         ii < (i + block_size_x + shiftx); ii++) {
                      int px = (ii + Nx) % Nx;
                      int py = (jj + Ny) % Ny;
                      int pz = (kk + Nz) % Nz;
                      AA[thread](count, e) = image(t, e, c)(px, py, pz);
                      count++;
                    }
                  }
                }
              }
            }
          }
        } break;

        case (2): {
          for (int c = 0; c < Ncoils; c++) {
            int count = 0;
            for (int e = 0; e < Ne; e++) {
              for (int t = 0; t < Nt; t++) {
                for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz);
                     kk++) {
                  for (int jj = (j + shifty); jj < (j + block_size_y + shifty);
                       jj++) {
                    for (int ii = (i + shiftx);
                         ii < (i + block_size_x + shiftx); ii++) {
                      int px = (ii + Nx) % Nx;
                      int py = (jj + Ny) % Ny;
                      int pz = (kk + Nz) % Nz;
                      AA[thread](count, c) = image(t, e, c)(px, py, pz);
                      count++;
                    }
                  }
                }
              }
            }
          }
        } break;
      }  // Switch Dim

      // --------------------------------------------------------
      // SVD Threshold
      // --------------------------------------------------------
      vec s;
      arma::svd(UU[thread], s, VV[thread], AA[thread], "std");

      for (vec::iterator miter = s.begin(); miter != s.end(); miter++) {
        (*miter) = max((*miter) - clear_alpha / block_iter, 0.0);
      }

      for (int pos = 0; pos < min(N, Np); pos++) {
        SS[thread](pos, pos) = s(pos);
      }

      // Reconstruct
      AA[thread] = UU[thread] * SS[thread] * trans(VV[thread]);

      // --------------------------------------------------------
      //  Block Scatter
      // --------------------------------------------------------
      switch (dim) {
        case (0): {
          for (int t = 0; t < Nt; t++) {
            int count = 0;
            for (int c = 0; c < Ncoils; c++) {
              for (int e = 0; e < Ne; e++) {
                for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz);
                     kk++) {
                  for (int jj = (j + shifty); jj < (j + block_size_y + shifty);
                       jj++) {
                    for (int ii = (i + shiftx);
                         ii < (i + block_size_x + shiftx); ii++) {
                      int px = (ii + Nx) % Nx;
                      int py = (jj + Ny) % Ny;
                      int pz = (kk + Nz) % Nz;
                      image(t, e, c)(px, py, pz) = AA[thread](count, t);
                      count++;
                    }
                  }
                }
              }
            }
          }
        } break;

        case (1): {
          for (int e = 0; e < Ne; e++) {
            int count = 0;
            for (int c = 0; c < Ncoils; c++) {
              for (int t = 0; t < Nt; t++) {
                for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz);
                     kk++) {
                  for (int jj = (j + shifty); jj < (j + block_size_y + shifty);
                       jj++) {
                    for (int ii = (i + shiftx);
                         ii < (i + block_size_x + shiftx); ii++) {
                      int px = (ii + Nx) % Nx;
                      int py = (jj + Ny) % Ny;
                      int pz = (kk + Nz) % Nz;
                      image(t, e, c)(px, py, pz) = AA[thread](count, e);
                      count++;
                    }
                  }
                }
              }
            }
          }
        } break;

        case (2): {
          for (int c = 0; c < Ncoils; c++) {
            int count = 0;
            for (int e = 0; e < Ne; e++) {
              for (int t = 0; t < Nt; t++) {
                for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz);
                     kk++) {
                  for (int jj = (j + shifty); jj < (j + block_size_y + shifty);
                       jj++) {
                    for (int ii = (i + shiftx);
                         ii < (i + block_size_x + shiftx); ii++) {
                      int px = (ii + Nx) % Nx;
                      int py = (jj + Ny) % Ny;
                      int pz = (kk + Nz) % Nz;
                      image(t, e, c)(px, py, pz) = AA[thread](count, c);
                      count++;
                    }
                  }
                }
              }
            }
          }
        } break;
      }  // Switch Dim

    }  // Block Loop

  }  // iteration

  delete[] AA;
  delete[] UU;
  delete[] SS;
  delete[] VV;
}

//-----------------------------------------------------
//  Clear Combine Call
//-----------------------------------------------------
void LOWRANKCOIL::combine(Array<Array<complex<float>, 3>, 3> &image,
                          Array<Array<complex<float>, 3>, 2> &im) {
  // Shorthand
  int Nt = image.extent(firstDim);
  int Ne = image.extent(secondDim);
  int Ncoils = image.extent(thirdDim);
  int Nx = image(0, 0).extent(firstDim);
  int Ny = image(0, 0).extent(secondDim);
  int Nz = image(0, 0).extent(thirdDim);

  int block_Nx = (int)(Nx / block_size_x);
  int block_Ny = (int)(Ny / block_size_y);
  int block_Nz = (int)(Nz / block_size_z);
  int total_blocks = block_Nx * block_Ny * block_Nz;

  int N = image.extent(thirdDim);
  int Np = image.numElements() * block_size_x * block_size_y * block_size_z / N;

  int shiftx = 0;
  int shifty = 0;
  int shiftz = 0;

#pragma omp parallel for
  for (int block = 0; block < total_blocks; block++) {
    // Nested parallelism workaround
    int i = (int)(block % block_Nx);
    int temp = (int)((float)block / (float)block_Nx);
    int j = (int)(temp % block_Ny);
    int k = (int)((float)temp / (float)(block_Ny));

    k *= block_size_z;
    j *= block_size_y;
    i *= block_size_x;

    cx_mat A;
    A.zeros(Np, N);  // Pixels x Coils

    cx_mat U;
    U.zeros(Np, Np);  // Pixels x Pixels

    fmat S;
    S.zeros(Np, N);  // Pixels x Coils

    cx_mat V;
    V.zeros(N, N);  // Coils x Coils

    //-----------------------------------------------------
    //   Block Gather
    //-----------------------------------------------------
    for (int c = 0; c < Ncoils; c++) {
      int count = 0;
      for (int e = 0; e < Ne; e++) {
        for (int t = 0; t < Nt; t++) {
          for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz); kk++) {
            for (int jj = (j + shifty); jj < (j + block_size_y + shifty); jj++) {
              for (int ii = (i + shiftx); ii < (i + block_size_x + shiftx); ii++) {
                int px = (ii + Nx) % Nx;
                int py = (jj + Ny) % Ny;
                int pz = (kk + Nz) % Nz;
                A(count, c) = image(t, e, c)(px, py, pz);
                count++;
              }
            }
          }
        }
      }
    }

    // --------------------------------------------------------
    // SVD Threshold
    // --------------------------------------------------------
    vec s;
    arma::svd(U, s, V, A, "std");

    // Grab Largest Value
    cx_mat A2 = A * V.cols(0, 0);

    // --------------------------------------------------------
    //  Block Scatter
    // --------------------------------------------------------

    int count = 0;
    for (int e = 0; e < Ne; e++) {
      for (int t = 0; t < Nt; t++) {
        for (int kk = (k + shiftz); kk < (k + block_size_z + shiftz); kk++) {
          for (int jj = (j + shifty); jj < (j + block_size_y + shifty); jj++) {
            for (int ii = (i + shiftx); ii < (i + block_size_x + shiftx); ii++) {
              int px = (ii + Nx) % Nx;
              int py = (jj + Ny) % Ny;
              int pz = (kk + Nz) % Nz;
              im(t, e)(px, py, pz) = A2(count, 0);
              count++;
            }
          }
        }
      }
    }

  } /*block*/
}
