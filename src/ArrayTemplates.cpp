#include "ArrayTemplates.hpp"

void NDarray::fftshift(Array<complex<float>, 3> &temp) {
  for (int k = 0; k < temp.extent(thirdDim); k++) {
    for (int j = 0; j < temp.extent(secondDim); j++) {
      for (int i = 0; i < temp.extent(firstDim); i++) {
        float mod = ((float)(2 * ((i + j + k) % 2) - 1));
        temp(i, j, k) *= mod;
      }
    }
  }
}
void NDarray::fftshift(Array<complex<float>, 4> &temp) {
  for (int t = 0; t < temp.extent(fourthDim); t++) {
    for (int k = 0; k < temp.extent(thirdDim); k++) {
      for (int j = 0; j < temp.extent(secondDim); j++) {
        for (int i = 0; i < temp.extent(firstDim); i++) {
          float mod = ((float)(2 * ((i + j + k + t) % 2) - 1));
          temp(i, j, k, t) *= mod;
        }
      }
    }
  }
}

void NDarray::nested_workaround(long index, int *N, int *idx, int total) {
  long tempi = index;
  for (int pos = 0; pos < total; pos++) {
    idx[pos] = tempi % N[pos];
    tempi = (tempi - idx[pos]) / N[pos];
  }
}

void NDarray::ifft(Array<complex<float>, 3> &temp) {
  // Shift to Center of K-space
  fftshift(temp);

  // Do FFT
  fftwf_complex *ptr = NULL;
  complex<float> *data = NULL;
  if (temp.isStorageContiguous()) {
    ptr = reinterpret_cast<fftwf_complex *>(temp.data());
  } else {
    cout << "Warning:: Doing FFT of Non-Contiguous Data. Will be slow!" << endl;
    data = new complex<float>[temp.numElements()];
    int pos = 0;
    for (int k = 0; k < temp.extent(thirdDim); k++) {
      for (int j = 0; j < temp.extent(secondDim); j++) {
        for (int i = 0; i < temp.extent(firstDim); i++) {
          data[pos] = temp(i, j, k);
          pos++;
        }
      }
    }
    ptr = reinterpret_cast<fftwf_complex *>(data);
  }

  fftwf_plan fft_plan = fftwf_plan_dft_3d(
      temp.length(thirdDim), temp.length(secondDim), temp.length(firstDim), ptr,
      ptr, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwf_execute(fft_plan);

  if (!temp.isStorageContiguous()) {
    int pos = 0;
    for (int k = 0; k < temp.extent(thirdDim); k++) {
      for (int j = 0; j < temp.extent(secondDim); j++) {
        for (int i = 0; i < temp.extent(firstDim); i++) {
          temp(i, j, k) = data[pos];
          pos++;
        }
      }
    }
    free(ptr);
  }

  // Fix Phase from Shift
  fftshift(temp);

  // Scale
  float scale =
      1.0 / sqrt((float)temp.extent(thirdDim) * (float)temp.extent(secondDim) *
                 (float)temp.extent(firstDim));
  temp *= scale;

  // Cleanup
  fftwf_destroy_plan(fft_plan);
}

void NDarray::ifft(Array<complex<float>, 4> &temp) {
  // Shift to Center of K-space
  fftshift(temp);

  // Do FFT
  fftwf_complex *ptr = NULL;
  complex<float> *data = NULL;
  if (temp.isStorageContiguous()) {
    ptr = reinterpret_cast<fftwf_complex *>(temp.data());
  } else {
    cout << "Warning:: Doing FFT of Non-Contiguous Data. Will be slow!" << endl;
    data = new complex<float>[temp.numElements()];
    int pos = 0;
    for (int t = 0; t < temp.extent(fourthDim); t++) {
      for (int k = 0; k < temp.extent(thirdDim); k++) {
        for (int j = 0; j < temp.extent(secondDim); j++) {
          for (int i = 0; i < temp.extent(firstDim); i++) {
            data[pos++] = temp(i, j, k, t);
          }
        }
      }
    }
    ptr = reinterpret_cast<fftwf_complex *>(data);
  }

  int n[4];
  n[0] = temp.length(fourthDim);
  n[1] = temp.length(thirdDim);
  n[2] = temp.length(secondDim);
  n[3] = temp.length(firstDim);
  fftwf_plan fft_plan =
      fftwf_plan_dft(4, n, ptr, ptr, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwf_execute(fft_plan);

  if (!temp.isStorageContiguous()) {
    int pos = 0;
    for (int t = 0; t < temp.extent(fourthDim); t++) {
      for (int k = 0; k < temp.extent(thirdDim); k++) {
        for (int j = 0; j < temp.extent(secondDim); j++) {
          for (int i = 0; i < temp.extent(firstDim); i++) {
            temp(i, j, k, t) = data[pos++];
          }
        }
      }
    }
    free(ptr);
  }

  // Fix Phase from Shift
  fftshift(temp);

  // Scale
  float scale =
      1.0 / sqrt((float)temp.extent(fourthDim) * (float)temp.extent(thirdDim) *
                 (float)temp.extent(secondDim) * (float)temp.extent(firstDim));
  temp *= scale;

  // Cleanup
  fftwf_destroy_plan(fft_plan);
}

void NDarray::fft(Array<complex<float>, 3> &temp) {
  // Shift to Center of K-space
  fftshift(temp);

  // Do FFT
  fftwf_complex *ptr = NULL;
  complex<float> *data = NULL;
  if (temp.isStorageContiguous()) {
    ptr = reinterpret_cast<fftwf_complex *>(temp.data());
  } else {
    cout << "Warning:: Doing FFT of Non-Contiguous Data. Will be slow!" << endl;
    data = new complex<float>[temp.numElements()];
    int pos = 0;
    for (int k = 0; k < temp.extent(thirdDim); k++) {
      for (int j = 0; j < temp.extent(secondDim); j++) {
        for (int i = 0; i < temp.extent(firstDim); i++) {
          data[pos] = temp(i, j, k);
          pos++;
        }
      }
    }
    ptr = reinterpret_cast<fftwf_complex *>(data);
  }

  fftwf_plan fft_plan = fftwf_plan_dft_3d(
      temp.length(thirdDim), temp.length(secondDim), temp.length(firstDim), ptr,
      ptr, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(fft_plan);

  if (!temp.isStorageContiguous()) {
    int pos = 0;
    for (int k = 0; k < temp.extent(thirdDim); k++) {
      for (int j = 0; j < temp.extent(secondDim); j++) {
        for (int i = 0; i < temp.extent(firstDim); i++) {
          temp(i, j, k) = data[pos];
          pos++;
        }
      }
    }
    free(ptr);
  }

  // Fix Phase from Shift
  fftshift(temp);

  // Scale
  float scale =
      1.0 / sqrt((float)temp.extent(thirdDim) * (float)temp.extent(secondDim) *
                 (float)temp.extent(firstDim));
  temp *= scale;

  // Cleanup
  fftwf_destroy_plan(fft_plan);
}

void NDarray::fft(Array<complex<float>, 4> &temp) {
  // Shift to Center of K-space
  fftshift(temp);

  // Do FFT
  fftwf_complex *ptr = NULL;
  complex<float> *data = NULL;
  if (temp.isStorageContiguous()) {
    ptr = reinterpret_cast<fftwf_complex *>(temp.data());
  } else {
    cout << "Warning:: Doing FFT of Non-Contiguous Data. Will be slow!" << endl;
    data = new complex<float>[temp.numElements()];
    int pos = 0;
    for (int t = 0; t < temp.extent(fourthDim); t++) {
      for (int k = 0; k < temp.extent(thirdDim); k++) {
        for (int j = 0; j < temp.extent(secondDim); j++) {
          for (int i = 0; i < temp.extent(firstDim); i++) {
            data[pos++] = temp(i, j, k, t);
          }
        }
      }
    }
    ptr = reinterpret_cast<fftwf_complex *>(data);
  }

  int n[4];
  n[0] = temp.length(fourthDim);
  n[1] = temp.length(thirdDim);
  n[2] = temp.length(secondDim);
  n[3] = temp.length(firstDim);
  fftwf_plan fft_plan =
      fftwf_plan_dft(4, n, ptr, ptr, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(fft_plan);

  if (!temp.isStorageContiguous()) {
    int pos = 0;
    for (int t = 0; t < temp.extent(fourthDim); t++) {
      for (int k = 0; k < temp.extent(thirdDim); k++) {
        for (int j = 0; j < temp.extent(secondDim); j++) {
          for (int i = 0; i < temp.extent(firstDim); i++) {
            temp(i, j, k, t) = data[pos++];
          }
        }
      }
    }
    free(ptr);
  }

  // Fix Phase from Shift
  fftshift(temp);

  // Scale
  float scale =
      1.0 / sqrt((float)temp.extent(fourthDim) * (float)temp.extent(thirdDim) *
                 (float)temp.extent(secondDim) * (float)temp.extent(firstDim));
  temp *= scale;

  // Cleanup
  fftwf_destroy_plan(fft_plan);
}

// FFT in only one dimension
void NDarray::fft(Array<complex<float>, 3> &temp, int dim) {
  fft3(temp, dim, FFTW_FORWARD, 1);
}

void NDarray::ifft(Array<complex<float>, 3> &temp, int dim) {
  fft3(temp, dim, FFTW_BACKWARD, 1);
}

double NDarray::Dmax(const Array<Array<double, 2>, 1> &A) {
  double val = 0;
  bool init = false;
  for (Array<Array<double, 2>, 1>::const_iterator miter = A.begin();
       miter != A.end(); miter++) {
    double temp = max(*miter);
    if (init == false) {
      val = temp;
    } else if (temp > val) {
      val = temp;
    }
  }
  return (val);
}

double NDarray::Dmin(const Array<Array<double, 2>, 1> &A) {
  double val = 0;
  bool init = false;
  for (Array<Array<double, 2>, 1>::const_iterator miter = A.begin();
       miter != A.end(); miter++) {
    double temp = min(*miter);
    if (init == false) {
      val = temp;
    } else if (temp < val) {
      val = temp;
    }
  }
  return (val);
}

void NDarray::fft3(Array<complex<float>, 3> &temp, int dim, int direction,
                   bool chop) {
  // Get Size
  int N;
  switch (dim) {
    case (0): {
      N = temp.length(firstDim);
    } break;
    case (1): {
      N = temp.length(secondDim);
    } break;
    case (2): {
      N = temp.length(thirdDim);
    } break;
    default: {
      cout << "Error trying to FFT a dimension that doesn't exist" << endl;
      exit(1);
    }
  }

  fftwf_plan plan;

#pragma omp critical
  {
    //  FFT planning is not thread safe
    // fftwf_init_threads();
    // fftwf_plan_with_nthreads(1);

    // Get a plan but never use
    complex<float> *data_temp = new complex<float>[N];
    fftwf_complex *ptr = reinterpret_cast<fftwf_complex *>(data_temp);
    plan = fftwf_plan_dft_1d(N, ptr, ptr, direction, FFTW_MEASURE);
    delete[] data_temp;
  }

  float scale = 1. / sqrt((float)N);

  switch (dim) {
    case (0): {
#pragma omp parallel for
      for (int k = 0; k < temp.extent(thirdDim); k++) {
        complex<float> *data = new complex<float>[N];
        fftwf_complex *data_ptr = reinterpret_cast<fftwf_complex *>(data);

        for (int j = 0; j < temp.extent(secondDim); j++) {
          // Copy
          if (chop == 1) {
            for (int i = 0; i < temp.extent(firstDim); i++) {
              data[i] = temp(i, j, k) * ((float)(2 * (i % 2) - 1));
            }
          } else {
            for (int i = 0; i < temp.extent(firstDim); i++) {
              data[i] = temp(i, j, k);
            }
          }

          // FFT
          fftwf_execute_dft(plan, data_ptr, data_ptr);

          // Copy Back
          if (chop) {
            for (int i = 0; i < temp.extent(firstDim); i++) {
              temp(i, j, k) = data[i] * ((float)(2 * (i % 2) - 1)) * scale;
            }
          } else {
            for (int i = 0; i < temp.extent(firstDim); i++) {
              temp(i, j, k) = data[i] * scale;
            }
          }
        }
        delete[] data;
      }
    } break;

    case (1): {
#pragma omp parallel for
      for (int k = 0; k < temp.extent(thirdDim); k++) {
        complex<float> *data = new complex<float>[N];
        fftwf_complex *data_ptr = reinterpret_cast<fftwf_complex *>(data);

        for (int i = 0; i < temp.extent(firstDim); i++) {
          // Copy
          if (chop == 1) {
            for (int j = 0; j < temp.extent(secondDim); j++) {
              data[j] = temp(i, j, k) * ((float)(2 * (j % 2) - 1));
            }
          } else {
            for (int j = 0; j < temp.extent(secondDim); j++) {
              data[j] = temp(i, j, k);
            }
          }
          // FFT
          fftwf_execute_dft(plan, data_ptr, data_ptr);

          // Copy Back
          if (chop) {
            for (int j = 0; j < temp.extent(secondDim); j++) {
              temp(i, j, k) = data[j] * ((float)(2 * (j % 2) - 1)) * scale;
            }
          } else {
            for (int j = 0; j < temp.extent(secondDim); j++) {
              temp(i, j, k) = data[j] * scale;
            }
          }
        }
        delete[] data;
      }

    } break;

    case (2): {
#pragma omp parallel for
      for (int j = 0; j < temp.extent(secondDim); j++) {
        complex<float> *data = new complex<float>[N];
        fftwf_complex *data_ptr = reinterpret_cast<fftwf_complex *>(data);

        for (int i = 0; i < temp.extent(firstDim); i++) {
          // Copy
          if (chop) {
            for (int k = 0; k < temp.extent(thirdDim); k++) {
              data[k] = temp(i, j, k) * ((float)(2 * (k % 2) - 1));
            }
          } else {
            for (int k = 0; k < temp.extent(thirdDim); k++) {
              data[k] = temp(i, j, k);
            }
          }

          // FFT
          fftwf_execute_dft(plan, data_ptr, data_ptr);

          // Copy Back
          if (chop) {
            for (int k = 0; k < temp.extent(thirdDim); k++) {
              temp(i, j, k) = data[k] * ((float)(2 * (k % 2) - 1)) * scale;
            }
          } else {
            for (int k = 0; k < temp.extent(thirdDim); k++) {
              temp(i, j, k) = data[k] * scale;
            }
          }
        }

        delete[] data;
      }

    } break;
  }  // Switch

  // Cleanup
  fftwf_destroy_plan(plan);
}



void NDarray::gaussian_filter(Array<float, 2> &temp, int fsize) {
  int Nx = temp.length(firstDim);
  float *filter_bank = new float[Nx];

  // cout << "Gaussian filtering" << endl;
  // cout << " X = " << temp.length(firstDim) << endl;
  // cout << " Y = " << temp.length(secondDim) << endl;
  // cout << " F = " << fsize << endl;

  for (int j = 0; j < temp.length(secondDim); j++) {
    /***Store Old**/
    for (int i = 0; i < Nx; i++) {
      filter_bank[i] = temp(i, j);
      temp(i, j) = 0.0;
    }

    /**DO FILTER*/
    for (int i = 0; i < Nx; i++) {
      int start = max(i - fsize * 3, 0);
      int stop = min(i + fsize * 3, Nx - 1);

      for (int ii = start; ii <= stop; ii++) {
        float rad = (float)ii - (float)i;
        float val =
            (float)exp(-(rad * rad) / (2.0 * (float)fsize * (float)fsize));
        temp(i, j) += (filter_bank[ii] * val);
      }
    }
  }

  delete[] filter_bank;
}

void NDarray::gaussian_filter(Array<complex<float>, 2> &temp, int fsize) {
  int Nx = temp.length(firstDim);
  complex<float> *filter_bank = new complex<float>[Nx];

  // cout << "Gaussian filtering" << endl;
  // cout << " X = " << temp.length(firstDim) << endl;
  // cout << " Y = " << temp.length(secondDim) << endl;
  // cout << " F = " << fsize << endl;

  for (int j = 0; j < temp.length(secondDim); j++) {
    /***Store Old**/
    for (int i = 0; i < Nx; i++) {
      filter_bank[i] = temp(i, j);
      temp(i, j) = complex<float>(0.0, 0.0);
    }

    /**DO FILTER*/
    for (int i = 0; i < Nx; i++) {
      int start = max(i - fsize * 3, 0);
      int stop = min(i + fsize * 3, Nx - 1);

      for (int ii = start; ii <= stop; ii++) {
        float rad = (float)ii - (float)i;
        float val =
            (float)exp(-(rad * rad) / (2.0 * (float)fsize * (float)fsize));
        temp(i, j) += (filter_bank[ii] * val);
      }
    }
  }

  delete[] filter_bank;
}
