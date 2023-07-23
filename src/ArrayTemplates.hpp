#ifndef hARRAYTEMPLATE
#define hARRAYTEMPLATE

// Switching to Blitz Based Arrays
#include <blitz/array.h>
#include <fftw3.h>
#include <omp.h>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

// It should be ok at populate namespace with complex
using std::complex;

namespace NDarray {

using namespace blitz;

// FFT Libraries complex Float
void fftshift(Array<complex<float>, 3> &temp);
void fftshift(Array<complex<float>, 4> &temp);
void fft(Array<complex<float>, 3> &temp);
void fft(Array<complex<float>, 4> &temp);
void ifft(Array<complex<float>, 3> &temp);
void ifft(Array<complex<float>, 4> &temp);

void fft(Array<complex<float>, 3> &temp, int);
void ifft(Array<complex<float>, 3> &temp, int);
void fft3(Array<complex<float>, 3> &temp, int, int, bool);

void gaussian_filter(Array<float, 2> &, int);
void gaussian_filter(Array<complex<float>, 2> &, int);
complex<float> conj_sum(Array<complex<float>, 3> P, Array<complex<float>, 3> R);
void gaussian_blur(Array<complex<float>, 3> &In, float sigX, float sigY, float sigZ);
void gaussian_blur(Array<float, 3> &In, float sigX, float sigY, float sigZ);

inline void endian_swap(int &x) {
  x = (x << 24 & 0xFF000000) | (x << 8 & 0x00FF0000) | (x >> 8 & 0x0000FF00) |
      (x >> 24 & 0x000000FF);
}

template <typename T, int N>
void ArrayRead(blitz::Array<T, N> &temp, const char *name) {
  FILE *fid;
  if ((fid = fopen(name, "r")) == NULL) {
    cout << "Array:Can't Open " << name << endl;
    cout << "Exiting" << endl;
    exit(1);
  } else {
    int j;
    if ((j = fread(temp.data(), sizeof(T), temp.numElements(), fid)) !=
        (int)(temp.numElements())) {
      cout << "Array3:Not enough data: only read " << j << "points of"
           << temp.numElements() << endl;
      exit(1);
    }
    fclose(fid);
  }
}

template <typename T, const int N_rank>
void ArrayWriteMagAppend(Array<complex<T>, N_rank> &temp, const char *name) {
  ofstream ofs(name, ios_base::binary | ios_base::app);

  for (typename Array<complex<T>, N_rank>::iterator miter = temp.begin();
       miter != temp.end(); miter++) {
    T val = abs(*miter);
    ofs.write((char *)&val, sizeof(T));
  }
}

template <typename T, const int N_rank>
void ArrayWriteAppend(Array<T, N_rank> &temp, const char *name) {
  ofstream ofs(name, ios_base::binary | ios_base::app);

  for (typename Array<T, N_rank>::iterator miter = temp.begin();
       miter != temp.end(); miter++) {
    T val = *miter;
    ofs.write((char *)&val, sizeof(T));
  }
}

template <typename T, const int N_rank>
void ArrayWriteAppendAsComplexFloat(Array<T, N_rank> &temp, const char *name) {
  ofstream ofs(name, ios_base::binary | ios_base::app);

  for (typename Array<T, N_rank>::iterator miter = temp.begin();
       miter != temp.end(); miter++) {
    complex<float> val = *miter;
    ofs.write((char *)&val, sizeof(complex<float>));
  }
}

template <typename T>
void ArrayWriteAppendZeros(int number, const char *name) {
  ofstream ofs(name, ios_base::binary | ios_base::app);

  for (int i = 0; i < number; i++) {
    T val = (T)0.0;  // Zero is zero for all types!
    ofs.write((char *)&val, sizeof(T));
  }
}

template <typename T, const int N_rank>
void ArrayWrite(Array<T, N_rank> &temp, const char *name) {
  remove(name);  // This speeds things up considerably
  ofstream ofs(name, ios_base::binary);

  for (typename Array<T, N_rank>::iterator miter = temp.begin();
       miter != temp.end(); miter++) {
    T val = *miter;
    ofs.write((char *)&val, sizeof(T));
  }
}

template <typename T, const int N_rank>
void ArrayWriteMag(Array<complex<T>, N_rank> &temp, const char *name) {
  remove(name);  // This speeds things up considerably
  ofstream ofs(name, ios_base::binary);
  for (typename Array<complex<T>, N_rank>::iterator miter = temp.begin();
       miter != temp.end(); miter++) {
    T val = abs(*miter);
    ofs.write((char *)&val, sizeof(T));
  }
}

template <typename T, const int N_rank>
double ArrayEnergy(Array<complex<T>, N_rank> &temp) {
  double EE = 0;
  for (typename Array<complex<T>, N_rank>::iterator miter = temp.begin();
       miter != temp.end(); miter++) {
    EE += (double)(norm(*miter));
  }
  return (EE);
}

template <typename T, const int N_rank, const int M_rank>
double ArrayEnergy(Array<Array<complex<T>, N_rank>, M_rank> &temp) {
  double EE = 0;
  for (typename Array<Array<complex<T>, N_rank>, M_rank>::iterator miter =
           temp.begin();
       miter != temp.end(); miter++) {
    EE += ArrayEnergy(*miter);
  }
  return (EE);
}

template <typename T, const int N_rank>
void ArrayWritePhase(Array<complex<T>, N_rank> &temp, const char *name) {
  ofstream ofs(name, ios_base::binary);

  for (typename Array<complex<T>, N_rank>::iterator miter = temp.begin();
       miter != temp.end(); miter++) {
    T val = arg(*miter);
    ofs.write((char *)&val, sizeof(T));
  }
}

template <typename T, const int N_rank>
void ArrayWritePhaseAppend(Array<complex<T>, N_rank> &temp, const char *name) {
  ofstream ofs(name, ios_base::binary | ios_base::app);

  for (typename Array<complex<T>, N_rank>::iterator miter = temp.begin();
       miter != temp.end(); miter++) {
    T val = arg(*miter);
    ofs.write((char *)&val, sizeof(T));
  }
}

template <typename T, const int N_rank, const int M_rank>
void WriteCFL(Array<Array<T, N_rank>, M_rank> &temp, const char *name,
              int pad_dim1 = 0) {
  // Create names for header and binary
  char name_hdr[2048];
  char name_bin[2048];

  sprintf(name_hdr, "%s.hdr", name);
  sprintf(name_bin, "%s.cfl", name);

  // Now determine the size of the ND array
  TinyVector<int, M_rank> Dim1 = temp.shape();

  // Loop over all sub-arrays to get the max
  typename Array<Array<T, N_rank>, M_rank>::iterator miter = temp.begin();

  TinyVector<int, N_rank> Dim2 = (*miter).shape();
  for (; miter != temp.end(); miter++) {
    TinyVector<int, N_rank> Dim2_temp = (*miter).shape();
    for (int i = 0; i < N_rank; i++) {
      Dim2(i) = max(Dim2(i), Dim2_temp(i));
    }
  }

  // Combine the Size to get the total size
  Array<int, 1> Dim(N_rank + M_rank + pad_dim1);
  int count = 0;
  for (int i = 0; i < pad_dim1; i++) {
    Dim(count) = 1;
    count++;
  }
  for (int i = 0; i < N_rank; i++) {
    Dim(count) = Dim2(i);
    count++;
  }
  for (int i = 0; i < M_rank; i++) {
    Dim(count) = Dim1(i);
    count++;
  }

  // Debug info (temp)
  cout << "File Name = " << name_bin << endl;
  cout << "Header Name = " << name_hdr << endl;
  cout << "Output Size = " << Dim << endl;

  // Write the header
  FILE *fid;
  fid = fopen(name_hdr, "w");
  fprintf(fid, "# Dimensions\n");
  for (int i = 0; i < M_rank + N_rank; i++) {
    fprintf(fid, "%d ", Dim(i));
  }
  if ((N_rank + M_rank) < 5) {
    for (int i = 0; i < (5 - (M_rank + N_rank)); i++) {
      fprintf(fid, "%d ", 1);
    }
  }
  fclose(fid);

  // Write the binary data
  remove(name_bin);
  int inner_size = product(Dim2);
  for (typename Array<Array<T, N_rank>, M_rank>::iterator miter = temp.begin();
       miter != temp.end(); miter++) {
    ArrayWriteAppendAsComplexFloat((*miter), name_bin);

    // Bart doesn't support container sizes. Just pad with zeros
    if ((*miter).numElements() < inner_size) {
      ArrayWriteAppendZeros<complex<float> >(
          inner_size - (*miter).numElements(), name_bin);
    }
  }
}

template <typename T, const int N_rank>
void WriteCFL(Array<T, N_rank> &temp, const char *name) {
  // Create names for header and binary
  char name_hdr[2048];
  char name_bin[2048];

  sprintf(name_hdr, "%s.hdr", name);
  sprintf(name_bin, "%s.cfl", name);

  // Now determine the size of the ND array
  TinyVector<int, N_rank> Dim = temp.shape();

  // Debug info (temp)
  cout << "File Name = " << name_bin << endl;
  cout << "Header Name = " << name_hdr << endl;
  cout << "Output Size = " << Dim << endl;

  // Write the header
  FILE *fid;
  fid = fopen(name_hdr, "w");
  fprintf(fid, "# Dimensions\n");
  for (int i = 0; i < N_rank; i++) {
    fprintf(fid, "%d ", Dim(i));
  }
  if ((N_rank) < 5) {
    for (int i = 0; i < (5 - (N_rank)); i++) {
      fprintf(fid, "%d ", 1);
    }
  }
  fclose(fid);

  // Write the binary data
  remove(name_bin);
  ArrayWriteAppendAsComplexFloat(temp, name_bin);
}

template <typename T, const int N_rank, const int M_rank>
void WriteCFL_triplet(Array<Array<T, N_rank>, M_rank> &temp,
                      Array<Array<T, N_rank>, M_rank> &temp2,
                      Array<Array<T, N_rank>, M_rank> &temp3,
                      const char *name) {
  // Create names for header and binary
  char name_hdr[2048];
  char name_bin[2048];

  sprintf(name_hdr, "%s.hdr", name);
  sprintf(name_bin, "%s.cfl", name);

  // Now determine the size of the ND array
  TinyVector<int, M_rank> Dim1 = temp.shape();

  // Loop over all sub-arrays to get the max
  typename Array<Array<T, N_rank>, M_rank>::iterator miter = temp.begin();

  TinyVector<int, N_rank> Dim2 = (*miter).shape();
  for (; miter != temp.end(); miter++) {
    TinyVector<int, N_rank> Dim2_temp = (*miter).shape();
    for (int i = 0; i < N_rank; i++) {
      Dim2(i) = max(Dim2(i), Dim2_temp(i));
    }
  }

  // Combine the Size to get the total size
  TinyVector<int, N_rank + M_rank> Dim;
  for (int i = 0; i < N_rank; i++) {
    Dim(i) = Dim2(i);
  }
  for (int i = 0; i < M_rank; i++) {
    Dim(i + N_rank) = Dim1(i);
  }

  TinyVector<int, N_rank + M_rank + 1> DimTriple;
  DimTriple(0) = 3;
  for (int i = 0; i < (N_rank + M_rank); i++) {
    DimTriple(i + 1) = Dim(i);
  }

  // Debug info (temp)
  cout << "File Name = " << name_bin << endl;
  cout << "Header Name = " << name_hdr << endl;
  cout << "Output Size = " << Dim << endl;

  // Write the header
  FILE *fid;
  fid = fopen(name_hdr, "w");
  fprintf(fid, "# Dimensions\n");
  for (int i = 0; i < M_rank + N_rank + 1; i++) {
    fprintf(fid, "%d ", DimTriple(i));
  }
  if ((N_rank + M_rank + 1) < 5) {
    for (int i = 0; i < (5 - (M_rank + N_rank + 1)); i++) {
      fprintf(fid, "%d ", 1);
    }
  }
  fclose(fid);

  // Write the binary data
  remove(name_bin);
  int inner_size = product(Dim2);
  ofstream ofs(name_bin, ios_base::binary | ios_base::app);

  // Three input arrays
  typename Array<Array<T, N_rank>, M_rank>::iterator miter1 = temp.begin();
  typename Array<Array<T, N_rank>, M_rank>::iterator miter2 = temp2.begin();
  typename Array<Array<T, N_rank>, M_rank>::iterator miter3 = temp3.begin();
  for (; (miter1 != temp.end()); miter1++, miter2++, miter3++) {
    // Pointers to the the inner array
    typename Array<T, N_rank>::iterator niter1 = (*miter1).begin();
    typename Array<T, N_rank>::iterator niter2 = (*miter2).begin();
    typename Array<T, N_rank>::iterator niter3 = (*miter3).begin();

    for (; (niter1 != (*miter1).end()); niter1++, niter2++, niter3++) {
      {
        complex<float> val = *niter1;
        ofs.write((char *)&val, sizeof(complex<float>));
      }
      {
        complex<float> val = *niter2;
        ofs.write((char *)&val, sizeof(complex<float>));
      }
      {
        complex<float> val = *niter3;
        ofs.write((char *)&val, sizeof(complex<float>));
      }
    }

    // Bart doesn't support container sizes. Just pad with zeros
    if ((*miter1).numElements() < inner_size) {
      ArrayWriteAppendZeros<complex<float> >(
          3 * (inner_size - (*miter1).numElements()), name_bin);
    }
  }
}

template <typename T>
Array<Array<T, 3>, 1> Alloc4DContainer(int x, int y, int z, int t) {
  Array<Array<T, 3>, 1> temp;
  temp.setStorage(ColumnMajorArray<1>());
  temp.resize(t);

  for (typename Array<Array<T, 3>, 1>::iterator miter = temp.begin();
       miter != temp.end(); miter++) {
    (*miter).setStorage(ColumnMajorArray<3>());
    (*miter).resize(x, y, z);
    (*miter) = (T)0;
  }
  return (temp);
}

template <typename T>
Array<Array<T, 3>, 3> Alloc6DContainer(int x, int y, int z, int d1, int d2,
                                       int d3) {
  Array<Array<T, 3>, 3> temp;
  temp.setStorage(ColumnMajorArray<3>());
  temp.resize(d1, d2, d3);

  for (typename Array<Array<T, 3>, 3>::iterator miter = temp.begin();
       miter != temp.end(); miter++) {
    (*miter).setStorage(ColumnMajorArray<3>());
    (*miter).resize(x, y, z);
    (*miter) = (T)0;
  }
  return (temp);
}

template <typename T>
Array<Array<T, 3>, 2> Alloc5DContainer(int x, int y, int z, int d1, int d2) {
  Array<Array<T, 3>, 2> temp;
  temp.setStorage(ColumnMajorArray<2>());
  temp.resize(d1, d2);

  for (int i = 0; i < d1; i++) {
    for (int j = 0; j < d2; j++) {
      temp(i, j).setStorage(ColumnMajorArray<3>());
      temp(i, j).resize(x, y, z);
      temp(i, j) = (T)0;
    }
  }

  /* This leads to errors for no real reason
  for( typename Array< Array<T,3>,2>::iterator miter=temp.begin();   miter
  !=temp.end(); miter++){
          (*miter).setStorage(ColumnMajorArray<3>());
          (*miter).resize(x,y,z);
          (*miter)= (T )0;
  }
  */
  return (temp);
}

double Dmax(const Array<Array<double, 2>, 1> &A);
double Dmin(const Array<Array<double, 2>, 1> &A);

void nested_workaround(long index, int *N, int *idx, int total);

}  // namespace NDarray

#endif