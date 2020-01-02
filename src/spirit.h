#pragma once

#include <iostream>

#include <fftw3.h>
#include <omp.h>
#include <armadillo>
#include <cmath>
#include "ArrayTemplates.hpp"
#include "tictoc.hpp"

enum { SP_SQUARE, SP_CIRCLE };
enum { SP_TIK, SP_TSVD };
enum { SP_COIL_PHASE, SP_SMOOTH };

void fftshift3(NDarray::Array<complex<float>, 3> &);
void ifft3(NDarray::Array<complex<float>, 3> &);

class SPIRIT {
 public:
  NDarray::Array<complex<float>, 5> k;
  NDarray::Array<complex<float>, 5> im;

  // enum flags
  int shape;

  int mapshrink;
  float mapthresh;

  // Size of the kernel (2*kr+1 for squares)
  int krx;
  int kry;
  int krz;

  float kr_f;
  float krx_f;
  float kry_f;
  float krz_f;
  float svd_thresh;
  int cr;
  int crx;
  int cry;
  int crz;

  int rcxres;
  int rcyres;
  int rczres;
  int debug;

  int ncoils;
  int nV;
  int phase_type;

  static void help_message(void);

  SPIRIT();
  void read_commandline(int numarg, char **pstring);
  void init(int, int, int, int);

  void generateEigenCoils(
      NDarray::Array<NDarray::Array<complex<float>, 3>, 1> &,
      NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &);
  void calibrate_ellipsoid(
      NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &);
  void prep();

  void getcoils(NDarray::Array<NDarray::Array<complex<float>, 3>, 1> &);
  void phase_correct(NDarray::Array<NDarray::Array<complex<float>, 3>, 1> &);
};
