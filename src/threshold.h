#ifndef hTHRESHLIB
#define hTHRESHLIB

#include <omp.h>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

#include "ArrayTemplates.hpp"
#include "wavelet3D.h"

// Thresholding methods
enum { TH_NONE,
       TH_FIXED,
       TH_FRACTION,
       TH_VISU,
       TH_BAYES,
       TH_SURE };
enum { NOISE_GLOBAL,
       NOISE_FRAME,
       NOISE_FIRST,
       NOISE_LAST };

class THRESHOLD {
 public:
  THRESHOLD();
  THRESHOLD(int numarg, const char **pstring);
  bool soft;   // soft/hard thresholding
  bool thapp;  // threshold approx band
  bool temporal;
  double thresh;
  bool group_complex;
  float imaginary_scale;
  int waveL;
  bool VERBOSE;
  float global_threshold;
  NDarray::Array<float, 5> subband_threshold;

  NDarray::Array<float, 2> independent_thresholds;
  bool independent;

  float noise;
  float noise_scale;
  int noise_est_type;
  int threshold_type;
  char th_type[20];
  char th_mode[20];
  static void help_message(void);

  // Get threshold value from coeffitients
  void update_threshold(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &Coef, WAVELET3D &wave, float);

  // Chooses the thresholding method
  void exec_threshold(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &Coef, WAVELET3D &wave);

  void setThresholdMethod(int ith_type);
  int getThresholdMethod();

  void setTemporalThresholding(bool flag);  // to skip the first frame in processing
  bool getTemporalThresholding();

 private:
  // Execute thresholding method
  void exec_multibandthreshold(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &Coef, WAVELET3D &wave);
  void thresholding(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &Coef, float value);
  void thresholding(NDarray::Array<complex<float>, 3> &Coef, float value);

  // Get Threshold
  float get_threshold(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &Coef, float);
  float get_threshold(NDarray::Array<complex<float>, 3> &Coef, float);
  float sure_cost(NDarray::Array<complex<float>, 3> &Coef, float thresh);

  void get_visuthreshold(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &Coef, WAVELET3D &wave);
  void get_bayesthreshold(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &Coef, WAVELET3D &wave);
  void get_surethreshold(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &Coef, WAVELET3D &wave);
  float get_bayesthreshold_subband(NDarray::Array<complex<float>, 3> &XX);
  float get_surethreshold_subband(NDarray::Array<complex<float>, 3> &XX);

  void robust_noise_estimate(NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &Coef, WAVELET3D &wave);
};

#endif
