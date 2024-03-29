#ifndef hGATINGLIB
#define hGATINGLIB

#include <omp.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <armadillo>

#include "ArrayTemplates.hpp"
#include "gridFFT.h"
#include "mri_data.h"
#include "tictoc.hpp"

// View sharing modes
#define VS_NONE 0
#define VS_SLIDING 1
#define VS_TORNADO 2

class GATING {
 public:
  enum ViewshareType { TORNADO,
                       NONE,
                       HIST_MODE };
  enum TornadoType { FLAT,
                     RADIAL,
                     VIPR };
  enum WeightType { ITERATIVE,
                    NON_ITERATIVE };
  enum FrameType { COMPOSITE,
                   TIME_FRAME };
  enum GateType { GATE_NONE,
                  RETRO_ECG,
                  ECG,
                  RESP,
                  TIME,
                  PREP };
  enum RespGateType { RESP_NONE,
                      RESP_THRESH,
                      RESP_WEIGHT,
                      RESP_HARD };

  GATING();
  GATING(int numarg, const char **pstring);
  void init(const MRI_DATA &data, int);
  void init_resp_gating(const MRI_DATA &data);
  void init_time_resolved(const MRI_DATA &data, int);
  int number_of_frames(void);

  // Tornado Filter Parameters
  int wdth_low;   // k=0 width
  int wdth_high;  // k=kmax width
  float kmax;
  TornadoType tornado_shape;  // kr^2 vs kr

  // Scaling for Waveform
  double scale_time;
  double offset_time;
  double actual_temporal_resolution;

  int recon_frames;

  // Control Gating Method
  ViewshareType vs_type;
  GateType gate_type;

  NDarray::Array<NDarray::Array<double, 2>, 1> gate_times;
  NDarray::Array<NDarray::Array<double, 2>, 1> resp_weight;

  // Control of Retrospective Respiratory Gating
  RespGateType resp_gate_type;
  int correct_resp_drift;
  float resp_gate_efficiency;
  float resp_gate_weight;

  // Respiratory Signal
  enum RespGateSignal { BELLOWS,
                        DC_DATA };
  RespGateSignal resp_gate_signal;
  float resp_sign;

  // Frame Centers
  NDarray::Array<double, 1> gate_frames;
  float gate_min_quantile;
  float gate_max_quantile;

  // Function Calls
  static void help_message(void);
  void weight_data(NDarray::Array<float, 3> &Tw, int e,
                   const NDarray::Array<float, 3> &kx,
                   const NDarray::Array<float, 3> &ky,
                   const NDarray::Array<float, 3> &kz, int t, WeightType,
                   FrameType);
  float temporal_resolution(void);
  void hist_weight(NDarray::Array<float, 3> &Tw, int e, int t);
  void tornado_weight(NDarray::Array<float, 3> &Tw, int e,
                      const NDarray::Array<float, 3> &kx,
                      const NDarray::Array<float, 3> &ky,
                      const NDarray::Array<float, 3> &kz, int t, WeightType);
  void filter_resp(const MRI_DATA &data);
  NDarray::Array<NDarray::Array<complex<float>, 2>, 1> combine_kspace_channels(const NDarray::Array<NDarray::Array<complex<float>, 2>, 2> &kdata_gating);

 private:
};

#endif
