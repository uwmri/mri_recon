#ifndef hL2REGLIB
#define hL2REGLIB

#include <iostream>

#include <omp.h>
#include <armadillo>
#include <cmath>
#include "ArrayTemplates.hpp"
#include "tictoc.hpp"

class L2REG {
 public:
  L2REG();
  L2REG(int numarg, const char **pstring);

  enum TransformType { NONE,
                       TV,
                       PHASE,
                       LOWRES };
  static void help_message(void);
  float lambda;
  float reg_scale;
  TransformType l2_type;
  NDarray::Array<complex<float>, 3> ZeroPad;
  int verbose;
  void regularize(NDarray::Array<complex<float>, 3> &,
                  NDarray::Array<complex<float>, 3> &);
  void set_scale(float, NDarray::Array<NDarray::Array<complex<float>, 3>, 2> &);
  void set_scale(float, NDarray::Array<complex<float>, 3> &);
};

#endif