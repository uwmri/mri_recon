#ifndef hVORDCF
#define hVORDCF

#include <omp.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <armadillo>

#include "ArrayTemplates.hpp"

/**
 *  @brief
 *  Class to perform Voroni Diagram of MRI Data
 *
 *  Usage Example:
 *  @code
 *  	include "voroni_dcf.h"

 *	@endcode
 */
class VORONOI_DCF {
 public:
  enum KShape { SPHERE, CYLINDER, CUBE };
  static void vor_dcf(NDarray::Array<float, 3> &, NDarray::Array<float, 3> &,
                      NDarray::Array<float, 3> &, NDarray::Array<float, 3> &,
                      KShape);

  static void vor_sphere(NDarray::Array<float, 1> &, NDarray::Array<float, 1> &,
                         NDarray::Array<float, 1> &,
                         NDarray::Array<float, 1> &);

 private:
};

#endif
