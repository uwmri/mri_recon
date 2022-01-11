#ifndef hPOLYFIT
#define hPOLYFIT

/**Code to Fit image data to a polynomial*/

#include <omp.h>
#include <armadillo>

// For Wrapper (this is the actual recon + NDarray)
#include "ArrayTemplates.hpp"

class POLYFIT {
 public:
  /*Code for fitting image to polynomial**/
  double poly3d(float x, float y, float z);
  double poly2d(float x, float y);

  /* Fitting */
  void poly_fitting3d(NDarray::Array<float, 3> &back_mag /*Binary Matrix*/,
                      blitz::Array<float, 3> &image, int);

  // Return image
  NDarray::Array<float, 3> image(void);

  /* Operations */
  void poly_subtract3d(NDarray::Array<float, 3> &);
  void poly_divide(NDarray::Array<float, 3> &);
  void poly_multiply(NDarray::Array<float, 3> &);
  void poly_add(NDarray::Array<float, 3> &);

 private:
  int poly_exists;
  int number;
  arma::vec alpha;
  arma::Col<int> px;
  arma::Col<int> py;
  arma::Col<int> pz;

  NDarray::TinyVector<int, 3> image_res;
  NDarray::TinyVector<int, 3> image_center;
};
#endif
