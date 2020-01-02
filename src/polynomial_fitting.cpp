#include "polynomial_fitting.h"
using namespace NDarray;
using namespace arma;

/*Code for fitting image to polynomial**/
double POLYFIT::poly3d(float x, float y, float z) {
  double value = 0;
  for (int pos = 0; pos < number; pos++) {
    value +=
        alpha(pos) * powf(x, px(pos)) * powf(y, py(pos)) * powf(z, pz(pos));
  }
  return value;
}

double POLYFIT::poly2d(float x, float y) {
  double value = 0;
  for (int pos = 0; pos < number; pos++) {
    value += alpha(pos) * powf(x, px(pos)) * powf(y, py(pos));
  }
  return value;
}

void POLYFIT::poly_fitting3d(Array<float, 3> &back_mag /*Binary Matrix*/,
                             Array<float, 3> &image, int fit_order) {
  number = 0;
  for (int x = 0; x <= fit_order; x++) {
    for (int y = 0; y <= fit_order; y++) {
      for (int z = 0; z <= fit_order; z++) {
        if ((x + y + z) <= fit_order) {
          number++;
        }
      }
    }
  }

  // Calculate Coef
  px = zeros<vec>(number);
  py = zeros<vec>(number);
  pz = zeros<vec>(number);
  int pos = 0;
  for (int x = 0; x <= fit_order; x++) {
    for (int y = 0; y <= fit_order; y++) {
      for (int z = 0; z <= fit_order; z++) {
        if ((x + y + z) <= fit_order) {
          px(pos) = (double)x;
          py(pos) = (double)y;
          pz(pos) = (double)z;
          pos++;
        }
      }
    }
  }

  /***2 Build Ah * A * r = Ah * r */
  arma::mat AA = zeros<mat>(number, number);
  arma::vec BB = zeros<vec>(number);
  arma::vec store = zeros<vec>(number);
  arma::vec val = zeros<vec>(number);

  double cx = (double)image.length(firstDim) / 2.0;
  double cy = (double)image.length(secondDim) / 2.0;
  double cz = (double)image.length(thirdDim) / 2.0;

  /*Perform SS and put into storage matrix*/
  for (int k = 0; k < image.length(thirdDim); k++) {
    for (int j = 0; j < image.length(secondDim); j++) {
      for (int i = 0; i < image.length(firstDim); i++) {
        if (back_mag(i, j, k) > 0) {
          double value = (double)image(i, j, k);
          double x = (double)i - cx;
          double y = (double)j - cy;
          double z = (double)k - cz;

          for (int pos = 0; pos < number; pos++) {
            store(pos) = powf(x, px(pos)) * powf(y, py(pos)) * powf(z, pz(pos));
          }

          // Matrix Allocation
          for (int im = 0; im < number; im++) {
            BB(im) += store(im) * value;
            for (int jm = 0; jm < number; jm++) {
              AA(im, jm) += store(im) * store(jm);
            }
          }
        }
      }
    }
  } /*x,y,z,thresh*/

  // Solve and assign
  alpha = zeros<vec>(number);
  if (number == 1) {
    alpha(0) = BB(0) / AA(0, 0);
  } else {
    arma::mat XX = solve(AA, BB);
    for (int index = 0; index < number; index++) {
      alpha(index) = XX(index);
    }
  }
}

void POLYFIT::poly_subtract3d(Array<float, 3> &image) {
  float cx = (float)image.length(firstDim) / 2.0;
  float cy = (float)image.length(secondDim) / 2.0;
  float cz = (float)image.length(thirdDim) / 2.0;
/***Do Correction*/
#pragma omp parallel for
  for (int k = 0; k < image.length(thirdDim); k++) {
    for (int j = 0; j < image.length(secondDim); j++) {
      for (int i = 0; i < image.length(firstDim); i++) {
        image(i, j, k) -= poly3d((float)i - cx, (float)j - cy, (float)k - cz);
      }
    }
  }
}
