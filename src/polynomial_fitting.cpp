#include "polynomial_fitting.h"
using namespace NDarray;
using namespace arma;

/*Code for fitting image to polynomial**/
double POLYFIT::poly3d(float x, float y, float z) {
  if (!poly_exists) {
    return (0.0);
  }

  double value = 0;
  x = (x - (double)image_center(0)) / (double)image_res(0);
  y = (y - (double)image_center(1)) / (double)image_res(1);
  z = (z - (double)image_center(2)) / (double)image_res(2);

  for (int pos = 0; pos < number; pos++) {
    value += alpha(pos) * pow(x, px(pos)) * pow(y, py(pos)) * pow(z, pz(pos));
  }
  return value;
}

double POLYFIT::poly2d(float x, float y) {
  if (!poly_exists) {
    return (0.0);
  }

  double value = 0;

  x = (x - (double)image_center(0)) / (double)image_res(0);
  y = (y - (double)image_center(1)) / (double)image_res(1);

  for (int pos = 0; pos < number; pos++) {
    value += alpha(pos) * powf(x, px(pos)) * powf(y, py(pos));
  }
  return value;
}

void POLYFIT::poly_fitting3d(Array<float, 3> &back_mag /*Binary Matrix*/, Array<float, 3> &image, int fit_order) {
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

  // Store in the image
  image_res = image.shape();
  image_center(0) = image_res(0) / 2;
  image_center(1) = image_res(1) / 2;
  image_center(2) = image_res(2) / 2;

  // Calculate Coef
  px = zeros<Col<int> >(number);
  py = zeros<Col<int> >(number);
  pz = zeros<Col<int> >(number);
  int pos = 0;
  for (int x = 0; x <= fit_order; x++) {
    for (int y = 0; y <= fit_order; y++) {
      for (int z = 0; z <= fit_order; z++) {
        if ((x + y + z) <= fit_order) {
          px(pos) = x;
          py(pos) = y;
          pz(pos) = z;
          pos++;
        }
      }
    }
  }

  cout << "Polynomial Coef" << std::endl;
  for (int pos = 0; pos < number; pos++) {
    cout << "Coefs" << px(pos) << "," << py(pos) << "," << pz(pos) << std::endl;
  }

  // Get the array
  arma::mat A = zeros<arma::mat>(number, number);
  arma::vec B = zeros<arma::vec>(number);

  /*Perform SS and put into storage matrix*/
  int background_fit_subsample = 1;
  for (int k = 0; k < image.length(thirdDim); k += background_fit_subsample) {
    for (int j = 0; j < image.length(secondDim); j += background_fit_subsample) {
      for (int i = 0; i < image.length(firstDim); i += background_fit_subsample) {
        if (back_mag(i, j, k) > 0) {
          // Get the array coordinates
          vec S = zeros<vec>(number);
          for (int pos = 0; pos < number; pos++) {
            double x = ((double)i - (double)image_center(0)) / (double)image_res(0);
            double y = ((double)j - (double)image_center(1)) / (double)image_res(1);
            double z = ((double)k - (double)image_center(2)) / (double)image_res(2);

            double val_x = pow(x, px(pos));
            double val_y = pow(y, py(pos));
            double val_z = pow(z, pz(pos));
            S(pos) = (val_x * val_y * val_z);
          }

          //cout << S << endl;

          // Add to array
          B += ((double)image(i, j, k) * S);
          A += S * S.t();
        }
      }
    }
  } /*x,y,z,thresh*/

  cout << "A" << A << endl;
  cout << "B" << B << endl;

  // Solve and assign
  alpha = zeros<vec>(number);
  if (number > 1) {
    cout << "Solving for Polynomial Coef " << endl
         << flush;
    arma::mat XX = arma::solve(A, B);
    arma::mat DIFF = A * XX - B;
    double error = arma::norm(DIFF, 1);
    std::cout << "Mean Error = " << error << std::endl;
    for (int index = 0; index < number; index++) {
      alpha(index) = XX(index);
      cout << "Poly Coef = " << alpha(index) << std::endl;
    }
  }

  poly_exists = true;
}

void POLYFIT::poly_subtract3d(Array<float, 3> &image) {
#pragma omp parallel for
  for (int k = 0; k < image.length(thirdDim); k++) {
    for (int j = 0; j < image.length(secondDim); j++) {
      for (int i = 0; i < image.length(firstDim); i++) {
        image(i, j, k) -= poly3d((float)i, (float)j, (float)k);
      }
    }
  }
}

void POLYFIT::poly_add(Array<float, 3> &image) {
#pragma omp parallel for
  for (int k = 0; k < image.length(thirdDim); k++) {
    for (int j = 0; j < image.length(secondDim); j++) {
      for (int i = 0; i < image.length(firstDim); i++) {
        image(i, j, k) += poly3d((float)i, (float)j, (float)k);
      }
    }
  }
}

void POLYFIT::poly_multiply(Array<float, 3> &image) {
#pragma omp parallel for
  for (int k = 0; k < image.length(thirdDim); k++) {
    for (int j = 0; j < image.length(secondDim); j++) {
      for (int i = 0; i < image.length(firstDim); i++) {
        image(i, j, k) *= poly3d((float)i, (float)j, (float)k);
      }
    }
  }
}

void POLYFIT::poly_divide(Array<float, 3> &image) {
#pragma omp parallel for
  for (int k = 0; k < image.length(thirdDim); k++) {
    for (int j = 0; j < image.length(secondDim); j++) {
      for (int i = 0; i < image.length(firstDim); i++) {
        image(i, j, k) /= poly3d((float)i, (float)j, (float)k);
      }
    }
  }
}

Array<float, 3> POLYFIT::image(void) {
  Array<float, 3> poly_image(image_res, ColumnMajorArray<3>());

#pragma omp parallel for
  for (int k = 0; k < image_res(thirdDim); k++) {
    for (int j = 0; j < image_res(secondDim); j++) {
      for (int i = 0; i < image_res(firstDim); i++) {
        poly_image(i, j, k) = poly3d((float)i, (float)j, (float)k);
      }
    }
  }

  return poly_image;
}
