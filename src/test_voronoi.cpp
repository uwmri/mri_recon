#include "recon_lib.h"
#include "voro++/voro++.hh"

using namespace NDarray;
using namespace std;
using namespace voro;

void test_radial_voronoi(void) {
  // Fake radial trajectory
  int xres = 256;
  int shots = 100;
  Array<float, 3> Kx(xres, shots, 1, ColumnMajorArray<3>());
  Array<float, 3> Ky(xres, shots, 1, ColumnMajorArray<3>());
  Array<float, 3> Kz(xres, shots, 1, ColumnMajorArray<3>());
  Array<float, 3> Kw(xres, shots, 1, ColumnMajorArray<3>());

  int count = 0;
  for (int shot = 0; shot < shots; shot++) {
    for (int x = 0; x < xres; x++) {
      float angle = (float)shot * (3.14159265358979323846 * (sqrt(5) - 1.0));
      float radius = ((float)x - (float)(xres / 2));
      Kx(x, shot, 0) = radius * cos(angle);
      Ky(x, shot, 0) = radius * sin(angle);
      Kz(x, shot, 0) = 0.0;
      Kw(x, shot, 0) = 0.0;
    }
  }

  VORONOI_DCF::vor_dcf(Kw, Kx, Ky, Kz, VORONOI_DCF::CYLINDER);
}

void test_sphere_voronoi(void) {
  int N = 10000;
  Array<float, 1> Kx(N);
  Array<float, 1> Ky(N);
  Array<float, 1> Kz(N);
  Array<float, 1> Kw(N);

  for (int pos = 0; pos < N; pos++) {
    double proj = (double)pos;

    /*  Just multiply. Note this will be off a bit due to rounding error
     * (1e-11)*/
    double mphi1 = (double)proj * 0.465571231876768;
    double mphi2 = (double)proj * 0.682327803828019;

    mphi1 = mphi1 - ((double)((int)mphi1));
    mphi2 = mphi2 - ((double)((int)mphi2));

    Kz(pos) = 2.0 * mphi1 - 1.0;
    Kx(pos) = sqrt(1 - Kz(pos) * Kz(pos)) * cos(2.0 * 3.14159265358979323846 * mphi2);
    Ky(pos) = sqrt(1 - Kz(pos) * Kz(pos)) * sin(2.0 * 3.14159265358979323846 * mphi2);
  }

  VORONOI_DCF::vor_sphere(Kw, Kx, Ky, Kz);
}

int main(void) {
  test_radial_voronoi();
  test_sphere_voronoi();

  return 0;
}
