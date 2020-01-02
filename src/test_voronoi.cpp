#include "recon_lib.h"

using namespace NDarray;
using namespace std;

int main(void) {
  int N = 100000;
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
    Kx(pos) =
        sqrt(1 - Kz(pos) * Kz(pos)) * cos(2.0 * 3.14159265358979323846 * mphi2);
    Ky(pos) =
        sqrt(1 - Kz(pos) * Kz(pos)) * sin(2.0 * 3.14159265358979323846 * mphi2);
  }

  VORONOI_DCF::vor_sphere(Kw, Kx, Ky, Kz);

  ArrayWrite(Kx, "Kx.dat");
  ArrayWrite(Ky, "Ky.dat");
  ArrayWrite(Kz, "Kz.dat");
  ArrayWrite(Kw, "Kw.dat");

  return 0;
}
