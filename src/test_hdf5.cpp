#include "recon_lib.h"

using namespace NDarray;
using namespace std;

int main(void) {
  // Write Section
  {
    HDF5 file = HDF5("TEST.h5", "w");

    Array<float, 3> TEMP(2, 3, 4, ColumnMajorArray<3>());
    Array<complex<float>, 3> TEMPC(4, 3, 2, ColumnMajorArray<3>());
    Array<short, 3> TEMPS(2, 3, 4, ColumnMajorArray<3>());

    for (int k = 0; k < TEMPC.length(thirdDim); k++) {
      for (int j = 0; j < TEMPC.length(secondDim); j++) {
        for (int i = 0; i < TEMPC.length(firstDim); i++) {
          TEMPC(i, j, k) = complex<float>(10 * k + j, i);
        }
      }
    }

    for (int k = 0; k < TEMP.length(thirdDim); k++) {
      for (int j = 0; j < TEMP.length(secondDim); j++) {
        for (int i = 0; i < TEMP.length(firstDim); i++) {
          TEMP(i, j, k) = 100 * k + 10 * j + i;
          TEMPS(i, j, k) = 100 * k + 10 * j + i - 2;
        }
      }
    }

    cout << "Exporting Attributes" << endl;
    int A = 10;
    float B = 3.14156;
    string C = "Super";

    file.AddH5Scaler("Kdata", "A", A);
    file.AddH5Scaler("Kdata", "B", B);
    file.AddH5Char("Kdata", "C", C.c_str());
    file.AddH5Array("Kdata", "D", TEMP);
    file.AddH5Array("Kdata", "E", TEMPC);
    file.AddH5Array("Kdata", "F", TEMPS);
  }

  // Write Section
  {
    HDF5 file = HDF5("TEST.h5", "r");

    cout << "Reading Attributes" << endl;
    int A;
    float B;
    std::string C;
    file.ReadH5Scaler("Kdata", "A", &A);
    cout << "A = " << A << endl
         << flush;

    file.ReadH5Scaler("Kdata", "B", &B);
    cout << "B = " << B << endl
         << flush;

    file.ReadH5Char("Kdata", "C", C);
    cout << "C = " << C << endl
         << flush;

    Array<float, 3> TEMP(3, 3, 3);
    file.ReadH5Array("Kdata", "D", TEMP);

    for (int k = 0; k < TEMP.length(thirdDim); k++) {
      for (int j = 0; j < TEMP.length(secondDim); j++) {
        for (int i = 0; i < TEMP.length(firstDim); i++) {
          cout << "k=" << k << ",j=" << j << ",i=" << i << " = "
               << TEMP(i, j, k) << endl;
        }
      }
    }

    Array<complex<float>, 3> TEMPC;
    file.ReadH5Array("Kdata", "E", TEMPC);

    for (int k = 0; k < TEMPC.length(thirdDim); k++) {
      for (int j = 0; j < TEMPC.length(secondDim); j++) {
        for (int i = 0; i < TEMPC.length(firstDim); i++) {
          cout << "k=" << k << ",j=" << j << ",i=" << i << " = "
               << TEMPC(i, j, k) << endl;
        }
      }
    }
    
    Array<short, 3> TEMPF;
    file.ReadH5Array("Kdata", "F", TEMPF);

    for (int k = 0; k < TEMPF.length(thirdDim); k++) {
      for (int j = 0; j < TEMPF.length(secondDim); j++) {
        for (int i = 0; i < TEMPF.length(firstDim); i++) {
          cout << "k=" << k << ",j=" << j << ",i=" << i << " = "
               << TEMPF(i, j, k) << endl;
        }
      }
    }
    
  }

  return 0;
}
