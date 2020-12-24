/*
        This contains Commandline Interface to the Recon
 */
#include "recon_lib.h"

using namespace NDarray;
using namespace std;

int main(int argc, const char **argv) {
  // ------------------------------------
  // Initialize Recon - reading the command line
  // ------------------------------------
  RECON recon(argc, argv);
  if (recon.threads != -1) {
    omp_set_num_threads(recon.threads);
  }

  // ------------------------------------
  // Read Data -
  // ------------------------------------

  MRI_DATA data;
  cout << "----Read Data-----" << endl;
  switch (recon.data_type) {
    case (RECON::BENCHMARK): {
      // This runs a system benchmark

    } break;

    case (RECON::SIMULATE): {
      // Use completely made up data
      data.Num_Encodings = 1;
      data.Num_Coils = 1;
      data.Num_Frames = 1;
      int Num_Proj = 10000;
      int Kmax = 64;
      int Xres = Kmax * 2 + 1;
      float dt = 4e-6;

      // This initializes the encode data which could normally be assigned
      data.init_memory(Xres, Num_Proj, 1);

      // 3D Radial
      data.trajectory_type(0) = MRI_DATA::NONCARTESIAN;
      data.trajectory_type(1) = MRI_DATA::NONCARTESIAN;
      data.trajectory_type(2) = MRI_DATA::NONCARTESIAN;
      data.dft_needed(0) = true;
      data.dft_needed(1) = true;
      data.dft_needed(2) = true;
      data.xres = 2 * Kmax;
      data.yres = 2 * Kmax;
      data.zres = 2 * Kmax;
      data.xfov = 1;
      data.yfov = 1;
      data.zfov = 1;

      // Generate
      std::cout << "Create " << std::endl
                << std::flush;
      float slr_freq = sqrt((float)Num_Proj * 3.14159);
      for (int proj = 0; proj < Num_Proj; proj++) {
        float z = (2.0 * (float)proj - Num_Proj + 1) / ((float)Num_Proj);
        float x = cos(slr_freq * asin(z)) * sqrt(1.0 - z * z);
        float y = sin(slr_freq * asin(z)) * sqrt(1.0 - z * z);

        data.time(0)(proj, 0) = (double)proj;
        for (int i = 0; i < Xres; i++) {
          float r = (float)i * 0.5;
          data.kx(0)(i, proj, 0) = x * r;
          data.ky(0)(i, proj, 0) = y * r;
          data.kz(0)(i, proj, 0) = z * r;
          data.kw(0)(i, proj, 0) = r * r;
          data.kt(0)(i, proj, 0) = i * dt;
        }
      }

      // Phantom Gridding Operator
      std::cout << "Precalc gridding for pahntom" << std::endl
                << std::flush;
      gridFFT phantom_gridding;
      phantom_gridding.kernel_type = KAISER_KERNEL;
      phantom_gridding.overgrid = 1.5;
      phantom_gridding.dwinX = 6;
      phantom_gridding.dwinY = 6;
      phantom_gridding.dwinZ = 6;
      phantom_gridding.precalc_gridding(2 * data.xres, 2 * data.yres, 2 * data.zres, data);

      // Update Image
      std::cout << "Precalc gridding for phantom" << std::endl
                << std::flush;

      // Now Inverse Grid
      std::cout << " Create Phantom Image" << std::endl
                << std::flush;
      Array<complex<float>, 3> IMAGE(2 * data.xres, 2 * data.yres, 2 * data.zres, ColumnMajorArray<3>());

      for (int k = 0; k < IMAGE.length(thirdDim); k++) {
        float z = ((float)k / (float)IMAGE.length(thirdDim)) - 0.5;

        for (int j = 0; j < IMAGE.length(secondDim); j++) {
          float y = ((float)j / (float)IMAGE.length(secondDim)) - 0.5;

          for (int i = 0; i < IMAGE.length(firstDim); i++) {
            float x = ((float)i / (float)IMAGE.length(firstDim)) - 0.5;

            float radius = sqrt(x * x + y * y + z * z);
            if (radius > 0.5 * 0.9) {
              IMAGE(i, j, k) = 0.0;
            } else {
              int p = (int)(radius * 10 / 0.9);
              if ((p % 2) == 0) {
                IMAGE(i, j, k) = complex<float>(1.0, 0.0);
              } else
                IMAGE(i, j, k) = complex<float>(0.0, 0.0);
            }
          }
        }
      }

      std::cout << "Max Image = " << max(abs(IMAGE)) << std::endl;

      Array<float, 3> TempWeight(data.kx(0).shape(), ColumnMajorArray<3>());
      TempWeight = 1.0;
      data.kdata(0, 0) = 0.0;  // Set to zero
      phantom_gridding.backward(IMAGE, data.kdata(0, 0), data.kx(0), data.ky(0), data.kz(0), TempWeight);

    } break;

    case (RECON::PSF): {
      // Use external Kx,Ky,Kz
      data.read_external_data(recon.filename);
      for (Array<Array<complex<float>, 3>, 2>::iterator miter = data.kdata.begin(); miter != data.kdata.end(); miter++) {
        *miter = complex<float>(1.0, 0.0);
      }
    } break;

    case (RECON::PHANTOM): {
      // Use external Kx,Ky,Kz
      data.read_external_data(recon.filename);
      for (Array<Array<complex<float>, 3>, 2>::iterator miter = data.kdata.begin(); miter != data.kdata.end(); miter++) {
        *miter = complex<float>(0.0, 0.0);
      }

      // Initialize Phantom
      cout << "Phantom " << endl;
      PHANTOM phantom;
      phantom.read_commandline(argc, argv);
      phantom.init(data.xres, data.yres, data.zres, recon.rcframes);

      // More accurate gridding for Phantom
      cout << "Grid " << endl;
      gridFFT phantom_gridding;
      phantom_gridding.kernel_type = KAISER_KERNEL;
      phantom_gridding.overgrid = 1.5;
      phantom_gridding.dwinX = 6;
      phantom_gridding.dwinY = 6;
      phantom_gridding.dwinZ = 6;
      phantom_gridding.precalc_gridding(phantom.IMAGE.length(firstDim),
                                        phantom.IMAGE.length(secondDim),
                                        phantom.IMAGE.length(thirdDim), data);

      GATING gate(argc, argv);
      // Weighting Array for Time coding
      Array<float, 3> TimeWeight(data.kx(0).shape(), ColumnMajorArray<3>());
      gate.init(data, recon.rcframes);

      /*-----------------------------
         Collect Data
       ------------------------------*/
      int e = 0;
      Array<float, 3> kxE = data.kx(e);
      Array<float, 3> kyE = data.ky(e);
      Array<float, 3> kzE = data.kz(e);
      Array<float, 3> kwE = data.kw(e);

      cout << "Num coils = " << data.Num_Coils << endl;
      for (int coil = 0; coil < data.Num_Coils; coil++) {
        cout << "Getting phantom" << coil << ":" << flush;
        cout << "Smap " << endl
             << flush;
        phantom.update_smap_biotsavart(coil, data.Num_Coils);

        // Update Image
        Array<complex<float>, 3> kdataE = data.kdata(e, coil);  // Get one encoding, one coil
        kdataE = 0;                                             // Zero data

        Array<complex<float>, 3> temp(kdataE.shape(), ColumnMajorArray<3>());

        for (int t = 0; t < recon.rcframes; t++) {
          // Get Image
          phantom.calc_image(t, recon.rcframes);

          // Weight Image
          TimeWeight = 1;
          gate.weight_data(TimeWeight, e, kxE, kyE, kzE, t, GATING::NON_ITERATIVE, GATING::TIME_FRAME);

          // Now Inverse Grid
          cout << " Inverse Grid :: " << t << endl;
          phantom_gridding.backward(phantom.IMAGE, kdataE, kxE, kyE, kzE, TimeWeight);
        }
      }

      // Add Noise
      phantom.add_noise(data.kdata);
      data.write_external_data("PhantomData.h5");
      phantom.write_matlab_truth_script("PhantomData/");

      cout << "Only generating data" << endl;
      exit(1);

    } break;

    default:
    case (RECON::EXTERNAL): {
      // Read in External Data Format
      data.read_external_data(recon.filename);
      data.scale_fov(1. / recon.zoom_x, 1. / recon.zoom_y, 1. / recon.zoom_z);
      data.demod_kdata(recon.demod_freq);
    } break;
  }

  // --------------------------------------------------
  // Code for recon (no PSD specific data/structures)
  // --------------------------------------------------
  recon.init_recon(argc, argv, data);
  Array<Array<complex<float>, 3>, 2> X = recon.reconstruct_all_frames(data);

  // ------------------------------------
  // Post Processing + Export
  // ------------------------------------

  HDF5 clearRAW("Images.h5", "w");
  for (int ee = 0; ee < X.length(firstDim); ee++) {
    for (int tt = 0; tt < X.length(secondDim); tt++) {
      char fname[80];
      sprintf(fname, "X_%03d_%03d.dat", ee, tt);
      Array<float, 3> IMAGE;
      IMAGE.setStorage(ColumnMajorArray<3>());
      IMAGE.resize(X(0, 0).shape());
      IMAGE = abs(X(tt, ee));
      clearRAW.AddH5Array("IMAGES", fname, IMAGE);
    }
  }

  HDF5 complexRAW("ComplexImages.h5", "w");
  for (int ee = 0; ee < X.length(firstDim); ee++) {
    for (int tt = 0; tt < X.length(secondDim); tt++) {
      char fname[80];
      sprintf(fname, "X_%03d_%03d.dat", ee, tt);
      complexRAW.AddH5Array("IMAGES", fname, X(tt, ee));
    }
  }

  return (0);
}
