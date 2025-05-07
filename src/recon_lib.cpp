
#include "recon_lib.h"
#include "io_templates.hpp"

using arma::cx_mat;
using arma::uvec;
using arma::vec;
using namespace NDarray;

// --------------------------------
//   Call for help without existing classes
//
void print_mri_recon_help(void) {
  cout << "----------------------------------------------" << endl;
  cout << "   MRI recon options " << endl;
  cout << "    this uses multiple modules with their own options" << endl;
  cout << "----------------------------------------------" << endl;

  const char *Oargv[] = {"prog", "-h", NULL};
  int Oargc = (int)(sizeof(Oargv) / sizeof(Oargv[0])) - 1;
  RECON recon;
  recon.parse_commandline(Oargc, Oargv);

  gridFFT::help_message();
  SPIRIT::help_message();
  THRESHOLD::help_message();
  PHANTOM::help_message();
  GATING::help_message();
  L2REG::help_message();
  LOWRANKCOIL::help_message();
}

// ----------------------
//  Basic constructor (no args)
// ----------------------
RECON::RECON(void) { set_defaults(); }

// ----------------------
//  Sets Default Recon Parameterrs
// ----------------------
void RECON::set_defaults(void) {
  // Help Message for recon
  recon_type = SOS;
  data_type = EXTERNAL;
  coil_combine_type = LOWRES;

  complex_diff = false;

  cs_spatial_transform = WAVELET;
  cs_temporal_transform = NONE;
  cs_encode_transform = NONE;

  zoom = 1.0;
  zoom_x = 1.0;
  zoom_y = 1.0;
  zoom_z = 1.0;

  demod_freq = 0.0;

  rcxres = -1;
  rcyres = -1;
  rczres = -1;
  rcframes = -1;
  rcencodes = 1;
  rc_frame_start = 0;
  parallel_coils = false;

  smap_res = 8;
  intensity_correction = false;
  intensity_correction_blurX = 20.0;
  intensity_correction_blurY = 20.0;
  intensity_correction_blurZ = 20.0;
  reset_dens = false;

  threads = -1;
  acc = 1;
  compress_coils = 0.0;
  compress_kr = 32;
  whiten = false;
  noise_scale_factor = 1.0;  //  2 doubles noise, 3 triples noise, ...
  export_smaps = false;
  debug_smaps = false;
  max_iter = 50;
  smap_use_all_encodes = false;
  smap_nex_encodes = false;
  smap_thresh = 0.0;
  smap_mask = SMAPMASK_NONE;
  smap_norm = SMAPNORM_SOS;
  body_coil_idx = -1;

  coil_rejection_flag = false;
  coil_rejection_radius = 0.5;
  coil_rejection_height = 0.5;
  coil_rejection_shape = 1;
  coil_rejection_thresh = 0.7;

  cycle_spins = 4;

  walsh_block_sizeX = 8;
  walsh_block_sizeY = 8;
  walsh_block_sizeZ = 8;

  // Gaussian Blur of smaps
  extra_blurX = 0.0;
  extra_blurY = 0.0;
  extra_blurZ = 0.0;

  prep_done = false;

  wavelet_levelsX = 4;
  wavelet_levelsY = 4;
  wavelet_levelsZ = 4;
  wavelet_type = WAVELET3D::WAVE_DB4;

  // Code to rotate
  phase_rotation = false;
  phase_rotation_sX = 20;
  phase_rotation_sY = 20;
  phase_rotation_sZ = 20;

  pregate_data_flag = false;

  dcf_type = SUPPLIED;
  dcf_iter = 20;
  dcf_dwin = 2.1;
  dcf_scale = 1.0;
  dcf_overgrid = 2.1;
  dcf_acc = 1.0;

  admm_gamma = 0.1;
  admm_max_iter = 20;
  admm_rho = 0.5;

  cauchy_update_number = 30;
  max_eigen_iterations = 30;
  iterative_step_type = STEP_MAXEIG;
  fast_maxeig = false;
  image_scale_normalization = false;
}

// ----------------------
//  Constructor with Command Line Read
// ----------------------

RECON::RECON(int numarg, const char **pstring) {
  set_defaults();

  // --------------------------
  // Help Messages for Commandline Inputs
  //   -Please add your own help message for new classes
  // --------------------------
  for (int pos = 0; pos < numarg; pos++) {
    if ((strcmp(pstring[pos], "-h") == 0) ||
        (strcmp(pstring[pos], "-help") == 0) ||
        (strcmp(pstring[pos], "--help") == 0)) {
      print_mri_recon_help();
      exit(0);
    }
  }

  // Get Input Parameters
  parse_commandline(numarg, pstring);
}

inline string help_str(string para, string help_string) {
  // Padded string to 25 for format
  string full_help_string(40 - para.length(), ' ');
  string para_str(para);
  full_help_string.insert(4, para_str);
  full_help_string.append(":");
  full_help_string.append(help_string);
  return (full_help_string);
}

inline string help_strIO(string para, string help_string) {
  // Padded string to 25 for format
  para.append(" []");
  string full_help_string(40 - para.length(), ' ');
  string para_str(para);
  full_help_string.insert(4, para_str);
  full_help_string.append(":");
  full_help_string.append(help_string);
  return (full_help_string);
}

// Better method without macros?
#define trig_flag(num, name, val, message)                         \
  help_message << help_str(string(name), string(message)) << endl; \
  for (int ii = 0; ii < numarg; ii++) {                            \
    if (strcmp(name, pstring[ii]) == 0) {                          \
      val = num;                                                   \
    }                                                              \
  }

#define float_flag(name, val, message)                               \
  help_message << help_strIO(string(name), string(message)) << endl; \
  for (int ii = 0; ii < numarg; ii++) {                              \
    if (strcmp(name, pstring[ii]) == 0) {                            \
      ii++;                                                          \
      val = atof(pstring[ii]);                                       \
    }                                                                \
  }

#define string_flag(name, val, message)                              \
  help_message << help_strIO(string(name), string(message)) << endl; \
  for (int ii = 0; ii < numarg; ii++) {                              \
    if (strcmp(name, pstring[ii]) == 0) {                            \
      ii++;                                                          \
      val = string(pstring[ii]);                                     \
    }                                                                \
  }

#define int_flag(name, val, message)                                 \
  help_message << help_strIO(string(name), string(message)) << endl; \
  for (int ii = 0; ii < numarg; ii++) {                              \
    if (strcmp(name, pstring[ii]) == 0) {                            \
      ii++;                                                          \
      val = atoi(pstring[ii]);                                       \
    }                                                                \
  }

void RECON::parse_commandline(int numarg, const char **pstring) {
  stringstream help_message;

  help_message << "----------------------------------------------" << endl;
  help_message << "   Recon :: Command line options " << endl;
  help_message << "----------------------------------------------" << endl;
  string_flag("-f", filename, "File name of MRI_Raw to load");

  // Source of data
  trig_flag(EXTERNAL, "-external_data", data_type, "load from .h5 file");
  trig_flag(PHANTOM, "-phantom", data_type, "create and use a phanton");
  trig_flag(SIMULATE, "-simulate", data_type, "simulate coordinates and data");
  trig_flag(PSF, "-psf", data_type, "load from .h5 but set data to ones");

  help_message << "----------------------------------------------" << endl;
  help_message << "  Recon :: Resolution and FOV control " << endl;
  help_message << "----------------------------------------------" << endl;

  // Reconstruction Geometry
  int_flag("-recon_xres", rcxres, "reconstructed resolution x");
  int_flag("-recon_yres", rcyres, "reconstructed resolution y");
  int_flag("-recon_zres", rczres, "reconstructed resolution z");
  int_flag("-rcframes", rcframes, "reconstructed frames");
  int_flag("-rc_frame_start", rc_frame_start, "offset to start reconstructing images");

  float_flag("-zoom", zoom, "global zoom in x/y/z");
  float_flag("-zoom_x", zoom_x, "zoom in x");
  float_flag("-zoom_y", zoom_y, "zoom in y");
  float_flag("-zoom_z", zoom_z, "zoom in z");

  help_message << "----------------------------------------------" << endl;
  help_message << "  Recon :: Recon types and options" << endl;
  help_message << "----------------------------------------------" << endl;

  // Type of Recons
  trig_flag(SOS, "-sos", recon_type, "coil combination using sum of squares");
  trig_flag(CG, "-isense", recon_type, "iterative sense using CG");
  trig_flag(PILS, "-pils", recon_type, "coil combination using derived sensitivity maps");
  trig_flag(IST, "-ist", recon_type, "gradient descent with iterative soft thresholding");
  trig_flag(FISTA, "-fista", recon_type, "fast iterative gradient descent with thresholding");
  trig_flag(CLEAR, "-clear", recon_type, "self calibrating with low rank constraint");
  trig_flag(CG, "-cg", recon_type, "same as isense");
  trig_flag(ADMM, "-admm", recon_type, "admm is not working");

  // Spatial Transforms
  help_message << help_strIO(string("-spatial_transform"), string("none/wavelet"))
               << endl;
  for (int ii = 0; ii < numarg; ii++) {
    if (strcmp("-spatial_transform", pstring[ii]) == 0) {
      ii++;
      if (ii == numarg) {
        cout << "Please provide spatial transform type..none/wavelet" << endl;
        exit(1);
      } else if (strcmp("none", pstring[ii]) == 0) {
        cs_spatial_transform = NONE;
      } else if (strcmp("wavelet", pstring[ii]) == 0) {
        cs_spatial_transform = WAVELET;
      } else {
        cout << "Please provide spatial transform type..none/wavelet" << endl;
        exit(1);
      }
    }
  }

  // Temporal Transforms
  help_message << help_strIO(string("-temporal_transform"), string("none/dft/diff/pca"))
               << endl;
  for (int ii = 0; ii < numarg; ii++) {
    if (strcmp("-temporal_transform", pstring[ii]) == 0) {
      ii++;
      if (ii == numarg) {
        cout << "Please provide temporal transform type..none/dft/diff/pca/composite_diff" << endl;
        exit(1);
      } else if (strcmp("none", pstring[ii]) == 0) {
        cs_temporal_transform = NONE;
      } else if (strcmp("wavelet", pstring[ii]) == 0) {
        cs_temporal_transform = WAVELET;
      } else if (strcmp("diff", pstring[ii]) == 0) {
        cs_temporal_transform = DIFF;
      } else if (strcmp("dft", pstring[ii]) == 0) {
        cs_temporal_transform = DFT;
      } else if (strcmp("pca", pstring[ii]) == 0) {
        cs_temporal_transform = PCA;
      } else if (strcmp("composite_diff", pstring[ii]) == 0) {
        cs_temporal_transform = COMPOSITE_DIFF;
      } else {
        cout << "Please provide temporal transform type..none/dft/diff/pca/composite_diff" << endl;
        exit(1);
      }
    }
  }

  // Temporal Transforms
  help_message << help_strIO(string("-encode_transform"), string("none/dft/diff/pca/composite_diff"))
               << endl;
  for (int ii = 0; ii < numarg; ii++) {
    if (strcmp("-encode_transform", pstring[ii]) == 0) {
      ii++;
      if (ii == numarg) {
        cout << "Please provide encode transform type..none/dft/diff/pca/composite_diff" << endl;
        exit(1);
      } else if (strcmp("none", pstring[ii]) == 0) {
        cs_temporal_transform = NONE;
      } else if (strcmp("wavelet", pstring[ii]) == 0) {
        cs_temporal_transform = WAVELET;
      } else if (strcmp("dft", pstring[ii]) == 0) {
        cs_temporal_transform = DFT;
      } else if (strcmp("pca", pstring[ii]) == 0) {
        cs_temporal_transform = PCA;
      } else if (strcmp("composite_diff", pstring[ii]) == 0) {
        cs_temporal_transform = COMPOSITE_DIFF;
      } else {
        cout << "Please provide encode transform type..none/dft/diff/pca/composite_diff" << endl;
        exit(1);
      }
    }
  }

  int_flag("-cauchy_update_number", cauchy_update_number, "recalc step size this many times before fixing");
  trig_flag(true, "-image_scale_normalization", image_scale_normalization, "scale images to be roughly 1");

  // Temporal Transforms
  help_message << help_strIO(string("-iterative_step_type"), string("cauchy/maxeig[default]"))
               << endl;
  for (int ii = 0; ii < numarg; ii++) {
    if (strcmp("-iterative_step_type", pstring[ii]) == 0) {
      ii++;
      if (ii == numarg) {
        cout << "Please provide iterative_step_type type for gradient descent cauchy/maxeig[default]" << endl;
        exit(1);
      } else if (strcmp("maxeig", pstring[ii]) == 0) {
        iterative_step_type = STEP_MAXEIG;
      } else if (strcmp("cauchy", pstring[ii]) == 0) {
        iterative_step_type = STEP_CAUCHY;
      } else {
        cout << "Please provide iterative_step_type type for gradient descent cauchy/maxeig[default]" << endl;
        exit(1);
      }
    }
  }

  int_flag("-wavelet_levelsX", wavelet_levelsX, "max number of wavelet levels x");
  int_flag("-wavelet_levelsY", wavelet_levelsY, "max number of wavelet levels x");
  int_flag("-wavelet_levelsZ", wavelet_levelsZ, "max number of wavelet levels x");

  // Temporal Transforms
  help_message
      << help_strIO(string("-wavelet_type"), string("db2/db4/db8/sym2/sym4/sym8/bo33"))
      << endl;
  for (int ii = 0; ii < numarg; ii++) {
    if (strcmp("-wavelet_type", pstring[ii]) == 0) {
      ii++;
      if (ii == numarg) {
        cout << "Please provide wavelet transform type [db2/db4/db8/sym2/sym4/sym8/bo33] " << endl;
        exit(1);
      } else if (strcmp("db2", pstring[ii]) == 0) {
        wavelet_type = WAVELET3D::WAVE_DB2;
      } else if (strcmp("db4", pstring[ii]) == 0) {
        wavelet_type = WAVELET3D::WAVE_DB4;
      } else if (strcmp("db8", pstring[ii]) == 0) {
        wavelet_type = WAVELET3D::WAVE_DB8;
      } else if (strcmp("sym2", pstring[ii]) == 0) {
        wavelet_type = WAVELET3D::WAVE_SYM2;
      } else if (strcmp("sym4", pstring[ii]) == 0) {
        wavelet_type = WAVELET3D::WAVE_SYM4;
      } else if (strcmp("sym8", pstring[ii]) == 0) {
        wavelet_type = WAVELET3D::WAVE_SYM8;
      } else if (strcmp("bo33", pstring[ii]) == 0) {
        wavelet_type = WAVELET3D::WAVE_BO33;
      } else {
        cout << "Please provide iterative_step_type type for gradient descent cauchy/maxeig[default]" << endl;
        exit(1);
      }
    }
  }

  // Iterations for IST
  int_flag("-max_iter", max_iter, "maximum number of iterations");
  int_flag("-cycle_spins", cycle_spins, "number of shifts to apply to images each iteration");

  float_flag("-admm_alpha", admm_alpha, "");
  float_flag("-admm_gamma", admm_gamma, "");
  int_flag("-admm_max_iter", admm_max_iter, "");
  float_flag("-admm_rho", admm_rho, "");

  help_message << "----------------------------------------------" << endl;
  help_message << "  Recon :: Data modifications" << endl;
  help_message << "----------------------------------------------" << endl;

  trig_flag(RECALC_DCF, "-recalc_dcf", dcf_type, "recalculate dcf using Pipe gridding method");
  trig_flag(RECALC_VOR, "-recalc_vor", dcf_type, "recalculate usin Voronoi diagram");
  int_flag("-dcf_iter", dcf_iter, "iterations for Pipe method");
  float_flag("-dcf_dwin", dcf_dwin, "window for Pipe method");
  float_flag("-dcf_scale", dcf_scale, "resolution scale to perfrom Pipe dcf compensation");
  float_flag("-dcf_overgrid", dcf_overgrid, "overgridding factor for dcf");
  float_flag("-dcf_acc", dcf_acc, "acceleration factor to assume for Pipe dcf");
  trig_flag(true, "-reset_dens", reset_dens, "set dcf to ones");
  float_flag("-demod", demod_freq, "frequency in Hz to demodulation the data");

  // Data modification
  int_flag("-acc", acc, "remove data with this acceleration integer");
  int_flag("-compress_coils", compress_coils, "compress to this many coils using PCA");
  float_flag("-compress_kr", compress_kr, "kspace radius of data to use for coil compression");
  trig_flag(true, "-whiten", whiten, "whiten the data using noise samples in MRI raw");
  float_flag("-noise_scale_factor", noise_scale_factor, "additive noise factor for simulating noisy data");
  trig_flag(true, "-complex_diff", complex_diff, "perfrom subtraction of first frame");

  help_message << "----------------------------------------------" << endl;
  help_message << "  Recon :: Coil sensitivity mapping" << endl;
  help_message << "----------------------------------------------" << endl;

  trig_flag(ESPIRIT, "-espirit", coil_combine_type, "Use ESPiRIT for coil sensitity estimate");
  trig_flag(WALSH, "-walsh", coil_combine_type, "Use blocked SVT (Walsh et al) for coil sensitity estimate");
  trig_flag(LOWRES, "-coil_lowres", coil_combine_type, "Use low resolution images for coil sensitity estimate");
  float_flag("-smap_res", smap_res, "resolution of images for coil sensitity mapping");
  trig_flag(true, "-export_smaps", export_smaps, "save sensitivity maps to .h5");
  trig_flag(true, "-debug_smaps", debug_smaps, "save debugging sensitity map images to .h5");
  trig_flag(true, "-intensity_correction", intensity_correction, "attempt to intensity correct the images");
  float_flag("-intensity_correction_blurX", intensity_correction_blurX, "Gaussian blurring in x for intensity correction");
  float_flag("-intensity_correction_blurY", intensity_correction_blurY, "Gaussian blurring in x for intensity correction");
  float_flag("-intensity_correction_blurZ", intensity_correction_blurZ, "Gaussian blurring in x for intensity correction");
  int_flag("-walsh_block_sizeX", walsh_block_sizeX, "size of block in x for Walsh method");
  int_flag("-walsh_block_sizeY", walsh_block_sizeY, "size of block in y for Walsh method");
  int_flag("-walsh_block_sizeZ", walsh_block_sizeZ, "size of block in z for Walsh method");
  float_flag("-extra_blurX", blurX, "Gaussian blurring in x for coil sensitivity mapping");
  float_flag("-extra_blurY", blurY, "Gaussian blurring in y for coil sensitivity mapping");
  float_flag("-extra_blurZ", blurZ, "Gaussian blurring in z for coil sensitivity mapping");
  trig_flag(true, "-phase_rotation", phase_rotation, "rotate the phase of the coil maps to match first coil");
  trig_flag(SMAPNORM_SOS, "-smap_norm_sos", smap_norm, "normalize the sensitivity maps to sum of squares");
  trig_flag(SMAPNORM_BODY, "-smap_norm_body", smap_norm, "normalize the sensitivity maps using body coil");
  int_flag("-body_coil_idx", body_coil_idx, "index of body coil for normalization");

  trig_flag(true, "-smap_use_all_encodes", smap_use_all_encodes, "use all encodes for mapping instead of just the first");
  trig_flag(true, "-smap_nex_encodes", smap_nex_encodes, "use all encodes for mapping but just add them");
  float_flag("-smap_thresh", smap_thresh, "fraction of max signal to threshold the images");
  trig_flag(true, "-coil_rejection", coil_rejection_flag, "turn on coil rejection");
  float_flag("-coil_rejection_radius", coil_rejection_radius, "fractional radius to reject coils outside of FOV");
  float_flag("-coil_rejection_height", coil_rejection_height, "height of cylinder to reject outside of FOV");
  float_flag("-coil_rejection_thresh", coil_rejection_thresh, "threshold to reject coils outside of FOV");
  int_flag("-coil_rejection_shape", coil_rejection_shape, "shape for coil rejection (0=sphere, 1=cylinder");

  help_message << help_strIO(string("-smap_mask"), string("none/circle/sphere"))
               << endl;
  for (int ii = 0; ii < numarg; ii++) {
    if (strcmp("-smap_mask", pstring[ii]) == 0) {
      ii++;
      if (ii == numarg) {
        cout << "Please provide sensitivty map..none/circle/sphere" << endl;
        exit(1);
      } else if (strcmp("none", pstring[ii]) == 0) {
        smap_mask = SMAPMASK_NONE;
      } else if (strcmp("circle", pstring[ii]) == 0) {
        smap_mask = SMAPMASK_CIRCLE;
      } else if (strcmp("sphere", pstring[ii]) == 0) {
        smap_mask = SMAPMASK_SPHERE;
      } else {
        cout << "Please provide sensitivty map..none/circle/sphere" << endl;
        exit(1);
      }
    }
  }

  help_message << help_strIO(string("-smap_skip_encode"), string("skip this encode for sensitivity mapping"))
               << endl;
  for (int ii = 0; ii < numarg; ii++) {
    if (strcmp("-smap_skip_encode", pstring[ii]) == 0) {
      ii++;
      smap_skip_encode.resizeAndPreserve(smap_skip_encode.numElements() + 1);
      smap_skip_encode(smap_skip_encode.numElements() - 1) = atoi(pstring[ii]);
    }
  }

  help_message << "----------------------------------------------" << endl;
  help_message << "  Recon :: Fine/compute control " << endl;
  help_message << "----------------------------------------------" << endl;

  trig_flag(true, "-parallel_coils", parallel_coils, "parallize over coils, uses more memory");
  int_flag("-threads", threads, "set number of threads");
  trig_flag(true, "-pregate_data", pregate_data_flag, "gate the data prior to reconstructing");
  int_flag("-max_eigen_iterations", max_eigen_iterations, "power iterations for max eigen calculations");
  trig_flag(true, "-fast_maxeig", fast_maxeig, "Calc step size on single frame");

  // First scan for help message
  for (int pos = 0; pos < numarg; pos++) {
    if (strcmp("-h", pstring[pos]) == 0) {
      cout << help_message.str();
    }
  }
}

void RECON::init_recon(int argc, const char **argv, MRI_DATA &data) {
  // Use resultion in the header
  rcxres = (rcxres == -1) ? (data.recon_res(0)) : (rcxres);
  rcyres = (rcyres == -1) ? (data.recon_res(1)) : (rcyres);
  rczres = (rczres == -1) ? (data.recon_res(2)) : (rczres);
  rcencodes = data.Num_Encodings;

  // Whiten
  if (whiten) {
    data.whiten();  // Requires noise samples inserted
    if (noise_scale_factor > 1.0) {
      data.add_noise(noise_scale_factor);
    }
  }

  // Option to compress coils
  if (compress_coils > 0) {
    data.coilcompress(compress_coils, compress_kr);
  }

  // Matlab like timer (openmp code base)
  tictoc T;

  // Setup Gridding + FFT Structure
  gridding.read_commandline(argc, argv);
  gridding.precalc_gridding(rczres, rcyres, rcxres, data);

  // Option to parallelize overs coils
  if (this->parallel_coils) {
    gridding_CoilThreaded.read_commandline(argc, argv);
    gridding_CoilThreaded.precalc_gridding(rczres, rcyres, rcxres, data);
  }

  // Recalculate the Density compensation
  switch (dcf_type) {
    default: {
      // Use supplied
    } break;

    case (RECALC_DCF): {
      dcf_calc(data);
    } break;

    case (RECALC_VOR): {
#ifdef USE_VORO
      if (data.trajectory_type(2) == MRI_DATA::NONCARTESIAN) {
        VORONOI_DCF::vor_dcf(data.kw(0), data.kx(0), data.ky(0), data.kz(0),
                             VORONOI_DCF::SPHERE);
      } else {
        VORONOI_DCF::vor_dcf(data.kw(0), data.kx(0), data.ky(0), data.kz(0),
                             VORONOI_DCF::CYLINDER);
      }
#else
      throw std::invalid_argument(
          "Invalid DCF Type [voronoi], code was not compiled with VORO++. "
          "Recompile with VORO++ to use.");
#endif
    } break;
  }

  // Calculate Sensitivity maps (using gridding struct)
  calc_sensitivity_maps(argc, argv, data);

  // Complex Difference
  if (complex_diff) {
    cout << "Doing Complex Diff" << endl;
    for (int coil = 0; coil < data.Num_Coils; coil++) {
      // Subtract off reference
      for (int e = 1; e < rcencodes; e++) {
        Array<complex<float>, 3> kdata1 = data.kdata(0, coil);
        Array<complex<float>, 3> kdata2 = data.kdata(e, coil);
        kdata2 -= kdata1;
      }

      // Rearrange Positions
      for (int e = 1; e < rcencodes; e++) {
        Array<complex<float>, 3> kdata1 = data.kdata(e - 1, coil);
        Array<complex<float>, 3> kdata2 = data.kdata(e, coil);
        kdata1 = kdata2;
      }
    }
    rcencodes -= 1;
    data.Num_Encodings -= 1;

    /* Need to resize Gating
    data.time.resizeAndPreserve(data.Num_Encodings);
    data.resp.resizeAndPreserve(data.Num_Encodings);
    data.prep.resizeAndPreserve(data.Num_Encodings);
    data.ecg.resizeAndPreserve(data.Num_Encodings);
    */
  }

  // -------------------------------------
  //	This handles all the gating, assuming mri_data physio data is populated
  // -------------------------------------
  gate = GATING(argc, argv);
  gate.init(data, rcframes);
  rcframes = gate.number_of_frames();
  std::cout << "Reconstructing to " << rcframes << std::endl;

  if (pregate_data_flag) {
    pregate_data(data);
  }

  // -------------------------------------
  //	This handles special preperations
  // -------------------------------------
  switch (recon_type) {
    // Non-Iterative Recons
    default:
    case (SOS):
    case (PILS): {
    } break;

    // Iterative Recons allow Regularization
    case (CG):
    case (IST):
    case (FISTA):
    case (CLEAR):
    case (ADMM): {
      // Setup 3D Wavelet
      wave = WAVELET3D(
          TinyVector<int, 3>(rcxres, rcyres, rczres),
          TinyVector<int, 3>(wavelet_levelsX, wavelet_levelsY, wavelet_levelsZ),
          wavelet_type);

      // Setup Soft Thresholding
      softthresh = THRESHOLD(argc, argv);

      // CLEAR
      if (recon_type == CLEAR) {
        lrankcoil = LOWRANKCOIL(argc, argv);
      }

      // Setup L2 Regularization Thresholding
      l2reg = L2REG(argc, argv);

      if (cs_temporal_transform == COMPOSITE_DIFF) {
        composite_image.free();
        composite_image.setStorage(ColumnMajorArray<3>());
        composite_image.resize(rcxres, rcyres, rczres);
        composite_image = complex<float>(0.0, 0.0);
      }
    } break;
  }

  // TEMPX
  lranktime = LOWRANKCOIL(argc, argv);

  // Signal that the recon is ready to perform the requested reconstruction
  prep_done = true;
}

/**
 * This function parses an MRI_DATA structure into time frames such that this is
 * not required during recon
 * @see setup()
 * @param data MRI_DATA to be transfromed
 */
void RECON::pregate_data(MRI_DATA &data) {
  // For now balloon memory
  MRI_DATA data2;

  data2.Num_Encodings = data.Num_Encodings * rcframes;
  data2.Num_Coils = data.Num_Coils;
  data2.Num_Frames = rcframes;

  // Resize the data structures but don't allocate sub-structure
  data2.kx.resize(data2.Num_Encodings);
  data2.ky.resize(data2.Num_Encodings);
  data2.kz.resize(data2.Num_Encodings);
  data2.kw.resize(data2.Num_Encodings);
  data2.kdata.setStorage(ColumnMajorArray<2>());
  data2.kdata.resize(data2.Num_Encodings, data2.Num_Coils);

  int count = 0;
  for (int e = 0; e < data.Num_Encodings; e++) {
#pragma omp parallel for
    for (int t = 0; t < rcframes; t++) {
      Array<float, 3> Kweight;
      Kweight.setStorage(ColumnMajorArray<3>());
      Kweight.resize(data.kx(e).shape());

      // First count th number of points required
      Kweight = 1;
      gate.weight_data(Kweight, e, data.kx(e), data.ky(e), data.kz(e), t, GATING::NON_ITERATIVE, GATING::TIME_FRAME);

      // Count the number of frames
      int number_of_points = 0;
      for (Array<float, 3>::iterator miter = Kweight.begin(); miter != Kweight.end(); miter++) {
        if (*miter > 0) {
          number_of_points++;
        }
      }
      float Kw_Scale = (1e6 / ((float)number_of_points));
#pragma omp critical
      {
        cout << "Frame " << t
             << ", Total number of point = " << number_of_points
             << " ,Kw_scale = " << Kw_Scale << endl;
      }

      data2.kx(count + t).setStorage(ColumnMajorArray<3>());
      data2.kx(count + t).resize(number_of_points, 1, 1);

      data2.ky(count + t).setStorage(ColumnMajorArray<3>());
      data2.ky(count + t).resize(number_of_points, 1, 1);

      data2.kz(count + t).setStorage(ColumnMajorArray<3>());
      data2.kz(count + t).resize(number_of_points, 1, 1);

      data2.kw(count + t).setStorage(ColumnMajorArray<3>());
      data2.kw(count + t).resize(number_of_points, 1, 1);

      for (int coil = 0; coil < data2.Num_Coils; coil++) {
        data2.kdata(count + t, coil).setStorage(ColumnMajorArray<3>());
        data2.kdata(count + t, coil).resize(number_of_points, 1, 1);
      }

      int point_number = 0;
      for (int k = 0; k < Kweight.length(thirdDim); k++) {
        for (int j = 0; j < Kweight.length(secondDim); j++) {
          for (int i = 0; i < Kweight.length(firstDim); i++) {
            if (Kweight(i, j, k) > 0) {
              data2.kx(count + t)(point_number, 0, 0) = data.kx(e)(i, j, k);
              data2.ky(count + t)(point_number, 0, 0) = data.ky(e)(i, j, k);
              data2.kz(count + t)(point_number, 0, 0) = data.kz(e)(i, j, k);
              if (reset_dens) {
                data2.kw(count + t)(point_number, 0, 0) = 1.0;
              } else {
                data2.kw(count + t)(point_number, 0, 0) = Kw_Scale * data.kw(e)(i, j, k);
              }

              for (int coil = 0; coil < data2.Num_Coils; coil++) {
                data2.kdata(count + t, coil)(point_number, 0, 0) = data.kdata(e, coil)(i, j, k);
              }
              point_number++;
            }
          }
        }
      }
    }
    count += rcframes;
  }

  cout << "Swap" << endl
       << flush;
  cycleArrays(data.kx, data2.kx);
  cycleArrays(data.ky, data2.ky);
  cycleArrays(data.kz, data2.kz);
  cycleArrays(data.kw, data2.kw);
  cycleArrays(data.kdata, data2.kdata);

  data.Num_Encodings = data2.Num_Encodings;
  data.Num_Frames = data2.Num_Frames;

  // data.stats();

  cout << "Done gating data" << endl
       << flush;
}

int RECON::get_rcframes(void) {
  return (this->rcframes);
}

void vor_dcf2(Array<Array<float, 3>, 2> &Kw, Array<Array<float, 3>, 2> &Ky,
              Array<Array<float, 3>, 2> &Kz) {
  // Get total number of points
  int total_points =
      Kz.numElements() * Kz(0, 0).length(secondDim) * Kz(0, 0).length(thirdDim);

  cout << "Redoing with new DCF" << endl;
  Array<float, 3> ky(total_points, 1, 1, ColumnMajorArray<3>());
  Array<float, 3> kz(total_points, 1, 1, ColumnMajorArray<3>());
  Array<float, 3> kx(total_points, 1, 1, ColumnMajorArray<3>());
  kx = 1.0;
  Array<float, 3> kw(total_points, 1, 1, ColumnMajorArray<3>());

  cout << "Copy points " << endl;
  int count = 0;
  for (int e = 0; e < Ky.length(firstDim); e++) {
    for (int t = 0; t < Ky.length(secondDim); t++) {
      for (int jj = 0; jj < Ky(e, t).length(secondDim); jj++) {
        for (int kk = 0; kk < Ky(e, t).length(thirdDim); kk++) {
          ky(count, 0, 0) = Ky(e, t)(0, jj, kk);
          kz(count, 0, 0) = Kz(e, t)(0, jj, kk);
          count++;
        }
      }
    }
  }

#ifdef USE_VORO
  cout << "Compute Vor" << endl;
  VORONOI_DCF::vor_dcf(kw, ky, kz, kx, VORONOI_DCF::CUBE);
#else
  throw std::invalid_argument(
      "Invalid DCF Type [voronoi], code was not compiled with VORO++. "
      "Recompile with VORO++ to use.");
#endif

  ArrayWrite(kz, "NEW_Kz.dat");
  ArrayWrite(ky, "NEW_Ky.dat");
  ArrayWrite(kw, "NEW_Kw.dat");

  count = 0;
  for (int e = 0; e < Ky.length(firstDim); e++) {
    for (int t = 0; t < Ky.length(secondDim); t++) {
      for (int jj = 0; jj < Ky(e, t).length(secondDim); jj++) {
        for (int kk = 0; kk < Ky(e, t).length(thirdDim); kk++) {
          for (int ii = 0; ii < Ky(e, t).length(firstDim); ii++) {
            Kw(e, t)(ii, jj, kk) = kw(count, 0, 0);
          }
          count++;
        }
      }
    }
  }
}

Array<Array<complex<float>, 3>, 1> RECON::reconstruct_one_frame(MRI_DATA &data, int frame_number) {
  Array<Array<complex<float>, 3>, 2> XX = full_recon(data, Range(frame_number, frame_number), Range(0, 0), false);
  Array<Array<complex<float>, 3>, 1> X = XX(0, Range::all());
  return (X);
}

Array<Array<complex<float>, 3>, 2> RECON::reconstruct_all_frames(MRI_DATA &data) {
  Array<Array<complex<float>, 3>, 2> XX = full_recon(data, Range(rc_frame_start, rcframes - 1), Range(0, rcframes - 1 - rc_frame_start), false);
  return (XX);
}

Array<Array<complex<float>, 3>, 1> RECON::reconstruct_composite(MRI_DATA &data) {
  Array<Array<complex<float>, 3>, 2> XX = full_recon(data, Range(0, 0), Range(0, 0), true);
  Array<Array<complex<float>, 3>, 1> X = XX(0, Range::all());
  return (X);
}

void RECON::dcf_calc(MRI_DATA &data) {
  // Setup Gridding + FFT Structure
  DCFgridFFT dcf_gridding;
  dcf_gridding.kernel_type = DCFgridFFT::POLY_KERNEL;
  dcf_gridding.dwinX = dcf_dwin;
  dcf_gridding.dwinY = dcf_dwin;
  dcf_gridding.dwinZ = dcf_dwin;
  dcf_gridding.grid_x = dcf_overgrid;
  dcf_gridding.grid_y = dcf_overgrid;
  dcf_gridding.grid_z = dcf_overgrid;
  dcf_gridding.grid_in_x = 1;
  dcf_gridding.grid_in_y = 1;
  dcf_gridding.grid_in_z = 1;
  dcf_gridding.precalc_kernel();
  dcf_gridding.acc = dcf_acc;

  // dcf_gridding.precalc_gridding(rczres+16,rcyres+16,rcxres+16,data.trajectory_dims,data.trajectory_type);

  // Weighting Array for Time coding
  Array<float, 3> X(dcf_overgrid * rczres + 16,
                    dcf_overgrid * rcyres + 16,
                    dcf_overgrid * rcxres + 16, ColumnMajorArray<3>());

  for (int e = 0; e < rcencodes; e++) {
    Array<float, 3> Kweight(data.kx(e).shape(), ColumnMajorArray<3>());
    Array<float, 3> Kweight2(data.kx(e).shape(), ColumnMajorArray<3>());

    data.kx(e) *= dcf_scale;
    data.ky(e) *= dcf_scale;
    data.kz(e) *= dcf_scale;

    Kweight = 1;
    for (int iter = 0; iter < dcf_iter; iter++) {
      cout << "Iteration = " << iter << flush;

      X = 0;
      Kweight2 = 0.0;
      dcf_gridding.forward(X, Kweight, data.kx(e), data.ky(e), data.kz(e));    // Grid data to K-Space
      dcf_gridding.backward(X, Kweight2, data.kx(e), data.ky(e), data.kz(e));  // Degrid

      float kw_sum = 0.0;
      Array<float, 3>::iterator kw_iter = Kweight.begin();
      Array<float, 3>::iterator kw2_iter = Kweight2.begin();
      for (; (kw_iter != Kweight.end()) && (kw2_iter != Kweight2.end()); kw_iter++, kw2_iter++) {
        if (abs(*kw2_iter) < (1e-3)) {
          *kw_iter = 0.0;
        } else {
          *kw_iter = *kw_iter / *kw2_iter;
          kw_sum += abs(*kw_iter);
        }
      }

      cout << " Sum = " << kw_sum << endl;
    }

    dcf_gridding.scale_kw(Kweight, data.kx(e), data.ky(e), data.kz(e));
    data.kw(e) = abs(Kweight);

    data.kx(e) /= dcf_scale;
    data.ky(e) /= dcf_scale;
    data.kz(e) /= dcf_scale;
  }
  ArrayWrite(data.kw(0), "Kweight_DCF.dat");
}

Array<Array<complex<float>, 3>, 2> RECON::full_recon(MRI_DATA &data,
                                                     Range times,
                                                     Range times_store,
                                                     bool composite) {
  cout << "Starting Recon" << endl
       << flush;

  // Shorthand for Blitz++
  Range all = Range::all();

  // Matlab like timer (openmp code base)
  tictoc T;

  // Calculate shared structures
  if (!prep_done) {
    cout << "Error need to run prep before running recon (obj.init_recon(arg,argc, data). Exiting." << endl;
    exit(1);
  }

  /*----------------------------Main Recons---------------------------------------*/

  int Nt = times.length();
  GATING::FrameType frame_type;
  if (composite) {
    frame_type = GATING::COMPOSITE;
  } else {
    frame_type = GATING::TIME_FRAME;
  }

  // Final Image Solution
  cout << "Alloc Container for Solution "
       << "( " << rcxres << "," << rcyres << "," << rczres << ") x (" << Nt << "," << rcencodes << ")" << endl
       << flush;
  Array<Array<complex<float>, 3>, 2> X = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes);

  // Weighting Array for Time coding
  Array<float, 3> TimeWeight;
  TimeWeight.setStorage(ColumnMajorArray<3>());

  cout << "Full recon for " << rcencodes << " encodes, " << Nt << "frames "
       << endl
       << flush;

  switch (recon_type) {
    default:
    case (SOS):
    case (PILS): {
      for (int e = 0; e < rcencodes; e++) {
        for (int t = 0; t < Nt; t++) {
          int act_t = times(t);
          int store_t = times_store(t);
          int act_e = (pregate_data_flag) ? (e * Nt + act_t) : (e);

          cout << "Recon Encode" << e << " Frame " << t << endl;
          T.tic();

          // Temporal weighting
          if (pregate_data_flag) {
            TimeWeight.reference(data.kw(act_e));
          } else {
            TimeWeight.resize(data.kw(act_e).shape());
            TimeWeight = data.kw(act_e);
            gate.weight_data(TimeWeight, act_e, data.kx(act_e), data.ky(act_e), data.kz(act_e), act_t, GATING::NON_ITERATIVE, frame_type);
          }

          if (this->parallel_coils) {
            // K-space to Image
            Array<Array<complex<float>, 3>, 1> TEMP = data.kdata(act_e, Range::all());
            if (recon_type == PILS) {
              gridding_CoilThreaded.forward(X(store_t, e), smaps, TEMP,
                                            data.kx(act_e), data.ky(act_e), data.kz(act_e),
                                            TimeWeight);
            } else {
              gridding_CoilThreaded.forward_sos(X(store_t, e), TEMP,
                                                data.kx(act_e), data.ky(act_e), data.kz(act_e),
                                                TimeWeight);
            }
          } else {
            // cout << "\tForward Gridding Coil ";
            for (int coil = 0; coil < data.Num_Coils; coil++) {
              // K-space to Image
              if (recon_type == PILS) {
                gridding.forward(X(store_t, e), smaps(coil), data.kdata(act_e, coil),
                                 data.kx(act_e), data.ky(act_e), data.kz(act_e),
                                 TimeWeight);
              } else {
                gridding.forward_sos(X(store_t, e), data.kdata(act_e, coil),
                                     data.kx(act_e), data.ky(act_e), data.kz(act_e),
                                     TimeWeight);
              }
            }
          }
          cout << "Took " << T << endl;

          // Take Square Root for SOS
          if (recon_type == SOS) {
            X(store_t, e) = csqrt(X(store_t, e));
          }

          // Correct the intensity image
          correct_intensity_bias(X(store_t, e));
        }
      }

    } break;

    case (CLEAR): {
      cout << "----------------------------------------" << endl;
      cout << "\tStarting Clear based recon" << endl;
      cout << "----------------------------------------" << endl;
      typedef Array<complex<float>, 3> Complex3D;

      // Image and Residue (note now need to store coils as well)
      Array<Complex3D, 3> XX = Alloc6DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes, data.Num_Coils);
      Array<Complex3D, 3> RR = Alloc6DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes, data.Num_Coils);

      // Temp variable for E'ER
      Complex3D P(rcxres, rcyres, rczres, ColumnMajorArray<3>());

      cout << "Iterate" << endl;
      double error0 = 0.0;
      for (int iteration = 0; iteration < max_iter; iteration++) {
        tictoc iteration_timer;
        iteration_timer.tic();
        cout << "\nIteration = " << iteration << endl;

        //----------------------------------------------------
        //  First get Gradient Descent
        //----------------------------------------------------

        // Zero this for Cauchy set size
        complex<double> scale_RhP(0, 0);

        cout << "\tGradient Calculation" << endl;
        for (int e = 0; e < rcencodes; e++) {
          // Storage for (Ex-d)
          Complex3D diff_data(data.kx(e).shape(), ColumnMajorArray<3>());

          for (int t = 0; t < Nt; t++) {
            int act_t = times(t);
            int store_t = times_store(t);

            T.tic();

            // Get Sub-Arrays for Encoding
            Array<float, 3> kxE = data.kx(e);
            Array<float, 3> kyE = data.ky(e);
            Array<float, 3> kzE = data.kz(e);
            Array<float, 3> kwE = data.kw(e);

            // Temporal weighting

            TimeWeight.resize(kwE.shape());
            TimeWeight = kwE;
            gate.weight_data(TimeWeight, e, kxE, kyE, kzE, act_t, GATING::ITERATIVE, frame_type);

            for (int coil = 0; coil < data.Num_Coils; coil++) {
              T.tic();

              // Alloc Data
              Array<complex<float>, 3> diff_data(data.kx(e).shape(), ColumnMajorArray<3>());
              diff_data = 0;

              // Ex
              gridding.backward(XX(store_t, e, coil), diff_data, kxE, kyE, kzE, TimeWeight);

              // Ex-d
              Array<complex<float>, 3> kdataC = data.kdata(e, coil);
              diff_data -= kdataC;

              // E'(Ex-d)
              RR(t, e, coil) = complex<float>(0.0, 0.0);
              gridding.forward(RR(store_t, e, coil), diff_data, kxE, kyE, kzE, TimeWeight);

              // Now Get Scale
              P = 0;

              // EE'(Ex-d)
              diff_data = 0;
              gridding.backward(RR(store_t, e, coil), diff_data, kxE, kyE, kzE, TimeWeight);

              // E'EE'(Ex-d)
              gridding.forward(P, diff_data, kxE, kyE, kzE, TimeWeight);

              scale_RhP += conj_sum(P, RR(t, e, coil));

              cout << "Coil " << coil << " took " << T << endl;
            }  // Coils

            // cout << e << "," << t << "took " << T << "s" << endl;
          }  // Time
        }  // Encode

        // Get Scaling Factor R'P / R'R
        complex<double> scale_RhR = complex<double>(ArrayEnergy(RR), 0);

        // Error check
        if (iteration == 0) {
          error0 = abs(scale_RhR);
        }
        cout << "Residue Energy =" << scale_RhR << " ( " << (abs(scale_RhR) / error0) << " % )  " << endl;

        // Export R (across coils)
        Array<complex<float>, 2> Rslice = RR(0, 0, 0)(all, all, RR(0, 0, 0).length(2) / 2);
        ArrayWriteMag(Rslice, "R.dat");

        // Step in direction
        complex<double> scale_double = (scale_RhR / scale_RhP);
        complex<float> scale(real(scale_double), imag(scale_double));
        cout << "Scale = " << scale << endl;

        for (int coil = 0; coil < data.Num_Coils; coil++) {
          for (int e = 0; e < rcencodes; e++) {
            for (int t = 0; t < Nt; t++) {
              RR(t, e, coil) *= scale;
              XX(t, e, coil) -= RR(t, e, coil);
            }
          }
        }
        cout << "Took " << iteration_timer << " s " << endl;

        // Export X slice
        {
          Array<float, 2> Xslice(rcxres, rcyres, ColumnMajorArray<2>());

          Xslice = 0;
          for (int coil = 0; coil < data.Num_Coils; coil++) {
            Array<complex<float>, 2> Xtemp = XX(0, 0, coil)(all, all, RR(0, 0, 0).length(2) / 2);
            Xslice += norm(Xtemp);
            ArrayWriteMagAppend(Xtemp, "X_coils.dat");
          }
          Xslice = sqrt(Xslice);
          ArrayWriteAppend(Xslice, "X_mag.dat");
        }

        // ----------------------------------
        //   Soft Thresh
        // ----------------------------------
        {
          cout << "Soft thresh" << endl;
          Array<Complex3D, 2> X2 = XX(0, all, all);
          if (softthresh.getThresholdMethod() != TH_NONE) {
            L1_threshold(X2);
          }

          // Export X slice
          {
            Array<float, 2> Xslice;
            Xslice.setStorage(ColumnMajorArray<2>());
            Xslice.resize(rcxres, rcyres);
            Xslice = 0;
            for (int coil = 0; coil < data.Num_Coils; coil++) {
              Xslice +=
                  norm(XX(0, 0, coil)(all, all, RR(0, 0, 0).length(2) / 2));
            }
            Xslice = sqrt(Xslice);
            ArrayWriteAppend(Xslice, "X_mag.dat");
          }
        }

        // ---------------------------------
        //   Clear
        // ---------------------------------

        iteration_timer.tic();
        lrankcoil.update_threshold(XX, 2, iteration);
        cout << "Get thresh took " << iteration_timer << " s" << endl;

        iteration_timer.tic();
        lrankcoil.thresh(XX, 2);
        cout << "Thresh took " << iteration_timer << " s" << endl;

        // Export X slice
        {
          Array<float, 2> Xslice;
          Xslice.setStorage(ColumnMajorArray<2>());
          Xslice.resize(rcxres, rcyres);
          Xslice = 0;
          for (int coil = 0; coil < data.Num_Coils; coil++) {
            Xslice += norm(XX(0, 0, coil)(all, all, RR(0, 0, 0).length(2) / 2));
          }
          Xslice = sqrt(Xslice);
          ArrayWriteAppend(Xslice, "X_mag.dat");
        }
      }  // Iteration

      {
        HDF5 clearRAW("ClearRaw.h5", "w");

        for (int t = 0; t < XX.length(firstDim); t++) {
          for (int e = 0; e < XX.length(secondDim); e++) {
            for (int coil = 0; coil < data.Num_Coils; coil++) {
              char name[1024];
              sprintf(name, "IM_T%03d_E%03d_C%03d", t, e, coil);
              clearRAW.AddH5Array("IMAGES", name, XX(t, e, coil));
            }
          }
        }
      }

      lrankcoil.combine(XX, X);

    } break;

    case (ADMM): {
      // ------------------------------------
      // Alternating Direction of Multipliers
      //   0.5||Ex-d||2  + gamma*||y||1 + rho*||x-y||2 + L'*(x-y)

      // CG Structures
      Array<Array<complex<float>, 3>, 2> R = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes);
      Array<Array<complex<float>, 3>, 2> P = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes);
      Array<Array<complex<float>, 3>, 2> LHS = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes);

      // Augmented Lagrangian Structures
      Array<Array<complex<float>, 3>, 2> Y = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes);
      Array<Array<complex<float>, 3>, 2> L = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes);

      float cg_thresh = 1e-9;

      // Initialize
      // X = 0;
      // Y = 0;
      // L = 0;
      for (int e = 0; e < rcencodes; e++) {
        for (int t = 0; t < Nt; t++) {
          X(t, e) = complex<float>(0.0, 0.0);
          Y(t, e) = complex<float>(0.0, 0.0);
          L(t, e) = complex<float>(0.0, 0.0);
          LHS(t, e) = complex<float>(0.0, 0.0);
        }
      }

      // Iterate
      for (int admm_iteration = 0; admm_iteration < admm_max_iter; admm_iteration++) {
        //
        //  Solve Subproblem  1
        //  	0.5||Ex-d||2  + gamma*||y||1 + rho*||x-y||2 + L'*(x-y)
        //      where y and L are fixed paramaters. This reduces to minimizing
        // 		0.5||Ex-d||2  + rho*||x-y||2 + L'*(x-y)

        // (E'E + rho*I)*x = E'd + rho*y - L;

        // First Calculate E'd
        cout << "ADMM Inner :: LHS Calculation" << endl;
        for (int e = 0; e < rcencodes; e++) {
          for (int t = 0; t < Nt; t++) {
            R(t, e) = complex<float>(0.0, 0.0);

            int act_t = times(t);
            int store_t = times_store(t);

            // Temporal weighting
            if (pregate_data_flag) {
              TimeWeight.reference(data.kw(e));
            } else if (reset_dens) {
              TimeWeight.resize(data.kw(e).shape());
              TimeWeight = 1.0;
              gate.weight_data(TimeWeight, e, data.kx(e), data.ky(e),
                               data.kz(e), act_t, GATING::ITERATIVE,
                               frame_type);
            } else {
              TimeWeight.resize(data.kw(e).shape());
              TimeWeight = data.kw(e);
              gate.weight_data(TimeWeight, e, data.kx(e), data.ky(e),
                               data.kz(e), act_t, GATING::ITERATIVE,
                               frame_type);
            }

            // Images
            for (int coil = 0; coil < data.Num_Coils; coil++) {
              Array<complex<float>, 3> diff_data(data.kx(e).shape(),
                                                 ColumnMajorArray<3>());
              diff_data = complex<float>(0.0, 0.0);

              // Ex
              gridding.backward(X(store_t, e), smaps(coil), diff_data,
                                data.kx(e), data.ky(e), data.kz(e), TimeWeight);

              // d - Ex
              diff_data = data.kdata(e, coil) - diff_data;

              // E'd - E'Ex
              gridding.forward(R(store_t, e), smaps(coil), diff_data,
                               data.kx(e), data.ky(e), data.kz(e), TimeWeight);

            }  // Coils

            // E'd - ( E'Ex + rho*Ix)
            R(store_t, e) -= admm_rho * X(store_t, e);

            // RHS Values ( E'd + rho*Y - L ) - (E' + rho)*X
            R(store_t, e) += admm_rho * Y(store_t, e);
            R(store_t, e) -= L(store_t, e);

            // Initiialize P
            P(store_t, e) = R(store_t, e);
          }
        }

        // Now Iterate
        cout << "ADMM Inner :: Iterate " << endl;
        for (int cg_iteration = 0; cg_iteration < max_iter; cg_iteration++) {
          tictoc iteration_timer;
          iteration_timer.tic();

          // E'Ex
          for (int e = 0; e < rcencodes; e++) {
            for (int t = 0; t < Nt; t++) {
              int act_t = times(t);
              // int store_t = times_store(t);

              LHS(t, e) = complex<float>(0.0, 0.0);
              T.tic();

              // Temporal weighting
              if (pregate_data_flag) {
                TimeWeight.reference(data.kw(e));
              } else if (reset_dens) {
                TimeWeight.resize(data.kw(e).shape());
                TimeWeight = 1.0;
                gate.weight_data(TimeWeight, e, data.kx(e), data.ky(e),
                                 data.kz(e), act_t, GATING::ITERATIVE,
                                 frame_type);
              } else {
                TimeWeight.resize(data.kw(e).shape());
                TimeWeight = data.kw(e);
                gate.weight_data(TimeWeight, e, data.kx(e), data.ky(e),
                                 data.kz(e), act_t, GATING::ITERATIVE,
                                 frame_type);
              }

              for (int coil = 0; coil < data.Num_Coils; coil++) {
                // Storage for (Ex-d) - dynamic for variable size
                Array<complex<float>, 3> diff_data(data.kx(e).shape(),
                                                   ColumnMajorArray<3>());
                diff_data = complex<float>(0.0, 0.0);

                // E'Ex
                gridding.backward(P(t, e), smaps(coil), diff_data, data.kx(e),
                                  data.ky(e), data.kz(e), TimeWeight);
                gridding.forward(LHS(t, e), smaps(coil), diff_data, data.kx(e),
                                 data.ky(e), data.kz(e), TimeWeight);
              }  // Coils

              // Values
              LHS(t, e) += admm_rho * P(t, e);

            }  // t
          }  // e

          //----------------------------------------------------
          //  Now perform gradient update
          // ---------------------------------------------------

          complex<float> sum_R0_R0(0.0, 0.0);
          complex<float> sum_R_R(0.0, 0.0);
          complex<float> sum_P_LHS(0.0, 0.0);

          // Calc R'R and P'*LHS
          for (int e = 0; e < rcencodes; e++) {
            for (int t = 0; t < Nt; t++) {
              for (int k = 0; k < rczres; k++) {
                for (int j = 0; j < rcyres; j++) {
                  for (int i = 0; i < rcxres; i++) {
                    sum_R0_R0 += norm(R(t, e)(i, j, k));
                    sum_P_LHS += conj(P(t, e)(i, j, k)) * LHS(t, e)(i, j, k);
                  }
                }
              }
            }
          }
          complex<float> scale = sum_R0_R0 / sum_P_LHS;

          // Take step size
          for (int e = 0; e < rcencodes; e++) {
            for (int t = 0; t < Nt; t++) {
              for (int k = 0; k < rczres; k++) {
                for (int j = 0; j < rcyres; j++) {
                  for (int i = 0; i < rcxres; i++) {
                    X(t, e)
                    (i, j, k) += (scale * (P(t, e)(i, j, k)));
                    R(t, e)
                    (i, j, k) -= (scale * (LHS(t, e)(i, j, k)));
                    sum_R_R += norm(R(t, e)(i, j, k));
                  }
                }
              }
            }
          }

          cout << "       Sum R'R = " << sum_R_R << endl;
          complex<float> scale2 = sum_R_R / sum_R0_R0;

          // Take step size
          for (int e = 0; e < rcencodes; e++) {
            for (int t = 0; t < Nt; t++) {
              for (int k = 0; k < rczres; k++) {
                for (int j = 0; j < rcyres; j++) {
                  for (int i = 0; i < rcxres; i++) {
                    P(t, e)
                    (i, j, k) = R(t, e)(i, j, k) + (scale2 * P(t, e)(i, j, k));
                  }
                }
              }
            }
          }

          {
            Array<complex<float>, 2> Xslice =
                X(0, 0)(all, all, X(0, 0).length(2) / 2);
            ArrayWriteMagAppend(Xslice, "X_mag.dat");
          }

          if (abs(sum_R_R) < cg_thresh) {
            break;
          }

        }  // Inner CG Iterations

        //
        //  Solve Subproblem  2
        //  	0.5||Ex-d||2  + gamma*||y||1 + rho*||x-y||2 + L'*(x-y)
        //      where x and L are fixed paramaters. This reduces to
        //		gamma*||y||1 + rho*||x-y||2 + L'*(x-y)

        // Take step size
        float alpha = 0;
        for (int e = 0; e < rcencodes; e++) {
          for (int t = 0; t < Nt; t++) {
            for (int k = 0; k < rczres; k++) {
              for (int j = 0; j < rcyres; j++) {
                for (int i = 0; i < rcxres; i++) {
                  Y(t, e)
                  (i, j, k) = alpha * X(t, e)(i, j, k) +
                              (1 - alpha) * Y(t, e)(i, j, k) +
                              L(t, e)(i, j, k) / admm_rho;

                  // Also update dual
                  L(t, e)
                  (i, j, k);
                }
              }
            }
          }
        }

        {
          Array<complex<float>, 2> Xslice =
              Y(0, 0)(all, all, X(0, 0).length(2) / 2);
          ArrayWriteMagAppend(Xslice, "Y_mag.dat");
        }

        if (admm_iteration == 0) {
          admm_gamma *= max(abs(X(0, 0)));
        }
        softthresh.threshold_type = TH_FIXED;
        softthresh.thresh = admm_gamma / admm_rho;
        L1_threshold(Y);

        {
          Array<complex<float>, 2> Xslice =
              Y(0, 0)(all, all, X(0, 0).length(2) / 2);
          ArrayWriteMagAppend(Xslice, "Y_mag.dat");
        }

        //
        // Update the lagrangian multiplier
        //
        for (int e = 0; e < rcencodes; e++) {
          for (int t = 0; t < Nt; t++) {
            for (int k = 0; k < rczres; k++) {
              for (int j = 0; j < rcyres; j++) {
                for (int i = 0; i < rcxres; i++) {
                  L(t, e)
                  (i, j, k) = L(t, e)(i, j, k) +
                              admm_rho * (X(t, e)(i, j, k) - Y(t, e)(i, j, k));
                }
              }
            }
          }
        }

        {
          Array<complex<float>, 2> Xslice =
              L(0, 0)(all, all, X(0, 0).length(2) / 2);
          ArrayWriteMagAppend(Xslice, "L_mag.dat");
        }

      }  // iteration

    } break;

    case (CG): {
      // ------------------------------------
      // Conjugate Gradient Recon
      //   -Uses much more memory than gradient descent but is faster (both
      //   convergence + per iteration)
      // ------------------------------------

      // Structures
      Array<Array<complex<float>, 3>, 2> R = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes);
      Array<Array<complex<float>, 3>, 2> P = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes);
      Array<Array<complex<float>, 3>, 2> LHS = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes);

      // First Calculate E'd
      cout << "LHS Calculation" << endl;
      for (int e = 0; e < rcencodes; e++) {
        for (int t = 0; t < Nt; t++) {
          int act_t = times(t);
          int store_t = times_store(t);
          int act_e = (pregate_data_flag) ? (e * Nt + act_t) : (e);

          // Temporal weighting
          if (pregate_data_flag) {
            TimeWeight.reference(data.kw(act_e));
          } else if (reset_dens) {
            TimeWeight.resize(data.kw(act_e).shape());
            TimeWeight = 1.0;
            gate.weight_data(TimeWeight, act_e,
                             data.kx(act_e), data.ky(act_e), data.kz(act_e),
                             act_t, GATING::ITERATIVE, frame_type);
          } else {
            TimeWeight.resize(data.kw(act_e).shape());
            TimeWeight = data.kw(act_e);
            gate.weight_data(TimeWeight, act_e,
                             data.kx(act_e), data.ky(act_e), data.kz(act_e),
                             act_t, GATING::ITERATIVE, frame_type);
          }

          // E'd
          if (this->parallel_coils) {
            Array<Array<complex<float>, 3>, 1> TEMP_KDATA = data.kdata(act_e, Range::all());
            gridding_CoilThreaded.forward(R(store_t, e), smaps, TEMP_KDATA, data.kx(act_e), data.ky(act_e), data.kz(act_e), TimeWeight);
          } else {
            for (int coil = 0; coil < data.Num_Coils; coil++) {
              gridding.forward(R(store_t, e), smaps(coil), data.kdata(act_e, coil), data.kx(act_e), data.ky(act_e), data.kz(act_e), TimeWeight);
            }
          }

          // Initialize P
          P(store_t, e) = R(store_t, e);
        }
      }

      // Now Iterate
      cout << "Iterate" << endl;
      // double error0=0.0;
      for (int iteration = 0; iteration < max_iter; iteration++) {
        tictoc iteration_timer;
        iteration_timer.tic();

        // E'Ex
        for (int e = 0; e < rcencodes; e++) {
          for (int t = 0; t < Nt; t++) {
            int act_t = times(t);
            int store_t = times_store(t);
            int act_e = (pregate_data_flag) ? (e * Nt + act_t) : (e);

            LHS(store_t, e) = complex<float>(0.0, 0.0);
            T.tic();

            // Temporal weighting
            if (pregate_data_flag) {
              TimeWeight.reference(data.kw(act_e));
            } else if (reset_dens) {
              TimeWeight.resize(data.kw(act_e).shape());
              TimeWeight = 1.0;
              gate.weight_data(TimeWeight, act_e, data.kx(act_e), data.ky(act_e), data.kz(act_e), act_t, GATING::ITERATIVE, frame_type);
            } else {
              TimeWeight.resize(data.kw(act_e).shape());
              TimeWeight = data.kw(act_e);
              gate.weight_data(TimeWeight, act_e, data.kx(act_e), data.ky(act_e), data.kz(act_e), act_t, GATING::ITERATIVE, frame_type);
            }

            if (this->parallel_coils) {
              // K-space to Image
              Array<Array<complex<float>, 3>, 1> diff_data = Alloc4DContainer<complex<float> >(data.kx(act_e).length(firstDim), data.kx(act_e).length(secondDim), data.kx(act_e).length(thirdDim), data.Num_Coils);

              // E'Ex
              gridding_CoilThreaded.backward(P(store_t, e), smaps, diff_data, data.kx(act_e), data.ky(act_e), data.kz(act_e), TimeWeight);
              gridding_CoilThreaded.forward(LHS(store_t, e), smaps, diff_data, data.kx(act_e), data.ky(act_e), data.kz(act_e), TimeWeight);

            } else {
              for (int coil = 0; coil < data.Num_Coils; coil++) {
                // Storage for (Ex-d) - dynamic for variable size
                Array<complex<float>, 3> diff_data(data.kx(act_e).shape(), ColumnMajorArray<3>());
                diff_data = 0;

                // E'Ex
                gridding.backward(P(store_t, e), smaps(coil), diff_data, data.kx(act_e), data.ky(act_e), data.kz(act_e), TimeWeight);
                gridding.forward(LHS(store_t, e), smaps(coil), diff_data, data.kx(act_e), data.ky(act_e), data.kz(act_e), TimeWeight);
              }  // Coils
            }

            // L2 R'R
            if (iteration > 0) {
            }

            cout << "\r" << e << "," << t << "took " << T << "s" << flush;

          }  // t
        }  // e

        // Regularization
        if ((iteration == 0) && (l2reg.lambda > 0)) {
          float sum_P_P = 0.0;
          float sum_L_L = 0.0;
          for (int e = 0; e < LHS.length(secondDim); e++) {
            for (int t = 0; t < LHS.length(firstDim); t++) {
              sum_P_P += sum(norm(P(t, e)));
              sum_L_L += sum(norm(LHS(t, e)));
            }
          }

          l2reg.reg_scale = l2reg.lambda * sqrt(sum_L_L / sum_P_P);
          cout << "L2 Scale = " << l2reg.reg_scale << endl;
        }

        for (int e = 0; e < LHS.length(firstDim); e++) {
          for (int t = 0; t < LHS.length(secondDim); t++) {
            l2reg.regularize(LHS(t, e), P(t, e));
          }
        }

        //----------------------------------------------------
        //  Now perform gradient update
        // ---------------------------------------------------

        Array<complex<double>, 1> Asum_R0_R0(rczres);
        Array<complex<double>, 1> Asum_R_R(rczres);
        Array<complex<double>, 1> Asum_P_LHS(rczres);
        Asum_R0_R0 = complex<double>(0.0, 0.0);
        Asum_R_R = complex<double>(0.0, 0.0);
        Asum_P_LHS = complex<double>(0.0, 0.0);

        complex<double> sum_R0_R0(0.0, 0.0);
        complex<double> sum_R_R(0.0, 0.0);
        complex<double> sum_P_LHS(0.0, 0.0);

        // Calc R'R and P'*LHS
        for (int e = 0; e < rcencodes; e++) {
          for (int t = 0; t < Nt; t++) {
#pragma omp parallel for
            for (int k = 0; k < rczres; k++) {
              for (int j = 0; j < rcyres; j++) {
                for (int i = 0; i < rcxres; i++) {
                  Asum_R0_R0(k) += norm(R(t, e)(i, j, k));
                  Asum_P_LHS(k) += conj(P(t, e)(i, j, k)) * LHS(t, e)(i, j, k);
                }
              }
            }
          }
        }
        sum_R0_R0 = sum(Asum_R0_R0);
        sum_P_LHS = sum(Asum_P_LHS);
        complex<double> dscale = sum_R0_R0 / sum_P_LHS;
        complex<float> scale = complex<float>(real(dscale), imag(dscale));

        // Take step size
        for (int e = 0; e < rcencodes; e++) {
          for (int t = 0; t < Nt; t++) {
#pragma omp parallel for
            for (int k = 0; k < rczres; k++) {
              for (int j = 0; j < rcyres; j++) {
                for (int i = 0; i < rcxres; i++) {
                  X(t, e)
                  (i, j, k) += (scale * (P(t, e)(i, j, k)));
                  R(t, e)
                  (i, j, k) -= (scale * (LHS(t, e)(i, j, k)));
                  Asum_R_R(k) += norm(R(t, e)(i, j, k));
                }
              }
            }
          }
        }
        sum_R_R = sum(Asum_R_R);

        cout << "Sum R'R = " << sum_R_R << endl;
        complex<double> Dscale2 = sum_R_R / sum_R0_R0;
        complex<float> scale2 = complex<float>(real(Dscale2), imag(Dscale2));

        // Take step size
        for (int e = 0; e < rcencodes; e++) {
          for (int t = 0; t < Nt; t++) {
#pragma omp parallel for
            for (int k = 0; k < rczres; k++) {
              for (int j = 0; j < rcyres; j++) {
                for (int i = 0; i < rcxres; i++) {
                  P(t, e)
                  (i, j, k) = R(t, e)(i, j, k) + (scale2 * P(t, e)(i, j, k));
                }
              }
            }
          }
        }

        // Export X slice
        {
          Array<complex<float>, 2> Xslice = X(0, 0)(all, all, X(0, 0).length(2) / 2);
          ArrayWriteMagAppend(Xslice, "X_mag.dat");

          Array<complex<float>, 2> LHSslice = LHS(0, 0)(all, all, LHS(0, 0).length(2) / 2);
          ArrayWriteMagAppend(LHSslice, "LHS_mag.dat");

          Array<complex<float>, 2> Pslice = X(0, 0)(all, all, P(0, 0).length(2) / 2);
          ArrayWriteMagAppend(Pslice, "P_mag.dat");

          Array<complex<float>, 2> Rslice = R(0, 0)(all, all, R(0, 0).length(2) / 2);
          ArrayWriteMagAppend(Rslice, "R_mag.dat");

          if (X.numElements() > 1) {
            int count = 0;
            for (Array<Array<complex<float>, 3>, 2>::iterator miter = X.begin(); miter != X.end(); miter++) {
              Array<complex<float>, 2> Xf = (*miter)(all, all, X(0, 0).length(2) / 2);
              if (count == 0) {
                ArrayWriteMag(Xf, "X_frames.dat");
                ArrayWritePhase(Xf, "X_frames.dat.phase");
              } else {
                ArrayWriteMagAppend(Xf, "X_frames.dat");
                ArrayWritePhaseAppend(Xf, "X_frames.dat.phase");
              }
              count++;
            }
          }
        }

      }  // iteration

      // Correct for the intensity bias
      correct_intensity_bias(X);

    } break;

    case (IST):
    case (FISTA): {
      // ------------------------------------
      // Iterative Soft Thresholding  x(n+1)=  thresh(   x(n) - E*(Ex(n) - d)  )
      //  Designed to not use memory
      // Uses gradient descent x(n+1) = x(n) - ( R'R ) / ( R'E'E R) * Grad  [ R
      // = E'(Ex-d)]
      // ------------------------------------

      // Arrays for FISTA
      Array<Array<complex<float>, 3>, 2> X_old;
      if (recon_type == FISTA) {
        cout << "Alloc Fista Matrix" << endl;
        Array<Array<complex<float>, 3>, 2> Temp = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes);
        X_old.reference(Temp);
      }

      // Residue
      Array<Array<complex<float>, 3>, 2> R = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, Nt, rcencodes);

      // Temp variable for E'ER
      Array<complex<float>, 3> P(rcxres, rcyres, rczres, ColumnMajorArray<3>());

      // Class for gradient descent step size
      complex<float> step_size;

      //
      // Max Eigen calculations
      //

      if (iterative_step_type == STEP_MAXEIG) {
        // Initialize x to random number
        std::cout << "Max Eigen Calc for Step size" << std::endl;
        std::cout << "Seeding X" << std::endl;

        int max_encode = rcencodes;
        int max_frame = Nt;
        if (fast_maxeig) {
          max_encode = 1;
          max_frame = 1;
        }

        for (int e = 0; e < max_encode; e++) {
          for (int t = 0; t < max_frame; t++) {
#pragma omp parallel for
            for (int k = 0; k < rczres; k++) {
              for (int j = 0; j < rcyres; j++) {
                for (int i = 0; i < rcxres; i++) {
                  float tmp_r = std::rand() / RAND_MAX - 0.5;
                  float tmp_i = std::rand() / RAND_MAX - 0.5;
                  X(t, e)
                  (i, j, k) = complex<float>(tmp_r, tmp_i);
                }
              }
            }
          }
        }

        // Now use iterations through power method
        for (int iteration = 0; iteration < max_eigen_iterations; iteration++) {
          // Get Residue
          for (int e = 0; e < max_encode; e++) {
            for (int t = 0; t < max_frame; t++) {
              R(t, e) = 0.0;
            }
          }

          // Get EH*E*x
          for (int e = 0; e < max_encode; e++) {
            for (int t = 0; t < max_frame; t++) {
              int act_t = times(t);
              int store_t = times_store(t);
              int act_e = (pregate_data_flag) ? (e * Nt + act_t) : (e);

              // Get Sub-Arrays for Encoding
              Array<float, 3> kxE = data.kx(act_e);
              Array<float, 3> kyE = data.ky(act_e);
              Array<float, 3> kzE = data.kz(act_e);
              Array<float, 3> kwE = data.kw(act_e);

              // Temporal weighting
              if (pregate_data_flag) {
                TimeWeight.reference(kwE);
              } else if (reset_dens) {
                TimeWeight.resize(kwE.shape());
                TimeWeight = 1.0;
                Array<float, 3>::iterator titer = TimeWeight.begin();
                Array<float, 3>::iterator kiter = kwE.begin();
                for (; (titer != TimeWeight.end()) && (kiter != kwE.end()); titer++, kiter++) {
                  if ((*kiter) < 0.1) {
                    *titer = *kiter;
                  } else {
                    *titer = 0.1;
                  }
                }
                gate.weight_data(TimeWeight, e, kxE, kyE, kzE, act_t, GATING::ITERATIVE, frame_type);
              } else {
                TimeWeight.resize(kwE.shape());
                TimeWeight = kwE;
                gate.weight_data(TimeWeight, e, kxE, kwE, kzE, act_t, GATING::ITERATIVE, frame_type);
              }

              // Differences (Ex)
              Array<complex<float>, 3> diff_data(kxE.shape(), ColumnMajorArray<3>());
              Array<Array<complex<float>, 3>, 1> diff_data_all;

              if (this->parallel_coils) {
                // K-space to Image
                Array<Array<complex<float>, 3>, 1> TEMP_KDATA = data.kdata(act_e, Range::all());
                Array<Array<complex<float>, 3>, 1> TEMP_POINT = Alloc4DContainer<complex<float> >(kxE.length(firstDim), kxE.length(secondDim), kxE.length(thirdDim), data.Num_Coils);
                diff_data_all.reference(TEMP_POINT);

                // Ex
                gridding_CoilThreaded.backward(X(store_t, e), smaps, diff_data_all, kxE, kyE, kzE, TimeWeight);

                // E'Ex
                gridding_CoilThreaded.forward(R(store_t, e), smaps, diff_data_all, kxE, kyE, kzE, TimeWeight);
              } else {
                for (int coil = 0; coil < data.Num_Coils; coil++) {
                  // Ex
                  gridding.backward(X(store_t, e), smaps(coil), diff_data, kxE, kyE, kzE, TimeWeight);

                  // E'Ex
                  gridding.forward(R(store_t, e), smaps(coil), diff_data, kxE, kyE, kzE, TimeWeight);
                }  // Coils
              }
            }
          }

          // Power Methods
          double max_eig = 0.0;
          for (int e = 0; e < max_encode; e++) {
            for (int t = 0; t < max_frame; t++) {
              max_eig += ArrayEnergy(R(t, e));
            }
          }
          max_eig = sqrt(max_eig);

          complex<float> scale = complex<float>(1. / max_eig, 0.0);
          for (int e = 0; e < max_encode; e++) {
            for (int t = 0; t < max_frame; t++) {
              X(t, e) = scale * R(t, e);
            }
          }
          cout << "Iter = " << iteration << " Max Eig = " << max_eig << endl;
          step_size = complex<float>(1. / max_eig, 0.0);
        }

        // Scale kweights instead of step size
        if (this->image_scale_normalization) {
          for (Array<Array<float, 3>, 1>::iterator riter = data.kw.begin(); riter != data.kw.end(); riter++) {
            *riter *= abs(step_size);
          }
          step_size = complex<float>(1.0, 0.0);
        }
      }

      // Reset X after max_eigen calculation
      for (Array<Array<complex<float>, 3>, 2>::iterator riter = X.begin(); riter != X.end(); riter++) {
        *riter = 0;
      }

      cout << "Iterate" << endl;
      double error0 = 0.0;
      for (int iteration = 0; iteration < max_iter; iteration++) {
        tictoc iteration_timer;
        iteration_timer.tic();
        cout << "\nIteration = " << iteration << endl;

        // Zero this for Cauchy set size
        complex<double> scale_RhP(0, 0);

        // Get Residue
        for (Array<Array<complex<float>, 3>, 2>::iterator riter = R.begin(); riter != R.end(); riter++) {
          *riter = 0;
        }

        // Copy FISTA update
        if (recon_type == FISTA) {
          fista_update(X, X_old, iteration);
        }

        cout << "\tGradient Calculation" << endl;

        for (int e = 0; e < rcencodes; e++) {
          for (int t = 0; t < Nt; t++) {
            int act_t = times(t);
            int store_t = times_store(t);
            int act_e = (pregate_data_flag) ? (e * Nt + act_t) : (e);

            // Get Sub-Arrays for Encoding
            Array<float, 3> kxE = data.kx(act_e);
            Array<float, 3> kyE = data.ky(act_e);
            Array<float, 3> kzE = data.kz(act_e);
            Array<float, 3> kwE = data.kw(act_e);

            T.tic();

            // Temporal weighting
            if (pregate_data_flag) {
              TimeWeight.reference(kwE);
            } else if (reset_dens) {
              TimeWeight.resize(kwE.shape());
              TimeWeight = 1.0;
              Array<float, 3>::iterator titer = TimeWeight.begin();
              Array<float, 3>::iterator kiter = kwE.begin();
              for (; (titer != TimeWeight.end()) && (kiter != kwE.end()); titer++, kiter++) {
                if ((*kiter) < 0.1) {
                  *titer = *kiter;
                } else {
                  *titer = 0.1;
                }
              }
              gate.weight_data(TimeWeight, e, kxE, kyE, kzE, act_t, GATING::ITERATIVE, frame_type);
            } else {
              TimeWeight.resize(kwE.shape());
              TimeWeight = kwE;
              gate.weight_data(TimeWeight, e, kxE, kwE, kzE, act_t, GATING::ITERATIVE, frame_type);
            }

            // Differences (Ex-d)
            Array<complex<float>, 3> diff_data(kxE.shape(), ColumnMajorArray<3>());
            Array<Array<complex<float>, 3>, 1> diff_data_all;

            if (this->parallel_coils) {
              // K-space to Image
              Array<Array<complex<float>, 3>, 1> TEMP_KDATA = data.kdata(act_e, Range::all());
              Array<Array<complex<float>, 3>, 1> TEMP_POINT = Alloc4DContainer<complex<float> >(kxE.length(firstDim), kxE.length(secondDim), kxE.length(thirdDim), data.Num_Coils);
              diff_data_all.reference(TEMP_POINT);

              // Ex - d
              gridding_CoilThreaded.backward_residual(X(store_t, e), smaps, diff_data_all, kxE, kyE, kzE, TimeWeight, TEMP_KDATA);

              // E'(Ex-d)
              gridding_CoilThreaded.forward(R(store_t, e), smaps, diff_data_all, kxE, kyE, kzE, TimeWeight);
            } else {
              for (int coil = 0; coil < data.Num_Coils; coil++) {
                // Ex - d
                gridding.backward_residual(X(store_t, e), smaps(coil), diff_data, kxE, kyE, kzE, TimeWeight, data.kdata(act_e, coil));

                // E'(Ex-d)
                gridding.forward(R(store_t, e), smaps(coil), diff_data, kxE, kyE, kzE, TimeWeight);

              }  // Coils
            }

            // L2
            if (iteration > 0) {
              l2reg.regularize(R(store_t, e), X(store_t, e));
            }

            // Now Get Scale factor (for Cauchy-Step Size)
            if ((iteration < this->cauchy_update_number) && (this->iterative_step_type == STEP_CAUCHY)) {
              P = 0;
              if (this->parallel_coils) {
                // EE'(Ex-d)
                gridding_CoilThreaded.backward(R(store_t, e), smaps, diff_data_all, kxE, kyE, kzE, TimeWeight);

                // E'EE'(Ex-d)
                gridding_CoilThreaded.forward(P, smaps, diff_data_all, kxE, kyE, kzE, TimeWeight);
              } else {
                for (int coil = 0; coil < data.Num_Coils; coil++) {
                  // EE'(Ex-d)
                  gridding.backward(R(store_t, e), smaps(coil), diff_data, kxE, kyE, kzE, TimeWeight);

                  // E'EE'(Ex-d)
                  gridding.forward(P, smaps(coil), diff_data, kxE, kyE, kzE, TimeWeight);
                }  // Coils
              }

              // TV of Image
              if (iteration > 0) {
                l2reg.regularize(P, R(store_t, e));
              }

              P *= conj(R(store_t, e));

              for (Array<complex<float>, 3>::iterator riter = P.begin(); riter != P.end(); riter++) {
                scale_RhP += complex<double>(real(*riter), imag(*riter));
              }
            }
            cout << "\r" << e << "," << t << "took " << T << "s" << flush;
          }  // Time
        }  // Encode

        // Scale the image to set aproximate value to 1
        if (this->image_scale_normalization && (iteration == 0)) {
          // Get max value
          float max_image_value = 0.0;
          for (Array<Array<complex<float>, 3>, 2>::iterator riter = R.begin(); riter != R.end(); riter++) {
            float max_current = max(abs(*riter));
            max_image_value = (max_current > max_image_value) ? (max_current) : (max_image_value);
          }
          complex<float> image_scale(1. / max_image_value, 0.0);
          std::cout << "Scaling Image to max ~1, Max is " << max_image_value << " scale to " << image_scale << std::endl;

          // Scale image
          for (Array<Array<complex<float>, 3>, 2>::iterator riter = R.begin(); riter != R.end(); riter++) {
            *riter *= image_scale;
          }

          // Scale data
          for (Array<Array<complex<float>, 3>, 2>::iterator riter = data.kdata.begin(); riter != data.kdata.end(); riter++) {
            *riter *= image_scale;
          }
        }

        // Get Scaling Factor R'P / R'R
        cout << endl
             << "Calc residue" << endl
             << flush;
        complex<double> scale_RhR = 0.0;
        for (Array<Array<complex<float>, 3>, 2>::iterator riter = R.begin(); riter != R.end(); riter++) {
          scale_RhR += complex<double>(ArrayEnergy(*riter), 0.0);
        }

        // Error check
        if (iteration == 0) {
          error0 = abs(scale_RhR);
        }
        cout << "Residue Energy =" << scale_RhR << " ( " << (abs(scale_RhR) / error0) << " % )  " << endl;
        cout << "RhP = " << scale_RhP << endl;

        // Export R
        export_slice(R(0, 0), "R.dat");

        // Step in direction
        if ((iteration < this->cauchy_update_number) && (this->iterative_step_type == STEP_CAUCHY)) {
          complex<double> scale_double = (scale_RhR / scale_RhP);
          step_size = complex<float>(real(scale_double), imag(scale_double));
        }

        cout << "Scale = " << step_size << endl
             << flush;
        for (int e = 0; e < rcencodes; e++) {
          for (int t = 0; t < Nt; t++) {
            R(t, e) *= step_size;
            X(t, e) -= R(t, e);
          }
        }
        cout << "Took " << iteration_timer << " s " << endl;

        // Error check
        if (iteration == 0) {
          l2reg.set_scale((float)abs(scale_RhR), X);
          cout << "L2 Scale = " << l2reg.reg_scale << endl;
        }

        // Export X slice
        export_slice(X(0, 0), "X_mag.dat");

        if (lranktime.clear_alpha_time > 0.0) {
          lranktime.update_threshold(X, 0, iteration);
          lranktime.thresh(X, 0);
          export_slice(X(0, 0), "X_mag.dat");
        }

        if (lranktime.clear_alpha_encode > 0.0) {
          lranktime.update_threshold(X, 1, iteration);
          lranktime.thresh(X, 1);
          export_slice(X(0, 0), "X_mag.dat");
        }

        if (X.numElements() > 1) {
          int count = 0;
          for (Array<Array<complex<float>, 3>, 2>::iterator miter = X.begin(); miter != X.end(); miter++) {
            Array<complex<float>, 2> Xf = (*miter)(all, all, X(0, 0).length(2) / 2);
            if (count == 0) {
              ArrayWriteMag(Xf, "X_frames.dat");
              ArrayWritePhase(Xf, "X_frames.dat.phase");
            } else {
              ArrayWriteMagAppend(Xf, "X_frames.dat");
              ArrayWritePhaseAppend(Xf, "X_frames.dat.phase");
            }
            count++;
          }
        }

        // ------------------------------------
        // Soft thresholding operation (need to add transform control)
        // ------------------------------------
        if (softthresh.getThresholdMethod() != TH_NONE) {
          correct_intensity_bias(X);
          L1_threshold(X);
          apply_intensity_bias(X);
        }
        export_slice(X(0, 0), "X_mag.dat");

      }  // Iteration

      // Correct
      correct_intensity_bias(X);

    } break;

  }  // Recon Type

  cout << "Recon was completed successfully " << endl;
  return (X);
}

void RECON::fista_update(Array<Array<complex<float>, 3>, 2> &X, Array<Array<complex<float>, 3>, 2> &X_old, int iteration) {
  float A = 1.0 + (iteration - 1.0) / (iteration + 1.0);
  float B = -(iteration - 1.0) / (iteration + 1.0);

  int Nx = X_old(0).length(firstDim);
  int Ny = X_old(0).length(secondDim);
  int Nz = X_old(0).length(thirdDim);
  int Ne = X_old.length(secondDim);
  int Nt = X_old.length(firstDim);

  cout << "fista_update: Matrix Size = " << Nx << " x " << Ny << " x " << Nz
       << " x " << Nt << " x " << Ne << endl;

  // Get Update
  for (int e = 0; e < Ne; e++) {
    for (int t = 0; t < Nt; t++) {
      Array<complex<float>, 3> XX = X(t, e);
      Array<complex<float>, 3> XX_old = X_old(t, e);

      for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
          for (int i = 0; i < Nx; i++) {
            // cout << "Pos = " << i << "," << j << "," << k << "," << t << ","
            // << e << endl;
            complex<float> Xn0 = XX_old(i, j, k);
            complex<float> Xn1 = XX(i, j, k);
            XX(i, j, k) = A * Xn1 + B * Xn0;
            XX_old(i, j, k) = Xn1;
          }
        }
      }  // Spatial
    }
  }
  cout << "Done with FISTA" << endl
       << flush;
}

double RECON::kspace_residual(MRI_DATA &data) {
  // Matlab like timer (openmp code base)
  tictoc T;

  // Setup Gridding + FFT Structure
  gridFFT Kgridding;
  Kgridding.overgrid = 1.0;
  Kgridding.kernel_type = TRIANGLE_KERNEL;
  Kgridding.dwinX = 1.0;
  Kgridding.dwinY = 1.0;
  Kgridding.dwinZ = 1.0;
  Kgridding.precalc_gridding(32, 32, 32, data);

  // Class for gradient descent step size
  complex<float> step_size;

  double residual = 0.0;
  for (int e = 0; e < data.kx.length(firstDim); e++) {
    // Get Sub-Arrays for Encoding
    Array<float, 3> kxE = data.kx(e);
    Array<float, 3> kyE = data.ky(e);
    Array<float, 3> kzE = data.kz(e);
    Array<float, 3> kwE = data.kw(e);

    // Differences (Ex-d)
    Array<complex<float>, 3> diff_data(kxE.shape(), ColumnMajorArray<3>());

    for (int coil = 0; coil < data.Num_Coils; coil++) {
      // Grid data
      T.tic();
      Kgridding.k3d_grid = 0;  // Zero K-Space
      Kgridding.chop_grid_forward(data.kdata(e, coil), kxE, kyE, kzE,
                                  kwE);  // Grid data to K-Space
      // cout << "Forward took " << T << endl;

      // Inverse
      T.tic();
      Array<complex<float>, 3> temp;
      Kgridding.chop_grid_backward(diff_data, kxE, kyE, kzE, kwE, temp, false);
      // cout << "Backward took " << T << endl;

      T.tic();

      // Residual
      complex<double> DhD(0.0, 0.0);
      complex<double> DhR(0., 0.);

#pragma omp parallel for
      for (int j = 0; j < diff_data.length(secondDim); j++) {
        complex<double> DhD_temp(0.0, 0.0);
        complex<double> DhR_temp(0.0, 0.0);

        for (int k = 0; k < diff_data.length(thirdDim); k++) {
          for (int i = 0; i < diff_data.length(firstDim); i++) {
            // Energy in D'*D
            DhD_temp += data.kdata(e, coil)(i, j, k) * conj(data.kdata(e, coil)(i, j, k));

            // Energy in D'*R
            DhR_temp += conj(data.kdata(e, coil)(i, j, k)) * diff_data(i, j, k);
          }
        }

// Prevent Race conditions in multi-threaded
#pragma omp critical
        {
          DhD += DhD_temp;
          DhR += DhR_temp;
        }
      }

      // Scale D*'D/(D'*R)
      complex<double> scale = DhD / DhR;
      complex<float> scaleF(real(scale), imag(scale));

#pragma omp parallel for
      for (int j = 0; j < diff_data.length(secondDim); j++) {
        double residual_temp = 0.0;
        for (int k = 0; k < diff_data.length(thirdDim); k++) {
          for (int i = 0; i < diff_data.length(firstDim); i++) {
            // Energy in D'*D
            residual_temp += norm(data.kdata(e, coil)(i, j, k) - scaleF * diff_data(i, j, k));
          }
        }

// Prevent Race conditions in multi-threaded
#pragma omp critical
        {
          residual += residual_temp;
        }
      }
      // cout << "Residual took " << T << endl;

    }  // Coils
  }  // Encode

  return (residual);
}

void RECON::export_slice(Array<complex<float>, 3> &temp, const char *fname) {
  int slice = (int)(temp.length(2) / 2);
  Array<complex<float>, 2> Xslice = temp(Range::all(), Range::all(), slice);
  ArrayWriteMagAppend(Xslice, fname);
}

void RECON::transform_in_encode(Array<Array<complex<float>, 3>, 2> &X, TransformDirection direction) {
  switch (cs_encode_transform) {
    case (DFT): {
      cout << "DFT in Encode" << endl;
      if (direction == FORWARD) {
        TRANSFORMS::fft_e(X);
      } else {
        TRANSFORMS::ifft_e(X);
      }
    } break;

    case (DIFF): {
      cout << "DIFF in Encode" << endl;
      if (direction == FORWARD) {
        TRANSFORMS::ediff(X);
      } else {
        TRANSFORMS::inv_ediff(X);
      }
    } break;

    case (PCA): {
      cout << "PCA in Encode" << endl;
      if (direction == FORWARD) {
        sparse_transform.eigen(X, 1, TRANSFORMS::FORWARD);
      } else {
        sparse_transform.eigen(X, 1, TRANSFORMS::BACKWARD);
      }
    } break;

    case (WAVELET): {
      cout << "WAVELET in Encode" << rcencodes << endl;
      if (direction == FORWARD) {
        TRANSFORMS::ewave(X);
      } else {
        TRANSFORMS::inv_ewave(X);
      }
    } break;

    default: {
    } break;
  }
}

void RECON::transform_in_time(Array<Array<complex<float>, 3>, 2> &X,
                              TransformDirection direction) {
  switch (cs_temporal_transform) {
    case (DFT): {
      cout << "DFT in Time" << endl;
      if (direction == FORWARD) {
        TRANSFORMS::fft_t(X);
      } else {
        TRANSFORMS::ifft_t(X);
      }
    } break;

    case (DIFF): {
      cout << "DIFF in Time" << endl;
      if (direction == FORWARD) {
        TRANSFORMS::tdiff(X);
      } else {
        TRANSFORMS::inv_tdiff(X);
      }
    } break;

    case (WAVELET): {
      cout << "WAVELET in Time" << endl;
      if (direction == FORWARD) {
        TRANSFORMS::twave(X);
      } else {
        TRANSFORMS::inv_twave(X);
      }
    } break;

    case (PCA): {
      cout << "PCA in Time" << endl;
      if (direction == FORWARD) {
        sparse_transform.eigen(X, 0, TRANSFORMS::FORWARD);
      } else {
        sparse_transform.eigen(X, 0, TRANSFORMS::BACKWARD);
      }
    } break;

    case (COMPOSITE_DIFF): {
      cout << "Composite Diff" << endl;
      for (Array<Array<complex<float>, 3>, 2>::iterator miter = X.begin(); miter != X.end(); miter++) {
        if (direction == FORWARD) {
          *miter -= composite_image;
        } else {
          *miter += composite_image;
        }
      }
    } break;

    default: {
    } break;
  }
}

void RECON::L1_threshold(Array<Array<complex<float>, 3>, 2> &X) {
  // Use a matrix to rotate into low resolution phase
  Array<complex<float>, 3> Phase(rcxres, rcyres, rczres, ColumnMajorArray<3>());
  if (phase_rotation) {
    cout << "Rotating to zero phase plane" << endl;
    Phase = complex<float>(0.0, 0.0);
    for (Array<Array<complex<float>, 3>, 2>::iterator miter = X.begin(); miter != X.end(); miter++) {
      Phase += (*miter);
    }

    // FFT to K-Space
    ifft(Phase);

    // Set Values to zero
    int cx = rcxres / 2;
    int cy = rcyres / 2;
    int cz = rczres / 2;

    // Calc windows
    Array<float, 1> Wx(Phase.length(firstDim));
    for (int i = 0; i < Phase.length(firstDim); i++) {
      float rx = abs(i - cx);
      Wx(i) = 1. / (1.0 + exp((rx - phase_rotation_sX) / 6));
    }

    Array<float, 1> Wy(Phase.length(secondDim));
    for (int i = 0; i < Phase.length(secondDim); i++) {
      float ry = abs(i - cy);
      Wy(i) = 1. / (1.0 + exp((ry - phase_rotation_sY) / 6));
    }

    Array<float, 1> Wz(Phase.length(thirdDim));
    for (int i = 0; i < Phase.length(thirdDim); i++) {
      float rz = abs(i - cz);
      Wz(i) = 1. / (1.0 + exp((rz - phase_rotation_sZ) / 6));
    }

    for (int k = 0; k < Phase.length(thirdDim); k++) {
      for (int j = 0; j < Phase.length(secondDim); j++) {
        for (int i = 0; i < Phase.length(firstDim); i++) {
          Phase(i, j, k) *= (Wx(i) * Wy(j) * Wz(k));
        }
      }
    }

    // FFT Back
    fft(Phase);

    // Normalize
    for (int k = 0; k < Phase.length(thirdDim); k++) {
      for (int j = 0; j < Phase.length(secondDim); j++) {
        for (int i = 0; i < Phase.length(firstDim); i++) {
          Phase(i, j, k) /= abs(Phase(i, j, k));
        }
      }
    }

    ArrayWritePhase(Phase, "PhaseCorr.dat");

    for (Array<Array<complex<float>, 3>, 2>::iterator miter = X.begin();
         miter != X.end(); miter++) {
      *miter *= conj(Phase);
    }
  }

  int actual_cycle_spins = cycle_spins;
  if (cs_spatial_transform != WAVELET) {
    actual_cycle_spins = 1;
  }

  transform_in_encode(X, FORWARD);
  transform_in_time(X, FORWARD);

  if (cs_spatial_transform == WAVELET) {
    cout << "WAVLET in Space" << endl;
    for (int cycle_spin = 0; cycle_spin < actual_cycle_spins; cycle_spin++) {
      cout << "    Spin " << cycle_spin << endl;
      wave.random_shift();
      for (Array<Array<complex<float>, 3>, 2>::iterator miter = X.begin();
           miter != X.end(); miter++) {
        wave.forward(*miter);
      }

      // Only update thresh once
      if (cycle_spin == 0) {
        softthresh.update_threshold(X, wave, 1. / (float)actual_cycle_spins);
      }
      softthresh.exec_threshold(X, wave);

      for (Array<Array<complex<float>, 3>, 2>::iterator miter = X.begin();
           miter != X.end(); miter++) {
        wave.backward(*miter);
      }
    }  // Cycle Spinning
  } else {
    softthresh.update_threshold(X, wave, 1.0);
    softthresh.exec_threshold(X, wave);

  }  // Wavelet

  transform_in_time(X, BACKWARD);
  transform_in_encode(X, BACKWARD);

  // Use a matrix to rotate into low resolution phase
  if (phase_rotation) {
    cout << "Rotating to non-zero phase plane" << endl;
    for (Array<Array<complex<float>, 3>, 2>::iterator miter = X.begin();
         miter != X.end(); miter++) {
      *miter *= Phase;
    }
  }
}

inline float sqr(float x) { return (x * x); }

void RECON::autofov(MRI_DATA &data, AutoFovMode automode, float autofov_thresh) {
  std::cout << "AUTOFOV :: Reconstruct images" << std::endl;
  std::cout << "AUTOFOV :: Native size " << data.xres << " x " << data.yres << " x " << data.zres << std::endl;

  HDF5 AutoFovDebug;
  AutoFovDebug = HDF5("DebugAutoFOV.h5", "w");

  gridFFT autofov_grid;
  autofov_grid.precalc_gridding(data.zres, data.yres, data.xres, data);

  Array<float, 3> sos_image;
  sos_image.setStorage(ColumnMajorArray<3>());
  sos_image.resize(autofov_grid.image.shape());
  sos_image = 0.0;

  Array<complex<float>, 3> complex_image;
  complex_image.setStorage(ColumnMajorArray<3>());
  complex_image.resize(autofov_grid.image.shape());

  for (int coil = 0; coil < data.Num_Coils; coil++) {
    complex_image = complex<float>(0.0, 0.0);
    for (int e = 0; e < 1; e++) {
      cout << "AUTOFOV :: Coil " << coil << "Encode " << e << endl;

      // Simple gridding
      complex_image = complex<float>(0.0, 0.0);
      autofov_grid.forward(complex_image, data.kdata(e, coil), data.kx(e), data.ky(e), data.kz(e), data.kw(e));
      sos_image += norm(complex_image);
    }
  }
  sos_image = sqrt(sos_image);
  cout << "AUTOFOV :: Shape" << sos_image.shape() << endl;

  // Blur in place
  AutoFovDebug.AddH5Array("Images", "SOS_PreBlur", sos_image);
  gaussian_blur(sos_image, 5.0, 5.0, 5.0);

  // Temp write to disk (remove once tested)
  AutoFovDebug.AddH5Array("Images", "SOS", sos_image);
  // ArrayWrite(sos_image, "AutoFov.dat");

  // Now we find the bounding box using a threshold
  float thresh = autofov_thresh * max(sos_image);
  int max_x = 0;
  int max_y = 0;
  int max_z = 0;
  int min_x = sos_image.length(firstDim);
  int min_y = sos_image.length(secondDim);
  int min_z = sos_image.length(thirdDim);
  for (int k = 0; k < sos_image.length(thirdDim); k++) {
    for (int j = 0; j < sos_image.length(secondDim); j++) {
      for (int i = 0; i < sos_image.length(firstDim); i++) {
        // Check for spherical mask
        if (automode == AUTOFOVSPHERE) {
          float rad = pow(2.0 * (float)i / (float)sos_image.length(firstDim) - 1.0, 2);
          rad += pow(2.0 * (float)j / (float)sos_image.length(secondDim) - 1.0, 2);
          rad += pow(2.0 * (float)k / (float)sos_image.length(thirdDim) - 1.0, 2);
          if (rad > 1.0) {
            continue;
          }
        } else if (automode == AUTOFOVCYLINDER) {
          float rad = pow(2.0 * (float)i / (float)sos_image.length(firstDim) - 1.0, 2);
          rad += pow(2.0 * (float)j / (float)sos_image.length(secondDim) - 1.0, 2);
          if (rad > 1.0) {
            continue;
          }
        }

        // Add the coordinates if greater than the threshold
        if (sos_image(i, j, k) > thresh) {
          max_x = max(i, max_x);
          max_y = max(j, max_y);
          max_z = max(k, max_z);
          min_x = min(i, min_x);
          min_y = min(j, min_y);
          min_z = min(k, min_z);
        }
      }
    }
  }

  int new_rcxres = 2 * max(sos_image.length(firstDim) / 2 - min_x, max_x - sos_image.length(firstDim) / 2);
  int new_rcyres = 2 * max(sos_image.length(secondDim) / 2 - min_y, max_y - sos_image.length(secondDim) / 2);
  int new_rczres = 2 * max(sos_image.length(thirdDim) / 2 - min_z, max_z - sos_image.length(thirdDim) / 2);

  // Make a multiple of 16 for FFT and block operators
  new_rcxres = 16 * (int)std::ceil((float)new_rcxres / 16.0);
  new_rcyres = 16 * (int)std::ceil((float)new_rcyres / 16.0);
  new_rczres = 16 * (int)std::ceil((float)new_rczres / 16.0);

  // Don't increase the resolution due to rounding
  new_rcxres = std::min(new_rcxres, data.xres);
  new_rcyres = std::min(new_rcyres, data.yres);
  new_rczres = std::min(new_rczres, data.zres);

  std::cout << "AUTOFOV :: X  [ " << min_x << " " << max_x << " ] -> recon to " << new_rcxres << std::endl;
  std::cout << "AUTOFOV :: Y  [ " << min_y << " " << max_y << " ] -> recon to " << new_rcyres << std::endl;
  std::cout << "AUTOFOV :: Z  [ " << min_z << " " << max_z << " ] -> recon to " << new_rczres << std::endl;

  float new_zoom_x = (float)sos_image.length(firstDim) / (float)new_rcxres;
  float new_zoom_y = (float)sos_image.length(secondDim) / (float)new_rcyres;
  float new_zoom_z = (float)sos_image.length(thirdDim) / (float)new_rczres;

  // Scale the FOV
  data.scale_fov(1. / new_zoom_x, 1. / new_zoom_y, 1. / new_zoom_z);

  // Set the resolution
  data.recon_res(0) = new_rcxres;
  data.recon_res(1) = new_rcyres;
  data.recon_res(2) = new_rczres;
}

void RECON::calc_sensitivity_maps(int argc, const char **argv, MRI_DATA &data) {
  // ------------------------------------
  //  Get coil sensitivity map ( move into function)
  // ------------------------------------

  HDF5 SmapsDebug;
  if (debug_smaps) {
    cout << "Exporting Smaps debuging to DebugSenseMaps.h5" << endl;
    SmapsDebug = HDF5("DebugSenseMaps.h5", "w");
  }
  switch (recon_type) {
    case (CLEAR): {
      if (intensity_correction) {
        // Clear doesn't need anything, but might need to handle sensitivity
        // correction
        for (int e = 0; e < 1; e++) {
          for (int coil = 0; coil < data.Num_Coils; coil++) {
            cout << "Coil = " << coil << " encode = " << e << endl;

            // gridding.forward(data.kdata(e,coil),data.kx(e),data.ky(e),data.kz(e),data.kw(e));
            IntensityCorrection += norm(gridding.image);
          }
        }
        IntensityCorrection = sqrt(IntensityCorrection);

        // Gaussian blur
        Array<float, 3> Blur;
        Blur.setStorage(ColumnMajorArray<3>());
        Blur.resize(rcxres, rcyres, rczres);
        normalized_gaussian_blur(IntensityCorrection, Blur,
                                 intensity_correction_blurX, intensity_correction_blurY, intensity_correction_blurZ);

        IntensityCorrection = Blur;
      }

    } break;

    case (ADMM):
    case (PILS):
    case (IST):
    case (FISTA):
    case (CG): {
      // Allocate Storage for Map	and zero
      cout << "Allocate Sense Maps" << endl
           << flush;
      {
        Array<Array<complex<float>, 3>, 1> temp = Alloc4DContainer<complex<float> >(rcxres, rcyres, rczres, data.Num_Coils);
        smaps.reference(temp);
      }

      // Break out of this loop if there is only one coil
      if (data.Num_Coils == 1) {
        cout << "One coil - assuming senitivity map is 1" << endl
             << flush;
        smaps(0) = complex<float>(1.0, 0.0);
        break;
      }

      // ------------------------------------
      cout << "Recon Low Resolution Images" << endl
           << flush;

      // Multiple Encode Smaps
      Array<Array<complex<float>, 3>, 2> image_store;
      int smap_num_encodes = 1;
      if (smap_use_all_encodes) {
        cout << "Using all encodes for coil maps" << endl;
        smap_num_encodes = data.Num_Encodings - this->smap_skip_encode.length(firstDim);
        Array<Array<complex<float>, 3>, 2> temp = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, smap_num_encodes, data.Num_Coils);
        image_store.reference(temp);
      } else if (coil_combine_type == WALSH) {
        cout << "Using one encodes for coil maps" << endl;
        Array<Array<complex<float>, 3>, 2> temp = Alloc5DContainer<complex<float> >(rcxres, rcyres, rczres, 1, data.Num_Coils);
        image_store.reference(temp);
      } else {
        // Point to Smaps
        image_store.setStorage(ColumnMajorArray<2>());
        image_store.resize(1, data.Num_Coils);
        for (int coil = 0; coil < data.Num_Coils; coil++) {
          image_store(0, coil).reference(smaps(coil));
        }
      }

      // Low Pass filtering for Sensitivity Map
      if (coil_combine_type != ESPIRIT) {
        gridding.k_rad = smap_res;
      }

      if (smap_use_all_encodes) {
        int smap_encode_count = 0;
        for (int e = 0; e < data.Num_Encodings; e++) {
          // Skip if in array
          bool skip_this_encode = false;
          for (int i = 0; i < (int)this->smap_skip_encode.numElements(); i++) {
            skip_this_encode |= (e == this->smap_skip_encode(i));
          }

          if (skip_this_encode) {
            std::cout << "Skipping Encode" << e << "for coil sensitivity mapping" << std::endl;
            continue;
          }

          for (int coil = 0; coil < data.Num_Coils; coil++) {
            cout << "Coil " << coil << "Encode " << e << endl;
            // Simple gridding
            gridding.forward(image_store(smap_encode_count, coil),
                             data.kdata(e, coil), data.kx(e), data.ky(e),
                             data.kz(e), data.kw(e));

            // Gaussian blur
            gaussian_blur(image_store(smap_encode_count, coil),
                          extra_blurX,
                          extra_blurY,
                          extra_blurZ);  // TEMP
          }
          smap_encode_count++;
        }
      } else {
        for (int coil = 0; coil < data.Num_Coils; coil++) {
          image_store(0, coil) = complex<float>(0.0, 0.0);

          int used_encodes = smap_nex_encodes ? data.Num_Encodings : 1;
          for (int e = 0; e < used_encodes; e++) {
            cout << "Coil " << coil << "Encode " << e << endl;
            // Simple gridding
            gridding.forward(image_store(0, coil), data.kdata(e, coil),
                             data.kx(e), data.ky(e), data.kz(e), data.kw(e));
          }
          // Gaussian blur
          gaussian_blur(image_store(0, coil), extra_blurX, extra_blurY,
                        extra_blurZ);  // TEMP
        }
      }
      gridding.k_rad = 9999;

      if (debug_smaps) {
        for (int encode = 0; encode < image_store.length(firstDim); encode++) {
          for (int coil = 0; coil < data.Num_Coils; coil++) {
            char name[512];
            sprintf(name, "Encode_%03d_Coil_%03d", encode, coil);
            SmapsDebug.AddH5Array("Images", name, image_store(encode, coil));

            Array<float, 3> Imag(image_store(encode, coil).shape(), ColumnMajorArray<3>());
            Imag = abs(image_store(encode, coil));
            SmapsDebug.AddH5Array("MagImages", name, Imag);
          }
        }
      }

      if (intensity_correction) {
        intensity_correct(this->IntensityCorrection, image_store);
      }

      if (this->coil_rejection_flag) {
        float cx = smaps(0).length(firstDim) / 2;
        float cy = smaps(0).length(secondDim) / 2;
        float cz = smaps(0).length(thirdDim) / 2;

        Array<double, 1> SignalFractionIn(data.Num_Coils);
        Array<double, 1> SignalFractionTotal(data.Num_Coils);

        SignalFractionIn = 0.0;
        SignalFractionTotal = 0.0;

        // Get the coils contribution to the image in the center
        switch (this->coil_rejection_shape) {
          default:
          case (0): {
            for (int coil = 0; coil < data.Num_Coils; coil++) {
              for (int e = 0; e < image_store.length(firstDim); e++) {
                for (int k = 0; k < smaps(0).length(thirdDim); k++) {
                  for (int j = 0; j < smaps(0).length(secondDim); j++) {
                    for (int i = 0; i < smaps(0).length(firstDim); i++) {
                      float x = ((float)i - cx) / (2.0 * cx);
                      float y = ((float)j - cy) / (2.0 * cy);
                      float z = ((float)k - cz) / (2.0 * cz);

                      // Check the radius
                      if ((x * x + y * y + z * z) < (this->coil_rejection_radius * this->coil_rejection_radius)) {
                        SignalFractionIn(coil) += abs(image_store(e, coil)(i, j, k));
                      }
                      SignalFractionTotal(coil) += abs(image_store(e, coil)(i, j, k));
                    }  // i
                  }  // j
                }  // k
              }  // e
            }  // coil

          } break;

          case (1): {
            for (int coil = 0; coil < data.Num_Coils; coil++) {
              for (int e = 0; e < image_store.length(firstDim); e++) {
                for (int k = 0; k < smaps(0).length(thirdDim); k++) {
                  for (int j = 0; j < smaps(0).length(secondDim); j++) {
                    for (int i = 0; i < smaps(0).length(firstDim); i++) {
                      // Check the radius
                      float z = ((float)k - cz) / (2.0 * cz);
                      float y = ((float)j - cy) / (2.0 * cy);
                      float x = ((float)i - cx) / (2.0 * cx);

                      // Check the radius
                      if (abs(z) < this->coil_rejection_height) {
                        if ((y * y + x * x) < this->coil_rejection_radius * this->coil_rejection_radius) {
                          SignalFractionIn(coil) += abs(image_store(e, coil)(i, j, k));
                        }
                      }
                      SignalFractionTotal(coil) += abs(image_store(e, coil)(i, j, k));

                    }  // i
                  }  // j
                }  // k
              }  // e
            }  // coil

          } break;

        }  // shape switch

        // Normalize
        for (int coil = 0; coil < data.Num_Coils; coil++) {
          cout << "Coil = " << coil << "Signal In=" << SignalFractionIn(coil) << " Signal Out=" << SignalFractionTotal(coil);
          SignalFractionIn(coil) /= SignalFractionTotal(coil);
          cout << " Fraction=" << SignalFractionIn(coil) << endl;
        }

        for (int coil = 0; coil < data.Num_Coils; coil++) {
          std::cout << "Coil = " << coil << "Signal Fraction =" << SignalFractionIn(coil) << std::endl;
          if (SignalFractionIn(coil) < this->coil_rejection_thresh) {
            std::cout << "    REJECTING COIL " << coil << std::endl;
            for (int e = 0; e < image_store.length(firstDim); e++) {
              image_store(e, coil) = complex<float>(0.0, 0.0);
            }
          }
        }  // coil
      }  // coil rejection

      // Spirit Code
      switch (coil_combine_type) {
        case (ESPIRIT): {
          // Espirit (k-space kernal Eigen method)
          cout << "eSPIRIT Based Maps" << endl;
          SPIRIT S;
          S.read_commandline(argc, argv);
          S.init(rcxres, rcyres, rczres, data.Num_Coils);
          S.generateEigenCoils(smaps, image_store);
        } break;

        case (WALSH): {
          // Image space eigen Method
          eigen_coils(smaps, image_store);

          switch (smap_norm) {
            case (SMAPNORM_SOS): {
              sos_normalize(smaps);
            } break;

            case (SMAPNORM_BODY): {
              body_coil_normalize(smaps, body_coil_idx);
            } break;
          }

        } break;

        case (LOWRES): {  // E-spirit Code

          // Sos Normalization
          cout << "Normalize Coils" << endl;
#pragma omp parallel for
          for (int k = 0; k < smaps(0).length(thirdDim); k++) {
            for (int j = 0; j < smaps(0).length(secondDim); j++) {
              for (int i = 0; i < smaps(0).length(firstDim); i++) {
                for (int coil = 0; coil < data.Num_Coils; coil++) {
                  complex<float> s(0.0, 0.0);
                  for (int e = 0; e < smap_num_encodes; e++) {
                    s += image_store(e, coil)(i, j, k);
                  }
                  smaps(coil)(i, j, k) = s;
                }
              }
            }
          }

          switch (smap_norm) {
            case (SMAPNORM_SOS): {
              sos_normalize(smaps);
            } break;

            case (SMAPNORM_BODY): {
              body_coil_normalize(smaps, body_coil_idx);
            } break;
          }

        } break;
      }  // Normalization

      // Export
      if (export_smaps || debug_smaps) {
        HDF5 SmapsOut = HDF5("SenseMaps.h5", "w");
        for (int coil = 0; coil < smaps.length(firstDim); coil++) {
          char name[256];
          sprintf(name, "SenseMaps_%d", coil);
          SmapsOut.AddH5Array("Maps", name, smaps(coil));
        }
      }

    } break;

    case (SOS): {
      // Allocate Storage for Map	and zero
      cout << "Allocate Sense Maps (SOS)" << endl
           << flush;
      {
        Array<Array<complex<float>, 3>, 1> temp = Alloc4DContainer<complex<float> >(rcxres, rcyres, rczres, data.Num_Coils);
        smaps.reference(temp);
      }

      for (int coil = 0; coil < smaps.length(firstDim); coil++) {
        smaps(coil) = complex<float>(1.0, 0.0);
      }

    } break;
  }

  if (smap_mask != SMAPMASK_NONE) {
    for (int coil = 0; coil < data.Num_Coils; coil++) {
#pragma omp parallel for
      for (int k = 0; k < smaps(0).length(thirdDim); k++) {
        for (int j = 0; j < smaps(0).length(secondDim); j++) {
          for (int i = 0; i < smaps(0).length(firstDim); i++) {
            float r = 0.0;
            switch (smap_mask) {
              case (SMAPMASK_CIRCLE): {
                float x = 2.0 * (i - 0.5 * ((float)smaps(0).length(firstDim))) / ((float)smaps(0).length(firstDim));
                float y = 2.0 * (j - 0.5 * ((float)smaps(0).length(secondDim))) / ((float)smaps(0).length(secondDim));
                r = sqrt(x * x + y * y);

              } break;

              case (SMAPMASK_SPHERE): {
                float x = 2.0 * (i - 0.5 * ((float)smaps(0).length(firstDim))) / ((float)smaps(0).length(firstDim));
                float y = 2.0 * (j - 0.5 * ((float)smaps(0).length(secondDim))) / ((float)smaps(0).length(secondDim));
                float z = 2.0 * (k - 0.5 * ((float)smaps(0).length(thirdDim))) / ((float)smaps(0).length(thirdDim));
                r = sqrt(x * x + y * y + z * z);
              } break;

              default: {
                r = 0.0;
              } break;
            }

            if (r > 1) {
              smaps(coil)(i, j, k) = 0.0;
            }
          }
        }
      }
    }
  }
}

void RECON::body_coil_normalize(Array<Array<complex<float>, 3>, 1> &A, int body_coil_idx) {
  if (body_coil_idx == -1) {
    cout << "Using last coil for body coil normalization" << endl;
    body_coil_idx = A.length(firstDim) - 1;
  }

#pragma omp parallel for
  for (int k = 0; k < A(0).length(thirdDim); k++) {
    for (int j = 0; j < A(0).length(secondDim); j++) {
      for (int i = 0; i < A(0).length(firstDim); i++) {
        // Signal from the body coil
        complex<float> body_coil_phase = polar(1.0f, -arg(A(body_coil_idx)(i, j, k)));

        // The sum of squares
        float sos = 0.0;
        for (int coil = 0; coil < A.length(firstDim); coil++) {
          sos += norm(A(coil)(i, j, k));
        }
        sos = sqrt(sos);

        // Normalize by SOS but use phase from body coil
        for (int coil = 0; coil < A.length(firstDim); coil++) {
          A(coil)
          (i, j, k) *= body_coil_phase / sos;
        }
      }
    }
  }  // x,y,z
}

void RECON::sos_normalize(Array<Array<complex<float>, 3>, 1> &A) {
  if (smap_thresh == 0.0) {
#pragma omp parallel for
    for (int k = 0; k < A(0).length(thirdDim); k++) {
      for (int j = 0; j < A(0).length(secondDim); j++) {
        for (int i = 0; i < A(0).length(firstDim); i++) {
          // The sum of squares
          float sos = 0.0;
          for (int coil = 0; coil < A.length(firstDim); coil++) {
            sos += norm(A(coil)(i, j, k));
          }
          sos = sqrt(sos);

          // Normalize
          for (int coil = 0; coil < A.length(firstDim); coil++) {
            A(coil)
            (i, j, k) /= sos;
          }
        }
      }
    }  // x,y,z

  } else {
    // Get the threshold
    Array<float, 3> IC(rcxres, rcyres, rczres, ColumnMajorArray<3>());
    IC = 0.0;
    for (int coil = 0; coil < A.length(firstDim); coil++) {
#pragma omp parallel for
      for (int k = 0; k < A(0).length(thirdDim); k++) {
        for (int j = 0; j < A(0).length(secondDim); j++) {
          for (int i = 0; i < A(0).length(firstDim); i++) {
            IC(i, j, k) += norm(A(coil)(i, j, k));
          }
        }
      }
    }
    IC = sqrt(IC);
    float max_IC = max(IC);
    cout << "Max IC = " << max_IC << endl;
    cout << "Thresh of " << max_IC * smap_thresh << endl;

    // Create a binary mask
    Array<int, 3> MASK(rcxres, rcyres, rczres, ColumnMajorArray<3>());
    MASK = 0;
    for (int k = 0; k < A(0).length(thirdDim); k++) {
      for (int j = 0; j < A(0).length(secondDim); j++) {
        for (int i = 0; i < A(0).length(firstDim); i++) {
          if (IC(i, j, k) < smap_thresh * max_IC) {
            MASK(i, j, k) = 0.0;
          } else {
            MASK(i, j, k) = 1.0;
          }
        }
      }
    }
    // ArrayWrite(MASK,"PreFill.dat");

    // Try to remove the holes
    for (int k = 0; k < MASK.length(thirdDim); k++) {
      for (int i = 0; i < MASK.length(firstDim); i++) {
        int leading_edge = -1;
        for (int j = 0; j < MASK.length(secondDim); j++) {
          if (MASK(i, j, k) > 0) {
            leading_edge = j;
            break;
          }
        }

        int trailing_edge = -1;
        for (int j = MASK.length(secondDim) - 1; j >= 0; j--) {
          if (MASK(i, j, k) > 0) {
            trailing_edge = j;
            break;
          }
        }

        if (leading_edge > -1) {
          for (int j = leading_edge; j <= trailing_edge; j++) {
            MASK(i, j, k) = 1;
          }
        }
      }
    }
    // ArrayWrite(MASK,"PostFill.dat");

    // Now normalize
    for (int coil = 0; coil < A.length(firstDim); coil++) {
#pragma omp parallel for
      for (int k = 0; k < A(0).length(thirdDim); k++) {
        for (int j = 0; j < A(0).length(secondDim); j++) {
          for (int i = 0; i < A(0).length(firstDim); i++) {
            if (MASK(i, j, k) > 0) {
              A(coil)
              (i, j, k) /= IC(i, j, k);
            } else {
              A(coil)
              (i, j, k) = complex<float>(0.0, 0.0);
            }
          }
        }
      }
    }
  }
}

void RECON::apply_intensity_bias(Array<Array<complex<float>, 3>, 2> &images) {
  if (intensity_correction) {
    for (Array<Array<complex<float>, 3>, 2>::iterator miter = images.begin(); miter != images.end(); miter++) {
      *miter /= IntensityCorrection;
    }
  }
}

void RECON::apply_intensity_bias(Array<Array<complex<float>, 3>, 1> &images) {
  if (intensity_correction) {
    for (Array<Array<complex<float>, 3>, 1>::iterator miter = images.begin(); miter != images.end(); miter++) {
      *miter /= IntensityCorrection;
    }
  }
}

void RECON::apply_intensity_bias(Array<complex<float>, 3> &images) {
  if (intensity_correction) {
    images /= IntensityCorrection;
  }
}

void RECON::correct_intensity_bias(Array<Array<complex<float>, 3>, 2> &images) {
  if (intensity_correction) {
    for (Array<Array<complex<float>, 3>, 2>::iterator miter = images.begin(); miter != images.end(); miter++) {
      *miter *= IntensityCorrection;
    }
  }
}

void RECON::correct_intensity_bias(Array<Array<complex<float>, 3>, 1> &images) {
  if (intensity_correction) {
    for (Array<Array<complex<float>, 3>, 1>::iterator miter = images.begin(); miter != images.end(); miter++) {
      *miter *= IntensityCorrection;
    }
  }
}

void RECON::correct_intensity_bias(Array<complex<float>, 3> &images) {
  if (intensity_correction) {
    images *= IntensityCorrection;
  }
}

void RECON::intensity_correct(Array<float, 3> &IC, Array<Array<complex<float>, 3>, 2> &smaps) {
  cout << "Correcting Intensity" << endl;

  // Allocate Storage for correction
  IC.free();
  IC.setStorage(ColumnMajorArray<3>());
  IC.resize(rcxres, rcyres, rczres);
  IC = 0.0;

  // Collect a Sum of Squares image
  for (int encode = 0; encode < smaps.length(firstDim); encode++) {
    for (int coil = 0; coil < smaps.length(secondDim); coil++) {
      IC += norm(smaps(encode, coil));
    }
  }

  // Take a square root
  for (Array<float, 3>::iterator miter = IC.begin(); miter != IC.end(); miter++) {
    (*miter) = sqrtf(*miter);
  }

  // In this case a body coil image is in the volume
  float zero_protection = 0.0001;
  if (smap_norm == SMAPNORM_BODY) {
    if (body_coil_idx == -1) {
      cout << "Using last coil for body coil normalization" << endl;
      body_coil_idx = smaps.length(secondDim) - 1;
    }

    // sos / body coil
    IC /= abs(smaps(0, body_coil_idx));

    zero_protection = 1e-6;

  } else {
    IC /= max(IC);

    // Threshold that image - to reduce error with air
    float max_sos = max(IC);
    float ic_threshold = 0.05 * max_sos;
    for (Array<float, 3>::iterator miter = IC.begin(); miter != IC.end(); miter++) {
      if ((*miter) < ic_threshold) {
        (*miter) = 0.0;
      }
    }
  }

  // Gaussian blur
  Array<float, 3> Blur;
  Blur.setStorage(ColumnMajorArray<3>());
  Blur.resize(rcxres, rcyres, rczres);
  normalized_gaussian_blur(IC, Blur, intensity_correction_blurX, intensity_correction_blurY, intensity_correction_blurZ);

  // Calc intensity correction
  float max_blur = max(Blur);
  for (int k = 0; k < rczres; k++) {
    for (int j = 0; j < rcyres; j++) {
      for (int i = 0; i < rcxres; i++) {
        IC(i, j, k) = Blur(i, j, k) / (Blur(i, j, k) * Blur(i, j, k) + zero_protection * max_blur * max_blur);
      }
    }
  }

  ArrayWrite(IC, "IntensityCorrection.dat");

  return;
}

void RECON::normalized_gaussian_blur(const Array<float, 3> &In, Array<float, 3> &Out,
                                     float sigmax, float sigmay, float sigmaz) {
  int dwinX = 3 * (int)sigmax;
  int dwinY = 3 * (int)sigmay;
  int dwinZ = 3 * (int)sigmaz;

  Array<float, 3> Normalization;
  Normalization.setStorage(ColumnMajorArray<3>());
  Normalization.resize(rcxres, rcyres, rczres);
  Normalization = 0.0;

  // Set to zero
  Out = 0.0;

  // Kernel to reduce calls to exp
  Array<float, 1> kernx(2 * dwinX + 1);
  for (int t = 0; t < (2 * dwinX + 1); t++) {
    kernx(t) = exp(-sqr((float)t - dwinX) / (2.0 * sqr(sigmax)));
  }

  // Gaussian Blur in X
  cout << "Blur in X" << endl;
#pragma omp parallel for
  for (int k = 0; k < rczres; k++) {
    for (int j = 0; j < rcyres; j++) {
      for (int i = 0; i < rcxres; i++) {
        int ks = 0;
        int sx = i - dwinX;
        if (sx < 0) {
          ks = -sx;
          sx = 0;
        }
        int ex = min(i + dwinX, rcxres);

        for (int ii = sx; ii < ex; ii++) {
          Out(i, j, k) += kernx(ks) * In(ii, j, k);
          Normalization(i, j, k) += kernx(ks) * (float)(In(ii, j, k) > 0);
          ks++;
        }
      }
    }
  }

  // Kernel to reduce calls to exp
  Array<float, 1> kerny(2 * dwinY + 1);
  for (int t = 0; t < (2 * dwinY + 1); t++) {
    kerny(t) = exp(-sqr((float)t - dwinY) / (2.0 * sqr(sigmay)));
  }

  // Gaussian Blur in Y
  cout << "Blur in Y" << endl;
#pragma omp parallel for
  for (int k = 0; k < rczres; k++) {
    for (int i = 0; i < rcxres; i++) {
      // Copy
      float *Line = new float[rcyres];
      float *NLine = new float[rcyres];
      for (int j = 0; j < rcyres; j++) {
        Line[j] = Out(i, j, k);
        NLine[j] = Normalization(i, j, k);
      }

      // Convolve
      for (int j = 0; j < rcyres; j++) {
        int ks = 0;
        int sy = j - dwinY;
        if (sy < 0) {
          ks = -sy;  // Offset;
          sy = 0;
        }
        int ey = min(j + dwinY + 1, rcyres);

        for (int jj = sy; jj < ey; jj++) {
          Out(i, j, k) += kerny(ks) * Line[jj];
          Normalization(i, j, k) += kerny(ks) * NLine[jj];  // Already a count
          ks++;
        }
      }

      delete[] Line;
      delete[] NLine;
    }
  }

  // Kernel to reduce calls to exp
  Array<float, 1> kernz(2 * dwinZ + 1);
  for (int t = 0; t < (2 * dwinZ + 1); t++) {
    kernz(t) = exp(-sqr((float)t - dwinZ) / (2.0 * sqr(sigmaz)));
  }

  // Gaussian Blur in Z
  cout << "Blur in Z" << endl;
#pragma omp parallel for
  for (int j = 0; j < rcyres; j++) {
    for (int i = 0; i < rcxres; i++) {
      // Copy
      float *Line = new float[rczres];
      float *NLine = new float[rczres];
      for (int k = 0; k < rczres; k++) {
        Line[k] = Out(i, j, k);
        NLine[k] = Normalization(i, j, k);
      }

      // Convolve
      for (int k = 0; k < rczres; k++) {
        int ks = 0;
        int sz = k - dwinZ;
        if (sz < 0) {
          ks = -sz;  // Offset;
          sz = 0;
        }
        int ez = min(k + dwinZ + 1, rczres);

        for (int kk = sz; kk < ez; kk++) {
          Out(i, j, k) += kernz(ks) * Line[kk];
          Normalization(i, j, k) += kernz(ks) * NLine[kk];
          ks++;
        }

        if (Normalization(i, j, k) > 0) {
          Out(i, j, k) /= Normalization(i, j, k);
        }
      }

      delete[] Line;
      delete[] NLine;
    }
  }
}

/**
 * Generate coils from SPIRiT kernel.
 */
void RECON::eigen_coils(Array<Array<complex<float>, 3>, 1> &smaps,
                        Array<Array<complex<float>, 3>, 2> &image) {
  // Shorthand
  int Nencodes = image.extent(firstDim);
  int Ncoils = image.extent(secondDim);
  int Nx = image(0).extent(firstDim);
  int Ny = image(0).extent(secondDim);
  int Nz = image(0).extent(thirdDim);

  // Get coil with max signal
  cout << "Getting maximum coil" << endl;
  int ref_coil = 0;
  float max_sig = 0;
  for (int coil = 0; coil < Ncoils; coil++) {
    float sig = 0;
    for (int e = 0; e < Nencodes; e++) {
      sig += sum(norm(image(e, coil)));
    }

    if (sig > max_sig) {
      ref_coil = coil;
      max_sig = sig;
    }
  }
  cout << "Reference coil = " << ref_coil << endl;

  int block_size_x = walsh_block_sizeX;
  int block_size_y = walsh_block_sizeY;
  int block_size_z = walsh_block_sizeZ;

  // Blocks shouldn't be larger than the dimension
  block_size_x = (block_size_x > Nx) ? (Nx) : (block_size_x);
  block_size_y = (block_size_y > Ny) ? (Ny) : (block_size_y);
  block_size_z = (block_size_z > Nz) ? (Nz) : (block_size_z);

  int block_hsize_x = block_size_x / 2;
  int block_hsize_y = block_size_y / 2;
  int block_hsize_z = block_size_z / 2;

  // Get the block stride
  int block_stride_x = 1;
  int block_stride_y = 1;
  int block_stride_z = 1;

  // Size of eigen block
  int Np = block_size_x * block_size_y * block_size_z * Nencodes;

  cout << "Eigen Coil ( " << Nx << " x " << Ny << " x " << Nz << " x " << Ncoils << " x " << Nencodes << endl;
  cout << "Actual Block Size = " << block_size_x << " x " << block_size_y << " x " << block_size_z << endl;
  cout << "Getting Low Rank threshold (N=" << Ncoils << ")(Np = " << Np << ")" << endl;
  cout << "Block stride = ( " << block_stride_x << "," << block_stride_y << "," << block_stride_z << ")" << endl;

  int block_Nx = (int)(Nx / block_stride_x);
  int block_Ny = (int)(Ny / block_stride_y);
  int block_Nz = (int)(Nz / block_stride_z);

  // Center the block
  int block_offset_x = block_stride_x / 2;
  int block_offset_y = block_stride_y / 2;
  int block_offset_z = block_stride_z / 2;

  int total_blocks = block_Nx * block_Ny * block_Nz;
  cout << "Total Block Number" << total_blocks << " ( " << block_Nx << "," << block_Ny << "," << block_Nz << ")" << endl;

  // This is for nested workaround in old versions of openmp
  int *N = new int[3];
  N[0] = block_Nx;
  N[1] = block_Ny;
  N[2] = block_Nz;

  tictoc T;
  T.tic();

#pragma omp parallel for
  for (int block = 0; block < total_blocks; block++) {
    // cout << "Block " << block << " of " << total_blocks << endl;

    // Get the actual position
    int *ID = new int[3];
    nested_workaround(block, N, ID, 3);
    int i = ID[0] * block_stride_x + block_offset_x;
    int j = ID[1] * block_stride_y + block_offset_y;
    int k = ID[2] * block_stride_z + block_offset_z;
    delete[] ID;

    //-----------------------------------------------------
    //   Gather Block coordinates
    //-----------------------------------------------------
    int kstart = k - block_hsize_z;
    int kstop = k + block_hsize_z;

    if (kstart < 0) {
      kstart = 0;
      kstop = block_size_z;
    }

    if (kstop > Nz) {
      kstop = Nz;
      kstart = Nz - block_size_z;
    }

    int jstart = j - block_hsize_y;
    int jstop = j + block_hsize_y;
    if (jstart < 0) {
      jstart = 0;
      jstop = block_size_y;
    }

    if (jstop > Ny) {
      jstop = Ny;
      jstart = Ny - block_size_y;
    }

    int istart = i - block_hsize_x;
    int istop = i + block_hsize_x;
    if (istart < 0) {
      istart = 0;
      istop = block_size_x;
    }

    if (istop > Nx) {
      istop = Nx;
      istart = Nx - block_size_x;
    }

    //-----------------------------------------------------
    //   Scatter Block coordinates
    //----------------------------------------------------
    int scatter_istart = i;
    int scatter_istop = i + block_stride_x;

    int scatter_jstart = j;
    int scatter_jstop = j + block_stride_y;

    int scatter_kstart = k;
    int scatter_kstop = k + block_stride_z;

    //-----------------------------------------------------
    //   Collect a Block
    //-----------------------------------------------------

    arma::cx_fmat R;
    R.zeros(Ncoils, Np);
    for (int c = 0; c < Ncoils; c++) {
      int count = 0;
      for (int e = 0; e < Nencodes; e++) {
        for (int kk = kstart; kk < kstop; kk++) {
          for (int jj = jstart; jj < jstop; jj++) {
            for (int ii = istart; ii < istop; ii++) {
              R(c, count) = image(e, c)(ii, jj, kk);
              count++;
            }
          }
        }
      }
    }

    // SVD to Get Eigen Vector
    arma::cx_fmat U;
    arma::cx_fmat V;
    arma::fvec s;
    arma::svd_econ(U, s, V, R, "left");

    // Scatter the sensitivity
    for (int c = 0; c < Ncoils; c++) {
      // Sensitivty for the map
      complex<float> temp = sqrt(s(0)) * U(c, 0) * polar<float>(1.0, -arg(U(ref_coil, 0)));

      // Scatter to the block
      for (int kk = scatter_kstart; kk < scatter_kstop; kk++) {
        for (int jj = scatter_jstart; jj < scatter_jstop; jj++) {
          for (int ii = scatter_istart; ii < scatter_istop; ii++) {
            smaps(c)(ii, jj, kk) = temp;
          }
        }
      }
    }

  }  // Block (threaded)
  cout << "(eigen coil took " << T << ")" << endl;
}
