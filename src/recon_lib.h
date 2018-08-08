#pragma once 

// Standard Libraries
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstring>
#include <string>
#include <complex>
#include <omp.h>

// External Libraries
#include <armadillo>

// Local Libraries
#include "ArrayTemplates.hpp"
#include "wavelet3D.h"	
#include "temporal_diff.h"	
#include "gridFFT.h"
#include "smsEncode.h"
#include "DCFgridFFT.h"
#include "mri_data.h"
#include "threshold.h"
#include "spirit.h"
#include "phantom.h"
#include "gating.h"
#include "tictoc.hpp"
#include "gating.h"
#include "clear.h"
#include "l2reg.h"
#include "hdf5_interface.h"
#include "polynomial_fitting.h"
#include "voronoi_dcf.h"

class RECON{
	public:
	  // Types of Recons
	  enum ReconType { SOS, PILS, CG, IST, FISTA, CLEAR, ADMM};
	   
	  // Data Types
	  enum DataType { EXTERNAL, SIMULATE, PSF, PHANTOM, BENCHMARK};

	  // Coil Combine Type
	  enum CoilCombineType { LOWRES, ESPIRIT, WALSH };
	  
	  // Enum Transform Types
	  enum TransformType {NONE,WAVELET,DIFF,DFT,PCA,COMPOSITE_DIFF};
	  enum TransformDirection {FORWARD,BACKWARD};
	  
	  // Enum Sensitivity Maps
	  enum SmapMaskType {SMAPMASK_NONE,SMAPMASK_CIRCLE,SMAPMASK_SPHERE};
	  SmapMaskType smap_mask;
	 	  	 
	  // Recon Flags  
	  ReconType recon_type;
	  DataType data_type;
	  CoilCombineType coil_combine_type;
	  
	  // CS Flags
	  TransformType cs_spatial_transform;
	  TransformType cs_temporal_transform;
	  TransformType cs_encode_transform;
	  
	  // Reconstruction modules
	  WAVELET3D wave;
	  THRESHOLD softthresh;
	  LOWRANKCOIL lrankcoil;
	  LOWRANKCOIL lranktime;
	  L2REG l2reg;
	  TRANSFORMS sparse_transform;  
	  	  
	  bool complex_diff;
	  
	  float zoom;
	  float zoom_x;
	  float zoom_y;
	  float zoom_z;
	  
	  int threads;
	  
	  // Reconstructed Images
	  int rcxres;
	  int rcyres;
	  int rczres;
	  int rcframes;
	  int rcencodes;
  	  
	  int rc_frame_start;
	  int step_update_frequency;
	  
	  // Code to rotate into low resolution images
	  bool phase_rotation;
  	  int phase_rotation_sX;
	  int phase_rotation_sY;
	  int phase_rotation_sZ;
	  		
  	  // Gridding operators
	  gridFFT gridding;
	  	  
	  // Arrays for calculations
	  NDarray::Array< NDarray::Array< complex<float>,3>,1 >smaps; // Array of arrays
	  NDarray::Array< complex<float>,3>composite_image;
	  NDarray::Array< float, 3 >IntensityCorrection;
	  GATING gate;
	  	  
	  // Density calcs
	  enum DcfType { SUPPLIED, RECALC_VOR, RECALC_DCF};
	  DcfType dcf_type;
	  int dcf_iter;
	  float dcf_dwin;
	  float dcf_scale;
	  float dcf_overgrid;
	  float dcf_acc;
	  
	  int acc;
	  float compress_coils;
	  bool whiten;
      bool coil_rejection_flag;
      float coil_rejection_radius;
      int coil_rejection_shape;

	  float smap_res;
	  float smap_thresh;
	  bool intensity_correction;
	  char filename[1024];
	  int cycle_spins;
	  bool iterative_smaps;
	  int walsh_block_sizeX;
	  int walsh_block_sizeY;
	  int walsh_block_sizeZ;
	  float extra_blurX;
	  float extra_blurY;
	  float extra_blurZ;
	  
	  float blurX;
	  float blurY;
	  float blurZ;
	  
	  bool reset_dens;
	  
	  
	  bool smap_use_all_encodes;
	  bool smap_nex_encodes;
	  
	  int wavelet_levelsX; 
	  int wavelet_levelsY; 
	  int wavelet_levelsZ; 
	  
	  float admm_alpha;
	  float admm_rho;
	  float admm_gamma;
	  int admm_max_iter;	  
		  
	  int max_iter;
	  int export_smaps;
	  bool prep_done;
	  bool pregate_data_flag;
	  
	  NDarray::Array< NDarray::Array< complex<float>,3>, 2 >  test_sms( MRI_DATA&,MRI_DATA&,int numarg,char **pstring);
	  
	  RECON(void); 	  
	  RECON(int numarg, char **pstring); 
	  static void help_message(void);
	  void parse_external_header(MRI_DATA &data);
	  void set_defaults(void);
	  void parse_commandline(int numarg, char **pstring);
	  void init_recon(int argc, char **argv, MRI_DATA& data );
	  void calc_sensitivity_maps( int argc, char **argv, MRI_DATA& data);
	
	  void L1_threshold( NDarray::Array< NDarray::Array< complex<float>,3>, 2 >&);
	  void transform_in_time( NDarray::Array< NDarray::Array< complex<float>,3>, 2 >&, TransformDirection);
	  void transform_in_encode( NDarray::Array< NDarray::Array< complex<float>,3>, 2 >&, TransformDirection);
	  
	  static double kspace_residual( MRI_DATA& data);
	  
	  NDarray::Array< NDarray::Array< complex<float>,3>, 2 > full_recon( MRI_DATA& data, NDarray::Range, NDarray::Range,bool);
	  NDarray::Array< NDarray::Array< complex<float>,3>, 1 > reconstruct_one_frame( MRI_DATA& data, int);
	  NDarray::Array< NDarray::Array< complex<float>,3>, 2 > reconstruct_all_frames( MRI_DATA& data);
	  NDarray::Array< NDarray::Array< complex<float>,3>, 1 > reconstruct_composite( MRI_DATA& data);
	  void eigen_coils( NDarray::Array< NDarray::Array< complex<float>,3 >,1 > &smaps,NDarray::Array< NDarray::Array< complex<float>,3 >,2 > &image);
	  void dcf_calc( MRI_DATA& data);
	  void dcf_calc( MRI_DATA& data, GATING& gate);
	  void gaussian_blur(  NDarray::Array< complex<float>, 3> & In, float,float,float);
	  void normalized_gaussian_blur( const  NDarray::Array< float, 3> & In,  NDarray::Array< float, 3> & out, float sigma);
	  void intensity_correct( NDarray::Array<float,3> &IC, NDarray::Array< NDarray::Array< complex<float>,3>,1 >&smaps );
	  void pregate_data( MRI_DATA&);
	  void sos_normalize(  NDarray::Array< NDarray::Array< complex<float>,3>, 1 >&);
		
	  void export_slice( NDarray::Array< complex<float>,3> &temp, const char * fname);
	private:	
		
};
