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

class RECON{
	public:
	  // Types of Recons
	  enum ReconType { SOS, PILS, CG, IST, FISTA, CLEAR};
	   
	  // Data Types
	  enum DataType { PFILE, EXTERNAL, SIMULATE, PSF, PHANTOM };

	  // Coil Combine Type
	  enum CoilCombineType { LOWRES, ESPIRIT, WALSH };
	  
	  // Enum Transform Types
	  enum TransformType {NONE,WAVELET,DIFF,DFT,PCA,COMPOSITE_DIFF};
	  
	  enum TransformDirection {FORWARD,BACKWARD};
	  	 
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
	  	  
	  bool complex_diff;
	  
	  float zoom;
	  float zoom_x;
	  float zoom_y;
	  float zoom_z;
	  
	  // Reconstructed Images
	  int rcxres;
	  int rcyres;
	  int rczres;
	  int rcframes;
	  int rcencodes;
  
  	  // Gridding operators
	  int grid_workers;
	  int grid_threads;
	  NDarray::Array< gridFFT,1> gridding;
	  
	  // Arrays for calculations
	  NDarray::Array< NDarray::Array< complex<float>,3>,1 >smaps; // Array of arrays
	  NDarray::Array< complex<float>,3>composite_image;
	  NDarray::Array< float, 3 >IntensityCorrection;
	  GATING gate;
	  	  
	  // Density calcs
	  bool recalc_dcf;
	  int dcf_iter;
	  float dcf_dwin;
	  float dcf_scale;
	  
	  int acc;
	  float compress_coils;
	  bool whiten;
	  float smap_res;
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
	  
	  int wavelet_levelsX; 
	  int wavelet_levelsY; 
	  int wavelet_levelsZ; 
	  	  
		  
	  int max_iter;
	  int export_smaps;
	  bool prep_done;
	  bool pregate_data_flag;
	  
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
	  
	  NDarray::Array< NDarray::Array< complex<float>,3>, 2 >full_recon( MRI_DATA& data, NDarray::Range, NDarray::Range,bool);
	  NDarray::Array< NDarray::Array< complex<float>,3>, 1 >reconstruct_one_frame( MRI_DATA& data, int);
	  NDarray::Array< NDarray::Array< complex<float>,3>, 2 >reconstruct_all_frames( MRI_DATA& data);
	  NDarray::Array< NDarray::Array< complex<float>,3>, 1 >reconstruct_composite( MRI_DATA& data);
	  void eigen_coils( NDarray::Array< NDarray::Array< complex<float>,3 >,1 > &image);
	  void dcf_calc( MRI_DATA& data);
	  void dcf_calc( MRI_DATA& data, GATING& gate);
	  void gaussian_blur(  NDarray::Array< complex<float>, 3> & In, float,float,float);
	  void normalized_gaussian_blur( const  NDarray::Array< float, 3> & In,  NDarray::Array< float, 3> & out, float sigma);
	  void intensity_correct( NDarray::Array<float,3> &IC, NDarray::Array< NDarray::Array< complex<float>,3>,1 >&smaps );
	  void pregate_data( MRI_DATA&);
	  void single_coil_cg( NDarray::Array< complex<float>,3> & X, NDarray::Array< complex<float>,3> &kdata, NDarray::Array<float,3> &kx, NDarray::Array<float,3> &ky, NDarray::Array<float,3> &kz, NDarray::Array<float,3> &kw);

	  void export_slice( NDarray::Array< complex<float>,3> &temp, const char * fname);
	private:	
		
};

