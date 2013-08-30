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
#include "tictoc.cpp"
#include "gating.h"
#include "clear.h"
#include "l2reg.h"

class RECON{
	public:
	  // Types of Recons
	  enum ReconType { SOS, PILS, CG, IST, FISTA, CLEAR};
	   
	  // Data Types
	  enum DataType { PFILE, EXTERNAL, SIMULATE, PSF, PHANTOM };

	  // Coil Combine Type
	  enum CoilCombineType { LOWRES, ESPIRIT };
	  
	  // Enum Transform Types
	  enum TransformType {NONE,WAVELET,DIFF,DFT,PCA};
	  	 
	  // Recon Flags  
	  ReconType recon_type;
	  DataType data_type;
	  CoilCombineType coil_combine_type;
	  
	  // CS Flags
	  TransformType cs_spatial_transform;
	  TransformType cs_temporal_transform;
	  TransformType cs_encode_transform;
	  	  
	  int numrecv;
	  float zero_fill;
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
  
  	  // Store for later
      gridFFT gridding;
	  NDarray::Array< NDarray::Array< complex<float>,3>,1 >smaps; // Array of arrays
	  GATING gate;
	  
	  int acc;
	  float compress_coils;
	  bool whiten;
	  float lp_frac;
	  float lp_sig;
      float smap_res;
	  char filename[1024];
	  
	  int max_iter;
	  int export_smaps;
	  bool prep_done;
	  
	  RECON(void); 	  
	  RECON(int numarg, char **pstring); 
	  static void help_message(void);
	  void parse_external_header(MRI_DATA &data);
	  void set_defaults(void);
	  void parse_commandline(int numarg, char **pstring);
	  void init_recon(int argc, char **argv, MRI_DATA& data );
	  
	  NDarray::Array< NDarray::Array< complex<float>,3>, 2 >reconstruction( int argc, char **argv, MRI_DATA& data);
	private:	
		
};

