#pragma once 

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstring>
#include <string>
#include <complex>
using namespace std;
#include <omp.h>
#include <armadillo>
using arma::cx_mat;
using arma::vec;
using arma::uvec;

// Local Libraries
#include "ArrayTemplates.hpp"
#include "wavelet3D.h"	
#include "temporal_diff.h"	
#include "gridFFT.h"
// #include "ge_pfile_lib.h"
#include "mri_data.h"
#include "threshold.h"
#include "spirit.h"
#include "phantom.h"
#include "tictoc.cpp"

// Types of Recons
enum { RECON_SOS, RECON_PILS, RECON_CG, RECON_IST, RECON_FISTA };

// Data Types
enum { PFILE_DATA, EXTERNAL_DATA, SIMULATE_DATA, PHANTOM_DATA };

// Coil Combine Type
enum { COIL_LOWRES, COIL_ESPIRIT };

class RECON{
	public:
	  int recon_type;
	  int data_type;
	  
	  int numrecv;
	  float zero_fill;
	  bool complex_diff;
	  
	  float zoom;
	  float zoom_x;
	  float zoom_y;
	  float zoom_z;
	  
	  int rcxres;
	  int rcyres;
	  int rczres;
	  int rcframes;
	  int rcencodes;
	  
	  float acq_bw;
	  int xres;
	  int num_readouts;
	  int num_coils;
	  int num_slices;
	  int ss_2d;
	  int multi_echo;
	    
	  int acc;
	  float compress_coils;
	  
	  float lp_frac;
	  float lp_sig;
      float smap_res;
	  char filename[1024];
	  
	  int max_iter;
	  
	  int coil_combine_type;
	  int export_smaps;
	  
	  RECON(void); 	  
	  RECON(int numarg, char **pstring); 
	  static void help_message(void);
	  void parse_external_header(void);
	  void set_defaults(void);
	  void parse_commandline(int numarg, char **pstring);
	  
	  Array< complex<float>, 5 >reconstruction( int argc, char **argv, MRI_DATA& data);
	private:	
		
};

