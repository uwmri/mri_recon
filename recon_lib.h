#ifndef hRECONLIB
#define hRECONLIB

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
#include "ArrayTemplates.cpp"

#include <armadillo>
using arma::cx_mat;
using arma::vec;
using arma::uvec;


#define THREADS 32

// Types of Recons
#define RECON_SOS 0
#define RECON_PILS 1
#define RECON_CG 2
#define RECON_IST 3
#define RECON_FISTA 4

// Data Types
#define RECON_PFILE 0
#define RECON_EXTERNAL 1

class RECON{
	public:
	  int recon_type;
	  int data_type;
	  
	  int numrecv;
	  float zero_fill;
	  
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
	  
	  
	  RECON(void); 	  
	  RECON(int numarg, char **pstring); 
	  static void help_message(void);
	  void parse_external_header(void);
	  void set_defaults(void);
	  void parse_commandline(int numarg, char **pstring);
	private:	
		
};


#endif
