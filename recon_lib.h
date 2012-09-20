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
	  
	  int frames; 
	  	  
	  RECON(int numarg, char **pstring); 
	  void parse_external_header(void);
	private:	
		
};


class MRI_DATA{
	public:
		/*Raw Data - 3rd Dimension is for encoding*/
		array4D<float> kx;
		array4D<float> ky;
		array4D<float> kz;
		array4D<float> kw;
		array5D< complex<float> > kdata;
		int Num_Encodings;
		int Num_Readouts;
		int Num_Slices;
		int Num_Pts;
		int Num_Coils;
		void read_external_data(char *folder,int NumRecv,int Ne,int Ns,int Npr,int Nx);
		void undersample(int);
		void coilcompress(float);
		MRI_DATA( MRI_DATA *);
		MRI_DATA( void );
	private:	
		
};


#endif
