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
	  int ss_2d;
	  int multi_echo;
	  
	  int acc;
	  float compress_coils;
	  
	  float lp_frac;
	  float lp_sig;
    float smap_res;
	  char filename[1024];
	  
	  int frames; 
	  
    int sp_maps;
    
	  RECON(int numarg, char **pstring); 
	  void parse_external_header(void);
	private:	
		
};


class MRI_DATA{
	public:
		/*Raw Data - 3rd Dimension is for encoding*/
		array3D<float> kx;
		array3D<float> ky;
		array3D<float> kz;
		array3D<float> kw;
		array4D< complex<float> > kdata;
		int Num_Encodings;
		int Num_Readouts;
		int Num_Pts;
		int Num_Coils;
		void read_external_data(char *folder,int NumRecv,int Ne,int Npr,int Nx);
		void undersample(int);
		void coilcompress(float);
		MRI_DATA( MRI_DATA *);
		MRI_DATA( void );
	private:	
		
};


#endif
