#ifndef hgridFFTLIB
#define hgridFFTLIB

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string.h>
#include <complex>
using namespace std;
#include <omp.h>
#include "ArrayTemplates.cpp"
#include <fftw3.h>

#define THREADS 32

// Kernel Types
#define TRIANGLE 0
#define KAISER   1 

#define PI 3.14159265359

class gridFFT{
	public:
		int threads;
		
  		array3D< complex<float> >k3d_grid; 	/*Actual Gridding Space*/
  		array3D< complex<float> >k3dref;   	/*Complex Storage Space*/
  		array3D< float >image_mag; /*Magnitude Storage Space*/
  
  		/*Overgridding Factor*/   
  		float grid_x;
  		float grid_y;
  		float grid_z;
  
  		/*Overgridding Crop values*/
  		int og_sx;
  		int og_sy;
  		int og_sz;
  		int og_ex;
  		int og_ey;
  		int og_ez;
		
		// Matrix Size
		int Nx;
		int Ny;
		int Nz;
		
		// Grid Size
		int Sx;
		int Sy;
		int Sz;
		  
  		/*Kaiser Bessel Beta - Calculated*/
  		float betaX;
  		float betaY;
  		float betaZ;
  
  		/**Gridding Variables*/
  		float *winx;
  		float *winy;
  		float *winz;
  
  		float dwinX;
  		float dwinY;
  		float dwinZ;
  
  		float grid_modX;
  		float grid_modY;
  		float grid_modZ;
  
  		float *grid_filterX;
  		float *grid_filterY;
  		float *grid_filterZ;
  		fftwf_plan fft_plan;
		fftwf_plan ifft_plan;
		int kernel_type;
		float overgrid;
		
		float k_rad; 
		
		gridFFT();
		~gridFFT();
		
		void alloc_grid();
		void read_commandline(int numarg, char **pstring);
		void precalc_gridding(int Nz,int Ny,int Nx,int directions);
		void deapp_chop();
		void forward( complex<float> *data, float *kx, float *ky, float *kz, float *kw,int);
		void backward( complex<float> *data, float *kx, float *ky, float *kz, float *kw,int);
		
		array3D< complex<float> >return_array( void);
		
		void grid_forward( complex<float> *data, float *kx, float *ky, float *kz, float *kw,int);
		void grid_backward( complex<float> *data, float *kx, float *ky, float *kz, float *kw,int);
		float bessi0(float);
		void plan_fft( void );
		void deapp_chop_crop(void);
		void icrop_deapp_chop(void);
		void chop(void);
		
	private:	
		
};

#endif

