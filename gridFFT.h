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

// Kernel Types
#define TRIANGLE 0
#define KAISER   1 

#define PI 3.14159265359

class gridFFT{
	public:
		int threads;
		
  		array3D< complex<float> >k3d_grid; 	/*Actual Gridding Space*/
  		array3D< complex<float> >image;   	/*Complex Storage Space*/
  		array3D< float >image_mag; /*Magnitude Storage Space*/
  
  		// Controls for phase encode / 2D directions
		int fft_in_x;
  		int fft_in_y;
  		int fft_in_z;
  		int grid_in_x;
  		int grid_in_y;
  		int grid_in_z;
  				
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
		int kernel_type;
		float overgrid;
		  
  		/*Kaiser Bessel Beta - Calculated*/
  		float betaX;
  		float betaY;
  		float betaZ;
  
  		// Discrete Gridding Kernel Variables
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
  		
		// FFT
  		fftwf_plan fft_plan;
		fftwf_plan ifft_plan;
		
		float k_rad; 
		
		int time_grid;
		gridFFT();
		~gridFFT();
		
		void alloc_grid();
		void read_commandline(int numarg, char **pstring);
		void precalc_gridding(int Nz,int Ny,int Nx,int directions);
		void deapp_chop();
		void forward( array3D< complex<float> >data, array3D< float >kx, array3D< float >ky, array3D< float >kz, array3D< float >kw);
		void backward( array3D< complex<float> >data, array3D< float >kx, array3D< float >ky, array3D< float >kz, array3D< float >kw);
		
		array3D< complex<float> >return_array( void);
		
		void chop_grid_forward( complex<float> *data, float *kx, float *ky, float *kz, float *kw,int);
		void chop_grid_backward( complex<float> *data, float *kx, float *ky, float *kz, float *kw,int);
		float bessi0(float);
		void plan_fft( void );
		void deapp_chop_crop(void);
		void icrop_deapp_chop(void);
		void chop(void);
		
	private:	
		
};

#endif

