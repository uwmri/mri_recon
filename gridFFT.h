#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string.h>
#include <complex>
using namespace std;
#include <omp.h>
#include "ArrayTemplates.hpp"
#include <fftw3.h>
#include <sys/unistd.h>

// Kernel Types
enum {TRIANGLE_KERNEL, KAISER_KERNEL, SINC_KERNEL};

#ifndef PI
#define PI 3.14159265359
#endif

class gridFFT{
	public:
		int threads;
		
  		Array< complex<float>,3 >k3d_grid; 	/*Actual Gridding Space*/
  		Array< complex<float>,3 >image;   	/*Complex Storage Space*/
  		
		// Controls for phase encode / 2D directions
		int fft_in_x;
  		int fft_in_y;
  		int fft_in_z;
  		int grid_in_x;
  		int grid_in_y;
  		int grid_in_z;
  		
		// Grid Testing for Double 
		int double_grid;
				
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
		
		// Main Calls with and without Sensitivity maps
		void forward( Array< complex<float>,3 >&X,const Array< complex<float>,3 >&smap,const Array< complex<float>,3 >&data, const Array< float,3 >&kx,const Array< float,3 >&ky, const Array< float,3 >&kz, const Array< float,3 >&kw);
		void forward( Array< complex<float>,3 >&X,const Array< complex<float>,3 >&data, const Array< float,3 >&kx,const Array< float,3 >&ky, const Array< float,3 >&kz, const Array< float,3 >&kw);
		
		void backward( const Array< complex<float>,3 >&X,const Array< complex<float>,3 >&smap,Array< complex<float>,3 >&data, const Array< float,3 >&kx, const Array< float,3 >&ky, const Array< float,3 >&kz, const Array< float,3 >&kw);
		void backward( const Array< complex<float>,3 >&X,Array< complex<float>,3 >&data, const Array< float,3 >&kx, const Array< float,3 >&ky, const Array< float,3 >&kz, const Array< float,3 >&kw);
						
		void chop_grid_forward( const Array< complex<float>,3 >&data, const Array< float,3 >&kx, const Array< float,3 >&ky, const Array< float,3 >&kz, const Array< float,3 >&kw);
		void chop_grid_backward( Array< complex<float>,3 >&data, const Array< float,3 >&kx, const Array< float,3 >&ky, const Array< float,3 >&kz, const Array< float,3 >&kw);
		float bessi0(float);
		void plan_fft( void );
		
		// Copy Gridding to Image 
		void forward_image_copy(Array< complex<float>,3 >&X);
		void forward_image_copy(Array< complex<float>,3 >&X,const Array< complex<float>,3 >&smap);
		
		// Copy Image to Gridding
		void backward_image_copy(const Array< complex<float>,3 >&X);
		void backward_image_copy(const Array< complex<float>,3 >&X,const Array< complex<float>,3 >&smap);
		
		
		static void help_message(void);
		
	private:	
		
};


