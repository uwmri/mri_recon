#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string.h>
#include <complex>
#include <omp.h>
#include "ArrayTemplates.hpp"
#include <fftw3.h>
#include <sys/unistd.h>
#include "mri_data.h"

// Kernel Types
enum {TRIANGLE_KERNEL, KAISER_KERNEL, SINC_KERNEL};

#ifndef PI
#define PI 3.14159265359
#endif

class gridFFT{
	public:
		int threads;
		
  		NDarray::Array< complex<float>,3 >k3d_grid; 	/*Actual Gridding Space*/
  		NDarray::Array< complex<float>,3 >image;   	/*Complex Storage Space*/
  		
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
		
		// arma::matrix Size
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
  		NDarray::Array<float,1> winx;
  		NDarray::Array<float,1> winy;
  		NDarray::Array<float,1> winz;
  		float dwinX;
  		float dwinY;
  		float dwinZ;
  		float grid_modX;
  		float grid_modY;
  		float grid_modZ;
  		NDarray::Array<float,1> grid_filterX;
  		NDarray::Array<float,1> grid_filterY;
  		NDarray::Array<float,1> grid_filterZ;
  		
		// FFT
  		fftwf_plan fft_plan;
		fftwf_plan ifft_plan;
		
		float k_rad; 
		int time_grid;
		
		gridFFT();
		~gridFFT();
		
		void alloc_grid();
		void read_commandline(int numarg, char **pstring);
		void precalc_gridding(int Nz,int Ny,int Nx,TrajDim trajectory_dims,TrajType trajectory_type);
		void precalc_kernel(void);
		void deapp_chop();
		void do_fft( void);
		void do_ifft( void );
		
		// Main Calls with and without Sensitivity maps
		void forward( NDarray::Array< complex<float>,3 >&X,const NDarray::Array< complex<float>,3 >&smap,const NDarray::Array< complex<float>,3 >&data, const NDarray::Array< float,3 >&kx,const NDarray::Array< float,3 >&ky, const NDarray::Array< float,3 >&kz, const NDarray::Array< float,3 >&kw);
		void forward( NDarray::Array< complex<float>,3 >&X,const NDarray::Array< float,3 >&smap,const NDarray::Array< complex<float>,3 >&data, const NDarray::Array< float,3 >&kx,const NDarray::Array< float,3 >&ky, const NDarray::Array< float,3 >&kz, const NDarray::Array< float,3 >&kw);
		void forward( NDarray::Array< complex<float>,3 >&X,const NDarray::Array< complex<float>,3 >&data, const NDarray::Array< float,3 >&kx,const NDarray::Array< float,3 >&ky, const NDarray::Array< float,3 >&kz, const NDarray::Array< float,3 >&kw);
		
		void backward( const NDarray::Array< complex<float>,3 >&X,const NDarray::Array< complex<float>,3 >&smap,NDarray::Array< complex<float>,3 >&data, const NDarray::Array< float,3 >&kx, const NDarray::Array< float,3 >&ky, const NDarray::Array< float,3 >&kz, const NDarray::Array< float,3 >&kw);
		void backward( const NDarray::Array< complex<float>,3 >&X,const NDarray::Array< float,3 >&smap,NDarray::Array< complex<float>,3 >&data, const NDarray::Array< float,3 >&kx, const NDarray::Array< float,3 >&ky, const NDarray::Array< float,3 >&kz, const NDarray::Array< float,3 >&kw);
		void backward( const NDarray::Array< complex<float>,3 >&X,NDarray::Array< complex<float>,3 >&data, const NDarray::Array< float,3 >&kx, const NDarray::Array< float,3 >&ky, const NDarray::Array< float,3 >&kz, const NDarray::Array< float,3 >&kw);
						
		void chop_grid_forward( const NDarray::Array< complex<float>,3 >&data, const NDarray::Array< float,3 >&kx, const NDarray::Array< float,3 >&ky, const NDarray::Array< float,3 >&kz, const NDarray::Array< float,3 >&kw);
		void chop_grid_backward( NDarray::Array< complex<float>,3 >&data, const NDarray::Array< float,3 >&kx, const NDarray::Array< float,3 >&ky, const NDarray::Array< float,3 >&kz, const NDarray::Array< float,3 >&kw);
		float bessi0(float);
		void plan_fft( void );
		
		// grid for dcf 
		void grid_forward( NDarray::Array< float, 3 >&X, const NDarray::Array< float,3 >&data, const NDarray::Array< float,3 >&kx,const NDarray::Array< float,3 >&ky, const NDarray::Array< float,3 >&kz); 
		void grid_backward( const NDarray::Array< float, 3 >&X, NDarray::Array< float,3 >&data, const NDarray::Array< float,3 >&kx,const NDarray::Array< float,3 >&ky, const NDarray::Array< float,3 >&kz); 

		// Copy Gridding to Image 
		void forward_image_copy(NDarray::Array< complex<float>,3 >&X);
		void forward_image_copy(NDarray::Array< complex<float>,3 >&X,const NDarray::Array< complex<float>,3 >&smap);
		void forward_image_copy(NDarray::Array< complex<float>,3 >&X,const NDarray::Array< float,3 >&smap);
		
		
		// Copy Image to Gridding
		void backward_image_copy(const NDarray::Array< complex<float>,3 >&X);
		void backward_image_copy(const NDarray::Array< complex<float>,3 >&X,const NDarray::Array< complex<float>,3 >&smap);
		void backward_image_copy(const NDarray::Array< complex<float>,3 >&X,const NDarray::Array< float,3 >&smap);
		
		static void help_message(void);
		
	private:	
		
};


