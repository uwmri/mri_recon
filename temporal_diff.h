#ifndef hTDIFFLIB
#define hTDIFFLIB

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <complex>
#include <omp.h>
#include "ArrayTemplates.hpp"
#include <fftw3.h>

 
#include <armadillo>


class TRANSFORMS{
	public:
	 TRANSFORMS();
	 
	 enum TransformDirection {FORWARD,BACKWARD};
	 
	 /* Low res*/
	 void eigen( NDarray::Array< NDarray::Array< complex<float>,3>, 2>&temp,int,int);
	 
	 static void get_difference_transform( arma::cx_mat & A, arma::cx_mat & Ai, int N);
	 static void get_wavelet_transform( arma::cx_mat & A, arma::cx_mat & Ai, int N);
	 	  
	 static void fft_e( NDarray::Array< NDarray::Array< complex<float>,3>, 2>&temp);	 
	 static void ifft_e(NDarray::Array< NDarray::Array< complex<float>,3>, 2>&temp);	 
	 
	 static void fft_t(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 static void ifft_t(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 
	 static void ediff(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 static void inv_ediff(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 
	 static void ewave(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 static void inv_ewave(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 
	 static void tdiff(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 	 
	 static void inv_tdiff(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 
	 static void twave(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 	 
	 static void inv_twave(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 	 
	 static void multiply_in_time( NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp,arma::cx_mat);
	 static void multiply_in_encode( NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp,arma::cx_mat);
	 
	private:	
	 /* Eigen Matrices*/
	 arma::cx_mat E;
	 arma::cx_mat Ei;
	 bool reinit_pca;
	 int pca_count;	 	 	
};



#endif

