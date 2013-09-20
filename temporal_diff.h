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


class TDIFF{
	public:
	 
	 TDIFF();
	 TDIFF( NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);
	 TDIFF( int, int);
	 	
	 arma::cx_mat AWt;	// forward in time
	 arma::cx_mat AWIt;	// inverse in time 
	 
	 arma::cx_mat At;	// forward in time
	 arma::cx_mat AIt;	// inverse in time
	 
	 arma::cx_mat Ae; 	// forward in encode
	 arma::cx_mat AIe; 	// inverse in encode
	 
	 int Nt;
	 int Ne;
	 
	 
	 void fft_e( NDarray::Array< NDarray::Array< complex<float>,3>, 2>&temp);	 
	 void ifft_e(NDarray::Array< NDarray::Array< complex<float>,3>, 2>&temp);	 
	 
	 void fft_t(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 void ifft_t(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 
	 
	 void ediff(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 void inv_ediff(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 
	 void tdiff(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 	 
	 void inv_tdiff(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 
	 void twave(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 	 
	 void inv_twave(NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp);	 
	 
	 
	 void mat_multiply( NDarray::Array< NDarray::Array< complex<float>,3>,2>&temp,arma::cx_mat);
	 
	private:	
		
};



#endif

