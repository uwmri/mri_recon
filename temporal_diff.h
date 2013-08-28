#ifndef hTDIFFLIB
#define hTDIFFLIB

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <complex>
using namespace std;
#include <omp.h>
#include "ArrayTemplates.hpp"
#include <fftw3.h>

#include <armadillo>
using arma::cx_mat;
using arma::vec;
using arma::uvec;
using arma::cx_vec;


class TDIFF{
	public:
	 
	 TDIFF( Array< Array< complex<float>,3>,2>&temp);	 
	 cx_mat At;		// forward in time
	 cx_mat AIt;	// inverse in time
	 cx_mat Ae; 	// forward in encode
	 cx_mat AIe; 	// inverse in encoee
	 int Nt;
	 int Ne;
	 
	 
	 void fft_e( Array< Array< complex<float>,3>, 2>&temp);	 
	 void ifft_e(Array< Array< complex<float>,3>, 2>&temp);	 
	 
	 void fft_t(Array< Array< complex<float>,3>,2>&temp);	 
	 void ifft_t(Array< Array< complex<float>,3>,2>&temp);	 
	 
	 
	 void ediff(Array< Array< complex<float>,3>,2>&temp);	 
	 void inv_ediff(Array< Array< complex<float>,3>,2>&temp);	 
	 
	 void tdiff(Array< Array< complex<float>,3>,2>&temp);	 	 
	 void inv_tdiff(Array< Array< complex<float>,3>,2>&temp);	 
	 
	private:	
		
};



#endif

