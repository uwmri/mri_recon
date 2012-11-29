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
	 
	 TDIFF( Array< complex<float>,5 >&temp);
	 cx_mat At;		// forward in time
	 cx_mat AIt;	// inverse in time
	 cx_mat Ae; 	// forward in encode
	 cx_mat AIe; 	// inverse in encoee
	 int Nt;
	 int Ne;
	 
	 
	 void fft_e(Array< complex<float>,5>&temp);	 
	 void ifft_e(Array< complex<float>,5>&temp);	 
	 
	 void fft_t(Array< complex<float>,5>&temp);	 
	 void ifft_t(Array< complex<float>,5>&temp);	 
	 
	 
	 void ediff(Array< complex<float>,5>&temp);	 
	 void inv_ediff(Array< complex<float>,5>&temp);	 
	 
	 void tdiff(Array< complex<float>,5>&temp);	 
	 void inv_tdiff(Array< complex<float>,5>&temp);	 
	 
	private:	
		
};



#endif

