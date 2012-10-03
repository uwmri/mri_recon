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
#include "ArrayTemplates.cpp"
#include <fftw3.h>

#include <armadillo>
using arma::cx_mat;
using arma::vec;
using arma::uvec;
using arma::cx_vec;


class TDIFF{
	public:
	 
	 TDIFF( array5D< complex<float> >temp);
	 cx_mat At;		// forward in time
	 cx_mat AIt;	// inverse in time
	 cx_mat Ae; 	// forward in encode
	 cx_mat AIe; 	// inverse in encoee
	 int Nt;
	 int Ne;
	 
	 
	 void fft_e(array5D< complex<float> >temp);	 
	 void ifft_e(array5D< complex<float> >temp);	 
	 
	 void fft_t(array5D< complex<float> >temp);	 
	 void ifft_t(array5D< complex<float> >temp);	 
	 
	 
	 void forward(array5D< complex<float> >temp);
	 void backward(array5D< complex<float> >temp);
	 
	 void tdiff(array5D< complex<float> >temp);
	 void inv_tdiff(array5D< complex<float> >temp);
	 
	private:	
		
};



#endif

