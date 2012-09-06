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
	 
	 void forward(array5D< complex<float> >temp);
	 void backward(array5D< complex<float> >temp);
	private:	
		
};



#endif

