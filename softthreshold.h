#ifndef hSOFTTHRESHLIB
#define hSOFTTHRESHLIB

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

class SOFTTHRESHOLD{
	public:
		// Declare Functions/variables here		
		SOFTTHRESHOLD(int numarg, char **pstring);
		double thresh;
		void soft_threshold( Array< complex<float>,5 >&Coef);
		void get_threshold( const Array< complex<float>,5 >&Coef);
		float threshold;
		
		// Iterative soft thresholding code
		void fista_update( Array< complex<float>,5 >&X,Array< complex<float>,5 >&X_old,int iteration);
		
	private:	
		
};


#endif

