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
		float thresh;
		void soft_threshold( array5D< complex<float> >Coef);
		void get_threshold( array5D< complex<float> >Coef);
		float threshold;
		
	private:	
		
};


#endif

