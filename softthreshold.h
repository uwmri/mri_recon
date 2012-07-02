#ifndef hSOFTTHRESHLIB
#define hSOFTTHRESHLIB

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <string.h>
#include <complex.h>
using namespace std;
#include <omp.h>
#include "ArrayTemplates.cpp"

class SOFTTHRESHOLD{
	public:
		// Declare Functions/variables here		
		SOFTTHRESHOLD(int numarg, char **pstring);
		float thresh;
		void hard_threshold( array5D< complex<float> >Coef);
		
	private:	
		
};


#endif

