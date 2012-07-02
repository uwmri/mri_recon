#ifndef hMRIDATA
#define hMRIDATA

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


class MRI_DATA{
	public:
		/*Raw Data*/
		array5D kx;
		array5D ky;
		array5D kz;
		array5D kw;
		
		int dims[5];	
		
	

	private:	
		
};

#endif
