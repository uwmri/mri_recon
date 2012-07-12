#ifndef hPCVIPRLIB
#define hPCVIPRLIB

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
#include "ge_pfile_lib.h"

class PCVIPR{
	public:
		int nproj;
		void interpret_pfile_header(PFILE pfile);
		void get_kspace_trajectory(PFILE pfile);
	private:	
		

};


#endif

