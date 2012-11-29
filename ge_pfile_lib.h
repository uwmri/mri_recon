#ifndef hGELIB
#define hGELIB

/***ESE version control - current is default*/
#ifndef hUW_LIB
#define hUW_LIB
#ifdef ESE11_RECON
   #include <uw_lnx9.h>
#endif 
#ifdef ESE12_RECON 
	#include <uw_lnx11.h> 
#endif
#ifdef ESE14_RECON 
	#include <uw_lnx14.h> 
#endif
#ifdef ESE15_RECON 
	#include <uw_lnx15.h> 
#endif
#ifdef ESE20_RECON 
	#include "uw_lnx20.h" 
#endif
#ifdef ESE21_RECON 
	#include <uw_lnx21.h> 
#endif
#ifdef ESE22_RECON 
	#include "GE_LIB/idbm22.h"
#endif
#ifdef ESE16_RECON 
	#include <uw_lnx16.h> 
#endif
#endif

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <complex>
#include <string.h>
#include <omp.h>
#include "ArrayTemplates.hpp"
using namespace std;

class PFILE{
	public:
		void read_header(char filename[]);		
		char filename[1024];
		
		RDB_HEADER_REC rdbhead;
		RDB_DATA_ACQ_TAB acq_tab;
		EXAM examhead;
		SERIES serieshead;
		IMAGE imagehead;
		
		void read_data(int coil);
		array3D< complex<int> >RawData;
	private:	
		
};


#endif


