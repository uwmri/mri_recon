#pragma once
#ifndef hGEPFILE
#define GEPFILE

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <complex>
#include <string.h>
#include <omp.h>

// For Wrapper (this is the actual recon and also contains Array library)
#include <recon_lib.h>

// Function Switches for Include
#ifndef hUW_LIB
#define hUW_LIB

#if defined(ESE11_RECON)
   #include <uw_lnx9.h>
#elif defined(ESE12_RECON)
   #include <uw_lnx11.h> 
#elif defined(ESE14_RECON) 
	#include <uw_lnx14.h> 
#elif defined(ESE15_RECON) 
	#include <uw_lnx15.h> 
#elif defined(ESE16_RECON)	
	#include <uw_lnx16.h> 
#elif defined(ESE20_RECON) 
	#include "uw_lnx20.h" 
#elif defined(ESE21_RECON) 
	#include <uw_lnx21.h> 
#elif defined(ESE22_RECON)
	#include <uw_lnx22.h> 
#elif defined(ESE23_RECON)
	#include <uw_lnx23.h> 
#elif defined(ESE24_RECON)
	#include <uw_lnx24.h> 
#endif

#endif

/*************************************
scan_table 
  Structure holds image header info. Points
  defined by ux,uy,uz the "unit" vectors which
  have magnitude equal to spacing and sx,sy,sz 
  the start position 
**************************************/	
typedef struct { 
	float ix;
	float iy;
	float iz;
	float jx;
	float jy;
	float jz;
	float kx;
	float ky;
	float kz;
	float sx;
	float sy;
	float sz;
}scan_table;


// Tools For Reading in P-File
class PFILE{
	public:
		
		// Define Type
		void read_header(char filename[]);		
		char filename[1024];
				
		RDB_HEADER_REC rdbhead;
		RDB_DATA_ACQ_TAB acq_tab;
		EXAM examhead;
		SERIES serieshead;
		IMAGE imagehead;
		
		// For scan table logic;
		scan_table logical_tbl;
		scan_table physical_tbl;
		
		void read_data(complex<float> *dat,int coil, int offset, int number, int stride);
		void read_data(complex<float> *dat,int coil, int offset, int number, int stride, int );
		void write_dicom(NDarray::Array< complex<float>,3> &X, int lx, const string series_description, float max_all,float,float,float);
		void write_dicom(NDarray::Array< float,3> &X, int lx, const string series_description, float max_all,float,float,float);
		void calc_unit_vectors( int,int,int,float,float,float);
		
		
	private:	
		
};

#endif


