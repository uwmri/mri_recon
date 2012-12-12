#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstring>
#include <cmath>
#include <string>
#include <complex>
using namespace std;
#include <omp.h>
#include "ArrayTemplates.hpp"
#include <armadillo>

class MRI_DATA{
	public:
		/*Raw Data - 3rd Dimension is for encoding*/
		Array<float,4> kx;
		Array<float,4> ky;
		Array<float,4> kz;
		Array<float,4> kw;
		Array<float,4> kt;	  // TE Time (s)
		Array<float,4> times;// Time of readout for gating etc (s)
		Array< complex<float>,5> kdata;
		int Num_Encodings;
		int Num_Readouts;
		int Num_Slices;
		int Num_Pts;
		int Num_Coils;
		void read_external_data(char *folder,int NumRecv,int Ne,int Ns,int Npr,int Nx,int);
		void undersample(int);
		void coilcompress(float);
		
		// Pointer 
		void set_kx(float *,int e);
		void set_ky(float *,int e);
		void set_kz(float *,int e);
		void set_kw(float *,int e);
		void set_kdata(complex<float> *,int e,int coil);
		
		
		void init_memory();
		MRI_DATA( MRI_DATA *);
		MRI_DATA( void );
	private:	
		
};

