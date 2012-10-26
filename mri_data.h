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
#include "ArrayTemplates.cpp"

#include <armadillo>
using arma::cx_mat;
using arma::vec;
using arma::uvec;

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
		void read_external_data(char *folder,int NumRecv,int Ne,int Ns,int Npr,int Nx);
		void undersample(int);
		void coilcompress(float);
		MRI_DATA( MRI_DATA *);
		MRI_DATA( void );
	private:	
		
};

