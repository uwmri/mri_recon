#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <complex>
using namespace std;
#include <omp.h>
#include "ArrayTemplates.hpp"
#include <armadillo>
#include "gridFFT.h"

#ifndef PI
#define PI 3.14159265359
#endif

enum {FRACTAL_PHANTOM, SHEPP_PHANTOM};

class PHANTOM{
	public:
		int Nx;
		int Ny;
		int Nz;
		int Nt;
		void init(int,int,int);
		int phantom_type;
		float phantom_noise;
		
		PHANTOM(void);
		
		float over_res;
		Array< complex<float>,3>IMAGE;
		Array< complex<float>,3>SMAP;
		Array< complex<float>,3>TOA;
		
		void add_phase(void);
		void fractal3D(int Nx, int Ny, int Nz);
		void update_smap(int,int);
		void help_message(void);
		void add_noise( Array<complex<float>,5>&kdata);
		void read_commandline(int numarg, char **pstring);
	private:	
		
};


