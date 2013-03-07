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
#include "tictoc.cpp"

#ifndef PI
#define PI 3.14159265359
#endif

enum PhantomType {FRACTAL, SHEPP, PSF, EXTERNAL};

class PHANTOM{
	public:
		
		enum PerfType { ASL,ELLIPSE};
		
		// Size of Phantom
		int Nx;
		int Ny;
		int Nz;
		int Nt;
		float over_res; // Create higher resolution
		
		// Noise Factor
		float phantom_noise;
		
		// Type of Phantom
		PhantomType phantom_type;
		
		// Fractal Input
		int fractal_pts;
		
		char *external_phantom_name;
		
		// Constructor,I/O, and Init				
		PHANTOM(void);
		void init(int,int,int);
		static void help_message(void);
		void read_commandline(int numarg, char **pstring);
		
		// Arrays for holding Images
		Array< complex<float>,3>IMAGE;
		Array< float,4>FUZZY;
		Array< float,4>FUZZYT;
		Array< complex<float>,3>SMAP;
		Array< float,3>TOA;
		
		// Functions
		void add_phase(void);
		void add_noise( Array<complex<float>,5>&kdata);
		void fractal3D(int Nx, int Ny, int Nz);
		void fractal3D_new(int Nx, int Ny, int Nz);
		void update_smap_biotsavart(int,int);
		void calc_image(int,int);
		Array<int,3> synthetic_perfusion(int xs, int ys, int zs, PerfType ptype);
	private:	
		bool debug;
};


