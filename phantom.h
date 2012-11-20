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

enum {FRACTAL_PHANTOM, SHEPP_PHANTOM};

class PHANTOM{
	public:
		int Nx;
		int Ny;
		int Nz;
		int Nt;
		PHANTOM(int,int,int);
		
		Array< float,3 > fractal3D(int Nx, int Ny, int Nz);
	private:	
		
};


