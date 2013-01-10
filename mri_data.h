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

// Data Types
enum TrajDim { THREED, TWOD };
enum TrajType { CARTESIAN, NONCARTESIAN};

class MRI_DATA{
	public:
		// Raw Data
		//  readout x phase encode x slice x encoding x coil 
		
		// Non-Cartesian Trajectory
		Array<float,4> kx;
		Array<float,4> ky;
		Array<float,4> kz;
		Array<float,4> kw;
		Array<float,4> kt;	  // TE Time (s)
		Array<float,4> times; // Time of readout for gating etc (s)
		Array< complex<float>,5> kdata;
		
		// Sizes
		int Num_Encodings;
		int Num_Readouts;
		int Num_Slices;
		int Num_Pts;
		int Num_Coils;
		
		// Native Resolution
		int xres;
		int yres;
		int zres;
						
		// 2D/3D Cartesian/Non-Cartesian
		TrajDim trajectory_dims;
		TrajType trajectory_type;
						
		// Data Operations (move?)
		void undersample(int);
		void coilcompress(float);
		
		// Initialization Filling Operations				
		void init_memory();
		void read_external_data(char *folder,int);
		MRI_DATA( MRI_DATA *);
		MRI_DATA( void );
	private:	
		
};

