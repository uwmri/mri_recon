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
#include <sys/stat.h>

// Data Types
enum TrajDim { THREED, TWOD };
enum TrajType { CARTESIAN, NONCARTESIAN};

class MRI_DATA{
	public:
		// Raw Data
		//  readout x phase encode x slice x encoding x coil 
		
		// Sizes for shorthand
		int Num_Encodings;
		int Num_Readouts;
		int Num_Slices;
		int Num_Pts;
		int Num_Coils;
		
		// Non-Cartesian Trajectory
		Array<float,4> kx;	// Fov = 1 unit, delta k =1
		Array<float,4> ky;
		Array<float,4> kz;
		Array<float,4> kw;
		Array<float,4> kt;	  // TE Time (s)
		Array< complex<float>,5> kdata;
		
		// Data for Noise samples 
		Array< complex<float>,2> noise_samples; // data for noise samples
				
		//Physiologic Data for gating 
		Array< float,3>ecg;		// Distance from ECG in MS
		Array< float,3>resp;	// Respiratory signal from bellows or navigator
		Array< float,3>time;	// Acquisition Time 
		Array< float,3>prep;	// Time since a prep event (for example inversion)
				
		// Native Resolution
		int xres;
		int yres;
		int zres;
						
		// 2D/3D Cartesian/Non-Cartesian
		TrajDim trajectory_dims;
		TrajType trajectory_type;
		
		//Temp
		char gate_name[1024];
		
		// Data Operations (move?)
		void undersample(int);
		void coilcompress(float);
		void whiten();
		
		// Initialization Filling Operations				
		void init_memory();
		void init_noise_samples(int);
		void read_external_data(const char *folder,int);
		void write_external_data(const char *folder);
		void parse_external_header(const char *filename);
		void load_pcvipr_gating_file(const char *full_filename); //Temp
		void data_stats(void);
				
		MRI_DATA( MRI_DATA *);
		MRI_DATA( void );
	private:	
		
};

