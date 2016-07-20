#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstring>
#include <cmath>
#include <string>
#include <complex>
#include <omp.h>
#include "ArrayTemplates.hpp"
#define ARMA_FAKE_GCC
#include <armadillo>
#include <sys/stat.h>
#include "hdf5_interface.h"

// Data Types
enum TrajDim { THREED, TWOD };
enum TrajType { CARTESIAN, NONCARTESIAN, THREEDNONCARTESIAN};
enum SmsType { SMSoff, SMSon};

class MRI_DATA{
	public:
		
		//
		//	Variables
		// 
		
		// Sizes for shorthand
		int Num_Encodings;
		int Num_Readouts;
		int Num_Slices;
		int Num_Pts;
		int Num_Coils;
		
		// Special code required for simultaneous multi-slice 
		int sms_factor;
		NDarray::Array< NDarray::Array<float,3>,2> z;	  // Z coordinate for multi-slice with overlapped slices
				
		// Non-Cartesian Trajectory
		NDarray::Array< NDarray::Array<float,3>,1> kx;	// Fov = 1 unit, delta k =1
		NDarray::Array< NDarray::Array<float,3>,1> ky;
		NDarray::Array< NDarray::Array<float,3>,1> kz;
		NDarray::Array< NDarray::Array<float,3>,1> kw;
		NDarray::Array< NDarray::Array<float,3>,1> kt;	  // TE Time (s)
		NDarray::Array< NDarray::Array<std::complex<float>,3>,2> kdata;
		
		// Data for Noise samples 
		NDarray::Array< std::complex<float>,2> noise_samples; // data for noise samples
				
		//Physiologic Data for gating 
		NDarray::Array< double,3>ecg;	// Distance from ECG in MS
		NDarray::Array< double,3>resp;	// Respiratory signal from bellows or navigator
		NDarray::Array< double,3>time;	// Acquisition Time 
		NDarray::Array< double,3>prep;	// Time since a prep event (for example inversion)
		NDarray::Array< complex<float>,5>kdata_gating;	// Repeated sample for gating, need to be the same for each data point, all coils
						
		// Native Resolution
		int xres;
		int yres;
		int zres;
		int tres; 
		
		// Image Size in each dimension
		float zfov;
		float yfov;
		float xfov;
		
		// 2D/3D Cartesian/Non-Cartesian
		TrajDim trajectory_dims;
		TrajType trajectory_type;
		SmsType sms_type;
		
		
		//
		//	Functions
		// 
		
				
		//Temp
		char gate_name[1024];
		
		// Data Operations (move?)
		void undersample(int);
		void coilcompress(float);
		void whiten();
		void demod_kdata( float);
		
		// Initialization Filling Operations				
		void clone_attributes( MRI_DATA &);
		void init_memory();
		void init_gating_kdata(int);
		void init_noise_samples(int);
		
		void write_external_data(const char *fname);
		void read_external_data( const char *fname);
		
		void load_pcvipr_gating_file(const char *full_filename); //Temp
		void stats(void);
		void dump_stats(const std::string, const NDarray::Array< NDarray::Array<float,3>,1> & in);
		void dump_stats(const std::string, const NDarray::Array< NDarray::Array<complex<float>,3>,2> & in);
		
		void scale_fov( float,float,float);
		
		MRI_DATA subframe( int,int,int);
					
		MRI_DATA( MRI_DATA *);
		MRI_DATA( void );
	private:	
		
};

