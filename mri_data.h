#pragma once

// System Libraries
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstring>
#include <cmath>
#include <string>
#include <complex>
#include <omp.h>
#include <sys/stat.h>

// Matrix Library
 
#include <armadillo>

// Local Libraries
#include "ArrayTemplates.hpp"
#include "hdf5_interface.h"


class MRI_DATA{
	public:
				
		// 
		//	Data Descriptors
		// 
		enum TrajType { CARTESIAN, NONCARTESIAN};
		enum SmsType { SMSoff, SMSon};
		
		//
		//	Variables
		// 
		
		// This is the container size
		int Num_Frames;
		int Num_Encodings;
		int Num_Coils;
		
		// Special code required for simultaneous multi-slice 
		int sms_factor;
		NDarray::Array< NDarray::Array<float,3>,2> z;	  // Z coordinate for multi-slice with overlapped slices
		
		// Non-Cartesian Trajectory
		//  note FOV = 1, delta k=1 for unit
		NDarray::Array< NDarray::Array<float,3>,1> kx;	
		NDarray::Array< NDarray::Array<float,3>,1> ky;
		NDarray::Array< NDarray::Array<float,3>,1> kz;
		NDarray::Array< NDarray::Array<float,3>,1> kw;
		NDarray::Array< NDarray::Array<float,3>,1> kt;	 
		NDarray::Array< NDarray::Array<std::complex<float>,3>,2> kdata;
		
		// Data for Noise samples 
		NDarray::Array< std::complex<float>,2> noise_samples; // data for noise samples
				
		//Physiologic Data for gating 
		NDarray::Array< NDarray::Array<double,2>,1>ecg;	// Distance from ECG in MS
		NDarray::Array< NDarray::Array<double,2>,1>resp;	// Respiratory signal from bellows or navigator
		NDarray::Array< NDarray::Array<double,2>,1>time;	// Acquisition Time 
		NDarray::Array< NDarray::Array<double,2>,1>prep;	// Time since a prep event (for example inversion)
		NDarray::Array< NDarray::Array<complex<float>,2>,2>kdata_gating; // Repeated sample for gating, need to be the same for each data point, all coils
							
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
		NDarray::TinyVector<TrajType,3> trajectory_type;
		NDarray::TinyVector<bool,3> dft_needed;
		SmsType sms_type;
						
		//
		//	Functions
		// 
		
		// Constructors
		MRI_DATA subframe( int,int,int);
		MRI_DATA( MRI_DATA *);
		MRI_DATA( void );
		
		// Combiners 
		void convert_encodes_to_coils(void);
				
		// Data Operations (move?)
		void coilcompress(float);
		void whiten();
		void demod_kdata( float);
		void scale_fov( float,float,float);
		
		// Initialization Filling Operations				
		void clone_attributes( MRI_DATA &);
		void init_memory();
		void init_memory(int,int,int);
		void init_gating_kdata(int);
		void init_noise_samples(int);
		void init_encode(int,int,int,int);
					
		// HDF5 Data 
		void write_external_data(const char *fname);
		void read_external_data( const char *fname);
		
		// Print to Stdout
		void stats(void);
		void dump_stats(const std::string, const NDarray::Array< NDarray::Array<float,3>,1> & in);
		void dump_stats(const std::string, const NDarray::Array< NDarray::Array<complex<float>,3>,2> & in);
		
	private:	
		
};

