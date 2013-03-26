#pragma once

#include <iostream>
#include <armadillo>
#include <math.h>
#include <fftw3.h>
#include <omp.h>
#include "tictoc.cpp"
#include "ArrayTemplates.hpp"


class LOWRANKCOIL {
  public: 
  
  	// Size of Block
	int block_size_x;
	int block_size_y;
	int block_size_z;
	
	// Max Singular Value
	float smax;
		
	// Regularization Factor
	float clear_alpha;
	
	// Control 
	int debug;
	
	// Help Message
	static void help_message(void);
	
	// Functions
    LOWRANKCOIL(int numarg, char **pstring);
	void update_threshold( Array< complex<float>,4 > &image);
    void thresh( Array< complex<float>,4 > &image,Array< complex<float>,4 > &image);
    void combine( Array< complex<float>,4 > &image, Array< complex<float>,3> &out);
};

