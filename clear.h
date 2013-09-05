#pragma once

#include <iostream>
#include <armadillo>
#include <cmath>
#include <omp.h>
#include "tictoc.hpp"
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
	LOWRANKCOIL();
    LOWRANKCOIL(int numarg, char **pstring);
	
	void update_threshold( NDarray::Array< NDarray::Array<complex<float>,3>, 3 > &image);
    void thresh( NDarray::Array< NDarray::Array<complex<float>,3>, 3 > &image,NDarray::Array< NDarray::Array<complex<float>,3>, 3 > &image2);
    void combine( NDarray::Array< NDarray::Array<complex<float>,3>, 3 > &image, NDarray::Array< NDarray::Array<complex<float>,3>, 2 > &out);
};

