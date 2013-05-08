#pragma once

#include <iostream>
#include <armadillo>
#include <math.h>
#include <fftw3.h>
#include <omp.h>
#include "tictoc.cpp"
#include "ArrayTemplates.hpp"


class L2REG {
  public: 
  
    L2REG(int numarg, char **pstring);
	
  	enum TransformType { NONE, TV, PHASE};
	static void help_message(void);
  	float lambda;
	float reg_scale;
	TransformType l2_type;
	int verbose;		
	void regularize( Array<complex<float>,3>&,Array<complex<float>,3>&);
	void set_scale(float, Array<complex<float>,5>&);
	
};

