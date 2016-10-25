#pragma once

#include <iostream>
 
#include <armadillo>
#include <cmath>
#include <omp.h>
#include "tictoc.hpp"
#include "ArrayTemplates.hpp"


class L2REG {
  public: 
  
    L2REG();
    L2REG(int numarg, char **pstring);
	
  	enum TransformType { NONE, TV, PHASE, LOWRES};
	static void help_message(void);
  	float lambda;
	float reg_scale;
	TransformType l2_type;
	NDarray::Array< complex<float>,3> ZeroPad;
	int verbose;		
	void regularize( NDarray::Array<complex<float>,3>&,NDarray::Array<complex<float>,3>&);
	void set_scale(float, NDarray::Array< NDarray::Array<complex<float>,3>,2>&);
	void set_scale(float, NDarray::Array<complex<float>,3>&);
	
};

