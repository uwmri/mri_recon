#pragma once 

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <complex>
#include <omp.h>
#include <armadillo>
#include "wavelet3D.h"
#include "ArrayTemplates.hpp"
#include "mri_data.h"
#include "gridFFT.h"
#include "tictoc.cpp"
#include "io_templates.cpp"

// View sharing modes
#define VS_NONE 0
#define VS_SLIDING 1
#define VS_TORNADO 2

/*Don't define namespace in .h 
using namespace std;
*/

class GATING {

    public:
        
		enum ViewshareType { TORNADO, NONE, HIST_MODE };
		enum TornadoType { FLAT, RADIAL, VIPR};
        
		GATING(int numarg,char **pstring);
		void init( Array<float,4>&times,int,const Array<float,4>&,const Array<float,4>&,const Array<float,4>&);
		
		// Filter Parameters
        int wdth_low;	// k=0 width
		int wdth_high; 	// k=kmax width
		float kmax;
		TornadoType tornado_shape; // kr^2 vs kr                         
        
		// Scaling
		float scale_time;
		float offset_time;
		
		//Control Gating Method		
		ViewshareType type;	
		
		// Frame Centers	
		float *gate_frames;
		
		// Function Calls		
		static void help_message(void);
        void weight_data(Array<float,3>&Tw, Array<float,3>&times, const Array<float,3> &kx, const Array<float,3> &ky,const Array<float,3> &kz,int t);

	   	void hist_weight( Array<float,3>&Tw,Array<float,3>&times, int t);
		void tornado_weight(Array<float,3>&Tw, Array<float,3>&times, const Array<float,3> &kx, const Array<float,3> &ky,const Array<float,3> &kz,int t);

    private:

};

