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
        enum WeightType { ITERATIVE,NON_ITERATIVE};
		enum GateType{ ECG,RESP,TIME,PREP,RESP_GATE}; 
				
		GATING(int numarg,char **pstring);
		void init( const MRI_DATA &data,int);
		
		// Filter Parameters
        int wdth_low;	// k=0 width
		int wdth_high; 	// k=kmax width
		float kmax;
		TornadoType tornado_shape; // kr^2 vs kr                         
        
		// Scaling
		float scale_time;
		float offset_time;
		
		//Control Gating Method		
		ViewshareType vs_type;	
		GateType gate_type;
		Array< float, 3>gate_times;
		
		// Frame Centers	
		float *gate_frames;
		
		int correct_resp_drift;
		
		// Function Calls		
		static void help_message(void);
        void weight_data(Array<float,3>&Tw, int e, const Array<float,3> &kx, const Array<float,3> &ky,const Array<float,3> &kz,int t,WeightType);

	   	void hist_weight( Array<float,3>&Tw, int e, int t);
		void tornado_weight(Array<float,3>&Tw, int e, const Array<float,3> &kx, const Array<float,3> &ky,const Array<float,3> &kz,int t,WeightType);
		void filter_resp(  const MRI_DATA &data );
    private:

};

