#pragma once 

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <complex>
#include <omp.h>

#define ARMA_FAKE_GCC
#include <armadillo>
#include "ArrayTemplates.hpp"
#include "gridFFT.h"
#include "mri_data.h"
#include "tictoc.hpp"

// View sharing modes
#define VS_NONE 0
#define VS_SLIDING 1
#define VS_TORNADO 2


class GATING {

    public:
        
		enum ViewshareType { TORNADO, NONE, HIST_MODE };
		enum TornadoType { FLAT, RADIAL, VIPR};
		enum WeightType { ITERATIVE,NON_ITERATIVE};
		enum FrameType { COMPOSITE, TIME_FRAME};
		enum GateType{ GATE_NONE,RETRO_ECG,ECG,RESP,TIME,PREP}; 
		enum RespGateType{RESP_NONE,RESP_THRESH,RESP_WEIGHT};
		
		
		GATING();						
		GATING(int numarg,char **pstring);
		void init( const MRI_DATA &data,int *);
		void init_resp_gating(const MRI_DATA &data);
		void init_time_resolved(const MRI_DATA &data,int *);
		
		// Tornado Filter Parameters
		int wdth_low;	// k=0 width
		int wdth_high; 	// k=kmax width
		float kmax;
		TornadoType tornado_shape; // kr^2 vs kr                         
        
		// Scaling for Waveform
		double scale_time;
		double offset_time;
		double actual_temporal_resolution;
		
		//Control Gating Method		
		ViewshareType vs_type;	
		GateType gate_type;
		
		NDarray::Array< double, 3>gate_times;
		NDarray::Array< float, 3>resp_weight;

		// Control of Retrospective Respiratory Gating
		RespGateType resp_gate_type;
		int correct_resp_drift;
		float resp_gate_efficiency;

		// Respiratory Signal 
		enum RespGateSignal{BELLOWS,DC_DATA};
		RespGateSignal resp_gate_signal;
		float resp_sign; 

		// Frame Centers	
		double *gate_frames;
		
		// Function Calls		
		static void help_message(void);
		void weight_data(NDarray::Array<float,3>&Tw, int e, const NDarray::Array<float,3> &kx, const NDarray::Array<float,3> &ky,const NDarray::Array<float,3>&kz,int t,WeightType, FrameType );
		float temporal_resolution(void);
	   	void hist_weight( NDarray::Array<float,3>&Tw, int e, int t);
		void tornado_weight(NDarray::Array<float,3>&Tw, int e, const NDarray::Array<float,3> &kx, const NDarray::Array<float,3> &ky,const NDarray::Array<float,3> &kz,int t,WeightType);
		void filter_resp(  const MRI_DATA &data );
		NDarray::Array< complex<float>,3> combine_kspace_channels(  const NDarray::Array< complex<float>,5> &kdata_gating );
    private:

};

