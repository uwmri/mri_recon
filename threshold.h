#ifndef hTHRESHLIB
#define hTHRESHLIB

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string.h>
#include <complex>
#include <omp.h>

#include "ArrayTemplates.hpp"
#include "wavelet3D.h" 

// Thresholding methods
enum { TH_NONE,TH_FRACTION,TH_VISU,TH_BAYES,TH_SURE};

class THRESHOLD{
        public:
        
		THRESHOLD();
		THRESHOLD(int numarg, char **pstring);
        bool soft;				// soft/hard thresholding
        bool thapp;				// threshold approx band
        bool temporal;
        double thresh;
        int waveL;
		bool VERBOSE;
        float global_threshold;
        NDarray::Array<float,5> subband_threshold;
		
		float noise;
		float noise_scale;
		int threshold_type;
	  	char th_type[20];
	  	char th_mode[20];
		static void help_message(void);
		
		
        // Get threshold value from coeffitients
        void update_threshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef, WAVELET3D &wave, float);
				
		// Chooses the thresholding method
		void exec_threshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef, WAVELET3D &wave);

		
		void setThresholdMethod(int ith_type);
		int getThresholdMethod();
		
		void setTemporalThresholding(bool flag);	// to skip the first frame in processing
		bool getTemporalThresholding();

        // Iterative soft thresholding code
        void fista_update(NDarray::Array<NDarray::Array< complex<float>,3>,2>&X,NDarray::Array< NDarray::Array< complex<float>,3>,2 >&X_old,int iteration);
                
        private:        
        
		// Execute thresholding method
		void exec_multibandthreshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef, WAVELET3D &wave );
        void thresholding(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef, float value);
        void thresholding( NDarray::Array< complex<float>,3>&Coef, float value);

		// Get Threshold
		float get_threshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef, float);
		float get_threshold(NDarray::Array< complex<float>,3> &Coef, float);
		float sure_cost( NDarray::Array< complex<float>,3>&Coef, float thresh);
		
        void get_visuthreshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef, WAVELET3D &wave);
        void get_bayesthreshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef, WAVELET3D &wave);
		void get_surethreshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef, WAVELET3D &wave);
		float get_bayesthreshold_subband( NDarray::Array< complex<float>,3> &XX);
		float get_surethreshold_subband( NDarray::Array< complex<float>,3> &XX);
       
		void robust_noise_estimate(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef, WAVELET3D &wave);

		        
};


#endif
