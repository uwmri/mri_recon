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
		int VERBOSE;
        float threshold;
        float noise;
		int threshold_type;
	  	char th_type[20];
	  	char th_mode[20];
		static void help_message(void);

        // Get threshold value from coeffitients
        void get_threshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef);
        void get_visuthreshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef);
        void get_bayesthreshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef);
        void get_surethreshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef);

		// Chooses the thresholding method
		void exec_threshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef);

		// Execute thresholding method
		void exec_visuthreshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef);
        void exec_multibandthreshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef);
        void exec_fractionthreshold(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef);
        void thresholding(NDarray::Array<NDarray::Array< complex<float>,3>,2>&Coef);

		void setThresholdMethod(int ith_type);
		int getThresholdMethod();
		
		void setTemporalThresholding(bool flag);	// to skip the first frame in processing
		bool getTemporalThresholding();

        // Iterative soft thresholding code
        void fista_update(NDarray::Array<NDarray::Array< complex<float>,3>,2>&X,NDarray::Array< NDarray::Array< complex<float>,3>,2 >&X_old,int iteration);
                
        private:        
                
};


#endif
