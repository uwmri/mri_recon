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
#include "io_templates.cpp"
using namespace std;

// Thresholding methods
enum { TH_NONE,TH_FRACTION,TH_VISU,TH_BAYES,TH_SURE};

class THRESHOLD{
        public:
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
                void get_threshold(Array< complex<float>,5>&Coef);
                void get_visuthreshold(Array< complex<float>,5>&Coef);
                void get_bayesthreshold(Array< complex<float>,5>&Coef);
                void get_surethreshold(Array< complex<float>,5>&Coef);

		// Chooses the thresholding method
		void exec_threshold(Array< complex<float>,5>&Coef);

		// Execute thresholding method
		void exec_visuthreshold(Array< complex<float>,5>&Coef);
                void exec_multibandthreshold(Array< complex<float>,5>&Coef);
                void exec_fractionthreshold(Array<complex<float>,5>&Coef);
                void thresholding(Array< complex<float>,5>&Coef);

		void setThresholdMethod(int ith_type);
		int getThresholdMethod();
		
		void setTemporalThresholding(bool flag);	// to skip the first frame in processing
		bool getTemporalThresholding();

                // Iterative soft thresholding code
                void fista_update(Array< complex<float>,5>&X,Array< complex<float>,5 >&X_old,int iteration);
                
        private:        
                
};


#endif
