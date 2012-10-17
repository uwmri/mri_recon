/************************************************
3D Wavlet Libraries for pcvipr.e

Initial Author: Kevin M. Johnson
Description: This code contains functions for 3D wavelet transform. I tried blitzwave 
  but ended up not liking the lack of updates for newer gcc compilers and lack of threads.

*************************************************/

#include "temporal_diff.h"


TDIFF::TDIFF( Array< complex<float>,5>&temp){	 
	
	Nt = temp.length(fourthDim);
	Ne = temp.length(fifthDim);
	At.zeros(Nt,Nt);
	Ae.zeros(Ne,Ne);
	
	// Construct Temporal Differences Matrix
	for(int i = 0; i < Ne; i++) {
    	Ae(0,i)=1.0 / sqrt(Ne);
	}
	for(int i = 1; i < Ne; i++) {
    	Ae(i,i-1)=1.0/sqrt(2.0);
		Ae(i,i)=-1.0/sqrt(2.0);
	}
	AIe = Ae.i();
	
	
	// Construct Temporal Differences Matrix
	for(int i = 0; i < Nt; i++) {
    	At(0,i)=1.0 / sqrt(Nt);
	}
	for(int i = 1; i < Nt; i++) {
    	At(i,i-1)=1.0/sqrt(2.0);
		At(i,i)=-1.0/sqrt(2.0);
	}
	AIt = At.i();
}


void TDIFF::tdiff( Array< complex<float>,5>&temp){	 
	
	for(int e=0; e<temp.extent(fifthDim); e++){
	#pragma omp parallel for 
	for(int k=0; k<temp.extent(thirdDim); k++){
		cx_vec s(Nt);
		cx_vec ss(Nt);
		for(int j=0; j<temp.extent(secondDim); j++){
			for(int i=0; i<temp.extent(firstDim); i++){	
				//Copy
				for(int t=0; t<temp.extent(fourthDim); t++){
					s(t)=temp(i,j,k,t,e);
				}
				
				// Transform
				ss=At*s;
				
				//Copy Back
				for(int t=0; t<temp.extent(fourthDim); t++){
					temp(i,j,k,t,e) =ss(t);
				}
		}
	}}
	}
}

void TDIFF::inv_tdiff( Array< complex<float>,5>&temp){	 
	
	for(int e=0; e<temp.extent(fifthDim); e++){
	#pragma omp parallel for 
	for(int k=0; k<temp.extent(thirdDim); k++){
		cx_vec s(Nt);
		cx_vec ss(Nt);
		for(int j=0; j<temp.extent(secondDim); j++){
			for(int i=0; i<temp.extent(firstDim); i++){	
				//Copy
				for(int t=0; t<temp.extent(fourthDim); t++){
					s(t)=temp(i,j,k,t,e);
				}
				
				// Transform
				ss=AIt*s;
				
				//Copy Back
				for(int t=0; t<temp.extent(fourthDim); t++){
					temp(i,j,k,t,e) =ss(t);
				}
		}
	}}
	}
}


void TDIFF::ediff( Array< complex<float>,5>&temp){	 
	
	for(int t=0; t<temp.extent(fourthDim); t++){
	#pragma omp parallel for 
	for(int k=0; k<temp.extent(thirdDim); k++){
		cx_vec s(Ne);
		cx_vec ss(Ne);
		for(int j=0; j<temp.extent(secondDim); j++){
			for(int i=0; i<temp.extent(firstDim); i++){	
				//Copy
				for(int e=0; e<temp.extent(fifthDim); e++){
					s(e)=temp(i,j,k,t,e);
				}
				
				// Transform
				ss=Ae*s;
				
				//Copy Back
				for(int e=0; e<temp.extent(fifthDim); e++){
					temp(i,j,k,t,e) =ss(e);
				}
		}
	}}
	}
}


void TDIFF::inv_ediff( Array< complex<float>,5>&temp){	 
	
	for(int t=0; t<temp.extent(fourthDim); t++){
	#pragma omp parallel for 
	for(int k=0; k<temp.extent(thirdDim); k++){
		cx_vec s(Ne);
		cx_vec ss(Ne);
		for(int j=0; j<temp.extent(secondDim); j++){
			for(int i=0; i<temp.extent(firstDim); i++){	
				//Copy
				for(int e=0; e<temp.extent(fifthDim); e++){
					s(e)=temp(i,j,k,t,e);
				}
				
				// Transform
				ss=AIe*s;
				
				//Copy Back
				for(int e=0; e<temp.extent(fifthDim); e++){
					temp(i,j,k,t,e) =ss(e);
				}
		}
	}}
	}
}

	 
void TDIFF::fft_t(Array< complex<float>,5>&temp){	 

	fftwf_plan p;
    
	complex<float> *in = new complex<float>[temp.length(fourthDim)];
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(1);
	p = fftwf_plan_dft_1d(temp.length(fourthDim), (fftwf_complex*)in, (fftwf_complex*)in, FFTW_FORWARD, FFTW_MEASURE);
    
	float fft_scale = 1/ sqrt(temp.length(fourthDim));
	
	for(int e=0; e<temp.extent(fifthDim); e++){
	for(int k=0; k<temp.extent(thirdDim); k++){
		for(int j=0; j<temp.extent(secondDim); j++){
			for(int i=0; i<temp.extent(firstDim); i++){	
				//Copy
				for(int t=0; t<temp.extent(fourthDim); t++){
					in[t]=temp(i,j,k,t,e);
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int t=0; t<temp.extent(fourthDim); t++){
					temp(i,j,k,t,e) =( fft_scale*in[t]);
				}
		}
	}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
} 

void TDIFF::ifft_t(Array< complex<float>,5>&temp){	 

	fftwf_plan p;
    
	complex<float> *in = new complex<float>[temp.length(fourthDim)];
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(1);
	p = fftwf_plan_dft_1d(temp.length(fourthDim), (fftwf_complex*)in, (fftwf_complex*)in, FFTW_BACKWARD, FFTW_MEASURE);
    
	float fft_scale = 1/ sqrt(temp.length(fourthDim));
	
	for(int e=0; e<temp.extent(fifthDim); e++){
	for(int k=0; k<temp.extent(thirdDim); k++){
		for(int j=0; j<temp.extent(secondDim); j++){
			for(int i=0; i<temp.extent(firstDim); i++){	
				//Copy
				for(int t=0; t<temp.extent(fourthDim); t++){
					in[t]=temp(i,j,k,t,e);
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int t=0; t<temp.extent(fourthDim); t++){
					temp(i,j,k,t,e) =( fft_scale*in[t]);
				}
		}
	}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
} 

void TDIFF::fft_e(Array< complex<float>,5>&temp){	 

	fftwf_plan p;
    
	complex<float> *in = new complex<float>[temp.length(fifthDim)];
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(1);
	p = fftwf_plan_dft_1d(temp.length(fifthDim), (fftwf_complex*)in, (fftwf_complex*)in, FFTW_FORWARD, FFTW_MEASURE);
    
	float fft_scale = 1/ sqrt(temp.length(fifthDim));
	
	for(int t=0; t<temp.extent(fourthDim); t++){
	for(int k=0; k<temp.extent(thirdDim); k++){
		for(int j=0; j<temp.extent(secondDim); j++){
			for(int i=0; i<temp.extent(firstDim); i++){	
				//Copy
				for(int e=0; e<temp.extent(fifthDim); e++){
					in[e]=temp(i,j,k,t,e);
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int e=0; e<temp.extent(fifthDim); e++){
					temp(i,j,k,t,e) =( fft_scale*in[e]);
				}
		}
	}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
} 

void TDIFF::ifft_e(Array< complex<float>,5>&temp){	 

	fftwf_plan p;
    
	complex<float> *in = new complex<float>[temp.length(fifthDim)];
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(1);
	p = fftwf_plan_dft_1d(temp.length(fifthDim), (fftwf_complex*)in, (fftwf_complex*)in, FFTW_BACKWARD, FFTW_MEASURE);
    
	float fft_scale = 1/ sqrt(temp.length(fifthDim));
	
	for(int t=0; t<temp.extent(fourthDim); t++){
	for(int k=0; k<temp.extent(thirdDim); k++){
		for(int j=0; j<temp.extent(secondDim); j++){
			for(int i=0; i<temp.extent(firstDim); i++){	
				//Copy
				for(int e=0; e<temp.extent(fifthDim); e++){
					in[e]=temp(i,j,k,t,e);
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int e=0; e<temp.extent(fifthDim); e++){
					temp(i,j,k,t,e) =( fft_scale*in[e]);
				}
		}
	}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
} 

