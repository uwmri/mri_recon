/************************************************
3D Wavlet Libraries for pcvipr.e

Initial Author: Kevin M. Johnson
Description: This code contains functions for 3D wavelet transform. I tried blitzwave 
  but ended up not liking the lack of updates for newer gcc compilers and lack of threads.

*************************************************/

#include "temporal_diff.h"


TDIFF::TDIFF(Array< Array< complex<float>,3>,2>&temp){ 
	
	Nt = temp.length(firstDim);
	Ne = temp.length(secondDim);
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


void TDIFF::tdiff( Array< Array< complex<float>,3>,2>&temp){
	
	
	for(int e=0; e<temp.extent(secondDim); e++){
	#pragma omp parallel for
	 
	for(int k=0; k<temp(0).extent(thirdDim); k++){
		cx_vec s(Nt);
		cx_vec ss(Nt);
		for(int j=0; j<temp(0).extent(secondDim); j++){
			for(int i=0; i<temp(0).extent(firstDim); i++){	
				//Copy
				for(int t=0; t<temp.extent(firstDim); t++){
					s(t)=temp(t,e)(i,j,k);
				}
				
				// Transform
				ss=At*s;
				
				//Copy Back
				for(int t=0; t<temp.extent(firstDim); t++){
					temp(t,e)(i,j,k) =ss(t);
				}
		}
	}}
	}
}

void TDIFF::inv_tdiff( Array< Array< complex<float>,3>,2>&temp){ 
	
	for(int e=0; e<temp.extent(secondDim); e++){
	#pragma omp parallel for
	 
	for(int k=0; k<temp(0).extent(thirdDim); k++){
		cx_vec s(Nt);
		cx_vec ss(Nt);
		for(int j=0; j<temp(0).extent(secondDim); j++){
			for(int i=0; i<temp(0).extent(firstDim); i++){	
				//Copy
				for(int t=0; t<temp.extent(firstDim); t++){
					s(t)=temp(t,e)(i,j,k);
				}
				
				// Transform
				ss=AIt*s;
				
				//Copy Back
				for(int t=0; t<temp.extent(firstDim); t++){
					temp(t,e)(i,j,k) =ss(t);
				}
		}
	}}
	}
}


void TDIFF::ediff( Array< Array< complex<float>,3>,2>&temp){ 
	
	for(int t=0; t<temp.extent(firstDim); t++){
	#pragma omp parallel for 
	for(int k=0; k<temp(0).extent(thirdDim); k++){
		cx_vec s(Ne);
		cx_vec ss(Ne);
		for(int j=0; j<temp(0).extent(secondDim); j++){
			for(int i=0; i<temp(0).extent(firstDim); i++){	
				//Copy
				for(int e=0; e<temp.extent(secondDim); e++){
					s(e)=temp(t,e)(i,j,k);
				}
				
				// Transform
				ss=Ae*s;
				
				//Copy Back
				for(int e=0; e<temp.extent(secondDim); e++){
					temp(t,e)(i,j,k) =ss(e);
				}
		}
	}}
	}
}


void TDIFF::inv_ediff( Array< Array< complex<float>,3>,2>&temp){ 
	
	for(int t=0; t<temp.extent(firstDim); t++){
	#pragma omp parallel for 
	for(int k=0; k<temp(0).extent(thirdDim); k++){
		cx_vec s(Ne);
		cx_vec ss(Ne);
		for(int j=0; j<temp(0).extent(secondDim); j++){
			for(int i=0; i<temp(0).extent(firstDim); i++){	
				//Copy
				for(int e=0; e<temp.extent(secondDim); e++){
					s(e)=temp(t,e)(i,j,k);
				}
				
				// Transform
				ss=AIe*s;
				
				//Copy Back
				for(int e=0; e<temp.extent(secondDim); e++){
					temp(t,e)(i,j,k) =ss(e);
				}
		}
	}}
	}
}


	 
void TDIFF::fft_t(Array< Array< complex<float>,3>,2>&temp){ 
	
	Nt = temp.length(firstDim);
	Ne = temp.length(secondDim);
	int Nx = temp(0).length(firstDim);
	int Ny = temp(0).length(secondDim);
	int Nz = temp(0).length(thirdDim);
		
	fftwf_plan p;
    
	complex<float> *in = new complex<float>[Nt];
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(1);
	p = fftwf_plan_dft_1d(Nt, (fftwf_complex*)in, (fftwf_complex*)in, FFTW_FORWARD, FFTW_MEASURE);
    
	float fft_scale = 1/ sqrt(Nt);
	
	for(int e=0; e< Ne; e++){
	for(int k=0; k< Nz; k++){
		for(int j=0; j< Ny; j++){
			for(int i=0; i<Nx; i++){	
				//Copy
				for(int t=0; t<Nt; t++){
					in[t]=temp(t,e)(i,j,k);
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int t=0; t<Nt; t++){
					temp(t,e)(i,j,k) =( fft_scale*in[t]);
				}
		}
	}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
} 

void TDIFF::ifft_t(Array< Array< complex<float>,3>,2>&temp){ 
	
	Nt = temp.length(firstDim);
	Ne = temp.length(secondDim);
	int Nx = temp(0).length(firstDim);
	int Ny = temp(0).length(secondDim);
	int Nz = temp(0).length(thirdDim);
		
	fftwf_plan p;
    
	complex<float> *in = new complex<float>[Nt];
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(1);
	p = fftwf_plan_dft_1d(Nt, (fftwf_complex*)in, (fftwf_complex*)in, FFTW_BACKWARD, FFTW_MEASURE);
    
	float fft_scale = 1/ sqrt(Nt);
	
	for(int e=0; e< Ne; e++){
	for(int k=0; k< Nz; k++){
		for(int j=0; j< Ny; j++){
			for(int i=0; i<Nx; i++){	
				//Copy
				for(int t=0; t<Nt; t++){
					in[t]=temp(t,e)(i,j,k);
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int t=0; t<Nt; t++){
					temp(t,e)(i,j,k) =( fft_scale*in[t]);
				}
		}
	}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
} 

void TDIFF::fft_e(Array< Array< complex<float>,3>,2>&temp){ 
	
	Nt = temp.length(firstDim);
	Ne = temp.length(secondDim);
	int Nx = temp(0).length(firstDim);
	int Ny = temp(0).length(secondDim);
	int Nz = temp(0).length(thirdDim);
		
	fftwf_plan p;
    
	complex<float> *in = new complex<float>[Ne];
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(1);
	p = fftwf_plan_dft_1d(Ne, (fftwf_complex*)in, (fftwf_complex*)in, FFTW_FORWARD, FFTW_MEASURE);
    
	float fft_scale = 1/ sqrt(Ne);
	
	for(int t=0; t<Nt; t++){
	for(int k=0; k< Nz; k++){
		for(int j=0; j< Ny; j++){
			for(int i=0; i<Nx; i++){	
				//Copy
				for(int e=0; e< Ne; e++){
					in[e]=temp(t,e)(i,j,k);
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int e=0; e< Ne; e++){
					temp(t,e)(i,j,k) =( fft_scale*in[e]);
				}
		}
	}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
} 



void TDIFF::ifft_e(Array< Array< complex<float>,3>,2>&temp){ 
	
	Nt = temp.length(firstDim);
	Ne = temp.length(secondDim);
	int Nx = temp(0).length(firstDim);
	int Ny = temp(0).length(secondDim);
	int Nz = temp(0).length(thirdDim);
		
	fftwf_plan p;
    
	complex<float> *in = new complex<float>[Ne];
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(1);
	p = fftwf_plan_dft_1d(Ne, (fftwf_complex*)in, (fftwf_complex*)in, FFTW_BACKWARD, FFTW_MEASURE);
    
	float fft_scale = 1/ sqrt(Ne);
	
	for(int t=0; t<Nt; t++){
	for(int k=0; k< Nz; k++){
		for(int j=0; j< Ny; j++){
			for(int i=0; i<Nx; i++){	
				//Copy
				for(int e=0; e< Ne; e++){
					in[e]=temp(t,e)(i,j,k);
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int e=0; e< Ne; e++){
					temp(t,e)(i,j,k) =( fft_scale*in[e]);
				}
		}
	}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
}
