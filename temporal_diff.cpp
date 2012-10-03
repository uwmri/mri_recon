/************************************************
3D Wavlet Libraries for pcvipr.e

Initial Author: Kevin M. Johnson
Description: This code contains functions for 3D wavelet transform. I tried blitzwave 
  but ended up not liking the lack of updates for newer gcc compilers and lack of threads.

*************************************************/

#include "temporal_diff.h"


TDIFF::TDIFF( array5D< complex<float> >temp){
	
	Nt = temp.Nt;
	Ne = temp.Ne;
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


void TDIFF::tdiff( array5D< complex<float> >temp){
	
	for(int e=0; e<temp.Ne; e++){
	#pragma omp parallel for 
	for(int k=0; k<temp.Nz; k++){
		cx_vec s(Nt);
		cx_vec ss(Nt);
		for(int j=0; j<temp.Ny; j++){
			for(int i=0; i<temp.Nx; i++){	
				//Copy
				for(int t=0; t<temp.Nt; t++){
					s(t)=temp[e][t][k][j][i];
				}
				
				// Transform
				ss=At*s;
				
				//Copy Back
				for(int t=0; t<temp.Nt; t++){
					temp[e][t][k][j][i] =ss(t);
				}
		}
	}}
	}
}

void TDIFF::inv_tdiff( array5D< complex<float> >temp){
	
	for(int e=0; e<temp.Ne; e++){
	#pragma omp parallel for 
	for(int k=0; k<temp.Nz; k++){
		cx_vec s(Nt);
		cx_vec ss(Nt);
		for(int j=0; j<temp.Ny; j++){
			for(int i=0; i<temp.Nx; i++){	
				//Copy
				for(int t=0; t<temp.Nt; t++){
					s(t)=temp[e][t][k][j][i];
				}
				
				// Transform
				ss=AIt*s;
				
				//Copy Back
				for(int t=0; t<temp.Nt; t++){
					temp[e][t][k][j][i] =ss(t);
				}
		}
	}}
	}
}


void TDIFF::forward( array5D< complex<float> >temp){
	for(int t=0; t<temp.Nt; t++){
	#pragma omp parallel for 
	for(int k=0; k<temp.Nz; k++){
		cx_vec s(Ne);
		cx_vec ss(Ne);
		for(int j=0; j<temp.Ny; j++){
			for(int i=0; i<temp.Nx; i++){	
				//Copy
				for(int e=0; e<temp.Ne; e++){
					s(e)=temp[e][t][k][j][i];
				}
				
				// Transform
				ss=Ae*s;
				
				//Copy Back
				for(int e=0; e<temp.Ne; e++){
					temp[e][t][k][j][i] =ss(e);
				}
		}
	}}
	}
}

void TDIFF::fft_t(array5D< complex<float> >temp){

	fftwf_plan p;
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(1);
	
	complex<float> *in = new complex<float>[temp.Nt];
    p = fftwf_plan_dft_1d(temp.Nt, (fftwf_complex*)in, (fftwf_complex*)in, FFTW_FORWARD, FFTW_MEASURE);
    
	float fft_scale = 1/ sqrt(temp.Nt);
	
	for(int e=0; e<temp.Ne; e++){
	for(int k=0; k<temp.Nz; k++){
		for(int j=0; j<temp.Ny; j++){
			for(int i=0; i<temp.Nx; i++){	
				//Copy
				for(int t=0; t<temp.Nt; t++){
					in[t]=temp[e][t][k][j][i];
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int t=0; t<temp.Nt; t++){
					temp[e][t][k][j][i] =( fft_scale*in[t]);
				}
		}
	}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
} 
	 

void TDIFF::ifft_t(array5D< complex<float> >temp){

	fftwf_plan p;
    
	complex<float> *in = new complex<float>[temp.Nt];
    
	fftwf_init_threads();
    fftwf_plan_with_nthreads(1);
	p = fftwf_plan_dft_1d(temp.Nt, (fftwf_complex*)in, (fftwf_complex*)in, FFTW_BACKWARD, FFTW_MEASURE);
    
	float fft_scale = 1/ sqrt(temp.Nt);
	
	for(int e=0; e<temp.Ne; e++){
	
	for(int k=0; k<temp.Nz; k++){
		for(int j=0; j<temp.Ny; j++){
			for(int i=0; i<temp.Nx; i++){	
				//Copy
				for(int t=0; t<temp.Nt; t++){
					in[t]=temp[e][t][k][j][i];
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int t=0; t<temp.Nt; t++){
					temp[e][t][k][j][i] =( fft_scale*in[t]);
				}
		}
	}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
} 
	
void TDIFF::fft_e(array5D< complex<float> >temp){

	fftwf_plan p;
    
	complex<float> *in = new complex<float>[temp.Ne];
    p = fftwf_plan_dft_1d(temp.Ne, (fftwf_complex*)in, (fftwf_complex*)in, FFTW_FORWARD, FFTW_ESTIMATE);
    
	float fft_scale = 1/ sqrt(temp.Ne);
	
	for(int t=0; t<temp.Nt; t++){
	for(int k=0; k<temp.Nz; k++){
		for(int j=0; j<temp.Ny; j++){
			for(int i=0; i<temp.Nx; i++){	
				//Copy
				for(int e=0; e<temp.Ne; e++){
					in[e]=temp[e][t][k][j][i];
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int e=0; e<temp.Ne; e++){
					temp[e][t][k][j][i] =( fft_scale*in[e]);
				}
		}
	}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
} 
	 

void TDIFF::ifft_e(array5D< complex<float> >temp){

	fftwf_plan p;
    
	complex<float> *in = new complex<float>[temp.Ne];
    p = fftwf_plan_dft_1d(temp.Ne, (fftwf_complex*)in, (fftwf_complex*)in, FFTW_BACKWARD, FFTW_ESTIMATE);
    float fft_scale = 1/ sqrt(temp.Ne);
	
	for(int t=0; t<temp.Nt; t++){
	for(int k=0; k<temp.Nz; k++){
		for(int j=0; j<temp.Ny; j++){
			for(int i=0; i<temp.Nx; i++){	
				//Copy
				for(int e=0; e<temp.Ne; e++){
					in[e]=temp[e][t][k][j][i];
				}
				fftwf_execute(p);
								
				//Copy Back
				for(int e=0; e<temp.Ne; e++){
					temp[e][t][k][j][i] =(fft_scale*in[e]);
				}
	}}}
	}
	
	fftwf_destroy_plan(p);
	delete [] in;
}

void TDIFF::backward( array5D< complex<float> >temp){
	
	for(int t=0; t<temp.Nt; t++){
	#pragma omp parallel for 
	for(int k=0; k<temp.Nz; k++){
		cx_vec s(Ne);
		cx_vec ss(Ne);
		for(int j=0; j<temp.Ny; j++){
			for(int i=0; i<temp.Nx; i++){	
				//Copy
				for(int e=0; e<temp.Ne; e++){
					s(e)=temp[e][t][k][j][i];
				}
				
				// Transform
				ss=AIe*s;
				
				//Copy Back
				for(int e=0; e<temp.Ne; e++){
					temp[e][t][k][j][i] =ss(e);
				}
		}
	}}
	}
}


