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


