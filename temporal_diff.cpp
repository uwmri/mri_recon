/************************************************
3D Wavlet Libraries for pcvipr.e

Initial Author: Kevin M. Johnson
Description: This code controls transfroms in the non-spatial dimensions

*************************************************/

#include "temporal_diff.h"

using namespace NDarray;
using arma::cx_mat;
using arma::cx_vec;

TRANSFORMS::TRANSFORMS(){
	reinit_pca = true;
	pca_count = 0;
}
  

void TRANSFORMS::get_difference_transform( cx_mat & A, cx_mat & Ai, int N){
	
	A.zeros(N,N);
	Ai.zeros(N,N);
	
	for(int i = 0; i < N; i++) {
    	A(0,i)=1.0;
	}
	for(int i = 1; i < N; i++) {
    	A(i,i-1)=1.0;
		A(i,i)= -1.0;
	}
	
	Ai = A.i();
	
	return;
}

void TRANSFORMS::get_wavelet_transform( cx_mat & A, cx_mat & Ai, int N){
	
	A.zeros(N,N);
	Ai.zeros(N,N);
	
	int levels = (int)log2( (double)N );
	//cout << "Levels = " << levels << endl;
	
	for( int level = 0; level < levels; level++){
		
		// R is the center of the matrix
		int R = (int)( 0.5*(float)N / pow((float)2.0,(float)level)); 
		int span = (int)pow(2.0f,(float)level);
	
		for(int i =0; i < R; i++){
			
			int offset = 2*span*i;
			
			// Zero the rows
			for(int j = 0; j < N; j++){
				A(i,j)  = complex<float>(0.0,0.0);
				A(i+R,j)= complex<float>(0.0,0.0);
			}
			
			// Set the average
			for(int j=offset; j< offset+2*span; j++){
				A(i,j) = 1.0;
			}
					
			// Set the difference
			for(int j=offset; j< offset+span; j++){
				A(i+R,j) = 1.0;
				A(i+R,j+span) = -1.0;
			}
		}
		//A.print("A-wavelet");
	}
	
	// Normalize
	for(int i =0; i<N; i++){
		float temp = 0.0;
		for(int j=0; j<N; j++){
			temp += norm( A(i,j));
		}
		temp = sqrt(temp);
		
		for(int j=0; j<N; j++){
			A(i,j)/= temp;
		}
	}
	//A.print("A-noramlize wavelet");	
	
	
	Ai = A.i();
	
	return;
}




void TRANSFORMS::tdiff( Array< Array< complex<float>,3>,2>&temp){
	 if(temp.length(firstDim)==1){
	 	return;
	 }
	 	 
	 cx_mat A;
	 cx_mat Ai;
	 get_difference_transform(A,Ai,temp.length(firstDim));
	 multiply_in_time(temp,A);
	 
	 return;
}

void TRANSFORMS::inv_tdiff( Array< Array< complex<float>,3>,2>&temp){
	 if(temp.length(firstDim)==1){
	 	return;
	 }
	 
	 cx_mat A;
	 cx_mat Ai;
	 get_difference_transform(A,Ai,temp.length(firstDim) );
	 multiply_in_time(temp,Ai);
	 
	 return;
}


void TRANSFORMS::twave( Array< Array< complex<float>,3>,2>&temp){
	 if(temp.length(firstDim)==1){
	 	return;
	 }
	 
	 cx_mat A;
	 cx_mat Ai;
	 get_wavelet_transform(A,Ai,temp.length(firstDim));
	 multiply_in_time(temp,A);
	 
	 return;
}

void TRANSFORMS::inv_twave( Array< Array< complex<float>,3>,2>&temp){
	 if(temp.length(firstDim)==1){
	 	return;
	 }
	 
	 cx_mat A;
	 cx_mat Ai;
	 get_wavelet_transform(A,Ai,temp.length(firstDim));
	 multiply_in_time(temp,Ai);
	 
	 return;
}







void TRANSFORMS::ediff( Array< Array< complex<float>,3>,2>&temp){
	 if(temp.length(secondDim)==1){
	 	return;
	 }
	 
	 cx_mat A;
	 cx_mat Ai;
	 get_difference_transform(A,Ai,temp.length(secondDim));
	 multiply_in_encode(temp,A);
	 
	 return;
}

void TRANSFORMS::inv_ediff( Array< Array< complex<float>,3>,2>&temp){
	 if(temp.length(secondDim)==1){
	 	return;
	 }
	 cx_mat A;
	 cx_mat Ai;
	 get_difference_transform(A,Ai,temp.length(secondDim) );
	 multiply_in_encode(temp,Ai);
	 return;
}


void TRANSFORMS::ewave( Array< Array< complex<float>,3>,2>&temp){
	 if(temp.length(secondDim)==1){
	 	return;
	 }
	 cx_mat A;
	 cx_mat Ai;
	 get_wavelet_transform(A,Ai,temp.length(secondDim));
	 multiply_in_encode(temp,A);
	 return;
}

void TRANSFORMS::inv_ewave( Array< Array< complex<float>,3>,2>&temp){
	 if(temp.length(secondDim)==1){
	 	return;
	 }
	 cx_mat A;
	 cx_mat Ai;
	 get_wavelet_transform(A,Ai,temp.length(secondDim));
	 multiply_in_encode(temp,Ai);
	 return;
}

//-----------------------------------------------------
//  Clear Threshold Call 
//-----------------------------------------------------
void TRANSFORMS::eigen( Array< Array<complex<float>,3>,2 > &image, int dim, int direction){
		
		
	pca_count++;
		
	if(reinit_pca && (direction==FORWARD) ){
		if(pca_count > 5){
			reinit_pca = false;
		}	
		
		// Shorthand
		int Nt =image.extent(firstDim);
		int Ne =image.extent(secondDim);
		int Nx =image(0,0).extent(firstDim);
		int Ny =image(0,0).extent(secondDim);
		int Nz =image(0,0).extent(thirdDim);
	
		int block_size_x = 8;
		int block_size_y = 8;
		int block_size_z = 8;
        
        block_size_x = ( block_size_x > Nx ) ? ( Nx ) : ( block_size_x);
        block_size_y = ( block_size_y > Ny ) ? ( Ny ) : ( block_size_y);
        block_size_z = ( block_size_z > Nz ) ? ( Nz ) : ( block_size_z);
                
		int block_Nx= (int)( Nx / block_size_x );
		int block_Ny= (int)( Ny / block_size_y );
		int block_Nz= (int)( Nz / block_size_z );
		int total_blocks = block_Nx*block_Ny*block_Nz;
	
		/* Get Image Size*/
		int N = image.extent(dim);
		int Np = total_blocks; 
	
		cout << " Learning Eigen Mat (N=" << N << ")(Np = " << Np << ")" << endl;
	
		// For nested parallelism fix
		int count =0;
		arma::Col<int> act_i(total_blocks);
		arma::Col<int> act_j(total_blocks);
		arma::Col<int> act_k(total_blocks);
		for(int i=0; i < block_Nx;i++){
		for(int j=0; j < block_Ny;j++){
		for(int k=0; k < block_Nz;k++){
			act_i(count) = i;
			act_j(count) = j;
			act_k(count) = k;
			count++;
		}}}
		
		// Storage Block
		arma::cx_mat A;
		A.zeros(Np,N); 
	
		cout << "Collecting Blocks" << endl;	
		#pragma omp parallel for 
		for( int block = 0; block < total_blocks; block++){
		
			// Nested parallelism workaround
			int i = act_i(block)*block_size_x;
			int j = act_j(block)*block_size_y;
			int k = act_k(block)*block_size_z;
					
			//-----------------------------------------------------
			//   Block Gather
			//-----------------------------------------------------
			
			switch(dim){
				case(0):{
					for(int t=0; t < Nt; t++){
						for(int e=0; e < Ne; e++){
						for(int kk=(k); kk < (k+block_size_z); kk++){
						for(int jj=(j); jj < (j+block_size_y); jj++){
						for(int ii=(i); ii < (i+block_size_x); ii++){
							int px = (ii + Nx)% Nx;
							int py = (jj + Ny)% Ny;
							int pz = (kk + Nz)% Nz;
							A(block,t) += image(t,e)(px,py,pz);
						}}}}
					}
				}break;
				
				case(1):{
					for(int e=0; e < Ne; e++){
						for(int t=0; t < Nt; t++){
						for(int kk=(k); kk < (k+block_size_z); kk++){
						for(int jj=(j); jj < (j+block_size_y); jj++){
						for(int ii=(i); ii < (i+block_size_x); ii++){
							int px = (ii + Nx)% Nx;
							int py = (jj + Ny)% Ny;
							int pz = (kk + Nz)% Nz;
							A(block,e) += image(t,e)(px,py,pz);
						}}}}
					}
				}break;
			}// Switch Dim	
		}// Block Loop
	
	
		cout << "Learning PCA " << endl;
		arma::cx_mat U;
		arma::cx_mat V;
		arma::vec s;
  		arma::svd_econ(U,s,V,A);
		
		cx_mat wV = diagmat(s)*V.t();
		E = wV;
		
		Ei= wV.i();
	
	}/* Reinit PCA*/
	
		
	if(dim==0){
		if(direction==FORWARD){
			multiply_in_time( image,E);
		}else{
			multiply_in_time( image,Ei);
		}
	}else{
		if(direction==FORWARD){
			multiply_in_encode( image,E);
		}else{
			multiply_in_encode( image,Ei);
		}
	}
	
	if(direction == FORWARD){
		if( image.numElements() > 1){
		 	int count=0;
			for( Array< Array<complex<float>,3>,2>::iterator miter=image.begin(); miter!=image.end(); miter++){
								
					Array<complex<float>,2>Xf=(*miter)(Range::all(),Range::all(),image(0,0).length(2)/2);
					if(count==0){
						ArrayWriteMag(Xf,"Xtrans_frames.dat");
					}else{
						ArrayWriteMagAppend(Xf,"Xtrans_frames.dat");
					}
					count++;
			}
		}
	}
	
}


void TRANSFORMS::multiply_in_time( Array< Array< complex<float>,3>,2>&temp, cx_mat A){
	int Nt = temp.length(firstDim);
	int Ne = temp.length(secondDim);
	int Nx = temp(0).length(firstDim);
	int Ny = temp(0).length(secondDim);
	int Nz = temp(0).length(thirdDim);
	
	for(int e=0; e< Ne; e++){
	
	#pragma omp parallel for
	for(int k=0; k< Nz; k++){
		cx_vec s(Nt);
		cx_vec ss(Nt);
		for(int j=0; j< Ny; j++){
			for(int i=0; i<Nx; i++){	
				
				//Copy
				for(int t=0; t<Nt; t++){
					s(t)=temp(t,e)(i,j,k);
				}
				
				// Transform
				ss=A*s;
				
				//Copy Back
				for(int t=0; t<Nt; t++){
					temp(t,e)(i,j,k) =ss(t);
				}
		}
	}}
	}
	return;
}


void TRANSFORMS::multiply_in_encode( Array< Array< complex<float>,3>,2>&temp, cx_mat A){ 
	
	int Nt = temp.length(firstDim);
	int Ne = temp.length(secondDim);
	int Nx = temp(0).length(firstDim);
	int Ny = temp(0).length(secondDim);
	int Nz = temp(0).length(thirdDim);
	
	for(int t=0; t<Nt; t++){
	#pragma omp parallel for 
	for(int k=0; k< Nz; k++){
		cx_vec s(Ne);
		cx_vec ss(Ne);
		for(int j=0; j< Ny; j++){
			for(int i=0; i< Nx; i++){	
				//Copy
				for(int e=0; e< Ne; e++){
					s(e)=temp(t,e)(i,j,k);
				}
				
				// Transform
				ss=A*s;
				
				//Copy Back
				for(int e=0; e< Ne; e++){
					temp(t,e)(i,j,k) =ss(e);
				}
		}
	}}
	}
	return;
}

	 
void TRANSFORMS::fft_t(Array< Array< complex<float>,3>,2>&temp){ 
	
	int Nt = temp.length(firstDim);
	int Ne = temp.length(secondDim);
	int Nx = temp(0).length(firstDim);
	int Ny = temp(0).length(secondDim);
	int Nz = temp(0).length(thirdDim);
		
	fftwf_plan p;
    
	complex<float> *in = new complex<float>[Nt];
    
	//fftwf_init_threads();
    	//fftwf_plan_with_nthreads(1);
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
	
	return;
} 

void TRANSFORMS::ifft_t(Array< Array< complex<float>,3>,2>&temp){ 
	
	int Nt = temp.length(firstDim);
	int Ne = temp.length(secondDim);
	int Nx = temp(0).length(firstDim);
	int Ny = temp(0).length(secondDim);
	int Nz = temp(0).length(thirdDim);
		
	fftwf_plan p;
    
	complex<float> *in = new complex<float>[Nt];
    
	//fftwf_init_threads();
    	//fftwf_plan_with_nthreads(1);
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
	
	return;
} 

void TRANSFORMS::fft_e(Array< Array< complex<float>,3>,2>&temp){ 
	
	int Nt = temp.length(firstDim);
	int Ne = temp.length(secondDim);
	int Nx = temp(0).length(firstDim);
	int Ny = temp(0).length(secondDim);
	int Nz = temp(0).length(thirdDim);
		
	fftwf_plan p;
    
	complex<float> *in = new complex<float>[Ne];
    
	//fftwf_init_threads();
    	//fftwf_plan_with_nthreads(1);
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
	return;
} 



void TRANSFORMS::ifft_e(Array< Array< complex<float>,3>,2>&temp){ 
	
	int Nt = temp.length(firstDim);
	int Ne = temp.length(secondDim);
	int Nx = temp(0).length(firstDim);
	int Ny = temp(0).length(secondDim);
	int Nz = temp(0).length(thirdDim);
		
	fftwf_plan p;
    
	complex<float> *in = new complex<float>[Ne];
    
	//fftwf_init_threads();
    	//fftwf_plan_with_nthreads(1);
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
	return;
}
