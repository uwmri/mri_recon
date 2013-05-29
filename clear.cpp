#include "clear.h"
#include "io_templates.cpp"

using arma::cx_fmat;
using arma::fmat;
using arma::fvec;
using arma::uvec;
using namespace std;


// ----------------------
// Help Message
// ----------------------
void LOWRANKCOIL::help_message(void){
	cout << "----------------------------------------------" << endl;
	cout << "  Control for CLEAR " << endl;
	cout << "----------------------------------------------" << endl;
	
	help_flag("-clear_nx []","kernel size in x");
	help_flag("-clear_ny []","kernel size in y");
	help_flag("-clear_ns []","kernel size in z");
	help_flag("-clear_alpha[]","regularization factor");
	
}

LOWRANKCOIL::LOWRANKCOIL(int numarg, char **pstring){
  
  // Kernel Size Defaults 
  block_size_x = 4;
  block_size_y = 4;
  block_size_z = 4;

  clear_alpha = 0.01;
  debug =0;
  smax = 1.0;

#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);

  for(int pos=0; pos < numarg; pos++){
  
  	if(strcmp("-h", pstring[pos] ) == 0) {
	  	
		int_flag("-clear_nx",block_size_x);
		int_flag("-clear_ny",block_size_y);
		int_flag("-clear_nz",block_size_z);
		float_flag("-clear_alpha",clear_alpha);
	}
  }
}    

//-----------------------------------------------------
//Clear Threshold Estimate
//
//  Just estimate the max singular value over the matrix  
//-----------------------------------------------------
void LOWRANKCOIL::update_threshold( Array< complex<float>,4 > &image){

	float *stemp = new float[image.extent(thirdDim)];
		
	#pragma omp parallel for
	for( int k=0; k < image.extent(thirdDim); k+=block_size_z){
		stemp[k] = 0.0;
		
		cx_fmat A;
		A.zeros(block_size_x*block_size_y*block_size_z,image.extent(fourthDim)); // Pixels x Coils
	
		fmat S;
		S.zeros(block_size_x*block_size_y*block_size_z,image.extent(fourthDim)); // Pixels x Coils
		
		for( int j=0; j < image.extent(secondDim); j+=block_size_y){
		for( int i=0; i < image.extent(firstDim); i+=block_size_x){
		
			int count = 0;
			for(int kk=k; kk < k+block_size_z; kk++){
			for(int jj=j; jj < j+block_size_y; jj++){
			for(int ii=i; ii < i+block_size_x; ii++){
				for(int coil=0; coil < image.extent(fourthDim); coil++){
					A(count,coil) = image(ii,jj,kk,coil);
				}
				count++;
			}}}
			
			// SVD
			fvec s;
  			arma::svd(s,A);
			
			stemp[k]= ( s(0) > stemp[k] ) ? ( s(0) ) : ( stemp[k] );
	}}}
	
	smax = stemp[0];
	for( int k=0; k < image.extent(thirdDim); k+=block_size_z){
		smax = ( smax > stemp[k] ) ? ( smax ) : ( stemp[k] );
	}
	
	delete [] stemp;

	cout << "Max singular value is " << smax << endl;

}

//-----------------------------------------------------
//  Clear Threshold Call 
//-----------------------------------------------------
void LOWRANKCOIL::thresh( Array< complex<float>,4 > &image,Array< complex<float>,4 > &temp){
	
	
	if(clear_alpha==0.0){
		return;
	}	
	
	// Copy to temp location
	temp = image;
	image = complex<float>(0.0,0.0);
	
	int Ncoils = image.extent(fourthDim);
	
	int shiftx;
	int shifty;
	int shiftz;
	
	int rand_shiftx = rand() % block_size_x - block_size_x/2;
	int rand_shifty = rand() % block_size_y - block_size_y/2;
	int rand_shiftz = rand() % block_size_z - block_size_z/2;
		
	int Nx = image.extent(firstDim);
	int Ny = image.extent(secondDim);
	int Nz = image.extent(thirdDim);
	
	int step_size_x = (int)ceil((float)block_size_x);
	int step_size_y = (int)ceil((float)block_size_y);
	int step_size_z = (int)ceil((float)block_size_z);
		
	//cout << "Block size = " << block_size_x << " x " << block_size_y << " x " << block_size_z << endl;
	//cout << "Shifts = " << shiftx << " x " << shifty << " x " << shiftz << endl;
	//cout << "Steps = " << step_size_x << " x " << step_size_y << " x " << step_size_z << endl;
	
	for(int offset=0; offset <2; offset++){
		if(offset ==0){
			shiftx = -(int)floor(block_size_x/2)+rand_shiftx;
			shifty = -(int)floor(block_size_y/2)+rand_shifty;
			shiftz = -(int)floor(block_size_z/2)+rand_shiftz;
		}else{
			shiftx = rand_shiftx;
			shifty = rand_shifty;
			shiftz = rand_shiftz;
		}
		
		
			
	#pragma omp parallel for
	for( int k=0; k < image.extent(thirdDim); k+=step_size_z){
		
		cx_fmat A;
		A.zeros(block_size_x*block_size_y*block_size_z,image.extent(fourthDim)); // Pixels x Coils
	
		cx_fmat U;
		U.zeros(block_size_x*block_size_y*block_size_z,block_size_x*block_size_y*block_size_z); // Pixels x Pixels
	    
		fmat S;
		S.zeros(block_size_x*block_size_y*block_size_z,image.extent(fourthDim)); // Pixels x Coils
	    		
		cx_fmat V;
		V.zeros(image.extent(fourthDim),image.extent(fourthDim)); // Coils x Coils
	    
		 		
		for( int j=0; j < image.extent(secondDim); j+=step_size_y){ 
		for( int i=0; i < image.extent(firstDim); i+=step_size_x){
		
			int count = 0;
			for(int kk=(k+shiftz); kk < (k+block_size_z+shiftz); kk++){
			for(int jj=(j+shifty); jj < (j+block_size_y+shifty); jj++){
			for(int ii=(i+shiftx); ii < (i+block_size_x+shiftx); ii++){
				
				int px = (ii + Nx)% Nx;
				int py = (jj + Ny)% Ny;
				int pz = (kk + Nz)% Nz;
				for(int coil=0; coil < image.extent(fourthDim); coil++){
					A(count,coil) = temp(px,py,pz,coil);
				}
				count++;
			}}}
			
			// SVD
			fvec s;
  			arma::svd(U,s,V,A);
								
			S.zeros();
			for(int pos =0; pos< Ncoils; pos++){
				S(pos,pos)=   max( s(pos) - smax*clear_alpha, 0.0f );
			}
												
			// Reconstruct 
			cx_fmat X = U*S*trans(V);
			
			// Copy Back in weighted Fashion
			count = 0;
			for(int kk=(k+shiftz); kk < (k+block_size_z+shiftz); kk++){
			for(int jj=(j+shifty); jj < (j+block_size_y+shifty); jj++){
			for(int ii=(i+shiftx); ii < (i+block_size_x+shiftx); ii++){
				int px = (ii + Nx)% Nx;
				int py = (jj + Ny)% Ny;
				int pz = (kk + Nz)% Nz;
				
				complex<float>w(0.5,0.0);
				for(int coil=0; coil < image.extent(fourthDim); coil++){
					image(px,py,pz,coil) += X(count,coil)*w;
				}
				count++;
			}}}

	}}}
	
	}
		  
		 		  
}

//-----------------------------------------------------
//  Clear Combine Call 
//-----------------------------------------------------
void LOWRANKCOIL::combine( Array< complex<float>,4 > &image,Array< complex<float>,3 > &temp){
	
	temp =complex<float>(0.0,0.0);
	
	int shiftx;
	int shifty;
	int shiftz;
	
	int rand_shiftx = 0.0;
	int rand_shifty = 0.0;
	int rand_shiftz = 0.0;
		
	int Nx = image.extent(firstDim);
	int Ny = image.extent(secondDim);
	int Nz = image.extent(thirdDim);
	
	int step_size_x = (int)ceil((float)block_size_x);
	int step_size_y = (int)ceil((float)block_size_y);
	int step_size_z = (int)ceil((float)block_size_z);
		
	//cout << "Block size = " << block_size_x << " x " << block_size_y << " x " << block_size_z << endl;
	//cout << "Shifts = " << shiftx << " x " << shifty << " x " << shiftz << endl;
	//cout << "Steps = " << step_size_x << " x " << step_size_y << " x " << step_size_z << endl;
	
	for(int offset=0; offset <1; offset++){
		if(offset ==0){
			shiftx = -(int)floor(block_size_x/2)+rand_shiftx;
			shifty = -(int)floor(block_size_y/2)+rand_shifty;
			shiftz = -(int)floor(block_size_z/2)+rand_shiftz;
		}else{
			shiftx = rand_shiftx;
			shifty = rand_shifty;
			shiftz = rand_shiftz;
		}
		
		
			
	#pragma omp parallel for
	for( int k=0; k < image.extent(thirdDim); k+=step_size_z){
		
		cx_fmat A;
		A.zeros(block_size_x*block_size_y*block_size_z,image.extent(fourthDim)); // Pixels x Coils
	
		cx_fmat U;
		U.zeros(block_size_x*block_size_y*block_size_z,block_size_x*block_size_y*block_size_z); // Pixels x Pixels
	    
		fmat S;
		S.zeros(block_size_x*block_size_y*block_size_z,image.extent(fourthDim)); // Pixels x Coils
	    		
		cx_fmat V;
		V.zeros(image.extent(fourthDim),image.extent(fourthDim)); // Coils x Coils
	    
		 		
		for( int j=0; j < image.extent(secondDim); j+=step_size_y){ 
		for( int i=0; i < image.extent(firstDim); i+=step_size_x){
		
			int count = 0;
			for(int kk=(k+shiftz); kk < (k+block_size_z+shiftz); kk++){
			for(int jj=(j+shifty); jj < (j+block_size_y+shifty); jj++){
			for(int ii=(i+shiftx); ii < (i+block_size_x+shiftx); ii++){
				
				int px = (ii + Nx)% Nx;
				int py = (jj + Ny)% Ny;
				int pz = (kk + Nz)% Nz;
				for(int coil=0; coil < image.extent(fourthDim); coil++){
					A(count,coil) = image(px,py,pz,coil);
				}
				count++;
			}}}
			
			// SVD
			fvec s;
  			arma::svd(U,s,V,A);
								
			// Grab Largest Value
			cx_fmat A2 = A*V.cols(0,0);
						
			// Copy Back in weighted Fashion
			count = 0;
			for(int kk=(k+shiftz); kk < (k+block_size_z+shiftz); kk++){
			for(int jj=(j+shifty); jj < (j+block_size_y+shifty); jj++){
			for(int ii=(i+shiftx); ii < (i+block_size_x+shiftx); ii++){
				int px = (ii + Nx)% Nx;
				int py = (jj + Ny)% Ny;
				int pz = (kk + Nz)% Nz;
				
				complex<float>w(1.0,0.0);
				temp(px,py,pz) += A2(count,0)*w;
				count++;
			}}}

	}}}
	
	}
		  
		 		  
}


