#include "l2reg.h"
#include "io_templates.hpp"

using arma::cx_fmat;
using arma::fmat;
using arma::fvec;
using arma::uvec;
using namespace std;
using namespace NDarray;


// ----------------------
// Help Message
// ----------------------
void L2REG::help_message(void){
	cout << "----------------------------------------------" << endl;
	cout << "  Control for L2 Regularization " << endl;
	cout << "----------------------------------------------" << endl;

	help_flag("-l2_lambda []","relative amount to regularize");
	help_flag("-l2_image","L2 of image");
	help_flag("-l2_tv","L2 of total variation");
	help_flag("-l2_phase","L2 of phase");
}

L2REG::L2REG(){
}

L2REG::L2REG(int numarg, char **pstring){
  
	verbose =0;
	lambda = 0.00;
	l2_type = NONE;
	reg_scale =0.0;
	
#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);

  for(int pos=0; pos < numarg; pos++){
  
  	if(strcmp("-h", pstring[pos] ) == 0) {
	  	float_flag("-l2_lambda",lambda);
		trig_flag(NONE,"-l2_image",l2_type);
		trig_flag(TV,"-l2_tv",l2_type);
		trig_flag(PHASE,"-l2_phase",l2_type);
		trig_flag(1,"-verbose",verbose);
		
	}
  }
}    

void L2REG::set_scale( float EhE_scale, Array< Array< complex<float>,3>,2>&X){
	
	if(lambda==0.0){
		return;
	}
	// Input should be energy of EhE X
	// Scales regularization by EhEx/TVx
	
	float transform_energy=0.0;
	switch(l2_type){
	
		case(NONE):{
	  
	  		#pragma omp parallel for reduction(+:transform_energy)
			for(int et=0; et < X.numElements(); et++){ 	  
				int t = et % X.length(firstDim);
				int e = et / X.length(firstDim);
								
				float temp = sum(norm(X(t,e)));
				transform_energy += temp; 
			}
		}break;
	
	
		case(TV):{
			
			for(Array< Array< complex<float>,3>,2>::iterator miter=X.begin(); miter !=X.end(); miter++){
			
	  			int Nx = (*miter).length(firstDim);
	  			int Ny = (*miter).length(secondDim);
	  			int Nz = (*miter).length(thirdDim);
	  			Array< complex<float>,3> Xref;
				Xref.reference(*miter);
	    
	  			for(int k =0; k < Xref.length(thirdDim);k++){
	  			for(int j =0; j < Xref.length(secondDim);j++){
	  			for(int i =0; i < Xref.length(firstDim);i++){
					complex<float>val = complex<float>(1.0/6.0,0.0)*( complex<float>(6.0,0.0)*Xref(i,j,k) - Xref((i+1)%Nx,j,k) - Xref((i+Nx-1)%Nx,j,k) - Xref(i,(j+1)%Ny,k) - Xref(i,(j+Ny-1)%Ny,k)  - Xref(i,j,(k+1)%Nz) - Xref(i,j,(k+Nz-1)%Nz));
					transform_energy += norm( val );
	  			}}}
	  	
	  		}	
	  
		}break;
	
		case(PHASE):{
	 		for(Array< Array< complex<float>,3>,2>::iterator miter=X.begin(); miter !=X.end(); miter++){
				transform_energy += sum( pow( imag(*miter),2));
	  		}
		}break;
	}
	cout << "Transform Energy = " << transform_energy << endl;
	cout << "EhE Energy = " << EhE_scale << endl;
		
	reg_scale = lambda*EhE_scale/transform_energy;
	
}



void L2REG::set_scale( float EhE_scale, Array< complex<float>,3> &X){
	
	if(lambda==0.0){
		return;
	}
	// Input should be energy of EhE X
	// Scales regularization by EhEx/TVx
	
	float transform_energy=0.0;
	switch(l2_type){
	
		case(NONE):{
	  		transform_energy = sum(norm(X)); 
			
		}break;
	
	
		case(TV):{
			
			int Nx = X.length(firstDim);
	  		int Ny = X.length(secondDim);
	  		int Nz = X.length(thirdDim);
	  		
				for(int k =0; k < X.length(thirdDim);k++){
	  			for(int j =0; j < X.length(secondDim);j++){
	  			for(int i =0; i < X.length(firstDim);i++){
					complex<float>val = complex<float>(1.0/6.0,0.0)*( complex<float>(6.0,0.0)*X(i,j,k) - X((i+1)%Nx,j,k) - X((i+Nx-1)%Nx,j,k) - X(i,(j+1)%Ny,k) - X(i,(j+Ny-1)%Ny,k)  - X(i,j,(k+1)%Nz) - X(i,j,(k+Nz-1)%Nz));
					transform_energy += norm( val );
	  			}}}
	  
		}break;
	
		case(PHASE):{
	 		transform_energy += sum( pow( imag(X),2));
	  	}break;
		
		case(LOWRES):{
			
		
			
		}break;
		
	}
	cout << "Transform Energy = " << transform_energy << endl;
	cout << "EhE Energy = " << EhE_scale << endl;
	
	
	reg_scale = lambda*sqrt( EhE_scale/transform_energy);
	cout << "Reg scale = " << reg_scale << endl;
}



//-----------------------------------------------------
// Regularize Function
//-----------------------------------------------------
void L2REG::regularize( Array< complex<float>,3 > &Rref, Array< complex<float>,3 > &Xref){
	if(lambda==0.0){
		return;
	}
	
	switch(l2_type){
	
	case(NONE):{
	  #pragma omp parallel for
	  for(int k =0; k < Xref.length(thirdDim);k++){
	  for(int j =0; j < Xref.length(secondDim);j++){
	  for(int i =0; i < Xref.length(firstDim);i++){
		  	Rref(i,j,k) += reg_scale*Xref(i,j,k);
	  }}}
	
	}break;
	
	
	case(TV):{
	
	  int Nx = Xref.length(firstDim);
	  int Ny = Xref.length(secondDim);
	  int Nz = Xref.length(thirdDim);
	  #pragma omp parallel for
	  for(int k =0; k < Xref.length(thirdDim);k++){
	  for(int j =0; j < Xref.length(secondDim);j++){
	  for(int i =0; i < Xref.length(firstDim);i++){
		  	Rref(i,j,k) += complex<float>(reg_scale/6.0,0.0)*( complex<float>(6.0,0.0)*Xref(i,j,k) - Xref((i+1)%Nx,j,k) - Xref((i+Nx-1)%Nx,j,k) - Xref(i,(j+1)%Ny,k) - Xref(i,(j+Ny-1)%Ny,k)  - Xref(i,j,(k+1)%Nz) - Xref(i,j,(k+Nz-1)%Nz));
	  }}}
	}break;
	
	case(PHASE):{
	  #pragma omp parallel for
	  for(int k =0; k < Xref.length(thirdDim);k++){
	  for(int j =0; j < Xref.length(secondDim);j++){
	  for(int i =0; i < Xref.length(firstDim);i++){
		  	Rref(i,j,k) += reg_scale*complex<float>( 0.0, imag(Xref(i,j,k)));
	  }}}
	}break;
	
	
	
	
	case(LOWRES):{
			
			
			int rcxres = Xref.length(firstDim);
			int rcyres = Xref.length(secondDim);
			int rczres = Xref.length(thirdDim);
				  
			// Zeropadded bluring
			if( ZeroPad.numElements() !=  ( (rcxres+64)*(rcyres+64)*(rczres+64) ) ){
				ZeroPad.setStorage( ColumnMajorArray<3>() );
				ZeroPad.resize( rcxres+64, rcyres+64, rczres+64 );
			}
			
			// Enforce Smoothness padded
			ZeroPad = complex<float>(0.0,0.0);  // Zero 
			for(int k=0; k < rczres; k++){
			for(int j=0; j < rcyres; j++){
			for(int i=0; i < rcxres; i++){
				ZeroPad(i+32,j+32,k+32)=Xref(i,j,k);
			}}}
						
			fft(ZeroPad);
			float cx = (float)ZeroPad.length(firstDim)/2.0;
			float cy = (float)ZeroPad.length(secondDim)/2.0;
			float cz = (float)ZeroPad.length(thirdDim)/2.0;
			
			#pragma omp parallel for
			for(int k=0; k < ZeroPad.length(thirdDim); k++){
			for(int j=0; j < ZeroPad.length(secondDim); j++){
			for(int i=0; i < ZeroPad.length(firstDim); i++){
				float rad = sqrt(   ((float)k-cz)*((float)k-cz)  + ((float)j-cy)*((float)j-cy) + ((float)i-cx)*((float)i-cx)  );
				rad /= 64.0;
				if(rad < 1){
					ZeroPad(i,j,k) *=  ( 0.355768 - 0.487396*cos( 3.14159265359*rad) + 0.144232*cos(2*3.14159265359*rad) - 0.012604*cos( 3*3.14159265359*rad) )  ;
				}
			}}}	
					
			ifft(ZeroPad);
			
			// Copy back
			for(int k=0; k < rczres; k++){
			for(int j=0; j < rcyres; j++){
			for(int i=0; i < rcxres; i++){
				Rref(i,j,k)+=reg_scale*ZeroPad(i+32,j+32,k+32);
			}}}
			
	}break;
	
	
	
	
	
	
	
	
	
	
	}
	
}




