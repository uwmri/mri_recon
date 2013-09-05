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
	lambda = 0.01;
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
	
	// Input should be energy of EhE X
	// Scales regularization by EhEx/TVx
	
	float transform_energy=0.0;
	switch(l2_type){
	
	case(NONE):{
	  
	  for(int e =0; e < X.length(secondDim);e++){
	  for(int t =0; t < X.length(firstDim);t++){
	  
	  for(int k =0; k < X(0).length(thirdDim);k++){
	  for(int j =0; j < X(0).length(secondDim);j++){
	  for(int i =0; i < X(0).length(firstDim);i++){
			transform_energy += norm( X(e,t)(i,j,k));
	  }}}
	  
	  }}
	
	}break;
	
	
	case(TV):{
	
	  int Nx = X(0).length(firstDim);
	  int Ny = X(0).length(secondDim);
	  int Nz = X(0).length(thirdDim);
	  
	  for(int e =0; e < X.length(secondDim);e++){
	  for(int t =0; t < X.length(firstDim);t++){
	  
	  Array<complex<float>, 3>Xref = X(t,e);
	  for(int k =0; k < Xref.length(thirdDim);k++){
	  for(int j =0; j < Xref.length(secondDim);j++){
	  for(int i =0; i < Xref.length(firstDim);i++){
			complex<float>val = complex<float>(1.0/6.0,0.0)*( complex<float>(6.0,0.0)*Xref(i,j,k) - Xref((i+1)%Nx,j,k) - Xref((i+Nx-1)%Nx,j,k) - Xref(i,(j+1)%Ny,k) - Xref(i,(j+Ny-1)%Ny,k)  - Xref(i,j,(k+1)%Nz) - Xref(i,j,(k+Nz-1)%Nz));
			
			transform_energy += norm( val );
	  }}}
	  
	  }}
	  
	}break;
	
	case(PHASE):{
	  for(int e =0; e < X.length(secondDim);e++){
	  for(int t =0; t < X.length(firstDim);t++){
	  
	  for(int k =0; k < X(0).length(thirdDim);k++){
	  for(int j =0; j < X(0).length(secondDim);j++){
	  for(int i =0; i < X(0).length(firstDim);i++){
			transform_energy += pow( imag(X(t,e)(i,j,k)),2);
	  }}}
	  
	  }}
	}break;
	
	}
	cout << "Transform Energy = " << transform_energy << endl;
	cout << "EhE Energy = " << EhE_scale << endl;
	
	
	reg_scale = lambda*EhE_scale/transform_energy;
}


//-----------------------------------------------------
// Regularize Function
//-----------------------------------------------------
void L2REG::regularize( Array< complex<float>,3 > &Rref, Array< complex<float>,3 > &Xref){
	
	switch(l2_type){
	
	case(NONE):{
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
	  for(int k =0; k < Xref.length(thirdDim);k++){
	  for(int j =0; j < Xref.length(secondDim);j++){
	  for(int i =0; i < Xref.length(firstDim);i++){
		  	Rref(i,j,k) += complex<float>(reg_scale/6.0,0.0)*( complex<float>(6.0,0.0)*Xref(i,j,k) - Xref((i+1)%Nx,j,k) - Xref((i+Nx-1)%Nx,j,k) - Xref(i,(j+1)%Ny,k) - Xref(i,(j+Ny-1)%Ny,k)  - Xref(i,j,(k+1)%Nz) - Xref(i,j,(k+Nz-1)%Nz));
	  }}}
	}break;
	
	case(PHASE):{
	  for(int k =0; k < Xref.length(thirdDim);k++){
	  for(int j =0; j < Xref.length(secondDim);j++){
	  for(int i =0; i < Xref.length(firstDim);i++){
		  	Rref(i,j,k) += reg_scale*complex<float>( 0.0, imag(Xref(i,j,k)));
	  }}}
	}break;
	
	}
	
}




