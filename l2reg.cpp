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
	
	}
	
}




