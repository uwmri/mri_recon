/************************************************
Soft Thresholding Library

Description: This code contains functions to perform thresholding operations


*************************************************/
#include "softthreshold.h"

/*----------------------------------------------
     Constructor  - read command line
 *----------------------------------------------*/ 

SOFTTHRESHOLD::SOFTTHRESHOLD( int numarg, char **pstring){
	thresh = 0.01;	// Default 

#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);

  for(int pos=0; pos < numarg; pos++){
  
  	if (strcmp("-h", pstring[pos] ) == 0) {
	  	printf("\n*********************************************\n");
	  	printf("Soft Threshold Control:\n");
	  	printf("*********************************************\n");
	  	
		printf("-thresh #           :persent max for thresh\n");
	 
		float_flag("-thresh",thresh);
		// trig_flag(KAISER,"-kaiser",kernel_type);
	}
   }    

}


/*----------------------------------------------
     Threshold
 *----------------------------------------------*/ 
void SOFTTHRESHOLD::hard_threshold(  array5D< complex<float> >Coef){
	
	complex <float> max_wave= Coef.max();
	cout << "Max Value " << max_wave << endl;
		
	float std_wave = abs(max_wave)*(thresh);
	for(int e=0; e< Coef.Ne;e++){
	for(int t=0; t< Coef.Nt;t++){
	for(int k=0; k< Coef.Nz; k++){
		for(int j=0; j< Coef.Ny; j++){
			for(int i=0; i< Coef.Nx; i++){	
					Coef[e][t][k][j][i] = Coef[e][t][k][j][i] *complex<float>(max(0.0f,abs(Coef[e][t][k][j][i])-std_wave)/abs(Coef[e][t][k][j][i])); 

	}}}}} 
}
