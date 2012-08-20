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
	
	// Get range
	float max_wave= abs( Coef.max() );
	float min_wave= abs( Coef.min() );
	cout << "Image Range:  " << min_wave << " to " << max_wave << endl;
	
	// Estimate histogram (sort)
	int points_found = Coef.Numel;
	int target = (int)( (1.0-thresh)*(float)Coef.Numel);
	cout << "Thresh = " << thresh << " target points " << target << endl;
	
	// Range of threshold vals
	int N = 64;
	float std_wave = max_wave;
	float *thresh_vals=new float[N];
	int *thresh_count=new int[N];
	thresh_vals[0]=max_wave;
	thresh_count[0]=0;
	for(int p=1; p < N; p++){
		thresh_count[p]=0;
		thresh_vals[p]=thresh_vals[p-1]/1.2;
	}
	
	for(int e=0; e< Coef.Ne;e++){
	for(int t=0; t< Coef.Nt;t++){
	
	#pragma omp parallel for
	for(int k=0; k< Coef.Nz; k++){
	for(int j=0; j< Coef.Ny; j++){
	for(int i=0; i< Coef.Nx; i++){	
		float v=abs( Coef[e][t][k][j][i]);
		for(int p=0; p < N; p++){
			if( v < thresh_vals[p]){
				#pragma omp atomic
				thresh_count[p]++;
			}
		}
	}}}}}
	
	//for(int p=0; p < N; p++){
	//	cout << "Val = " << thresh_vals[p] << " count = " << thresh_count[p] << endl;
	//}
		
	// Find best val
	int p =0;
	while( thresh_count[p] > target){
		p++;
	}
	std_wave = 	thresh_vals[p];
	cout << "Determined thresh = " << std_wave; 	
		
	// Apply theshold	
	int count = 0;
	for(int e=0; e< Coef.Ne;e++){
	for(int t=0; t< Coef.Nt;t++){
	for(int k=0; k< Coef.Nz; k++){
		for(int j=0; j< Coef.Ny; j++){
			for(int i=0; i< Coef.Nx; i++){	
					if( abs( Coef[e][t][k][j][i]) < std_wave){
						count++;
					}
					Coef[e][t][k][j][i] = Coef[e][t][k][j][i] *complex<float>(max(0.0f,abs(Coef[e][t][k][j][i])-std_wave)/abs(Coef[e][t][k][j][i])); 
	}}}}} 
	cout << "Removed " << (float)((float)count/(float)Coef.Numel) << " of the values" << endl;
}
