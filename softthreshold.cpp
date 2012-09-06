/************************************************
Soft Thresholding Library

Description: This code contains functions to perform thresholding operations


*************************************************/
#include "softthreshold.h"

/*----------------------------------------------
     Constructor  - read command line
 *----------------------------------------------*/ 

SOFTTHRESHOLD::SOFTTHRESHOLD( int numarg, char **pstring){
	thresh = 0.0;	// Default 
	threshold =0.0;
	
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
void SOFTTHRESHOLD::get_threshold(  array5D< complex<float> >Coef){
	
	// Get range
	float max_wave= abs( Coef.max() );
	float min_wave= abs( Coef.min() );
	cout << "Image Range:  " << min_wave << " to " << max_wave << endl;
	
	// Estimate histogram (sort)
	int points_found = Coef.Numel;
	int target = (int)( thresh*(double)Coef.Numel);
	int accuracy = (int)(0.001*(double)Coef.Numel);
	cout << "Thresh = " << thresh << " target points " << target << endl;
	
	if(target < 2){
		threshold=0.0;
		return;
	}
		
	// New alogorithm based on compartments
	threshold = 0.5*(max_wave - min_wave);
	float max_t = max_wave;
	float min_t = min_wave;
	int iter = 0;
	while( abs(points_found - target) > accuracy){
		
		// Update to midpoint of good compartment	
		if(iter > 0){
			if( points_found > target){
				// Too many points found - make smaller
				max_t = threshold;
			}else{
				// Too few points found - make larger
				min_t = threshold;
			}
			threshold = 0.5*(max_t + min_t);
			// cout << "New threshold = " << threshold << " Min= " << min_t << " Max= " << max_t << endl;
		}
		iter++;
		
		// Calculate Threshold
		points_found= 0;
		for(int e=0; e< Coef.Ne;e++){
		for(int t=0; t< Coef.Nt;t++){
		#pragma omp parallel for reduction(+:points_found)
		for(int k=0; k< Coef.Nz; k++){
		for(int j=0; j< Coef.Ny; j++){
		for(int i=0; i< Coef.Nx; i++){	
			float v=abs( Coef[e][t][k][j][i]);
			if( v < threshold){
				points_found++;
			}
		}
		}}}}
		// cout << "Threshold = " << threshold << " points found= " << points_found << " error = " <<  (fabs((float)points_found - (float)target)/(float)Coef.Numel) <<endl; 
	
	}
}
	

void SOFTTHRESHOLD::soft_threshold(  array5D< complex<float> >Coef){
	
	cout << "Determined thresh = " << threshold; 	
	// Apply theshold	
	int count = 0;
	for(int e=0; e< Coef.Ne;e++){
	for(int t=0; t< Coef.Nt;t++){
	#pragma omp parallel for reduction(+:count)
	for(int k=0; k< Coef.Nz; k++){
		for(int j=0; j< Coef.Ny; j++){
			for(int i=0; i< Coef.Nx; i++){	
					if( abs( Coef[e][t][k][j][i]) < threshold){
						count++;
						Coef[e][t][k][j][i] = complex<float>(0.0,0.0);
					}else{
						float theta = arg( Coef[e][t][k][j][i] );
						Coef[e][t][k][j][i] = (abs(Coef[e][t][k][j][i])-threshold)*complex<float>(cos(theta),sin(theta));
					}
	}}}}} 
	cout << "Removed " << (float)((float)count/(float)Coef.Numel) << " of the values" << endl;
}
