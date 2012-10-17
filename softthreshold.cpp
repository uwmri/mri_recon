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
void SOFTTHRESHOLD::get_threshold(  const Array< complex<float>,5>&Coef){
	
	// Get range
	float max_wave= max(abs(Coef));
	float min_wave= min(abs(Coef));
	cout << "Image Range:  " << min_wave << " to " << max_wave << endl;
	
	// Estimate histogram (sort)
	int points_found = Coef.numElements();
	int target = (int)( thresh*(double)Coef.numElements());
	int accuracy = (int)(0.001*(double)Coef.numElements());
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
		for(int e=0; e< Coef.extent(fifthDim);e++){
		for(int t=0; t< Coef.extent(fourthDim);t++){
		#pragma omp parallel for reduction(+:points_found)
		for(int k=0; k< Coef.extent(thirdDim); k++){
		for(int j=0; j< Coef.extent(secondDim); j++){
		for(int i=0; i< Coef.extent(firstDim); i++){	
			float v=abs( Coef(i,j,k,t,e));
			if( v < threshold){
				points_found++;
			}
		}
		}}}}
		// cout << "Threshold = " << threshold << " points found= " << points_found << " error = " <<  (fabs((float)points_found - (float)target)/(float)Coef.Numel) <<endl; 
	
	}
}
	

void SOFTTHRESHOLD::soft_threshold(  Array< complex<float>,5 >&Coef){
	
	cout << "Determined thresh = " << threshold; 	
	// Apply theshold	
	int count = 0;
		for(int e=0; e< Coef.extent(fifthDim);e++){
		for(int t=0; t< Coef.extent(fourthDim);t++){
		#pragma omp parallel for reduction(+:count)
		for(int k=0; k< Coef.extent(thirdDim); k++){
		for(int j=0; j< Coef.extent(secondDim); j++){
		for(int i=0; i< Coef.extent(firstDim); i++){	
					if( abs( Coef(i,j,k,t,e)) < threshold){
						count++;
						Coef(i,j,k,t,e) = complex<float>(0.0,0.0);
					}else{
						float theta = arg( Coef(i,j,k,t,e) );
						Coef(i,j,k,t,e) = (abs(Coef(i,j,k,t,e))-threshold)*complex<float>(cos(theta),sin(theta));
					}
	}}}}} 
	cout << "Removed " << (float)((float)count/(float)Coef.numElements()) << " of the values" << endl;
}

void SOFTTHRESHOLD::fista_update(  Array< complex<float>,5 >&X,Array< complex<float>,5 >&X_old,int iteration){
	  
	float A = 1.0 + (iteration - 1.0)/(iteration+1.0);
	float B =     - (iteration - 1.0)/(iteration+1.0);

	// Get Update
	for(int e=0; e< X.extent(fifthDim);e++){
		for(int t=0; t< X.extent(fourthDim);t++){
		for(int k=0; k< X.extent(thirdDim); k++){
		for(int j=0; j< X.extent(secondDim); j++){
		for(int i=0; i< X.extent(firstDim); i++){	
				complex<float>Xn0 = X_old(i,j,k,t,e);								
				complex<float>Xn1 = X(i,j,k,t,e);
				X(i,j,k,t,e) = A*Xn1 + B*Xn0;
				X_old(i,j,k,t,e) = Xn1;
		  
		
	}}}}}			  

}
	

