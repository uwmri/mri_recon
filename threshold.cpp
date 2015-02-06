/************************************************
  Thresholding Library

Description: This code contains functions to perform thresholding operations

Authors:
Kevin Johnson (kmjohnso3@wisc.edu)
Grzegorz Bauman (gbauman@wisc.edu)

Includes: 
Soft/hard thresholding
Fractional thresholding (based on percentage)
Visu thresholding
Bayes thresholding
SURE thresholding

 *************************************************/
#include "threshold.h"
#include "io_templates.hpp"
using namespace NDarray;

THRESHOLD::THRESHOLD() {

}

/*----------------------------------------------
  Constructor  - read command line
 *----------------------------------------------*/ 

THRESHOLD::THRESHOLD( int numarg, char **pstring) {

	soft = true;
	thapp = false;
	thresh = 0.0;   // Default 
	group_complex = true;
	global_threshold =0.0;
	waveL = 4; // Wavelet Levels need better code to handle
	VERBOSE = false;
	temporal = false;
	threshold_type = TH_NONE;
	noise_scale = 1.0;
	noise_est_type = NOISE_FRAME;
	
#define trig_flag(num,name,val) }else if(strcmp(name,pstring[pos]) == 0){ val=num;
#define bool_flag(name,val) }else if(strcmp(name,pstring[pos]) == 0){ val=true;
#define float_flag(name,val) }else if(strcmp(name,pstring[pos]) == 0){ pos++; val=atof(pstring[pos]);
#define int_flag(name,val) }else if(strcmp(name,pstring[pos]) == 0){ pos++; val=atoi(pstring[pos]);
#define char_flag(name,val)   }else if(strcmp(name,pstring[pos]) == 0){ pos++; strcpy(val,pstring[pos]); 
#define char_trig_flag(num,name,name2,val)   }else if(strcmp(name,pstring[pos]) == 0){ pos++;  if(strcmp(name2,pstring[pos])==0){ val=num; pos++;}


	for(int pos=0; pos < numarg; pos++){

		if (strcmp("-h", pstring[pos] ) == 0) {
			float_flag("-thresh",thresh);
			char_flag("-th",th_type);
			char_flag("-thmode",th_mode);
			bool_flag("-thlow",thapp);
			trig_flag(true,"-verbose",VERBOSE);
			float_flag("-noise_scale",noise_scale);
			
		}else if(strcmp("-thmode",pstring[pos]) == 0) {
			pos++;
			if( pos==numarg){
				cout << "Please provide threshold type..hard/soft" << endl;
				trig_flag(false,"soft",soft);
				trig_flag(true,"hard",soft);
			}
		}else if(strcmp("-threshold_type",pstring[pos]) == 0) {
			pos++;
			if( pos==numarg){
				cout << "Please provide threshold type..fraction/bayes/visu/sure/fixed" << endl;
				trig_flag(TH_FIXED,"fixed",threshold_type);
				trig_flag(TH_FRACTION,"fraction",threshold_type);
				trig_flag(TH_VISU,"visu",threshold_type);
				trig_flag(TH_BAYES,"bayes",threshold_type);
				trig_flag(TH_SURE,"sure",threshold_type);
			}
		}else if(strcmp("-th_group_complex []",pstring[pos]) == 0) {
			pos++;
			if (pos == numarg) {
				cout << "Please specify setting for grouped complex thresholding..on/off (default is on)" << endl;
				trig_flag(false,"off",group_complex);
				trig_flag(true,"on",group_complex);
			}
		}else if(strcmp("-noise_est_type",pstring[pos]) == 0) {
			pos++;
			if( pos==numarg){
				cout << "Please provide noise estimation type..global/frame/first/last" << endl;
				trig_flag(NOISE_GLOBAL,"global",noise_est_type);
				trig_flag(NOISE_FRAME,"frame",noise_est_type);
				trig_flag(NOISE_FIRST,"first",noise_est_type);
				trig_flag(NOISE_LAST,"last",noise_est_type);
			}
		}
	}    
}

void THRESHOLD::help_message() {
	cout << "----------------------------------------------" << endl;
	cout << "   Thresholding Control " << endl;
	cout << "----------------------------------------------" << endl;
	help_flag("-threshold_type []","thresholding method l1reg/fraction/visu/bayes/sure");
	help_flag("-thresh []","fraction of coefficient to be removed for fractional method (default=0.0)");
	help_flag("-noise_est_type []","method for noise estimation global(estimate using all frames)/frame(estimate for each frame)/last(estimate using last frame)/first(estimate using first frame) (default=frame)");
	help_flag("-noise_scale []","scale factor for noise estimate(s) (e.g. for zero-padding correction)");
	help_flag("-thmode []","soft/hard (default=soft)");
	help_flag("-th_group_complex []", "on/off --> on performs group soft thresholding of complex data. When off (default) re/im data is treated individually");
	help_flag("-thlow","thresholding of the lowest resolusion band");
}

// Executes a chosen thresholding method
void THRESHOLD::exec_threshold( Array<  Array< complex<float>,3>,  2>&Coef, WAVELET3D &wave){

	switch(threshold_type){
		case TH_FIXED:
		case TH_FRACTION:
		case TH_VISU:{
				     thresholding(Coef,  global_threshold);
		}break;
				
		case TH_BAYES:
		case TH_SURE: {
				     exec_multibandthreshold(Coef, wave);
		}break;
		
		case TH_NONE:
		default:{
				return;
			}
	}

}


/*--------------------------------------------------------------------
  Upodate threshold. Seperated for cycle spin
 *------------------------------------------------------------------*/ 
 void THRESHOLD::update_threshold( Array<  Array< complex<float>,3>,  2>&Coef, WAVELET3D &wave, float scale){

	switch(threshold_type){

	
		case TH_FIXED:{
			global_threshold = scale; // This is odd?
		}break;

		case TH_FRACTION:{
			global_threshold = get_threshold( Coef, thresh); 
			global_threshold *= scale; 
		}break;
		
		case TH_VISU: {
			get_visuthreshold( Coef, wave);  
			global_threshold *= scale;
		}break;
		
		case TH_BAYES:{
			get_bayesthreshold( Coef, wave);
			subband_threshold *= scale; // Cycle spin correction  
		}break;
		
		case TH_SURE: {
			get_surethreshold( Coef, wave);
			subband_threshold *= scale; // Cycle spin correction  
		}break;
		
		case TH_NONE:
		default:{
				return;
			}
	}

}

/*--------------------------------------------------------------------
  Calc threshold based on percentage removed in one or all bands
  used also to get median of Coef for thresh=0.5
 *------------------------------------------------------------------*/ 
float THRESHOLD::get_threshold(Array<Array< complex<float>,3>,2>&Coef, float fraction){
	
	float value;

	if(VERBOSE){
		 cout << "Calc Range" << endl << flush;
	}
	
	int Nz = Coef(0,0).length(thirdDim);
	
	// Get range
	float *max_wave_t = new float[Nz];
	float *min_wave_t = new float[Nz];
	for(int k=0; k< Nz; k++){
		max_wave_t[k] = 0.0;
		min_wave_t[k] = 0.0;
	}
		
	if(VERBOSE) {
		cout << " ------- Computing min and max coefficients " << endl;
	}
	for(int e=0; e< Coef.length(secondDim);e++){
		for(int t=0; t< Coef.length(firstDim);t++){
			Array< complex<float>,3> XX;
			XX.reference( Coef(t,e) ); 

			#pragma omp parallel for 
			for(int k=0; k< XX.extent(thirdDim); k++){
			for(int j=0; j< XX.extent(secondDim); j++){
			for(int i=0; i< XX.extent(firstDim); i++){    
			
				if (group_complex) {
					float v=abs( XX(i,j,k));
					max_wave_t[k] = ( max_wave_t[k] > v ) ? ( max_wave_t[k] ) : ( v );
					min_wave_t[k] = ( min_wave_t[k] < v ) ? ( min_wave_t[k] ) : ( v );
				}else {
					float v=abs( real(XX(i,j,k)));
					max_wave_t[k] = ( max_wave_t[k] > v ) ? ( max_wave_t[k] ) : ( v );
					min_wave_t[k] = ( min_wave_t[k] < v ) ? ( min_wave_t[k] ) : ( v );
					v=abs( imag(XX(i,j,k)));
					max_wave_t[k] = ( max_wave_t[k] > v ) ? ( max_wave_t[k] ) : ( v );
					min_wave_t[k] = ( min_wave_t[k] < v ) ? ( min_wave_t[k] ) : ( v );
				}
			}}}	// Spatial			
	}}
	
	float max_wave= max_wave_t[0];
	float min_wave= min_wave_t[0];
	
	for(int k=0; k< Nz; k++){
		max_wave = ( max_wave_t[k] > max_wave ) ? ( max_wave_t[k] ) : ( max_wave );
		min_wave = ( min_wave_t[k] < max_wave ) ? ( min_wave_t[k] ) : ( min_wave );	
	}
	
	
	
	if(VERBOSE) cout << "Image Range:  " << min_wave << " to " << max_wave << endl;

	// Estimate histogram (sort)
	int total_points = Coef.numElements()*Coef(0,0).numElements() * (group_complex ? 1 : 2); // twice as many points without complex grouping
	int points_found = total_points;
	int target = (int)( fraction*(double)total_points);
	int accuracy = (int)(0.001*(double)total_points);
	if(VERBOSE) cout << "Thresh = " << fraction << " target points " << target << endl;

	if(target < 2){
		return(0.0);
	}

	// New alogorithm based on compartments
	value  = 0.5*(max_wave - min_wave);
	float max_t = thresh*max_wave;
	float min_t = min_wave;
	int iter = 0;
	while( abs(points_found - target) > accuracy){

		// Update to midpoint of good compartment       
		if(iter > 0){
			if( points_found > target){
				// Too many points found - make smaller
				max_t = value;
			}else{
				// Too few points found - make larger
				min_t = value;
			}
			value = 0.5*(max_t + min_t);
			if(VERBOSE) cout << "Points = " << points_found << " New threshold = " << value << " Min= " << min_t << " Max= " << max_t << endl;
		}
		iter++;

		// Calculate Threshold
		points_found= 0;
		for(int e=0; e< Coef.extent(secondDim);e++){
		for(int t=0; t< Coef.extent(firstDim);t++){
			Array< complex<float>,3> XX = Coef(t,e); 

			#pragma omp parallel for reduction(+:points_found)			
			for(int k=0; k< XX.extent(thirdDim); k++){
			for(int j=0; j< XX.extent(secondDim); j++){
			for(int i=0; i< XX.extent(firstDim); i++){
				if (group_complex) {
					float v=abs( XX(i,j,k));
					if( v < value){
						points_found++;
					}
				} else {
					float v=abs(real(XX(i,j,k)));
					if( v < value){
						points_found++;
					}
					v = abs(imag(XX(i,j,k)));
					if( v < value){
						points_found++;
					}
				}
					
			}}}//Spatial
		}}
			
	}//While
	
	return(value);
}

/*---------------------------------------------------------------------------------------
  Calc threshold based on universal threshold method using noise variance (VisuShrink)
 *-------------------------------------------------------------------------------------*/ 

void THRESHOLD::get_visuthreshold(Array<Array< complex<float>,3>,2>&Coef, WAVELET3D &wave){

	// Populated noise variable
	
	switch(noise_est_type){
		case(NOISE_GLOBAL):{
					   robust_noise_estimate( Coef, wave);
					   if(VERBOSE) cout<<"Estimated noise: "<<noise<<endl;
				   }break;
		case(NOISE_FIRST):{
					   Array< Array<complex<float>,3>,2> Cf;
					   Cf.setStorage(ColumnMajorArray<2>());
					   Cf.resize(1,1);
					   Cf(0,0).reference( Coef(0,0)); 
					   
					   robust_noise_estimate( Cf, wave);
					   if(VERBOSE) cout<<"Estimated noise: "<<noise<<endl;
				   }break;
		case(NOISE_LAST):{
					   Array< Array<complex<float>,3>,2> Cf;
					   Cf.setStorage(ColumnMajorArray<2>());
					   Cf.resize(1,1);
					   Cf(0,0).reference(Coef(Coef.length(firstDim)-1,Coef.length(secondDim-1))); 
					   
					   robust_noise_estimate( Cf, wave);
					   if(VERBOSE) cout<<"Estimated noise: "<<noise<<endl;
				   }break;
		case(NOISE_FRAME):{
					  // This is an invalid choice for visu thresholding.  Default to global.
					  noise_est_type = NOISE_GLOBAL; 
					  robust_noise_estimate( Coef, wave);
					  if(VERBOSE) cout<<"Estimated noise: "<<noise<<endl;
				  }break;
		default:{
				cout << "Bad value for noise estimation type" << endl;
				exit(1);
			}
	}

	robust_noise_estimate(Coef,wave);

	// VisuShrink tends to overestimate the threshold especially for large number of samples
	// Reduce it with an empirical factor
	// The original formula is: sigma_noise * sqrt(2*N_samples)
	// May need a better scaling factor
	global_threshold = 1.0/8.0*noise*sqrt(2*log(Coef(0,0).numElements()*Coef.numElements() / 8.0 ) );
}


/*---------------------------------------------------------------------------------------
  Calc threshold based on Bayesian threshold method
  	-Modified to account for complex data. 
 *-------------------------------------------------------------------------------------*/ 

void THRESHOLD::get_bayesthreshold(Array<Array< complex<float>,3>,2>&Coef, WAVELET3D &wave){
	
	subband_threshold.setStorage(ColumnMajorArray<5>());
	subband_threshold.resize( wave.L[0]+1, wave.L[1]+1, wave.L[2]+1, Coef.length(firstDim), Coef.length(secondDim)); 
	subband_threshold = 0.0;		
	
	int Ne = Coef.length(secondDim);
	int Nt = Coef.length(firstDim);

	// Get Noise
	Array< Array<complex<float>,3>,2> Cf = Alloc5DContainer< complex<float> >(Coef(0,0).length(firstDim),Coef(0,0).length(secondDim),Coef(0,0).length(thirdDim),1,1);

	switch(noise_est_type){
		case(NOISE_GLOBAL):{
					   robust_noise_estimate( Coef, wave);
					   if(VERBOSE) cout<<"Estimated noise: "<<noise<<endl;
				   }break;
		case(NOISE_FIRST):{
					   Cf(0,0).resize(Coef(0,0).shape());
					   Cf(0,0) = Coef(0,0); 
					   robust_noise_estimate( Cf, wave);
					   if(VERBOSE) cout<<"Estimated noise: "<<noise<<endl;
				   }break;
		case(NOISE_LAST):{
					 Cf(0,0).resize(Coef(Nt-1,Ne-1).shape());
					 Cf(0,0) = Coef(Nt-1,Ne-1); 
					 robust_noise_estimate( Cf, wave);
					 if(VERBOSE) cout<<"Estimated noise: "<<noise<<endl;
				   }break;
		case(NOISE_FRAME):{
					  // compute noise for each frame/encode within the loop below
				  }break;
		default:{
				cout << "Bad value for noise estimation type" << endl;
				exit(1);
			}
	}
	
	// Threshold each level
		

	// Get Update
	for(int e=0; e< Ne;e++){
		for(int t=0; t< Nt;t++){
			
			if (noise_est_type == NOISE_FRAME) {
				Cf(0,0).resize(Coef(t,e).shape());
				Cf(0,0) = Coef(t,e); 
				robust_noise_estimate( Cf, wave);
				if(VERBOSE) cout<<"Estimated noise (" << t << "," << e << "): " << noise << endl;
			}

			#pragma omp parallel for
			for( int lz=0; lz < wave.L[2]+1; lz++){
			for( int ly=0; ly < wave.L[1]+1; ly++){
			for( int lx=0; lx < wave.L[0]+1; lx++){
			Range rx,ry,rz;
			wave.get_subband_range(rx,ry,rz,lx,ly,lz);	// Get Range

			// Get the reference volume
			Array<complex<float>,3> Subband;
			Subband.reference( Coef(t,e)(rx,ry,rz));
		
			// Get and store subband threshold	(	
			subband_threshold(lx,ly,lz,t,e) = get_bayesthreshold_subband(Subband);
			// cout << "Level(" <<lx<<","<<ly<<","<<lz<<") - thresh = " << subband_threshold(lx,ly,lz,0,0) << endl;
	}}}
	
	}}

}

void THRESHOLD::get_surethreshold(Array<Array< complex<float>,3>,2>&Coef, WAVELET3D &wave){
	
	subband_threshold.setStorage(ColumnMajorArray<5>());
	subband_threshold.resize( wave.L[0]+1, wave.L[1]+1, wave.L[2]+1, Coef.length(firstDim), Coef.length(secondDim)); 
	subband_threshold = 0.0;		
	
	// Get Noise
	robust_noise_estimate( Coef, wave);
	if(VERBOSE) cout<<"Estimated noise: "<<noise<<endl;
	
	// Threshold each level
	int Ne = Coef.length(secondDim);
	int Nt = Coef.length(firstDim);
		
	// Get Update
	for(int e=0; e< Ne;e++){
		for(int t=0; t< Nt;t++){

			#pragma omp parallel for
			for( int lz=0; lz < wave.L[2]+1; lz++){
			for( int ly=0; ly < wave.L[1]+1; ly++){
			for( int lx=0; lx < wave.L[0]+1; lx++){
				Range rx,ry,rz;
				wave.get_subband_range(rx,ry,rz,lx,ly,lz);	// Get Range

				// Get the reference volume
				Array<complex<float>,3> Subband;
				Subband.reference( Coef(t,e)(rx,ry,rz));
		
				// Get and store subband threshold	(	
				subband_threshold(lx,ly,lz,t,e) = get_surethreshold_subband(Subband);
				// cout << "Level(" <<lx<<","<<ly<<","<<lz<<") - thresh = " << subband_threshold(lx,ly,lz,0,0) << endl;
			}}}
		}
	}
}

void THRESHOLD::exec_multibandthreshold(Array<Array< complex<float>,3>,2>&Coef, WAVELET3D &wave){
	
	int Ne = Coef.length(secondDim);
	int Nt = Coef.length(firstDim);
		
	// Get Update
	for(int e=0; e< Ne;e++){
		for(int t=0; t< Nt;t++){
		
		#pragma omp parallel for
		for( int lz=0; lz < wave.L[2]+1; lz++){
		for( int ly=0; ly < wave.L[1]+1; ly++){
		for( int lx=0; lx < wave.L[0]+1; lx++){
			Range rx,ry,rz;
			wave.get_subband_range(rx,ry,rz,lx,ly,lz);	// Get Range

			// Get the reference volume
			Array<complex<float>,3> Subband;
			Subband.reference( Coef(t,e)(rx,ry,rz));
		
			// Get and store subband threshold	(	
			thresholding( Subband , subband_threshold(lx,ly,lz,t,e) );
		
		}}}
	
	}}

}


float THRESHOLD::get_bayesthreshold_subband( Array< complex<float>,3> &XX){
	
	float threshold;
	
	// Get Std 
	double N = (double)XX.numElements();
	double sumXX = 0.0; 
	double maxCoef = 0.0;
	complex<double> meanXX = 0.0;
    	double meanR,meanI = 0.0;
	meanR = (double)(sum(real(XX)))/N;
	meanI = (double)(sum(imag(XX)))/N;
	meanXX = complex<double>(meanR,meanI);

	// Standard Deviation of Image
			
	for(int k=0; k< XX.extent(thirdDim); k++){
		for(int j=0; j< XX.extent(secondDim); j++){
			for(int i=0; i< XX.extent(firstDim); i++){    
					complex<double> val = complex<double>(XX(i,j,k));
					maxCoef = max(abs(real(val)),maxCoef);
					maxCoef = max(abs(imag(val)),maxCoef);
					sumXX += norm( val - meanXX ); 
	}}}
		
	float wave_dev = (float)( sumXX/(2.0*N) );
	
	// if(VERBOSE) cout << "N  = " << N << endl;
	// if(VERBOSE) cout << "Wave Deviation = " << wave_dev << endl;
		
	// Now Estimate Threshold	
	float sigmaS = sqrt( max( wave_dev - powf(noise,2),(float)0.0)); 
	
	if(sigmaS == 0.0) { 
		threshold = maxCoef;
	} else {
		threshold = pow(noise,2)/sigmaS;
	}
	
	return(threshold);

}


/*---------------------------------------------------------------------------------------
  Calc threshold based on SURE threshold method (still doen't work)
 *-------------------------------------------------------------------------------------*/ 


float THRESHOLD::sure_cost( Array< complex<float>,3>&Coef, float thresh){

	float count = 0.0;
	float energy = 0.0;
	for( Array<complex<float>,3>::iterator miter=Coef.begin(); miter!=Coef.end(); miter++){
		count += (float)( abs(*miter) < thresh );
		energy += pow( min( abs(*miter), thresh), 2.0f);	
	}
	
	float N = Coef.numElements();
	// cout << "Noise = " << noise << ", N = " << N << ", Energy " << energy << " ,count=" << count << endl;
	return( N*noise*noise - 2.0*noise*noise*count + energy );

}


float THRESHOLD::get_surethreshold_subband( Array< complex<float>,3>&Coef){

	// Search for best 
	float min_risk = 0.0;
	float min_thresh = 0.0;
	
	int sureres=40;
	
	for(int pos=0; pos < sureres; pos++){
		float thresh =  1.0/8.0*noise*sqrt(2*log((float)Coef.numElements()))*((float)pos)/(((float)sureres)-1.0);		
		
		float risk = sure_cost(Coef, thresh);
		if( (risk < min_risk) || (pos==0) ){
			min_risk =risk;
			min_thresh = thresh;
		}else{
			// Should be convex
			break;		
		}
				
	}
	
	return(min_thresh);
		
}


void THRESHOLD::robust_noise_estimate(Array<Array< complex<float>,3>,2>&Coef, WAVELET3D &wave){
	
	// To get noise from HHH band
	Range rx,ry,rz;
	wave.get_subband_range(rx,ry,rz,0,0,0);
	
	int Ne = Coef.length(secondDim);
	int Nt = Coef.length(firstDim);

	// The (now commented) code block below may cause problems with multi-frame/encode data, especially when the number
	// of frames (amount of sampled data) varies significantly between frames/encodes (i.e. some are noisier than others).
/* *******************************************************************************************	
	// Get 3D volume
	Array<complex<float>,3> HHHref;
	HHHref.reference( Coef(Coef.length(firstDim)-1,Coef.length(secondDim)-1)(rx,ry,rz));
	

	// Convert to 5D for get_threshold
	Array< Array<complex<float>,3>,2> DH;
	DH.setStorage(ColumnMajorArray<2>());
	DH.resize(1,1);
	DH(0,0).reference( HHHref);
**********************************************************************************************/

	// Compute "global" noise from all provided frames and encodes (may be all or some subset)
	Array< complex<float>,3>HHHref = Coef(0,0)(rx,ry,rz);
	Array< Array<complex<float>,3>,2> DH = Alloc5DContainer< complex<float> >(HHHref.length(firstDim),HHHref.length(secondDim),HHHref.length(thirdDim),Nt,Ne);
	
	for(int e=0; e< Ne;e++){
		for(int t=0; t< Nt;t++){
			DH(t,e) = Coef(t,e)(rx,ry,rz);
		}
	}

	cout << "Size of array for noise estimation -- (nt,ne): " << DH.shape() << " --- (nx,ny,nz): " << DH(0,0).shape() << endl;
		
	// Noise Estimate = median/ 0.6745	 - noise scale is needed for zero filled data
	noise = get_threshold(DH,0.5)/0.6745*noise_scale;
}	 



// Option of thresh of appriximation wavelets or all wavelets, hard vs soft threshold
void THRESHOLD::thresholding(Array<Array< complex<float>,3>,2>&Coef, float value){

	// if(VERBOSE) cout << "Determined thresh = " << value <<endl;

	// Apply theshold       
	int count = 0;
	
	int Nx = Coef(0,0).length(firstDim);
	int Ny = Coef(0,0).length(secondDim);
	int Nz = Coef(0,0).length(thirdDim);
	int Ne = Coef.length(secondDim);
	int Nt = Coef.length(firstDim);
		
	// Get Update
	for(int e=0; e< Ne;e++){
		for(int t=0; t< Nt;t++){
			Array< complex<float>,3>XX = Coef(t,e);
			#pragma omp parallel for reduction(+:count)
			for(int k=0; k< Nz; k++){
				for(int j=0; j< Ny; j++){
					for(int i=0; i< Nx; i++){

						if (group_complex) {
							float cc = abs(XX(i,j,k));
							if( cc <= value ){
								XX(i,j,k) = complex<float>(0.0,0.0);
							}else{
								count++;
								XX(i,j,k) *= (1.0f - value/cc);
							}
						} else {
							float newR = max(0.0f,real(XX(i,j,k)) - value) - max(0.0f, -real(XX(i,j,k)) - value);
							float newI = max(0.0f,imag(XX(i,j,k)) - value) - max(0.0f, -imag(XX(i,j,k)) - value);
							if(soft==true){
								XX(i,j,k) = complex<float>(newR,newI);
							}else{
								newR = (newR == 0.0) ? 0.0 : real(XX(i,j,k));
								newI = (newI == 0.0) ? 0.0 : imag(XX(i,j,k));
								XX(i,j,k) = complex<float>(newR,newI);
							}


						}

				}}}
	}} 
	// if(VERBOSE) cout << "Removed " << (float)((float)count/(float)Coef.numElements()) << " of the values" << endl;
}


// Option of thresh of appriximation wavelets or all wavelets, hard vs soft threshold
void THRESHOLD::thresholding( Array< complex<float>,3> &XX, float value){

	//if(VERBOSE) cout << "Determined thresh = " << value <<endl;

	// Apply theshold       
	int count = 0;
	
	int Nx = XX.length(firstDim);
	int Ny = XX.length(secondDim);
	int Nz = XX.length(thirdDim);
		
	// Get Update
	#pragma omp parallel for reduction(+:count)
	for(int k=0; k< Nz; k++){
		for(int j=0; j< Ny; j++){
			for(int i=0; i< Nx; i++){
				
		
				if (group_complex) {
					float cc = abs(XX(i,j,k));
					if( cc <= value ){
						XX(i,j,k) = complex<float>(0.0,0.0);
					}else{
						count++;
						XX(i,j,k) *= (1.0f - value/cc);
					}
				} else {
					float newR = max(0.0f,real(XX(i,j,k)) - value) - max(0.0f, -real(XX(i,j,k)) - value);
					float newI = max(0.0f,imag(XX(i,j,k)) - value) - max(0.0f, -imag(XX(i,j,k)) - value);
					if(soft==true){
						XX(i,j,k) = complex<float>(newR,newI);
					}else{
						newR = (newR == 0.0) ? 0.0 : real(XX(i,j,k));
						newI = (newI == 0.0) ? 0.0 : imag(XX(i,j,k));
						XX(i,j,k) = complex<float>(newR,newI);
					}
				}
				 
			
		
	}}}
	//if(VERBOSE) cout << "Removed " << (float)((float)count/(float)XX.numElements()) << " of the values" << endl;
}




void THRESHOLD::setThresholdMethod(int type) {
	threshold_type = type;
}

int THRESHOLD::getThresholdMethod() {
	return threshold_type;
}

void THRESHOLD::setTemporalThresholding(bool flag) {
	temporal=flag;
}

bool THRESHOLD::getTemporalThresholding() {
	return temporal;
}

void THRESHOLD::fista_update(  Array<Array< complex<float>,3>,2>&X,Array<Array< complex<float>,3>,2>&X_old,int iteration){
	
	float A = 1.0 + (iteration - 1.0)/(iteration+1.0);
	float B =     - (iteration - 1.0)/(iteration+1.0);

	int Nx = X_old(0).length(firstDim);
	int Ny = X_old(0).length(secondDim);
	int Nz = X_old(0).length(thirdDim);
	int Ne = X_old.length(secondDim);
	int Nt = X_old.length(firstDim);
	
	cout << "fista_update: Matrix Size = " << Nx << " x " << Ny << " x " << Nz << " x " << Nt << " x " << Ne << endl;
	
	// Get Update
	for(int e=0; e< Ne;e++){
		for(int t=0; t< Nt;t++){
			Array< complex<float>,3>XX     = X(t,e);
			Array< complex<float>,3>XX_old = X_old(t,e);
			
			for(int k=0; k< Nz; k++){
				for(int j=0; j< Ny; j++){
					for(int i=0; i< Nx; i++){
						
						// cout << "Pos = " << i << "," << j << "," << k << "," << t << "," << e << endl;
	
						
						complex<float>Xn0 = XX_old(i,j,k);
						complex<float>Xn1 = XX(i,j,k);
						XX(i,j,k) = A*Xn1 + B*Xn0;
						XX_old(i,j,k) = Xn1;
			}}} // Spatial
	}}
	cout << "Done with FISTA" << endl << flush;

}

