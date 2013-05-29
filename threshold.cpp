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

/*----------------------------------------------
  Constructor  - read command line
 *----------------------------------------------*/ 

THRESHOLD::THRESHOLD( int numarg, char **pstring) {

	soft = true;
	thapp = false;
	thresh = 0.0;   // Default 
	threshold =0.0;
	waveL = 4; // Wavelet Levels need better code to handle
	VERBOSE = 0;
	temporal = false;
	threshold_type = TH_NONE;
	
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
			trig_flag(1,"-verbose",VERBOSE);
		
		}else if(strcmp("-thmode",pstring[pos]) == 0) {
			pos++;
			if( pos==numarg){
				cout << "Please provide threshold type..fraction/bayes/visu/sure" << endl;
				trig_flag(false,"soft",soft);
				trig_flag(true,"hard",soft);
			}
		}else if(strcmp("-threshold_type",pstring[pos]) == 0) {
			pos++;
			if( pos==numarg){
				cout << "Please provide threshold type..fraction/bayes/visu/sure" << endl;
				trig_flag(TH_FRACTION,"fraction",threshold_type);
				trig_flag(TH_VISU,"visu",threshold_type);
				trig_flag(TH_BAYES,"bayes",threshold_type);
				trig_flag(TH_SURE,"sure",threshold_type);
			}
		}
	}    
}

void THRESHOLD::help_message() {
	cout << "----------------------------------------------" << endl;
	cout << "   Thresholding Control " << endl;
	cout << "----------------------------------------------" << endl;
	help_flag("-thresh_type []","thresholding method fraction/visu/bayes/sure");
	help_flag("-thresh []","fraction of coefficient to be removed for fractional method (default=0.0)");
	help_flag("-thmode []","soft/hard (default=soft)");
	help_flag("-thlow","thresholding of the lowest resolusion band");
}

// Executes a chosen thresholding method
void THRESHOLD::exec_threshold(Array<complex<float>,5>&Coef) {

	switch(threshold_type){

		case TH_FRACTION:{
					 exec_fractionthreshold(Coef);
				 }break;
		case TH_VISU: {
				      exec_visuthreshold(Coef);
			      }break;
		case TH_BAYES:
		case TH_SURE: {
				      exec_multibandthreshold(Coef);
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
void THRESHOLD::get_threshold(Array< complex<float>,5>&Coef){

	// Get range
	float *max_wave_t = new float[Coef.extent(thirdDim)];
	float *min_wave_t = new float[Coef.extent(thirdDim)];
		
	for(int e=0; e< Coef.extent(fifthDim);e++){
		for(int t=0; t< Coef.extent(fourthDim);t++){
#pragma omp parallel for 
			for(int k=0; k< Coef.extent(thirdDim); k++){
			for(int j=0; j< Coef.extent(secondDim); j++){
			for(int i=0; i< Coef.extent(firstDim); i++){    
				float v=abs( Coef(i,j,k,t,e));
				max_wave_t[k] = ( max_wave_t[k] > v ) ? ( max_wave_t[k] ) : ( v );
				min_wave_t[k] = ( min_wave_t[k] < v ) ? ( min_wave_t[k] ) : ( v );
				
				
			}}}				
		}
	}
	
	float max_wave= max_wave_t[0];
	float min_wave= min_wave_t[0];
	
	for(int k=0; k< Coef.extent(thirdDim); k++){
		max_wave = ( max_wave_t[k] > max_wave ) ? ( max_wave_t[k] ) : ( max_wave );
		min_wave = ( min_wave_t[k] < max_wave ) ? ( min_wave_t[k] ) : ( min_wave );	
	}
	
	
	
	if(VERBOSE) cout << "Image Range:  " << min_wave << " to " << max_wave << endl;

	// Estimate histogram (sort)
	int points_found = Coef.numElements();
	int target = (int)( thresh*(double)Coef.numElements());
	int accuracy = (int)(0.001*(double)Coef.numElements());
	if(VERBOSE) cout << "Thresh = " << thresh << " target points " << target << endl;

	if(target < 2){
		threshold=0.0;
		return;
	}

	// New alogorithm based on compartments
	threshold = 0.5*(max_wave - min_wave);
	float max_t = thresh*max_wave;
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
			if(VERBOSE) cout << "Points = " << points_found << " New threshold = " << threshold << " Min= " << min_t << " Max= " << max_t << endl;
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
	}
}

/*---------------------------------------------------------------------------------------
  Calc threshold based on universal threshold method using noise variance (VisuShrink)
 *-------------------------------------------------------------------------------------*/ 

void THRESHOLD::get_visuthreshold(Array< complex<float>,5 >&Coef){

	int sx = Coef.extent(firstDim); 
	int sy = Coef.extent(secondDim);
	int sz = Coef.extent(thirdDim);
	int st = Coef.extent(fourthDim);
	int se = Coef.extent(fifthDim);

	// HHH band
	Array<complex<float>,5> Cf = Coef(Range(sx/2,sx-1),Range(sy/2,sy-1),Range(sz/2,sz-1),Range::all(),Range::all());

	// Get median for noise estimation, sigma_noise = median(Coef)/0.6745
	thresh=0.5;
	get_threshold(Cf);

	noise = threshold/0.6745;

	// VisuShrink tends to overestimate the threshold especially for large number of samples
	// Reduce it with an empirical factor
	// The original formula is: sigma_noise * sqrt(2*N_samples)
	// May need a better scaling factor

	threshold = 1.0/8.0*noise*sqrt(2*log((float)sx*(float)sy*(float)sz*(float)st*(float)se/8.0));
}

/*---------------------------------------------------------------------------------------
  Calc threshold based on Bayesian threshold method
  	-Modified to account for complex data. 
 *-------------------------------------------------------------------------------------*/ 

void THRESHOLD::get_bayesthreshold(Array< complex<float>,5 >&Coef){
	
	// Get Std 
	double sumXX = 0.0; 
	complex<double> sumX = 0.0;
		
	// Standard Deviation of Image
	for(int e=0; e< Coef.extent(fifthDim);e++){
		for(int t=0; t< Coef.extent(fourthDim);t++){
			for(int k=0; k< Coef.extent(thirdDim); k++){
				for(int j=0; j< Coef.extent(secondDim); j++){
					for(int i=0; i< Coef.extent(firstDim); i++){    
						complex<double>val = complex<double>(Coef(i,j,k,t,e));
						sumXX += pow(abs(val),2);
						sumX += val;
					}
	}}}}
	
	double N = (double)Coef.size();
	float wave_dev = (float)sqrt( sumXX/N - abs(sumX)/(N*N));
		
	// Now Estimate Threshold	
	float sigmaS = sqrt( max( wave_dev - powf(noise,2),(float)0.0)); 
	
	if(sigmaS<=0.0) { 
		threshold = max(abs(Coef));
	} else {
		sigmaS = sqrt(sigmaS);
		threshold = pow(noise,2)/sigmaS;
	}

}


/*---------------------------------------------------------------------------------------
  Calc threshold based on SURE threshold method (still doen't work)
 *-------------------------------------------------------------------------------------*/ 

void THRESHOLD::get_surethreshold(Array< complex<float>,5 >&Coef) {

	//	Formula: sigmaˆ2+1/n(sum(min(abs(x),thres)ˆ2))-2*sigmaˆ2/n*sum(abs(x) < thres)

	Array<float,1> thr;
	thr.resize(20);
	Array<float,1> thc;
	thc.resize(20);
	float tmp1, tmp2;

	thr[0]=0.0;
	for(int i=1;i<20;i++) {
		thr(i) = noise*sqrt(2.0*log((float)Coef.numElements()))*(((float)i)/19.0);			// that seems to be the problem, the values are too small
	}

	for(int i=0;i<20;i++) {
		tmp1 = 0.0;
		tmp2 = 0.0;

		for(int e=0; e<Coef.extent(fifthDim);e++){
			for(int t=0; t<Coef.extent(fourthDim);t++){
				for(int k=0; k<Coef.extent(thirdDim); k++){
					for(int j=0; j<Coef.extent(secondDim); j++){
						for(int m=0; m<Coef.extent(firstDim); m++){    
							if(abs(Coef(m,j,k,t,e))<thr(i)){					// this must be fixed, since it is never true
								tmp1+=(float)abs(Coef(m,j,k,t,e))/(float)Coef.numElements();	
								tmp2+=(float)pow(abs(Coef(m,j,k,t,e)),2)/(float)Coef.numElements();
							}
							if(abs(Coef(m,j,k,t,e))>=thr(i)) {
								tmp2+=(float)pow(thr(i),2)/(float)Coef.numElements();
							}
						}}}}}
		thc(i) = pow(noise,2) - 2.0*pow(noise,2)*tmp1 + tmp2;
	}

	threshold = min(thc);

}

void THRESHOLD::exec_multibandthreshold(Array< complex<float>,5 >&Coef) {

	int sx = Coef.extent(firstDim); 
	int sy = Coef.extent(secondDim);
	int sz = Coef.extent(thirdDim);
	int st = Coef.extent(fourthDim);
	int tt = 0;

	if(temporal==true) tt=1;

	if(threshold_type==TH_BAYES) {
		cout<<"Bayes Shrink ";
		if(soft==true)	{ cout<<"(soft)"<<endl; }
		else { cout<<"(hard)"<<endl; }
	} else if(threshold_type==TH_SURE) {
		cout<<"SURE Shrink ";
		if(soft==true)	{ cout<<"(soft)"<<endl; }
		else { cout<<"(hard)"<<endl; }
	} else {
		return;
	}

	// Calculation of subbands borders
	Array<int,3>xx;
	xx.setStorage(ColumnMajorArray<3>());
	xx.resize(waveL,8,2);
	Array<int,3>yy;
	yy.resize(waveL,8,2);
	yy.setStorage(ColumnMajorArray<3>());
	Array<int,3>zz;
	zz.resize(waveL,8,2);
	zz.setStorage(ColumnMajorArray<3>());

	// To get noise from HHH band
	Array<complex<float>,5> HHHref = Coef(Range(sx/2,sx-1),Range(sy/2,sy-1),Range(sz/2,sz-1),Range(tt,st-1),Range::all());
	thresh = 0.5;
	get_threshold(HHHref);
	noise = threshold/0.6745;

	if(VERBOSE) cout<<"Estimated noise: "<<noise<<endl;

	// Wavelet band borders (assumes 3D-5D data)
	int dx = sx, dy = sy, dz = sz;
	for(int l=0;l<waveL;l++) {
		dx=dx/2; dy=dy/2; dz=dz/2;
		xx(l,0,0)=dx; xx(l,0,1)=dx*2-1; yy(l,0,0)=0;  yy(l,0,1)=dy-1;   zz(l,0,0)=0;   zz(l,0,1)=dz-1;
		xx(l,1,0)=dx; xx(l,1,1)=dx*2-1; yy(l,1,0)=dy; yy(l,1,1)=dy*2-1; zz(l,1,0)=0;   zz(l,1,1)=dz-1;
		xx(l,2,0)=0;  xx(l,2,1)=dx-1;   yy(l,2,0)=dy; yy(l,2,1)=dy*2-1; zz(l,2,0)=0;   zz(l,2,1)=dz-1;
		xx(l,3,0)=dx; xx(l,3,1)=dx*2-1; yy(l,3,0)=0;  yy(l,3,1)=dy-1;   zz(l,3,0)=dz;  zz(l,3,1)=dz*2-1;
		xx(l,4,0)=dx; xx(l,4,1)=dx*2-1; yy(l,4,0)=dy; yy(l,4,1)=dy*2-1; zz(l,4,0)=dz;  zz(l,4,1)=dz*2-1;
		xx(l,5,0)=0;  xx(l,5,1)=dx-1;   yy(l,5,0)=dy; yy(l,5,1)=dy*2-1; zz(l,5,0)=dz;  zz(l,5,1)=dz*2-1;
		xx(l,6,0)=0;  xx(l,6,1)=dx-1;   yy(l,6,0)=0;  yy(l,6,1)=dy-1;   zz(l,6,0)=dz;  zz(l,6,1)=dz*2-1;
		// LLL
		if(l==(waveL-1)) {
			xx(l,7,0)=0; xx(l,7,1)=dx-1; yy(l,7,0)=0; yy(l,7,1)=dy-1; zz(l,7,0)=0; zz(l,7,1)=dz-1;
		}
	}

	// Get and execute thresholding
	for(int l=0;l<waveL;l++) {
		if(VERBOSE) cout<<"Decomposition level: "<<l+1<<endl;
		for(int k=0;k<=7;k++) {
			if( ( (k==7) && (l!=(waveL-1))) || ( (k==7) && (thapp==false)) ) {
				continue;
			}
			Array<complex<float>,5> Cfref = Coef(Range(xx(l,k,0),xx(l,k,1)),Range(yy(l,k,0),yy(l,k,1)),Range(zz(l,k,0),zz(l,k,1)),Range(tt,st-1),Range::all());
			switch(threshold_type) {
				case TH_BAYES:{
						      get_bayesthreshold(Cfref);
					      }break;
				case TH_SURE:{
						     get_surethreshold(Cfref);
					     }break;
			}
			thresholding(Cfref);
		}
	}

	// Do thresholding for all branches at once (for testing only)
	// Based on HHH1 branch only
	// Array<complex<float>,5> Cf = Coef(Range(sx/2,sx-1),Range(sy/2,sy-1),Range(sz/2,sz-1),Range::all(),Range::all());
	// get_bayesthreshold(Cf);
	// thresholding(Coef);
}

void THRESHOLD::exec_visuthreshold(Array< complex<float>,5>&Coef) {

	int tt = 0;
	int st = Coef.extent(fourthDim);

	if(temporal==true) tt=1;

	cout<<"Visu Shrink ";
	if(soft==true)	{ cout<<"(soft)"<<endl; }
	else { cout<<"(hard)"<<endl; }

	Array<complex<float>,5> Cf = Coef(Range::all(),Range::all(),Range::all(),Range(tt,st-1),Range::all());
	get_visuthreshold(Cf);
	thresholding(Cf);

}

void THRESHOLD::exec_fractionthreshold(Array< complex<float>,5>&Coef) {
	// This makes no sense ? if(temporal==true) tt=1;

	cout<<"Fractional Thresholding ";
	if(soft==true)	{ cout<<"(soft)"<<endl; }
	else { cout<<"(hard)"<<endl; }

	get_threshold(Coef);
	thresholding(Coef);

}

// Option of thresh of appriximation wavelets or all wavelets, hard vs soft threshold
void THRESHOLD::thresholding(Array< complex<float>,5 >&Coef){

	if(VERBOSE) cout << "Determined thresh = " << threshold<<endl;

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
						}else if(soft==true){
							float theta = arg(Coef(i,j,k,t,e));
							Coef(i,j,k,t,e) = (abs(Coef(i,j,k,t,e))-threshold)*complex<float>(cos(theta),sin(theta));
						}
	}}}}} 
	if(VERBOSE) cout << "Removed " << (float)((float)count/(float)Coef.numElements()) << " of the values" << endl;
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

void THRESHOLD::fista_update(  Array< complex<float>,5 >&X,Array< complex<float>,5 >&X_old,int iteration){

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

