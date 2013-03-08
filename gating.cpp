/************************************************
View sharing techniques

The class uses information from Times.dat file to create
a mask used during reconstruction of each frame in recon.cxx

Initial Author:
        Grzegorz Bauman (gbauman@wisc.edu)

Changelog: 
        Stan Kruger (sjkruger@wisc.edu) 130107
        tornado filter should now be more robust.  Bins should be much closer to the same number of projections, and can now overlap as desired depending on < "frames," "vs_a," and "vs_b" >
		
		2013-03-08  KMJ: Major changes. Renamed many variables to practical names. Fixed tornado filer for non equidistant spacing. Added 
		respiratory gating. Better commenting,etc.
	
Init:
        recon_binary -f data_header.txt -rcframes 32 -vs tornado -vs_a 1 -vs_b 5 -vs_shape 2
        GATING vs(argc,**argv);
        vs.createmask(TimeWeight,timesE,t);


*************************************************/

#include "gating.h"

// Setup of 
GATING::GATING( int numarg, char **pstring) {

        // Setting default values, configurable
        wdth_low  = 1;
        wdth_high = 4;
        
		type = NONE;
        tornado_shape = VIPR; // Kr^2 shape
		kmax = 128; // TEMP
		
        // Catch command line switches
#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);
#define char_flag(name,val)   }else if(strcmp(name,pstring[pos]) == 0){ pos++; strcpy(val,pstring[pos]);
 
        for(int pos=0; pos<numarg; pos++) {
	       
			if(strcmp("-viewshare_type",pstring[pos]) == 0) {
				pos++;
				if( pos==numarg){
					cout << "Please provide gating type (-h for usage)" << endl;
					exit(1);
				trig_flag(TORNADO,"tornado",type);
				trig_flag(HIST_MODE,"hist",type);
				trig_flag(NONE,"none",type);
				
				}else{
					cout << "Please provide temporal transform type..none/dft/diff/pca" << endl;
					exit(1);
				}
			int_flag("-vs_wdth_low",wdth_low);
            int_flag("-vs_wdth_high",wdth_high);
     		}
		}
}

void GATING::help_message() {
        cout << "----------------------------------------------" << endl;
        cout << "   View Sharing Control " << endl;
        cout << "----------------------------------------------" << endl;
        cout << "Usage:" << endl;
        help_flag("-viewshare_type []","view sharing method tornado/none/hist (defult=none)");
		help_flag("","  tornado = variable width filter in kr");
		help_flag("","  none = images at equal time intervals, no sharing between frames");
		help_flag("","  hist = images with equal data points");
		
		cout << "Filter parameters for tornado:" << endl;
        help_flag("-vs_wdth_low []","width in the center of k-space in frames");
        help_flag("-vs_wdth_high []","width in the periphery of k-space in frames");
        help_flag("-vs_vipr_tornado","3D radial tornado (default)");
		help_flag("-vs_radial_tornado","2D radial tornado");
}

void GATING::init( Array<float,4>&times,int frames,const Array<float,4> &kx, const Array<float,4> &ky,const Array<float,4> &kz){
	
	// Get Range
	float max_time =max(times);
	float min_time =min(times);
	
	// Rescale to Frames
	scale_time= frames/(max_time-min_time);
	offset_time = min_time;
	
	times -= min_time;
	times *= scale_time;
	
	// Set time points
	gate_frames = new float[frames];
	for(int i=0; i < frames; i++){
		gate_frames[i] = 0.5+(float)i;
	}
	cout << "Time Range :: " << min_time << " to " << max_time << endl;


	if(type==TORNADO && (wdth_low != wdth_high) ){
		switch(tornado_shape){
			case(VIPR):{
				kmax = max( kx*kx + ky*ky + kz*kz);
				kmax = sqrt(kmax);
			}break;
			
			case(RADIAL):{
				kmax = max( kx*kx + ky*ky);
				kmax = sqrt(kmax);
			}break;
			
			case(FLAT):{
				kmax =999;
			}
					}
		cout << "Kmax = " << kmax << endl;
	}
		
	if( type== HIST_MODE){
		cout << "Sorting Data into Histogram" << endl;
		// Use Aradillo Sort function
		arma::fvec time_sort(times.length(fourthDim)*times.length(thirdDim)*times.length(secondDim));
		
		// Copy into array
		int count = 0;
		for(int e=0; e< times.length(fourthDim); e++){
		 for(int slice=0; slice< times.length(thirdDim); slice++){
		  for(int view=0; view< times.length(secondDim); view++){
		   time_sort(count) = times(0,view,slice,e);
		   count++;
		}}}
		int Ncount = count;
		
		// Sort
		arma::uvec sort_temp = sort_index( time_sort );
		arma::uvec sort_idx = sort_index( sort_temp );
		
		// Now Split into frames
		count = 0;
		for(int e=0; e< times.length(fourthDim); e++){
		 for(int slice=0; slice< times.length(thirdDim); slice++){
		  for(int view=0; view< times.length(secondDim); view++){
		   int t_frame = (int)( (float)(sort_idx(count)*frames) / Ncount);  
		   
		   for(int i=0; i< times.length(firstDim); i++){
		   	   times(i,view,slice,e) = t_frame;
		   }
		   count++;
		}}}
	}	

}




void GATING::weight_data(Array<float,3>&Tw, Array<float,3>&times, const Array<float,3> &kx, const Array<float,3> &ky,const Array<float,3> &kz,int t){
        
	switch(type){
		case(TORNADO):{
			tornado_weight(Tw,times, kx, ky,kz, t);
     	}break;
		
		case(HIST_MODE):
		case(NONE):{
			hist_weight( Tw, times,t);
		}break;
	}
}



/*
Simple 1 to 1 frames. No sharing
*/

void GATING::hist_weight( Array<float,3>&Tw, Array<float,3>&times, int t){
	Tw *= ( floor(times)==t);
}

/* Tornado-like filter in k-space in rcframe units
  ________________________________________________
              /wdth_low\
            /           \
          /   | |        \
        /     | |         \
     /    c   | |   c      \
    <----------b------------>
	
	KMJ: Rewrote entirely. Otherwise would be wrong for all but center out without ramp sampling/variable density!
*/
void GATING::tornado_weight(Array<float,3>&Tw, Array<float,3>&times, const Array<float,3> &kx, const Array<float,3> &ky,const Array<float,3> &kz,int t){
        
	float current_time = gate_frames[t];
	
	for(int k=0; k<Tw.length(thirdDim); k++){
	for(int j=0; j<Tw.length(secondDim); j++){
	for(int i=0; i<Tw.length(firstDim); i++){
		
		float t_diff = abs( times(i,j,k) - current_time );
	
		// Get K-space Radius
		float kr=0;
		float k_power=0;
		switch(tornado_shape){
			
			case(FLAT):{
				kr= 0.0;
				k_power=0.0;
			}break;
			
			case(RADIAL):{
				kr = sqrt( kx(i,j,k)*kx(i,j,k) + ky(i,j,k)*ky(i,j,k));
				k_power=1.0;
			}break;
		
			case(VIPR):{
				kr = sqrt( kx(i,j,k)*kx(i,j,k) + ky(i,j,k)*ky(i,j,k) + kz(i,j,k)*kz(i,j,k));
				k_power=2.0;	
			}break;
		}
		
		
		float wdth = 0.5*(  (wdth_high - wdth_low)*pow(kr/kmax,k_power) + wdth_low);
		Tw(i,j,k) *= ( t_diff < wdth)? ( 1./wdth) : ( 0.0); 		
	}}}
}

