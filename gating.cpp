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
        
		vs_type = NONE;
        tornado_shape = VIPR; // Kr^2 shape
		kmax = 128; // TEMP
		gate_type = TIME;
		
        // Catch command line switches
#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);
#define char_flag(name,val)   }else if(strcmp(name,pstring[pos]) == 0){ pos++; strcpy(val,pstring[pos]);
 
        for(int pos=0; pos<numarg; pos++) {
	       
			if(strcmp("-viewshare_type",pstring[pos]) == 0) {
				pos++;
				if( pos==numarg){
					cout << "Please provide vieshare type (-h for usage)" << endl;
					exit(1);
				trig_flag(TORNADO,"tornado",vs_type);
				trig_flag(HIST_MODE,"hist",vs_type);
				trig_flag(NONE,"none",vs_type);
				
				}else{
					cout << "Please provide vieshare type (-h for usage)" << endl;
					exit(1);
				}
			}else if(strcmp("-gating_type",pstring[pos]) == 0) {
				pos++;
				if( pos==numarg){
					cout << "Please provide gating type (-h for usage)" << endl;
					exit(1);
				trig_flag(RESP,"resp",gate_type);
				trig_flag(ECG,"ecg",gate_type);
				trig_flag(TIME,"time",gate_type);
				trig_flag(PREP,"prep",gate_type);
								
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
		
		help_flag("-gating_type []","how to gate images");
		help_flag("","  resp = respiratory phases");
		help_flag("","  ecg = gate by cardiac");
		help_flag("","  time = bin by acquisition time");
		help_flag("","  prep = bin by time from prep pulses");
		
		cout << "Filter parameters for tornado:" << endl;
        help_flag("-vs_wdth_low []","width in the center of k-space in frames");
        help_flag("-vs_wdth_high []","width in the periphery of k-space in frames");
        help_flag("-vs_vipr_tornado","3D radial tornado (default)");
		help_flag("-vs_radial_tornado","2D radial tornado");
}

void GATING::init( const MRI_DATA& data,int frames){
	
	
	// Create Array and Fill with Base 
	gate_times.setStorage(ColumnMajorArray<3>());
	switch(gate_type){
		case(RESP):{
			gate_times.resize( data.resp.shape());				  
			gate_times = data.resp;
		}break;
		
		case(ECG):{
			gate_times.resize( data.ecg.shape());				  
			gate_times = data.ecg;
		}break;
		
		case(TIME):{
			gate_times.resize( data.time.shape());				  
			gate_times = data.time;
		}break;
		
		case(PREP):{
			gate_times.resize( data.time.shape());				  
			gate_times = data.prep;
		}break;
	}
				
	// Get Range
	float max_time =max(gate_times);
	float min_time =min(gate_times);
	
	// Rescale to Frames
	scale_time= frames/(max_time-min_time);
	offset_time = min_time;
	
	gate_times -= offset_time;
	gate_times *= scale_time;
	
	// Set time points
	gate_frames = new float[frames];
	for(int i=0; i < frames; i++){
		gate_frames[i] = 0.5+(float)i;
	}
	cout << "Time Range :: " << min_time << " to " << max_time << endl;

	if(vs_type==TORNADO && (wdth_low != wdth_high) ){
		switch(tornado_shape){
			case(VIPR):{
				kmax = max( data.kx*data.kx + data.ky*data.ky + data.kz*data.kz);
				kmax = sqrt(kmax);
			}break;
			
			case(RADIAL):{
				kmax = max( data.kx*data.kx + data.ky*data.ky );
				kmax = sqrt(kmax);
			}break;
			
			case(FLAT):{
				kmax =999;
			}
					}
		cout << "Kmax = " << kmax << endl;
	}
		
	if( vs_type== HIST_MODE){
		cout << "Sorting Data into Histogram" << endl;
		// Use Aradillo Sort function
		arma::fvec time_sort(gate_times.numElements());
		
		// Copy into array
		int count = 0;
		for(int e=0; e< gate_times.length(thirdDim); e++){
		 for(int slice=0; slice< gate_times.length(secondDim); slice++){
		  for(int view=0; view< gate_times.length(firstDim); view++){
		   time_sort(count) = gate_times(view,slice,e);
		   count++;
		}}}
		int Ncount = count;
		
		// Sort
		arma::uvec sort_temp = sort_index( time_sort );
		arma::uvec sort_idx = sort_index( sort_temp );
		
		// Now Split into frames
		count = 0;
		for(int e=0; e< gate_times.length(thirdDim); e++){
		 for(int slice=0; slice< gate_times.length(secondDim); slice++){
		  for(int view=0; view< gate_times.length(firstDim); view++){
		   int t_frame = (int)( (float)(sort_idx(count)*frames) / Ncount);  
		   gate_times(view,slice,e) = t_frame;
		   count++;
		}}}
	}	

}




void GATING::weight_data(Array<float,3>&Tw, int e, const Array<float,3> &kx, const Array<float,3> &ky,const Array<float,3> &kz,int t,WeightType w_type){
        
	switch(vs_type){
		case(TORNADO):{
			tornado_weight(Tw,e, kx, ky,kz, t,w_type);
     	}break;
		
		case(HIST_MODE):
		case(NONE):{
			hist_weight( Tw, e,t);
		}break;
	}
}



/*
Simple 1 to 1 frames. No sharing
*/

void GATING::hist_weight( Array<float,3>&Tw,int e, int t){
	
	for(int k=0; k<Tw.length(thirdDim); k++){
	for(int j=0; j<Tw.length(secondDim); j++){
	for(int i=0; i<Tw.length(firstDim); i++){
		// Get K-space Radius
		Tw(i,j,k) *= ( floor(gate_times(j,k,e) == t) ) ? ( 1.0 ) : ( 0.0 );
	}}}
}

/* Tornado-like filter in k-space in rcframe units
  ________________________________________________
              /wdth_low\
            /           \
          /              \
        /                 \
     /                     \
    <-----wdth_high-------->
	
	KMJ: Rewrote entirely. Otherwise would be wrong for all but center out without ramp sampling/variable density!
*/
void GATING::tornado_weight(Array<float,3>&Tw, int e, const Array<float,3> &kx, const Array<float,3> &ky,const Array<float,3> &kz,int t,WeightType w_type){
        
	float current_time = gate_frames[t];
	
	for(int k=0; k<Tw.length(thirdDim); k++){
	for(int j=0; j<Tw.length(secondDim); j++){
	for(int i=0; i<Tw.length(firstDim); i++){
		
		float t_diff = abs( gate_times(j,k,e) - current_time );
	
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
		if(w_type == ITERATIVE){
			Tw(i,j,k) *= ( t_diff < wdth)? ( 1.0) : ( 0.0); // Don't Divide
		}else{
			Tw(i,j,k) *= ( t_diff < wdth)? ( 1./wdth) : ( 0.0);
		}
				
	}}}
}

