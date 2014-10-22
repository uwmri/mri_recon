/************************************************
View sharing techniques

The class uses information from Times.dat file to create
a mask used during reconstruction of each frame in recon.cxx

Initial Author:
        Grzegorz Bauman (gbauman@wisc.edu)

Changelog: 
        Stan Kruger (sjkruger@wisc.edu) 130107
        tornado filter should now be more robust.  Bins should be much closer to the same number of projections, and can now overlap as desired depending on < "frames," "vs_a," and "vs_b" >
		init
		2013-03-08  KMJ: Major changes. Renamed many variables to practical names. Fixed tornado filer for non equidistant spacing. Added 
		respiratory gating. Better commenting,etc.
	
Init:
        recon_binary -f data_header.txt -rcframes 32 -vs tornado -vs_a 1 -vs_b 5 -vs_shape 2
        GATING vs(argc,**argv);
        vs.createmask(TimeWeight,timesE,t);


*************************************************/

#include "gating.h"
#include "io_templates.hpp"
using namespace NDarray;

GATING::GATING(){
}

// Setup of 
GATING::GATING( int numarg, char **pstring) {

        // Setting default values, configurable
        wdth_low  = 1;
        wdth_high = 4;
        
		vs_type = NONE;
        	tornado_shape = VIPR; // Kr^2 shape
		kmax = 128; // TEMP
		gate_type = GATE_NONE;;
		
		// Respiratory Efficiency
		correct_resp_drift = 0;
		resp_gate_efficiency = 0.5;
        	resp_gate_type = RESP_NONE;
		resp_gate_signal = BELLOWS;
		
		external_weights = 0;
		strcpy(external_weights_filename,"");

		// Catch command line switches
#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);
#define char_flag(name,val)   }else if(strcmp(name,pstring[pos]) == 0){ pos++; strcpy(val,pstring[pos]);
 
        for(int pos=0; pos<numarg; pos++) {
	       
			if(strcmp("-viewshare_type",pstring[pos]) == 0) {
				pos++;
				if( pos==numarg){
					cout << "Please provide vieshare type..none/hist/tornado (-h for usage)" << endl;
					exit(1);
				trig_flag(TORNADO,"tornado",vs_type);
				trig_flag(HIST_MODE,"hist",vs_type);
				trig_flag(NONE,"none",vs_type);
				
				}else{
					cout << "Please provide vieshare type..none/hist/tornado (-h for usage)" << endl;
					exit(1);
				}
			}else if(strcmp("-gating_type",pstring[pos]) == 0) {
				pos++;
				if( pos==numarg){
					cout << "Please provide gating type..ecg/retro_ecg/resp/time/prep (-h for usage)" << endl;
					exit(1);
				trig_flag(RETRO_ECG,"retro_ecg",gate_type);
				trig_flag(RESP,"resp",gate_type);
				trig_flag(ECG,"ecg",gate_type);
				trig_flag(TIME,"time",gate_type);
				trig_flag(PREP,"prep",gate_type);
								
				}else{
					cout << "Please provide gating type..ecg/retro_ecg/resp/time/prep (-h for usage)" << endl;
					exit(1);
				}
			}else if(strcmp("-resp_gate",pstring[pos]) == 0) {
				pos++;
				if( pos==numarg){
					cout << "Please provide respiratory gating type..thresh/weight (-h for usage)" << endl;
					exit(1);
				trig_flag(RESP_THRESH,"thresh",resp_gate_type);
				trig_flag(RESP_WEIGHT,"weight",resp_gate_type);
								
				}else{
					cout << "Please provide respiratory gating type..thresh/weight (-h for usage)" << endl;
					exit(1);
				}			
			}else if(strcmp("-resp_gate_signal",pstring[pos]) == 0) {
				pos++;
				if( pos==numarg){
					cout << "Please provide a data source for estimation of respiratory phase..bellows/internal (-h for usage)" << endl;
					exit(1);
				trig_flag(BELLOWS,"bellows",resp_gate_signal);
				trig_flag(DC_DATA ,"internal",resp_gate_signal);
								
				}else{
					cout << "Please provide a data source for estimation of respiratory phase..bellows/internal (-h for usage)" << endl;
					exit(1);
				}	
			int_flag("-vs_wdth_low",wdth_low);
            		int_flag("-vs_wdth_high",wdth_high);
     			trig_flag(1,"-correct_resp_drift",correct_resp_drift);
			float_flag("-resp_gate_efficiency",resp_gate_efficiency);
			trig_flag(1,"-external_weights",external_weights);
			char_flag("-external_weights_file",external_weights_filename);
			}
	}

	if (resp_gate_type == RESP_THRESH) {
		cout << "Using threshold based respiratory gating" << endl;
	} else if (resp_gate_type == RESP_WEIGHT) {
		cout << "Using (fuzzy) weight based respiratory gating" << endl;
	}

	if ((external_weights == 1) && (strcmp("",external_weights_filename) == 0)) {
		cout << "No external fuzzy weight file specified.  Specify with '-external_weights_file []' option (-h for usage)" << endl;
		exit(1);
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
		help_flag("","  retro_ecg = retrospective gate by cardiac");
		help_flag("","  time = bin by acquisition time");
		help_flag("","  prep = bin by time from prep pulses");
		
		help_flag("-resp_gate []","In addition to other gating, perform respiratory gating");
		help_flag("","  thresh = threshold values");
		help_flag("","  weight = downweight bad values (see Johnson et al. MRM 67(6):1600");

		help_flag("-resp_gate_signal","Specify source for the data used to estimate respiratory phase");
		help_flag("","  bellows = signal from respiratory bellows belt in gating file (default)");
		help_flag("","  internal = use dc/low spatial frequency information extracted from acquired data");
				
		cout << "Filter parameters for tornado:" << endl;
        help_flag("-vs_wdth_low []","width in the center of k-space in frames");
        help_flag("-vs_wdth_high []","width in the periphery of k-space in frames");
        help_flag("-vs_vipr_tornado","3D radial tornado (default)");
		help_flag("-vs_radial_tornado","2D radial tornado");
		
		cout << "Control for Resp Data" << endl;
		help_flag("-correct_resp_drift","Median filter with 10s interval");
		help_flag("-resp_gate_efficiency","Fraction of data to accept");
		
		cout << "Control for ECG Data" << endl;
		help_flag("-bad_ecg_filter","Filter Bad ECG Vals (>10,000ms)");
		
		cout << "Control for external fuzzy weighting file (more generic functionality similar to '-resp_gate weight' option above)" << endl;
		help_flag("-external_weights","Read fuzzy weights from external file");
		help_flag("-external_weights_file []","Specify path to file with fuzzy weights (see Johnson et. al MRM 67(6):1600");
		
}


/*----------------------------------------------
     Smooths the Resp and subtracts off to correct drift
 *----------------------------------------------*/ 
void GATING::filter_resp(  const MRI_DATA &data ){
	
		cout << "Time range = " << ( max(data.time)-min(data.time) ) << endl;
		int fsize = (int)( 5.0 / (  ( max(data.time)-min(data.time) ) / data.time.numElements() ) ); // 10s filter / delta time
	
		cout << "Sorting Gate Data by Acquisition Time" << endl;
		
		// Use Aradillo Sort function
		arma::fvec time(gate_times.numElements());
		arma::fvec resp(gate_times.numElements());
		arma::fvec time_linear_resp(gate_times.numElements());
		
		// Put into Matrix for Armadillo
		int count = 0;
		for(int e=0; e< gate_times.length(thirdDim); e++){
		 for(int slice=0; slice< gate_times.length(secondDim); slice++){
		  for(int view=0; view< gate_times.length(firstDim); view++){
		   time(count) = data.time(view,slice,e);
		   resp(count) = data.resp(view,slice,e);
		   count++;
		}}}
		
		// Sort
		arma::uvec idx = arma::sort_index(time); 
		
		// Copy Resp
		idx.save("Sorted.dat", arma::raw_ascii);
		for(int i=0; i< (int)gate_times.numElements(); i++){
			time_linear_resp( i )= resp( idx(i));
		}
		time_linear_resp.save("RSorted.dat",arma::raw_ascii);
		
		
		// Now Filter
		cout << "Filtering Resp Data by" << fsize << endl;
		for(int i=0; i< (int)gate_times.numElements(); i++){
			int start = i - fsize;
			int stop  = i + fsize;
			if(start < 0){
				stop  = 2*fsize;
				start = 0; 
			}
			
			if(stop >=  (int)gate_times.numElements()){
				stop  = gate_times.numElements() -1;
				start = gate_times.numElements() - 1 - 2*fsize;
			}
			
			float thresh = median(  time_linear_resp.rows( start,stop) );
			resp(idx(i)) -= thresh;
			
		}
		resp.save("RFiltered.dat",arma::raw_ascii);
		
		// Copy Resp
		for(int i=0; i< (int)gate_times.numElements(); i++){
			time_linear_resp( i )= resp( idx(i));
		}
		time_linear_resp.save("RSorted_Filter.dat",arma::raw_ascii);
						
		// Write Back to Blitz
		count = 0;
		for(int e=0; e< gate_times.length(thirdDim); e++){
		 for(int slice=0; slice< gate_times.length(secondDim); slice++){
		  for(int view=0; view< gate_times.length(firstDim); view++){
		   gate_times(view,slice,e)=resp(count);
		   count++;
		}}}
		
		
}
void GATING::init( const MRI_DATA& data,int frames){
	
	init_resp_gating(data,frames);
	init_time_resolved(data,frames);
}

float GATING::temporal_resolution(void){
	return( actual_temporal_resolution );
}

void GATING::init_resp_gating( const MRI_DATA& data,int frames){
	
	cout << "Initializing Respiratory Gating with " << frames << " frames " << endl;
	

	if (resp_gate_type != RESP_NONE) {
		
		switch (resp_gate_signal) {
			case(BELLOWS):{
				     cout << "Respiratory gating using bellows belt waveform" << endl;
						resp_weight.resize(data.resp.shape());
				      }break;
			case(DC_DATA):{

					cout << "Respiratory gating using internal DC waveform(s)" << endl;
						resp_weight.resize(data.resp.shape());

						int views_per_grid = 16;
						int Ngrids = (int)floor((float)data.Num_Readouts / (float)views_per_grid);
						Array< Array<complex<float>, 3>, 2 >Kdc = Alloc5DContainer< complex<float> >(4,4,4,Ngrids,data.Num_Coils);
						
						extract_dc_data(Kdc,data);
						for (int coil = 0; coil < data.Num_Coils; coil++) {
	
						//	for( Array<Array<complex<float>,3>,1>::iterator miter=TempK.begin(); miter!=TempK.end(); miter++){
						
							for (int g = 0; g < Ngrids; g++) {
								Array<complex<float>, 3> Kgrid = Kdc(g,coil);
								ArrayWriteMagAppend(Kgrid,"dcMagnitude.dat");
								ArrayWritePhaseAppend(Kgrid,"dcPhase.dat");
							}
						}

						/*  The question is now how to process/convert this data for the gating algorithm below?
						 *  For now, just use the vector magnitude over all coils of k = 0, then apply a moving
						 *  average to upsample to full temporal resolution					*/


						int view_start = 0;
						int view_stop = views_per_grid-1;

						cout << "views_per_grid = " << views_per_grid << endl;

						int kzero_x = Kdc(0,0).length(firstDim)/2 + 1;
						int kzero_y = Kdc(0,0).length(secondDim)/2 + 1;
						int kzero_z = Kdc(0,0).length(thirdDim)/2 + 1;

						cout << kzero_x << "  " << kzero_y << "   " << kzero_z << endl;
						Array<float,1>temp_weight(data.Num_Readouts, ColumnMajorArray<1>());
						temp_weight = 0.0;

						cout << "Compress low frq grids to metric (frob norm for now)" << endl;

						#pragma omp parallel for
						for (int grid = 0; grid < Ngrids; grid++){
							float dc_resp = 0.0;
							for (int coil = 0; coil < data.Num_Coils; coil++) {
								dc_resp += norm(Kdc(grid,coil)(kzero_x,kzero_y,kzero_z));
							}
							dc_resp = sqrt(dc_resp);
							temp_weight(Range(grid*views_per_grid,(grid+1)*views_per_grid-1)) = dc_resp;
						}

						cout << "Marker" << endl;
						/* Get the acquisition order indices */
						arma::fvec time(data.time.numElements());
						// Put into Matrix for Armadillo
						int count = 0;
						for(int e=0; e< data.time.length(thirdDim); e++){
							for(int slice=0; slice< data.time.length(secondDim); slice++){
								for(int view=0; view< data.time.length(firstDim); view++){
									time(count) = data.time(view,slice,e);
									count++;
								}}}

						// Sort
						arma::uvec idx = arma::sort_index(time); 

						cout << "Moving average upsampling" << endl;

						/*  Now compute moving average, with a window 2x the views per grid, store the result in the original order */

						for (int pos = 0; pos < data.Num_Readouts; pos++) {
							if ((pos - views_per_grid) < 0){
								int numel = pos+views_per_grid;
								data.resp(idx(pos),Range::all(),Range::all()) = sum(temp_weight(Range(0,numel-1)))/(float)numel;
							
							} else if ((pos + views_per_grid - 1) >= data.Num_Readouts) {
								int numel = views_per_grid + (data.Num_Readouts - pos);
								data.resp(idx(pos),Range::all(),Range::all()) = sum(temp_weight(Range(pos-views_per_grid,data.Num_Readouts-1)))/(float)numel;
							
							} else {
								data.resp(idx(pos),Range::all(),Range::all()) = sum(temp_weight(Range(pos-views_per_grid,pos+views_per_grid-1)))/(2.0*(float)views_per_grid);
							}
						}
					  	cout << "DC signal processing complete" << endl;
					}break;
			default:{
						cout << "Invalid option for respiratory gating waveform" << endl;
				}break;
		}

	}


	switch (resp_gate_type){

		case(RESP_THRESH):{

					  cout << "Copying Resp Waveform" << endl;
					  resp_weight.resize( data.resp.shape());				  
					  resp_weight = data.resp;

					  cout << "Time Sorting Data" << endl;

					  // Use Aradillo Sort function
					  arma::fvec time(resp_weight.numElements());
					  arma::fvec resp(resp_weight.numElements());
					  arma::fvec arma_resp_weight(resp_weight.numElements());

					  arma::fvec time_linear_resp(resp_weight.numElements());
					  arma::fvec time_sort_resp_weight(resp_weight.numElements());

					  // Put into Matrix for Armadillo
					  int count = 0;
					  for(int e=0; e< resp_weight.length(thirdDim); e++){
						  for(int slice=0; slice< resp_weight.length(secondDim); slice++){
							  for(int view=0; view< resp_weight.length(firstDim); view++){
								  time(count) = data.time(view,slice,e);
								  resp(count) = data.resp(view,slice,e);
								  count++;
							  }}}
					  time.save("Time.txt",arma::raw_ascii);
					  resp.save("Resp.txt",arma::raw_ascii);

					  // Sort
					  arma::uvec idx = arma::sort_index(time); 

					  // Copy Resp
					  idx.save("Sorted.dat", arma::raw_ascii);
					  for(int i=0; i< (int)resp_weight.numElements(); i++){
						  time_linear_resp( i )= resp( idx(i));
					  }
					  time_linear_resp.save("TimeResp.txt",arma::raw_ascii);


					  // Size of histogram
					  cout << "Time range = " << ( max(data.time)-min(data.time) ) << endl;
					  int fsize = (int)( 5.0 / (  ( max(data.time)-min(data.time) ) / data.time.numElements() ) ); // 10s filter / delta time


					  // Now Filter
					  cout << "Thresholding Data Frame Size = " << fsize << endl;
					  for(int i=0; i< (int)time_linear_resp.n_elem; i++){
						  int start = i - fsize;
						  int stop  = i + fsize;
						  if(start < 0){
							  stop  = 2*fsize;
							  start = 0; 
						  }

						  if(stop >=  (int)time_linear_resp.n_elem){
							  stop  = time_linear_resp.n_elem -1;
							  start = time_linear_resp.n_elem - 1 - 2*fsize;
						  }

						  arma::fvec temp = time_linear_resp.rows( start,stop);
						  arma::fvec temp2= sort(temp);
						  float thresh = temp2( (int)( (float)temp2.n_elem*( 1.0- resp_gate_efficiency )));

						  arma_resp_weight(idx(i))= ( time_linear_resp(i) > thresh ) ? ( 1.0 ) : ( 0.0);
						  time_sort_resp_weight(i ) =arma_resp_weight(idx(i));
					  }
					  time_sort_resp_weight.save("TimeWeight.txt",arma::raw_ascii);
					  arma_resp_weight.save("Weight.txt",arma::raw_ascii);


					  // Copy Back
					  
					  for(int e=0; e< resp_weight.length(thirdDim); e++){
						  for(int slice=0; slice< resp_weight.length(secondDim); slice++){
							  for(int view=0; view< resp_weight.length(firstDim); view++){
								  resp_weight(view,slice,e)=arma_resp_weight(count);
								  count++;
							  }}}


				  }break;
		case(RESP_WEIGHT):{


					  if (external_weights == 1) {
						  cout << "Reading fuzzy retrospective weights from file" << endl;
					  
						  ArrayRead(resp_weight,external_weights_filename);
					  
					  }else{

						  cout << "Copying Resp Waveform" << endl;
						  resp_weight = data.resp;


						  // Use Aradillo Sort function
						  arma::fvec time(resp_weight.numElements());
						  arma::fvec resp(resp_weight.numElements());
						  arma::fvec arma_resp_weight(resp_weight.numElements());

						  arma::fvec time_linear_resp(resp_weight.numElements());
						  arma::fvec time_sort_resp_weight(resp_weight.numElements());

						  // Put into Matrix for Armadillo
						  int count = 0;
						  for(int e=0; e< resp_weight.length(thirdDim); e++){
							  for(int slice=0; slice< resp_weight.length(secondDim); slice++){
								  for(int view=0; view< resp_weight.length(firstDim); view++){
									  time(count) = data.time(view,slice,e);
									  resp(count) = data.resp(view,slice,e);
									  count++;
								  }}}
						  time.save("Time.txt",arma::raw_ascii);
						  resp.save("Resp.txt",arma::raw_ascii);

						  // Sort
						  arma::uvec idx = arma::sort_index(time); 

						  // Copy Resp
						  idx.save("Sorted.dat", arma::raw_ascii);
						  for(int i=0; i< (int)resp_weight.numElements(); i++){
							  time_linear_resp( i )= resp( idx(i));
						  }
						  time_linear_resp.save("TimeResp.txt",arma::raw_ascii);


						  // Size of histogram
						  cout << "Time range = " << ( max(data.time)-min(data.time) ) << endl;
						  int fsize = (int)( 5.0 / (  ( max(data.time)-min(data.time) ) / data.time.numElements() ) ); // 10s filter / delta time


						  // Now Filter
						  cout << "Estimating median within " << resp_gate_efficiency*100 << "% efficiency window of first 10sec of data" << fsize << endl;
							
						  	  arma::fvec temp = time_linear_resp.rows( 100,100+2*fsize);
							  arma::fvec temp2= sort(temp);
							  float med_resp = temp2( (int)( (float)temp2.n_elem*( 1.0- resp_gate_efficiency/2.0 )));

							  float sigma = temp2(temp2.n_elem-1) - temp2((int)((float)temp2.n_elem*(1.0 - resp_gate_efficiency )));


							  //arma_resp_weight(idx(i))= ( time_linear_resp(i) > thresh ) ? ( 1.0 ) : ( 0.0);
							 
					  	for(int i=0; i< (int)time_linear_resp.n_elem; i++){
							  arma_resp_weight(idx(i))= ( 1.0 / (abs(med_resp - time_linear_resp(i)) + sigma));
							  time_sort_resp_weight(i ) =arma_resp_weight(idx(i));
						  }
						  time_sort_resp_weight.save("TimeWeight.txt",arma::raw_ascii);
						  arma_resp_weight.save("Weight.txt",arma::raw_ascii);


						  // Copy Back
						  count = 0;
						  for(int e=0; e< resp_weight.length(thirdDim); e++){
							  for(int slice=0; slice< resp_weight.length(secondDim); slice++){
								  for(int view=0; view< resp_weight.length(firstDim); view++){
									  resp_weight(view,slice,e)=arma_resp_weight(count);
									  count++;
								  }}}
					  }
					
				  }break;
		case(RESP_NONE):
		default:{
				return;
			}
	}


	ArrayWrite(resp_weight,"RespWeight.dat");
}

void GATING::init_time_resolved( const MRI_DATA& data,int frames){
	
	cout << "Initializing Time resolved for" << frames << " frames" << endl;
	
	// Don't run 
	if( frames < 2 ){
		gate_type = GATE_NONE;
		return;
	}
			
	// Create Array and Fill with Base 
	gate_times.setStorage(ColumnMajorArray<3>());
	switch(gate_type){
		case(RESP):{
			cout << "Using Resp gate" << endl;
			gate_times.resize( data.resp.shape());				  
			gate_times = data.resp;
			if( correct_resp_drift ==1){
				cout << "Correcting Drift" << endl;
				filter_resp( data );
			}
		}break;
		
		case(RETRO_ECG):
		case(ECG):{
			gate_times.resize( data.ecg.shape());				  
			gate_times = data.ecg;
		}break;
		
		case(TIME):{
			cout << "Using Time Resolved" << endl;
			gate_times.resize( data.time.shape());				  
			gate_times = data.time;
		}break;
		
		case(PREP):{
			gate_times.resize( data.time.shape());				  
			gate_times = data.prep;
		}break;
		
		default:{
			gate_type = GATE_NONE;
			return;
		}
	}
	
		
	// Get Range
	float max_time =max(gate_times);
	float min_time =min(gate_times);
	
	if(gate_type==RETRO_ECG){
		gate_times -= min_time;
		min_time = 0;
				
		// Use Median to set value
		arma::fvec temp(gate_times.numElements());
		int count=0;
		for( Array<float,3>::iterator miter=gate_times.begin(); miter!=gate_times.end(); miter++,count++){
			temp(count) = *miter;
		}
		max_time = 2.0*median(temp);
		cout << "Retro ECG::RR is estimated to be " << max_time << endl;
	}
		
	// Rescale to Frames
	scale_time= (frames)/(max_time-min_time)*(1+1e-9); // Extra factor is to map last point to < frames
	offset_time = min_time;
	
	// Temporal resolution
	actual_temporal_resolution = ( max_time -min_time ) / frames;
	
	cout << "Actual temporal resolution = " << actual_temporal_resolution << endl;
	cout << " Gate offset = " << offset_time << endl;
	cout << " Gate scale = " << scale_time << endl;
	
	gate_times -= offset_time;
	gate_times *= scale_time;
	
	cout << "Time Range :: " << min_time << " to " << max_time << endl;
	
	/* Histogram*/
	{
		arma::vec temp(frames);
		temp.fill(0);
		for( Array<float,3>::iterator miter=gate_times.begin(); miter!=gate_times.end(); miter++){
		
			int pos = (int)floor( *miter);
			if( (pos < frames) && (pos >= 0) ){
		 		temp(pos)++;
			}
		}
		
		// Export
		cout << "Values per frames" << endl;
		for(int i=0; i < frames; i++){
			cout << " Frame " << i << " ,count = " << temp(i) << endl;
		}
		
	}
	
	
	cout << "Setting up view share" << endl;
	switch(vs_type){

	case(TORNADO ):{
		// Set time points
		gate_frames = new float[frames];
		for(int i=0; i < frames; i++){
			gate_frames[i] = 0.5+(float)i;
		}
		switch(tornado_shape){
			case(VIPR):{
				kmax = max( data.kx(0)*data.kx(0) + data.ky(0)*data.ky(0) + data.kz(0)*data.kz(0));
				kmax = sqrt(kmax);
			}break;
			
			case(RADIAL):{
				kmax = max( data.kx(0)*data.kx(0) + data.ky(0)*data.ky(0) );
				kmax = sqrt(kmax);
			}break;
			
			case(FLAT):{
				kmax =999;
			}
					}
		cout << "Kmax = " << kmax << endl;
	}break;
	
	case(NONE):{
		gate_times = floor(gate_times);
		cout << "Max gate time = " << min(gate_times) << endl;
		cout << "Min gate time = " << max(gate_times) << endl;
	}break;
			
	case(HIST_MODE):{
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
	}break;	
  }//Switch
  
}


void GATING::extract_dc_data(Array<Array<complex<float>, 3>,2> &K_dc_data, const MRI_DATA &data){

	Range all = Range::all();

	/* Grid the low frequency kspace samples to a grid for easier analysis and handling  */
	int Ngrids = K_dc_data.length(firstDim);
	int views_per_grid = 16;

	// Setup Gridding + FFT Structure
	gridFFT dc_gridding;
	dc_gridding.kernel_type = KAISER_KERNEL;
	dc_gridding.dwinX = 2.5;
	dc_gridding.dwinY = 2.5;
	dc_gridding.dwinZ = 2.5;
	dc_gridding.overgrid = 2;
	dc_gridding.precalc_gridding(4,4,4,3);

	// Weighting Array for Time coding
	Array< float, 3 >Kweight(data.Num_Pts,views_per_grid,1,ColumnMajorArray<3>());
	Kweight = 0.0;
	Array< float, 3 > KX(Kweight.shape(),ColumnMajorArray<3>());
	Array< float, 3 > KY(Kweight.shape(),ColumnMajorArray<3>());
	Array< float, 3 > KZ(Kweight.shape(),ColumnMajorArray<3>());
	
	Array<Array<complex<float>, 3>, 1> kdataC = Alloc4DContainer< complex<float> >(data.Num_Pts,views_per_grid,data.kdata.length(thirdDim),data.Num_Coils);


	// Use Aradillo Sort function
	arma::fvec time(resp_weight.numElements());

	// Put into Matrix for Armadillo
	int count = 0;
	for(int e=0; e< resp_weight.length(thirdDim); e++){
		for(int slice=0; slice< resp_weight.length(secondDim); slice++){
			for(int view=0; view< resp_weight.length(firstDim); view++){
				time(count) = data.time(view,slice,e);
				count++;
			}}}
					
	time.save("Time.txt",arma::raw_ascii);

	// Sort indices (acquisition order))
	arma::uvec idx = arma::sort_index(time); 				  
	idx.save("Sorted.dat", arma::raw_ascii);
	
	cout << "DC gridding... ";
	/* Note that the exported low frequency kspace grids are stored in time sorted order */
	for (int d = 0; d < Ngrids; d++) {
		
		int view_offset = d*views_per_grid;
		Kweight = 0.0;

		for (int v = 0; v < views_per_grid; v++) {
			Kweight(all,v,all) = data.kw(0)(all,idx(view_offset + v),all);
			KX(all,v,all) = data.kx(0)(all,idx(view_offset + v),all);
			KY(all,v,all) = data.ky(0)(all,idx(view_offset + v),all);
			KZ(all,v,all) = data.kz(0)(all,idx(view_offset + v),all);
			for (int c = 0; c < data.Num_Coils; c++) {
				kdataC(c)(all,v,all) = data.kdata(0,c)(all,idx(view_offset + v),all);
			}
		}
	
		for (int c = 0; c < data.Num_Coils; c++) {

			Array<complex<float>, 3> Kref = K_dc_data(d,c);
			Array<complex<float>, 3> Kdat = kdataC(c);

			dc_gridding.forward(Kref,Kdat,KX,KY,KZ,Kweight);
			ifft(Kref);
		}
	}
	cout << "Complete" << endl;
}


void GATING::weight_data(Array<float,3>&Tw, int e, const Array<float,3> &kx, const Array<float,3> &ky,const Array<float,3> &kz,int t,WeightType w_type,FrameType comp_type){
    
	switch(resp_gate_type){
		
		case(RESP_WEIGHT):{
			for(int k=0; k<Tw.length(thirdDim); k++){
			for(int j=0; j<Tw.length(secondDim); j++){
			for(int i=0; i<Tw.length(firstDim); i++){
				Tw(i,j,k) *= resp_weight(j,k,e);
			}}}
		}break;
		
		case(RESP_THRESH):{
			for(int k=0; k<Tw.length(thirdDim); k++){
			for(int j=0; j<Tw.length(secondDim); j++){
			for(int i=0; i<Tw.length(firstDim); i++){
				Tw(i,j,k) *= resp_weight(j,k,e);
			}}}
		}break;
		
		default:{
			
		}
	}
	
	
	
	if( (gate_type!= GATE_NONE) && (comp_type !=COMPOSITE) ){

		switch(vs_type){
			case(TORNADO):	{
						tornado_weight(Tw,e,kx,ky,kz,t,w_type);
					}break;

			case(HIST_MODE):
			case(NONE):	{		 
						hist_weight(Tw,e,t);
					}break;
		}

	}

	// Normalize Weighting
	//float sum_Tw = sum(Tw);
	//cout << "Sum Time Weight = " << sum_Tw << endl;
	//Tw /= sum_Tw;
}


/*
Simple 1 to 1 frames. No sharing
*/

void GATING::hist_weight( Array<float,3>&Tw,int e, int t){
	
	for(int k=0; k<Tw.length(thirdDim); k++){
	for(int j=0; j<Tw.length(secondDim); j++){
	for(int i=0; i<Tw.length(firstDim); i++){
		// Get K-space Radius
		Tw(i,j,k) *= ( floor(gate_times(j,k,e)) == t ) ? ( 1.0 ) : ( 0.0 );
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
	ArrayWrite(Tw,"TimeWeight.dat");
}

