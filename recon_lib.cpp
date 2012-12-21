#include "recon_lib.h"
#include "io_templates.cpp"

// ----------------------
//  Basic constructor (no args)
// ----------------------
RECON::RECON(void){
	set_defaults();
}

// ----------------------
//  Sets Default Recon Parameters
// ----------------------
void RECON::set_defaults( void){
	// Help Message for recon
	recon_type = SOS; 
	data_type = EXTERNAL;
	coil_combine_type = LOWRES;
	
	complex_diff = false;
	
	numrecv = 1;
	zero_fill = 1.0;
	zoom = 1.0;
	zoom_x = 1.0;
	zoom_y = 1.0;
	zoom_z = 1.0;
	  
	rcxres=-1;
	rcyres=-1;
	rczres=-1;
	rcframes=1;
	rcencodes=1;
	num_slices =1;    
	lp_frac=1.0;
	smap_res=16;
	
	acc = 1;
	compress_coils = 0.0;
	export_smaps = 0;
	max_iter = 50;		
}

// ----------------------
//  Constructor with Command Line Read
// ----------------------

RECON::RECON(int numarg, char **pstring){
	set_defaults();	
	
	// --------------------------
	// Help Messages for Commandline Inputs
	//   -Please add your own help message for new classes
	// --------------------------
	for(int pos=0;pos< numarg;pos++){
		if( (strcmp(pstring[pos],"-h")==0) || (strcmp(pstring[pos],"-help")==0) || (strcmp(pstring[pos],"--help")==0)){
			help_message();	
			gridFFT::help_message();
			SPIRIT::help_message();	
			THRESHOLD::help_message();
			PHANTOM::help_message();	
			exit(0);
		}
	}
	
	// Get Input Parameters
	parse_commandline(numarg,pstring);
}



// ----------------------
// Help Message
// ----------------------
void RECON::help_message(void){
	cout << "----------------------------------------------" << endl;
	cout << "   Basic Recon Control " << endl;
	cout << "----------------------------------------------" << endl;
	cout << "Usage:" << endl;
	cout << "   recon_binary -f header.txt [flags]" << endl;
	
	cout << "Recon Size:" << endl;
	help_flag("-rcxres []","matrix size in x");
	help_flag("-rcyres []","matrix size in y");
	help_flag("-rczres []","matrix size in z");
	help_flag("-rcframes []","reconstructed temporal frames");
	help_flag("-zoom_x []","zoom factor in x");
	help_flag("-zoom_y []","zoom factor in x");
	help_flag("-zoom_z []","zoom factor in x");
	
	cout << "Recon Types:" << endl;
	help_flag("-sos","sum of squares");
	help_flag("-pils","pils (coil combine with low resolution images)");
	help_flag("-ist","iterative soft thresholding");
	help_flag("-fista","fast iterative soft thresholding");
	help_flag("-complex_diff","Subtract first encode");
	
	cout << "Iterative Recon Control:" << endl;
	help_flag("-max_iter []","max iterations for iterative recons");
    
	cout << "Coil Control:" << endl;
	help_flag("-espirit","use ESPIRIT to get coil sensitivies");
	help_flag("-coil_lowres","default,use low resolution images to get coil sensitivies");
	help_flag("-export_smaps","write sensitivity maps");
}

// --------------------
//  Read command line and set variables
// --------------------
void RECON::parse_commandline(int numarg, char **pstring){

#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);
#define char_flag(name,val)   }else if(strcmp(name,pstring[pos]) == 0){ pos++; strcpy(val,pstring[pos]);
  	  
  for(int pos=0; pos < numarg; pos++){
 	
  	if (strcmp("-h", pstring[pos] ) == 0) {
		char_flag("-f",filename);
		
		// Reconstruction Geometry
		int_flag("-rcxres",rcxres);
		int_flag("-rcyres",rcyres);
		int_flag("-rczres",rczres);
		int_flag("-rcframes",rcframes);
		
		float_flag("-zoom",zoom);
		float_flag("-zoom_x",zoom_x);
		float_flag("-zoom_y",zoom_y);
		float_flag("-zoom_z",zoom_z);
		
		// Type of Recons		
		trig_flag(SOS,"-sos",recon_type);
		trig_flag(CG,"-isense",recon_type);
		trig_flag(PILS,"-pils",recon_type);
		trig_flag(IST,"-ist",recon_type);
		trig_flag(FISTA,"-fista",recon_type);
		
		trig_flag(ESPIRIT,"-espirit",coil_combine_type);
		trig_flag(LOWRES,"-coil_lowres",coil_combine_type);
		trig_flag(1,"-export_smaps",export_smaps);
		
		// Source of data
		trig_flag(EXTERNAL,"-external_data",data_type);
		trig_flag(PFILE,"-pfile",data_type);
		trig_flag(PHANTOM,"-phantom",data_type);
		trig_flag(SIMULATE,"-simulate",data_type);
				
		// Data modification
		int_flag("-acc",acc);
		float_flag("-compress_coils",compress_coils);
		trig_flag(true,"-complex_diff",complex_diff);
		
		// Coil Combination + Resolution		
		float_flag("-lp_frac",lp_frac);
		float_flag("-smap_res",smap_res);
		
		// Iterations for IST
		int_flag("-max_iter",max_iter);
		
	}
  }
} 

//--------------------------------------------------
//  Read external header
//--------------------------------------------------

void RECON::parse_external_header(void){
	
	char parameter[80];
	float value;
	float value2;
	float value3;
	float value4;
	FILE *fid;
	char line[200];
	
	cout << "Reading External Header: " << endl;
	
	fid = fopen(filename,"r");
	while( fgets(line, sizeof(line),fid) != NULL ) {
    	if(sscanf(line,"%s\t%f\t%f\t%f\t%f",parameter,&value,&value2,&value3,&value4) < 2){
		}else if(strcmp("acq_bw",parameter) == 0){ acq_bw = value;
		}else if(strcmp("xres",parameter) == 0){ xres = (int)value;
		}else if(strcmp("numrecv",parameter) == 0){ num_coils = (int)value;
		}else if(strcmp("slices",parameter) == 0){ num_slices = (int)value;
		}else if(strcmp("2d_flag",parameter) == 0){ ss_2d = (int)value;
		}else if(strcmp("nproj",parameter) == 0){ num_readouts = (int)value;
		}else if(strcmp("rcxres",parameter) == 0){ 
			rcxres     = (rcxres == -1 ) ?  ( (int)value ) : ( rcxres);
		}else if(strcmp("rcyres",parameter) == 0){ 
			rcyres     = (rcyres == -1 ) ?  ( (int)value ) : ( rcyres);
		}else if(strcmp("rczres",parameter) == 0){ 
			rczres     = (rczres == -1 ) ?  ( (int)value ) : ( rczres);
		}else if(strcmp("multi_echo",parameter) == 0){ multi_echo = (int)value;
		}else if(strcmp("num_encodes",parameter) == 0){ rcencodes = (int)value;
		}
	}
	fclose(fid);
}



Array< complex<float>,5 > RECON::reconstruction( int argc, char **argv, MRI_DATA& data){
	
	
	// Option to compress coils
	if (compress_coils > 0){
		data.coilcompress(compress_coils);
		num_coils = data.Num_Coils;
	}

	// Turn of parallel processing for 2D due to thread overhead
	if(rczres ==1){
		omp_set_num_threads(1);
	}else{
		if( omp_get_max_threads() > 8){
			omp_set_num_threads(omp_get_max_threads()-2);
		}
	}
		
	// Shorthand for Blitz++
	Range all=Range::all();

	// Matlab like timer (openmp code base)
	tictoc T; 

	// Setup Gridding + FFT Structure
	gridFFT gridding;
	gridding.read_commandline(argc,argv);
	gridding.precalc_gridding(rczres,rcyres,rcxres,3);

	// ------------------------------------
	//  Get coil sensitivity map ( move into function)
	// ------------------------------------

	Array<complex<float>,4>smaps; 
	if( (recon_type != SOS) && (data.Num_Coils>1)){ 
		cout << "Getting Coil Sensitivities " << endl<< flush; 

		// Low Pass filtering for Sensitivity Map
		gridding.k_rad = 16;

		// Allocate Storage for Map	and zero	
		cout << "Allocate Sense Maps"  << endl << flush;
		smaps.setStorage( ColumnMajorArray<4>());
		smaps.resize(rcxres,rcyres,rczres,data.Num_Coils);
		smaps=0;
		
				
		cout << "Recon Low Resolution Images"  << endl<< flush; 
		for(int e=0; e< rcencodes;e++){
			Array< float,3 >kxE = data.kx(all,all,all,e); 
			Array< float,3 >kyE = data.ky(all,all,all,e); 
			Array< float,3 >kzE = data.kz(all,all,all,e); 
			Array< float,3 >kwE = data.kw(all,all,all,e); 
	
			for(int coil=0; coil< data.Num_Coils; coil++){
				cout << "Coil = " << coil  << " encode = " << e << endl;
				// Arrays 
				Array<complex<float>,3>kdataE = data.kdata(all,all,all,e,coil); 
				
				//Do Gridding
				Array< complex<float>,3>smapC =smaps(all,all,all,coil);			
				gridding.forward(smapC,kdataE,kxE,kyE,kzE,kwE);
			}
		}
		gridding.k_rad = 999;
		

		// Spirit Code
		if (coil_combine_type==ESPIRIT){
 			cout << "eSPIRIT Based Maps"  << endl; 
			SPIRIT S;
      	 	S.read_commandline(argc,argv);
      		S.init(rcxres,rcyres,rczres,data.Num_Coils);
		    S.generateEigenCoils(smaps);		  
		}else{ // E-spirit Code 
		
			// Sos Normalization 
			cout << "Normalize Coils" << endl;
			#pragma omp parallel for 
			for(int k=0; k<smaps.length(2); k++){
				for(int j=0; j<smaps.length(1); j++){
					for(int i=0; i<smaps.length(0); i++){
						float sos=0.0;
						for(int coil=0; coil< data.Num_Coils; coil++){
							sos+= norm(smaps(i,j,k,coil));
						}
						sos = 1./sqrtf(sos);
						for(int coil=0; coil< data.Num_Coils; coil++){
							smaps(i,j,k,coil) *= sos;
						}
			}}}

		} // Normalization

		// Export 
		if(export_smaps==1){
			cout << "Exporting Smaps" << endl;
			ArrayWrite(smaps,"SenseMaps.dat");
			
		}
	}else{
		// Allocate Storage for Map	and zero	
		cout << "Allocate Sense Maps:" <<   data.Num_Coils << endl << flush;
		smaps.setStorage( ColumnMajorArray<4>());
		smaps.resize(rcxres,rcyres,rczres,data.Num_Coils);
		smaps=complex<float>(1.0,0.0);
	
	}	
	

	// ------------------------------------
	//  If time resolved need to sort the times in to bins (need to move to function calls)
	// ------------------------------------

	if(rcframes>1){
		float max_time =max(data.times);
		float min_time =min(data.times);

		cout << "Time Range :: " << min_time << " to " << max_time << endl;

		float delta_t = (max_time - min_time ) / ((float)(rcframes));

		// Sort times into discrete frames
		int *frame_count = new int[rcframes];
		memset( (void *)frame_count,0,(size_t)((rcframes)*sizeof(int)));
		for(int e=0; e< data.times.length(3); e++){
			for(int k=0; k< data.times.length(2); k++){
				for(int j=0; j< data.times.length(1); j++){
					for(int i=0; i< data.times.length(0); i++){

						data.times(i,j,k,e) -= min_time;
						data.times(i,j,k,e) /= delta_t;
						int pos = (int)data.times(i,j,k,e);
						if(pos > (rcframes-1)){
							pos= (rcframes-1);
						}
						data.times(i,j,k,e) = (float)pos;

						frame_count[pos]++;
					}}}}

		for(int t =0; t<rcframes; t++){
			cout << "Frame " << t << " count " << frame_count[t] << endl;
		}
		delete [] frame_count;
	}
	
	/* ----- Temp For complex diff ----*/
	if(complex_diff){
		cout << "Doing Complex Diff" << endl;
		for(int coil=0; coil <data.Num_Coils; coil++){ 
			Array<complex<float>,3>kdata1 = data.kdata(all,all,all,0,coil); 
			Array<complex<float>,3>kdata2 = data.kdata(all,all,all,1,coil); 
			kdata1 -= kdata2;
		}
		rcencodes = 1; 
	}
	

	/*----------------------------Main Recons---------------------------------------*/	

	
	// Final Image Solution
	Array< complex<float>,5 >X(rcxres,rcyres,rczres,rcframes,rcencodes,ColumnMajorArray<5>());
	X=0;

	// Weighting Array for Time coding
	Array< float, 3 >TimeWeight(data.kx.length(0),data.kx.length(1),data.kx.length(2),ColumnMajorArray<3>());

	switch(recon_type){
		default:
		case(SOS):
		case(PILS):{

					 for(int e=0; e< rcencodes; e++){
						 for(int t=0; t< rcframes; t++){
							 cout << "Recon Encode" << e << " Frame " << t << endl;

							 // Get Sub-Arrays for Encoding (Blitz-Subarray reference, no memory copied)
							 Array< float,3 >kxE = data.kx(all,all,all,e); 
							 Array< float,3 >kyE = data.ky(all,all,all,e); 
							 Array< float,3 >kzE = data.kz(all,all,all,e); 
							 Array< float,3 >kwE = data.kw(all,all,all,e); 
							 Array< float,3 >timesE = data.times(all,all,all,e); 

							 // Temporal weighting (move to functions )
							 TimeWeight = kwE;
							 if(rcframes>1){
								 for(int k=0;k<timesE.extent(2);k++){
									 for(int j=0;j<timesE.extent(1);j++){
										 for(int i=0;i<timesE.extent(0);i++){
											 if( (timesE(i,j,k)) != (float)t){
												 TimeWeight(i,j,k)= 0.0;
											 }
										 }}}
							 }

							 cout << "\tForward Gridding Coil ";
							 for(int coil=0; coil< data.Num_Coils; coil++){
								 // Subarray for Data
								 Array<complex<float>,3>kdataE = data.kdata(all,all,all,e,coil); 

								 cout << coil << "," << flush;

								 // Image to add too
								 Array<complex<float>,3>xet = X(all,all,all,t,e);
								 
								 T.tic();
								 if(recon_type==PILS){
								 	 // adds to xet, does gridding + multiply by conj(smap)
									 Array< complex<float>,3>smapC =smaps(all,all,all,coil);			
									 gridding.forward(xet,smapC,kdataE,kxE,kyE,kzE,TimeWeight);
								 }else{
								 	 gridding.image = 0;
									 gridding.forward(gridding.image,kdataE,kxE,kyE,kzE,TimeWeight);
									 gridding.image *=conj(gridding.image);
									 xet += gridding.image;
								 }								 
								 cout << "Gridding took = " << T << endl;
							 }
							 cout << endl;
						 }
					 }

					 // Take Square Root for SOS	
					 if(recon_type==SOS){
						 X = csqrt(X);
					 }

				 }break;

		case(CG):{
				       // ------------------------------------
				       // Conjugate Gradient Recon not yet
				       // ------------------------------------



			       }


		case(IST):
		case(FISTA):{

					  // ------------------------------------
					  // Iterative Soft Thresholding  x(n+1)=  thresh(   x(n) - E*(Ex(n) - d)  )
					  //  Designed to not use memory
					  // Uses gradient descent x(n+1) = x(n) - ( R'R ) / ( R'E'E R) * Grad  [ R = E'(Ex-d)]
					  // ------------------------------------
					  
					  float reg_scale = 0.0;	
						
					  // Previous Array for FISTA
					  Array< complex<float>,5>X_old;
					  if( recon_type == FISTA){
						  X_old.setStorage(ColumnMajorArray<5>());
						  X_old.resize( X.shape());				  
						  X_old = 0.0;
					  }

					  // Residue 	
					  Array< complex<float>,5>R(X.shape(),ColumnMajorArray<5>());
					  R=0.0;

					  // Temp variable for E'ER 
					  Array< complex<float>,3 >P(rcxres,rcyres,rczres,ColumnMajorArray<3>());

					  // Storage for (Ex-d)
					  Array< complex<float>,3 >diff_data(data.kdata.length(0),data.kdata.length(1),data.kdata.length(2),ColumnMajorArray<3>());

					  // Setup 3D Wavelet
					  int dirs[3] = {4, 4, 4};
					  Array< complex<float>,3>Xref=X(all,all,all,0,0);
					  WAVELET3D wave(Xref,dirs,WAVE_DB4);

					  // Temporal differences or FFT
					  TDIFF tdiff(X);

					  // Setup Soft Thresholding
					  THRESHOLD softthresh(argc,argv);

					  cout << "Iterate" << endl;
					  double error0=0.0;
					  for(int iteration =0; iteration< max_iter; iteration++){

						  tictoc iteration_timer;
						  iteration_timer.tic();
						  cout << "\nIteration = " << iteration << endl;

						  // Update X based on FISTA 
						  if(recon_type==FISTA){
							  softthresh.fista_update(X,X_old,iteration);
						  }

						  // Zero this for Cauchy set size
						  complex<float>scale_RhP(0,0);						  

						  // Get Residue
						  R=0; 								
						  cout << "\tGradient Calculation" << endl;
						  for(int e=0; e< rcencodes; e++){
							  for(int t=0; t< rcframes; t++){
								  T.tic();

								  // Get Sub-Arrays for Encoding
								  Array< float,3 >kxE = data.kx(all,all,all,e); 
								  Array< float,3 >kyE = data.ky(all,all,all,e); 
								  Array< float,3 >kzE = data.kz(all,all,all,e); 
								  Array< float,3 >kwE = data.kw(all,all,all,e); 
								  Array< float,3 >timesE = data.times(all,all,all,e); 

								  // Temporal weighting
								  TimeWeight = kwE;
								  if(rcframes>1){
									  for(int k=0;k<timesE.extent(2);k++){
										  for(int j=0;j<timesE.extent(1);j++){
											  for(int i=0;i<timesE.extent(0);i++){
												  if( (timesE(i,j,k)) != (float)t){
													  TimeWeight(i,j,k)= 0.0;
												  }
											  }}}
								  }


								  Array<complex<float>,3>Xref=X(all,all,all,t,e);
								  Array<complex<float>,3>Rref=R(all,all,all,t,e);
								  for(int coil=0; coil< data.Num_Coils; coil++){
									  // Ex
									  diff_data=0;
									  Array<complex<float>,3>smapC=smaps(all,all,all,coil);
									  gridding.backward(Xref,smapC,diff_data,kxE,kyE,kzE,TimeWeight);

									  //Ex-d
									  Array< complex<float>,3>kdataC = data.kdata(all,all,all,e,coil); 
									  diff_data -= kdataC;

									  //E'(Ex-d)
									  gridding.forward( Rref,smapC,diff_data,kxE,kyE,kzE,TimeWeight);
								  }//Coils
								 
								  // TV of Image
								  int Nx = Xref.length(firstDim);
								  int Ny = Xref.length(secondDim);
								  int Nz = Xref.length(thirdDim);
								  for(int k =0; k < Xref.length(thirdDim);k++){
								  for(int j =0; j < Xref.length(secondDim);j++){
								  for(int i =0; i < Xref.length(firstDim);i++){
								  	Rref(i,j,k) += reg_scale*( complex<float>(6.0,0.0)*Xref(i,j,k) - Xref((i+1)%Nx,j,k) - Xref((i+Nx-1)%Nx,j,k) - Xref(i,(j+1)%Ny,k) - Xref(i,(j+Ny-1)%Ny,k)  - Xref(i,j,(k+1)%Nz) - Xref(i,j,(k+Nz-1)%Nz));
								  }}}
								    
								  //Now Get Scale factor (for Cauchy-Step Size)
								  P=0;
								  for(int coil=0; coil< data.Num_Coils; coil++){
									  // EE'(Ex-d)
									  diff_data=0;
									  Array<complex<float>,3>smapC=smaps(all,all,all,coil);
									  gridding.backward(Rref,smapC, diff_data,kxE,kyE,kzE,TimeWeight);

									  //E'EE'(Ex-d)
									  gridding.forward(P,smapC,diff_data,kxE,kyE,kzE,TimeWeight);
								  }//Coils
								  
								  // TV of Image
								  for(int k =0; k < Xref.length(thirdDim);k++){
								  for(int j =0; j < Xref.length(secondDim);j++){
								  for(int i =0; i < Xref.length(firstDim);i++){
								  	P(i,j,k) += reg_scale*( complex<float>(6.0,0.0)*Rref(i,j,k) - Rref((i+1)%Nx,j,k) - Rref((i+Nx-1)%Nx,j,k) - Rref(i,(j+1)%Ny,k) - Rref(i,(j+Ny-1)%Ny,k)  - Rref(i,j,(k+1)%Nz) - Rref(i,j,(k+Nz-1)%Nz));
								  }}}
								  
								  P*=conj(Rref);
								  
								  scale_RhP += sum(P); 
								  cout << "took " << T << "s" << endl;
							  }//Time
						  }//Encode

						  // Get Scaling Factor R'P / R'R 
						  complex<float>scale_RhR = complex<float>(ArrayEnergy(R),0);

						  // Error check
						  if(iteration==0){
							  error0 = abs(scale_RhR);
						  }
						  cout << "Residue Energy =" << scale_RhR << " ( " << (abs(scale_RhR)/error0) << " % )  " << endl;

						  // Export R		
						  Array<complex<float>,2>Rslice=R(all,all,R.length(2)/2,0,0);
						  ArrayWriteMag(Rslice,"R.dat");						  

						  // Step in direction
						  complex<float>scale = (scale_RhR/scale_RhP);
						  cout << "Scale = " << scale << endl;
						  R *= scale;
						  X -= R;
						  cout << "Took " << iteration_timer << " s " << endl;
						  
						  //Temp L2 of image
						  if(iteration == 1){
						  	reg_scale = 0.005*sqrt(abs(scale_RhR))/(float)R.size(); // Add sqrt of scale Ex-d
						  	cout << "Reg Scale = " << reg_scale << endl;
							cout << "Energy X = " << ArrayEnergy(X) << endl;
						  }
						  
						  // Export X slice
						  Array<complex<float>,2>Xslice=X(all,all,X.length(2)/2,0,0);
						  ArrayWriteMag(Xslice,"X_mag.dat");
						  ArrayWritePhase(Xslice,"X_phase.dat");

						  // ------------------------------------
						  // Soft thresholding operation (need to add transform control)
						  // ------------------------------------
						  if(iteration==0){
						  	for(int k=0; k<X.extent(thirdDim);k++){
							for(int j=0; j<X.extent(thirdDim);j++){
							for(int i=0; i<X.extent(thirdDim);i++){
							
							complex<float>s(0,0);							
							for(int t=0; t<X.extent(fourthDim);t++){
								s += X(i,j,k,t,0);
							}
							s /= (float)X.extent(fourthDim);
							for(int t=0; t<X.extent(fourthDim);t++){
								X(i,j,k,t,0)=s;
							}
														
							
							}}}
							
						  }else if(softthresh.getThresholdMethod() != TH_NONE){
							  cout << "Soft thresh" << endl;
							  tdiff.fft_t(X);
							  wave.random_shift();
							  wave.forward(X);	
							  softthresh.exec_threshold(X);
							  wave.backward(X);
							  tdiff.ifft_t(X);
						  }
  						  ArrayWriteMag(Xslice,"X_mag_post.dat");
	
					  }// Iteration			

				  }break;

	}//Recon Type


	cout << "Recon was completed successfully " << endl;	
	return(X);
}





