#include "recon_lib.h"
#include "io_templates.hpp"

using arma::cx_mat;
using arma::vec;
using arma::uvec;
using namespace NDarray;

// ----------------------
//  Basic constructor (no args)
// ----------------------
RECON::RECON(void){
	set_defaults();
}

// ----------------------
//  Sets Default Recon Parameterrs
// ----------------------
void RECON::set_defaults( void){
	// Help Message for recon
	recon_type = SOS; 
	data_type = EXTERNAL;
	coil_combine_type = LOWRES;
	
	complex_diff = false;
	
	cs_spatial_transform = WAVELET;
	cs_temporal_transform = NONE;
	cs_encode_transform = NONE;
	
	zoom = 1.0;
	zoom_x = 1.0;
	zoom_y = 1.0;
	zoom_z = 1.0;
	  
	rcxres=-1;
	rcyres=-1;
	rczres=-1;
	rcframes=1;
	rcencodes=1;
	
	smap_res=8;
	
	acc = 1;
	compress_coils = 0.0;
	whiten = false;
	export_smaps = 0;
	max_iter = 50;

	prep_done = false;
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
			GATING::help_message();
			L2REG::help_message();
			LOWRANKCOIL::help_message();
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
	help_flag("-clear","clear (low rank coil aproximation)");
	help_flag("-ist","iterative soft thresholding");
	help_flag("-fista","fast iterative soft thresholding");
	help_flag("-cg","conjugate gradients");
	help_flag("-complex_diff","Subtract first encode");
	
	cout << "Transforms for Compressed Sensing:" << endl;
	help_flag("-spatial_transform []","none/wavelet");
    help_flag("-temporal_transform []","none/diff/pca/dft/wavelet");
    help_flag("-encode_transform []","none/diff");
    	
	cout << "Iterative Recon Control:" << endl;
	help_flag("-max_iter []","max iterations for iterative recons");
    
	cout << "Coil Control:" << endl;
	help_flag("-espirit","use ESPIRIT to get coil sensitivies");
	help_flag("-coil_lowres","default,use low resolution images to get coil sensitivies");
	help_flag("-export_smaps","write sensitivity maps");
}

/** --------------------
 *  Doxygen test: Read command line and set variables
 * --------------------*/

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
		trig_flag(CLEAR,"-clear",recon_type);
		trig_flag(CG,"-cg",recon_type);
		
		// Spatial Transforms
		}else if(strcmp("-spatial_transform",pstring[pos]) == 0) {
			pos++;
			if( pos==numarg){
				cout << "Please provide spatial transform type..none/wavelet" << endl;
				exit(1);
				trig_flag(NONE,"none",cs_spatial_transform);
				trig_flag(WAVELET,"wavelet",cs_spatial_transform);
			}else{
				cout << "Please provide spatial transform type..none/wavelet" << endl;
				exit(1);
			}
		
		// Temporal Transforms
		}else if(strcmp("-temporal_transform",pstring[pos]) == 0) {
			pos++;
			if( pos==numarg){
				cout << "Please provide temporal transform type..none/dft/diff/pca" << endl;
				exit(1);
				trig_flag(NONE,"none",cs_temporal_transform);
				trig_flag(DIFF,"diff",cs_temporal_transform);
				trig_flag(DFT,"dft",cs_temporal_transform);
				trig_flag(PCA,"pca",cs_temporal_transform);
				trig_flag(WAVELET,"wavelet",cs_temporal_transform);
				trig_flag(COMPOSITE_DIFF,"composite",cs_temporal_transform);
				
			}else{
				cout << "Please provide temporal transform type..none/dft/diff/pca" << endl;
				exit(1);
			}
				
		// Encode Transforms
		}else if(strcmp("-encode_transform",pstring[pos]) == 0) {
			pos++;
			if( pos==numarg){
				cout << "Please provide encode transform type..none/dft/diff/pca" << endl;
				exit(1);
				trig_flag(NONE,"none",cs_encode_transform);
				trig_flag(DIFF,"diff",cs_encode_transform);
			}else{
				cout << "Please provide encode transform type..none/diff" << endl;
				exit(1);
			}
				
		// Coil Combination		
		trig_flag(ESPIRIT,"-espirit",coil_combine_type);
		trig_flag(LOWRES,"-coil_lowres",coil_combine_type);
		float_flag("-smap_res",smap_res);
		trig_flag(1,"-export_smaps",export_smaps);
		
		// Source of data
		trig_flag(EXTERNAL,"-external_data",data_type);
		trig_flag(PFILE,"-pfile",data_type);
		trig_flag(PHANTOM,"-phantom",data_type);
		trig_flag(SIMULATE,"-simulate",data_type);
		trig_flag(PSF,"-psf",data_type);
				
		// Data modification
		int_flag("-acc",acc);
		float_flag("-compress_coils",compress_coils);
		trig_flag(true,"-whiten",whiten);
		trig_flag(true,"-complex_diff",complex_diff);
		
		// Iterations for IST
		int_flag("-max_iter",max_iter);
		
	}
  }
} 

complex<float> conj_sum( Array<complex<float>,3>P,Array<complex<float>,3>R){
		complex<float>s(0,0);
		
		complex<float> *stemp = new complex<float>[R.length(thirdDim)];
				
		#pragma omp parallel for					
		for(int k =0; k< R.length(thirdDim); k++){
		stemp[k] = complex<float>(0.0,0.0);
		for(int j =0; j< R.length(secondDim); j++){
		for(int i =0; i< R.length(firstDim); i++){
			stemp[k]+= P(i,j,k)*conj( R(i,j,k));
		}}}
		
		for(int k =0; k< R.length(thirdDim); k++){
			s+= stemp[k];
		}
			
		delete [] stemp;					
		return(s);
}

void RECON::init_recon(int argc, char **argv, MRI_DATA& data ){
	
	// Use Native Resultion 
	rcxres = (rcxres == -1) ? ( data.xres ) : ( rcxres );
	rcyres = (rcyres == -1) ? ( data.yres ) : ( rcyres );
	rczres = (rczres == -1) ? ( data.zres ) : ( rczres );
	rcencodes = data.Num_Encodings;
		
	// Whiten
	if( whiten ){
		data.whiten(); // Requires noise samples inserted
	}
	
	// Option to compress coils
	if (compress_coils > 0){
		data.coilcompress(compress_coils);
	}
		
	// Shorthand for Blitz++
	Range all=Range::all();

	// Matlab like timer (openmp code base)
	tictoc T; 

	// Setup Gridding + FFT Structure
	gridding.read_commandline(argc,argv);
	gridding.precalc_gridding(rczres,rcyres,rcxres,3);
	
	// Calculate Sensitivity maps (using gridding struct)	
	calc_sensitivity_maps(  argc,argv, data);
	
	
	// Complex Difference
	if(complex_diff){
		cout << "Doing Complex Diff" << endl;
		for(int coil=0; coil <data.Num_Coils; coil++){ 
			
			// Subtract off reference
			for(int e=1; e< rcencodes; e++){
				Array<complex<float>,3>kdata1 = data.kdata(0,coil); 
				Array<complex<float>,3>kdata2 = data.kdata(e,coil); 
				kdata2 -= kdata1;
			}
			
			// Rearrange Positions
			for(int e=1; e< rcencodes; e++){
				Array<complex<float>,3>kdata1 = data.kdata(e-1,coil); 
				Array<complex<float>,3>kdata2 = data.kdata(e,coil); 
				kdata1 = kdata2;
			}
		}
		rcencodes -= 1; 
		
	}
		
	// -------------------------------------
	//	This handles all the gating, assuming mri_data physio data is populated 
	// -------------------------------------
	gate=GATING(argc,argv);
	gate.init( data,rcframes);
		
	
	// -------------------------------------
	//	This handles special preperations
	// -------------------------------------
	switch(recon_type){
		
		// Non-Iterative Recons 
		default:
		case(SOS):
		case(PILS):{
			
	
		}break;
				
		
		// Iterative Recons allow Regularization		
		case(CG):
		case(IST):
		case(FISTA):
		case(CLEAR):{
			
			// Setup 3D Wavelet
			wave = WAVELET3D( TinyVector<int,3>(rcxres,rcyres,rczres),TinyVector<int,3>(4,4,4),WAVELET3D::WAVE_DB4);

			// Temporal differences or FFT
			tdiff=TDIFF( rcframes,rcencodes );

			// Setup Soft Thresholding
			softthresh = THRESHOLD(argc,argv);
					
			// CLEAR  
			if(recon_type==CLEAR){
				lrankcoil=LOWRANKCOIL(argc,argv);
			}
			
			// Setup L2 Regularization Thresholding
			l2reg=L2REG(argc,argv);
			
			if(cs_temporal_transform == COMPOSITE_DIFF){
				composite_image.setStorage( ColumnMajorArray<3>());
				composite_image.resize( rcxres,rcyres,rczres);
				composite_image = complex<float>(0.0,0.0);
			}
		}break;
	}
		
	// Signal that the recon is ready to perform the requested reconstruction
	prep_done = true;
}

Array< Array<complex<float>,3 >,1 >RECON::reconstruct_one_frame( MRI_DATA& data, int frame_number){
	Array< Array<complex<float>,3 >,2 >XX = full_recon( data, Range(frame_number,frame_number),Range(0,0),false);
	Array< Array<complex<float>,3 >,1 >X = XX(0,Range::all());
	return(X);
}

Array< Array<complex<float>,3 >,2 >RECON::reconstruct_all_frames( MRI_DATA& data){
	return(full_recon( data, Range(0,rcframes-1),Range(0,rcframes-1),false));
}

Array< Array<complex<float>,3 >,1 >RECON::reconstruct_composite( MRI_DATA& data){
	Array< Array<complex<float>,3 >,2 >XX = full_recon( data, Range(0,0),Range(0,0),true);
	Array< Array<complex<float>,3 >,1 >X = XX(0,Range::all());
	return(X);
}


Array< Array<complex<float>,3 >,2 >RECON::full_recon( MRI_DATA& data, Range times, Range times_store, bool composite){
	
	// Shorthand for Blitz++
	Range all=Range::all();

	// Matlab like timer (openmp code base)
	tictoc T; 
	
	// Calculate shared structures
	if(!prep_done){
		cout << "Error need to run prep before running recon (obj.init_recon(arg, argc, data). Exiting." << endl;
		exit(1);
	}
	
	/*----------------------------Main Recons---------------------------------------*/	
	
	int Nt = times.length();
	GATING::FrameType frame_type;
	if( composite){
		frame_type=GATING::COMPOSITE;
	}else{
		frame_type=GATING::TIME_FRAME;
	}
	
	// Final Image Solution
	Array< Array<complex<float>,3>,2>X = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
	if( cs_temporal_transform==COMPOSITE_DIFF){
		for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
			(*miter)=composite_image;
		}
	}
	
	// Weighting Array for Time coding
	Array< float, 3 >TimeWeight(data.Num_Pts,data.Num_Readouts,data.Num_Slices,ColumnMajorArray<3>());

	switch(recon_type){
		default:
		case(SOS):
		case(PILS):{

					 for(int e=0; e< rcencodes; e++){
						 for(int t=0; t< Nt; t++){
							 
							 int act_t   = times(t);
							 int store_t = times_store(t);
							 
							 cout << "Recon Encode" << e << " Frame " << t << endl;

							 // Temporal weighting (move to functions )
							 TimeWeight = data.kw(e);
							 gate.weight_data( TimeWeight, e, data.kx(e),data.ky(e),data.kz(e),act_t,GATING::NON_ITERATIVE,frame_type);
   							 
							 cout << "\tForward Gridding Coil ";
							 for(int coil=0; coil< data.Num_Coils; coil++){
								 
								 cout << coil << "," << flush;
								 T.tic();
								 if(recon_type==PILS){
								 	 // adds to xet, does gridding + multiply by conj(smap)
									 gridding.forward(X(store_t,e),smaps(coil),data.kdata(e,coil),data.kx(e),data.ky(e),data.kz(e),TimeWeight);
								 }else{
								 	 gridding.image = 0;
									 gridding.forward(gridding.image,data.kdata(e,coil),data.kx(e),data.ky(e),data.kz(e),TimeWeight);
									 gridding.image *=conj(gridding.image);
									 X(store_t,e) += gridding.image;
								 }								 
								 cout << "Gridding took = " << T << endl;
							 }
						
							 // Take Square Root for SOS	
					 		 if(recon_type==SOS){
							 	X(store_t,e) = csqrt(X(store_t,e));
							 }
							 cout << endl;
						 }
					 }

				 }break;

		case(CLEAR):{
					
					cout << "----------------------------------------" << endl;
					cout << "\tStarting Clear based recon" << endl;
					cout << "----------------------------------------" << endl;
					typedef Array<complex<float>,3> Complex3D;
					
					
					// Image and Residue (note now need to store coils as well)
					Array< Complex3D, 3 >XX = Alloc6DContainer<complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes,data.Num_Coils);
					Array< Complex3D, 3 >RR = Alloc6DContainer<complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes,data.Num_Coils);
				
					// Temp variable for E'ER 
					Complex3D P(rcxres,rcyres,rczres,ColumnMajorArray<3>());

					// Storage for (Ex-d)
					Complex3D diff_data(data.Num_Pts,data.Num_Readouts,data.Num_Slices,ColumnMajorArray<3>());
					
					cout << "Iterate" << endl;
					double error0=0.0;
					for(int iteration =0; iteration< max_iter; iteration++){

						  tictoc iteration_timer;
						  iteration_timer.tic();
						  cout << "\nIteration = " << iteration << endl;

						  //----------------------------------------------------
						  //  First get Gradient Descent 
						  //----------------------------------------------------
						  						  
						  // Zero this for Cauchy set size
						  complex<float>scale_RhP(0,0);						  
						  
						  cout << "\tGradient Calculation" << endl;
						  for(int e=0; e< rcencodes; e++){
							  for(int t=0; t< Nt; t++){
							 	  int act_t   = times(t);
							 	  int store_t = times_store(t);
								  
								  T.tic();

								  // Get Sub-Arrays for Encoding
								  Array< float,3 >kxE = data.kx(e); 
								  Array< float,3 >kyE = data.ky(e); 
								  Array< float,3 >kzE = data.kz(e); 
								  Array< float,3 >kwE = data.kw(e); 
								  
								  // Temporal weighting
								  TimeWeight = kwE;
								  gate.weight_data( TimeWeight,e, kxE, kyE,kzE,act_t,GATING::ITERATIVE, frame_type);
   							 	  
								  for(int coil=0; coil< data.Num_Coils; coil++){
									  T.tic();
									  
									  // Ex
									  diff_data=0;
									  gridding.backward(XX(store_t,e,coil),diff_data,kxE,kyE,kzE,TimeWeight);
									  									  
									  //Ex-d
									  Array< complex<float>,3>kdataC = data.kdata(e,coil); 
									  diff_data -= kdataC;
									  
									  //E'(Ex-d)
									  RR(t,e,coil)=complex<float>(0.0,0.0);
									  gridding.forward( RR(store_t,e,coil),diff_data,kxE,kyE,kzE,TimeWeight);
									  
									  // Now Get Scale
									  P=0;
									   							  
									  // EE'(Ex-d)
									  diff_data=0;
									  gridding.backward(RR(store_t,e,coil), diff_data,kxE,kyE,kzE,TimeWeight);
									  
									  //E'EE'(Ex-d)
									  gridding.forward(P,diff_data,kxE,kyE,kzE,TimeWeight);
								  	  
									  scale_RhP += conj_sum(P,RR(t,e,coil)); 
									  
									  cout << "Coil " << coil << " took " << T << endl;
								  }//Coils
								  
								  // cout << e << "," << t << "took " << T << "s" << endl;
							  }//Time
						  }//Encode

						  // Get Scaling Factor R'P / R'R 
						  complex<float>scale_RhR = complex<float>(ArrayEnergy(RR),0);

						  // Error check
						  if(iteration==0){
							  error0 = abs(scale_RhR);
						  }
						  cout << "Residue Energy =" << scale_RhR << " ( " << (abs(scale_RhR)/error0) << " % )  " << endl;

						  // Export R (across coils)	
						  Array<complex<float>,2>Rslice=RR(0,0,0)(all,all,RR(0,0,0).length(2)/2);
						  ArrayWriteMag(Rslice,"R.dat");						  

						  // Step in direction
						  complex<float>scale = (scale_RhR/scale_RhP);
						  cout << "Scale = " << scale << endl;
						  
						  for(int e=0; e< rcencodes; e++){
							  for(int t=0; t< Nt; t++){
						  		RR(t,e) *= scale;
						  		XX(t,e) -= RR(t,e);
						  }}
						  cout << "Took " << iteration_timer << " s " << endl;
												  
						  // Export X slice
						  Array<complex<float>,2>Xslice=XX(0,0,0)(all,all,RR(0,0,0).length(2)/2);
						  ArrayWriteMag(Xslice,"X_mag.dat");
						  ArrayWritePhase(Xslice,"X_phase.dat");
						  
						 					  
						  
						  // ---------------------------------
						  //   Clear 
						  // ---------------------------------
						 
						  iteration_timer.tic();
						  lrankcoil.update_threshold( XX);
						  cout << "Get thresh took " << iteration_timer << " s" << endl;
						  
						  iteration_timer.tic();
						  lrankcoil.thresh(XX ,RR );
						  cout << "Thresh took " << iteration_timer << " s" << endl;
						  
						  // -- Note missing thresholding ---
						  // ToDo
							
					  }// Iteration	
				      
					  lrankcoil.combine(XX,X);
					  
		
		}break;
		
		case(CG):{
				 // ------------------------------------
				 // Conjugate Gradient Recon 
				 //   -Uses much more memory than gradient descent but is faster (both convergence + per iteration)
				 // ------------------------------------
				
				 // Structures  	
				 Array< Array< complex<float>,3>, 2>R = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
				 for( Array< Array<complex<float>,3>,2>::iterator miter =R.begin(); miter != R.end(); miter++){
				 	(*miter)= complex<float>(0.0,0.0);
				 }
				 
				 Array< Array< complex<float>,3>, 2>P = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
				 for( Array< Array<complex<float>,3>,2>::iterator miter =P.begin(); miter != P.end(); miter++){
				 	(*miter)= complex<float>(0.0,0.0);
				 }
				 
				 Array< Array< complex<float>,3>, 2>LHS = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
				 				 
				 // Storage for (Ex-d)
				 Array< complex<float>,3 >diff_data(data.Num_Pts,data.Num_Readouts,data.Num_Slices,ColumnMajorArray<3>());

				 			 
				 // First Calculate E'd
				 cout << "LHS Calculation" << endl;
				 for(int e=0; e< rcencodes; e++){
					for(int t=0; t< Nt; t++){
						int act_t = times(t);
						int store_t = times_store(t);
						
						// Temporal weighting
						TimeWeight = -data.kw(e); // Note negative --> Residue
						gate.weight_data( TimeWeight,e,data.kx(e),data.ky(e),data.kz(e),act_t,GATING::ITERATIVE, frame_type);
   						
						// E'd
						for(int coil=0; coil< data.Num_Coils; coil++){
							gridding.forward( R(store_t,e),smaps(coil),data.kdata(e,coil),data.kx(e),data.ky(e),data.kz(e),TimeWeight);
						} 
						 
				 }}
				 
				 // Initiialize P
				 P  = R;	
														 
				 // Now Iterate
				 cout << "Iterate" << endl;
				 double error0=0.0;
				 for(int iteration =0; iteration< max_iter; iteration++){

					  tictoc iteration_timer;
					  iteration_timer.tic();
					  
					  for(int e=0; e< rcencodes; e++){
						  for(int t=0; t< Nt; t++){
							  int act_t = times(t);
							  int store_t = times_store(t);
							  
							  LHS(t,e ) = complex<float>(0.0,0.0);
							  T.tic();
 
							  // Temporal weighting
							  TimeWeight = data.kw(e);
							  gate.weight_data( TimeWeight,e,data.kx(e),data.ky(e),data.kz(e),act_t,GATING::ITERATIVE, frame_type);
   							  	  
							  for(int coil=0; coil< data.Num_Coils; coil++){
								  // E'Ex
								  diff_data=0;
								  gridding.backward( P(store_t,e),smaps(coil),diff_data,data.kx(e),data.ky(e),data.kz(e),TimeWeight);
								  gridding.forward( LHS(store_t,e),smaps(coil),diff_data,data.kx(e),data.ky(e),data.kz(e),TimeWeight);
							  }//Coils
						}// t
					}// e
						  
					//----------------------------------------------------
					//  Now perform gradient update
					// ---------------------------------------------------
					
					complex< float> sum_R0_R0(0.0,0.0);
					complex< float> sum_R_R(0.0,0.0);
					complex< float> sum_P_LHS(0.0,0.0);
					
					// Calc R'R and P'*LHS
					for(int e=0; e< rcencodes; e++){
					 for(int t=0; t< Nt; t++){
						for(int k=0; k < rczres; k++){
						for(int j=0; j < rcyres; j++){
						for(int i=0; i < rcxres; i++){
							sum_R0_R0 += norm( R(t,e)(i,j,k) );
							sum_P_LHS += conj( P(t,e)(i,j,k) )*LHS(t,e)(i,j,k);
						}}}
					}}	
					complex< float> scale = sum_R0_R0 / sum_P_LHS; 
					
										
					// Take step size
					for(int e=0; e< rcencodes; e++){
					 for(int t=0; t< Nt; t++){
						for(int k=0; k < rczres; k++){
						for(int j=0; j < rcyres; j++){
						for(int i=0; i < rcxres; i++){
							X(t,e)(i,j,k) += ( scale* (  P(t,e)(i,j,k)) );
							R(t,e)(i,j,k) -= ( scale* (LHS(t,e)(i,j,k)) );
							sum_R_R += norm( R(t,e)(i,j,k) );
						}}}
					}}
										
					cout << "Sum R'R = " << sum_R_R << endl;
					complex< float> scale2 = sum_R_R / sum_R0_R0; 
					
					// Take step size
					for(int e=0; e< rcencodes; e++){
					 for(int t=0; t< Nt; t++){
						for(int k=0; k < rczres; k++){
						for(int j=0; j < rcyres; j++){
						for(int i=0; i < rcxres; i++){
							P(t,e)(i,j,k) = R(t,e)(i,j,k) + ( scale2*P(t,e)(i,j,k) );
						}}}
					}}
					
					// Export X slice
					{
						  Array<complex<float>,2>Xslice=X(0,0)(all,all,X(0,0).length(2)/2);
						  ArrayWriteMagAppend(Xslice,"X_mag.dat");
						  
						  Array<complex<float>,2>LHSslice=LHS(0,0)(all,all,LHS(0,0).length(2)/2);
						  ArrayWriteMagAppend(LHSslice,"LHS_mag.dat");
						  
						  Array<complex<float>,2>Pslice=X(0,0)(all,all,P(0,0).length(2)/2);
						  ArrayWriteMagAppend(Pslice,"P_mag.dat");
						  
						  Array<complex<float>,2>Rslice=R(0,0)(all,all,R(0,0).length(2)/2);
						  ArrayWriteMagAppend(Rslice,"R_mag.dat");
					}
					
				}//iteration	
					
										



			   }break;


		case(IST):
		case(FISTA):{
	
					 
					  // ------------------------------------
					  // Iterative Soft Thresholding  x(n+1)=  thresh(   x(n) - E*(Ex(n) - d)  )
					  //  Designed to not use memory
					  // Uses gradient descent x(n+1) = x(n) - ( R'R ) / ( R'E'E R) * Grad  [ R = E'(Ex-d)]
					  // ------------------------------------
						
					  // Previous Array for FISTA
					  
					  Array< Array< complex<float>,3>, 2>X_old;
					  if( recon_type == FISTA){
						  	cout << "Alloc Fista Matrix" << endl;
							Array< Array< complex<float>,3>, 2>Temp =Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
							X_old.reference(Temp);
					  }

					  // Residue 	
					  Array< Array< complex<float>,3>, 2>R = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
	
					  // Temp variable for E'ER 
					  Array< complex<float>,3 >P(rcxres,rcyres,rczres,ColumnMajorArray<3>());

					  // Storage for (Ex-d)
					  Array< complex<float>,3 >diff_data(data.Num_Pts,data.Num_Readouts,data.Num_Slices,ColumnMajorArray<3>());


					  cout << "Iterate" << endl;
					  double error0=0.0;
					  for(int iteration =0; iteration< max_iter; iteration++){

						  tictoc iteration_timer;
						  iteration_timer.tic();
						  cout << "\nIteration = " << iteration << endl;

						  // Update X based on FISTA 
						  if(recon_type==FISTA){
							  cout << "Fista update" << endl << flush;
							  softthresh.fista_update(X, X_old, iteration);
						  }

						  // Zero this for Cauchy set size
						  complex<float>scale_RhP(0,0);						  

						  // Get Residue
						  for( Array< Array<complex<float>,3>,2>::iterator riter =R.begin(); riter != R.end(); riter++){
						  		*riter=0;
						  }
						  
						  cout << "\tGradient Calculation" << endl;
						  for(int e=0; e< rcencodes; e++){
							  for(int t=0; t< Nt; t++){
							  	int act_t = times(t);
								int store_t = times_store(t);
							  
								  T.tic();

								  // Get Sub-Arrays for Encoding
								  Array< float,3 >kxE = data.kx(e); 
								  Array< float,3 >kyE = data.ky(e); 
								  Array< float,3 >kzE = data.kz(e); 
								  Array< float,3 >kwE = data.kw(e); 
								  
								  // Temporal weighting
								  TimeWeight = kwE;
								  gate.weight_data( TimeWeight,e, kxE, kyE,kzE,act_t,GATING::ITERATIVE, frame_type);
   							 	  TimeWeight /= sum(TimeWeight);
								  
								  for(int coil=0; coil< data.Num_Coils; coil++){
									  // Ex
									  diff_data=0;
									  gridding.backward(X(store_t,e),smaps(coil),diff_data,kxE,kyE,kzE,TimeWeight);

									  //Ex-d
									  Array< complex<float>,3>kdataC = data.kdata(e,coil); 
									  diff_data -= kdataC;

									  //E'(Ex-d)
									  gridding.forward( R(store_t,e),smaps(coil),diff_data,kxE,kyE,kzE,TimeWeight);
								  }//Coils
								 
								  // L2 
								  if(iteration > 0){
								  	l2reg.regularize(R(store_t,e),X(store_t,e) );
								  }
								    
								  //Now Get Scale factor (for Cauchy-Step Size)
								  P=0;
								  for(int coil=0; coil< data.Num_Coils; coil++){
									  // EE'(Ex-d)
									  diff_data=0;
									  gridding.backward(R(store_t,e),smaps(coil), diff_data,kxE,kyE,kzE,TimeWeight);

									  //E'EE'(Ex-d)
									  gridding.forward(P,smaps(coil),diff_data,kxE,kyE,kzE,TimeWeight);
								  }//Coils
								  
								  // TV of Image
								  if(iteration > 0){
								  	l2reg.regularize(P,R(store_t,e));
								  }
								 
								  P*=conj(R(store_t,e));
								  
								  scale_RhP += sum(P); 
								  cout << e << "," << t << "took " << T << "s" << endl;
							  }//Time
						  }//Encode
						  
						  // Get Scaling Factor R'P / R'R 
						  cout << "Calc residue" << endl << flush;
						  complex<float>scale_RhR = 0.0;
						  for( Array< Array<complex<float>,3>,2>::iterator riter =R.begin(); riter != R.end(); riter++){
						  		scale_RhR += complex<float>( ArrayEnergy( *riter ), 0.0);
						  }
						  
						  // Error check
						  if(iteration==1){
							  cout << "L2 set scale " << endl << flush;
						  	  error0 = abs(scale_RhR);
							  l2reg.set_scale(error0,X);
						  }
						  cout << "Residue Energy =" << scale_RhR << " ( " << (abs(scale_RhR)/error0) << " % )  " << endl;

						  // Export R
						  		
						  Array<complex<float>,2>Rslice=R(0,0)(all,all,R(0,0).length(2)/2);
						  ArrayWriteMag(Rslice,"R.dat");						  

						  // Step in direction
						  complex<float>scale = (scale_RhR/scale_RhP);
						  cout << "Scale = " << scale << endl << flush;
						  for(int e=0; e< rcencodes; e++){
							  for(int t=0; t< Nt; t++){
						  		R(t,e) *= scale;
						  		X(t,e) -= R(t,e);
							}
						  }
						  cout << "Took " << iteration_timer << " s " << endl;
												  
						  // Export X slice
						  Array<complex<float>,2>Xslice=X(0,0)(all,all,X(0,0).length(2)/2);
						  ArrayWriteMagAppend(Xslice,"X_mag.dat");
						  
						  if( Nt > 1){
						  	ArrayWriteMag(Xslice,"X_frames.dat");
							ArrayWritePhase(Xslice,"X_frames.dat.phase");
							for(int t=1; t< Nt; t++){
								Array<complex<float>,2>Xf=X(t,0)(all,all,X(0,0).length(2)/2);
								ArrayWriteMagAppend(Xf,"X_frames.dat");
								ArrayWritePhaseAppend(Xf,"X_frames.dat.phase");
							}
						  }
						  
						  // ------------------------------------
						  // Soft thresholding operation (need to add transform control)
						  // ------------------------------------
						  if(softthresh.getThresholdMethod() != TH_NONE){
							  L1_threshold(X);
						  }
  						  ArrayWriteMagAppend(Xslice,"X_mag.dat");
						  
					  }// Iteration			

				  }break;

	}//Recon Type


	cout << "Recon was completed successfully " << endl;	
	return(X);
}


void RECON::L1_threshold( Array< Array< complex<float>,3>, 2>&X){
	
	cout << "Soft thresh" << endl;
	switch(cs_temporal_transform){
		case(DFT):{	
			cout << "DFT in Time" << endl;
			tdiff.fft_t(X); 
		}break;
								
		case(DIFF):{
			cout << "DIFF in Time" << endl;
			tdiff.tdiff(X); 
		}break;
								
		case(WAVELET):{
			cout << "WAVELET in Time" << endl;
			tdiff.twave(X); 
		}break;
		
		case(COMPOSITE_DIFF):{
				cout << "Composite Diff" << endl;
				for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
						*miter -= composite_image;
				}
		}break;
								
		default:{
		}break;
	  }
		
	  switch(cs_spatial_transform){
			case(WAVELET):{
				cout << "Wavelet in Space" << endl;
				wave.random_shift();
				for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
						wave.forward(*miter);
									}
				}break;
			
				default:{
								
				}break;
	 }
							  
						/*	  
							  if( Nt > 1){
						  		ArrayWriteMag(Xslice,"X_frames.dat.td");
								ArrayWritePhase(Xslice,"X_frames.dat.phase.td");
								for(int t=1; t< Nt; t++){
									Array<complex<float>,2>Xf=X(t,0)(all,all,X(0,0).length(2)/2);
									ArrayWriteMagAppend(Xf,"X_frames.dat.td");
									ArrayWritePhaseAppend(Xf,"X_frames.dat.phase.td");
						}*/
		 
		ArrayWriteMag(X(0,0),"PostTransform.dat");
							  
		softthresh.exec_threshold(X);
							  
		switch(cs_spatial_transform){
							  	case(WAVELET):{
							  		for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
										wave.backward(*miter);
									}
								}break;
								default:{
								
								}break;
		 }
							  
		switch(cs_temporal_transform){
							  	case(DFT):{ 
									tdiff.ifft_t(X); 
								}break;
								
								case(DIFF):{
									cout << "DIFF in Time" << endl;
									tdiff.inv_tdiff(X); 
								}break;
								
								case(WAVELET):{
									cout << "WAVE in Time" << endl;
									tdiff.inv_twave(X); 
								}break;
								
								case(COMPOSITE_DIFF):{
									cout << "Composite Diff" << endl;
									for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
										*miter += composite_image;
									}
								}break;																
								
								default:{
								
								}break;
		}


}



void RECON::calc_sensitivity_maps( int argc, char **argv, MRI_DATA& data){

	// ------------------------------------
	//  Get coil sensitivity map ( move into function)
	// ------------------------------------
	Range all = Range::all();	
	if( (recon_type != SOS) && (recon_type != CLEAR) && (data.Num_Coils>1) ){ 
		cout << "Getting Coil Sensitivities " << endl<< flush; 

		// Low Pass filtering for Sensitivity Map
		if(coil_combine_type!=ESPIRIT){
			gridding.k_rad = smap_res;
		}
		
		// Allocate Storage for Map	and zero	
		cout << "Allocate Sense Maps"  << endl << flush;
		smaps.setStorage( ColumnMajorArray<1>());
		smaps.resize( data.Num_Coils );
		for(int coil=0; coil< smaps.length(firstDim); coil++){
			smaps(coil).setStorage( ColumnMajorArray<3>());
			smaps(coil).resize(rcxres,rcyres,rczres,data.Num_Coils);
			smaps(coil)=0;
		}	
						
		cout << "Recon Low Resolution Images"  << endl<< flush; 
		for(int e=0; e< 1;e++){
			Array< float,3 >kxE = data.kx(e); 
			Array< float,3 >kyE = data.ky(e); 
			Array< float,3 >kzE = data.kz(e); 
			Array< float,3 >kwE = data.kw(e); 
	
			for(int coil=0; coil< data.Num_Coils; coil++){
				cout << "Coil = " << coil  << " encode = " << e << endl;
				// Arrays 
				Array<complex<float>,3>kdataE = data.kdata(e,coil); 
				
				//Do Gridding
				Array< complex<float>,3>smapC =smaps(coil);			
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
			for(int k=0; k<smaps(0).length(thirdDim); k++){
				for(int j=0; j<smaps(0).length(secondDim); j++){
					for(int i=0; i<smaps(0).length(firstDim); i++){
						float sos=0.0;
						for(int coil=0; coil< data.Num_Coils; coil++){
							sos+= norm(smaps(coil)(i,j,k));
						}
						sos = 1./sqrtf(sos);
						for(int coil=0; coil< data.Num_Coils; coil++){
							smaps(coil)(i,j,k) *= sos;
						}
			}}}

		} // Normalization

		// Export 
		if(export_smaps==1){
			cout << "Exporting Smaps" << endl;
			for(int coil=0; coil< smaps.length(firstDim); coil++){
				char name[256];
				sprintf(name,"SenseMaps_%2d.dat",coil);
				ArrayWrite(smaps(coil),name);
			}
		}
	}else if(recon_type != CLEAR){
		// Allocate Storage for Map	and zero	
		cout << "Allocate Sense Maps"  << endl << flush;
		smaps = Alloc4DContainer< complex<float> >(rcxres,rcyres,rczres,data.Num_Coils);
		for(int coil=0; coil< smaps.length(firstDim); coil++){
			smaps(coil)=complex<float>(1.0,0.0);
		}	
	}

} 








