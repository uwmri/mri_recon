/* 

   Recon Code for Cartesian and Non-Cartesian Data
   (in process of reorganizing for template + 5D recons)



 */

#include "wavelet3D.h"	
#include "temporal_diff.h"	
#include "gridFFT.h"
// #include "ge_pfile_lib.h"
#include "recon_lib.h"
#include "mri_data.h"
#include "softthreshold.h"
#include "ArrayTemplates.cpp"
#include "tictoc.cpp"
using namespace std;

Array< complex<float>, 5 >reconstruction( int argc, char **argv, MRI_DATA data,RECON recon);

int main(int argc, char **argv){

	// --------------------------
	// Check for help message and output help
	// --------------------------
	for(int pos=0;pos<argc;pos++){
		if( (strcmp(argv[pos],"-h")==0) || (strcmp(argv[pos],"-help")==0) || (strcmp(argv[pos],"--help")==0)){
			RECON::help_message();	
			exit(0);
		}
	}
	
	// ------------------------------------
	// Setup Recon
	// ------------------------------------

	RECON recon(argc,argv);
	MRI_DATA data;

	cout << "----Read Data-----" << endl;	
	if(recon.data_type==RECON_PFILE){	
		// Read in P-File (doesn't work)
		//PFILE pfile;
		//pfile.read_header(recon.filename);
		//pfile.read_data(0);
	}else{
		// Read in External Data Format
		recon.parse_external_header();
		data.read_external_data("./",recon.num_coils,recon.rcencodes,recon.num_slices,recon.num_readouts,recon.xres);
	}

	cout << "----Geometry Modification-----" << endl;	
	// Geometry Modification by Recon
	data.kx *= ((float)(1.0/recon.zoom_x));
	data.ky *= ((float)(1.0/recon.zoom_y));
	data.kz *= ((float)(1.0/recon.zoom_z));

	if (recon.acc > 1){
		data.undersample(recon.acc);
	}

	if (recon.compress_coils > 0){
		data.coilcompress(recon.compress_coils);
	}

	// Turn of parallel processing for 2D due to thread overhead
	if(recon.rczres ==1){
		omp_set_num_threads(1);
	}

	// --------------------------------------------------
	// Code for recon (no PSD specific data/structures)
	// --------------------------------------------------
	Array< complex<float>,5 >X = reconstruction(argc,argv,data,recon);

	// ------------------------------------
	// Post Processing + Export
	// ------------------------------------

	//Export Binary Images for Now
	for(int ee=0; ee<recon.rcencodes; ee++){
		for(int tt=0; tt<recon.rcframes; tt++){
			char fname[80];
			sprintf(fname,"X_%d_%d.dat.complex",ee,tt);
			Array< complex<float>,3>Xref = X(Range::all(),Range::all(),Range::all(),tt,ee);
			ArrayWrite( Xref,fname);
			sprintf(fname,"X_%d_%d.dat",ee,tt);
			ArrayWriteMag( Xref,fname);
	}}


	return(0);
}

Array< complex<float>,5 >reconstruction( int argc, char **argv, MRI_DATA data,RECON recon){
	// Shorthand for Blitz++
	Range all=Range::all();
	
	// Matlab like timer (openmp code base)
	tictoc T; 

	// Setup Gridding + FFT Structure
	gridFFT gridding;
	gridding.read_commandline(argc,argv);
	gridding.precalc_gridding(recon.rczres,recon.rcyres,recon.rcxres,3);

	// ------------------------------------
	//  Get coil sensitivity map ( move into function)
	// ------------------------------------

	Array<complex<float>,4 >smaps;
	if( (recon.recon_type != RECON_SOS) && (data.Num_Coils>1)){ 
		cout << "Getting Coil Sensitivities " << endl;

		// Low Pass filtering for Sensitivity Map
		gridding.k_rad = 8;

		// Allocate Storage for Map	and zero	
		smaps.setStorage( ColumnMajorArray<4>()); // Hopefully temporary
		smaps.resize(recon.rcxres,recon.rcyres,recon.rczres,data.Num_Coils);
		smaps = 0;

		for(int coil=0; coil< data.Num_Coils; coil++){
			int e =0;

			// Blitz Referencing is a bit wordy
			Array<complex<float>,3>kdataE = data.kdata(all,all,all,e,coil); 
			Array< float,3 >kxE = data.kx(all,all,all,e); 
			Array< float,3 >kyE = data.ky(all,all,all,e); 
			Array< float,3 >kzE = data.kz(all,all,all,e); 
			Array< float,3 >kwE = data.kw(all,all,all,e); 
			
			//Do Gridding			
			gridding.forward( kdataE,kxE,kyE,kzE,kwE);
			
			//Add to 			
			Array< complex<float>,3>SmapC = smaps(all,all,all,coil);	
			SmapC += gridding.image;
		}

		// Restore Full Resolution
		gridding.k_rad = 9999;

		// E-Spirit Code in seperate branch -- need to talk to Michael Loecher about merging


		// Sos Normalization S(coil)= I(coil)/sum(I^2,coils)
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


		// Export need to add flag "-export_smaps"
		if(1==0){
			ArrayWrite(smaps,"SenseMaps.dat");
		}
	}else if(recon.recon_type != RECON_SOS){
		// Allocate Storage for Map	and zero (simplifies IST/PILS code)	
		smaps.setStorage( ColumnMajorArray<4>()); // Hopefully temporary
		smaps.resize(recon.rcxres,recon.rcyres,recon.rczres,data.Num_Coils);
		smaps = complex<float>(1.0,0);
	}


	// ------------------------------------
	//  If time resolved need to sort the times in to bins (need to move to function calls)
	// ------------------------------------

	if(recon.rcframes>1){
		float max_time =max(data.times);
		float min_time =min(data.times);

		cout << "Time Range :: " << min_time << " to " << max_time << endl;

		float delta_t = (max_time - min_time ) / ((float)(recon.rcframes));

		// Sort times into discrete frames
		int *frame_count = new int[recon.rcframes];
		memset( (void *)frame_count,0,(size_t)((recon.rcframes)*sizeof(int)));
		for(int e=0; e< data.times.length(3); e++){
			for(int k=0; k< data.times.length(2); k++){
				for(int j=0; j< data.times.length(1); j++){
					for(int i=0; i< data.times.length(0); i++){

						data.times(i,j,k,e) -= min_time;
						data.times(i,j,k,e) /= delta_t;
						int pos = (int)data.times(i,j,k,e);
						if(pos > (recon.rcframes-1)){
							pos= (recon.rcframes-1);
						}
						data.times(i,j,k,e) = (float)pos;

						frame_count[pos]++;
					}}}}

		for(int t =0; t<recon.rcframes; t++){
			cout << "Frame " << t << " count " << frame_count[t] << endl;
		}
		delete [] frame_count;
	}


	/*----------------------------Main Recons---------------------------------------*/	

	// Final Image Solution
	Array< complex<float>,5 >X(recon.rcxres,recon.rcyres,recon.rczres,recon.rcframes,recon.rcencodes,ColumnMajorArray<5>());
	X=0;

	// Weighting Array for Time coding
	Array< float, 3 >TimeWeight(data.kx.length(0),data.kx.length(1),data.kx.length(2),ColumnMajorArray<3>());

	switch(recon.recon_type){
		default:
		case(RECON_SOS):
		case(RECON_PILS):{

					 for(int e=0; e< recon.rcencodes; e++){
						 for(int t=0; t< recon.rcframes; t++){
							 cout << "Recon Encode" << e << " Frame " << t << endl;

							 // Get Sub-Arrays for Encoding (Blitz-Subarray reference, no memory copied)
							 Array< float,3 >kxE = data.kx(all,all,all,e); 
							 Array< float,3 >kyE = data.ky(all,all,all,e); 
							 Array< float,3 >kzE = data.kz(all,all,all,e); 
							 Array< float,3 >kwE = data.kw(all,all,all,e); 
							 Array< float,3 >timesE = data.times(all,all,all,e); 

							 // Temporal weighting (move to functions )
							 TimeWeight = kwE;
							 if(recon.rcframes>1){
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
								 
								 //Gridding/FFT
								 T.tic();
								 gridding.forward( kdataE,kxE,kyE,kzE,TimeWeight);
								 cout << "\tGridding took = " << T << endl;
								 
								 // Multiply by Sensitivity Map
								 T.tic();
								 if(recon.recon_type==RECON_PILS){
									 Array<complex<float>,3>smapC = smaps(all,all,all,coil);
									 gridding.image *=conj(smapC);
								 }else{
									 gridding.image *=conj(gridding.image);
								 }								 
								 cout << "Conj took = " << T << endl;
								 
								 // Add to accumulated image
								 T.tic();
								 Array<complex<float>,3>xet = X(all,all,all,t,e);
								 xet += gridding.image;
								 cout << "\tAdd to X = " << T << endl;
							 }
							 cout << endl;
						 }
					 }

					 // Take Square Root for SOS	
					 if(recon.recon_type==RECON_SOS){
						 X = csqrt(X);
					 }

				 }break;

		case(RECON_CG):{
				       // ------------------------------------
				       // Conjugate Gradient Recon not yet
				       // ------------------------------------



			       }


		case(RECON_IST):
		case(RECON_FISTA):{

					  // ------------------------------------
					  // Iterative Soft Thresholding  x(n+1)=  thresh(   x(n) - E*(Ex(n) - d)  )
					  //  Designed to not use memory
					  // Uses gradient descent x(n+1) = x(n) - ( R'R ) / ( R'E'E R) * Grad  [ R = E'(Ex-d)]
					  // ------------------------------------

					  // Previous Array for FISTA
					  Array< complex<float>,5>X_old;
					  if( recon.recon_type == RECON_FISTA){
						  X_old.setStorage(ColumnMajorArray<5>());
						  X_old.resize( X.shape());				  
						  X_old = 0.0;
					  }

					  // Residue 	
					  Array< complex<float>,5>R(X.shape(),ColumnMajorArray<5>());
					  R=0.0;

					  // Temp variable for E'ER 
					  Array< complex<float>,3 >P(recon.rcxres,recon.rcyres,recon.rczres,ColumnMajorArray<3>());

					  // Storage for (Ex-d)
					  Array< complex<float>,3 >diff_data(data.kdata.length(0),data.kdata.length(1),data.kdata.length(2),ColumnMajorArray<3>());

					  // Setup 3D Wavelet
					  int dirs[3] = {4, 4, 4};
					  Array< complex<float>,3>Xref=X(all,all,all,0,0);
					  WAVELET3D wave(Xref,dirs,WAVE_DB4);

					  // Temporal differences or FFT
					  TDIFF tdiff(X);

					  // Setup Soft Thresholding
					  SOFTTHRESHOLD softthresh(argc,argv);

					  cout << "Iterate" << endl;
					  double error0=0.0;
					  for(int iteration =0; iteration< recon.max_iter; iteration++){

						  tictoc iteration_timer;
						  iteration_timer.tic();
						  cout << "\nIteration = " << iteration << endl;

						  // Update X based on FISTA 
						  if(recon.recon_type==RECON_FISTA){
						  	  softthresh.fista_update(X,X_old,iteration);
						  }
						  
						  // Zero this for Cauchy set size
						  complex<float>scale_RhP(0,0);						  

						  // Get Residue
						  R=0; 								
						  cout << "\tGradient Calculation" << endl;
						  for(int e=0; e< recon.rcencodes; e++){
							  for(int t=0; t< recon.rcframes; t++){
								  
								  // Get Sub-Arrays for Encoding
								  Array< float,3 >kxE = data.kx(all,all,all,e); 
								  Array< float,3 >kyE = data.ky(all,all,all,e); 
								  Array< float,3 >kzE = data.kz(all,all,all,e); 
								  Array< float,3 >kwE = data.kw(all,all,all,e); 
								  Array< float,3 >timesE = data.times(all,all,all,e); 

								  // Temporal weighting
								  TimeWeight = kwE;
								  if(recon.rcframes>1){
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
									  gridding.image = Xref;
									  gridding.image*= smapC;
									  gridding.backward( diff_data,kxE,kyE,kzE,TimeWeight);

									  //Ex-d
									  Array< complex<float>,3>kdataC = data.kdata(all,all,all,e,coil); 
									  diff_data -= kdataC;

									  //E'(Ex-d)
									  gridding.forward( diff_data,kxE,kyE,kzE,TimeWeight);
									  gridding.image*=conj(smapC);
									  Rref += gridding.image;
								  }//Coils


								  //Now Get Scale factor (for Cauchy-Step Size)
								  P=0;
								  for(int coil=0; coil< data.Num_Coils; coil++){
									  // EE'(Ex-d)
									  diff_data=0;
									  Array<complex<float>,3>smapC=smaps(all,all,all,coil);
									  gridding.image = Rref;
									  gridding.image*= smapC;
									  gridding.backward( diff_data,kxE,kyE,kzE,TimeWeight);

									  //E'EE'(Ex-d)
									  gridding.forward( diff_data,kxE,kyE,kzE,TimeWeight);
									  gridding.image*=conj(smapC);
									  P += gridding.image;
								  }//Coils
								  P*=conj(Rref);
								  scale_RhP += sum(P); 
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

						  // Export X slice
						  Array<complex<float>,2>Xslice=X(all,all,X.length(2)/2,0,0);
						  ArrayWriteMag(Xslice,"X.dat");

						  // ------------------------------------
						  // Soft thresholding operation (need to add transform control)
						  // ------------------------------------

						  if(softthresh.thresh > 0.0){
							  tdiff.fft_t(X);
							  cout << "Wavelet " << endl;
							  wave.random_shift();
							  wave.forward(X);	

							  softthresh.get_threshold(X);
							  softthresh.soft_threshold(X);
							  wave.backward(X);
							  cout << "Wavelet Done" << endl;
							  tdiff.ifft_t(X);

						  }

					  }// Iteration			

				  }break;

	}//Recon Type


	cout << "Recon was completed successfully " << endl;	
	return(X);
}



