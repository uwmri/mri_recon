/* 

   Recon Code for Cartesian and Non-Cartesian Data
   (in process of reorganizing for template + 4D recons)


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


double gettime(void);
Array< complex<float>, 5 >reconstruction( int argc, char **argv, MRI_DATA data,RECON recon);

int main(int argc, char **argv){

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

	// Turn of parallel processing for 2D
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

	//Export Complex Images for Now
	for(int ee=0; ee<recon.rcencodes; ee++){
		for(int tt=0; tt<recon.rcframes; tt++){
			char fname[80];
			sprintf(fname,"X_%d_%d.dat.complex",ee,tt);
			Array< complex<float>,3>Xref = X(Range::all(),Range::all(),Range::all(),tt,ee);
			sprintf(fname,"X_%d_%d.dat",ee,tt);
			ArrayWriteMag( Xref,fname);
	}}


	return(0);
}

Array< complex<float>,5 >reconstruction( int argc, char **argv, MRI_DATA data,RECON recon){


	// Setup Gridding + FFT Structure
	gridFFT gridding;
	gridding.read_commandline(argc,argv);
	gridding.precalc_gridding(recon.rczres,recon.rcyres,recon.rcxres,3);

	// ------------------------------------
	//  Get coil sensitivity map
	// ------------------------------------

	Array<complex<float>,4 >smaps;
	if( (recon.recon_type != RECON_SOS) && (data.Num_Coils>1)){ 
		cout << "Getting Coil Sensitivities " << endl;

		gridding.k_rad = 8;
		smaps.setStorage( ColumnMajorArray<4>());
		smaps.resize(recon.rcxres,recon.rcyres,recon.rczres,data.Num_Coils);
		smaps = 0;
		
		for(int coil=0; coil< data.Num_Coils; coil++){
			// Blitz Referencing is a bit wordy
			Array< complex<float>,3>SmapC = smaps(Range::all(),Range::all(),Range::all(),coil);	
			int e =0;
			Array<complex<float>,3>kdataE = data.kdata(Range::all(),Range::all(),Range::all(),e,coil); 
			Array< float,3 >kxE = data.kx(Range::all(),Range::all(),Range::all(),e); 
			Array< float,3 >kyE = data.ky(Range::all(),Range::all(),Range::all(),e); 
			Array< float,3 >kzE = data.kz(Range::all(),Range::all(),Range::all(),e); 
			Array< float,3 >kwE = data.kw(Range::all(),Range::all(),Range::all(),e); 
			gridding.forward( kdataE,kxE,kyE,kzE,kwE);
			SmapC += gridding.image;
		}
		gridding.k_rad = 9999;

		// Spirit Code?


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


		// Export
		if(1==0){
			ArrayWrite(smaps,"SenseMaps.dat");
		}
	}

	// ------------------------------------
	//  If time resolved need to sort the times in to bins
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
	
	Array< float, 3 >TimeWeight(data.kx.length(0),data.kx.length(1),data.kx.length(2),ColumnMajorArray<3>());
	
	switch(recon.recon_type){
		default:

#ifdef LKJKL
		case(RECON_SOS):{

					// ------------------------------------
					// Sum of Squares Recon (only for testing)
					// ------------------------------------

					X.zero();
					for(int e=0; e< recon.rcencodes; e++){
						for(int t=0; t< recon.rcframes; t++){
							cout << "\tForward Gridding Coil ";
							for(int coil=0; coil< data.Num_Coils; coil++){
								cout << coil << "," << flush;
								gridding.forward( data.kdata[coil][e],data.kx[e],data.ky[e],data.kz[e],data.kw[e]);
								gridding.image.conjugate_multiply(gridding.image);
								X[e][t] += gridding.image;
							}
							cout << endl;
						}
					}

				}break;
#endif
		case(RECON_PILS):{

					 // ------------------------------------
					 // PILS
					 // ------------------------------------

					 tictoc T;
					 Range all=Range::all();
					 for(int e=0; e< recon.rcencodes; e++){
						 for(int t=0; t< recon.rcframes; t++){
							 cout << "Recon Encode" << e << " Frame " << t << endl;
							 
							 // Get Sub-Arrays for Encoding
							 Array< float,3 >kxE = data.kx(all,all,all,e); 
							 Array< float,3 >kyE = data.ky(all,all,all,e); 
							 Array< float,3 >kzE = data.kz(all,all,all,e); 
							 Array< float,3 >kwE = data.kw(all,all,all,e); 
							 Array< float,3 >timesE = data.times(all,all,all,e); 
							
							 // Temporal weighting
							 TimeWeight = kwE;
							 if(recon.rcframes>1){
								  cout << "Set Times" << endl << flush;
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
								 Array<complex<float>,3>kdataE = data.kdata(all,all,all,e,coil); 
							
								 cout << coil << "," << flush;
								 T.tic();
								 gridding.forward( kdataE,kxE,kyE,kzE,TimeWeight);
								 cout << "\tGridding took = " << T << endl;
								
								 T.tic();
								 Array<complex<float>,3>smapC = smaps(all,all,all,coil);
								 gridding.image*=conj(smapC);
								 cout << "Conj took = " << T << endl;

								 T.tic();
								 Array<complex<float>,3>xet = X(all,all,all,t,e);
								 xet += gridding.image;
								 cout << "\tAdd to X = " << T << endl;
							 }
							 cout << endl;
						 }
					 }


				 }break;

		case(RECON_CG):{
				       // ------------------------------------
				       // Conjugate Gradient Recon
				       // ------------------------------------



			       }


		case(RECON_IST):
		case(RECON_FISTA):{

#ifdef JLLKJL
					  // ------------------------------------
					  // Iterative Soft Thresholding  x(n+1)=  thresh(   x(n) - E*(Ex(n) - d)  )
					  //  Designed to not use memory
					  // Uses gradient descent x(n+1) = x(n) - ( R'R ) / ( R'E'E R) * Grad  [ R = E'(Ex-d)]
					  // ------------------------------------

					  array5D< complex<float> >X_old;
					  if( recon.recon_type == RECON_FISTA){
						  X_old.alloc(recon.rcencodes,recon.rcframes,recon.rczres,recon.rcyres,recon.rcxres);
						  X_old.zero();
					  }

					  // Residue 	
					  array5D< complex<float> >R;
					  R.alloc(recon.rcencodes,recon.rcframes,recon.rczres,recon.rcyres,recon.rcxres);
					  R.zero();

					  // Temp variable for E'ER 
					  array3D< complex<float> >P;
					  P.alloc(recon.rczres,recon.rcyres,recon.rcxres);
					  P.zero();

					  // Storage for (Ex-d)
					  array5D< complex<float> >diff_data;
					  diff_data.samesize(&data.kdata);

					  // Setup 3D Wavelet
					  int dirs[3] = {4, 4, 4};
					  //WAVELET3D wave(&X,dirs,WAVE_DB4);

					  // Temporal differences or FFT
					  TDIFF tdiff(X);

					  // Setup Soft Thresholding
					  SOFTTHRESHOLD softthresh(argc,argv);

					  cout << "Iterate" << endl;
					  double error0=0.0;
					  for(int iteration =0; iteration< recon.max_iter; iteration++){

						  double start = gettime();
						  cout << "\nIteration = " << iteration << endl;

						  // Update X based on FISTA 
						  if(recon.recon_type==RECON_FISTA){
							  softthresh.fista_update(X,X_old,iteration);
						  }

						  // Ex 								
						  cout << "\tInverse Gridding " << endl;
						  diff_data.zero();
						  for(int e=0; e< recon.rcencodes; e++){
							  for(int t=0; t< recon.rcframes; t++){

								  //Temporal weighting
								  TimeWeight=data.kw[e];
								  if(recon.rcframes>1){
									  for(int pos=0; pos< TimeWeight.Numel;pos++){
										  if( (data.times[e](pos)) != (float)t){
											  TimeWeight(pos)= 0.0;
								  }}}
								  
								  for(int coil=0; coil< data.Num_Coils; coil++){
									  gridding.image.multi_equal( X[e][t], smaps[coil]);
									  gridding.backward( diff_data[coil][e],data.kx[e],data.ky[e],data.kz[e],TimeWeight);
								  }
							  }
						  }

						  // Ex-d
						  diff_data -= data.kdata;


						  // E'(Ex-d)
						  R.zero();
						  cout << "\tForward Gridding Residue " << endl;
						  for(int e=0; e< recon.rcencodes; e++){
							  for(int t=0; t< recon.rcframes; t++){
								  //Temporal weighting
								  TimeWeight=data.kw[e];
								  if(recon.rcframes>1){
									  for(int pos=0; pos< TimeWeight.Numel;pos++){
										  if( (data.times[e](pos)) != (float)t){
											  TimeWeight(pos)= 0.0;
								  }}}
								  
								  for(int coil=0; coil< data.Num_Coils; coil++){
									  gridding.forward( diff_data[coil][e],data.kx[e],data.ky[e],data.kz[e],TimeWeight);
									  gridding.image.conjugate_multiply(smaps[coil]);
									  R[e][t] += gridding.image;
								  }
							  }
						  }

						  // Get Scaling Factor
						  complex<float>scale_RhP(0,0);
						  cout << "\tInverse Gridding Residue " << endl;
						  for(int e=0; e< recon.rcencodes; e++){
							  for(int t=0; t< recon.rcframes; t++){
								  //Temporal weighting
								  TimeWeight=data.kw[e];
								  if(recon.rcframes>1){
									  for(int pos=0; pos< TimeWeight.Numel;pos++){
										  if( (data.times[e](pos)) != (float)t){
											  TimeWeight(pos)= 0.0;
								  }}}
								  	
								  // Compute Ex of Residue EE'(Ex-d)
								  diff_data.zero();
								  for(int coil=0; coil< data.Num_Coils; coil++){
									  gridding.image.multi_equal( R[e][t], smaps[coil]);
									  gridding.backward( diff_data[coil][e],data.kx[e],data.ky[e],data.kz[e],TimeWeight);
								  }
								  
								  // Compute the Residue  E'E ( E'(Ex-d))
								  P.zero();
								  for(int coil=0; coil< data.Num_Coils; coil++){
									  gridding.forward( diff_data[coil][e],data.kx[e],data.ky[e],data.kz[e],TimeWeight);
									  gridding.image.conjugate_multiply(smaps[coil]);
									  P += gridding.image;
								  }
								  
								  // Compute Scale  scale =  [ (Ex-d)'(Ex-d) ] / [ (Ex-d)' EE' (Ex-d) ]
								  //						= (r'r)/( r'EE'r) where r is k-space residue
								  P.conjugate_multiply(R[e][t]);
								  scale_RhP += P.sum();

							  }// Time Frame
						  }// Encoding
						  complex<float>scale_RhR = complex<float>(R.energy(),0);

						  // Error check
						  if(iteration==0){
							  error0 = abs( scale_RhR);
						  }
						  cout << "Residue Energy =" << scale_RhR << " ( " << (abs(scale_RhR)/error0) << " % )  " << endl;

						  // Export X		
						  R[0][0].write_mag("R.dat",R.Nz/2,"a+"); // Just one slice from one encoding

						  // Step in direction
						  complex<float>scale = (scale_RhR/scale_RhP);
						  cout << "Scale = " << scale << endl;
						  R *= scale;
						  X -= R;
						  cout << "Took " << (gettime()-start) << " s " << endl;

						  // ------------------------------------
						  // Soft thresholding operation
						  // ------------------------------------
						  
						  if(softthresh.thresh > 0.0){

							  //tdiff.forward(X);	
							  tdiff.fft_t(X);
							  cout << "Wavelet " << endl;
							  //wave.random_shift();
							  //wave.forward();	
							  
							  // Export X		
						  	for(int ee=0; ee<recon.rcencodes; ee++){
							  for(int t=0; t< recon.rcframes; t++){
								  X[ee][t].write_mag("X.dat",X.Nz/2,"a+"); // Just one slice from one encoding
								  X[ee][t].write_phase("X_Phase.dat",X.Nz/2,"a+"); // Just one slice from one encoding
						  }}
							  							 
							  //if(iteration==1){
							  softthresh.get_threshold(X);
							  //  }
							  softthresh.soft_threshold(X);
							  //wave.backward();
							  cout << "Wavelet Done" << endl;
							  tdiff.ifft_t(X);
							  //tdiff.backward(X);	
						  }
		 				 
						  // Export X		
						  for(int ee=0; ee<recon.rcencodes; ee++){
							  for(int t=0; t< recon.rcframes; t++){
								  X[ee][t].write_mag("X.dat",X.Nz/2,"a+"); // Just one slice from one encoding
								  X[ee][t].write_phase("X_Phase.dat",X.Nz/2,"a+"); // Just one slice from one encoding
						  }}	


					  }// Iteration			

#endif

				  }break;

	}//Recon Type


	cout << "Recon was completed successfully " << endl;	
	return(X);
}



