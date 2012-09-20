/* 

   Recon Code for Cartesian and Non-Cartesian Data
   (in process of reorganizing for template + 4D recons)

 
 */

#include "wavelet3D.h"	
#include "temporal_diff.h"	
#include "gridFFT.h"
#include "ge_pfile_lib.h"
#include "recon_lib.h"
#include "softthreshold.h"
#include "ArrayTemplates.cpp"
#include "tictoc.cpp"
using namespace std;

double gettime(void);
array5D< complex<float> >reconstruction( int argc, char **argv, MRI_DATA data,RECON recon);

int main(int argc, char **argv){

	// ------------------------------------
	// Setup Recon
	// ------------------------------------

	RECON recon(argc,argv);
	MRI_DATA data;
	
		
	cout << "----Read Data-----" << endl;	
	if(recon.data_type==RECON_PFILE){	
		// Read in P-File (doesn't work)
		PFILE pfile;
		pfile.read_header(recon.filename);
		pfile.read_data(0);
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
	array5D< complex<float> >X = reconstruction(argc,argv,data,recon);
	
	// ------------------------------------
	// Post Processing + Export
	// ------------------------------------
	
	// Export Complex Images for Now
	for(int ee=0; ee<recon.rcencodes; ee++){
			char fname[80];
			sprintf(fname,"X_%d.dat",ee);
			X[ee][0].write_mag(fname);
					
			sprintf(fname,"X_%d.complex.dat",ee);
			X[ee][0].write(fname);
	}


	return(0);
}




array5D< complex<float> >reconstruction( int argc, char **argv, MRI_DATA data,RECON recon){


	// Setup Gridding + FFT Structure
	gridFFT gridding;
	gridding.read_commandline(argc,argv);
	gridding.precalc_gridding(recon.rczres,recon.rcyres,recon.rcxres,3);

	// ------------------------------------
	//  Get coil sensitivity map
	// ------------------------------------

	array4D< complex<float> >smaps;
	if( recon.recon_type != RECON_SOS){ 
		cout << "Getting Coil Sensitivities " << endl;

		
		gridding.k_rad = 16;
		smaps.alloc(data.Num_Coils,recon.rczres,recon.rcyres,recon.rcxres);
		smaps.zero();
		for(int coil=0; coil< data.Num_Coils; coil++){
				gridding.forward( data.kdata[coil][0],data.kx[0],data.ky[0],data.kz[0],data.kw[0]);
				smaps[coil] =( gridding.return_array() );
		}
		gridding.k_rad = 9999;

		// Spirit Code?


		// Sos Normalization 
		cout << "Normalize Coils" << endl;
#pragma omp parallel for 
		for(int k=0; k<smaps[0].Nz; k++){
			for(int j=0; j<smaps[0].Ny; j++){
				for(int i=0; i<smaps[0].Nx; i++){
					float sos=0.0;
					for(int coil=0; coil< data.Num_Coils; coil++){
						sos+= norm(smaps[coil][k][j][i]);
					}
					sos = 1./sqrtf(sos);
					for(int coil=0; coil< data.Num_Coils; coil++){
						smaps[coil][k][j][i] *= sos;
					}

		}}}
		
		
		// Export
		if(0==1){
			for(int coil=0; coil< data.Num_Coils; coil++){
				char fname[80];
				sprintf(fname,"SMAP_C%d.dat",coil);
				smaps[coil].write(fname);
			}
		}
	}
	

	/*----------------------------Main Recons---------------------------------------*/	
	
	// Final Image Solution
	array5D< complex<float> >X;
	X.alloc(recon.rcencodes,recon.rcframes,recon.rczres,recon.rcyres,recon.rcxres);
	X.zero();
					
	switch(recon.recon_type){
		default:
		case(RECON_SOS):{

					// ------------------------------------
					// Sum of Squares Recon (not for PC)
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

		case(RECON_PILS):{

					 // ------------------------------------
					 // PILS
					 // ------------------------------------
					 
					 tictoc T;
					 X.zero();
					 for(int e=0; e< recon.rcencodes; e++){
						for(int t=0; t< recon.rcframes; t++){
							cout << "\tForward Gridding Coil ";
							for(int coil=0; coil< data.Num_Coils; coil++){
								cout << coil << "," << flush;
								T.tic();
								gridding.forward( data.kdata[coil][e],data.kx[e],data.ky[e],data.kz[e],data.kw[e]);
								cout << "\tGridding took = " << T << endl;
								
								T.tic();
								gridding.image.conjugate_multiply(smaps[coil]);
								cout << "\tConj Multiple took = " << T << endl;
								
								T.tic();
								X[e][t] += gridding.image;
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
					int dirs[5] = {4, 4, 4};
					WAVELET3D wave(&X,dirs,WAVE_DB4);
					
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
							X[0][0].write_mag("X.dat",X.Nz/2,"a+"); // Just one slice from one encoding
							
							float A = 1.0 + (iteration - 1.0)/(iteration+1.0);
							float B =     - (iteration - 1.0)/(iteration+1.0);
																					
							// Get Update
							for(int jj=0; jj<X[0][0].Numel; jj++){
								complex<float>Xn0 = X_old[0][0](jj);								
								complex<float>Xn1 = X[0][0](jj);
								X[0][0](jj) = A*Xn1 + B*Xn0;
								X_old[0][0](jj) = Xn1;
							}
							
							cout << "FISTA: A = " <<  A << " B = " << B << endl;
							X[0][0].write_mag("X.dat",X.Nz/2,"a+"); // Just one slice from one encoding
						}
												
						// Ex 								
						diff_data.zero();
						for(int e=0; e< recon.rcencodes; e++){
							for(int t=0; t< recon.rcframes; t++){
								cout << "\tInverse Gridding Coil ";
								for(int coil=0; coil< data.Num_Coils; coil++){
									cout << coil << "," << flush;
									gridding.image.multi_equal( X[e][t], smaps[coil]);
									gridding.backward( diff_data[coil][e],data.kx[e],data.ky[e],data.kz[e],data.kw[e]);
								}
								cout << endl;
							}
						}
						
						// Ex-d
						diff_data -= data.kdata;
						
						cout << "Energy in difference = " << diff_data.energy() << endl;
						
						
						// E'(Ex-d)
						R.zero();
						for(int e=0; e< recon.rcencodes; e++){
							for(int t=0; t< recon.rcframes; t++){
								cout << "\tForward Gridding Coil ";
								for(int coil=0; coil< data.Num_Coils; coil++){
									cout << coil << "," << flush;
									gridding.forward( diff_data[coil][e],data.kx[e],data.ky[e],data.kz[e],data.kw[e]);
									gridding.image.conjugate_multiply(smaps[coil]);
									R[e][t] += gridding.image;
								}
								cout << endl;
							}
						}
						
						// Get Scaling Factor
						complex<float>scale_RhP(0,0);
						for(int e=0; e< recon.rcencodes; e++){
							for(int t=0; t< recon.rcframes; t++){
						
								// Compute Ex of Residue EE'(Ex-d)
								cout << "\tInverse Gridding Coil ";
								diff_data.zero();
								for(int coil=0; coil< data.Num_Coils; coil++){
									cout << coil << "," << flush;
									gridding.image.multi_equal( R[e][t], smaps[coil]);
									gridding.backward( diff_data[coil][e],data.kx[e],data.ky[e],data.kz[e],data.kw[e]);
								}
								cout << endl;

								// Compute the Residue  E'E ( E'(Ex-d))
								cout << "\tForward Gridding Coil ";
								P.zero();
								for(int coil=0; coil< data.Num_Coils; coil++){
									cout << coil << "," << flush;
									gridding.forward( diff_data[coil][e],data.kx[e],data.ky[e],data.kz[e],data.kw[e]);
									gridding.image.conjugate_multiply(smaps[coil]);
									P += gridding.image;
								}
								cout << endl;

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
						complex<float>scale = scale_RhR/scale_RhP;
						cout << "Scale = " << scale << endl;
						R *= scale;
						X -= R;
						cout << "Took " << (gettime()-start) << " s " << endl;

						// ------------------------------------
					  	// Soft thresholding operation
					  	// ------------------------------------
						
						// Export X		
						for(int ee=0; ee<recon.rcencodes; ee++){
							X[ee][0].write_mag("X.dat",X.Nz/2,"a+"); // Just one slice from one encoding
							X[ee][0].write_phase("X_Phase.dat",X.Nz/2,"a+"); // Just one slice from one encoding
						}
						
						if(softthresh.thresh > 0.0){
							
							tdiff.fft_e(X);
							cout << "Wavelet " << endl;
							wave.random_shift();
							wave.forward();	
						
							if(iteration==1){
								softthresh.get_threshold(X);
							}
							softthresh.soft_threshold(X);
							wave.backward();
							cout << "Wavelet Done" << endl;
							tdiff.ifft_e(X);
						}
						
						
					}// Iteration			
					

				}break;

	}//Recon Type

	
	cout << "Recon was completed successfully " << endl;	
	return(X);
}


