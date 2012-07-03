/* 

   Recon Code for Cartesian and Non-Cartesian Data
   (in process of reorganizing for template + 4D recons)

 
 */

#include "wavelet3D.h"			 
#include "gridFFT.h"
#include "ge_pfile_lib.h"
#include "recon_lib.h"
#include "softthreshold.h"
#include "ArrayTemplates.cpp"
using namespace std;

double gettime(void);

int main(int argc, char **argv){

	// ------------------------------------
	// Setup Recon
	// ------------------------------------

	RECON recon(argc,argv);
	PFILE pfile;
	MRI_DATA data;
	
	cout << "----Read Data-----" << endl;	
	if(recon.data_type==RECON_PFILE){	
		// Read in P-File
		pfile.read_header(recon.filename);
		pfile.read_data(0);
	}else{
		// Read in External Data Format
		recon.parse_external_header();
		data.read_external_data("./",recon.num_coils,recon.rcencodes,recon.num_readouts,recon.xres);
	}

	cout << "----Geometry Modification-----" << endl;	
	// Geometry Modification by Recon
	data.kx *= ((float)(1.0/recon.zoom_x));
	data.ky *= ((float)(1.0/recon.zoom_y));
	data.kz *= ((float)(1.0/recon.zoom_z));

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

		int e=0;
		gridding.k_rad = 24;
		smaps.alloc(data.Num_Coils,recon.rczres,recon.rcyres,recon.rcxres);
		for(int coil=0; coil< data.Num_Coils; coil++){
			gridding.forward( data.kdata[coil][e][0],data.kx[e][0],data.ky[e][0],data.kz[e][0],data.kw[e][0],data.Num_Pts*data.Num_Readouts);
			smaps[coil]= ( gridding.return_array() );
		}
		gridding.k_rad = 99999;

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
	}
	
	

	/*----------------------------Main Recons---------------------------------------*/	

	switch(recon.recon_type){
		default:
		case(RECON_SOS):{

					// ------------------------------------
					// Sum of Squares Recon 
					// ------------------------------------
					array3D< float >sos;
					sos.alloc(recon.rczres,recon.rcyres,recon.rcxres);
					int e=0;
					for(int coil=0; coil< data.Num_Coils; coil++){
						cout << "Gridding Coil " << coil << endl;
						gridding.forward( data.kdata[coil][e][0],data.kx[e][0],data.ky[e][0],data.kz[e][0],data.kw[e][0],data.Num_Pts*data.Num_Readouts);
						for(int ii=0; ii< sos.Numel; ii++){
							sos(ii) += norm( gridding.image(ii));	// adds square
						}
					}
					sos.sqrt();
					sos.write("Sos.dat");
					
					}break;

		case(RECON_PILS):{

					 // ------------------------------------
					 // PILS
					 // ------------------------------------
					 array3D< complex<float> >pils;
					 pils.alloc(recon.rczres,recon.rcyres,recon.rcxres);
					 int e =0;
					 for(int coil=0; coil< data.Num_Coils; coil++){
						 cout << "Gridding Coil " << coil << endl;
						 gridding.forward( data.kdata[coil][e][0],data.kx[e][0],data.ky[e][0],data.kz[e][0],data.kw[e][0],data.Num_Pts*data.Num_Readouts);
						 gridding.image.conjugate_multiply(smaps[coil]);
						 pils += ( gridding.return_array() );
					 }
					 pils.write("PILS.dat");

				 }break;

		case(RECON_CG):{
				  // ------------------------------------
				  // Conjugate Gradient Recon
				  // ------------------------------------



			       }

		case(RECON_IST):{

					// ------------------------------------
					// Iterative Soft Thresholding  x(n+1)=  thresh(   x(n) - E*(Ex(n) - d)  )
					//  Designed to not use memory
					// Uses gradient descent x(n+1) = x(n) - ( R'R ) / ( R'E'E R) * Grad  [ R = E'(Ex-d)]
					// ------------------------------------

					// Final Image Solution
					array5D< complex<float> >X;
					X.alloc(recon.rcencodes,recon.rcframes,recon.rczres,recon.rcyres,recon.rcxres);
					X.zero();

					// Residue 	
					array5D< complex<float> >R;
					R.alloc(recon.rcencodes,recon.rcframes,recon.rczres,recon.rcyres,recon.rcxres);
					R.zero();

					// Temp variable for E'ER 
					array3D< complex<float> >P;
					P.alloc(recon.rczres,recon.rcyres,recon.rcxres);
					P.zero();

					// Storage for (Ex-d)
					array4D< complex<float> >diff_data;
					diff_data.samesize(&data.kdata);
				
					// Setup 5D Wavelet
					int dirs[5] = {4, 4, 4, 0, 3};
					int waves[5] = {WAVE_DB4, WAVE_DB4, WAVE_DB4, WAVE_DB2, WAVE_DB2};
					WAVELET3D wave(&X,dirs,waves);
					
					// Setup Soft Thresholding
					SOFTTHRESHOLD softthresh(argc,argv);
										
					cout << "Iterate" << endl;
					double error0=0.0;
					for(int iteration =0; iteration<50; iteration++){

						double start = gettime();
						cout << "\nIteration = " << iteration << endl;
						diff_data.zero();
						
						// Ex 								
						for(int e=0; e< recon.rcencodes; e++){
							for(int t=0; t< recon.rcframes; t++){
								cout << "\tInverse Gridding Coil ";
								for(int coil=0; coil< data.Num_Coils; coil++){
									cout << coil << "," << flush;
									gridding.image  = X[e][t];
									gridding.image *=smaps[coil];
									gridding.backward( diff_data[coil][e][0],data.kx[e][0],data.ky[e][0],data.kz[e][0],data.kw[e][0],data.Num_Pts*data.Num_Readouts);
								}
								cout << endl;
							}
						}
						
						// Ex-d
						diff_data -= data.kdata;
						
						// E'(Ex-d)
						R.zero();
						for(int e=0; e< recon.rcencodes; e++){
							for(int t=0; t< recon.rcframes; t++){
								cout << "\tForward Gridding Coil ";
								for(int coil=0; coil< data.Num_Coils; coil++){
									cout << coil << "," << flush;
									gridding.forward( diff_data[coil][e][0],data.kx[e][0],data.ky[e][0],data.kz[e][0],data.kw[e][0],data.Num_Pts*data.Num_Readouts);
									gridding.image.conjugate_multiply(smaps[coil]);
									R[e][t] += ( gridding.return_array() );
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
									gridding.image  = R[e][t];
									gridding.image *=smaps[coil];
									gridding.backward( diff_data[coil][e][0],data.kx[e][0],data.ky[e][0],data.kz[e][0],data.kw[e][0],data.Num_Pts*data.Num_Readouts);
								}
								cout << endl;

								// Compute the Residue  E'E ( E'(Ex-d))
								cout << "\tForward Gridding Coil ";
								P.zero();
								for(int coil=0; coil< data.Num_Coils; coil++){
									cout << coil << "," << flush;
									gridding.forward( diff_data[coil][e][0],data.kx[e][0],data.ky[e][0],data.kz[e][0],data.kw[e][0],data.Num_Pts*data.Num_Readouts);
									gridding.image.conjugate_multiply(smaps[coil]);
									P += ( gridding.return_array() );
								}
								cout << endl;

								// Compute Scale  scale =  [ (Ex-d)'(Ex-d) ] / [ (Ex-d)' EE' (Ex-d) ]
								//						= (r'r)/( r'EE'r) where r is k-space residue
								P.conjugate_multiply(R[e][t]);
								scale_RhP += P.sum();

							}// Time Frame
						}// Encoding
						complex<float>scale_RhR = complex<float>(R.energy(),0);
						if(iteration==0){
							error0 = abs( scale_RhR);
						}
						cout << "Residue Energy =" << scale_RhR << " ( " << (abs(scale_RhR)/error0) << " % )  " << endl;
						
						// Step 
						complex<float>scale = scale_RhR/scale_RhP;
						cout << "Scale = " << scale << endl;
						R *= scale;
						X -= R;
						cout << "Took " << (gettime()-start) << " s " << endl;

						// X	
						FILE *fid;
						for(int ee=0; ee<recon.rcencodes; ee++){
							fid = fopen("X.dat","a+");
							fwrite( X[ee][0][X.Nz/2][0],X.Nx*X.Ny, sizeof(complex<float>),fid);
							fclose(fid);
						}

						// ------------------------------------
					  	// Soft thresholding operation
					  	// ------------------------------------
						cout << "Wavelet " << endl;
						wave.random_shift();
						wave.forward();	
						
						// Export X		
						for(int ee=0; ee<recon.rcencodes;ee++){
							fid = fopen("X.dat","a+");
							fwrite( X[ee][0][X.Nz/2][0],X.Nx*X.Ny, sizeof(complex<float>),fid);
							fclose(fid);
						}

						softthresh.hard_threshold(X);
						
						wave.backward();
						cout << "Wavelet Done" << endl;
						
						// Export X		
						for(int ee=0; ee<recon.rcencodes;ee++){
							fid = fopen("X.dat","a+");
							fwrite( X[ee][0][X.Nz/2][0],X.Nx*X.Ny, sizeof(complex<float>),fid);
							fclose(fid);
						}

					}// Iteration			
					
					// Export Complex Images for Now
					for(int ee=0; ee<recon.rcencodes; ee++){
						char fname[80];
						sprintf(fname,"IST_%d.dat",ee);
						X[ee][0].write(fname);
					}
				}break;

	}//Recon Type



	// ------------------------------------
	// Post Processing + Export
	// ------------------------------------



	cout << "Recon was completed successfully " << endl;	
	return(0);
}


