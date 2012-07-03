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



//------------OLD RECON -------------------------------//




#ifdef NNNN


// Old Code








#include "iterative_lib.h" 
#include "master_lib.h" 
#include "matrix_lib.h" 
#include "pcvipr_gradwarp.h"
#include "io_lib.h"
#include "polynomial_fitting.h"


int main(int argc, char **argv)
{

	double start_time=gettime();
	double map_time=0;
	double recon_time;

	printf("(Start time = %.3f)",start_time - start_time );

	/*GE Header Variables*/
	RDB_HEADER_REC rdbhead;
	RDB_DATA_ACQ_TAB acq_tab;
	EXAM examhead;
	SERIES serieshead;
	IMAGE imagehead;

	/*Phase Storage for Phase Contrast*/
	mtx3d_short **theta=0;		/*Current Phase Data*/
	mtx3d_short **theta_comp=0;	/*Composite Phase Data*/

	/*Regridding (K-space)*/
	cmtx3d *k3d_grid=0; 			/*Overgrid k-space*/
	cmtx3d *image=0;				/*Kspace Bipolar Up */
	cmtx3d *k3dvel=0;				/*Kspace Bipolar down*/
	cmtx3d *k3d_s=0;				/*Kspace Phase Store*/
	cmtx3d *k3d_warp=0;			/*Complex grad warp*/
	kspace_map *kmap=0;			/*holds all info for recon*/

	/*Image Space*/
	mtx3d *time_mip=0;			/* Mip in time*/
	mtx3d *image_comp=0;			/* Composite CD Image*/
	mtx3d *image_mag=0;			/* Composite Mag Image*/
	mtx3d *image_warp=0;         	/* Storage for grad warp*/ 

	/*Raw Data Stores*/
	cmtx *dataset;				/* Time frame data*/

	/*IO Routines*/
	FILE *fp;						/*Output 1*/
	char fname[80];				/*File Name*/
	char fname_vx[80];			/*File Name*/
	char fname_vy[80];			/*File Name*/
	char fname_vz[80];			/*File Name*/
	int im_seno;					/*for dicom output*/
	float max_signal;				/*Max on float image*/
	float scale_factor;			/*Scaling to right bit resolution*/

	/*Temp Variables*/
	int  i, j, k;					/*Positions in 3-orthogonal dir*/
	int recv_cnt=0;				/*Reciever number*/
	int vd=0;						/*Velocity Direction*/
	int n_ph=0;
	int time_frame;
	int vdt;
	float mag_scale = -1;  
	char output_filename[1024];
	int vdgrid=0;
	float *kspace_radius=0;

	/***Coil Sensitivity***/
	float mag;

	/**OFFRESONANCE CORRECTIONS*/ 
	mtx3d_short  *FREQ=0;
	float demod_freq=0;
	float demod_arg=0;

	/**Phase Correction*/
	float min_mag;
	float max_cd;
	float cd;
	float vx, vy, vz;
	float vel;

	/*For Timing*/
	double grid_time;
	double fft_time;
	double igrid_time;
	double ifft_time;
	double read_time;
	double tornado_time;
	double combine_time;
	double mfi_time;
	double deapp_time; 
	double maxwell_time; 
	double L1_regularize_time; 
	double L2_regularize_time; 

	/*Background phase fitting*/
	POLYFIT_PARAMS *polyfit_vx=0;
	POLYFIT_PARAMS *polyfit_vy=0;
	POLYFIT_PARAMS *polyfit_vz=0;

	/**Structures for Libraries defaults + allocation in construct function**/
	ISENSE_DATA *isense_data = construct_isense_data();
	GRIDDING_INFO *grd_info = construct_gridding_info();
	ARB_TRAJ *arb_traj = construct_arb_traj();
	PCFLAG *flag = construct_pcflag();			
	TORNADO_INFO *tornado_info = construct_tornado_flags();

	/***********************************************************************
	  Input
	  -Read inputs for command line
	  -Get parameters from headers
	  -Convert parameters for recon
	 ***********************************************************************/
	printf("\nCommand Line Args: \n");
	for(int apos=0; apos<argc; apos++){
		printf("%s ",argv[apos]);
	}
	printf("\n\n");

	isense_read_commandline(argc,argv, isense_data);			
	tornado_read_commandline(argc,argv, tornado_info);	
	grid_read_commandline(argc,argv, grd_info);					
	master_read_commandline(argc, argv, flag,&im_seno);

	/*Read Header*/
	if(flag->external_data==0){
		read_raw_ge_header(flag->rawfilename, &rdbhead, &acq_tab, &examhead, &serieshead, &imagehead,flag);
		setup_flags( flag, &rdbhead, &imagehead, &serieshead);
		export_scan_info(flag, &examhead, &serieshead, &imagehead, &acq_tab, &rdbhead);
	}else{
		setup_external( flag );
	}
	isense_data->threads  = flag->threads;
	tornado_info->threads = 1;

	printf("Done Reading Command Line\n");	

	/************************************************************************
	  Gridding Prep
	  -Get Gridding deapp window + convolution kernel
	 ***********************************************************************/
	printdbg("Setup Recon Basics\n");
	setup_recon_type( flag, isense_data);
	setup_gridding_lookup_table( grd_info,flag->rczres,flag->rcyres,flag->rcxres,flag->trajectory_type);
	setup_deapp_window( grd_info,flag->rczres,flag->rcyres,flag->rcxres);
	get_offcenter_shift( flag, &imagehead, &serieshead, &rdbhead);
	if(flag->gating!=GATE_NONE){
		setup_gating( flag, &rdbhead);
		setup_default_tornado_filter(tornado_info,flag->time_phase[1] -flag->time_phase[0], flag->time_phase[flag->frames-1] -flag->time_phase[0],flag->res/2);
	}

	/***********************************************************************
	  Trajectory
	  -Get K-space trajectory
	  -Setup density weighting
	 ***********************************************************************/
	printf("Calculate Trajectory\n");
	if(flag->kmap_memory_save==0){
		kmap = kspace_map_alloc(flag->xres, flag->nproj);
		flag->kmap = kmap;

		if(flag->external_data==1){
			read_external_kmap( flag,0);
			kspace_radius = get_kspace_radius( flag );
			get_weighted_echoes( flag);	
		}
	}

	if(flag->external_data ==0){
		if(arb_traj->arb_fov==1){
			get_arb_traj(arb_traj);
		}
		calc_trajectory(flag,&rdbhead, flag->xres, kmap);
		kspace_radius = get_kspace_radius( flag );

		if( (flag->maxwell==1) && (flag->recon_type==PHASE_DIFF) ){
			printf("Getting Maxwell Coeficients from Trajectory\n");
			get_maxwell_coef( flag);
		}

		printf("Exporting Flags for Matlab\n");
		export_flags( flag, examhead.ex_no);
		if(flag->export_matlab_header == 1){
			exit(1);
		}
	}

	/***********************************************************************
	 * Dicom Init - Done Here to Allow Use for Grad Warp                               
	 ***********************************************************************/
#ifdef UW_DICOM 	
	printf("Initialize Dicom Header\n");  
	DATA_INFO DI;

	/*TEMP*/
	strcpy(serieshead.prtcl,"");


	if( (flag->output_type != DAT_ONLY) || (flag->dat_to_dicom==1) ){

		initialize_dicom_header( flag, &DI, &rdbhead,&imagehead, &serieshead, &examhead);

		if(flag->dat_to_dicom==1){
			printf("External Dat to Dicom\n");
			flag->mag_out = 1;
			flag->cd_out = 0;
			image_mag = (mtx3d *)mtx3d_alloc(flag->rczres,flag->rcyres,flag->rcxres);
			fp = fopen(flag->dat_name,"r");
			fread(image_mag->m[0][0],sizeof(float),flag->rcxres*flag->rcyres*flag->rczres,fp);
			fclose(fp);

			max_signal = maxsig(image_mag);
			printf("Mag Max Signal %f\n",max_signal);
			scale_factor =( DICOM_MAX/max_signal ); 
			scale_mtx3d( image_mag, scale_factor);
			export_dicom(image_mag,image_comp,&DI,flag,&rdbhead,&imagehead,&serieshead,&examhead,0);
			exit(1);
		}
	}
#endif 

	/***********************************************************************
	  Allocate matrices required for all stages                                
	 ***********************************************************************/
	dataset=(cmtx *)cmtx_alloc(flag->xres,flag->nproj);
	if(flag->preload_data==1){
		flag->preloaded_data = ( cmtx **)malloc( sizeof(cmtx *)*flag->numrecv);
		for (recv_cnt=flag->startrecv; recv_cnt <= flag->endrecv; recv_cnt++){
			printf("Preload Data for Coil %d\n",recv_cnt);
			flag->preloaded_data[recv_cnt] = (cmtx *)cmtx_alloc(flag->xres,flag->nproj);
			read_comp(flag->preloaded_data[recv_cnt],recv_cnt,flag,0,flag->demod_freq,1);
			apply_offisocenter_data(flag->preloaded_data[recv_cnt],kmap,flag,flag->xshift,flag->yshift,flag->zshift,vd);
		}}

		if(  (isense_data->isense_flag==1) && (flag->gating!=GATE_NONE)  ){
			tornado_info->store_tornado_weight = 0;
			tornado_info->tornado_weight = (mtx *)mtx_alloc(flag->xres,flag->nproj);
		}
		image_comp= (mtx3d *)mtx3d_alloc(flag->rczres,flag->rcyres,flag->rcxres);
		image_mag = (mtx3d *)mtx3d_alloc(flag->rczres,flag->rcyres,flag->rcxres);
		grd_info->image_mag = image_mag;

		/*Calculates off-resonance map*/
		if(flag->correct_offres){
			FREQ = get_offres_map( grd_info,flag, dataset);
		}
		/*****************************************************************
		  Allocation moved here to reduce off-resonance memory use
		 ******************************************************************/

		printf("Matrix(og): %d(%f) x %d(%f) x %d(%f)\n",flag->rczres,grd_info->grid_z,flag->rcyres,grd_info->grid_y,flag->rcxres,grd_info->grid_x );
		printf("Og Starts: x:%d y:%d z:%d\n",grd_info->og_sx,grd_info->og_sy,grd_info->og_sz);
		printf("Og Stops: x:%d y:%d z:%d\n",grd_info->og_ex,grd_info->og_ey,grd_info->og_ez);

		/*Gridding Matrix + Place to Store*/
		k3d_grid =  (cmtx3d *)cmtx3d_alloc((int)(grd_info->grid_z * flag->rczres),(int)(grd_info->grid_y * flag->rcyres),(int)(grd_info->grid_x *flag->rcxres));
		image   = (cmtx3d *)cmtx3d_alloc(flag->rczres,flag->rcyres,flag->rcxres);	
		grd_info->k3d_grid = k3d_grid;
		grd_info->image   = image;

		if(flag->grad_warp){
			k3d_warp = (cmtx3d *)cmtx3d_alloc(flag->rczres,flag->rcyres,flag->rcxres);	
		}

		/**Allocate Iterative Sense Structures VERY MEMORY INTENSIVE!!!*/
		if(isense_data->isense_flag==1){
			allocate_isense_structs(isense_data,flag->rczres,flag->rcyres,flag->rcxres,flag->xres,flag->nproj,flag->numrecv);
		}

		/**Phase Contrast Arrays*/
		if(flag->recon_type == PHASE_DIFF){
			/**In iterative sense this is handled later*/
			k3dvel  = (cmtx3d *)cmtx3d_alloc(flag->rczres,flag->rcyres,flag->rcxres);
			k3d_s   = (cmtx3d *)cmtx3d_alloc(flag->rczres,flag->rcyres,flag->rcxres);

			int number_of_theta = ( flag->nvd > 3 ) ? ( flag->nvd-1) : ( 3 );
			theta = ( mtx3d_short **)malloc( sizeof(mtx3d_short *)*( number_of_theta));
			for(i=0; i< number_of_theta;i++){
				theta[i] = (mtx3d_short *)mtx3d_short_alloc(flag->rczres,flag->rcyres,flag->rcxres);
			}
		}

		/**For Time Resolved Mipping*/
		if(flag->time_mip == 1){
			time_mip = (mtx3d *)mtx3d_alloc(flag->rczres,flag->rcyres,flag->rcxres);
			reinitimagespace(time_mip);
		}

		/**Initialize Fourier transform*/
		init_recon_struct(&(k3d_grid->fftw_plan_forward), k3d_grid->xs, k3d_grid->ys, k3d_grid->zs,flag->threads, k3d_grid ,FFTW_FORWARD,flag->ss_2d);
		if(isense_data->isense_flag==1){
			init_recon_struct(&(k3d_grid->fftw_plan_backward), k3d_grid->xs, k3d_grid->ys, k3d_grid->zs,flag->threads, k3d_grid ,FFTW_BACKWARD,flag->ss_2d);
		}


		/*****************************************************************
		  Prep Things for Recon
		 ******************************************************************/

		/*The mag phase will always be zero*/
		printf("( Start Flow = %.3f)\n", gettime()-start_time );
		printf("********************************\n");
		printf("**     Standard Recon         **\n");
		printf("********************************\n");

		/*Setting up multi-frequency reconstruction*/
		setup_mfi_parameters( flag, FREQ);

		if(isense_data->isense_flag){

			reinit_isense_structs( isense_data);	

			printf("*****Getting Coil Sensitivity******\n");
			if ( flag->phantom==1){
				/**Collects data for digital phantom*/
				setup_phantom(grd_info,flag,isense_data,image,k3d_grid);
			}else{
				get_coil_sensitivity( flag, grd_info,dataset,isense_data,grd_info);
			}

			if(isense_data->foccus_flag==1){
				load_foccus_image( isense_data,flag->rczres,flag->rcyres,flag->rcxres);
			}

			printf("Getting E^H Scale Factor\n");
			reinit_isense_structs( isense_data);
			if( isense_data->cg_max_iter > 0){
				get_cg_dens(isense_data,flag);
				setup_isense_scale(grd_info,flag,isense_data,image,k3d_grid,kmap);
			}
		} 

		/****************************************************************
		 *** The iterative CG is getting complicated.
		 ***
		 *** Seperating into 2 Recons
		 ***   1) Standard Gridding w/wo MFI 
		 ***   2) Iterative CG Recon
		 *****************************************************************/
		for( time_frame = flag->fr_start; time_frame< flag->end_frame; time_frame++){ 

			printf("Recon Time Frame %d\n",time_frame);

			/********REINIT EVERYTHING*********/
			reinitimagespace(image_mag);
			reinitimagespace(image_comp);

			if(isense_data->isense_flag == 0){
				for (vdt = flag->start_vd; vdt <  flag->stop_vd; vdt++) {
					printf("\t( Start vd %d = %.3f)\n",vdt, gettime()-start_time );
					if(flag->external_data==1){
						read_external_kmap( flag,vdt);
					}
					if(flag->recon_type == PHASE_DIFF){
						reinitkspace(k3d_s);
					}

					if( flag->multi_echo == 2 && flag->recon_type == CSI){
						vd = ( vdt < flag->nvd ) ? ( vdt ) : ( vdt-flag->nvd);
						flag->echo_start = ( vdt < flag->nvd ) ? ( 0 ) : ( 1 );
						flag->echo_stop = ( vdt < flag->nvd ) ? ( 0 ) : ( 1 );
						printf("VDT %d VD %d\n",vdt,vd);	
					}else{
						vd = vdt;
					}

					for (recv_cnt=flag->startrecv; recv_cnt <= flag->endrecv; recv_cnt++){

						/***Keep track of times ****/
						grid_time = igrid_time = ifft_time = fft_time = read_time = combine_time = mfi_time =  deapp_time = maxwell_time =  tornado_time=0;

						/**************************************
						 *    Coil Skippping                  *
						 **************************************/
						if(flag->skipcoil[recv_cnt]){ 
							printf("Skipping coil %d \n",recv_cnt);
							continue;
						}

						/**************************************
						 ***   Reference Reconstruction        *
						 ***************************************/
						reinitkspace(image);	
						for(int dem=flag->start_demod; dem <= flag->stop_demod; dem++){

							recon_time += -gettime();
							/**************************************
							 ***   Setup for Demodulation          *
							 ***************************************/
							printf("\t\t( Dem=%d Coil=%d t= %.3f)\n",dem,recv_cnt,gettime()-start_time );
							demod_freq = ( flag->correct_offres == 1) ? ( flag->interp_freqs[dem] ) : ( flag->demod_freq );

							/************************************
							 ***  Base Recon                     *
							 *************************************/

							/***READ or DEGRID***/
							printdbg("\t\t\t Reading Data\n");
							vdgrid = (flag->recon_type == PHASE_DIFF) ? ( 0 ) : ( vd );

							/**just read data*/
							read_time+= -gettime();
							if(flag->preload_data==0){
								read_comp(dataset,recv_cnt,flag,vdgrid,demod_freq,1);
								apply_offisocenter_data(dataset,kmap,flag,flag->xshift,flag->yshift,flag->zshift,vd);
							}else{
								copy_cmtx( flag->preloaded_data[recv_cnt],dataset);
							}
							read_time+= gettime();

							tornado_time+= -gettime();
							printdbg("\t\t\t Tornado Weighting\n");
							/**Weighting of Data**/				
							if(time_frame >= 0){
								torweight_data(dataset,flag->proj_position,tornado_info,kspace_radius,flag->time_phase[time_frame],vdgrid,1,0);
							}
							tornado_time+= gettime();

							/***GRID to Cartesian***/
							printdbg("\t\t\t Gridding Data\n");
							grid_time+= -gettime();
							reinitkspace(k3d_grid);
							grid_data(grd_info,(cmtx *)dataset, (kspace_map *)kmap, 0, flag,isense_data,tornado_info,(cmtx3d*)k3d_grid,flag->xshift,flag->yshift,flag->zshift,DENS_NORMAL,vdgrid,recv_cnt,GRID_FORWARD);
							grid_time+= gettime();

							/***Write out data if needed*/	
							if( flag->rawon == 3){
								sprintf(fname,"KSPACE_c%d_v%d.dat",recv_cnt,vdt);
								export_cmtx3d(k3d_grid,fname);
							}

							printdbg("\t\t\t Doing FFT\n");
							fft_time+= -gettime();
							fft_3d(&(k3d_grid->fftw_plan_forward),k3d_grid,flag->threads);
							fft_time+= gettime();

							printdbg("\t\t\t Deappodization\n");
							deapp_time += -gettime();
							deapp3d_grid(k3d_grid,grd_info,flag->threads);
							deapp_time +=  gettime();

							/***Write out data if needed*/	
							if(flag->rawon==3){
								sprintf(fname,"OG_IMAGE.dat");
								export_cmtx3d(k3d_grid,fname);
							}

							printdbg("\t\t\t Doing MFI or CROP\n");
							int mfi_type = (flag->correct_offres==1) ? ( MFI ) : ( CROP );
							mfi_time += -gettime();
							multi_freq_interp( grd_info, flag, isense_data, k3d_grid, image, FREQ, demod_arg, dem,mfi_type,recv_cnt);
							mfi_time +=  gettime();
						}

						/**Done With Main Recon Correct Maxwell Phase*/
						if(flag->maxwell==1 && ( flag->recon_type == PHASE_DIFF  || flag->recon_type==CSI)){
							maxwell_time += -gettime();
							maxwell_correction(vdgrid,image,flag,&imagehead,&serieshead,&acq_tab);
							maxwell_time += gettime();
						}

						/***Write out data if needed*/	
						if( (flag->recon_type == CSI) || flag->rawon == 3){
							if( (mag_scale == -1.0) && (flag->scale_complex_images==1)){
								mag_scale = 4000.0/maxsig_cmtx3d( image);
								printf("Scaling Complex data by %f\n",mag_scale);
								scale_cmtx3d( image, mag_scale);
							}else{
								scale_cmtx3d( image, mag_scale);
							}			
							sprintf(fname,"IMAGE_c%d_v%d.dat",recv_cnt,vdt);
							export_cmtx3d(image,fname);
						}

						/************************************
						 ***  Flow Recon                     *
						 *************************************/
						if(flag->recon_type == PHASE_DIFF){
							reinitkspace(k3dvel);
							for(int dem=flag->start_demod; dem <= flag->stop_demod; dem++){
								recon_time += -gettime();

								/**************************************
								 ***   Setup for Demodulation          *
								 ***************************************/
								printf("\t\t( Dem=%d Coil=%d t= %.3f)\n",dem,recv_cnt,gettime()-start_time );
								demod_freq = ( flag->correct_offres == 1) ? ( flag->interp_freqs[dem] ) : ( flag->demod_freq );

								/*read data change for DualVenc*/						
								read_time += -gettime();
								read_comp(dataset,recv_cnt,flag, vd/*temp*/,demod_freq,1);
								apply_offisocenter_data(dataset,kmap,flag,flag->xshift,flag->yshift,flag->zshift,vd);
								vdgrid = vd;
								read_time += gettime();

								tornado_time += -gettime();
								if(time_frame >= 0){
									torweight_data(dataset,flag->proj_position,tornado_info,kspace_radius,flag->time_phase[time_frame],vdgrid,1,0);
								}
								tornado_time +=  gettime();


								reinitkspace(k3d_grid);
								grid_time += -gettime();
								grid_data(grd_info,(cmtx *)dataset, (kspace_map *)kmap, 0, flag, isense_data,tornado_info,(cmtx3d*)k3d_grid,flag->xshift,flag->yshift,flag->zshift,0 /*Offres flag*/,vdgrid,recv_cnt,GRID_FORWARD);
								grid_time += gettime();

								fft_time += -gettime();
								fft_3d(&(k3d_grid->fftw_plan_forward),k3d_grid,flag->threads);
								recon_time += gettime();
								fft_time += gettime();

								deapp_time += -gettime();
								deapp3d_grid(k3d_grid,grd_info,flag->threads);
								deapp_time +=  gettime();

								int mfi_type = ( flag->correct_offres == 1) ? ( MFI ) : ( CROP );
								mfi_time += -gettime();
								multi_freq_interp( grd_info,flag, isense_data, k3d_grid, k3dvel, FREQ, demod_arg, dem,mfi_type,recv_cnt);
								mfi_time +=  gettime();

							}		
							/***Maxwell phase correction*/
							if(flag->maxwell){
								maxwell_time += -gettime();
								maxwell_correction(vdgrid,k3dvel,flag,&imagehead,&serieshead,&acq_tab);
								maxwell_time += gettime();
							} 	
						}

						/**************************************
						 ***  Coil Combination                 *
						 ***************************************/
						printdbg("\t\t\t Coil Combination\n");
						combine_time = -gettime();
						coil_combine(flag,isense_data,image, k3dvel,k3d_s,image_mag,recv_cnt, FORWARD);
						combine_time += gettime();

						printdbg("\t\t( TIMES:\n");
						printdbg("\t\t\t read      %.3f \n",read_time);
						printdbg("\t\t\t tornado   %.3f \n",tornado_time);
						printdbg("\t\t\t grid      %.3f \n",grid_time);
						printdbg("\t\t\t fft       %.3f \n",fft_time);
						printdbg("\t\t\t igrid     %.3f \n",igrid_time);
						printdbg("\t\t\t ifft      %.3f \n",ifft_time);
						printdbg("\t\t\t combine   %.3f \n",combine_time);
						printdbg("\t\t\t mfi       %.3f \n",mfi_time );
						printdbg("\t\t\t coil sens %.3f \n",map_time);
						printdbg("\t\t\t deapp     %.3f \n",deapp_time);
						printdbg("\t\t\t maxwell   %.3f \n",maxwell_time);
						if(flag->verbose){ print_memusage();}

					} /* receiver loop ends */

					/******************************
					 **  Updates for Phase Contrast
					 **   DualVenc 
					 *******************************/
					if(flag->recon_type == PHASE_DIFF){

						if(flag->grad_warp){
							printf("(Grad Warp = %.3f)\n", gettime()-start_time ); 
							perform_grad_warp_complex(flag,k3d_warp,k3d_s,&imagehead,&serieshead,&acq_tab);
						}

						printf("\t\t\t( Pre-Phase = %.3f)\n", gettime()-start_time );
						image_phase( (cmtx3d*)k3d_s, (mtx3d_short *)theta[vd-1], flag);

						printf("\t\t\t( Post-Phase = %.3f)\n", gettime()-start_time );
					}
					/*End Phase difference updates*/

					printf("\t\t( End vd %d = %.3f)\n",vd, gettime()-start_time );

				}/*Velocity Direction done*/

			}else{

				/*********************************************************
				 ***  CG  Looping                                        **
				 **********************************************************/

				for (vdt = flag->start_vd; vdt <  flag->stop_vd; vdt++) {
					printf("\t( Start vd %d = %.3f)\n",vdt, gettime()-start_time );

					if(flag->external_data==1){
						read_external_kmap( flag,vdt);
					}

					if(flag->recon_type == PHASE_DIFF){
						reinitkspace(k3d_s);
					}

					if( flag->multi_echo == 2 && flag->recon_type == CSI){
						vd = ( vdt < flag->nvd ) ? ( vdt ) : ( vdt-flag->nvd);
						flag->echo_start = ( vdt < flag->nvd ) ? ( 0 ) : ( 1 );
						flag->echo_stop = ( vdt < flag->nvd ) ? ( 0 ) : ( 1 );
						printf("VDT %d VD %d\n",vdt,vd);	
					}else{
						vd = vdt;
					}

					int iteration_number = 1;
					reinit_isense_structs( isense_data);

					if(isense_data->L1_diff_flag==1){
						FILE *fid = fopen(isense_data->L1_diff_filename,"r");
						fread(isense_data->Image_L1_diff->m[0][0],sizeof(complex_u),flag->rcxres*flag->rcyres*flag->rczres,fid);
						fclose(fid);
					}
					iteration_number = isense_data->cg_max_iter+1;
					isense_data->max_kdata = 0.0;
					isense_data->scale_kdata = 0.0;
					isense_data->X0_scaled = 0;
					isense_data->initial_error=0.0; 

					/****************************************************************************
					 ****  Get Initial Estimate and Initialize
					 *****************************************************************************/ 

					/**Iterative Sense Looping - First Pass Stup RHS*/
					printf("Initialize Image\n");
					if(isense_data->initialize_type == INIT_ZERO){
						isense_data->scale_kdata=1.0;
						isense_data->X0_scaled = 1;
					}else if(isense_data->initialize_type == INIT_ONES){
						for(int k=0; k < isense_data->X0->zs; k++){
							for(int j=0; j < isense_data->X0->ys; j++){
								for(int i=0; i < isense_data->X0->xs; i++){
									isense_data->X0->m[k][j][i].i = 1.0;
									isense_data->X0->m[k][j][i].q = 0.0;
								}}}
					}else if(isense_data->initialize_type != INIT_EXTERNAL){
						tornado_info->store_tornado_weight = 0;
						for (recv_cnt=flag->startrecv; recv_cnt <= flag->endrecv; recv_cnt++){
							/**   Coil Skippping **/
							if(flag->skipcoil[recv_cnt]){ 
								printf("Skipping coil %d \n",recv_cnt);
								continue;
							}

							if(flag->preload_data==0){  
								read_comp(dataset,recv_cnt,flag,vd,flag->demod_freq,1);
								apply_offisocenter_data(dataset,kmap,flag,flag->xshift,flag->yshift,flag->zshift,vd);
							}else{
								copy_cmtx( flag->preloaded_data[recv_cnt],dataset);
							}

							if( (time_frame >= 0) && (isense_data->initialize_type != INIT_COMPOSITE) ){
								float old_wdth_low = 0;
								float old_wdth_high = 0;
								if( isense_data->initialize_type == INIT_SLIDING_WINDOW){
									old_wdth_low = tornado_info->wdth_low;
									old_wdth_high= tornado_info->wdth_high;
									tornado_info->wdth_low  *= isense_data->sliding_window_frac;
									tornado_info->wdth_high *= isense_data->sliding_window_frac;
								}
								torweight_data(dataset,flag->proj_position,tornado_info,kspace_radius,flag->time_phase[time_frame],vd,1,0);
								if( isense_data->initialize_type == INIT_SLIDING_WINDOW){
									tornado_info->wdth_low  = old_wdth_low;
									tornado_info->wdth_high = old_wdth_high;
								}
							}

							reinitkspace(k3d_grid);
							grid_data(grd_info,(cmtx *)dataset, (kspace_map *)kmap, 0, flag, isense_data,tornado_info,(cmtx3d*)k3d_grid,flag->xshift,flag->yshift,flag->zshift,DENS_NORMAL,vd,recv_cnt,GRID_FORWARD);
							fermi_3d( k3d_grid, flag->fermi_r, flag->fermi_w,flag->threads,flag->zoom_x/grd_info->grid_x,flag->zoom_y/grd_info->grid_y,flag->zoom_z/grd_info->grid_z);
							fft_3d(&(k3d_grid->fftw_plan_forward),k3d_grid,flag->threads);
							deapp3d_grid(k3d_grid,grd_info,flag->threads);
							multi_freq_interp(grd_info, flag, isense_data,k3d_grid, image, FREQ, demod_arg,0,CROP,recv_cnt);
							coil_combine(flag,isense_data,image, k3dvel,isense_data->X0,image_mag,recv_cnt, ISENSE);
						} /* receiver loop ends */
					}else{
						printf("Reading External Image\n");
						fp=fopen(isense_data->init_external_name,"r");
						fread(isense_data->X0->m[0][0],sizeof(complex_u),isense_data->X0->xs*isense_data->X0->ys*isense_data->X0->zs,fp);
						fclose(fp);
						isense_data->scale_kdata=1.0;
						isense_data->X0_scaled = 1;
					}

					if(flag->rawon == 11){
						export_isense_structs(isense_data);
					}	

					/***************************************************************
					 *** Data Weighting for CG Weighted Least Squares
					 ***************************************************************/ 		
					if( isense_data->cg_max_iter > 0){

						printf("******Getting CG Dens****************\n");
						get_cg_dens(isense_data,flag);

						tornado_info->store_tornado_weight=1;
						tornado_info->tornado_weight_exists=0;
						if(time_frame >= 0){
							torweight_data(dataset,flag->proj_position,tornado_info,kspace_radius,flag->time_phase[time_frame],vd,0 /*No density comp*/,1);
							scale_mtx_by_mtx( isense_data->cg_dens, tornado_info->tornado_weight,0);
							/* export_mtx(isense_data->cg_dens,"CG_DENS.dat");*/
						}
						tornado_time+= gettime();

					}

					/*****************************************************************************
					 *** Now iterative recon	 
					 ******************************************************************************/ 

					int irls_iterations = ( (isense_data->L1_image_flag==1) || (isense_data->L1_wave_flag==1) || (isense_data->L1_diff_flag==1) || (isense_data->L1_phase_flag==1)|| (isense_data->L1_tv_flag==1)) ? ( isense_data->irls_max_iter ) : ( 1 );

					printf("Iterative Solver Image\n");

					/**Iterative Sense Looping - First Pass Setup RHS*/
					for(int irls_iter =0; irls_iter < irls_iterations; irls_iter++){
						isense_L1_update_reweighting(isense_data,0);

						for(int isense_iter= 0; isense_iter< isense_data->cg_max_iter; isense_iter++){
							double iteration_total = -gettime();

							if(isense_iter%isense_data->isense_reinit ==0){
								reinitkspace(isense_data->P0);
								reinitkspace(isense_data->R0);
							}else{
								reinitkspace(isense_data->LHS);
							}

							/***Keep track of times ****/
							grid_time = 0;
							igrid_time = 0;
							ifft_time = 0;
							fft_time = 0;
							read_time = 0;
							combine_time = 0;
							mfi_time = 0;
							deapp_time =0;
							maxwell_time = 0;
							tornado_time=0;
							L1_regularize_time = 0;
							L2_regularize_time = 0;
							double CG_time = 0;

							for (recv_cnt=flag->startrecv; recv_cnt <= flag->endrecv; recv_cnt++){


								/**    Coil Skippping      **/
								if(flag->skipcoil[recv_cnt]){ 
									printf("Skipping coil %d \n",recv_cnt);
									continue;
								}

								recon_time += -gettime();

								/**Forward Model*/
								mfi_time += -gettime();
								reinitkspace(k3d_grid);
								if(isense_iter%isense_data->isense_reinit==0){
									multi_freq_interp( grd_info,flag, isense_data,k3d_grid, isense_data->X0, FREQ, demod_arg,0,ISENSE,recv_cnt);
								}else{
									multi_freq_interp( grd_info,flag, isense_data,k3d_grid, isense_data->P0, FREQ, demod_arg,0,ISENSE,recv_cnt);
								}
								mfi_time += gettime();

								deapp_time += -gettime();
								deapp3d_grid(k3d_grid,grd_info,flag->threads);
								deapp_time +=  gettime();

								ifft_time += -gettime();
								fft_3d_backwards(&(k3d_grid->fftw_plan_backward),k3d_grid,flag->threads);
								ifft_time +=  gettime();

								igrid_time+= -gettime();
								grid_data(grd_info,isense_data->dataset,kmap, 0, flag, isense_data,tornado_info,k3d_grid,0,0,0,0 /*Offres flag*/,vd,recv_cnt,GRID_BACKWARD);
								scale_cmtx_by_mtx( isense_data->dataset,isense_data->cg_dens,flag->threads);
								igrid_time+=  gettime();

								/*Residue Calculation*/
								read_time+= -gettime();
								if(isense_iter%isense_data->isense_reinit == 0){
									if(flag->preload_data==0){
										read_comp(dataset,recv_cnt,flag,vd,flag->demod_freq,1);
										apply_offisocenter_data(dataset,kmap,flag,flag->xshift,flag->yshift,flag->zshift,vd);
									}else{
										copy_cmtx( flag->preloaded_data[recv_cnt],dataset);
									}
									scale_cmtx_by_mtx( dataset,isense_data->cg_dens,flag->threads);

									if(isense_data->scale_kdata==0.0){
										isense_data->scale_kdata = energy_cmtx(dataset) / energy_cmtx(isense_data->dataset);
										printf("Kdata Scale %f\n",isense_data->scale_kdata);
										scale_cmtx3d(isense_data->X0,isense_data->scale_kdata);
										scale_cmtx(isense_data->dataset,isense_data->scale_kdata);
										isense_data->X0_scaled = 1;
									}
									// float pre_energy = energy_cmtx(dataset);
									isense_subtract_kdata( dataset, isense_data->dataset, isense_data);
									// printf("Energy Difference %f percent (%e vs %e)\n",100.0*energy_cmtx(dataset)/energy_cmtx(isense_data->dataset),pre_energy,energy_cmtx(isense_data->dataset));

								}else{
									/*Swap Pointers */
									cmtx *temp_dataset;
									temp_dataset = dataset;
									dataset = isense_data->dataset;
									isense_data->dataset = temp_dataset;
								}
								read_time+= gettime();


								/***GRID to Cartesian***/
								grid_time+= -gettime();
								reinitkspace(k3d_grid);
								grid_data(grd_info,(cmtx *)dataset, (kspace_map *)kmap, 0, flag,isense_data,tornado_info,(cmtx3d*)k3d_grid,0,0,0,DENS_NONE,vd,recv_cnt,GRID_FORWARD);
								fermi_3d( k3d_grid, flag->fermi_r, flag->fermi_w,flag->threads,flag->zoom_x/grd_info->grid_x,flag->zoom_y/grd_info->grid_y,flag->zoom_z/grd_info->grid_z);
								grid_time+= gettime();

								fft_time+= -gettime();
								fft_3d(&(k3d_grid->fftw_plan_forward),k3d_grid,flag->threads);
								fft_time+= gettime();

								/* No need to deappodize twice*/	
								deapp_time += -gettime();
								deapp3d_grid(k3d_grid,grd_info,flag->threads);
								deapp_time +=  gettime();

								mfi_time += -gettime();
								multi_freq_interp( grd_info,flag,isense_data,k3d_grid, image, FREQ, demod_arg,0,CROP,recv_cnt);
								mfi_time +=  gettime();


								/**************************************
								 ***  Coil Combination                 *
								 ***************************************/
								if(isense_iter%isense_data->isense_reinit ==0){
									/**Collect RHS E'd */
									combine_time = -gettime();
									coil_combine(flag, isense_data,image, k3dvel,isense_data->R0,image_mag,recv_cnt, ISENSE);
									combine_time += gettime();
								}else{
									combine_time = -gettime();
									coil_combine(flag, isense_data,image, k3dvel,isense_data->LHS,image_mag,recv_cnt, ISENSE);
									combine_time += gettime();
								}
							} /* receiver loop ends */


							if(isense_iter%isense_data->isense_reinit ==0){
								L2_regularize_time = -gettime();
								/**Collect Residue - Initialize */
								isense_L2_regularization( isense_data,1);
								L2_regularize_time += gettime();

								L1_regularize_time = -gettime();
								isense_L1_regularization( isense_data,1);
								L1_regularize_time += gettime();

								copy_cmtx3d( isense_data->R0, isense_data->P0);
								isense_data->scale_kdata = 1.0;
								if(isense_data->initial_error==0){
									isense_data->initial_error = energy_cmtx3d( isense_data->R0);
								}
							}else{
								L2_regularize_time = -gettime();
								isense_L2_regularization( isense_data,0);
								L2_regularize_time += gettime();

								L1_regularize_time = -gettime();
								isense_L1_regularization( isense_data,0);
								L1_regularize_time += gettime();

								CG_time = -gettime();
								isense_cg_update( isense_data);
								CG_time +=  gettime();
								isense_data->current_error = energy_cmtx3d( isense_data->R0);
								printf("\tIteration %d :: Current Error = %e\n",isense_iter,isense_data->current_error/isense_data->initial_error);
							}

							printdbg("\t\t( TIMES:\n");
							printdbg("\t\t\t read      %.3f \n",read_time);
							printdbg("\t\t\t tornado   %.3f \n",tornado_time);
							printdbg("\t\t\t grid      %.3f \n",grid_time);
							printdbg("\t\t\t fft       %.3f \n",fft_time);
							printdbg("\t\t\t igrid     %.3f \n",igrid_time);
							printdbg("\t\t\t ifft      %.3f \n",ifft_time);
							printdbg("\t\t\t combine   %.3f \n",combine_time);
							printdbg("\t\t\t mfi       %.3f \n",mfi_time );
							printdbg("\t\t\t coil sens %.3f \n",map_time);
							printdbg("\t\t\t deapp     %.3f \n",deapp_time);
							printdbg("\t\t\t maxwell   %.3f \n",maxwell_time);
							printdbg("\t\t\t L1 Reg    %.3f \n",L1_regularize_time);
							printdbg("\t\t\t L2 Reg    %.3f \n",L2_regularize_time);
							printdbg("\t\t\t CG        %.3f \n",CG_time);
							printdbg("\t\t\t Total     %.3f \n",iteration_total+gettime());


							// if(flag->verbose){ print_memusage();}
							if(flag->rawon == 11){
								export_isense_structs(isense_data);
							}
						}/*ISENSE ITER*/
						isense_data->X0 =complex_wavelet_thresh( isense_data->X0,0,0,0);



					}/*IRLS ITER*/

					if(isense_data->normalize_coils ==0){
						renormalize_image( isense_data->X0,image_mag,isense_data,flag->startrecv,flag->endrecv);
					}

					if(isense_data->foccus_flag==1){
						scale_cmtx3d_by_cmtx3d( isense_data->X0,isense_data->constraining_image,flag->threads);	
					}


					/******************************
					 **  Updates for Phase Contrast
					 **   DualVenc 
					 *******************************/

					/*GradWarp*/	
					if(flag->grad_warp){
						printf("(Grad Warp = %.3f)\n", gettime()-start_time ); 
						perform_grad_warp_complex(flag,k3d_warp,isense_data->X0,&imagehead,&serieshead,&acq_tab);
					}

					if( flag->export_complex_images){
						if(time_frame==-1){
							sprintf(fname,"Complex_Image_Comp_V_%03d.dat",vd);
						}else{
							sprintf(fname,"Complex_Image_Frame_%03d_V_%03d.dat",time_frame,vd);
						}
						export_cmtx3d(isense_data->X0,fname);
					}


					if(flag->recon_type == PHASE_DIFF){

						if(flag->maxwell){
							maxwell_correction(vd,isense_data->X0,flag,&imagehead,&serieshead,&acq_tab);
						} 	

						/**Special Updates for Iterative SENSE*/
						if(vd ==0 ){
							copy_cmtx3d(isense_data->X0,k3dvel);
							continue;
						}else{
							copy_cmtx3d(isense_data->X0,image);
							coil_combine(flag, isense_data,k3dvel, image,k3d_s,image_mag,recv_cnt, FORWARD);
						}

						printf("\t\t\t( Pre-Phase = %.3f)\n", gettime()-start_time );
						image_phase( (cmtx3d*)k3d_s, (mtx3d_short *)theta[vd-1], flag);

						printf("\t\t\t( Post-Phase = %.3f)\n", gettime()-start_time );
						if(flag->grad_warp ){	
							k3d_s = k3d_warp; /* FFTW requires static address*/
						}

					}else if(flag->recon_type== CSI){
						/***Write out data if needed*/	
						if(flag->maxwell){
							maxwell_correction(vd,isense_data->X0,flag,&imagehead,&serieshead,&acq_tab);
						} 	

						if(mag_scale == -1.0){
							mag_scale = 4000.0/maxsig_cmtx3d( isense_data->X0);
							printf("Scaling Complex data by %f\n",mag_scale);
						}			
						scale_cmtx3d(  isense_data->X0, mag_scale);
						sprintf(fname,"IMAGE_c%d_v%d.dat",0,vdt);
						export_cmtx3d( isense_data->X0,fname);

					}else{
						isense_store_solution( isense_data->X0, image_mag);
					}
					/*End Phase difference updates*/

					printf("\t\t( End vd %d = %.3f)\n",vd, gettime()-start_time );
				}/*Velocity Direction done*/

			}

			/***Correct Square of Signal*/
			coil_combine(flag, isense_data,image, k3dvel,k3d_s,image_mag,-1,BACKWARD);


			/******GRAD WARP TEMP*****/
			if(flag->grad_warp==1 && isense_data->isense_flag==0 ){
				printf("(Before Scale = %.3f)\n", gettime()-start_time ); 
				image_warp = (mtx3d *)mtx3d_alloc(flag->rczres,flag->rcyres,flag->rcxres);	
				perform_grad_warp(flag,image_warp,image_mag,&imagehead,&serieshead,&acq_tab);
			}

			/***************************************************
			 *** Raw Phase Output                               *
			 ****************************************************/	
			if(flag->rawon==1 && flag->recon_type == PHASE_DIFF){  
				/*Write Data*/
				printf("(Before Write = %.3f)\n", gettime()-start_time );

				if(time_frame== -1){
					for(vd=0; vd<flag->nvd-1;vd++){
						sprintf(fname,"theta_%d.dat",vd);
						export_mtx3d_short( theta[vd],fname);
					}
				}else{
					for(vd=0; vd<flag->nvd-1;vd++){
						sprintf(fname,"%stheta_%03d_vd_%d.dat",flag->out_folder,time_frame,vd);
						export_mtx3d_short( theta[vd],fname);
					}
				}
				printf("(After Write = %.3f)\n", gettime()-start_time );
				printf("\t( End Phase %d = %.3f)\n",n_ph, gettime()-start_time );

			}/***Vel mats****/

			/***************************************************
			 *** Velocity Calculations                          *
			 ****************************************************/


			/*Do Balanced Recon - DualVenc */
			printf("(Before Vel Det = %.3f)\n", gettime()-start_time );
			if(flag->recon_type == PHASE_DIFF){

				if(flag->composite_antialiasing==1 && time_frame==-1){
					int number_of_theta = ( flag->nvd > 3 ) ? ( flag->nvd-1) : ( 3 );
					theta_comp = ( mtx3d_short **)malloc( sizeof(mtx3d_short *)*( number_of_theta));
					for(i=0; i< number_of_theta;i++){
						theta_comp[i] = (mtx3d_short *)mtx3d_short_alloc(flag->rczres,flag->rcyres,flag->rcxres);
						memcpy( theta_comp[i]->m[0][0],theta[i]->m[0][0],(size_t)( sizeof(short)*flag->rcxres*flag->rcyres*flag->rczres));
					}
				}

				if(flag->grad_warp){
					generalized_pc(theta,theta_comp,flag); 
				}else{
					balanced_recon(theta,theta_comp,time_frame,flag); 
				}
			}

			printf("(After Vel Det  = %.3f)\n", gettime()-start_time );

			//if(flag->image_intensity_correction){
			//	image_intensity_correction( image_mag, flag->image_intensity_correction_order, flag->image_intensity_correction_,flag->image_intensity_correction_external);
			//}

			max_signal = maxsig(image_comp);
			if( (flag->recon_type == PHASE_DIFF) || (flag->output_type != DAT_ONLY) ){ 
				max_signal = maxsig(image_mag);
				printf("Mag Max Signal %f\n",max_signal);
				scale_factor = (flag->output_type == DAT_ONLY ) ? ( DAT_MAX/max_signal ) : ( DICOM_MAX/max_signal ); 
				scale_mtx3d( image_mag, scale_factor);
			}else if( flag->output_scale != -1){
				if(flag->output_scale_factor == -1){ 
					flag->output_scale_factor = flag->output_scale / maxsig(image_mag);
					printf("Mag Max Signal %f\n",max_signal);
				}	
				scale_mtx3d( image_mag, flag->output_scale_factor);
			}


			/***************************************************
			 *** Eddy current Phase Correction                  *
			 ****************************************************/
			if(flag->recon_type == PHASE_DIFF){    
				float cx = (float)flag->rcxres / 2.0;
				float cy = (float)flag->rcyres / 2.0;
				float cz = (float)flag->rczres / 2.0;

				if( (( flag->phase_correct==1) && (time_frame == -1)) || ( (flag->phase_correct==1) && (flag->phase_fit_dim != 3))  ){

					if(flag->rawon){  
						for(vd=0; vd<flag->nvd-1;vd++){
							sprintf(fname,"theta_%d_pre.dat",vd);
							export_mtx3d_short( theta[vd],fname);
						}	
					}/***Vel mats****/


					/***CD double check*/
					max_signal = maxsig(image_mag);
					min_mag= 0.15 * max_signal;
					max_cd = 0.30 * max_signal;

					for (k=0; k< flag->rczres; k++){
						for (j=0; j< flag->rcyres; j++){ 
							for (i=0; i< flag->rcxres; i++){
								mag = image_mag->m[k][j][i];

								vx=(float)theta[0]->m[k][j][i] * 2.0/flag->vencx;
								if(flag->flow_dir > 0){
									vy = 0;
									vz = 0;
								}else{
									vy=(float)theta[1]->m[k][j][i] * 2.0/flag->vency;
									vz=(float)theta[2]->m[k][j][i] * 2.0/flag->vencz;
								}
								vel = sqrtf( vx*vx + vy*vy + vz*vz);

								if(vel > 1.0){
									cd = fabs( image_mag->m[k][j][i]);
								}else{
									cd = fabs( image_mag->m[k][j][i]*sinf( PI/2 * vel));
								}
								if( (mag > min_mag) && (cd < max_cd) ){
									image_comp->m[k][j][i]=1.0;
								}else{
									image_comp->m[k][j][i]=0.0; 
								}}}}

								sprintf(fname,"segment_image.dat");
								export_mtx3d(image_comp,fname);

								if(flag->phase_fit_dim==3){	

									polyfit_vx = poly_fitting3d(image_comp, theta[0],flag->fit_order);
									polyfit_vy = poly_fitting3d(image_comp, theta[1],flag->fit_order);
									polyfit_vz = poly_fitting3d(image_comp, theta[2],flag->fit_order);

									poly_subtract3d(theta[0],polyfit_vx);
									poly_subtract3d(theta[1],polyfit_vy);
									poly_subtract3d(theta[2],polyfit_vz);
								}else{
									mtx *image_comp_slice; 
									mtx *theta1_slice;
									mtx *theta2_slice;
									mtx *theta3_slice;

									image_comp_slice=(mtx *)mtx_alloc(flag->rcyres,flag->rcxres);
									theta1_slice=(mtx *)mtx_alloc(flag->rcyres,flag->rcxres);
									theta2_slice=(mtx *)mtx_alloc(flag->rcyres,flag->rcxres);
									theta3_slice=(mtx *)mtx_alloc(flag->rcyres,flag->rcxres);


									for (k=0; k< flag->rczres; k++){
										printf("Fitting Slice %d\n",k);
										for (j=0; j< flag->rcyres; j++){ 
											for (i=0; i< flag->rcxres; i++){
												image_comp_slice->m[j][i]= image_comp->m[k][j][i];
												theta1_slice->m[j][i] = (float)theta[0]->m[k][j][i];
												theta2_slice->m[j][i] = (float)theta[1]->m[k][j][i];
												theta3_slice->m[j][i] = (float)theta[2]->m[k][j][i];
											}}	

										poly_fitting2d(image_comp_slice,theta1_slice,polyfit_vx);
										poly_fitting2d(image_comp_slice,theta1_slice,polyfit_vy);
										poly_fitting2d(image_comp_slice,theta1_slice,polyfit_vz);

										for (j=0; j< flag->rcyres; j++){ 
											for (i=0; i< flag->rcxres; i++){
												image_comp_slice->m[j][i]=image_comp->m[k][j][i];
												theta[0]->m[k][j][i] -= (short)poly2d(polyfit_vx,(float)i - cx, (float)j - cy);
												theta[1]->m[k][j][i] -= (short)poly2d(polyfit_vy,(float)i - cx, (float)j - cy);
												theta[2]->m[k][j][i] -= (short)poly2d(polyfit_vz,(float)i - cx, (float)j - cy);
											}}	

									}
								}

				}else if(flag->phase_correct ){

					/***Do Correction*/
					for (k=0; k< flag->rczres; k++){
						for (j=0; j< flag->rcyres; j++){ 
							for (i=0; i< flag->rcxres; i++){
								theta[0]->m[k][j][i] -= (short)poly3d(polyfit_vx,(float)i - cx, (float)j - cy, (float)k - cz);
								theta[1]->m[k][j][i] -= (short)poly3d(polyfit_vy,(float)i - cx, (float)j - cy, (float)k - cz);
								theta[2]->m[k][j][i] -= (short)poly3d(polyfit_vz,(float)i - cx, (float)j - cy, (float)k - cz);
							}}}
				}/*Phase correct*/
			}		

			/***************************************************
			 *** Phase Anti-Aliasing /Velocity Calcs              *
			 ****************************************************/

			if(flag->recon_type == PHASE_DIFF){	
				/****Common Phase code***/
				if(flag->anti_aliasing){
					printf("(Before 26 = %.3f)\n", gettime()-start_time );
					phase_comparison26(theta[0], flag->filtersteps_large, flag->rcxres, flag->rczres, flag, image_comp, max_signal,flag->vencx); 
					phase_comparison26(theta[1], flag->filtersteps_large, flag->rcxres, flag->rczres, flag, image_comp, max_signal,flag->vency); 
					phase_comparison26(theta[2], flag->filtersteps_large, flag->rcxres, flag->rczres, flag, image_comp, max_signal,flag->vencz); 

					printf("(After  26 = %.3f)\n", gettime()-start_time );

					printf("(Before 6 = %.3f)\n", gettime()-start_time );
					phase_comparison6(theta[0], flag->filtersteps_small, flag->rcxres,flag->rczres, flag, image_comp, max_signal,flag->vencx);
					phase_comparison6(theta[1], flag->filtersteps_small, flag->rcxres, flag->rczres, flag, image_comp, max_signal,flag->vency);
					phase_comparison6(theta[2], flag->filtersteps_small, flag->rcxres, flag->rczres, flag, image_comp, max_signal,flag->vencz);
					printf("(After  6 = %.3f)\n", gettime()-start_time );
				} 
			}

			/************************************************
			 **  Image Calculations                         **
			 *************************************************/


			if(flag->recon_type == PHASE_DIFF){
				printf("(Make CD Start = %.3f)\n", gettime()-start_time );
				make_cd( theta,image_mag,image_comp,flag);
				printf("(Make CD Done = %.3f)\n", gettime()-start_time );
			}

			if(flag->time_mip==1 && time_frame >= 0 ){

				for (k=0; k< flag->rczres; k++){
					for (j=0; j< flag->rcyres; j++){ 
						for (i=0; i< flag->rcxres; i++){
							time_mip->m[k][j][i] += image_comp->m[k][j][i]*image_comp->m[k][j][i];
							/*( time_mip->m[k][j][i] > image_comp->m[k][j][i] ) ? ( time_mip->m[k][j][i] ) : ( image_comp->m[k][j][i] );*/
						}}}

			}/*time mip*/



			if(flag->output_type != DICOM_ONLY){
				printf("Output Folder %s\n",flag->out_folder);
				if(flag->recon_type == PHASE_DIFF && flag->ss_2d==1){
					sprintf(fname_vx,"%scard_vd_1.dat",flag->out_folder);
					sprintf(fname_vy,"%scard_vd_2.dat",flag->out_folder);
					sprintf(fname_vz,"%scard_vd_3.dat",flag->out_folder);
					export_dat_2D_short(theta[0],fname_vx);
					export_dat_2D_short(theta[1],fname_vy);
					export_dat_2D_short(theta[2],fname_vz);
				}else if(flag->recon_type == PHASE_DIFF){
					if(time_frame == -1){
						sprintf(fname_vx,"%scomp_vd_1.dat",flag->out_folder);
						sprintf(fname_vy,"%scomp_vd_2.dat",flag->out_folder);
						sprintf(fname_vz,"%scomp_vd_3.dat",flag->out_folder);
					}else{
						sprintf(fname_vx,"%sph_%03d_vd_1.dat",flag->out_folder,time_frame);
						sprintf(fname_vy,"%sph_%03d_vd_2.dat",flag->out_folder,time_frame);
						sprintf(fname_vz,"%sph_%03d_vd_3.dat",flag->out_folder,time_frame);
					}
					export_dat_3D_short( theta[0],fname_vx);
					export_dat_3D_short( theta[1],fname_vy);
					export_dat_3D_short( theta[2],fname_vz);
				}

				if(flag->recon_type != PHASE_DIFF){
					if(flag->ss_2d){
						sprintf(fname,"%scard_mag.dat",flag->out_folder);
						export_dat_2D_float(image_mag,fname);
					}else{
						if(time_frame == -1){
							sprintf(fname,"%sMAG.dat",flag->out_folder);
						}else{
							sprintf(fname,"%sph_%03d_mag.dat",flag->out_folder,time_frame);
						}
						printf("%s\n",fname);
						export_dat_3D_float(image_mag,fname);
					}
				}else{
					convert_mtx3d_to_mtx3d_short( image_comp,theta[2]);
					convert_mtx3d_to_mtx3d_short( image_mag,theta[1]);

					if(flag->ss_2d){
						sprintf(fname,"%scard_cd.dat",flag->out_folder);
						export_dat_2D_short(theta[2],fname);
						sprintf(fname,"%scard_mag.dat",flag->out_folder);
						export_dat_2D_short(theta[1],fname);

					}else{
						if(time_frame == -1){
							sprintf(fname,"%sCD.dat",flag->out_folder);
						}else{
							sprintf(fname,"%sph_%03d_cd.dat",flag->out_folder,time_frame);
						}
						export_dat_3D_short(theta[2],fname);

						if(time_frame == -1){
							sprintf(fname,"%sMAG.dat",flag->out_folder);
						}else{
							sprintf(fname,"%sph_%03d_mag.dat",flag->out_folder,time_frame);
						}
						export_dat_3D_short(theta[1],fname);
					}
				}

			}

#ifdef UW_DICOM
			if( flag->output_type == DAT_AND_DICOM ){
				if( time_frame == -1){
					export_dicom(image_mag,image_comp,&DI,flag,&rdbhead,&imagehead,&serieshead,&examhead,time_frame);
				}
			}else if( flag->output_type == DICOM_ONLY){
				export_dicom(image_mag,image_comp,&DI,flag,&rdbhead,&imagehead,&serieshead,&examhead,time_frame);
			}
#endif

		}/***End cardiac time***/


		/***Export DICOM**/


		/******************TIME MIP OUTPUT****************************/

		if(flag->time_mip==1 && flag->gating!=GATE_NONE){

			for (k=0; k< flag->rczres; k++){
				for (j=0; j< flag->rcyres; j++){ 
					for (i=0; i< flag->rcxres; i++){
						time_mip->m[k][j][i] = sqrt(time_mip->m[k][j][i]);
					}}}

			if( flag->fout==0){
				max_signal = maxsig(time_mip);
				printf("Mag Max Signal %f\n",max_signal);
				scale_factor = (flag->output_type ==DAT_ONLY) ? ( DAT_MAX/max_signal ) : ( DICOM_MAX/max_signal ); 

				for (k=0; k< flag->rczres; k++){
					for (j=0; j< flag->rcyres; j++){ 
						for (i=0; i< flag->rcxres; i++){
							time_mip->m[k][j][i] *= scale_factor;
						}}}
			}



			if(flag->output_type != DICOM_ONLY){
				printf("(Before Write = %.3f)\n", gettime()-start_time );
				sprintf(fname,"%stime_mip.dat",flag->out_folder);
				fp=fopen(fname,"w");
				fwrite( time_mip->m[0][0],flag->rcxres*flag->rcyres*flag->rczres, sizeof(float), fp);
				fclose(fp);
			}

#ifdef UW_DICOM
			if(flag->output_type != DAT_ONLY){
				sprintf(DI.series_description,"SRC:PCVIPR TIME MIP");
				imagehead.im_seno    = flag->lx + 49;
				DI.new_series_number = flag->lx + 49;
				DI.base_series_number= flag->lx + 49;
				short *signa_image;
				for (k = 0; k < flag->rczres; k++){
					DI.image_number = k+1;
					strcpy(DI.outfilename, "");

					if (procGEImage(time_mip->m[k][0], &signa_image,&DI,&rdbhead,&examhead, 
								&serieshead, &imagehead) != 0) {
						printf("Error processing magnitude image data.\n");
						exit(1);
					}
					sprintf(output_filename, "%s%s",flag->out_folder,DI.outfilename);
					if(writeDCM(output_filename, rdbhead, examhead, serieshead,
								imagehead, signa_image, 0, 0, 0) != 0) {
						printf("Error writing magnitude image file.\n");
						exit (1);
					}
				}
			}
#endif
		}/*time mip*/


		printf("(End time = %.3f)\n", gettime()-start_time );

}

#endif

