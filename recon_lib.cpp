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
	
	external_dcf = false;
	recalc_dcf =false;
	dcf_iter = 20;
	dcf_dwin = 5.5;
	dcf_scale = 1.0;
	
	smap_res=8;
	intensity_correction = false;
	iterative_smaps = false;
		
	acc = 1;
	compress_coils = 0.0;
	whiten = false;
	export_smaps = 0;
	max_iter = 50;
	
	cycle_spins = 4;

	walsh_block_sizeX = 8;
	walsh_block_sizeY = 8;
	walsh_block_sizeZ = 8;
	
	// Gaussian Blur of smaps
	extra_blurX = 0.0;
	extra_blurY = 0.0;
	extra_blurZ = 0.0;
	
	prep_done = false;
	
	wavelet_levelsX=4; 
	wavelet_levelsY=4; 
	wavelet_levelsZ=4; 
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
	help_flag("-walsh","Use eigen vector aproach to estimate coil sensitivities");
	help_flag("-coil_lowres","default,use low resolution images to get coil sensitivies");
	help_flag("-export_smaps","write sensitivity maps");
	help_flag("-intensity_correction","polynomial fit for intensity correction");
	
	cout << "Options  Before Recon" << endl;
	help_flag("-recalc_dcf","Use iterative calc for density");
	help_flag("-dcf_iter","Iterations for DCF");
	help_flag("-dcf_dwin []","Size of kernel (default 5.5 -radius)");
	help_flag("-dcf_scale []","Scaling of matrix");
	help_flag("-external_dcf","Use (precomputed) DCF from external file");
	help_flag("-dcf_file []","Filename for external DCF usage");
	
	gridFFT::help_message();
	SPIRIT::help_message();	
	THRESHOLD::help_message();
	PHANTOM::help_message();
	GATING::help_message();
	L2REG::help_message();
	LOWRANKCOIL::help_message();
	
	
	
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
		
		trig_flag(true,"-recalc_dcf",recalc_dcf);
		int_flag("-dcf_iter",dcf_iter);
		float_flag("-dcf_dwin",dcf_dwin);
		float_flag("-dcf_scale",dcf_scale);
		
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
		int_flag("-wavelet_levelsX",wavelet_levelsX);
		int_flag("-wavelet_levelsY",wavelet_levelsY);
		int_flag("-wavelet_levelsZ",wavelet_levelsZ);
						
		// Coil Combination		
		trig_flag(ESPIRIT,"-espirit",coil_combine_type);
		trig_flag(WALSH,"-walsh",coil_combine_type);
		trig_flag(LOWRES,"-coil_lowres",coil_combine_type);
		float_flag("-smap_res",smap_res);
		trig_flag(1,"-export_smaps",export_smaps);
		trig_flag(true,"-intensity_correction",intensity_correction);
		trig_flag(true,"-iterative_smaps",iterative_smaps);
		int_flag("-walsh_block_sizeX",walsh_block_sizeX);
		int_flag("-walsh_block_sizeY",walsh_block_sizeY);
		int_flag("-walsh_block_sizeZ",walsh_block_sizeZ);
		float_flag("-extra_blurX",extra_blurX);
		float_flag("-extra_blurY",extra_blurY);
		float_flag("-extra_blurZ",extra_blurZ);
						
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
		int_flag("-cycle_spins",cycle_spins);

		trig_flag(true,"-external_dcf",external_dcf);
		char_flag("-dcf_file",dcffilename);
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
	gridding.precalc_gridding(rczres,rcyres,rcxres,data.trajectory_dims,data.trajectory_type);
	
	//
	if(recalc_dcf && (rcframes == 1)){
		dcf_calc(data);
	} else if (external_dcf) {
		ArrayRead(data.kw(0),dcffilename);
		for (int e = 1; e < rcencodes; e++) {
			data.kw(e) = data.kw(0);
		}
	}
	
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
		
	
	if(recalc_dcf && (rcframes > 1)){
		dcf_calc(data,gate);
	}


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
			wave = WAVELET3D( TinyVector<int,3>(rcxres,rcyres,rczres),TinyVector<int,3>(wavelet_levelsX,wavelet_levelsY,wavelet_levelsZ),WAVELET3D::WAVE_DB4);

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
	
	// TEMPX
	lranktime=LOWRANKCOIL(argc,argv);
			
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


void RECON::dcf_calc( MRI_DATA& data){
	
	// Setup Gridding + FFT Structure
	gridFFT dcf_gridding;
	dcf_gridding.kernel_type = KAISER_KERNEL;
	dcf_gridding.dwinX = dcf_dwin;
	dcf_gridding.dwinY = dcf_dwin;
	dcf_gridding.dwinZ = dcf_dwin;
	dcf_gridding.grid_x = 3.2;
	dcf_gridding.grid_y = 3.2;
	dcf_gridding.grid_z = 3.2;
	dcf_gridding.grid_in_x = 1;
	dcf_gridding.grid_in_y = 1;
	dcf_gridding.grid_in_z = 1;
	dcf_gridding.precalc_kernel(); // data.trajectory_dims,data.trajectory_type);
	/*dcf_gridding.grid_x = 1;
	dcf_gridding.grid_y = 1;
	dcf_gridding.grid_z = 1;
		*/
			
	/* Take square root of kernel*/
	for( int pos=0; pos< dcf_gridding.grid_filterX.length(firstDim); pos++){
		dcf_gridding.grid_filterX(pos) = sqrt( abs(dcf_gridding.grid_filterX(pos)))*( (dcf_gridding.grid_filterX(pos) > 0.0) ?( 1.0 ) : ( -1.0 )); 
	}
	
	for( int pos=0; pos< dcf_gridding.grid_filterY.length(firstDim); pos++){
		dcf_gridding.grid_filterY(pos) = sqrt( abs(dcf_gridding.grid_filterY(pos)))*( (dcf_gridding.grid_filterY(pos) > 0.0) ?( 1.0 ) : ( -1.0 )); 
	}
	
	for( int pos=0; pos< dcf_gridding.grid_filterZ.length(firstDim); pos++){
		dcf_gridding.grid_filterZ(pos) = sqrt( abs(dcf_gridding.grid_filterZ(pos)))*( (dcf_gridding.grid_filterZ(pos) > 0.0) ?( 1.0 ) : ( -1.0 )); 
	}
				
	 // Weighting Array for Time coding
	Array< float, 3 >Kweight(data.Num_Pts,data.Num_Readouts,data.Num_Slices,ColumnMajorArray<3>());
	Array< float, 3 >Kweight2(data.Num_Pts,data.Num_Readouts,data.Num_Slices,ColumnMajorArray<3>());
	
	// Array to grid to
	Array< float, 3 >Xtemp(4*rcxres,4*rcyres,4*rczres,ColumnMajorArray<3>());
	cout << "Size Xtemp = " << Xtemp.length(firstDim) << " x " << Xtemp.length(secondDim) << " x " << Xtemp.length(thirdDim) << endl;
	
	for(int e=0; e< rcencodes; e++){
		
		data.kx(e)*= dcf_scale;
		data.ky(e)*= dcf_scale;
		data.kz(e)*= dcf_scale;
		
		Kweight = 1;
		for(int iter=0; iter < dcf_iter; iter++){
			cout << "Iteration = " << iter << flush;			 
						
			Xtemp=0;
			dcf_gridding.grid_forward(Xtemp,Kweight,data.kx(e),data.ky(e),data.kz(e));
			dcf_gridding.grid_backward(Xtemp,Kweight2,data.kx(e),data.ky(e),data.kz(e));
			Kweight = Kweight / ( Kweight2);
			cout << " Sum = " << sum(Kweight) << endl;
		}
		data.kw(e) = Kweight;
		
		
		data.kx(e)/= dcf_scale;
		data.ky(e)/= dcf_scale;
		data.kz(e)/= dcf_scale;
		
		
	}
	ArrayWrite(data.kw(0),"Kweight_DCF.dat");
	
}

void RECON::dcf_calc( MRI_DATA& data, GATING& gate){
	
	// Setup Gridding + FFT Structure
	gridFFT dcf_gridding;
	dcf_gridding.kernel_type = KAISER_KERNEL;
	dcf_gridding.dwinX = dcf_dwin;
	dcf_gridding.dwinY = dcf_dwin;
	dcf_gridding.dwinZ = dcf_dwin;
	dcf_gridding.grid_x = 3.2;
	dcf_gridding.grid_y = 3.2;
	dcf_gridding.grid_z = 3.2;
	dcf_gridding.grid_in_x = 1;
	dcf_gridding.grid_in_y = 1;
	dcf_gridding.grid_in_z = 1;
	dcf_gridding.precalc_kernel(); 
			
	 // Weighting Array for Time coding
	Array< float, 3 >Kweight(data.Num_Pts,data.Num_Readouts,data.Num_Slices,ColumnMajorArray<3>());
	Array< float, 3 >Kweight2(data.Num_Pts,data.Num_Readouts,data.Num_Slices,ColumnMajorArray<3>());
	
	// Array to grid to
	Array< float, 3 >Xtemp(2*rcxres,2*rcyres,2*rczres,ColumnMajorArray<3>());

	for(int e=0; e< rcencodes; e++){
		data.kx(e)*= dcf_scale;
		data.ky(e)*= dcf_scale;
		data.kz(e)*= dcf_scale;

		for(int t=0; t < rcframes; t++){

			Kweight = 1;					 
			gate.weight_data( Kweight, e, data.kx(e),data.ky(e),data.kz(e),t,GATING::ITERATIVE,GATING::TIME_FRAME);

			for(int iter=0; iter < dcf_iter; iter++){
				cout << "Frame " << t << ", Iteration = " << iter << endl;			 

				Xtemp=0;
				dcf_gridding.grid_forward(Xtemp,Kweight,data.kx(e),data.ky(e),data.kz(e));
				dcf_gridding.grid_backward(Xtemp,Kweight2,data.kx(e),data.ky(e),data.kz(e));

				Kweight = Kweight / Kweight2;
			}


			cout << "size of data.kw = " << data.kw(e).shape() << endl;
			for(int s=0; s< Kweight.length(thirdDim); s++){
				for(int v=0; v< Kweight.length(secondDim); v++){
					for(int r=0; r< Kweight.length(firstDim); r++){
						float wt = Kweight(r,v,s);			  	
						if (abs(wt) > 0.0) {
							data.kw(e)(s,v,r) = wt;
						}
					}
				}
			}

		}
		data.kx(e)/= dcf_scale;
		data.ky(e)/= dcf_scale;
		data.kz(e)/= dcf_scale;
		
		
	}
	ArrayWrite(data.kw(0),"Kweight_DCF.dat");
	
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
						  complex<double>scale_RhP(0,0);						  
						  
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
									  if(intensity_correction){
									  	gridding.backward(XX(store_t,e,coil),IntensityCorrection,diff_data,kxE,kyE,kzE,TimeWeight);
									  }else{
									  	gridding.backward(XX(store_t,e,coil),diff_data,kxE,kyE,kzE,TimeWeight);
									  }
									  									  									  									  
									  //Ex-d
									  Array< complex<float>,3>kdataC = data.kdata(e,coil); 
									  diff_data -= kdataC;
									  
									  //E'(Ex-d)
									  RR(t,e,coil)=complex<float>(0.0,0.0);
									  if(intensity_correction){
									  	gridding.forward( RR(store_t,e,coil),IntensityCorrection,diff_data,kxE,kyE,kzE,TimeWeight);
									  }else{
									  	gridding.forward( RR(store_t,e,coil),diff_data,kxE,kyE,kzE,TimeWeight);
									  }
									  
									  // Now Get Scale
									  P=0;
									   							  
									  // EE'(Ex-d)
									  diff_data=0;
									  if(intensity_correction){
									  	gridding.backward(RR(store_t,e,coil),IntensityCorrection, diff_data,kxE,kyE,kzE,TimeWeight);
									  }else{
									  	gridding.backward(RR(store_t,e,coil), diff_data,kxE,kyE,kzE,TimeWeight);
									  }
									  
									  //E'EE'(Ex-d)
									  if(intensity_correction){
									  	gridding.forward(P,IntensityCorrection,diff_data,kxE,kyE,kzE,TimeWeight);
								  	  }else{
									  	gridding.forward(P,diff_data,kxE,kyE,kzE,TimeWeight);
								  	  }
									  scale_RhP += conj_sum(P,RR(t,e,coil)); 
									  
									  cout << "Coil " << coil << " took " << T << endl;
								  }//Coils
								  
								  // cout << e << "," << t << "took " << T << "s" << endl;
							  }//Time
						  }//Encode

						  // Get Scaling Factor R'P / R'R 
						  complex<double>scale_RhR = complex<double>(ArrayEnergy(RR),0);

						  // Error check
						  if(iteration==0){
							  error0 = abs(scale_RhR);
						  }
						  cout << "Residue Energy =" << scale_RhR << " ( " << (abs(scale_RhR)/error0) << " % )  " << endl;

						  // Export R (across coils)	
						  Array<complex<float>,2>Rslice=RR(0,0,0)(all,all,RR(0,0,0).length(2)/2);
						  ArrayWriteMag(Rslice,"R.dat");						  
						  
						  
						  // Step in direction
						  complex<double>scale_double = (scale_RhR/scale_RhP);
						  complex<float>scale( real(scale_double),imag(scale_double) );
						  cout << "Scale = " << scale << endl;
						  
						  for(int coil=0; coil< data.Num_Coils; coil++){
						     for(int e=0; e< rcencodes; e++){
							  for(int t=0; t< Nt; t++){
						  		RR(t,e,coil) *= scale;
						  		XX(t,e,coil) -= RR(t,e,coil);
						  }}}
						  cout << "Took " << iteration_timer << " s " << endl;
												  
						  // Export X slice
						  {
						  	Array< float,2>Xslice;
							Xslice.setStorage(ColumnMajorArray<2>());
							Xslice.resize( rcxres,rcyres);
							Xslice =0;
							for(int coil=0; coil< data.Num_Coils; coil++){
																
								Xslice += norm( XX(0,0,coil)(all,all,RR(0,0,0).length(2)/2) );
						  	}
							Xslice = sqrt( Xslice );
							ArrayWriteAppend(Xslice,"X_mag.dat");
						  }
						  
						  				  
						  // ----------------------------------
						  //   Soft Thresh
						  // ----------------------------------
						  {
						  	cout << "Soft thresh" << endl;
							Array< Complex3D, 2 >X2 = XX(0,all,all);
						  	if(softthresh.getThresholdMethod() != TH_NONE){
							  L1_threshold(X2);
						  	}
						  
						  	 // Export X slice
						  	{	
						  	Array< float,2>Xslice;
							Xslice.setStorage(ColumnMajorArray<2>());
							Xslice.resize( rcxres,rcyres);
							Xslice =0;
							for(int coil=0; coil< data.Num_Coils; coil++){
								Xslice += norm( XX(0,0,coil)(all,all,RR(0,0,0).length(2)/2) );
						  	}
							Xslice = sqrt( Xslice );
							ArrayWriteAppend(Xslice,"X_mag.dat");
						  	}
						  }
						  
						  
						  // ---------------------------------
						  //   Clear 
						  // ---------------------------------
						 
						  iteration_timer.tic();
						  lrankcoil.update_threshold( XX,2);
						  cout << "Get thresh took " << iteration_timer << " s" << endl;
						  
						  iteration_timer.tic();
						  lrankcoil.thresh(XX,2);
						  cout << "Thresh took " << iteration_timer << " s" << endl;
						  
						  // Export X slice
						  {
						  	Array< float,2>Xslice;
							Xslice.setStorage(ColumnMajorArray<2>());
							Xslice.resize( rcxres,rcyres);
							Xslice =0;
							for(int coil=0; coil< data.Num_Coils; coil++){
								Xslice +=   norm( XX(0,0,coil)(all,all,RR(0,0,0).length(2)/2) );
						  	}
							Xslice = sqrt( Xslice );
							ArrayWriteAppend(Xslice,"X_mag.dat");
						  }
						  
						  
							
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
				 
				 // Initial Guess
				 
				 			 
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
				 // double error0=0.0;
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
						  complex<double>scale_RhP(0,0);						  

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
								  

								  if (data.Num_Coils > 1) {
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
								  }else{ 
									  // Ex
									  diff_data=0;
									  gridding.backward(X(store_t,e),diff_data,kxE,kyE,kzE,TimeWeight);

									  //Ex-d
									  Array< complex<float>,3>kdataC = data.kdata(e,0); 
									  diff_data -= kdataC;

									  //E'(Ex-d)
									  gridding.forward( R(store_t,e),diff_data,kxE,kyE,kzE,TimeWeight);

								  }
								  
								  
								  // L2 
								  if(iteration > 0){
								  	l2reg.regularize(R(store_t,e),X(store_t,e) );
								  }
								    
								  //Now Get Scale factor (for Cauchy-Step Size)
								  P=0;

								  if (data.Num_Coils > 1) {
									  for(int coil=0; coil< data.Num_Coils; coil++){
										  // EE'(Ex-d)
										  diff_data=0;
										  gridding.backward(R(store_t,e),smaps(coil), diff_data,kxE,kyE,kzE,TimeWeight);

										  //E'EE'(Ex-d)
										  gridding.forward(P,smaps(coil),diff_data,kxE,kyE,kzE,TimeWeight);
									  }//Coils
								  }else{
									  // EE'(Ex-d)
									  diff_data=0;
									  gridding.backward(R(store_t,e), diff_data,kxE,kyE,kzE,TimeWeight);

									  //E'EE'(Ex-d)
									  gridding.forward(P,diff_data,kxE,kyE,kzE,TimeWeight);
								  }
								  
								  // TV of Image
								  if(iteration > 0){
								  	l2reg.regularize(P,R(store_t,e));
								  }
								 
								  P*=conj(R(store_t,e));
								  
								  for( Array<complex<float>,3>::iterator riter=P.begin(); riter != P.end(); riter++){
								    	scale_RhP += complex< double>( real(*riter),imag(*riter)); 
								  }
								  cout << e << "," << t << "took " << T << "s" << endl;
							  }//Time
						  }//Encode
						  
						  // Get Scaling Factor R'P / R'R 
						  cout << "Calc residue" << endl << flush;
						  complex<double>scale_RhR = 0.0;
						  for( Array< Array<complex<float>,3>,2>::iterator riter =R.begin(); riter != R.end(); riter++){
						  		scale_RhR += complex<double>( ArrayEnergy( *riter ), 0.0);
						  }
						  
						  // Error check
						  if(iteration==0){
							  cout << "L2 set scale " << endl << flush;
						  	  error0 = abs(scale_RhR);
							  l2reg.set_scale(error0,X);
						  }
						  cout << "Residue Energy =" << scale_RhR << " ( " << (abs(scale_RhR)/error0) << " % )  " << endl;
						  cout << "RhP = " << scale_RhP << endl;
						  
						  // Export R
						  Array<complex<float>,2>Rslice=R(0,0)(all,all,R(0,0).length(2)/2);
						  ArrayWriteMagAppend(Rslice,"R.dat");						  

						  // Step in direction
						  complex<double>scale_double = (scale_RhR/scale_RhP);
						  complex<float>scale( real(scale_double),imag(scale_double) );
						  
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
						  
						  /* TEMP
						  lranktime.update_threshold(X,0);
						  lranktime.thresh(X,0);
						  ArrayWriteMagAppend(Xslice,"X_mag.dat");
						  */
						  
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
	
	
	int actual_cycle_spins=cycle_spins;
	if(cs_spatial_transform != WAVELET){
		actual_cycle_spins = 1;
	}
	
	
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
	
	
	if(cs_spatial_transform==WAVELET){
		cout << "WAVLET in Space" << endl;
		for( int cycle_spin = 0; cycle_spin < actual_cycle_spins; cycle_spin++){
	  			cout << "    Spin " << cycle_spin << endl;
				wave.random_shift();
				for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
						wave.forward(*miter);
				}
							  
				// Only update thresh once 
				if( cycle_spin==0){
					softthresh.update_threshold(X,wave, 1./(float)actual_cycle_spins );
				}
				softthresh.exec_threshold(X,wave);
				
				for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
					wave.backward(*miter);
				}
		}// Cycle Spinning
	}else{
		softthresh.update_threshold(X,wave,1.0);
		softthresh.exec_threshold(X,wave);

	}// Wavelet
						  
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

inline float sqr ( float x){ return(x*x);}

void RECON::calc_sensitivity_maps( int argc, char **argv, MRI_DATA& data){

	// ------------------------------------
	//  Get coil sensitivity map ( move into function)
	// ------------------------------------
	Range all = Range::all();	
	
	switch( recon_type){
	
		case(CLEAR):{
			
			if(intensity_correction){
			
			
				// Clear doesn't need anything, but might need to handle sensitivity correction
				IntensityCorrection.setStorage( ColumnMajorArray<3>());
				IntensityCorrection.resize(rcxres,rcyres,rczres);
				IntensityCorrection = 0.0;
									
				for(int e=0; e< 1;e++){
					for(int coil=0; coil< data.Num_Coils; coil++){
						cout << "Coil = " << coil  << " encode = " << e << endl;
					
						gridding.image = 0;
						gridding.forward(gridding.image,data.kdata(e,coil),data.kx(e),data.ky(e),data.kz(e),data.kw(e));
						IntensityCorrection += norm( gridding.image);
					}
				}
				IntensityCorrection = sqrt(IntensityCorrection);
			
				// Gaussian blur
				Array< float , 3 > Blur;
				Blur.setStorage( ColumnMajorArray<3>() );
				Blur.resize(rcxres,rcyres,rczres);
				normalized_gaussian_blur( IntensityCorrection, Blur, 20.0);
			
				IntensityCorrection = Blur;
			
			}
				
		}break;
		
		case(PILS):
		case(IST):
		case(FISTA):
		case(CG):{
			// Allocate Storage for Map	and zero	
			cout << "Allocate Sense Maps"  << endl << flush;
			smaps.setStorage( ColumnMajorArray<1>());
			smaps.resize( data.Num_Coils );
			for(int coil=0; coil< smaps.length(firstDim); coil++){
				smaps(coil).setStorage( ColumnMajorArray<3>());
				smaps(coil).resize(rcxres,rcyres,rczres,data.Num_Coils);
				smaps(coil)=0;
			} 
			
			if(data.Num_Coils ==1){
				smaps(0)=complex<float>(1.0,0.0);
				break;
			}
						
			cout << "Recon Low Resolution Images"  << endl<< flush; 
			
			
			// Recon the maps
			for(int e=0; e< 1;e++){
				for(int coil=0; coil< data.Num_Coils; coil++){
					
					cout << "Coil = " << coil  << " encode = " << e << endl;
					if( !iterative_smaps){
						
						// Low Pass filtering for Sensitivity Map
						if(coil_combine_type!=ESPIRIT){
							gridding.k_rad = smap_res;
						}
						
						// Simple gridding
						gridding.forward( smaps(coil),  data.kdata(e,coil), data.kx(e), data.ky(e), data.kz(e) ,data.kw(e) );
						
						gaussian_blur(smaps(coil),extra_blurX,extra_blurY,extra_blurZ); // TEMP						
																								
					}else{
						// Regularized 
						single_coil_cg( smaps(coil),  data.kdata(e,coil), data.kx(e), data.ky(e), data.kz(e) ,data.kw(e) );
					}
				}
			}
			gridding.k_rad = 9999;
			
			Array< float , 3 > IC;
			if( intensity_correction ){
				intensity_correct( IC, smaps);
			}
		
			// Spirit Code
			switch(coil_combine_type){
		
				case(ESPIRIT):{
 					// Espirit (k-space kernal Eigen method)
					cout << "eSPIRIT Based Maps"  << endl; 
					SPIRIT S;
      	 				S.read_commandline(argc,argv);
      					S.init(rcxres,rcyres,rczres,data.Num_Coils);
		    			S.generateEigenCoils(smaps);		  
				}break;
		
				case(WALSH):{
					// Image space eigen Method
					eigen_coils(smaps);
				}break;
		
				case(LOWRES):{ // E-spirit Code 
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
					
				}break;
			} // Normalization
	
			if( intensity_correction ){
				cout << "Intensity correcting maps" << endl << flush;
				for(int coil=0; coil< data.Num_Coils; coil++){
					#pragma omp parallel for 
					for(int k=0; k<smaps(0).length(thirdDim); k++){
					for(int j=0; j<smaps(0).length(secondDim); j++){
					for(int i=0; i<smaps(0).length(firstDim); i++){
						smaps(coil)(i,j,k) *= IC(i,j,k);
					}}}
				}
				ArrayWrite(IC,"Intensity.dat");
			}
		
		
			// Export 
			if(export_smaps==1){
				cout << "Exporting Smaps" << endl;
				for(int coil=0; coil< smaps.length(firstDim); coil++){
					char name[256];
					sprintf(name,"SenseMaps_%2d.dat",coil);
					ArrayWrite(smaps(coil),name);
				}
			}
		
		}break;
		
		
		case(SOS):{
			// Allocate Storage for Map	and zero	
			cout << "Allocate Sense Maps"  << endl << flush;
			smaps = Alloc4DContainer< complex<float> >(rcxres,rcyres,rczres,data.Num_Coils);
			for(int coil=0; coil< smaps.length(firstDim); coil++){
				smaps(coil)=complex<float>(1.0,0.0);
			}	
		
		}break;
	}
} 


void RECON::single_coil_cg( Array< complex<float>,3> & X, Array< complex<float>,3> &kdata, Array<float,3> &kx, Array<float,3> &ky, Array<float,3> &kz, Array<float,3> &kw){

	// ------------------------------------
	// Conjugate Gradient Recon 
    //   -Uses much more memory than gradient descent but is faster (both convergence + per iteration)
	// ------------------------------------
	
	L2REG l2reg_coil;
	l2reg_coil.lambda  = 1.0;
	l2reg_coil.l2_type = L2REG::LOWRES;
				
	 // Structures  	
	Array< complex<float>,3> R( rcxres, rcyres, rczres,ColumnMajorArray<3>() );
	
	Array< complex<float>,3> P( rcxres, rcyres, rczres,ColumnMajorArray<3>() );
	P = complex<float>(0.0,0.0);

	Array< complex<float>,3> LHS( rcxres, rcyres, rczres,ColumnMajorArray<3>() );
	LHS = complex<float>(0.0,0.0);
				
	// Zeropadded bluring
	Array< complex<float>,3> ZeroPad( rcxres+64, rcyres+64, rczres+64,ColumnMajorArray<3>() );
	
	// Storage for (Ex-d)
	Array< complex<float>,3 >diff_data;
	diff_data.setStorage( ColumnMajorArray<3>() );
	diff_data.resize( kdata.shape());
	
	// Get LHS
	R = complex<float>(0.0,0.0);
	gridding.forward( R,kdata,kx,ky,kz,kw);	
	R = -R;
	P = R;	
		
			 
	// Now Iterate
	cout << "Iterate" << endl;
	for(int iteration =0; iteration< 40; iteration++){

			// E'Ex
			diff_data= complex<float>(0.0,0.0);
			LHS = complex<float>(0.0,0.0);
			gridding.backward( P ,diff_data,kx,ky,kz,kw);
			gridding.forward( LHS ,diff_data,kx,ky,kz,kw);
			
			// Get scale
			if(iteration==0){
				float sum_P_P = sum(norm(P));
				float sum_L_L = sum(norm(LHS));
				l2reg_coil.reg_scale = l2reg_coil.lambda*sqrt(sum_L_L/sum_P_P);
				cout << "L2 Scale = " << l2reg_coil.reg_scale << endl;
			}
			l2reg_coil.regularize(LHS,P);
								
			complex< float> sum_R0_R0(0.0,0.0);
			complex< float> sum_R_R(0.0,0.0);
			complex< float> sum_P_LHS(0.0,0.0);
			
			// Get scale
			for(int k=0; k < rczres; k++){
				for(int j=0; j < rcyres; j++){
						for(int i=0; i < rcxres; i++){
							sum_R0_R0 += norm( R(i,j,k) );
							sum_P_LHS += conj( P(i,j,k) )*LHS(i,j,k);
			}}}
			complex< float> scale = sum_R0_R0 / sum_P_LHS; 
			
			// Step
			for(int k=0; k < rczres; k++){
			for(int j=0; j < rcyres; j++){
			for(int i=0; i < rcxres; i++){
					X(i,j,k) += ( scale* (  P(i,j,k)) );
					R(i,j,k) -= ( scale* (LHS(i,j,k)) );
					sum_R_R += norm( R(i,j,k) );
			}}}
			
			cout << "Sum R'R = " << sum_R_R << endl;
			complex< float> scale2 = sum_R_R / sum_R0_R0;
			
			
			for(int k=0; k < rczres; k++){
			for(int j=0; j < rcyres; j++){
			for(int i=0; i < rcxres; i++){
					P(i,j,k) = R(i,j,k) + ( scale2*P(i,j,k) );
			}}}
			
			// Export X slice
			Array<complex<float>,2>Xslice=X(Range::all(),Range::all(),X.length(2)/2);
			ArrayWriteMagAppend(Xslice,"CG_mag.dat");
		
		}//iteration	
		

					
}

void  RECON::intensity_correct( Array< float, 3> & IC, Array< Array< complex<float>,3>,1 > &smaps ){

	cout << "Correcting Intensity" << endl;
	
	IC.setStorage( ColumnMajorArray<3>() );
	IC.resize(rcxres,rcyres,rczres);
			
	// Collect a Sum of Squares image
	IC = 0.0;
	for(int coil=0; coil< smaps.length(firstDim); coil++){
		IC += norm( smaps(coil ));
	}
	for( Array< float,3>::iterator miter=IC.begin(); miter!=IC.end(); miter++){
		(*miter) = sqrtf( *miter );
	}
	IC /= max(IC);
			
	// Threshold that image - to reduce error with air
	float max_sos = max(IC);
	float ic_threshold = 0.05 * max_sos;
	for( Array< float,3>::iterator miter=IC.begin(); miter!=IC.end(); miter++){
		if(  (*miter) < ic_threshold){
		 	(*miter) = 0.0;
		}
	}
			
	// Gaussian blur
	Array< float , 3 > Blur;
	Blur.setStorage( ColumnMajorArray<3>() );
	Blur.resize(rcxres,rcyres,rczres);
	normalized_gaussian_blur( IC, Blur, 20.0);
									
	// Calc intensity correction
	float max_blur = max(Blur);
	for(int k=0; k < rczres; k++){
		for(int j=0; j < rcyres; j++){
			for(int i=0; i < rcxres; i++){
				switch(recon_type){
					case(SOS):
					case(PILS):{
						IC(i,j,k) = Blur(i,j,k) / ( Blur(i,j,k)*Blur(i,j,k) + 0.01*max_blur*max_blur);
					}break;
					
					default:{
						IC(i,j,k) = ( Blur(i,j,k) > 0.05*max_blur) ? ( Blur(i,j,k) ) : ( 0.0);
					}break;
				}
	}}}
	return;
}
	
void RECON::gaussian_blur( Array< complex<float> , 3> & In, float sigmaX, float sigmaY, float sigmaZ){
	
	// Extent of kernel
	int dwinX = 3*(int)sigmaX;
	int dwinY = 3*(int)sigmaY;
	int dwinZ = 3*(int)sigmaZ;
	
	if( sigmaX > 0){
		// Kernel to reduce calls to exp
		Array< float,1> kern(2*dwinX+1);	
		for(int t=0; t< (2*dwinX +1); t++){
			kern(t) = exp( -sqr( (float)t - dwinX ) / (2.0*sqr(sigmaX))); 
		}				
			
		// Gaussian Blur in X
		cout << "Blur in X" << endl;
		#pragma omp parallel for
		for(int k=0; k < rczres; k++){
		for(int j=0; j < rcyres; j++){
		
		Array< complex<float>, 1>TEMP( rcxres);
		for(int i=0; i < rcxres; i++){
			TEMP(i) = In(i,j,k);
			In(i,j,k) = complex<float>(0.0,0.0);
		}
		
		for(int i=0; i < rcxres; i++){
			int sx = max( i - dwinX, 0);	
			int ex = min( i + dwinX, rcxres);
			int ks = sx - (i-dwinX); 
						
			for(int ii=sx; ii<ex; ii++){
				In(i,j,k) += kern(ks)*TEMP(ii);
				ks++;
			}
		}}}
	}
	
	if( sigmaY > 0){
		// Kernel to reduce calls to exp
		Array< float,1> kern(2*dwinY+1);	
		for(int t=0; t< (2*dwinY +1); t++){
			kern(t) = exp( -sqr( (float)t - dwinY ) / (2.0*sqr(sigmaY))); 
		}				
			
		// Gaussian Blur in Y
		cout << "Blur in Y" << endl;
		#pragma omp parallel for
		for(int k=0; k < rczres; k++){
		for(int i=0; i < rcxres; i++){
				
		Array< complex<float>, 1>TEMP( rcyres);
		for(int j=0; j < rcyres; j++){
			TEMP(j) = In(i,j,k);
			In(i,j,k) = complex<float>(0.0,0.0);
		}
		
		for(int j=0; j < rcyres; j++){
			int sy = max( j - dwinY, 0);	
			int ey = min( j + dwinY, rcyres);
			int ks = sy - (j-dwinY);
				
			for(int ii=sy; ii<ey; ii++){
				In(i,j,k) += kern(ks)*TEMP(ii);
				ks++;
			}
		}
		
		}}
	}	
	
	if( sigmaZ > 0){
		// Kernel to reduce calls to exp
		Array< float,1> kern(2*dwinZ+1);	
		for(int t=0; t< (2*dwinZ +1); t++){
			kern(t) = exp( -sqr( (float)t - dwinZ ) / (2.0*sqr(sigmaZ))); 
		}				
			
		// Gaussian Blur in Y
		cout << "Blur in Z" << endl;
		#pragma omp parallel for
		for(int j=0; j < rcyres; j++){
		for(int i=0; i < rcxres; i++){
				
		Array< complex<float>, 1>TEMP( rczres);
		for(int k=0; k < rczres; k++){
			TEMP(k) = In(i,j,k);
			In(i,j,k) = complex<float>(0.0,0.0);
		}
		
		for(int k=0; k < rczres; k++){
			int sz =  max( k - dwinZ, 0);	
			int ez =  min( k + dwinZ, rczres);
			int ks = sz - (k-dwinZ);
			
			for(int ii=sz; ii<ez; ii++){
				In(i,j,k) += kern(ks)*TEMP(ii);
				ks++;
			}
		}
		
		}}
	}
}	



void RECON::normalized_gaussian_blur( const Array< float, 3> & In, Array< float, 3> & Out, float sigma){

	int dwinX = 3*(int)sigma;
	int dwinY = 3*(int)sigma;
	int dwinZ = 3*(int)sigma;
	
	Array< float,3> Normalization;
	Normalization.setStorage( ColumnMajorArray<3>() );
	Normalization.resize(rcxres,rcyres,rczres);
	Normalization = 0.0;
	
	// Set to zero
	Out = 0.0;
	
	// Kernel to reduce calls to exp
	float *kern = new float [2*dwinX+1];	
	for(int t=0; t< (2*dwinX +1); t++){
		kern[t] = exp( -sqr( (float)t - dwinX ) / (2.0*sqr(sigma))); 
		// cout << "Kern (" << t << ") = " << kern[t] << endl;
	}				
			
	// Gaussian Blur in X
	cout << "Blur in X" << endl;
	#pragma omp parallel for
	for(int k=0; k < rczres; k++){
	for(int j=0; j < rcyres; j++){
	for(int i=0; i < rcxres; i++){
		int ks = 0;
		int sx =  i - dwinX;	
		if(sx < 0){
			ks = -sx;
			sx = 0;
		}		
		int ex = min( i + dwinX, rcxres);
		
		for(int ii=sx; ii<ex; ii++){
			Out(i,j,k) += kern[ks]*In(ii,j,k);
			Normalization(i,j,k) += kern[ks]*(float)( In(ii,j,k) > 0 );
			ks++;
		}
	}}}			
	
	// Gaussian Blur in Y
	cout << "Blur in Y" << endl;
	#pragma omp parallel for
	for(int k=0; k < rczres; k++){
	for(int i=0; i < rcxres; i++){
		// Copy
		float *Line = new float[rcyres];
		float *NLine = new float[rcyres];
		for(int j=0; j < rcyres; j++){
			Line[j]= Out(i,j,k);
			NLine[j] = Normalization(i,j,k);
		}
	
		// Convolve
		for(int j=0; j < rcyres; j++){
		int ks=0;
		int sy = j - dwinY;
		if(sy<0){
		 	ks = -sy; // Offset;
		 	sy = 0;
		}			
		int ey = min( j + dwinY+1, rcyres);
		
		for(int jj=sy; jj<ey; jj++){
			Out(i,j,k) += kern[ks]*Line[jj];
			Normalization(i,j,k) += kern[ks]*NLine[jj];  // Already a count
			ks++;
		}
		}
		
		delete [] Line;
		delete [] NLine;
		
	}}			
	
	// Gaussian Blur in Z
	cout << "Blur in Z" << endl;
	#pragma omp parallel for
	for(int j=0; j < rcyres; j++){
	for(int i=0; i < rcxres; i++){
		// Copy
		float *Line = new float[rczres];
		float *NLine = new float[rczres];
		for(int k=0; k < rczres; k++){
			Line[k]= Out(i,j,k);
			NLine[k] = Normalization(i,j,k);
		}
	
		// Convolve
		for(int k=0; k < rczres; k++){
			int ks=0;
			int sz = k - dwinZ;
			if(sz<0){
		 		ks = -sz; // Offset;
		 		sz = 0;
			}			
			int ez = min( k + dwinZ+1, rczres);
		
			for(int kk=sz; kk<ez; kk++){
				Out(i,j,k) += kern[ks]*Line[kk];
				Normalization(i,j,k) += kern[ks]*NLine[kk];
				ks++;
			}
			
			if( Normalization(i,j,k) > 0){
				Out(i,j,k)/= Normalization(i,j,k);
			}
			
		}
		
		delete [] Line;
		delete [] NLine;
	}}
	
}


/**
 * Generate coils from SPIRiT kernel.
 */
void RECON::eigen_coils( Array< Array< complex<float>,3 >,1 > &image)
{
    	
	typedef Array<complex<float>,3> Complex3D;
	
	// Shorthand
	int Ncoils = image.extent(firstDim);
	int Nx =image(0).extent(firstDim);
	int Ny =image(0).extent(secondDim);
	int Nz =image(0).extent(thirdDim);
	
	int block_size_x = walsh_block_sizeX;
	int block_size_y = walsh_block_sizeY;
	int block_size_z = walsh_block_sizeZ;
		
	// Blocks shouldn't be larger than the dimension
	block_size_x = ( block_size_x > Nx) ? ( Nx ) : ( block_size_x );
	block_size_y = ( block_size_y > Ny) ? ( Ny ) : ( block_size_y );
	block_size_z = ( block_size_z > Nz) ? ( Nz ) : ( block_size_z );
	
	int Np = block_size_x*block_size_y*block_size_z;
	
	cout << "Actual Block Size = " << block_size_x << " x " << block_size_y << " x " << block_size_z << endl;	
	cout << "Getting Low Rank threshold (N=" << Ncoils << ")(Np = " << Np << ")" << endl;
	
	int block_Nx= (int)( Nx / block_size_x );
	int block_Ny= (int)( Ny / block_size_y );
	int block_Nz= (int)( Nz / block_size_z );
	
	int total_blocks = block_Nx * block_Ny * block_Nz;
	cout << "Total Block Size" << total_blocks << " ( " << block_Nx << "," << block_Ny << "," << block_Nz << ")" << endl;
	
		
	#pragma omp parallel for
	for( int block = 0; block < total_blocks; block++){
		
		// Nested parallelism workaround
		int i = (int)( block % block_Nx);
		
		int temp = (int)( (float)block / (float)block_Nx);
		int j = (int)(        temp % block_Ny );
		int k = (int)( (float)temp / (float)(block_Ny) );
		
		k*= block_size_z;
		j*= block_size_y;
		i*= block_size_x;
				
		arma::cx_fmat A;
		A.zeros(Np,Ncoils); // Pixels x Coils
	
		arma::cx_fmat U;
		U.zeros(Np,Np); // Pixels x Pixels
	    
		arma::fmat S;
		S.zeros(Np,Ncoils); // Pixels x Coils
	    		
		arma::cx_fmat V;
		V.zeros(Ncoils,Ncoils); // Coils x Coils
		
		//-----------------------------------------------------
		//   Block section
		//-----------------------------------------------------
		
		for(int c=0; c< Ncoils; c++){
			int count=0;
			for(int kk=k; kk < k+block_size_z; kk++){
			for(int jj=j; jj < j+block_size_y; jj++){
			for(int ii=i; ii < i+block_size_x; ii++){
				A(count,c) = image(c)(ii,jj,kk);
				count++;
			}}}
		}
		
				
		// SVD
		arma::fvec s;
  		arma::svd(U,s,V,A);
		
		//  V.print("V");
		arma::cx_fvec sens = V.col(0);
		
		for(int c=0; c< Ncoils; c++){
			int count=0;
			for(int kk=k; kk < k+block_size_z; kk++){
			for(int jj=j; jj < j+block_size_y; jj++){
			for(int ii=i; ii < i+block_size_x; ii++){
				image(c)(ii,jj,kk) = conj( sens(c,0));
				count++;
			}}}
		}
		
	}// Block (threaded)
}







