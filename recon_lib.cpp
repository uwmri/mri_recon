
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
	rcframes=-1;
	rcencodes=1;
	rc_frame_start = 0;
	
	
	smap_res=8;
	intensity_correction = false;
	iterative_smaps = false;
	reset_dens=false;
	
	threads = -1;	
	acc = 1;
	compress_coils = 0.0;
	whiten = false;
	export_smaps = 0;
	max_iter = 50;
	smap_use_all_encodes = false;
	smap_nex_encodes = false;
	smap_thresh = 0.0;
	smap_mask = SMAPMASK_NONE;
	
    coil_rejection_flag = false;
    coil_rejection_radius = 0.5;
    coil_rejection_shape = 1; 

	cycle_spins = 4;

	walsh_block_sizeX = 8;
	walsh_block_sizeY = 8;
	walsh_block_sizeZ = 8;
	
	// Gaussian Blur of smaps
	extra_blurX = 0.0;
	extra_blurY = 0.0;
	extra_blurZ = 0.0;
	
	// Gaussian Blur of smaps
	blurX = 0.0;
	blurY = 0.0;
	blurZ = 0.0;
	
	prep_done = false;
	
	wavelet_levelsX=4; 
	wavelet_levelsY=4; 
	wavelet_levelsZ=4;
	
	// Code to rotate
	phase_rotation = false;
	phase_rotation_sX = 20;
	phase_rotation_sY = 20;
	phase_rotation_sZ = 20;
		
	pregate_data_flag = false;
	
	dcf_type = SUPPLIED;
	dcf_iter = 20;
	dcf_dwin = 2.1;
	dcf_scale = 1.0;
	dcf_overgrid = 2.1;
	dcf_acc = 1.0;
	
	admm_gamma = 0.1;
	admm_max_iter = 20;
	admm_rho = 0.5;
	
	step_update_frequency = 1;
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
	help_flag("-zoom_x []","zoom factor in x (external data)");
	help_flag("-zoom_y []","zoom factor in y (external data)");
	help_flag("-zoom_z []","zoom factor in z (external data)");
	
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
	help_flag("-dcf_acc []","Increase size of kernel at edge of k-space by this amount");
	/* Shouldn't be here. It's inherently external to this code
	help_flag("-external_dcf","Use (precomputed) DCF from external file");
	help_flag("-dcf_file []","Filename for external DCF usage");
	*/
	
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
		int_flag("-rc_frame_start",rc_frame_start);
		
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
		trig_flag(ADMM,"-admm",recon_type);
		
		trig_flag(RECALC_DCF,"-recalc_dcf",dcf_type);
		trig_flag(RECALC_VOR,"-recalc_vor",dcf_type);
		int_flag("-dcf_iter",dcf_iter);
		float_flag("-dcf_dwin",dcf_dwin);
		float_flag("-dcf_scale",dcf_scale);
		float_flag("-dcf_overgrid",dcf_overgrid);
		float_flag("-dcf_acc",dcf_acc);
		trig_flag(true,"-reset_dens",reset_dens);
		
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
				trig_flag(DFT,"dft",cs_encode_transform);
				trig_flag(PCA,"pca",cs_encode_transform);
				trig_flag(WAVELET,"wavelet",cs_encode_transform);
			}else{
				cout << "Please provide encode transform type..none/dft/wavelet/diff" << endl;
				exit(1);
			}
		int_flag("-wavelet_levelsX",wavelet_levelsX);
		int_flag("-wavelet_levelsY",wavelet_levelsY);
		int_flag("-wavelet_levelsZ",wavelet_levelsZ);
		float_flag("-blurX",blurX);
		float_flag("-blurY",blurY);
		float_flag("-blurZ",blurZ);
		trig_flag(true,"-phase_rotation",phase_rotation);
						
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
		trig_flag(true,"-smap_use_all_encodes",smap_use_all_encodes);
		trig_flag(true,"-smap_nex_encodes",smap_nex_encodes);
		float_flag("-smap_thresh",smap_thresh);
        trig_flag(true,"-coil_rejection",coil_rejection_flag);
        float_flag("-coil_rejection_radius", coil_rejection_radius);
        int_flag("-coil_rejection_shape", coil_rejection_shape);
		int_flag("-step_update_frequency",step_update_frequency);
	
		// Encode Transforms
		}else if(strcmp("-smap_mask",pstring[pos]) == 0) {
			pos++;
			if( pos==numarg){
				cout << "Please provide sensitivty map..none/circle/sphere" << endl;
				exit(1);
				trig_flag(SMAPMASK_NONE,"none",smap_mask);
				trig_flag(SMAPMASK_CIRCLE,"circle",smap_mask);
				trig_flag(SMAPMASK_SPHERE,"sphere",smap_mask);
			}else{
				cout << "Please provide sensitivty map..none/circle/sphere" << endl;
				exit(1);
			}
								
		// Source of data
		trig_flag(EXTERNAL,"-external_data",data_type);
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
		
		float_flag("-admm_alpha",admm_alpha);
		float_flag("-admm_gamma",admm_gamma);
		int_flag("-admm_max_iter",admm_max_iter);
		float_flag("-admm_rho",admm_rho);
		
		trig_flag(true,"-pregate_data",pregate_data_flag);
		
		int_flag("-threads",threads);
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

	// Matlab like timer (openmp code base)
	tictoc T; 
	
	// Setup Gridding + FFT Structure
	gridding.read_commandline(argc,argv);
	gridding.precalc_gridding(rczres,rcyres,rcxres,data);
			
	// Recalculate the Density compensation
	switch(dcf_type){
		
		default:{
			// Use supplied
		}break;
		
		case(RECALC_DCF):{
			dcf_calc(data);
		}break;
		
		case(RECALC_VOR):{
			if( data.trajectory_type(2) == MRI_DATA::NONCARTESIAN){
				VORONOI_DCF::vor_dcf(data.kw(0),data.kx(0),data.ky(0),data.kz(0),VORONOI_DCF::SPHERE);
			}else{
				VORONOI_DCF::vor_dcf(data.kw(0),data.kx(0),data.ky(0),data.kz(0),VORONOI_DCF::CYLINDER);
			}
		}break;
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
		data.Num_Encodings -= 1;
		
		
		/* Need to resize Gating
		data.time.resizeAndPreserve(data.Num_Encodings); 
		data.resp.resizeAndPreserve(data.Num_Encodings); 
		data.prep.resizeAndPreserve(data.Num_Encodings); 
		data.ecg.resizeAndPreserve(data.Num_Encodings); 
		*/
	}
		
	// -------------------------------------
	//	This handles all the gating, assuming mri_data physio data is populated 
	// -------------------------------------
	gate=GATING(argc,argv);
	gate.init( data,&rcframes);
	
	if(pregate_data_flag){
		pregate_data( data ); 	
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
		case(CLEAR):
		case(ADMM):{
			
			// Setup 3D Wavelet
			wave = WAVELET3D( TinyVector<int,3>(rcxres,rcyres,rczres),TinyVector<int,3>(wavelet_levelsX,wavelet_levelsY,wavelet_levelsZ),WAVELET3D::WAVE_DB4);

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


/**
 * This function parses an MRI_DATA structure into time frames such that this is not required during recon
 * @see setup()
 * @param data MRI_DATA to be transfromed
 */
void RECON::pregate_data( MRI_DATA &data){
	
	//data.stats();
	
	// For now balloon memory
	MRI_DATA data2;
	
	data2.Num_Encodings = data.Num_Encodings * rcframes;
	data2.Num_Coils = data.Num_Coils;
	
	// Resize the data structures but don't allocate sub-structure
	data2.kx.resize( data2.Num_Encodings);
	data2.ky.resize( data2.Num_Encodings);
	data2.kz.resize( data2.Num_Encodings);
	data2.kw.resize( data2.Num_Encodings);
	data2.kdata.setStorage( ColumnMajorArray<2>());
	data2.kdata.resize(data2.Num_Encodings,data2.Num_Coils);
		
	int count = 0;
	for( int e=0; e<data.Num_Encodings; e++){
		
		Array< float, 3 >Kweight;
		Kweight.setStorage( ColumnMajorArray<3>());
		Kweight.resize( data.kx(e).shape());
			
		for( int t= 0; t < rcframes; t++){
		
			// First count th number of points required
			Kweight = 1;					 
			gate.weight_data( Kweight, e, data.kx(e),data.ky(e),data.kz(e),t,GATING::NON_ITERATIVE,GATING::TIME_FRAME);
			
			// Count the number of frames
			int number_of_points=0;
			for( Array< float,3>::iterator miter=Kweight.begin(); miter!=Kweight.end(); miter++){
				if( *miter > 0){
					number_of_points++;
				}
			}
			float Kw_Scale = (1e6/((float)number_of_points));
			cout << "Frame " << t << ", Total number of point = " << number_of_points << "Kw_scale = " << Kw_Scale << endl;
			
			data2.kx(count).setStorage( ColumnMajorArray<3>());
			data2.kx(count).resize(number_of_points,1,1);
			
			data2.ky(count).setStorage( ColumnMajorArray<3>());
			data2.ky(count).resize(number_of_points,1,1);
			
			data2.kz(count).setStorage( ColumnMajorArray<3>());
			data2.kz(count).resize(number_of_points,1,1);
			
			data2.kw(count).setStorage( ColumnMajorArray<3>());
			data2.kw(count).resize(number_of_points,1,1);
			
			for(int coil=0; coil < data2.Num_Coils; coil++){
				data2.kdata(count,coil).setStorage( ColumnMajorArray<3>());
				data2.kdata(count,coil).resize(number_of_points,1,1);
			}
			
			int point_number=0;
			for( int k=0; k< Kweight.length(thirdDim); k++){
				for( int j=0; j< Kweight.length(secondDim); j++){
					for( int i=0; i< Kweight.length(firstDim); i++){
						
						if( Kweight(i,j,k) > 0){
					
							data2.kx(count)(point_number,0,0) = data.kx(e)(i,j,k);
							data2.ky(count)(point_number,0,0) = data.ky(e)(i,j,k);
							data2.kz(count)(point_number,0,0) = data.kz(e)(i,j,k);
							if(reset_dens){
								data2.kw(count)(point_number,0,0) = 1.0;
							}else{
								data2.kw(count)(point_number,0,0) = Kw_Scale*data.kw(e)(i,j,k);
							}
							
							for(int coil=0; coil<data2.Num_Coils; coil++){
								data2.kdata(count,coil)(point_number,0,0) = data.kdata(e,coil)(i,j,k);
							}
							point_number++;
						}
			}}}	
			
			count++;		
		}
		
	}
	
	cout << "Swap" << endl << flush;
	cycleArrays( data.kx, data2.kx);
	cycleArrays( data.ky, data2.ky);
	cycleArrays( data.kz, data2.kz);
	cycleArrays( data.kw, data2.kw);
	cycleArrays( data.kdata, data2.kdata);
	
	data.Num_Encodings = data2.Num_Encodings;
	// data.stats();
		
	cout << "Done gating data" << endl << flush;
}	



void vor_dcf2(Array< Array<float,3 >,2 >&Kw,Array< Array<float,3 >,2 > &Ky,Array< Array<float,3 >,2 > &Kz){
	
	
	// Get total number of points
	int total_points = Kz.numElements()*Kz(0,0).length(secondDim)*Kz(0,0).length(thirdDim);
	
	cout << "Redoing with new DCF" << endl;
	Array< float, 3> ky(total_points,1,1,ColumnMajorArray<3>());
	Array< float, 3> kz(total_points,1,1,ColumnMajorArray<3>());
	Array< float, 3> kx(total_points,1,1,ColumnMajorArray<3>());
	kx = 1.0;
	Array< float, 3> kw(total_points,1,1,ColumnMajorArray<3>());	
	
	cout << "Copy points " << endl;
	int count =0;
	for( int e=0; e < Ky.length(firstDim); e++){
		for( int t=0; t < Ky.length(secondDim); t++){
			
			
			for(int jj = 0; jj < Ky(e,t).length(secondDim); jj++){
				for(int kk = 0; kk < Ky(e,t).length(thirdDim); kk++){
			
					ky(count,0,0) = Ky(e,t)(0,jj,kk);
					kz(count,0,0) = Kz(e,t)(0,jj,kk);
					count++;
			}}
		}
	}
	
	cout << "Compute Vor" << endl;
	VORONOI_DCF::vor_dcf( kw,ky,kz,kx,VORONOI_DCF::CUBE);	
	
	ArrayWrite(kz,"NEW_Kz.dat");
	ArrayWrite(ky,"NEW_Ky.dat");
	ArrayWrite(kw,"NEW_Kw.dat");
	
	count =0;
	for( int e=0; e < Ky.length(firstDim); e++){
		for( int t=0; t < Ky.length(secondDim); t++){
			
			
			for(int jj = 0; jj < Ky(e,t).length(secondDim); jj++){
				for(int kk = 0; kk < Ky(e,t).length(thirdDim); kk++){
					
					for(int ii=0; ii< Ky(e,t).length(firstDim); ii++){
						Kw(e,t)(ii,jj,kk) = kw(count,0,0);
					}
					count++;
			}}
		}
	}


}


Array< Array<complex<float>,3 >,2 > RECON::test_sms( MRI_DATA& data_cal, MRI_DATA& data, int argc, char **argv){
	
	
	// Input should be:
	//  data_cal     - A dataset without sms used to calc sensitvity maps
	//  data	 	 - A dataset with sms reconed in dynamic fashion
	//  passthrough for arguments
	
	int flex = 1;
	int cal_frames = data_cal.Num_Encodings/(flex+1);
	int sms_frames = data.Num_Encodings/(flex+1);
	int cal_encodes = flex + 1;
	int sms_encodes = flex + 1;
	int Ncoils = data.Num_Coils;
	int zpad = 8;  // Pad the slice to account for edge effect.
	int sms_factor = data.sms_factor;
		
	cout <<"Starting SMS " << endl;
	cout <<"Flex = " << flex << endl;
	cout <<"Input Frames (cal) = " << data_cal.Num_Encodings << endl;
	cout <<"Input Frames (dynamic) = " << data.Num_Encodings << endl;
	cout <<"Cal Frames = " << cal_frames << endl;
	cout <<"Sms Frames = " << sms_frames << endl;
	cout <<"Sms Factor = " << sms_factor << endl;
		
	// Use Native Resultion 
	rcxres = (rcxres == -1) ? ( data.xres ) : ( rcxres );
	rcyres = (rcyres == -1) ? ( data.yres ) : ( rcyres );
	rczres = (rczres == -1) ? ( data.zres + zpad ) : ( rczres + zpad);
	
	typedef Array<complex<float>,3> Complex3D;
	typedef Array< float,3> Float3D;
					
	// Allocate Sensiticty Maps
	Array< Array< complex<float>, 3>,1> temp =  Alloc4DContainer< complex<float> >(rcxres,rcyres,rczres,Ncoils);
	smaps.reference(temp);
		
	cout << "Allocate Sense Maps"  << endl << flush;
	if( Ncoils==1){
		smaps(0) = complex<float>(1.0,0.0);
	}else{
		
		int Num_Pts = data_cal.kx(0).length(firstDim);
		int Num_Readouts = data_cal.kx(0).length(secondDim);
		int Num_Slices = data_cal.kx(0).length(thirdDim);
				
		// Get the Calibration
		cout << "Copying Cal Data" << endl;
		Array< Complex3D,3> CalData = Alloc6DContainer< complex<float> >(Num_Pts,Num_Readouts,Num_Slices,cal_frames,cal_encodes,Ncoils);
		Array< Float3D,2> CalKx = Alloc5DContainer< float >(Num_Pts,Num_Readouts,Num_Slices,cal_frames,cal_encodes);
		Array< Float3D,2> CalKy = Alloc5DContainer< float >(Num_Pts,Num_Readouts,Num_Slices,cal_frames,cal_encodes); 
		Array< Float3D,2> CalKz = Alloc5DContainer< float >(Num_Pts,Num_Readouts,Num_Slices,cal_frames,cal_encodes);
		Array< Float3D,2> CalKw = Alloc5DContainer< float >(Num_Pts,Num_Readouts,Num_Slices,cal_frames,cal_encodes);
		Array< Float3D,3> CalZ  = Alloc6DContainer< float >(Num_Pts,Num_Readouts,Num_Slices,cal_frames,cal_encodes,1);
		
		for( int e=0; e < cal_encodes; e++){
			for( int t=0; t < cal_frames; t++){
				CalKx(t,e) = data_cal.kx(t*cal_encodes + e);
				CalKy(t,e) = data_cal.ky(t*cal_encodes + e);
				CalKz(t,e) = data_cal.kz(t*cal_encodes + e);
				CalKw(t,e) = data_cal.kw(t*cal_encodes + e);
				CalZ(t,e,0)= data_cal.z(t*cal_encodes + e,0);
				
				for(int coil =0; coil < Ncoils; coil++){
					CalData(t,e,coil) = data_cal.kdata(t*cal_encodes + e,coil);
					data_cal.kdata(t*cal_encodes + e,coil).resize(1,1,1);
				}
			}
		}
		
		// This calculates the weighting based on voronoi diagram
		vor_dcf2(CalKw,CalKy,CalKz);
				
		//  This is just to get sensitivty maps 
		smsEncode smsCgrid;
		smsCgrid.read_commandline(argc,argv);
		smsCgrid.sms_factor = 1;
		smsCgrid.precalc_gridding(rczres,rcyres,rcxres,cal_frames,cal_encodes,1,data);
	
		cout << "Alloc Image" << endl << flush;
		Array< Array< complex<float>,3>,2> image_store =  Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,cal_encodes,Ncoils);
			
		// Compute Forward Transform
		cout << "Create X" << endl << flush;
		Array< Array<complex<float>,3>,2> xdf = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,cal_frames,cal_encodes);
		
		Complex3D smaptest(rcxres,rcyres,rczres,ColumnMajorArray<3>());
		smaptest = complex<float>(1.0,0.0);
				
		// Now grid each frame
		for( int coil =0; coil < Ncoils; coil++){
			
			// Zero the storage
			for( int e= 0; e< cal_encodes; e++){
				for( int t =0; t<cal_frames; t++){
					xdf(t,e) = complex<float>(0.0,0.0);
				}
			}
			
			// Grid the data
			Array< Complex3D ,2> temp = CalData(Range::all(),Range::all(),coil);			 
			smsCgrid.forward( xdf, smaptest, temp,CalKx,CalKy,CalKz,CalKw,CalZ);
			
			// Sum over the coils
			for(int e = 0; e< cal_encodes; e++){
				image_store(e,coil) = complex<float>(0.0,0.0);
				for(int t=0; t<cal_frames; t++){
					image_store(e,coil) += xdf(t,e); 
				}
				
				{
					char name[1024];
					sprintf(name,"PreNormSmap_Encode_%d_Coil_%03d.dat",e,coil);
					ArrayWrite(image_store(e,coil),name);
				}
			}
		}
		
		float max_val = 0.0;
		for( Array< Array<complex<float>,3>,2>::iterator miter=image_store.begin(); miter!=image_store.end(); miter++){
			gaussian_blur(*miter,blurX,blurY,blurZ); // TEMP		
			float temp = max( abs(*miter));
			max_val += temp*temp;
		}
		
		// Spirit Code
		switch(coil_combine_type){
		
			case(ESPIRIT):{
 				// Espirit (k-space kernal Eigen method)
				cout << "eSPIRIT Based Maps"  << endl; 
				SPIRIT S;
      			S.read_commandline(argc,argv);
      			S.init(rcxres,rcyres,rczres,data.Num_Coils);
				S.generateEigenCoils(smaps, image_store);
			}break;
		
			case(WALSH):{
				// Image space eigen Method
				eigen_coils(smaps, image_store);
			}break;
		
			case(LOWRES):{
				
				
				float scale = 1/sqrt(max_val);
				cout << "Scale" << scale << endl;				
				for(int coil = 0; coil < Ncoils; coil++){
					smaps(coil) = scale*image_store(cal_encodes-1,coil);
				}
				sos_normalize(smaps);
			}
		}
	
		cout << "Writing Smaps " << endl << flush;
		for( int coil =0; coil < Ncoils; coil++){
			char name[1024];
			sprintf(name,"Smap%03d.dat",coil);
			ArrayWrite(smaps(coil),name);
		}
	}
	
	int Num_Pts = data.kx(0).length(firstDim);
	int Num_Readouts = data.kx(0).length(secondDim);
	int Num_Slices = data.kx(0).length(thirdDim);
	
	// Get the Dynamic Data
	cout << "Alloc Dynamic Data" << endl;
	Array< Complex3D,3> SmsData = Alloc6DContainer< complex<float> >(Num_Pts,Num_Readouts,Num_Slices,sms_frames,sms_encodes,Ncoils);
	Array< Float3D,2> SmsKx = Alloc5DContainer< float >(Num_Pts,Num_Readouts,Num_Slices,sms_frames,sms_encodes);
	Array< Float3D,2> SmsKy = Alloc5DContainer< float >(Num_Pts,Num_Readouts,Num_Slices,sms_frames,sms_encodes); 
	Array< Float3D,2> SmsKz = Alloc5DContainer< float >(Num_Pts,Num_Readouts,Num_Slices,sms_frames,sms_encodes);
	Array< Float3D,2> SmsKw = Alloc5DContainer< float >(Num_Pts,Num_Readouts,Num_Slices,sms_frames,sms_encodes);
	Array< Float3D,3> SmsZ  = Alloc6DContainer< float >(Num_Pts,Num_Readouts,Num_Slices,sms_frames,sms_encodes,sms_factor);
	
	cout << "Copy Dynamic Data (sms_factor = " << data.sms_factor << ", coils = " << Ncoils << ")" << endl;
	for( int e=0; e < sms_encodes; e++){
		for( int t=0; t < sms_frames; t++){
			int offset = t*sms_encodes + e;
		
			SmsKx(t,e) = data.kx(offset);
			SmsKy(t,e) = data.ky(offset);
			SmsKz(t,e) = data.kz(offset);
			SmsKw(t,e) = data.kw(offset);
			
			for(int sms_pos=0; sms_pos < data.sms_factor; sms_pos++){			
				SmsZ(t,e,sms_pos) = data.z(offset,sms_pos);
			}
			
			for(int coil =0; coil < Ncoils; coil++){
				SmsData(t,e,coil) = data.kdata(offset,coil);
				data.kdata(offset,coil).resize(1,1,1);
			}
		}
	}

	cout << "Prep Thresholds" << endl << flush;
	{
			// Setup 3D Wavelet
			wave = WAVELET3D( TinyVector<int,3>(rcxres,rcyres,rczres),TinyVector<int,3>(wavelet_levelsX,wavelet_levelsY,wavelet_levelsZ),WAVELET3D::WAVE_DB4);

			// Setup Soft Thresholding
			softthresh = THRESHOLD(argc,argv);
	
			// TEMPX
			lranktime=LOWRANKCOIL(argc,argv);
	
	}
				
	// Signal that the recon is ready to perform the requested reconstruction
	prep_done = true;
	
	// Setup Gridding + FFT Structure
	smsEncode smsgrid;
	smsgrid.read_commandline(argc,argv);
	smsgrid.sms_factor = 1;
	smsgrid.precalc_gridding(rczres,rcyres,rcxres,sms_frames,sms_encodes,2,data);
	int rcframes = sms_frames + smsgrid.sms_factor -1;
	cout << "Recon for " << rcframes << " frames " << endl;
	
	Array< Array<complex<float>,3>,2> X = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,rcframes,sms_encodes);
	
	// Compute Forward Transform
	cout << "Create X" <<endl << flush;
	Array< Complex3D,2>R = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,rcframes,sms_encodes);
	Array< Complex3D,2>P = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,rcframes,sms_encodes);
	Array< Complex3D,2> diff_data = Alloc5DContainer< complex<float> >(SmsKx(0).length(firstDim),SmsKx(0).length(secondDim),SmsKx(0).length(thirdDim),sms_frames,sms_encodes);
	
	{ 
	
		tictoc timer;
		
		cout << "Iterate" << endl;
		double error0=0.0;
		for(int iteration =0; iteration< max_iter; iteration++){

			tictoc iteration_timer;
			iteration_timer.tic();
			cout << "\nIteration = " << iteration << endl;

			// Set R to Zero
			cout << "\tZero R " << flush;
			timer.tic();
			for( int t=0; t< rcframes; t++){
				for(int e =0; e< sms_encodes; e++){
					R(t,e) = complex<float>(0.0,0.0);
				}
			}
			cout << "(took " << timer << ")" << endl << flush;
		
									  
			cout << "\tGradient Calculation" << flush;
			for(int coil=0; coil< Ncoils; coil++){
				
				cout << "\t\tZeroing" << flush;
				timer.tic();
				for( int e =0; e < sms_encodes; e++){
					for( int t=0; t< sms_frames; t++){
						diff_data(t,e) = complex<float>(0.0,0.0);
				}}				
				cout << "(took " << timer << ")" << endl << flush;
		
				// Ex
				if(iteration > 0){
					cout << "\t\tBackwards" << flush;
					timer.tic();
					smsgrid.backward( X, smaps(coil), diff_data, SmsKx, SmsKy, SmsKz, SmsKw,SmsZ);
					cout << "(took " << timer << ")" << endl << flush;
				}
				
				//Ex-d
				cout << "\t\tSubtract" << flush;
					cout << "(took " << timer << ")" << endl << flush;
				
				for( int e =0; e < sms_encodes; e++){
					for( int t=0; t< sms_frames; t++){
						diff_data(t,e) -= SmsData(t,e,coil);
				}}
				cout << "(took " << timer << ")" << endl << flush;
				
				
				// E'Ex
				cout << "\t\tForward" << flush;
				timer.tic();
				smsgrid.forward( R, smaps(coil), diff_data, SmsKx, SmsKy, SmsKz, SmsKw,SmsZ);
				cout << "(took " << timer << ")" << endl << flush;
				
			}//Coils
													  
			
			// Set P to Zero
			cout << "\tZero P " << flush;
			timer.tic();
			for( int t=0; t< rcframes; t++){
				for(int e =0; e< sms_encodes; e++){
					P(t,e) = complex<float>(0.0,0.0);
				}
			}
			cout << "(took " << timer << ")" << endl << flush;
				
				
			for(int coil=0; coil< Ncoils; coil++){
				
				cout << "\t\tZeroing"  << flush;
				timer.tic();
				for( int e =0; e < sms_encodes; e++){
					for( int t=0; t< sms_frames; t++){
						diff_data(t,e) = complex<float>(0.0,0.0);
				}}		
				cout << "(took " << timer << ")" << endl << flush;
										  
				// EE'(Ex-d)
				cout << "\t\tBackwards"  << flush;
				timer.tic();
				smsgrid.backward( R, smaps(coil), diff_data, SmsKx, SmsKy, SmsKz, SmsKw,SmsZ);
				cout << "(took " << timer << ")" << endl << flush;
				
				
				//E'EE'(Ex-d)
				cout << "\t\tForward"  << flush;
				timer.tic();
				smsgrid.forward( P, smaps(coil), diff_data, SmsKx, SmsKy, SmsKz, SmsKw,SmsZ);
				cout << "(took " << timer << ")" << endl << flush;
				
			}//Coils
			
			for(int t =0; t< rcframes; t++){
				for( int e =0; e < sms_encodes; e++){
					export_slice( R(t,e), "R_mag.dat");
					export_slice( P(t,e), "P_mag.dat");
				}
			}
			
			// Now get the scale factors
			complex<double>scale_RhP(0.0,0.0);
			complex<double>scale_RhR(0.0,0.0);
			
			
			cout << "Calc Scales" << flush;
			timer.tic();
			for( int t=0; t< rcframes; t++){
				for(int e =0; e< sms_encodes; e++){
				
					P(t,e)*= conj(R(t,e));
					for( Array<complex<float>,3>::iterator riter=P(t,e).begin(); riter != P(t,e).end(); riter++){
				   		scale_RhP += complex< double>( real(*riter),imag(*riter)); 
					}
				
					scale_RhR += complex<double>( ArrayEnergy( R(t,e)), 0.0);
				}
			}
			cout << "(took " << timer << ")" << endl << flush;
				
			
			if(iteration ==0){
				error0 = abs(scale_RhR);
			}
			cout << "Error = " << abs(scale_RhR)/error0 << endl;
			
			// Step in direction
			complex<double>scale_double = (scale_RhR/scale_RhP);
			complex<float>scale( real(scale_double),imag(scale_double) );
			
			cout << "RhR = " << scale_RhR;
			cout << "RhP = " << scale_RhP;
			cout << "Scale = " << scale << endl << flush;
			for( int t=0; t< rcframes; t++){
				for(int e =0; e< sms_encodes; e++){
					R(t,e) *= scale;
					X(t,e) -= R(t,e);
			}}
			cout << "Took " << iteration_timer << " s " << endl;
												  
			// Export X slice
			for(int t =0; t< rcframes; t++){
				for( int e =0; e < sms_encodes; e++){
				export_slice( X(t,e), "X_mag.dat");
			}}
			
			if(softthresh.getThresholdMethod() != TH_NONE){
				L1_threshold(X);
			}
			
			// Export X slice
			for(int t =0; t< rcframes; t++){
				for( int e =0; e < sms_encodes; e++){
				export_slice( X(t,e), "X_mag.dat");
			}}
			
			
		}// Iteration			

	}
	
	return(X);

}


Array< Array<complex<float>,3 >,1 > RECON::reconstruct_one_frame( MRI_DATA& data, int frame_number){
	Array< Array<complex<float>,3 >,2 >XX = full_recon( data, Range(frame_number,frame_number),Range(0,0),false);
	Array< Array<complex<float>,3 >,1 >X = XX(0,Range::all());
	return(X);
}

Array< Array<complex<float>,3 >,2 > RECON::reconstruct_all_frames( MRI_DATA& data){
	Array< Array<complex<float>,3 >,2 >XX = full_recon( data, Range(rc_frame_start,rcframes-1),Range(0,rcframes-1-rc_frame_start),false);
	return(XX);
}

Array< Array<complex<float>,3 >,1 > RECON::reconstruct_composite( MRI_DATA& data){
	Array< Array<complex<float>,3 >,2 >XX = full_recon( data, Range(0,0),Range(0,0),true);
	Array< Array<complex<float>,3 >,1 >X = XX(0,Range::all());
	return(X);
}


void RECON::dcf_calc( MRI_DATA& data){

	
	// Setup Gridding + FFT Structure
	DCFgridFFT dcf_gridding;
	dcf_gridding.kernel_type = DCFgridFFT::POLY_KERNEL;
	dcf_gridding.dwinX = dcf_dwin;
	dcf_gridding.dwinY = dcf_dwin;
	dcf_gridding.dwinZ = dcf_dwin;
	dcf_gridding.grid_x = dcf_overgrid;
	dcf_gridding.grid_y = dcf_overgrid;
	dcf_gridding.grid_z = dcf_overgrid;
	dcf_gridding.grid_in_x = 1;
	dcf_gridding.grid_in_y = 1;
	dcf_gridding.grid_in_z = 1;
	dcf_gridding.precalc_kernel();
	dcf_gridding.acc = dcf_acc;
	
	//dcf_gridding.precalc_gridding(rczres+16,rcyres+16,rcxres+16,data.trajectory_dims,data.trajectory_type);
  
	 // Weighting Array for Time coding
	Array<float,3> X(dcf_overgrid*rczres+16,dcf_overgrid*rcyres+16,dcf_overgrid*rcxres+16,ColumnMajorArray<3>());
	
	for(int e=0; e< rcencodes; e++){
		
		Array< float, 3 >Kweight(data.kx(e).shape(),ColumnMajorArray<3>());
		Array< float, 3 >Kweight2(data.kx(e).shape(),ColumnMajorArray<3>());
	
		data.kx(e)*= dcf_scale;
		data.ky(e)*= dcf_scale;
		data.kz(e)*= dcf_scale;
		
		Kweight = 1;
		for(int iter=0; iter < dcf_iter; iter++){
			cout << "Iteration = " << iter << flush;			 
						
			X=0;
			Kweight2 =0.0;
			dcf_gridding.forward(X,Kweight,data.kx(e),data.ky(e),data.kz(e)); // Grid data to K-Space
			dcf_gridding.backward(X,Kweight2,data.kx(e),data.ky(e),data.kz(e)); // Degrid
			
			float kw_sum = 0.0;
			Array< float,3>::iterator kw_iter=Kweight.begin();
			Array< float,3>::iterator kw2_iter=Kweight2.begin();
			for(  ;(kw_iter!=Kweight.end()) && (kw2_iter!=Kweight2.end()); kw_iter++,kw2_iter++){
				if( abs(*kw2_iter) < (1e-3) ){
					*kw_iter = 0.0; 
				}else{
					*kw_iter = *kw_iter / * kw2_iter;
					kw_sum += abs(*kw_iter);
				}
			} 
						
			cout << " Sum = " << kw_sum << endl;
		}
		
		dcf_gridding.scale_kw( Kweight, data.kx(e),data.ky(e),data.kz(e));
		data.kw(e) = abs(Kweight);
				
		data.kx(e)/= dcf_scale;
		data.ky(e)/= dcf_scale;
		data.kz(e)/= dcf_scale;
		
		
	}
	ArrayWrite(data.kw(0),"Kweight_DCF.dat");

}


Array< Array<complex<float>,3 >,2 >RECON::full_recon( MRI_DATA& data, Range times, Range times_store, bool composite){
	
	cout << "Starting Recon" << endl << flush;
			
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
	cout << "Alloc Container for Solution ( " << rcxres << "," << rcyres << "," << rczres << ") x (" << Nt << "," << rcencodes << ")" << endl << flush;
	Array< Array<complex<float>,3>,2>X = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
	
	// Weighting Array for Time coding
	Array< float, 3 >TimeWeight;
	TimeWeight.setStorage( ColumnMajorArray<3>() );

	cout << "Full recon for " << rcencodes << " encodes, " << Nt << "frames " << endl << flush;
	
	switch(recon_type){
		default:
		case(SOS):
		case(PILS):{
			
			
			for(int e=0; e< rcencodes; e++){
				for(int t=0; t< Nt; t++){
					int act_t   = times(t);
					int store_t = times_store(t);
					int act_e = ( pregate_data_flag) ? ( e*Nt + act_t ) : ( e);
									 
					cout << "Recon Encode" << e << " Frame " << t << endl;
					T.tic();
					
					// Temporal weighting 
					if(pregate_data_flag){
					 	TimeWeight.reference( data.kw(act_e) );
					}else{
					 	TimeWeight.resize(data.kw(act_e).shape());
						TimeWeight = data.kw(act_e);
					 	gate.weight_data( TimeWeight, act_e, data.kx(act_e),data.ky(act_e),data.kz(act_e),act_t,GATING::NON_ITERATIVE,frame_type);
   					}
					
					// cout << "\tForward Gridding Coil ";
					for(int coil=0; coil< data.Num_Coils; coil++){
						// K-space to Image
						if(recon_type==PILS){
							gridding.forward(X(store_t,e),smaps(coil),data.kdata(act_e,coil),data.kx(act_e),data.ky(act_e),data.kz(act_e),TimeWeight);
						}else{
							gridding.forward_sos(X(store_t,e),data.kdata(act_e,coil),data.kx(act_e),data.ky(act_e),data.kz(act_e),TimeWeight);
						}
						
					}
					cout << "Took " << T << endl;	
					
					// Take Square Root for SOS	
					if(recon_type==SOS){
				 		X(store_t,e) = csqrt(X(store_t,e));
					}
					
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
							  
							  // Storage for (Ex-d)
							  Complex3D diff_data(data.kx(e).shape(),ColumnMajorArray<3>());
					
	
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
								  
								  TimeWeight.resize(kwE.shape());
								  TimeWeight = kwE;
								  gate.weight_data( TimeWeight,e, kxE, kyE,kzE,act_t,GATING::ITERATIVE, frame_type);
   							 	  
								  for(int coil=0; coil< data.Num_Coils; coil++){
										  
									  T.tic();
									  
									  // Alloc Data
									  Array< complex<float>,3 >diff_data( data.kx(e).shape(),ColumnMajorArray<3>());
								      diff_data=0;
								  
									  // Ex
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
						  	Array< float,2>Xslice( rcxres,rcyres,ColumnMajorArray<2>());
							
							Xslice =0;
							for(int coil=0; coil< data.Num_Coils; coil++){
								Array< complex<float>,2> Xtemp  = XX(0,0,coil)(all,all,RR(0,0,0).length(2)/2);
								Xslice += norm( Xtemp );
						  		ArrayWriteMagAppend( Xtemp,"X_coils.dat");
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
						  lrankcoil.update_threshold( XX,2, iteration);
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
				      
					  
					  {
					  	HDF5 clearRAW("ClearRaw.h5","w");
						
						for(int t =0; t< XX.length(firstDim); t++){
						for(int e =0; e < XX.length(secondDim); e++){
						for(int coil=0; coil< data.Num_Coils; coil++){
							char name[1024];
							sprintf(name,"IM_T%03d_E%03d_C%03d",t,e,coil);
							clearRAW.AddH5Array( "IMAGES",name,XX(t,e,coil));	
						
						}}}
					  }
					  
					  lrankcoil.combine(XX,X);

		}break;
		
		case(ADMM):{

				 // ------------------------------------
				 // Alternating Direction of Multipliers
				 //   0.5||Ex-d||2  + gamma*||y||1 + rho*||x-y||2 + L'*(x-y)  
				
				 // CG Structures  	
				 Array< Array< complex<float>,3>, 2>R = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
				 Array< Array< complex<float>,3>, 2>P = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
				 Array< Array< complex<float>,3>, 2>LHS = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
				 
				 // Augmented Lagrangian Structures
				 Array< Array< complex<float>,3>, 2>Y = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
				 Array< Array< complex<float>,3>, 2>L = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
				 
				
				 float cg_thresh = 1e-9;
				 
				 // Initialize
				 //X = 0;
				 //Y = 0;
				 //L = 0;
				 for(int e=0; e< rcencodes; e++){
					for(int t=0; t< Nt; t++){
						X(t,e) = complex<float>(0.0,0.0);
						Y(t,e) = complex<float>(0.0,0.0);
						L(t,e) = complex<float>(0.0,0.0);
						LHS(t,e) = complex<float>(0.0,0.0);
				 }}		
							
				 
				 // Iterate
				 for(int admm_iteration =0; admm_iteration< admm_max_iter; admm_iteration++){
				 				 				 
				 		
						//   
						//  Solve Subproblem  1  
						//  	0.5||Ex-d||2  + gamma*||y||1 + rho*||x-y||2 + L'*(x-y)
						//      where y and L are fixed paramaters. This reduces to minimizing
						// 		0.5||Ex-d||2  + rho*||x-y||2 + L'*(x-y)	
						
						// (E'E + rho*I)*x = E'd + rho*y - L; 
							
						// First Calculate E'd
				 		cout << "ADMM Inner :: LHS Calculation" << endl;
						for(int e=0; e< rcencodes; e++){
							for(int t=0; t< Nt; t++){
								
								R(t,e) = complex<float>(0.0,0.0);
								
								int act_t = times(t);
								int store_t = times_store(t);
																
								// Temporal weighting 
								if(pregate_data_flag){
					 				TimeWeight.reference( data.kw(e) );
								}else if(reset_dens){
								    TimeWeight.resize(data.kw(e).shape());
									TimeWeight = 1.0;								  
									gate.weight_data( TimeWeight, e, data.kx(e),data.ky(e),data.kz(e),act_t,GATING::ITERATIVE,frame_type);
   								}else{
					 				TimeWeight.resize(data.kw(e).shape());
									TimeWeight = data.kw(e);
					 				gate.weight_data( TimeWeight, e, data.kx(e),data.ky(e),data.kz(e),act_t,GATING::ITERATIVE,frame_type);
   								}
						 								
								// Images
								for(int coil=0; coil< data.Num_Coils; coil++){
									 
									Array< complex<float>,3 >diff_data( data.kx(e).shape(),ColumnMajorArray<3>());
									diff_data = complex<float>(0.0,0.0);
									
									// Ex
									gridding.backward(X(store_t,e),smaps(coil),diff_data,data.kx(e),data.ky(e),data.kz(e),TimeWeight);

									// d - Ex
									diff_data = data.kdata(e,coil) - diff_data;

									// E'd - E'Ex
									gridding.forward(R(store_t,e),smaps(coil),diff_data,data.kx(e),data.ky(e),data.kz(e),TimeWeight);
										  
								 }//Coils
								 
								 // E'd - ( E'Ex + rho*Ix) 
								 R(store_t,e) -= admm_rho*X(store_t,e);
								 
								 
								 // RHS Values ( E'd + rho*Y - L ) - (E' + rho)*X 
								 R(store_t,e) += admm_rho*Y(store_t,e);
								 R(store_t,e) -= L(store_t,e);
								 
								 // Initiialize P
								 P(store_t,e) = R(store_t,e);
							}
				 		}
						
																				 
				 		// Now Iterate
						cout << "ADMM Inner :: Iterate " << endl;
						for(int cg_iteration =0; cg_iteration< max_iter; cg_iteration++){

							tictoc iteration_timer;
							iteration_timer.tic();
					  
							// E'Ex
							for(int e=0; e< rcencodes; e++){
								for(int t=0; t< Nt; t++){
									int act_t = times(t);
									//int store_t = times_store(t);
							  
									LHS(t,e ) = complex<float>(0.0,0.0);
									T.tic();
 
									// Temporal weighting 
									if(pregate_data_flag){
					 					TimeWeight.reference( data.kw(e) );
									}else if(reset_dens){
								        TimeWeight.resize(data.kw(e).shape());
										TimeWeight = 1.0;								  
										gate.weight_data( TimeWeight, e, data.kx(e),data.ky(e),data.kz(e),act_t,GATING::ITERATIVE,frame_type);
   									}else{
					 					TimeWeight.resize(data.kw(e).shape());
										TimeWeight = data.kw(e);
					 					gate.weight_data( TimeWeight, e, data.kx(e),data.ky(e),data.kz(e),act_t,GATING::ITERATIVE,frame_type);
   									}
							  
							  		for(int coil=0; coil< data.Num_Coils; coil++){
										  
										// Storage for (Ex-d) - dynamic for variable size
			  		  					Array< complex<float>,3 >diff_data( data.kx(e).shape(),ColumnMajorArray<3>());
										diff_data=complex<float>(0.0,0.0);

										// E'Ex
										gridding.backward( P(t,e),smaps(coil),diff_data,data.kx(e),data.ky(e),data.kz(e),TimeWeight);
										gridding.forward( LHS(t,e),smaps(coil),diff_data,data.kx(e),data.ky(e),data.kz(e),TimeWeight);
									}//Coils
									
									// Values
									LHS(t,e) += admm_rho*P(t,e);
							
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
										
							cout << "       Sum R'R = " << sum_R_R << endl;
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
							
														
							{
								Array<complex<float>,2>Xslice=X(0,0)(all,all,X(0,0).length(2)/2);
						  		ArrayWriteMagAppend(Xslice,"X_mag.dat");
							}
							
							if( abs(sum_R_R) < cg_thresh){
								break;
							}
							
						}// Inner CG Iterations
							
							
							
						//   
						//  Solve Subproblem  2  
						//  	0.5||Ex-d||2  + gamma*||y||1 + rho*||x-y||2 + L'*(x-y)
						//      where x and L are fixed paramaters. This reduces to
						//		gamma*||y||1 + rho*||x-y||2 + L'*(x-y)  				
						
						
						// Take step size
						float alpha = 0;
						for(int e=0; e< rcencodes; e++){
					 		for(int t=0; t< Nt; t++){
								for(int k=0; k < rczres; k++){
									for(int j=0; j < rcyres; j++){
										for(int i=0; i < rcxres; i++){
											Y(t,e)(i,j,k) = alpha*X(t,e)(i,j,k) + (1-alpha)*Y(t,e)(i,j,k) + L(t,e)(i,j,k)/admm_rho;
											
											// Also update dual
											L(t,e)(i,j,k);
											
								}}}
						}}
						
						
						
						{
							Array<complex<float>,2>Xslice=Y(0,0)(all,all,X(0,0).length(2)/2);
							ArrayWriteMagAppend(Xslice,"Y_mag.dat");
						}
						
						
						
						if(admm_iteration==0){
							admm_gamma *= max(abs(X(0,0)));
						}
						softthresh.threshold_type = TH_FIXED;
						softthresh.thresh = admm_gamma / admm_rho;
						L1_threshold(Y);
						
						{
							Array<complex<float>,2>Xslice=Y(0,0)(all,all,X(0,0).length(2)/2);
							ArrayWriteMagAppend(Xslice,"Y_mag.dat");
						}
								
						// 
						// Update the lagrangian multiplier
						// 
						for(int e=0; e< rcencodes; e++){
					 		for(int t=0; t< Nt; t++){
								for(int k=0; k < rczres; k++){
									for(int j=0; j < rcyres; j++){
										for(int i=0; i < rcxres; i++){
												L(t,e)(i,j,k) = L(t,e)(i,j,k) +  admm_rho*( X(t,e)(i,j,k) - Y(t,e)(i,j,k));
								}}}
						}}	
						
						{
							Array<complex<float>,2>Xslice=L(0,0)(all,all,X(0,0).length(2)/2);
							ArrayWriteMagAppend(Xslice,"L_mag.dat");
						}
						
			}//iteration	

		}break;

		case(CG):{

				 // ------------------------------------
				 // Conjugate Gradient Recon 
				 //   -Uses much more memory than gradient descent but is faster (both convergence + per iteration)
				 // ------------------------------------
				
				 // Structures  	
				 Array< Array< complex<float>,3>, 2>R = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
				 Array< Array< complex<float>,3>, 2>P = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
				 Array< Array< complex<float>,3>, 2>LHS = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
				 				 
				 				 			 
				 // First Calculate E'd
				 cout << "LHS Calculation" << endl;
				 for(int e=0; e< rcencodes; e++){
					for(int t=0; t< Nt; t++){
						int act_t = times(t);
						int store_t = times_store(t);
						int act_e = ( pregate_data_flag) ? ( e*Nt + act_t ) : ( e);
								
						// Temporal weighting 
						if(pregate_data_flag){
					 		TimeWeight.reference( data.kw(act_e) );
						}else if(reset_dens){
						    TimeWeight.resize(data.kw(act_e).shape());
							TimeWeight = 1.0;								  
							gate.weight_data( TimeWeight, act_e, data.kx(act_e),data.ky(act_e),data.kz(act_e),act_t,GATING::ITERATIVE,frame_type);
   						}else{
					 		TimeWeight.resize(data.kw(act_e).shape());
							TimeWeight = data.kw(act_e);
					 		gate.weight_data( TimeWeight, act_e, data.kx(act_e),data.ky(act_e),data.kz(act_e),act_t,GATING::ITERATIVE,frame_type);
   						}
						
						// E'd
						for(int coil=0; coil< data.Num_Coils; coil++){
							gridding.forward( R(store_t,e),smaps(coil),data.kdata(act_e,coil),data.kx(act_e),data.ky(act_e),data.kz(e),TimeWeight);
						} 
						
				 	}
				 }
				 
				 // Initiialize P
				 P  = R;	
														 
				 // Now Iterate
				 cout << "Iterate" << endl;
				 // double error0=0.0;
				 for(int iteration =0; iteration< max_iter; iteration++){

					  tictoc iteration_timer;
					  iteration_timer.tic();
					  
					  // E'Ex
					  for(int e=0; e< rcencodes; e++){
						  for(int t=0; t< Nt; t++){
							  int act_t = times(t);
							  int store_t = times_store(t);
							  int act_e = ( pregate_data_flag) ? ( e*Nt + act_t ) : ( e);
						
							  LHS(t,e ) = complex<float>(0.0,0.0);
							  T.tic();
 
							  // Temporal weighting 
							  if(pregate_data_flag){
					 			TimeWeight.reference( data.kw(act_e) );
							  }else if(reset_dens){	
								TimeWeight.resize(data.kw(act_e).shape());
								TimeWeight = 1.0;								  
								gate.weight_data( TimeWeight, act_e, data.kx(act_e),data.ky(act_e),data.kz(act_e),act_t,GATING::ITERATIVE,frame_type);
   							  }else{
					 			TimeWeight.resize(data.kw(act_e).shape());
								TimeWeight = data.kw(act_e);
					 			gate.weight_data( TimeWeight, act_e, data.kx(act_e),data.ky(act_e),data.kz(act_e),act_t,GATING::ITERATIVE,frame_type);
   							  }
							  
							  for(int coil=0; coil< data.Num_Coils; coil++){
										  
								  // Storage for (Ex-d) - dynamic for variable size
			  		  			  Array< complex<float>,3 >diff_data( data.kx(act_e).shape(),ColumnMajorArray<3>());
								  diff_data=0;

								  // E'Ex
								  gridding.backward( P(store_t,e),smaps(coil),diff_data,data.kx(act_e),data.ky(act_e),data.kz(act_e),TimeWeight);
								  gridding.forward( LHS(store_t,e),smaps(coil),diff_data,data.kx(act_e),data.ky(act_e),data.kz(act_e),TimeWeight);
							  }//Coils
							  
							  // L2 R'R 
							  if(iteration> 0){
							  	
							  }
							  
						}// t
					}// e
					
					
					// Regularization
					if( (iteration==0) && (l2reg.lambda>0) ){
						float sum_P_P = 0.0; 
						float sum_L_L = 0.0; 
						for(int e=0; e< LHS.length(secondDim); e++){
						  for(int t=0; t< LHS.length(firstDim); t++){
							sum_P_P += sum(norm(P(t,e)));
							sum_L_L += sum(norm(LHS(t,e)));
						}}
						
						l2reg.reg_scale = l2reg.lambda*sqrt(sum_L_L/sum_P_P);
						cout << "L2 Scale = " << l2reg.reg_scale << endl;
					}
					
					for(int e=0; e< LHS.length(firstDim); e++){
						  for(int t=0; t< LHS.length(secondDim); t++){
							l2reg.regularize(LHS(t,e),P(t,e));
						}
					}
											  
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
						  
						  
						  if( X.numElements() > 1){
						  	int count=0;
							for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
								
								Array<complex<float>,2>Xf=(*miter)(all,all,X(0,0).length(2)/2);
								if(count==0){
									ArrayWriteMag(Xf,"X_frames.dat");
									ArrayWritePhase(Xf,"X_frames.dat.phase");
								}else{
									ArrayWriteMagAppend(Xf,"X_frames.dat");
									ArrayWritePhaseAppend(Xf,"X_frames.dat.phase");
								}
								count++;
							}
						  }
						  
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
						
					  // Residue 	
					  Array< Array< complex<float>,3>, 2>R = Alloc5DContainer< complex<float> >(rcxres,rcyres,rczres,Nt,rcencodes);
	
					  // Temp variable for E'ER 
					  Array< complex<float>,3 >P(rcxres,rcyres,rczres,ColumnMajorArray<3>());

					  // Class for gradient descent step size
					  complex<float>step_size; 
					    
					  cout << "Iterate" << endl;
					  double error0=0.0;
					  for(int iteration =0; iteration< max_iter; iteration++){

						  tictoc iteration_timer;
						  iteration_timer.tic();
						  cout << "\nIteration = " << iteration << endl;

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
							  	int act_e = ( pregate_data_flag) ? ( e*Nt + act_t ) : ( e);
									
								// Get Sub-Arrays for Encoding
							  	Array< float,3 >kxE = data.kx(act_e); 
							  	Array< float,3 >kyE = data.ky(act_e); 
							  	Array< float,3 >kzE = data.kz(act_e);
							  	Array< float,3 >kwE = data.kw(act_e);
							  	
								T.tic();
																  								  
								// Temporal weighting 
								if(pregate_data_flag){
					 				TimeWeight.reference( kwE );
								}else if(reset_dens){
									TimeWeight.resize(kwE.shape());
									TimeWeight = 1.0;
									Array< float,3>::iterator titer=TimeWeight.begin();
									Array< float,3>::iterator kiter=kwE.begin();
									for( ; (titer!=TimeWeight.end()) && (kiter!=kwE.end()); titer++,kiter++){
										if( (*kiter) < 0.1 ){
											*titer = *kiter;
										}else{
											*titer = 0.1;
										}	
									} 
									gate.weight_data( TimeWeight, e,kxE,kyE,kzE,act_t,GATING::ITERATIVE,frame_type);
   								}else{
					 				TimeWeight.resize(kwE.shape());
									TimeWeight = kwE;
					 				gate.weight_data( TimeWeight, e, kxE,kwE,kzE,act_t,GATING::ITERATIVE,frame_type);
   								}
								 								 								  
								// Differences (Ex-d)
								Array< complex<float>,3 >diff_data( kxE.shape(),ColumnMajorArray<3>());
								
								for(int coil=0; coil< data.Num_Coils; coil++){
									
									// Ex - d
									gridding.backward_residual(X(store_t,e),smaps(coil),diff_data,kxE,kyE,kzE,TimeWeight,data.kdata(act_e,coil));

									//E'(Ex-d)
									gridding.forward(R(store_t,e),smaps(coil),diff_data,kxE,kyE,kzE,TimeWeight);
										  
								}//Coils
													  
								  // L2 
								if(iteration > 0){
								  	l2reg.regularize(R(store_t,e),X(store_t,e) );
								}
								   									
								// Now Get Scale factor (for Cauchy-Step Size)
								if( (iteration % this->step_update_frequency) == 0){
								 
									P=0;
									for(int coil=0; coil< data.Num_Coils; coil++){
										  
										// EE'(Ex-d)
										gridding.backward(R(store_t,e),smaps(coil), diff_data,kxE,kyE,kzE,TimeWeight);

										//E'EE'(Ex-d)
										gridding.forward(            P,smaps(coil),diff_data,kxE,kyE,kzE,TimeWeight);
									}//Coils
								  
								  
									// TV of Image
									if(iteration > 0){
										l2reg.regularize(P,R(store_t,e));
									}
								 
									P*=conj(R(store_t,e));
								  
									for( Array<complex<float>,3>::iterator riter=P.begin(); riter != P.end(); riter++){
										scale_RhP += complex< double>( real(*riter),imag(*riter)); 
									}
								}
								cout << "\r" << e << "," << t << "took " << T << "s" << flush;
							  }//Time
						  }//Encode
						  
						  // Get Scaling Factor R'P / R'R 
						  cout << endl << "Calc residue" << endl << flush;
						  complex<double>scale_RhR = 0.0;
						  for( Array< Array<complex<float>,3>,2>::iterator riter =R.begin(); riter != R.end(); riter++){
						  		scale_RhR += complex<double>( ArrayEnergy( *riter ), 0.0);
						  }
						  
						  // Error check
						  if(iteration==0){
							  error0 = abs(scale_RhR);
						  }
						  cout << "Residue Energy =" << scale_RhR << " ( " << (abs(scale_RhR)/error0) << " % )  " << endl;
						  cout << "RhP = " << scale_RhP << endl;
						  
						  // Export R
						  export_slice( R(0,0),"R.dat");
						  
						  // Step in direction
						  if( (iteration % this->step_update_frequency) == 0){
						  	complex<double>scale_double = (scale_RhR/scale_RhP);
						  	step_size = complex<float>( real(scale_double),imag(scale_double) );
						  } 
						  
						  cout << "Scale = " << step_size << endl << flush;
						  for(int e=0; e< rcencodes; e++){
							  for(int t=0; t< Nt; t++){
						  		R(t,e) *= step_size;
						  		X(t,e) -= R(t,e);
							}
						  }
						  cout << "Took " << iteration_timer << " s " << endl;
						
						
						   // Error check
						  if(iteration==0){
							  l2reg.set_scale( (float)abs(scale_RhR), X);
							  cout << "L2 Scale = " << l2reg.reg_scale << endl;
						  }
						  						  
						  // Export X slice
						  export_slice( X(0,0), "X_mag.dat");
						  
						  if(lranktime.clear_alpha_time > 0.0){
						  	lranktime.update_threshold(X,0,iteration);
						  	lranktime.thresh(X,0);
						  	export_slice( X(0,0), "X_mag.dat");
						  }
						  
						  if(lranktime.clear_alpha_encode > 0.0){
						  	lranktime.update_threshold(X,1,iteration);
						  	lranktime.thresh(X,1);
						  	export_slice( X(0,0), "X_mag.dat");
						  }
						  
						  if( X.numElements() > 1){
						  	int count=0;
							for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
								
								Array<complex<float>,2>Xf=(*miter)(all,all,X(0,0).length(2)/2);
								if(count==0){
									ArrayWriteMag(Xf,"X_frames.dat");
									ArrayWritePhase(Xf,"X_frames.dat.phase");
								}else{
									ArrayWriteMagAppend(Xf,"X_frames.dat");
									ArrayWritePhaseAppend(Xf,"X_frames.dat.phase");
								}
								count++;
							}
						  }
							
						  // ------------------------------------
						  // Soft thresholding operation (need to add transform control)
						  // ------------------------------------
						  if(softthresh.getThresholdMethod() != TH_NONE){
							  L1_threshold(X);
						  }
  						  export_slice( X(0,0), "X_mag.dat");
						  
					  }// Iteration			

				  }break;

	}//Recon Type


	cout << "Recon was completed successfully " << endl;	
	return(X);
}


double RECON::kspace_residual( MRI_DATA& data){
	
	// Matlab like timer (openmp code base)
	tictoc T; 
	
	// Setup Gridding + FFT Structure
	gridFFT Kgridding;
	Kgridding.overgrid = 1.0;
	Kgridding.kernel_type = TRIANGLE_KERNEL;
	Kgridding.dwinX = 1.0;
	Kgridding.dwinY = 1.0;
	Kgridding.dwinZ = 1.0;
	Kgridding.precalc_gridding(32,32,32,data);

	// Class for gradient descent step size
	complex<float>step_size; 
					    
	double residual = 0.0;
	for(int e=0; e< data.kx.length(firstDim); e++){
	  					  	
						
		// Get Sub-Arrays for Encoding
		Array< float,3 >kxE = data.kx(e); 
		Array< float,3 >kyE = data.ky(e); 
		Array< float,3 >kzE = data.kz(e);
		Array< float,3 >kwE = data.kw(e);
								 								 								  
		// Differences (Ex-d)
		Array< complex<float>,3 >diff_data( kxE.shape(),ColumnMajorArray<3>());
		
		for(int coil=0; coil< data.Num_Coils; coil++){
			
			// Grid data
			T.tic();
			Kgridding.k3d_grid=0; // Zero K-Space
			Kgridding.chop_grid_forward(data.kdata(e,coil),kxE,kyE,kzE,kwE); // Grid data to K-Space
			//cout << "Forward took " << T << endl;
			
			// Inverse
			T.tic(); 
			Array<complex<float>,3> temp;
			Kgridding.chop_grid_backward(diff_data,kxE,kyE,kzE,kwE,temp,false);
			//cout << "Backward took " << T << endl;
			
			T.tic();
			
			// Residual
			complex<double> DhD(0.0,0.0);
			complex<double> DhR(0.,0.);
	
			#pragma omp parallel for
			for(int j=0; j < diff_data.length(secondDim); j++){
			
				complex<double> DhD_temp(0.0,0.0);
				complex<double> DhR_temp(0.0,0.0);
			
				for(int k=0; k < diff_data.length(thirdDim); k++){
				for(int i=0; i < diff_data.length(firstDim); i++){
				
					// Energy in D'*D			
					DhD_temp +=  data.kdata(e,coil)(i,j,k)*conj( data.kdata(e,coil)(i,j,k) );
			
					// Energy in D'*R
					DhR_temp +=  conj( data.kdata(e,coil)(i,j,k))*diff_data(i,j,k);
				}}
			
				// Prevent Race conditions in multi-threaded
				#pragma omp critical
				{
					DhD +=DhD_temp;
					DhR +=DhR_temp;
				}
			}			
			
			// Scale D*'D/(D'*R)		
			complex<double> scale = DhD/DhR;		
			complex<float> scaleF( real(scale),imag(scale));		
			
			#pragma omp parallel for
			for(int j=0; j < diff_data.length(secondDim); j++){
				
				double residual_temp=0.0;
				for(int k=0; k < diff_data.length(thirdDim); k++){
				for(int i=0; i < diff_data.length(firstDim); i++){
				
					// Energy in D'*D			
					residual_temp +=  norm( data.kdata(e,coil)(i,j,k)-scaleF*diff_data(i,j,k));
				}}
			
				// Prevent Race conditions in multi-threaded
				#pragma omp critical
				{
					residual += residual_temp;
				}
			}			
			//cout << "Residual took " << T << endl;
			
		}//Coils
	}// Encode
													  
	
	return(residual);
}


void RECON::export_slice( Array< complex<float>,3> &temp, const char *fname){
	
	int slice = (int)(  temp.length(2)/2 );
	Array<complex<float>,2>Xslice=temp(Range::all(),Range::all(),slice);
	ArrayWriteMagAppend(Xslice,fname);
}



void RECON::transform_in_encode( Array< Array< complex<float>,3>, 2>&X, TransformDirection direction){

	
	switch(cs_encode_transform){
		case(DFT):{	
			cout << "DFT in Encode" << endl;
			if(direction==FORWARD){
				TRANSFORMS::fft_e(X); 
			}else{
				TRANSFORMS::ifft_e(X); 
			}
		}break;
								
		case(DIFF):{
			cout << "DIFF in Encode" << endl;
			if(direction==FORWARD){
				TRANSFORMS::ediff(X); 
			}else{
				TRANSFORMS::inv_ediff(X); 
			}
		}break;
		
		case(PCA):{
			cout << "PCA in Encode" << endl;
			if(direction==FORWARD){
				sparse_transform.eigen(X,1,TRANSFORMS::FORWARD); 
			}else{
				sparse_transform.eigen(X,1,TRANSFORMS::BACKWARD); 
			}
		}break;
								
		case(WAVELET):{
			cout << "WAVELET in Encode" << rcencodes << endl;
			if(direction==FORWARD){
				TRANSFORMS::ewave(X); 
			}else{
				TRANSFORMS::inv_ewave(X); 
			}
		}break;
		
		default:{
		}break;
	}
}

void RECON::transform_in_time( Array< Array< complex<float>,3>, 2>&X, TransformDirection direction){

	
	switch(cs_temporal_transform){
		case(DFT):{	
			cout << "DFT in Time" << endl;
			if(direction==FORWARD){
				TRANSFORMS::fft_t(X); 
			}else{
				TRANSFORMS::ifft_t(X); 
			}
		}break;
								
		case(DIFF):{
			cout << "DIFF in Time" << endl;
			if(direction==FORWARD){
				TRANSFORMS::tdiff(X); 
			}else{
				TRANSFORMS::inv_tdiff(X); 
			}
		}break;
								
		case(WAVELET):{
			cout << "WAVELET in Time" << endl;
			if(direction==FORWARD){
				TRANSFORMS::twave(X); 
			}else{
				TRANSFORMS::inv_twave(X); 
			}
		}break;
		
		case(PCA):{
			cout << "PCA in Time" << endl;
			if(direction==FORWARD){
				sparse_transform.eigen(X,0,TRANSFORMS::FORWARD); 
			}else{
				sparse_transform.eigen(X,0,TRANSFORMS::BACKWARD); 
			}
		}break;
		
		
		case(COMPOSITE_DIFF):{
				cout << "Composite Diff" << endl;
				for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
						if(direction==FORWARD){
							*miter -= composite_image;
						}else{
							*miter += composite_image;
						}
				}
		}break;
								
		default:{
		}break;
	}
}



void RECON::L1_threshold( Array< Array< complex<float>,3>, 2>&X){
	
	// Use a matrix to rotate into low resolution phase
	Array<complex<float>,3>Phase(rcxres,rcyres,rczres, ColumnMajorArray<3>());
	if(phase_rotation){
		cout << "Rotating to zero phase plane" << endl;
		Phase = complex<float>(0.0,0.0);
		for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
			Phase += (*miter);
		}
		
		// FFT to K-Space
		ifft(Phase);		
		
		// Set Values to zero
		int cx = rcxres /2;
		int cy = rcyres /2;
		int cz = rczres /2;
		
		// Calc windows
		Array< float,1> Wx(Phase.length(firstDim));
		for( int i=0; i< Phase.length(firstDim); i++){
			float rx = abs( i - cx);
			Wx(i) = 1./( 1.0 + exp( (rx - phase_rotation_sX ) / 6));
		}
		
		
		Array< float,1> Wy(Phase.length(secondDim));
		for( int i=0; i< Phase.length(secondDim); i++){
			float ry = abs( i - cy);
			Wy(i) = 1./( 1.0 + exp( (ry - phase_rotation_sY ) / 6));
		}
				
		Array< float,1> Wz(Phase.length(thirdDim));
		for( int i=0; i< Phase.length(thirdDim); i++){
			float rz = abs( i - cz);
			Wz(i) = 1./( 1.0 + exp( (rz - phase_rotation_sZ ) / 6));
		}
					
		
		for(int k=0; k < Phase.length(thirdDim); k++){
		for(int j=0; j < Phase.length(secondDim); j++){
		for(int i=0; i < Phase.length(firstDim); i++){
			Phase(i,j,k) *= (Wx(i)*Wy(j)*Wz(k));
		}}}
		
		// FFT Back
		fft(Phase);
		
		// Normalize
		for(int k=0; k < Phase.length(thirdDim); k++){
		for(int j=0; j < Phase.length(secondDim); j++){
		for(int i=0; i < Phase.length(firstDim); i++){
			Phase(i,j,k) /= abs( Phase(i,j,k));
		}}}
		
		ArrayWritePhase(Phase,"PhaseCorr.dat");
		
		for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
			*miter *= conj(Phase);
		}
	}
	
	int actual_cycle_spins=cycle_spins;
	if(cs_spatial_transform != WAVELET){
		actual_cycle_spins = 1;
	}
	
	transform_in_encode( X, FORWARD);
	transform_in_time( X , FORWARD);
	
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
	
	transform_in_time( X , BACKWARD);
	transform_in_encode( X, BACKWARD);

	// Use a matrix to rotate into low resolution phase
	if(phase_rotation){
		cout << "Rotating to non-zero phase plane" << endl;
		for( Array< Array<complex<float>,3>,2>::iterator miter=X.begin(); miter!=X.end(); miter++){
			*miter *= Phase;
		}
	}
}

inline float sqr ( float x){ return(x*x);}

void RECON::calc_sensitivity_maps( int argc, char **argv, MRI_DATA& data){

	// ------------------------------------
	//  Get coil sensitivity map ( move into function)
	// ------------------------------------
	
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
					
						//gridding.forward(data.kdata(e,coil),data.kx(e),data.ky(e),data.kz(e),data.kw(e));
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
		
		case(ADMM):
		case(PILS):
		case(IST):
		case(FISTA):
		case(CG):{
			
			// Allocate Storage for Map	and zero	
			cout << "Allocate Sense Maps"  << endl << flush;
			{
				Array< Array< complex<float>, 3>,1> temp =  Alloc4DContainer< complex<float> >(rcxres,rcyres,rczres,data.Num_Coils);
				smaps.reference(temp);
			}

			if(data.Num_Coils ==1){
				smaps(0)=complex<float>(1.0,0.0);
				break;
			}
			 
			cout << "Recon Low Resolution Images"  << endl<< flush; 
			
			// Multiple Encode Smaps
			Array< Array<complex<float>,3>,2> image_store;
			int smap_num_encodes=1;
			if( smap_use_all_encodes){
				cout << "Using all encodes for coil maps" << endl;
				Array< Array<complex<float>,3>,2> temp = Alloc5DContainer< complex<float> >(rcxres, rcyres, rczres, data.Num_Encodings,data.Num_Coils);
				image_store.reference( temp);
				smap_num_encodes = data.Num_Encodings;
			}else if( coil_combine_type == WALSH){	
				cout << "Using one encodes for coil maps" << endl;
				Array< Array<complex<float>,3>,2> temp = Alloc5DContainer< complex<float> >(rcxres, rcyres, rczres, 1,data.Num_Coils);
				image_store.reference( temp);
			}else{
				// Point to Smaps
				image_store.setStorage( ColumnMajorArray<2>() );
				image_store.resize(1,data.Num_Coils);
				for( int coil = 0; coil < data.Num_Coils; coil++){
					image_store(0,coil).reference( smaps(coil) );
				}			
			}
			 
			// Low Pass filtering for Sensitivity Map
			if(coil_combine_type!=ESPIRIT){
				gridding.k_rad = smap_res;
			}
			
			if( smap_use_all_encodes){
			
				for(int e=0; e< smap_num_encodes;e++){
					for(int coil=0; coil< data.Num_Coils; coil++){
						cout << "Coil " << coil << "Encode " << e << endl;	
						// Simple gridding
						gridding.forward( image_store(e,coil),data.kdata(e,coil),data.kx(e), data.ky(e), data.kz(e) ,data.kw(e) );
					
						// Gaussian blur	
						gaussian_blur(image_store(e,coil),extra_blurX,extra_blurY,extra_blurZ); // TEMP						
					}
				}
			}else{
				for(int coil=0; coil< data.Num_Coils; coil++){
					image_store(0,coil) =complex<float>(0.0,0.0);
					
					int used_encodes = smap_nex_encodes ?  data.Num_Encodings : 1;
					for( int e = 0; e < used_encodes; e++){
						cout << "Coil " << coil << "Encode " << e << endl;	
						// Simple gridding
						gridding.forward( image_store(0,coil),data.kdata(e,coil),data.kx(e), data.ky(e), data.kz(e) ,data.kw(e) );
					}
					// Gaussian blur	
					gaussian_blur(image_store(0,coil),extra_blurX,extra_blurY,extra_blurZ); // TEMP	
				}
			}
			gridding.k_rad = 9999;
			
			Array< float , 3 > IC;
			if( intensity_correction ){
				intensity_correct( IC, smaps);
			}
		
            if (this->coil_rejection_flag) {
                
                for (int coil = 0; coil< data.Num_Coils; coil++) {
                    double sumI=0;
                    double sumIX = 0;
                    double sumIY = 0;
                    double sumIZ = 0;
                    double cx = 0.5*(double)smaps(0).length(firstDim) - 0.5;
                    double cy = 0.5*(double)smaps(0).length(secondDim) - 0.5;
                    double cz = 0.5*(double)smaps(0).length(thirdDim) - 0.5;
                    
                    // Get the coil center
                    for (int e = 0; e < image_store.length(firstDim); e++) {
                        for (int k = 0; k<smaps(0).length(thirdDim); k++) {
                            for (int j = 0; j<smaps(0).length(secondDim); j++) {
                                for (int i = 0; i<smaps(0).length(firstDim); i++) {
                                    double val = norm( image_store(e,coil)(i, j, k) );
                                    val *= val;
                                    sumI += val;
                                    sumIX+= (val)*( (double)i - cx);
                                    sumIY+= (val)*( (double)j - cy);
                                    sumIZ+= (val)*( (double)k - cz);
                                }
                            }
                        }
                    }

                    double coil_cx = 2.0*sumIX / sumI / (double)smaps(0).length(firstDim);
                    double coil_cy = 2.0*sumIY / sumI / (double)smaps(0).length(secondDim);
                    double coil_cz = 2.0*sumIZ / sumI / (double)smaps(0).length(thirdDim);

                    std::cout << "Coil = " << coil << "center=" << coil_cx << "," << coil_cy << "," << coil_cz << std::endl;
                    
                    double coil_radius;
                    switch (this->coil_rejection_shape) {
                        default:
                        case(0): {
                            coil_radius = sqrt(coil_cx*coil_cx + coil_cy*coil_cy + coil_cz*coil_cz);
                        }break;

                        case(1): {
                            coil_radius = sqrt(coil_cx*coil_cx + coil_cy*coil_cy );
                        }break;
                    }
                                      
                    if (coil_radius > this->coil_rejection_radius ) {
                        std::cout << "REJECTING COIL " << coil << std::endl;
                        for (int e = 0; e < image_store.length(firstDim); e++) {
                            image_store(e, coil) = complex<float>(0.0,0.0);
                        }
                    }

                }
            }




			// Spirit Code
			switch(coil_combine_type){
		
				case(ESPIRIT):{
 					// Espirit (k-space kernal Eigen method)
					cout << "eSPIRIT Based Maps"  << endl; 
					SPIRIT S;
      	 				S.read_commandline(argc,argv);
      					S.init(rcxres,rcyres,rczres,data.Num_Coils);
		    			S.generateEigenCoils(smaps, image_store);
				}break;
		
				case(WALSH):{
					// Image space eigen Method
					eigen_coils(smaps, image_store);
				}break;
		
				case(LOWRES):{ // E-spirit Code 
					
					// Sos Normalization 
					cout << "Normalize Coils" << endl;
					#pragma omp parallel for 
					for(int k=0; k<smaps(0).length(thirdDim); k++){
					for(int j=0; j<smaps(0).length(secondDim); j++){
					for(int i=0; i<smaps(0).length(firstDim); i++){
						for(int coil=0; coil< data.Num_Coils; coil++){
							complex<float>s(0.0,0.0);
							for(int e=0; e< smap_num_encodes;e++){
								s += image_store(e,coil)(i,j,k);
							}
							smaps(coil)(i,j,k) = s;
						}
					}}}
					sos_normalize(smaps);
					
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
				//ArrayWrite(IC,"Intensity.dat");
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
			{
				Array< Array< complex<float>, 3>,1> temp =  Alloc4DContainer< complex<float> >(rcxres,rcyres,rczres,data.Num_Coils);
				smaps.reference(temp);
			}

			for(int coil=0; coil< smaps.length(firstDim); coil++){
				smaps(coil)=complex<float>(1.0,0.0);
			}	
		
		}break;
	}
	
	
	if(smap_mask != SMAPMASK_NONE){
	
		for(int coil=0; coil< data.Num_Coils; coil++){
			#pragma omp parallel for 
			for(int k=0; k<smaps(0).length(thirdDim); k++){
				for(int j=0; j<smaps(0).length(secondDim); j++){
					for(int i=0; i<smaps(0).length(firstDim); i++){
						
						float r = 0.0;
						switch(smap_mask){
							case(SMAPMASK_CIRCLE):{
								float x = 2.0*( i - 0.5*((float)smaps(0).length(firstDim)) )/((float)smaps(0).length(firstDim)); 
								float y = 2.0*( j - 0.5*((float)smaps(0).length(secondDim)) )/((float)smaps(0).length(secondDim));
								r = sqrt( x*x + y*y ); 
														
							}break;
							
							case(SMAPMASK_SPHERE):{
								float x = 2.0*( i - 0.5*((float)smaps(0).length(firstDim)) )/((float)smaps(0).length(firstDim)); 
								float y = 2.0*( j - 0.5*((float)smaps(0).length(secondDim)) )/((float)smaps(0).length(secondDim));
								float z = 2.0*( k - 0.5*((float)smaps(0).length(thirdDim)) )/((float)smaps(0).length(thirdDim));
								r = sqrt( x*x + y*y + z*z ); 
							}break;
							
							default:{
								r = 0.0;
							}break;
						}
							
						if( r > 1){
							smaps(coil)(i,j,k) = 0.0;
						}
			}}}
		}
	}

}


void RECON::sos_normalize( Array< Array<complex<float>,3>,1>&A){

	if(smap_thresh==0.0){
		
		#pragma omp parallel for 
		for(int k=0; k<A(0).length(thirdDim); k++){
			for(int j=0; j<A(0).length(secondDim); j++){
				for(int i=0; i<A(0).length(firstDim); i++){
					
					// The sum of squares
					float sos = 0.0;				
					for(int coil=0; coil< A.length(firstDim); coil++){
						sos +=  norm(A(coil)(i,j,k));
					}
					sos = sqrt(sos);
					
					// Normalize
					for(int coil=0; coil< A.length(firstDim); coil++){
						A(coil)(i,j,k) /= sos;
					}
		}}} // x,y,z
	
	}else{
	
	// Get the threshold
	Array<float,3>IC(rcxres,rcyres,rczres, ColumnMajorArray<3>());
	IC = 0.0;
	for(int coil=0; coil< A.length(firstDim); coil++){
		#pragma omp parallel for 
		for(int k=0; k<A(0).length(thirdDim); k++){
			for(int j=0; j<A(0).length(secondDim); j++){
				for(int i=0; i<A(0).length(firstDim); i++){
					IC(i,j,k) += norm(A(coil)(i,j,k));

		}}}		
	}
	IC = sqrt(IC);
	float max_IC = max(IC);
	cout << "Max IC = " << max_IC << endl;
	cout << "Thresh of " << max_IC*smap_thresh << endl;
			
	// Create a binary mask
	Array<int,3>MASK(rcxres,rcyres,rczres, ColumnMajorArray<3>());
	MASK = 0;
	for(int k=0; k<A(0).length(thirdDim); k++){
		for(int j=0; j<A(0).length(secondDim); j++){
			for(int i=0; i<A(0).length(firstDim); i++){
				if( IC(i,j,k) < smap_thresh*max_IC ){
					MASK(i,j,k) = 0.0;
				}else{
					MASK(i,j,k) = 1.0;
				}
	}}}
	//ArrayWrite(MASK,"PreFill.dat");
	
	// Try to remove the holes
	for(int k=0; k< MASK.length(thirdDim); k++){
		for(int i=0; i< MASK.length(firstDim); i++){
			
			int leading_edge= -1;
			for(int j=0; j < MASK.length(secondDim);  j++){
				if( MASK(i,j,k) > 0 ){
					leading_edge = j;
					break;
				}
			}
						
			int trailing_edge= -1;
			for(int j=MASK.length(secondDim)-1; j>=0; j--){
				if( MASK(i,j,k) > 0 ){
					trailing_edge = j;
					break;
				}
			}
			
			if( leading_edge > -1){
				for(int j=leading_edge; j <= trailing_edge; j++){
					MASK(i,j,k) = 1;
				}
			}
	}}
	//ArrayWrite(MASK,"PostFill.dat");
	
	// Now normalize
	for(int coil=0; coil< A.length(firstDim); coil++){
		#pragma omp parallel for 
		for(int k=0; k<A(0).length(thirdDim); k++){
			for(int j=0; j<A(0).length(secondDim); j++){
				for(int i=0; i<A(0).length(firstDim); i++){
					if( MASK(i,j,k) > 0){
						A(coil)(i,j,k) /= IC(i,j,k);
					}else{
						A(coil)(i,j,k) = complex<float>(0.0,0.0);
					}
		}}}
	}
	
	}
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
	int dwinX = (int)( 5*sigmaX);
	int dwinY = (int)( 5*sigmaY);
	int dwinZ = (int)( 5*sigmaZ);
	
	if( sigmaX > 0){
		// Kernel to reduce calls to exp
		Array< float,1> kern(2*dwinX+1);	
		for(int t=0; t< (2*dwinX +1); t++){
			kern(t) = exp( -sqr( (float)t - dwinX ) / (2.0*sqr(sigmaX))); 
		}				
		kern /= sum(kern);
			
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
		kern /= sum(kern);
			
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
		kern /= sum(kern);
			
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
void RECON::eigen_coils( Array< Array< complex<float>,3 >,1 > &smaps, Array< Array< complex<float>,3 >,2 > &image)
{
    		
	// Shorthand
	int Nencodes = image.extent(firstDim);
	int Ncoils = image.extent(secondDim);
	int Nx =image(0).extent(firstDim);
	int Ny =image(0).extent(secondDim);
	int Nz =image(0).extent(thirdDim);
	
	// Get coil with max signal
	cout << "Getting maximum coil" << endl;
	int ref_coil = 0;
	float max_sig = 0;
	for( int coil =0; coil < Ncoils; coil++){
		
		float sig = 0;
		for(int e=0; e< Nencodes; e++){
			sig += sum(norm(image(e,coil)));
		}
		
		if( sig > max_sig){
			ref_coil = coil;
			max_sig = sig;
		}
	}
	cout << "Reference coil = " << ref_coil << endl;
	
	int block_size_x = walsh_block_sizeX;
	int block_size_y = walsh_block_sizeY;
	int block_size_z = walsh_block_sizeZ;
		
	// Blocks shouldn't be larger than the dimension
	block_size_x = ( block_size_x > Nx) ? ( Nx ) : ( block_size_x );
	block_size_y = ( block_size_y > Ny) ? ( Ny ) : ( block_size_y );
	block_size_z = ( block_size_z > Nz) ? ( Nz ) : ( block_size_z );
	
	int block_hsize_x = block_size_x/2;
	int block_hsize_y = block_size_y/2;
	int block_hsize_z = block_size_z/2;
			
	int Np = block_size_x*block_size_y*block_size_z*Nencodes;
	
	cout << "Eigen Coil ( " << Nx << " x " << Ny << " x " << Nz << " x " << Ncoils << " x " << Nencodes << endl;
	cout << "Actual Block Size = " << block_size_x << " x " << block_size_y << " x " << block_size_z << endl;	
	cout << "Getting Low Rank threshold (N=" << Ncoils << ")(Np = " << Np << ")" << endl;
	
	int block_Nx= (int)( Nx / block_size_x );
	int block_Ny= (int)( Ny / block_size_y );
	int block_Nz= (int)( Nz / block_size_z );
	
	int total_blocks = block_Nx * block_Ny * block_Nz;
	cout << "Total Block Size" << total_blocks << " ( " << block_Nx << "," << block_Ny << "," << block_Nz << ")" << endl;
	
	// This code is for accelerating the calculation
	int EaccX = 1;
	int EaccY = 1;
	int EaccZ = 1;
	
	cout << "Eigen Acceleration ( " << EaccX << " x " << EaccY  << " x " << EaccZ << " ) " << endl;
	
	int NNx = max( Nx/EaccX,1);
	int NNy = max( Ny/EaccY,1);
	int NNz = max( Nz/EaccZ,1);
	
	cout << "Eigen Size ( " << NNx << " x " << NNy  << " x " << NNz << " ) " << endl;
				
	int *N = new int[3];
	N[0] = NNx;
	N[1] = NNy;
	N[2] = NNz;
	
	tictoc T;
	T.tic(); 
	
	#pragma omp parallel for
	for( int block = 0; block < (NNx*NNy*NNz); block++){
		
		//cout << "Block " << block << " of " << total_blocks << endl;
		
		// Get the actual position
		int *I = new int[3];
		nested_workaround(block,N,I,3);
		int i = I[0]*EaccX;
		int j = I[1]*EaccY;
		int k = I[2]*EaccZ;
		delete [] I;
		
		//-----------------------------------------------------
		//   Block coordinates
		//-----------------------------------------------------
		int kstart = k - block_hsize_z;
		int kstop = k - block_hsize_z + block_size_z;
		if( kstart < 0){
			kstart =0;
			kstart = k + block_size_z;
		}
		
		if(kstop > Nz){
			kstop = Nz;
			kstart = Nz - block_size_z;
		}
		
		
		int jstart = j - block_hsize_y;
		int jstop = j - block_hsize_y + block_size_y;
		if( jstart < 0){
			jstart =0;
			jstart = j + block_size_y;
		}
		
		if(jstop > Ny){
			jstop = Ny;
			jstart = Ny - block_size_y;
		}
		
		int istart = i - block_hsize_x;
		int istop = i - block_hsize_x + block_size_x;
		if( istart < 0){
			istart =0;
			istart = i + block_size_x;
		}
		
		if(istop > Nx){
			istop = Nx;
			istart = Nx - block_size_x;
		}
		
		//-----------------------------------------------------
		//   Collect a Block 
		//-----------------------------------------------------
		
		arma::cx_mat R;
		R.zeros(Ncoils,Np);
		for(int c=0; c< Ncoils; c++){
			int count = 0;
			for(int e=0; e< Nencodes; e++){
				for(int kk=kstart; kk < kstop; kk++){
					for(int jj=jstart; jj < jstop; jj++){
						for(int ii=istart; ii < istop; ii++){
							R(c,count) = image(e,c)(ii,jj,kk);
							count++;
			}}}}
		}
		
		// SVD to Get Eigen Vector
		arma::cx_mat U;
		arma::cx_mat V;
		arma::vec s;
  		arma::svd_econ(U, s, V, R,"left");
								
		for(int c=0; c< Ncoils; c++){
			complex<double>temp = sqrt(s(0))*U(c,0)*polar<double>(1.0,-arg(U(ref_coil,0)));
			smaps(c)(i,j,k) = complex<float>((float)real(temp),(float)imag(temp));
		}
		
	}// Block (threaded)
	cout << "(eigen coil took " << T << ")" << endl << flush;
	
	//-----------------------------------------------------
	//   Interpolate for acclerated calculation
	//-----------------------------------------------------
	T.tic(); 
	if( (EaccX != 1) || (EaccY != 1) || (EaccZ != 1) ){
		cout << "Interpolating Smaps" << endl; 
		
		// Actually do the convolution
		#pragma omp parallel for
		for( int k=0; k < Nz; k++){
		for( int j=0; j < Ny; j++){
		for( int i=0; i < Nx; i++){
				
			// Skip points already computed
			if( ((i%EaccX)==0) && ((j%EaccY)==0) && ((k%EaccZ)==0) ){
				continue;
			}
				
			// Grab nearest 8 points
			int i0 = i - i%EaccX;
			int j0 = j - j%EaccY;
			int k0 = k - k%EaccZ;
			
			int i1 = (i + EaccX + Nx)% Nx;
			int j1 = (j + EaccY + Ny)% Ny;
			int k1 = (k + EaccZ + Nz)% Nz;
						
			// Create the matrix
			arma::cx_mat C;
			C.zeros(Ncoils,8);
			for( int c=0; c < Ncoils; c++){
				C(c,0) = smaps(c)(i0,j0,k0);
				C(c,1) = smaps(c)(i1,j0,k0);
				C(c,2) = smaps(c)(i0,j0,k0);
				C(c,3) = smaps(c)(i1,j0,k0);
				C(c,4) = smaps(c)(i0,j1,k1);
				C(c,5) = smaps(c)(i1,j1,k1);
				C(c,6) = smaps(c)(i0,j1,k1);
				C(c,7) = smaps(c)(i1,j1,k1);
			}	
			
			// Phase the matrix via SVD
			arma::cx_mat U;
			arma::cx_mat V;
			arma::vec s;
  			arma::svd_econ(U, s, V, C,"left");	
			
			// Copy back
			for(int c=0; c< Ncoils; c++){
				smaps(c)(i,j,k) = complex<float>((float)real(C(c,0)),(float)imag(C(c,0)));
			}
			
		}}} // i,j,k
		cout << "Interpolation took " << T << endl;
	}
	
	// Normalize the maps
	#pragma omp parallel for
	for( int k=0; k < Nz; k++){
	for( int j=0; j < Ny; j++){
	for( int i=0; i < Nx; i++){
		
		double SS = 0;
		for(int c=0; c< Ncoils; c++){
			SS += norm( smaps(c)(i,j,k));
		}
		SS = 1./sqrt(SS);
		
		for(int c=0; c< Ncoils; c++){
			smaps(c)(i,j,k) *= SS;
		}
	}}}
}






