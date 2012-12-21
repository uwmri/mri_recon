/* 
   	This contains Commandline Interface to the Recon
 */
#include "recon_lib.h"

int main(int argc, char **argv){
	
	// ------------------------------------
	// Initialize Recon - reading the command line
	// ------------------------------------
	RECON recon(argc,argv);
		
	// ------------------------------------
	// Read Data - 
	// ------------------------------------
		
	MRI_DATA data;
	cout << "----Read Data-----" << endl;	
	switch(recon.data_type){
		case(RECON::PFILE):{	
			// Read in P-File (doesn't work)
			//PFILE pfile;
			//pfile.read_header(recon.filename);
			//pfile.read_data(0);
		}break;
		
		case(RECON::SIMULATE):{
			// Use completely made up data
			
		}break;
		
		case(RECON::PHANTOM):{
			// Use external Kx,Ky,Kz 
			recon.parse_external_header();
			data.read_external_data("./",recon.num_coils,recon.rcencodes,recon.num_slices,recon.num_readouts,recon.xres,0);
			
			// Initialize Phantom
			PHANTOM phantom;
			phantom.read_commandline(argc,argv);  
			phantom.init(recon.rcxres,recon.rcyres,recon.rczres);
			
			// Accurate gridding for Phantom
			gridFFT phantom_gridding;
			phantom_gridding.kernel_type = KAISER_KERNEL;
			phantom_gridding.overgrid = 1.25;
			phantom_gridding.dwinX = 6;
			phantom_gridding.dwinY = 6;
			phantom_gridding.dwinZ = 6;
			phantom_gridding.precalc_gridding(phantom.IMAGE.length(firstDim),phantom.IMAGE.length(secondDim),phantom.IMAGE.length(thirdDim),3);
			
			//Collect Data
			int e= 0;
			Range all=Range::all();
			Array< float,3 >kxE = data.kx(all,all,all,e); 
			Array< float,3 >kyE = data.ky(all,all,all,e); 
			Array< float,3 >kzE = data.kz(all,all,all,e); 
			Array< float,3 >kwE = data.kw(all,all,all,e); 
			for(int coil =0; coil < recon.num_coils; coil++){
				cout << "Getting phantom" << coil << flush;
				cout << "Smap " << flush;
				phantom.update_smap_biotsavart(coil,recon.num_coils);				
				Array<complex<float>,3>kdataE = data.kdata(all,all,all,e,coil); 
				kdataE=0;			
				cout << " Grid " << flush;
				phantom_gridding.backward(phantom.IMAGE,phantom.SMAP,kdataE,kxE,kyE,kzE,kwE);
				cout << " Done " << endl;
			}
			
			// Add Noise
			phantom.add_noise( data.kdata );
			
			
		}break;
		
		default:
		case(RECON::EXTERNAL):{
			// Read in External Data Format
			recon.parse_external_header();
			data.read_external_data("./",recon.num_coils,recon.rcencodes,recon.num_slices,recon.num_readouts,recon.xres,1);
		}break;
	}

	cout << "----Geometry Modification-----" << endl;	
	// Geometry Modification by Recon
	data.kx *= ((float)(1.0/recon.zoom_x));
	data.ky *= ((float)(1.0/recon.zoom_y));
	data.kz *= ((float)(1.0/recon.zoom_z));
	
	if (recon.acc > 1){
		data.undersample(recon.acc);
	}

	// --------------------------------------------------
	// Code for recon (no PSD specific data/structures)
	// --------------------------------------------------
	Array< complex<float>,5 >X = recon.reconstruction(argc,argv,data);
	
	// ------------------------------------
	// Post Processing + Export
	// ------------------------------------
	
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


