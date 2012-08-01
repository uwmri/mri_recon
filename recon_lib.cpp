
#include "recon_lib.h"


RECON::RECON(int numarg, char **pstring){
	
	// Default Values
	recon_type = RECON_SOS;
	data_type = RECON_PFILE;
	
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
	    
	lp_frac=1.0;
	smap_res=16;
	
	frames = 1;
	
	acc = 1;
	compress_coils = 0.0;
	


#define trig_flag(num,name,val)   }else if(strcmp(name,pstring[pos]) == 0){ val = num; 
#define float_flag(name,val)  }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atof(pstring[pos]); 
#define int_flag(name,val)    }else if(strcmp(name,pstring[pos]) == 0){ pos++; val = atoi(pstring[pos]);
#define char_flag(name,val)   }else if(strcmp(name,pstring[pos]) == 0){ pos++; strcpy(val,pstring[pos]);

  for(int pos=0; pos < numarg; pos++){
  
  	if (strcmp("-h", pstring[pos] ) == 0) {
	  	printf("\n*********************************************\n");
	  	printf(" Basic Recon Control:\n");
	  	printf("*********************************************\n");
	  
		char_flag("-f",filename);
		
		// Reconstruction Geometry
		int_flag("-rcxres",rcxres);
		int_flag("-rcyres",rcyres);
		int_flag("-rczres",rczres);
		float_flag("-zoom",zoom);
		float_flag("-zoom_x",zoom_x);
		float_flag("-zoom_y",zoom_y);
		float_flag("-zoom_z",zoom_z);
		
		// Type of Recons		
		trig_flag(RECON_SOS,"-sos",recon_type);
		trig_flag(RECON_CG,"-isense",recon_type);
		trig_flag(RECON_PILS,"-pils",recon_type);
		trig_flag(RECON_IST,"-ist",recon_type);
		
		// Source of data
		trig_flag(RECON_EXTERNAL,"-external_data",data_type);
		trig_flag(RECON_PFILE,"-pfile",data_type);
		
		// Data modification
		int_flag("-acc",acc);
		float_flag("-compress_coils",compress_coils);
		
		// Coil Combination + Resolution		
		float_flag("-lp_frac",lp_frac);
		float_flag("-smap_res",smap_res);
		
		// Time Resolved Flags
		int_flag("-frames",frames);
		
		// Spirit Coils
		trig_flag(1,"-sp_maps",sp_maps);
		
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



//---------------------------------------------------
//    This function duplicates an esting MRI dataset
//---------------------------------------------------

// Constructer for MRI data type
MRI_DATA::MRI_DATA( MRI_DATA *base_data){
	
	// Copy Size
	Num_Encodings = base_data->Num_Encodings;
	Num_Readouts = base_data->Num_Readouts;
	Num_Pts = base_data->Num_Pts;
	Num_Coils = base_data->Num_Coils;
	
	// Allocate Memory and Copy Values
	kx.alloc(Num_Encodings,Num_Readouts,Num_Pts);  
	kx = base_data->kx;
	
	ky.alloc(Num_Encodings,Num_Readouts,Num_Pts);  
	ky = base_data->ky;
	
	kz.alloc(Num_Encodings,Num_Readouts,Num_Pts);  
	kz = base_data->kz;
	
	kw.alloc(Num_Encodings,Num_Readouts,Num_Pts);  
	kw = base_data->kw;
	
	cout << "4D Copy " << endl;
	kdata.alloc(Num_Coils,Num_Encodings,Num_Readouts,Num_Pts);
	//kdata = base_data.kdata;
}

MRI_DATA::MRI_DATA( void){
	Num_Encodings = -1;
	Num_Readouts = -1;
	Num_Pts = -1;
	Num_Coils = -1;
}

//---------------------------------------------------
//    This function allocates and reads all data into memory
//---------------------------------------------------

void MRI_DATA::read_external_data( char *folder,int coils,int Ne,int Npr,int Nx){

	FILE *fid;
	char fname[1024];
	
	Num_Encodings = Ne;
	Num_Readouts = Npr;
	Num_Pts = Nx;
	Num_Coils = coils;
	
	cout << "Data size= " << Num_Coils << " coils x " << Num_Encodings << " encodings x "<<  Num_Readouts << " readouts x" << Num_Pts << " pts" << endl;
	cout << "Alloc Kx " << endl;
	kx.alloc(Num_Encodings,Num_Readouts,Num_Pts);
	cout << "Read  Kx " << endl;
	for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKMAPX_VD_%d.dat",folder,e);
			if( (fid=fopen(fname,"r")) == NULL){
				cout << "Can't Open " << fname << endl;
				cout << "Exiting" << endl;
				exit(1);
			}else{	
				int j;
				if( (j=fread(kx[e][0],sizeof(float),Num_Readouts*Num_Pts,fid)) != (Num_Readouts*Num_Pts)){
					cout << "Not enough data: only read " << j << "points" << endl;
				}
				fclose(fid);
			}
	}
		
	cout << "Alloc Ky " << endl;
	ky.alloc(Num_Encodings,Num_Readouts,Num_Pts);
	cout << "Read  Ky " << endl;
	for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKMAPY_VD_%d.dat",folder,e);
			if( (fid=fopen(fname,"r")) == NULL){
				cout << "Can't Open " << fname << endl;
				cout << "Exiting" << endl;
				exit(1);
			}else{	
				if( (int)fread(ky[e][0],sizeof(float),Num_Readouts*Num_Pts,fid) != (Num_Readouts*Num_Pts) ){
					cout << "Can't Read Ky" << endl;
					exit(1);
				}
				fclose(fid);
			}
	}

	cout << "Alloc Kz " << endl;
	kz.alloc(Num_Encodings,Num_Readouts,Num_Pts);
	cout << "Read  Kz " << endl;
	for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKMAPZ_VD_%d.dat",folder,e);
			if( (fid=fopen(fname,"r")) == NULL){
				cout << "Can't Open " << fname << endl;
				cout << "Assume 2D" << endl;
			}else{	
				if( (int)fread(kz[e][0],sizeof(float),Num_Readouts*Num_Pts,fid) != (Num_Readouts*Num_Pts)){
					cout << "Can't Read Kz" << endl;
					exit(1);
				}
				fclose(fid);
			}
	}

	cout << "Alloc Kw " << endl;
	kw.alloc(Num_Encodings,Num_Readouts,Num_Pts);
	cout << "Read  Kw " << endl;
	for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKWEIGHT_VD_%d.dat",folder,e);
			if( (fid=fopen(fname,"r")) == NULL){
				cout << "Can't Open " << fname << endl;
				cout << "Trying Other Name" << endl;
				sprintf(fname,"%sKWEIGHT.dat",folder);
				if((fid=fopen(fname,"r")) == NULL){	
					cout << "Can't Open " << fname << " either " << endl;
					exit(1);
				}
			}
			
			if(fid!=NULL){
				cout << "Read Kw "  << endl;
				if( (int)fread(kw[e][0],sizeof(float),Num_Readouts*Num_Pts,fid) != (Num_Readouts*Num_Pts)){
					cout << "Can't Read Kw" << endl;
					exit(1);
				}
				fclose(fid);
			}
	}
	cout << "Max Kx = " << kx.max() << endl;
	cout << "Max Ky = " << ky.max() << endl;
	cout << "Max Kz = " << kz.max() << endl;
	cout << "Completed Reading Kx,Ky,Kz" << endl;
	cout << "Reading Kdata" << endl;
	
	kdata.alloc(Num_Coils,Num_Encodings,Num_Readouts,Num_Pts);
	
	for(int c=0; c<Num_Coils;c++){
		cout << "Reading CoilL: " << c << endl;
		for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKDATA_C%02d_VD_%d.dat",folder,c,e);
			if( (fid=fopen(fname,"r")) == NULL){
				cout << "Can't Open " << fname << endl;
				cout << "Trying without encoding" << endl;
				sprintf(fname,"%sKDATA_C%02d.dat",folder,c);
				if((fid=fopen(fname,"r")) == NULL){	
					cout << "Can't Open " << fname << " either " << endl;
					exit(1);
				}
			}
			
			if(fid!=NULL){	
				if( (int)fread(kdata[c][e][0],sizeof( complex<float>),Num_Readouts*Num_Pts,fid) != (Num_Readouts*Num_Pts)){
					cout << "Error can read data for coil " << c << endl;
					exit(1);
				}
				fclose(fid);
			}
	}}
	
	
}

/** Undersample the data by deleting argument 'us' readouts
*
* TODO:
* -Need to check if the num_readouts value in the RECON class needs to change
* -Needs to handle non radial trajectories as they are implemented
*/
void MRI_DATA::undersample(int us){
  if (us>1) { // Check argument in case we want to always run this function
  
  printf("Retrospectively undersampling by: %d (%d -> %d)\n", us, Num_Readouts, Num_Readouts/us); 
    
  Num_Readouts = Num_Readouts/us;
  
  // Undersampled arrays	
	array3D<float> kx_us;
	array3D<float> ky_us;
	array3D<float> kz_us;
	array3D<float> kw_us;
	array4D< complex<float> > kdata_us;
	
	// Allocate Memory and Copy Values
	kx_us.alloc(Num_Encodings,Num_Readouts,Num_Pts); 
	ky_us.alloc(Num_Encodings,Num_Readouts,Num_Pts);  
	kz_us.alloc(Num_Encodings,Num_Readouts,Num_Pts);  
	kw_us.alloc(Num_Encodings,Num_Readouts,Num_Pts);  
	kdata_us.alloc(Num_Coils,Num_Encodings,Num_Readouts,Num_Pts);
	
	for(int e = 0; e < Num_Encodings; e++) {
	  for(int r = 0; r < Num_Readouts; r++) {
	    for(int p = 0; p < Num_Pts; p++) {
	       kx_us[e][r][p] = kx[e][r*us][p];
	       ky_us[e][r][p] = ky[e][r*us][p];   
	       kz_us[e][r][p] = kz[e][r*us][p];   
	       kw_us[e][r][p] = kw[e][r*us][p];   
	}}}
	
	for(int c = 0; c < Num_Coils; c++) {
	  for(int e = 0; e < Num_Encodings; e++) {
      for(int r = 0; r < Num_Readouts; r++) {
        for(int p = 0; p < Num_Pts; p++) {
           kdata_us[c][e][r][p] = kdata[c][e][r*us][p];
	}}}}
	
	// Free the full arrays and then point them to the undersampled ones
	kx.freeArray();
	ky.freeArray();
	kz.freeArray();
	kw.freeArray();
	kdata.freeArray();
	
	kx.point_to_3D(&kx_us);
	ky.point_to_3D(&ky_us);
	kz.point_to_3D(&kz_us);
	kw.point_to_3D(&kw_us);
	kdata.point_to_4D(&kdata_us);
	
	}
}

/** Coil compress data with a cutoff of thresh*max(SV)
*
*/
void MRI_DATA::coilcompress(float thresh)
{ 
  cout << "Coil compression . . . " << flush;
  cx_mat A;
  A.zeros(kdata.Nx*kdata.Ny, kdata.Nt);
  
  for(int i = 0; i < kdata.Nt; i++) {
    for(int j = 0; j < kdata.Nx*kdata.Ny; j++) {
      A(j,i) = kdata[i][0][0][j];
  }}
  
  cx_mat U;
  vec s;
  cx_mat V;
  
  arma::svd_econ(U,s,V,A);
  
  s = s/s(0);
  
  //cout << s << endl << endl;
  
  uvec cc_find = arma::find(s > thresh, 1, "last");
  int cc_ncoils = cc_find(0)+1; // New number of coils
  
  cout << "(" << Num_Coils << " coils -> " << cc_ncoils << " coils) . . . " << flush;
  
  Num_Coils = cc_ncoils;
  
  A = A*V.cols(0,cc_ncoils-1);
  
  array4D< complex<float> > kdata_cc;
  kdata_cc.alloc(Num_Coils,Num_Encodings,Num_Readouts,Num_Pts);
  
  for(int i = 0; i < kdata_cc.Nt; i++) {
    for(int j = 0; j < kdata_cc.Nx*kdata_cc.Ny; j++) {
      kdata_cc[i][0][0][j] = A(j,i);
  }}
  
  kdata.freeArray();
	kdata.point_to_4D(&kdata_cc);
  
  cout << "done" << endl;
}


