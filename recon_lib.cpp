
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
	
	sp_maps = 0;
	sp_mapshrink = 4;
	sp_mapthresh = 0.0;
	

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
		
		// Coil Combination + Resolution		
		float_flag("-lp_frac",lp_frac);
		float_flag("-smap_res",smap_res);
		
		// Time Resolved Flags
		int_flag("-frames",frames);
		
		// Spirit Coils
		trig_flag(1,"-sp_maps",sp_maps);
		int_flag("-sp_mapshrink",sp_mapshrink);
		float_flag("-sp_mapthresh",sp_mapthresh);
		
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


