//---------------------------------------------------
//    This function duplicates an existing MRI dataset
//---------------------------------------------------
#include "mri_data.h"

using arma::cx_fmat;
using arma::fvec;

void MRI_DATA::data_stats(void){
	
	
	if( kx.numElements()==0){
		cout << "Kspace does not exist" << endl;
	}else{
		cout << "Range Kx = " << min(kx) << " to " << max(kx) << endl;
		cout << "Range Ky = " << min(ky) << " to " << max(ky) << endl;
		cout << "Range Kz = " << min(kz) << " to " << max(kz) << endl;
		cout << "Range Kw = " << min(kw) << " to " << max(kw) << endl;
	}
	
	if( kdata.numElements()==0){
		cout << "Kdata does not exist yet" << endl;
	}else{
		cout << "Range Kdata = " << min(abs(kdata)) << " to " << max(abs(kdata)) << endl;
	}
	
	// gating
	if( ecg.numElements()==0){
		cout << "Physiologic data does not exist yet" << endl;
	}else{
		cout << "Range ECG = " << min(ecg) << " to " << max(ecg) << endl;
		cout << "Range RESP = " << min(resp) << " to " << max(resp) << endl;
		cout << "Range TIME = " << min(time) << " to " << max(time) << endl;
		cout << "Range PREP = " << min(prep) << " to " << max(prep) << endl;
	}
}


//--------------------------------------------------
//  Read external header
//--------------------------------------------------

void MRI_DATA::parse_external_header(char *filename){
	
	char parameter[1024];
	char cvalue[1024];
	float value;
	FILE *fid;
	char line[1024];
	
	cout << "Reading External Header: " << endl;
	
	fid = fopen(filename,"r");
	while( fgets(line, sizeof(line),fid) != NULL ) {
    	// Float Reads
		if(sscanf(line,"%s\t%f",parameter,&value) < 2){
		}else if(strcmp("acq_bw",parameter) == 0){ 
		}else if(strcmp("xres",parameter) == 0){  Num_Pts = (int)value;
		}else if(strcmp("numrecv",parameter) == 0){ Num_Coils = (int)value;
		}else if(strcmp("slices",parameter) == 0){ Num_Slices = (int)value;
		}else if(strcmp("2d_flag",parameter) == 0){ 
			if( (int)value ==1){
				trajectory_dims = TWOD;
			}else{
				trajectory_dims = THREED;
			}		
		}else if(strcmp("nproj",parameter) == 0){  Num_Readouts = (int)value;
		}else if(strcmp("rcxres",parameter) == 0){ 
			xres = (int)value;
		}else if(strcmp("rcyres",parameter) == 0){ 
			yres = (int)value;
		}else if(strcmp("rczres",parameter) == 0){ 
			zres = (int)value;
		}else if(strcmp("num_encodes",parameter) == 0){ 
			Num_Encodings = (int)value;
		}
		
		// Text Reads
		if(sscanf(line,"%s\t%s",parameter,cvalue) < 2){
		}else if(strcmp("gate_name",parameter) == 0){ 
			strcpy(gate_name,cvalue);
			printf("Gate name = %s\n",gate_name);
		}
	}
	fclose(fid);
}


// Constructer for MRI data type
MRI_DATA::MRI_DATA( MRI_DATA *base_data){
	
	// Copy Size
	Num_Encodings = base_data->Num_Encodings;
	Num_Readouts = base_data->Num_Readouts;
	Num_Pts = base_data->Num_Pts;
	Num_Coils = base_data->Num_Coils;
	Num_Slices = base_data->Num_Slices;
	
	// Allocate Memory and Copy Values
	kx.setStorage( ColumnMajorArray<4>());
	kx.resize(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings);
	kx = base_data->kx;
	
	ky.setStorage( ColumnMajorArray<4>());
	ky.resize( kx.shape());  
	ky = base_data->ky;
	
	kz.setStorage( ColumnMajorArray<4>());
	kz.resize( kx.shape());  
	kz = base_data->kz;
	
	kw.setStorage( ColumnMajorArray<4>());
	kw.resize( kx.shape());  
	kw = base_data->kw;
	
	cout << "5D Copy " << endl;
	kdata.setStorage( ColumnMajorArray<5>());
	kdata.resize( Num_Pts,Num_Readouts,Num_Slices,Num_Encodings,Num_Coils);
	kdata = base_data->kdata;
}


// Constructer for MRI data type
void MRI_DATA::init_memory(void){
	
	// Allocate Memory and Copy Values
	kx.setStorage( ColumnMajorArray<4>());
	kx.resize(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings);
	
	ky.setStorage( ColumnMajorArray<4>());
	ky.resize( kx.shape());  
	
	kz.setStorage( ColumnMajorArray<4>());
	kz.resize( kx.shape());  
	
	kw.setStorage( ColumnMajorArray<4>());
	kw.resize( kx.shape());  
	
	kdata.setStorage( ColumnMajorArray<5>());
	kdata.resize( Num_Pts,Num_Readouts,Num_Slices,Num_Encodings,Num_Coils);

	// Timesd
	time.setStorage( ColumnMajorArray<3>());
	time.resize(Num_Readouts,Num_Slices,Num_Encodings);

	ecg.setStorage( ColumnMajorArray<3>());
	ecg.resize(Num_Readouts,Num_Slices,Num_Encodings);

	resp.setStorage( ColumnMajorArray<3>());
	resp.resize(Num_Readouts,Num_Slices,Num_Encodings);

	prep.setStorage( ColumnMajorArray<3>());
	prep.resize(Num_Readouts,Num_Slices,Num_Encodings);

}


MRI_DATA::MRI_DATA( void){
	Num_Encodings = -1;
	Num_Readouts = -1;
	Num_Pts = -1;
	Num_Coils = -1;
	Num_Slices =  1; 
}



//---------------------------------------------------
//    This function allocates and reads all data into memory
//---------------------------------------------------

void MRI_DATA::read_external_data( char *folder, int read_kdata){

	FILE *fid;
	char fname[1024];
	
	Range all = Range::all();
	
	cout << "Data size= " << Num_Coils << " coils x " << Num_Encodings << " encodings x "<< Num_Slices<< " slices x "<<  Num_Readouts << " readouts x" << Num_Pts << " pts" << endl;
	
	cout << "Alloc Kx " << endl;
	kx.setStorage( ColumnMajorArray<4>());
	kx.resize(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings); 
	cout << "Read  Kx " << endl;
	for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKMAPX_VD_%d.dat",folder,e);
			if( (fid=fopen(fname,"r")) == NULL){
				cout << "Can't Open " << fname << endl;
				cout << "Exiting" << endl;
				exit(1);
			}else{	
				Array<float,3>KxRef = kx(all,all,all,e);
				ArrayRead(KxRef,fname);
				fclose(fid);
			}
	}
		
	cout << "Alloc Ky " << endl;
	ky.setStorage( ColumnMajorArray<4>());
	ky.resize(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings); 
	cout << "Read  Ky " << endl;
	for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKMAPY_VD_%d.dat",folder,e);
			if( (fid=fopen(fname,"r")) == NULL){
				cout << "Can't Open " << fname << endl;
				cout << "Exiting" << endl;
				exit(1);
			}else{	
				Array<float,3>KyRef = ky(all,all,all,e);
				ArrayRead(KyRef,fname);
				fclose(fid);
			}
	}

	cout << "Alloc Kz " << endl;
	kz.setStorage( ColumnMajorArray<4>());
	kz.resize(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings); 
	cout << "Read  Kz " << endl;
	for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKMAPZ_VD_%d.dat",folder,e);
			if( (fid=fopen(fname,"r")) == NULL){
				cout << "Can't Open " << fname << endl;
				cout << "Assume 2D" << endl;
			}else{	
				Array<float,3>KzRef = kz(all,all,all,e);
				ArrayRead(KzRef,fname);
				fclose(fid);
			}
	}
	
	cout << "Alloc Kw " << endl;
	kw.setStorage( ColumnMajorArray<4>());
	kw.resize(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings); 
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
				Array<float,3>KwRef = kw(all,all,all,e);
				ArrayRead(KwRef,fname);
				fclose(fid);
			}
	}
	cout << "Completed Reading Kspace Sampling" << endl;
	
	// Now Read Physiologic Data
	MRI_DATA::load_pcvipr_gating_file(gate_name); // TEMP	
	data_stats();
		
	
	cout << "Reading Kdata" << endl;
	kdata.setStorage( ColumnMajorArray<5>());
	kdata.resize(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings,Num_Coils); 
	
	
	if(read_kdata==1){
	
	for(int c=0; c<Num_Coils;c++){
		cout << "Reading CoilL: " << c << endl;
		for(int e=0; e<Num_Encodings; e++){
		
			while(1==1){
			 
			sprintf(fname,"%sKDATA_C%02d_VD_%d.dat",folder,c,e);
			if( (fid=fopen(fname,"r")) != NULL){
				cout << "\tUsing name " << fname << endl;
				break;
			}
			
			sprintf(fname,"%sKDATA_C%02d_V%02d.dat",folder,c,e);
			if( (fid=fopen(fname,"r")) != NULL){
				cout << "\tUsing name " << fname << endl;
				break;
			}
			
			sprintf(fname,"%sKDATA_C%02d.dat",folder,c);
			if( (fid=fopen(fname,"r")) != NULL){
				cout << "\tUsing name " << fname << endl;
				break;
			}
			
			cout << "Can't Open " << fname << " or other forms " << endl;
			exit(1);
			
			}
			
			Array< complex<float>,3>kdataC= kdata(Range::all(),Range::all(),Range::all(),e,c);
			if(fid!=NULL){	
				ArrayRead(kdataC,fname);
			}
	}}
	
	}// Read Kdata 
}



//---------------------------------------------------
//  Temporary Function to Write Data ( will be replaced by ismrmd ) 
//---------------------------------------------------

void MRI_DATA::write_external_data( char *folder){

	char fname[1024];
	mkdir(folder,0777);
	
	Range all = Range::all();
	
	cout << "Data size= " << Num_Coils << " coils x " << Num_Encodings << " encodings x "<< Num_Slices<< " slices x "<<  Num_Readouts << " readouts x" << Num_Pts << " pts" << endl;
	
	cout << "Write Kx " << endl;
	for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKMAPX_VD_%d.dat",folder,e);
			Array<float,3>KxRef = kx(all,all,all,e);
			ArrayWrite(KxRef,fname);
	}
	
	cout << "Write Ky " << endl;
	for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKMAPY_VD_%d.dat",folder,e);
			Array<float,3>KyRef = ky(all,all,all,e);
			ArrayWrite(KyRef,fname);
	}
	
	cout << "Write Kz " << endl;
	for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKMAPZ_VD_%d.dat",folder,e);
			Array<float,3>KzRef = kz(all,all,all,e);
			ArrayWrite(KzRef,fname);
	}		
	
	cout << "Write Kw " << endl;
	for(int e=0; e<Num_Encodings; e++){
			sprintf(fname,"%sKWEIGHT_VD_%d.dat",folder,e);
			Array<float,3>KwRef = kw(all,all,all,e);
			ArrayWrite(KwRef,fname);
	}		
	
	
	cout << "Writing Kdata" << endl;
	for(int c=0; c<Num_Coils;c++){
		cout << "Write Coil: " << c << endl;
		for(int e=0; e<Num_Encodings; e++){
		
			sprintf(fname,"%sKDATA_C%02d_VD_%d.dat",folder,c,e);
			Array< complex<float>,3>kdataC= kdata(Range::all(),Range::all(),Range::all(),e,c);
			ArrayWrite(kdataC,fname);
			
	}}
	
	
	// Write ECG file
	if( ecg.numElements()==0){
		cout << "Physiologic data does not exist yet" << endl;
	}else{
		sprintf(fname,"%sgating_track",folder);
		FILE *fid = fopen(fname,"w");
		for( int i=0; i< ecg.numElements(); i++){
			int temp = ( (int)(1e3*ecg(i)));
			endian_swap(temp);
			fwrite(&temp,1,sizeof(int),fid);
		}
		
		for( int i=0; i< ecg.numElements(); i++){
			int temp = ( (int)(resp(i)));
			endian_swap(temp);
			fwrite(&temp,1,sizeof(int),fid);
		}
		
		for( int i=0; i< ecg.numElements(); i++){
			int temp = ( (int)(1e6*time(i)));
			endian_swap(temp);
			fwrite(&temp,1,sizeof(int),fid);
		}
		
		for( int i=0; i< ecg.numElements(); i++){
			int temp = ( (int)(1e6*prep(i)));
			endian_swap(temp);
			fwrite(&temp,1,sizeof(int),fid);
		}
		fclose(fid);
		
	}
		
	// Write a Header
	sprintf(fname,"%sHeader.txt",folder);
	ofstream header;
  	header.open (fname);
  	header << "version " << 1 << endl;
	header << "numrecv " << Num_Coils << endl;
	header << "xres " << Num_Pts << endl;
	header << "nproj " << Num_Readouts << endl;
	header << "rcxres " << xres << endl;
	header << "rcyres " << yres << endl;
	header << "rczres " << zres << endl;
	header << "slices " << Num_Slices << endl;
	header << "num_encodes " << Num_Encodings << endl;	
	header << "2d_flag " << trajectory_dims << endl;
	header << "gate_name " << "gating_track" << endl;
	header.close();
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
    
  	int Num_Readouts_us = Num_Readouts/us;
  
  	// Undersampled arrays use resize to ensure contigous memory
	
	cout<< "\t Kx" << endl;
	Array<float,4>temp_us(Num_Pts, Num_Readouts, Num_Slices, Num_Encodings);
	temp_us = kx;
	Array<float,4> temp = temp_us(Range(fromStart,toEnd,1),Range(fromStart,toEnd,us),Range(fromStart,toEnd,1),Range(fromStart,toEnd,1));
	kx.resize(Num_Pts,Num_Readouts_us,Num_Slices,Num_Encodings);
	kx = temp;
	
	cout<< "\t Ky" << endl;
	temp_us = ky;
	temp = temp_us(Range(fromStart,toEnd,1),Range(fromStart,toEnd,us),Range(fromStart,toEnd,1),Range(fromStart,toEnd,1));
	ky.resize(kx.shape());
	ky = temp;
	
	cout<< "\t Kz" << endl;
	temp_us = kz;
	temp = temp_us(Range(fromStart,toEnd,1),Range(fromStart,toEnd,us),Range(fromStart,toEnd,1),Range(fromStart,toEnd,1));
	kz.resize(kx.shape());
	kz = temp;
	
	cout<< "\t Kw" << endl;
	temp_us = kw;
	temp = temp_us(Range(fromStart,toEnd,1),Range(fromStart,toEnd,us),Range(fromStart,toEnd,1),Range(fromStart,toEnd,1));
	kw.resize(kx.shape());
	kw = temp;
	
	cout<< "\t Kdata" << endl;
	Array<complex<float>,5>temp_data_us(kdata.shape());
	temp_data_us = kdata;
	Array<complex<float>,5> temp_data = temp_data_us(Range(fromStart,toEnd,1),Range(fromStart,toEnd,us),Range(fromStart,toEnd,1),Range(fromStart,toEnd,1),Range(fromStart,toEnd,1));
	kdata.resize(Num_Pts,Num_Readouts_us,Num_Slices,Num_Encodings,Num_Coils);
	kdata = temp_data;
	
	}
}


/** Coil compress data with a cutoff of thresh*max(SV)
*
*/
void MRI_DATA::coilcompress(float thresh)
{ 

  cout << "Coil compression . . .(thresh = " << thresh << " )" << flush;
  cx_fmat A;
  A.zeros(kdata.length(fourthDim)*kdata.length(thirdDim)*kdata.length(secondDim)*kdata.length(firstDim), kdata.length(fifthDim));
  
  cout << "Copy to Array" << endl << flush; 
  for(int coil = 0; coil < kdata.length(fifthDim); coil++) {
   	for(int e =0; e< kdata.length(fourthDim); e++){
	for(int k =0; k< kdata.length(thirdDim); k++){
	for(int j =0; j< kdata.length(secondDim); j++){
	for(int i =0; i< kdata.length(firstDim); i++){
		int offset = e*kdata.length(thirdDim)*kdata.length(secondDim)*kdata.length(firstDim)+  k*kdata.length(secondDim)*kdata.length(firstDim) + j*kdata.length(firstDim) + i;
		A(offset,coil) = kdata(i,j,k,e,coil);
  	}}}}
  }
  
  cout << "Size A " << A.n_rows << " x " << A.n_cols << endl;
  
  cout << "Svd" << endl << flush; 
  cx_fmat U;
  fvec s;
  cx_fmat V;
  arma::svd_econ(U,s,V,A,'r');
  s = s/s(0);
  
  //cout << s << endl << endl;
  //uvec cc_find = arma::find(s > thresh, 1, "last");
  //int cc_ncoils = cc_find(0)+1; // New number of coils
  int cc_ncoils = (int)thresh;
  
  cout << "(" << Num_Coils << " coils -> " << cc_ncoils << " coils) . . . " << flush;
  
  Num_Coils = cc_ncoils;
  
  A = A*V.cols(0,cc_ncoils-1);


  cout << "Copy back to Array" << endl; 
  kdata.resizeAndPreserve(kdata.length(firstDim),kdata.length(1),kdata.length(2),kdata.length(3),Num_Coils);
  for(int coil = 0; coil < kdata.length(fifthDim); coil++) {
   	for(int e =0; e< kdata.length(fourthDim); e++){
	for(int k =0; k< kdata.length(thirdDim); k++){
	for(int j =0; j< kdata.length(secondDim); j++){
	for(int i =0; i< kdata.length(firstDim); i++){
		int offset = e*kdata.length(thirdDim)*kdata.length(secondDim)*kdata.length(firstDim)+  k*kdata.length(secondDim)*kdata.length(firstDim) + j*kdata.length(firstDim) + i;
		kdata(i,j,k,e,coil) = A(offset,coil); 
  }}}}}
  

  cout << "done" << endl;
}


/*----------------------------------------------
     Loads the Gating file into memory (temp)
 *----------------------------------------------*/ 
void MRI_DATA::load_pcvipr_gating_file(char *full_filename){

	// Open the file
	FILE *fid;
	if ((fid = fopen (full_filename, "r")) == NULL){  
		printf("**********************************************");
		printf("\nERROR!!! Could not open gating file %s\n",full_filename);
      	printf("**********************************************");
		return;
	}
	
	// Get Size
	fseek (fid, 0, SEEK_END);
    int size=ftell(fid);
	int pts = size / 4/ sizeof(int);
	
	// To lazy to rewrite
	int *ecg_vals = new int[pts];
	int *resp_vals = new int[pts];
	int *acquisition_time = new int[pts];
	int *prep_time = new int[pts];
	rewind(fid);
	
	if( pts!=(int)fread(ecg_vals,sizeof(int),pts,fid)){ printf("Error Reading Tracking\n"); exit(1);}
	if( pts!=(int)fread(resp_vals,sizeof(int),pts,fid)){ printf("Error Reading Tracking\n"); exit(1);}
	if( pts!=(int)fread(acquisition_time,sizeof(int),pts,fid)){ printf("Error Reading Tracking\n"); exit(1);}
	if( pts!=(int)fread(prep_time,sizeof(int),pts,fid)){ printf("Error Reading Tracking\n"); exit(1);}
	fclose(fid);
	
	for(int pos=0;pos<pts;pos++){
		endian_swap(ecg_vals[pos]);
		endian_swap(resp_vals[pos]);
		endian_swap(acquisition_time[pos]);
		endian_swap(prep_time[pos]);
	}
	
	// Alloc new Arrays
	ecg.setStorage( ColumnMajorArray<3>());
	ecg.resize(Num_Readouts,Num_Slices,Num_Encodings);
	
	resp.setStorage( ColumnMajorArray<3>());
	resp.resize(Num_Readouts,Num_Slices,Num_Encodings);
	
	prep.setStorage( ColumnMajorArray<3>());
	prep.resize(Num_Readouts,Num_Slices,Num_Encodings);
	
	time.setStorage( ColumnMajorArray<3>());
	time.resize(Num_Readouts,Num_Slices,Num_Encodings);
	
	// Copy
	for(int e=0; e< Num_Encodings; e++){
	for(int slice=0; slice< Num_Slices; slice++){
	for(int view=0; view< Num_Readouts; view++){
		ecg(view,slice,e)= 1e-3*(float)ecg_vals[slice*Num_Readouts + view];
		resp(view,slice,e)=(float)resp_vals[slice*Num_Readouts + view];
		time(view,slice,e)=1e-6*(float)acquisition_time[slice*Num_Readouts + view];
		prep(view,slice,e)=1e-6*(float)prep_time[slice*Num_Readouts + view];
	}}}
}


#ifdef SUPERMAN

//  This is from old recon (temp)

/************************************************
Gating Libraries for pcvipr.e

Initial Author: Kevin M. Johnson
Description: This code contains functions for gating

*************************************************/

#include "gating_lib.h"

/*----------------------------------------------
     This cleans up memory
 *----------------------------------------------*/ 
GATING_ARRAY::~GATING_ARRAY(){
	/*
	delete [] ecg_vals;
	delete [] resp_vals;
	delete [] acquisition_time;
	delete [] prep_time;				
	delete [] store_order;*/
}

/*----------------------------------------------
     This converts the floating point user variable to
	 an actual time.
 *----------------------------------------------*/ 
void GATING_ARRAY::parse_filename(float ftime){
	
	printf("Ftime = %e\n",ftime);
	
	int itime=0;
	memcpy(&itime,&ftime,sizeof(int));
	
	printf("Itime = %d\n",itime);
	
	int gate_year 	= (int)( itime 					/( 12*32*24*60*60)) + 2000;
   	int gate_month 	= (int)( itime%(12*32*24*60*60) /(    32*24*60*60)) + 1;
   	int gate_day 	= (int)( itime%(32*24*60*60)	/(       24*60*60));
   	int gate_hour 	= (int)( itime%(24*60*60) 		/(          60*60));
   	int gate_minute = (int)( itime%(60*60) 			/(             60));
   	int gate_second = (int)( itime%(60) 			/(              1));
	
	/* TempTo Deal with PSD Missing Output
	 gate_year 	= 2012;
   	 gate_month 	= 1;
   	 gate_day 	= 5;
   	 gate_hour 	= 16;
   	 gate_minute = 20;
   	 gate_second = 31;
	*/
		
	sprintf(filename,"Gating_Track_y%d_m%d_d%d_h%d_m%02d_s%02d.pcvipr_track",gate_year,gate_month,gate_day,gate_hour,gate_minute,gate_second);
	printf("gating Name= %s\n",filename);
}

/*----------------------------------------------
     In this older format, there are different file
	  names for everything.
 *----------------------------------------------*/ 
void GATING_ARRAY::parse_filename_old(float ecg_file_id){
	sprintf(ecg_filename,"%s%06d",ECG_TRACK_FILE,(int)ecg_file_id);
	sprintf(resp_filename,"%s%06d",RESP_TRACK_FILE,(int)ecg_file_id);
	sprintf(prep_filename,"%s%06d",PREP_TRACK_FILE,(int)ecg_file_id);
}



inline void endian_swap( int& x)
{
    x = ( x<<24 & 0xFF000000) |
        ( x<<8  & 0x00FF0000) |
        ( x>>8  & 0x0000FF00) |
        ( x>>24 & 0x000000FF);
}

/*----------------------------------------------
     Loads the Gating file into memory
 *----------------------------------------------*/ 
void GATING_ARRAY::load_gating_file(char *data_folder){

	char full_filename[1024];
	sprintf(full_filename,"%s%s",data_folder,filename);
	
	
	FILE *fid;
		
	if ((fid = fopen (full_filename, "r")) == NULL){  
		printf("**********************************************");
		printf("\nERROR!!! Could not open file %s\n",full_filename);
      	printf("**********************************************");
		exit(1);
	}
	
	fseek (fid, 0, SEEK_END);
    int size=ftell(fid);
	pts = size / 4/ sizeof(int);
	
	printf("Points = %d\n",pts);	
		
	ecg_vals = new int[pts];
	resp_vals = new int[pts];
	acquisition_time = new int[pts];
	prep_time = new int[pts];
	rewind(fid);
	
	if( pts!=(int)fread(ecg_vals,sizeof(int),pts,fid)){ printf("Error Reading Tracking\n"); exit(1);}
	if( pts!=(int)fread(resp_vals,sizeof(int),pts,fid)){ printf("Error Reading Tracking\n"); exit(1);}
	if( pts!=(int)fread(acquisition_time,sizeof(int),pts,fid)){ printf("Error Reading Tracking\n"); exit(1);}
	if( pts!=(int)fread(prep_time,sizeof(int),pts,fid)){ printf("Error Reading Tracking\n"); exit(1);}
	fclose(fid);
	
	for(int pos=0;pos<pts;pos++){
		endian_swap(ecg_vals[pos]);
		endian_swap(resp_vals[pos]);
		endian_swap(acquisition_time[pos]);
		endian_swap(prep_time[pos]);
	}
	
	store_order = new int[pts];
  	for(int pos=0; pos <pts;pos++){
  		store_order[pos]=pos;
  	}
	
	gating_stats();
}

/*----------------------------------------------
     Loads the Gating file into memory for old files
 *----------------------------------------------*/ 
void GATING_ARRAY::load_gating_file_old(char *data_folder, int nt){

	char full_filename[1024];
	
	pts = nt;
	cout << "Pts= " << pts << "\n";
		
	FILE *fid;
	
	sprintf(full_filename,"%s%s",data_folder,ecg_filename);
	ecg_vals = new int[pts];
	if ((fid = fopen (full_filename, "r")) == NULL){  
		printf("***   WARNING - NO ECG FILE              **\n");
		memset(ecg_vals,0,(size_t)(sizeof(int)*pts));
		
	}else{
		if( pts!=(int)fread(ecg_vals,sizeof(int),pts,fid)){ printf("Error Reading Tracking\n"); exit(1);}
		fclose(fid);
		for(int pos=0;pos<pts;pos++){
			endian_swap(ecg_vals[pos]);
		}
	}
	
	sprintf(full_filename,"%s%s",data_folder,resp_filename);
	resp_vals = new int[pts];
	if ((fid = fopen (full_filename, "r")) == NULL){  
		printf("***   WARNING - NO RESP FILE             **\n");
		memset(resp_vals,0,(size_t)(sizeof(int)*pts));
		
	}else{
		if( pts!=(int)fread(resp_vals,sizeof(int),pts,fid)){ printf("Error Reading Tracking\n"); exit(1);}
		fclose(fid);
		for(int pos=0;pos<pts;pos++){
			endian_swap(resp_vals[pos]);
		}
	}


	sprintf(full_filename,"%s%s",data_folder,prep_filename);
	prep_time = new int[pts];
	if ((fid = fopen (full_filename, "r")) == NULL){  
		printf("***   WARNING - NO PREP FILE             **\n");
		memset(prep_time,0,(size_t)(sizeof(int)*pts));
		
	}else{
		if( pts!=(int)fread(prep_time,sizeof(int),pts,fid)){ printf("Error Reading Tracking\n"); exit(1);}
		fclose(fid);
		for(int pos=0;pos<pts;pos++){
			endian_swap(prep_time[pos]);
		}
	}
	
	acquisition_time = new int[pts];
	memset(acquisition_time,0,(size_t)(sizeof(int)*pts));	


 	store_order = new int[pts];
  	for(int pos=0; pos <pts;pos++){
  		store_order[pos]=pos;
  	}

	printf("Old File Acquisition Time Not in File\n");	
}

/*----------------------------------------------
     Smooths the Resp and subtracts off to correct drift
 *----------------------------------------------*/ 
void GATING_ARRAY::correct_resp_drift(int fsize){
	
	int *filtered_resp = new int[pts];
	
	for(int pos=0;pos< pts; pos++){
		int start =  ( pos > fsize ) ? ( pos-fsize ) : ( 0 );
		int stop =  ( pos > pts - fsize) ? (pts ) : ( pos + fsize);
		
		long int count = 0;
		long int val = 0;
		for(int fpos=start; fpos<stop; fpos++){
			val += resp_vals[fpos];
			count++;
		}
		filtered_resp[pos] = val / count;
	}
	
	for(int pos=0;pos< pts; pos++){
		resp_vals[pos] -= filtered_resp[pos] - 4095;
	}
	delete [] filtered_resp;
}

void GATING_ARRAY::export_gating_info(char *name){
	FILE *fid=fopen(name,"w");
	fwrite(ecg_vals , pts , sizeof(int) , fid);
	fwrite(resp_vals , pts , sizeof(int) , fid);
	fwrite(acquisition_time , pts , sizeof(int) , fid);
	fwrite(prep_time , pts , sizeof(int) , fid);
	fwrite(store_order , pts , sizeof(int) , fid);
	fclose(fid);
}

void GATING_ARRAY::bitreverse(int interleaves, unsigned int *bitreversal)
{
  unsigned int i;
  unsigned int j;
  unsigned int tmp;
  unsigned int num_bits;
  unsigned int diff_element;
  unsigned int diff;
  unsigned int *tmp_bitrev;
  unsigned int curr_high;

  if (interleaves == 1){
    bitreversal[0] = 0;
    return;
  }
  
  tmp_bitrev = (unsigned int *)malloc(interleaves * sizeof(unsigned int));
  for (i=0; i< interleaves;i++){
    tmp_bitrev[i] = 0;
  }
  num_bits = 1;
  tmp = interleaves - 1;
  /* find out how many bits in interleaves */
  while ((tmp = (tmp>>1)) > 0)
    num_bits++;
  
  
  /* Code to reverse bit order*/
  for (j = 0; j < interleaves; j++){
    tmp = j;
    for (i=0; i< num_bits; i++){
	  tmp_bitrev[j] += pow(2,num_bits - i - 1) * ( (tmp >> i)  & 0x00000001);
    }
  }


  /* In case, interleaves is not a power of 2 sort into order */
  bitreversal[0] = tmp_bitrev[0];
  curr_high = 0;

  for (j = 1; j < interleaves; j++){
    tmp = interleaves +1;
    diff_element = 0;
    for (i = 1; i < interleaves; i++)
      if (tmp_bitrev[i] > curr_high) {
	diff = tmp_bitrev[i] - curr_high;
	if (diff < tmp){
	  diff_element = i;
	  tmp = diff;
	}
      }
    curr_high = tmp_bitrev[diff_element];
    bitreversal[diff_element] = j;
  }
  
  free(tmp_bitrev);
}
	

#define K_VIPR 0
#define K_SOS 1
#define K_CART 2


void GATING_ARRAY::determine_acquisition_order(int trajectory_type,int angio_core,float tr,int slices,int inter,int pr_order){

	cout << "TR = " << tr << "\n";

	if( angio_core){
		cout << "Using Angio Core View Ordering\n";		

		unsigned int *bitreversal = (unsigned int *)malloc((inter + 1) * sizeof(unsigned int));
		bitreverse(inter, bitreversal);
		
		float time=0.0;
			
		/******** Get time of acquisition ***********/
		if(trajectory_type == K_VIPR){
			int subproj = pts/inter;
			cout << "Subproj=" << subproj <<"\n";
			cout << "Inter=" << inter <<"\n";
			
			for(int pos=0; pos< inter; pos++){
 				for( int sb=0; sb< subproj; sb++){
       				int proj = sb*inter + bitreversal[pos];
					acquisition_time[proj] = (int)time;
					time += tr;
			}}		
			
	  	}else if(pr_order !=9){
				
			for(int pos=0; pos< inter; pos++){
 				for( int sb=0; sb< pts/inter/slices; sb++){
					int proj = sb*inter + bitreversal[pos]; 											
					for(int pe = 0; pe < slices; pe++){
				   		acquisition_time[ (int)(pe * pts/inter + proj) ] = (int)time;
					}/*pe*/
					time += tr*slices;
			 	}
			}/*proj*/
		}else if(trajectory_type == K_SOS){
			cout << "Golden Angle=" << endl;
			
			for(int pos=0; pos< pts/slices; pos++){
 				for(int pe = 0; pe < slices; pe++){
					acquisition_time[ (int)(pe * pts/slices + pos) ] = (int)time;
				}/*pe*/
				time += tr*slices;
			}/*proj*/
		}/*K_SOS*/
	}else{
  		cout << "Using Cardiac Ordering\n";		
  		cout << "Pts=" << pts << "\n";
		cout << "Inter=" << inter <<"\n";		

		float time=0.0;
		unsigned int *bitreversal = (unsigned int *)malloc((inter + 1) * sizeof(unsigned int));
		bitreverse(inter, bitreversal);
		
		int subproj = pts/inter;
		cout << "Subproj=" << subproj <<"\n";		
		
		if(trajectory_type==K_VIPR){
			switch(pr_order){
			
			case(2):{
			for(int si=0; si< subproj; si++){
 	 	 		for(int j=0; j< inter; j++){
					int proj = si*inter + bitreversal[j];
					
					if( si < 5){
						cout << "Pos =" << (si*inter +j) << "Proj=" << proj << " Bitrev = " << bitreversal[j] << "\n";   					
					}
					acquisition_time[proj] = (int)time;
					time += tr;
			}}
			}break;
			
			case(4):{
			
			unsigned int *sub_order = (unsigned int *)malloc((subproj + 1) * sizeof(unsigned int));
			bitreverse(inter, bitreversal);
			bitreverse(subproj, sub_order);
		
			for( int pos=0; pos < pts; pos++){
				int inter_index = pos%inter; /* Interleave / Time Frame */
				int inter_pos = (int)( (float)pos/(float)inter); /* Projection within time frame*/
				int sub_pos = inter_index + inter_pos;
				
				int proj = inter*sub_order[sub_pos%subproj] + bitreversal[inter_index];
				acquisition_time[proj] = (int)time;
				time += tr;
				
				
			}
			}break;				
			
			}
		}else{
			int nhb = (pts/inter);
			for(int slice =0; slice < slices; slice++){
				for(int hb =0; hb < nhb; hb++){
					for( int j=0; j<inter; j++){
						int proj = slice*pts/slices + hb*inter + bitreversal[j];
						acquisition_time[proj]= (int)time;				
						time += tr;
					
			}}}
		}
	}
}


template <typename T>
T array_max(T *a, int l) {
  T result=a[0]; 
  for(int pos=0; pos< l; pos++){
  	result = ( a[pos] > result ) ? ( a[pos] ) : ( result); 
  }
  
  return (result);
}

template <typename T>
T array_min(T *a, int l) {
  
  T result=a[0]; 
  for(int pos=0; pos< l; pos++){
  	result = ( a[pos] < result ) ? ( a[pos] ) : ( result); 
  }
  
  return (result);
}

void GATING_ARRAY::check_badecg(void){
	int bad_ecg_max = 5000; /*5s is a long time for a heart not to beat*/
	int bad_ecg_min = 0; /*Can't have negative time*/
	
	for(int pos=0;pos< pts; pos++){
		if( ( ecg_vals[pos] > bad_ecg_max ) || ( ecg_vals[pos] <  bad_ecg_min )){
			// cout << "Found a bad gate with value (ms)" <<  ecg_vals[pos] << "\n";
			ecg_vals[pos]= NULL_GATE;
		}
	}
}

void GATING_ARRAY::gating_stats(void){
	
	printf("---------Gating Stats--------------\n");
	
	max_ecg= array_max<int>(ecg_vals,pts);
	min_ecg= array_min<int>(ecg_vals,pts);
	printf("ECG Range (ms) = %d - %d\n",min_ecg,max_ecg);
	
	max_resp= array_max<int>(resp_vals,pts);
	min_resp= array_min<int>(resp_vals,pts);
	printf("Resp Range = %d - %d\n",min_resp,max_resp);

	max_acquisition= array_max<int>(acquisition_time,pts);
	min_acquisition= array_min<int>(acquisition_time,pts);
	printf("Acquisition Range (s) = %f - %f\n",(float)min_acquisition/1e6,(float)max_acquisition/1e6);

	max_prep= array_max<int>(prep_time,pts);
	min_prep= array_min<int>(prep_time,pts);
	printf("Prep Range (ms) = %f - %f\n",(float)min_prep/1e3,(float)max_prep/1e3);
}


void BottomUpMerge(int *A, int iLeft, int iRight, int iEnd, int *B,int *OLD,int *OLD_B)
{
  int i0 = iLeft;
  int i1 = iRight;
  int j;
 
  /* while there are elements in the left or right lists */
  for (j = iLeft; j < iEnd; j++)
    {
      /* if left list head exists and is <= existing right list head */
      if (i0 < iRight && (i1 >= iEnd || A[i0] <= A[i1]))
        {
          B[j] = A[i0];
		  OLD_B[j] = OLD[i0];
          i0 = i0 + 1;
        }
      else
        {
          B[j] = A[i1];
		  OLD_B[j] = OLD[i1];
          i1 = i1 + 1;
        }
    }
}

/* array A[] has the items to sort; array B[] is a work array */
void GATING_ARRAY::sort(int *order){
  
  int width;
  int *B = new int[pts];
  
  /*For Tracking New Order*/
  int *A= new int[pts];
  int *OLD = new int[pts];
  int *OLD_B = new int[pts];
  memset(OLD_B,0,(size_t)(sizeof(int)*pts));
  for(int pos=0; pos <pts;pos++){
  	OLD[pos]=pos;
	A[pos] = order[pos];
  }
       
  
  /* each 1-element run in A is already "sorted". */
 
  /* Make successively longer sorted runs of length 2, 4, 8, 16... until whole array is sorted */
  for (width = 1; width < pts; width = 2 * width){
      
	  int i;
 
      /* array A is full of runs of length width */
      for (i = 0; i < pts; i = i + 2 * width)
        {
          /* merge two runs: A[i:i+width-1] and A[i+width:i+2*width-1] to B[] */
          /*  or copy A[i:n-1] to B[] ( if(i+width >= n) ) */
          BottomUpMerge(A, i, min(i+width, pts), min(i+2*width, pts),B,OLD,OLD_B);
        }
 
      /* now work array B is full of runs of length 2*width */
      /* copy array B to array A for next iteration */
      /*   a more efficient implementation would swap the roles of A and B */
      memcpy(A,B,(size_t)(sizeof(int)*pts));
      memcpy(OLD,OLD_B,(size_t)(sizeof(int)*pts));
      
	  /* now array A is full of runs of length 2*width */
    }
	
	
	/*Now Rearrange All the Values*/
	memcpy(B,ecg_vals,(size_t)(sizeof(int)*pts));
	for( int pos=0;pos < pts; pos++){
		ecg_vals[pos] = B[OLD[pos]]; 
	}	
	
	memcpy(B,resp_vals,(size_t)(sizeof(int)*pts));
	for( int pos=0;pos < pts; pos++){
		resp_vals[pos] = B[OLD[pos]]; 
	}	
	
	memcpy(B,acquisition_time,(size_t)(sizeof(int)*pts));
	for( int pos=0;pos < pts; pos++){
		acquisition_time[pos] = B[OLD[pos]]; 
	}	
	
	memcpy(B,prep_time,(size_t)(sizeof(int)*pts));
	for( int pos=0;pos < pts; pos++){
		acquisition_time[pos] = B[OLD[pos]]; 
	}	
	
	memcpy(B,store_order,(size_t)(sizeof(int)*pts));
	for( int pos=0;pos < pts; pos++){
		store_order[pos] = B[OLD[pos]]; 
	}	
	
	delete [] A;		
	delete [] B;
	delete [] OLD_B;
	delete [] OLD;
	
}


void GATING_ARRAY::shift_ecg(int shift){
	cout << "Shifting ECG by " << shift << "\n";
	for(int pos=0; pos< pts; pos++){
		
		int new_pos = pos + shift;
		if( (new_pos< pts) & (new_pos > 0) ){
			ecg_vals[pos] = ecg_vals[new_pos];
		}else{
			ecg_vals[pos] = NULL_GATE;
		}
	}
} 

#endif
















