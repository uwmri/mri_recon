//---------------------------------------------------
//    This function duplicates an existing MRI dataset
//---------------------------------------------------
#include "mri_data.h"

using arma::cx_fmat;
using arma::fvec;
using namespace NDarray;

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

void MRI_DATA::parse_external_header(const char *filename){
	
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


// Data for Whitening
void MRI_DATA::init_noise_samples(int total_samples){
	noise_samples.setStorage( ColumnMajorArray<2>());
	noise_samples.resize(total_samples,Num_Coils);
	noise_samples = complex<float>(0,0);
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

void MRI_DATA::read_external_data( const char *folder, int read_kdata){

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

void MRI_DATA::write_external_data( const char *folder){

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
		for( int i=0; i< (int)ecg.numElements(); i++){
			int temp = ( (int)(1e3*ecg(i)));
			endian_swap(temp);
			fwrite(&temp,1,sizeof(int),fid);
		}
		
		for( int i=0; i< (int)ecg.numElements(); i++){
			int temp = ( (int)(resp(i)));
			endian_swap(temp);
			fwrite(&temp,1,sizeof(int),fid);
		}
		
		for( int i=0; i< (int)ecg.numElements(); i++){
			int temp = ( (int)(1e6*time(i)));
			endian_swap(temp);
			fwrite(&temp,1,sizeof(int),fid);
		}
		
		for( int i=0; i< (int)ecg.numElements(); i++){
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



//--------------------------------------------------
//  Whiten Data
//--------------------------------------------------

void MRI_DATA::whiten(void){
	
	if( noise_samples.numElements()==0){
		cout << "Noise Samples do not exist:: Can't whiten data" << endl;
		return;
	}		
	
	cout << "Noise Samples : " << noise_samples.length(firstDim) << endl;
	cout << "Noise Pre-Whitening" << endl << flush;
  	
	// Copy into matrix
	arma::cx_mat NoiseData = arma::randu<arma::cx_mat>(noise_samples.length(secondDim),noise_samples.length(firstDim));
  	for(int coil=0; coil < Num_Coils; coil++){
	for(int i=0; i< noise_samples.length(firstDim); i++){
		NoiseData(coil,i) = noise_samples(i,coil);
	}}
	
	
	cout << "Calc Cov" << endl << flush;
	arma::cx_mat CV =  NoiseData*NoiseData.t();
	CV.save("CovMatrix.dat",arma::raw_binary);
		
	cout << "Whiten" << endl;
	arma::cx_mat V = chol(CV);
	arma::cx_mat VT = V.t(); 
	arma::cx_mat Decorr = VT.i();
			
	// Test Whitening
	arma::cx_mat W = NoiseData;
	arma::cx_mat temp = arma::randu<arma::cx_mat>(Num_Coils);
	for(int i =0; i< noise_samples.length(firstDim); i++){
		
		for(int coil=0; coil < Num_Coils; coil++){
			temp(coil,0)=W(coil,i);
		}
		arma::cx_mat temp2 = Decorr*temp;
		
		for(int coil=0; coil < Num_Coils; coil++){
			W(coil,i)=temp2(coil,0);
		}
	}
			
	arma::cx_mat CV_POST = W*W.t();
	CV_POST.save("CovMatrixPost.dat", arma::raw_binary);
	
	
	// Now Whiten Actual Data
	{
	cout << "Whiten all data" << endl;
	
	
	for(int e =0; e< kdata.length(fourthDim); e++){
	for(int k =0; k< kdata.length(thirdDim); k++){
	#pragma omp parallel for
	for(int j =0; j< kdata.length(secondDim); j++){
	for(int i =0; i< kdata.length(firstDim); i++){
		
		// Copy Sample
		arma::cx_mat temp(Num_Coils,1);
		for(int coil = 0; coil < kdata.length(fifthDim); coil++) {
   			temp(coil,0) = kdata(i,j,k,e,coil);
		}
		
		arma::cx_mat temp2 = Decorr*temp;
		
		// Copy Back
		for(int coil = 0; coil < kdata.length(fifthDim); coil++) {
   			 kdata(i,j,k,e,coil) = temp2(coil,0);
		}
	}}}}
	
	}/*actual whitening*/
	
}


/*----------------------------------------------
     Loads the Gating file into memory (temp)
 *----------------------------------------------*/ 
void MRI_DATA::load_pcvipr_gating_file(const char *full_filename){

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



