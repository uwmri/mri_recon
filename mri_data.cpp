//---------------------------------------------------
//    This function duplicates an existing MRI dataset
//---------------------------------------------------
#include "mri_data.h"

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
	//kdata = base_data.kdata;
}

MRI_DATA::MRI_DATA( void){
	Num_Encodings = -1;
	Num_Readouts = -1;
	Num_Pts = -1;
	Num_Coils = -1;
	Num_Slices = -1; // Default for 3D Non-Cartesian
}

//---------------------------------------------------
//    This function allocates and reads all data into memory
//---------------------------------------------------

void MRI_DATA::read_external_data( char *folder,int coils,int Ne,int Ns,int Npr,int Nx){

	FILE *fid;
	char fname[1024];
	
	Num_Encodings = Ne;
	Num_Readouts = Npr;
	Num_Pts = Nx;
	Num_Slices = Ns;
	Num_Coils = coils;
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
	
	cout << "Alloc Times " << endl;
	times.setStorage( ColumnMajorArray<4>());
	times.resize(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings); 
	cout << "Read  Times " << endl;
	for(int e=0; e<Num_Encodings; e++){
			
			while(1==1){
			 
			sprintf(fname,"%sTimes_VD_%d.dat",folder,e);
			if( (fid=fopen(fname,"r")) != NULL){
				cout << "\tUsing name " << fname << endl;
				break;
			}
			
			sprintf(fname,"%sTimes.dat",folder);
			if( (fid=fopen(fname,"r")) != NULL){
				cout << "\tUsing name same name for all encodings " << fname << endl;
				break;
			}
			
			cout << "Can't Open " << fname << " or other forms " << endl;
			break;
			}
			
			if(fid!=NULL){	
				Array<float,3>TimesRef = times(all,all,all,e);
				ArrayRead(TimesRef,fname);
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
	cout << "Max Kx = " << max(kx) << endl;
	cout << "Max Ky = " << max(ky) << endl;
	cout << "Max Kz = " << max(kz) << endl;
	cout << "Completed Reading Kx,Ky,Kz" << endl;
	cout << "Reading Kdata" << endl;
	
	kdata.setStorage( ColumnMajorArray<5>());
	kdata.resize(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings,Num_Coils); 
	
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
  cx_mat A;
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
  cx_mat U;
  vec s;
  cx_mat V;
  
  arma::svd_econ(U,s,V,A);
  
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

