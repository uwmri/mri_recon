//---------------------------------------------------
//    This function duplicates an existing MRI dataset
//---------------------------------------------------
#include "mri_data.h"

using arma::cx_fmat;
using arma::fvec;
using namespace NDarray;

void MRI_DATA::dump_stats(const string name, const Array< Array<float,3>,1> & in){
	cout << name << endl;
	cout << "\tContainer Size = " << in.length(firstDim) << endl;
	cout << "\tElement size = " << in(0).length(firstDim) << " x " << in(0).length(secondDim) << " x " << in(0).length(thirdDim) << endl;
	float high = 0;
	float low = 9999;
	for( Array< Array<float,3>,1> ::const_iterator miter=in.begin();   miter !=in.end(); miter++){
		high = max(high, max(*miter) );
		low = min(low, min(*miter));
	}
	cout << "\tRange = " << low << " to " << high << endl;
}

void MRI_DATA::dump_stats(const string name, const Array< Array<complex<float>,3>,2> & in){
	cout << name << endl;
	cout << "\tContainer Size = " << in.length(firstDim) << " x " << in.length(secondDim) << endl;
	cout << "\tElement size = " << in(0).length(firstDim) << " x " << in(0).length(secondDim) << " x " << in(0).length(thirdDim) << endl;
	float high = 0;
	float low = 9999;
	for( Array< Array< complex<float>,3>,2>::const_iterator miter=in.begin();   miter !=in.end(); miter++){
		float temp = max(abs(*miter));
		high = max(high,temp);
		
		temp = min(abs(*miter));
		low = min(low, temp);
	}
	cout << "\tRange = " << low << " to " << high << endl;
}

void MRI_DATA::scale_fov(float scale_x, float scale_y, float scale_z){
	
	if( (scale_x ==1.0) && (scale_y==1.0) && (scale_z==1.0) ){
		return;
	}
	
	cout << "Scaling the fov by" << scale_x << "," << scale_y << "," << scale_z << endl;
	
	// Multiply kspace by inverse
	for( Array< Array<float,3>,1>::iterator miter=kx.begin();   miter !=kx.end(); miter++){
		(*miter) *= ( 1./scale_x);
	}

	for( Array< Array<float,3>,1>::iterator miter=ky.begin();   miter !=ky.end(); miter++){
		(*miter) *= ( 1./scale_y);
	}
	
	for( Array< Array<float,3>,1>::iterator miter=kz.begin();   miter !=kz.end(); miter++){
		(*miter) *= ( 1./scale_z);
	}
		
	// Multiple image by scale
	zfov *= scale_z;
	yfov *= scale_y;
	xfov *= scale_x;
	
	if( sms_type == SMSon){
		for( Array< Array<float,3>,2>::iterator miter=z.begin();   miter !=z.end(); miter++){
			(*miter) *= scale_z;
		}
	}

}

MRI_DATA MRI_DATA::subframe( int eStart, int eStop, int eStride ){

	MRI_DATA data2;
	data2.clone_attributes(*this);
	
	// Copy info
	data2.Num_Encodings = floor( ( 1 + eStop - eStart)/eStride);
	data2.init_memory();
	
	int count = 0;
	for( int ee = eStart; ee<= eStop; ee+= eStride){
	
		data2.kx(count) = this->kx(ee);
		data2.ky(count) = this->ky(ee);
		data2.kz(count) = this->kz(ee);
		data2.kw(count) = this->kw(ee);
		data2.kt(count) = this->kt(ee);
		
		// SMS need z positions
		if(data2.sms_type == SMSon){
			for(int sms_pos=0; sms_pos < sms_factor; sms_pos++){
				data2.z(sms_pos,count) = this->z(sms_pos,ee);
			}
		}
			
		// Copy the data	
		for( int coil=0; coil< data2.Num_Coils; coil++){
			data2.kdata(count,coil) = this->kdata(ee,coil);
		}
		count = count+1;
	}

	return(data2);
}

void MRI_DATA::stats(void){

	if( kx.numElements()==0){
		cout << "Kspace does not exist" << endl;
	}else{
		dump_stats("Kx",kx);
		dump_stats("Ky",ky);
		dump_stats("Kz",kz);
		dump_stats("Kw",kw);
	}
	
	
	if( kdata.numElements()==0){
		cout << "Kdata does not exist yet" << endl;
	}else{
		dump_stats("Kdata",kdata);
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

// Data for Whitening
void MRI_DATA::init_noise_samples(int total_samples){
	noise_samples.setStorage( ColumnMajorArray<2>());
	noise_samples.resize(total_samples,Num_Coils);
	noise_samples = complex<float>(0,0);
}

void MRI_DATA::clone_attributes( MRI_DATA &data){

	Num_Encodings= data.Num_Encodings;;
	Num_Readouts=data.Num_Readouts;
	Num_Slices=data.Num_Slices;
	Num_Pts= data.Num_Pts;
	Num_Coils=data.Num_Coils;
	
	trajectory_dims=data.trajectory_dims;
	trajectory_type=data.trajectory_type;
	
	xres = data.xres;
	yres = data.yres;
	zres = data.zres;
	tres = data.tres;
	
	sms_factor = data.sms_factor;
	sms_type = data.sms_type;
	
	zfov = data.xfov;
	yfov = data.yfov;
	zfov = data.zfov;
}



// Constructer for MRI data type
void MRI_DATA::init_memory(void){
	
	
	cout << "Container Size = " << Num_Encodings << " x " << Num_Coils << endl;
	cout <<"3D Size " << Num_Pts << " x " << Num_Readouts << " x " << Num_Slices << endl;

	{
		Array< Array< float,3>,1> temp = Alloc4DContainer<float>(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings);
		kx.reference(temp);
	}

	{
		Array< Array< float,3>,1> temp = Alloc4DContainer<float>(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings);
		ky.reference(temp);
	}

	{
		Array< Array< float,3>,1> temp = Alloc4DContainer<float>(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings);
		kz.reference(temp);
	}

	{
		Array< Array< float,3>,2> temp = Alloc5DContainer<float>(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings,sms_factor);
		z.reference(temp);
	}

	{
		Array< Array< float,3>,1> temp = Alloc4DContainer<float>(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings);
		kw.reference(temp);
	}
	
	{
		Array< Array< float,3>,1> temp = Alloc4DContainer<float>(Num_Pts,Num_Readouts,Num_Slices,Num_Encodings);
		kt.reference(temp);
	}
		
	{
		Array< Array< complex<float>,3>,2> temp = Alloc5DContainer<complex<float> >( Num_Pts,Num_Readouts,Num_Slices,Num_Encodings,Num_Coils);
		kdata.reference(temp);
	} 
	
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

void MRI_DATA::init_gating_kdata(int gating_samples){
	kdata_gating.setStorage( ColumnMajorArray<5>());
	kdata_gating.resize(gating_samples,Num_Readouts,Num_Slices,Num_Encodings,Num_Coils);
}


MRI_DATA::MRI_DATA( void){
	Num_Encodings = -1;
	Num_Readouts = -1;
	Num_Pts = -1;
	Num_Coils = -1;
	Num_Slices =  1; 
	trajectory_dims = THREED;
	trajectory_type = THREEDNONCARTESIAN;
	sms_type = SMSoff;
	sms_factor = 1;
}

//---------------------------------------------------
//    This function allocates and reads all data into memory
//---------------------------------------------------

void MRI_DATA::demod_kdata( float demod){
	for(int c=0; c<Num_Coils;c++){
		for(int e=0; e<Num_Encodings; e++){
	 		for(int slice =0; slice< Num_Slices; slice++){
  				#pragma omp parallel for
  				for(int readout =0; readout< Num_Readouts; readout++){
  					for(int i =0; i< Num_Pts; i++){
  						kdata(e,c)(i,readout,slice) *= polar<float>(1.0,demod*2.0*arma::datum::pi*kt(e)(i,readout,slice));
				}}
	}}}
}

//---------------------------------------------------
//  Temporary Function to Write Data ( will be replaced by ismrmd ) 
//---------------------------------------------------

void MRI_DATA::write_external_data( const char *fname){
		
	HDF5 file = HDF5(fname,"w");
	
	cout << "Exporting Attributes" << endl;
	
	// Add dimensions
	file.AddH5Scaler("Kdata","Num_Encodings",Num_Encodings);
	file.AddH5Scaler("Kdata","Num_Readouts",Num_Readouts);
	file.AddH5Scaler("Kdata","Num_Slices",Num_Slices);
	file.AddH5Scaler("Kdata","Num_Pts",Num_Pts);
	file.AddH5Scaler("Kdata","Num_Coils",Num_Coils);
	
	// 2D/3D Cartesian/Non-Cartesian
	int temp1 = (int)trajectory_dims;
	file.AddH5Scaler("Kdata","trajectory_dims",temp1);
	
	int temp2 = (int)trajectory_type;
	file.AddH5Scaler("Kdata","trajectory_type",temp2);
	
	for(int encode = 0; encode < kdata.length(firstDim); encode++){
		
		cout << "Exporting " << encode << endl;
		
		{
			stringstream ss;
  			ss << "KT_E" << encode;	
			string s = ss.str();
			file.AddH5Array( "Kdata",s.c_str(),kt(encode));	
		}
		
	
		{
			stringstream ss;
  			ss << "KX_E" << encode;	
			string s = ss.str();
			file.AddH5Array( "Kdata",s.c_str(),kx(encode));	
		}	
				
		{
			stringstream ss;
  			ss << "KY_E" << encode;	
			string s = ss.str();
			file.AddH5Array( "Kdata",s.c_str(),ky(encode));	
		}	
			
		{
			stringstream ss;
  			ss << "KZ_E" << encode;	
			string s = ss.str();
			file.AddH5Array( "Kdata",s.c_str(),kz(encode));	
		}		
	
		if(sms_type == SMSon){
			for( int sms_pos=0; sms_pos < sms_factor; sms_pos++){
				stringstream ss;
  				ss << "Z_E" << encode << "_S" << sms_pos;	
				string s = ss.str();
				file.AddH5Array( "Kdata",s.c_str(),z(encode,sms_pos));
			}	
		}
	
		{
			stringstream ss;
  			ss << "KW_E" << encode;	
			string s = ss.str();
			file.AddH5Array( "Kdata",s.c_str(),kw(encode));	
		}	
	
		for(int coil = 0; coil < kdata.length(secondDim); coil++){
			{
				stringstream ss;
  				ss << "KData_E" << encode << "_C" << coil;	
				string s = ss.str();
				file.AddH5Array( "Kdata",s.c_str(),kdata(encode,coil));	
			}	
			
		}
	}
	
	// Noise Samples
	if( noise_samples.numElements()!=0){
		file.AddH5Array( "Kdata","Noise",noise_samples);	
	}
	
	// Gating
	if(ecg.numElements() != 0){
		file.AddH5Array( "Gating","ecg",ecg);	
		file.AddH5Array( "Gating","resp",resp);	
		file.AddH5Array( "Gating","prep",prep);	
		file.AddH5Array( "Gating","time",time);	
	}
	
	if(kdata_gating.numElements() != 0){
		file.AddH5Array( "Gating","kdata_gating",kdata_gating);	
	}
}

void MRI_DATA::read_external_data( const char *fname){
		
	HDF5 file = HDF5(fname,"r");
	
	cout << "Reading External File " << fname  << endl;
	
	// Read atrributes common to all
	file.ReadH5Scaler("Kdata","Num_Encodings",&Num_Encodings);
	file.ReadH5Scaler("Kdata","Num_Readouts",&Num_Readouts);
	file.ReadH5Scaler("Kdata","Num_Slices",&Num_Slices);
	file.ReadH5Scaler("Kdata","Num_Pts",&Num_Pts);
	file.ReadH5Scaler("Kdata","Num_Coils",&Num_Coils);
	
	init_memory():
		
	// 2D/3D Cartesian/Non-Cartesian
	int temp1;
	file.ReadH5Scaler("Kdata","trajectory_dims",&temp1);
	trajectory_dims = static_cast<TrajDim>(temp1);
	
	file.ReadH5Scaler("Kdata","trajectory_type",&temp1);
	trajectory_type = static_cast<TrajType>(temp1);
		
	for(int encode = 0; encode < kdata.length(firstDim); encode++){
		
		cout << "Importing " << encode << endl;
		
		{
			stringstream ss;
  			ss << "KT_E" << encode;	
			string s = ss.str();
			file.ReadH5Array( "Kdata", s, kt(encode));	
		}
		
	
		{
			stringstream ss;
  			ss << "KX_E" << encode;	
			string s = ss.str();
			file.ReadH5Array( "Kdata",s, kx(encode));	
		}	
				
		{
			stringstream ss;
  			ss << "KY_E" << encode;	
			string s = ss.str();
			file.ReadH5Array( "Kdata",s, ky(encode));	
		}	
			
		{
			stringstream ss;
  			ss << "KZ_E" << encode;	
			string s = ss.str();
			file.ReadH5Array( "Kdata",s, kz(encode));	
		}		
	
		if(sms_type == SMSon){
			for( int sms_pos=0; sms_pos < sms_factor; sms_pos++){
				stringstream ss;
  				ss << "Z_E" << encode << "_S" << sms_pos;	
				string s = ss.str();
				file.ReadH5Array( "Kdata",s,z(encode,sms_pos));
			}	
		}
	
		{
			stringstream ss;
  			ss << "KW_E" << encode;	
			string s = ss.str();
			file.ReadH5Array( "Kdata",s,kw(encode));	
		}	
	
		for(int coil = 0; coil < kdata.length(secondDim); coil++){
			{
				stringstream ss;
  				ss << "KData_E" << encode << "_C" << coil;	
				string s = ss.str();
				file.ReadH5Array( "Kdata",s, kdata(encode,coil));	
			}	
			
		}
	}
	
	// Noise Samples
	file.ReadH5Array( "Kdata","Noise",noise_samples);	
	
	// Gating
	file.ReadH5Array( "Gating","ecg",ecg);	
	file.ReadH5Array( "Gating","resp",resp);	
	file.ReadH5Array( "Gating","prep",prep);	
	file.ReadH5Array( "Gating","time",time);	
	
	file.ReadH5Array( "Gating","kdata_gating",kdata_gating);	
	
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
	for(int e = 0; e< Num_Encodings; e++){
		
		cout<< "\t Kx" << endl;
		{
		Array<float,3> temp_us = kx(e);
		Array<float,3> temp = temp_us(Range(fromStart,toEnd,1),Range(fromStart,toEnd,us),Range(fromStart,toEnd,1));
		kx(e).resize(Num_Pts,Num_Readouts_us,Num_Slices);
		kx(e) = temp;
		}
		
		{	
		cout<< "\t Ky" << endl;
		Array<float,3> temp_us = ky(e);
		Array<float,3> temp = temp_us(Range(fromStart,toEnd,1),Range(fromStart,toEnd,us),Range(fromStart,toEnd,1));
		ky(e).resize(Num_Pts,Num_Readouts_us,Num_Slices);
		ky(e) = temp;
		}
		
		{
		cout<< "\t Kz" << endl;
		Array<float,3> temp_us = kz(e);
		Array<float,3> temp = temp_us(Range(fromStart,toEnd,1),Range(fromStart,toEnd,us),Range(fromStart,toEnd,1));
		kz(e).resize(Num_Pts,Num_Readouts_us,Num_Slices);
		kz(e) = temp;
		}
		
		
		{
		cout<< "\t Kw" << endl;
		Array<float,3> temp_us = kw(e);
		Array<float,3> temp = temp_us(Range(fromStart,toEnd,1),Range(fromStart,toEnd,us),Range(fromStart,toEnd,1));
		kw(e).resize(Num_Pts,Num_Readouts_us,Num_Slices);
		kw(e) = temp;
		}
		
		for( int coil=0; coil < Num_Coils; coil++){
			{
			cout<< "\t Kdata" << endl;
			Array<complex<float>,3> temp_data_us = kdata(e,coil);
			Array<complex<float>,3> temp_data = temp_data_us(Range(fromStart,toEnd,1),Range(fromStart,toEnd,us),Range::all());
			kdata(e,coil).resize(Num_Pts,Num_Readouts_us,Num_Slices);
			kdata(e,coil) = temp_data;
			}
		}
	
	}
	}
}


/** Coil compress data with a cutoff of thresh*max(SV)
*
*/
void MRI_DATA::coilcompress(float thresh)
{ 

  int Num_Pixels = Num_Encodings*kdata(0,0).numElements();
      
  arma::cx_fmat all_data;
  all_data.zeros(Num_Pixels,Num_Coils);
  cout << "Num_pixels = " << Num_Pixels << endl;
  
  cout << "Collect Data" << endl;
  for(int coil=0; coil< Num_Coils; coil++){
  
  	int pos = 0;
  	for(int e =0; e< Num_Encodings; e++){
  	for(int slice =0; slice< Num_Slices; slice++){
  	for(int readout =0; readout< Num_Readouts; readout++){
  	for(int i =0; i< Num_Pts; i++){
  		all_data(pos,coil) = kdata(e,coil)(i,readout,slice);
		pos++;
  	}}}}
	cout << "Copied pixels = " << pos << endl;
  }
  
  cout << "SVD " << endl << flush;
  arma::fvec s;
  arma::cx_fmat U;
  arma::cx_fmat V;
  arma::svd_econ(U,s,V,all_data); 
  s = s/s(0);
  
  s.print("S");
  
  arma::cx_fmat VV = V.cols(0,(int)thresh-1);
  VV.print("V");    
  cout << "Rotate " << endl << flush;  
  all_data = all_data*VV;
  
  cout << "Resize to " << thresh << endl << flush;
  Array< Array<complex<float>,3>,2>kdata2(kdata.length(firstDim),(int)thresh,ColumnMajorArray<2>());
  for(int e =0; e< Num_Encodings; e++){
  for(int c =0; c< thresh; c++){
  	kdata2(e,c).reference(kdata(e,c));
  }}
  kdata.reference(kdata2); 
  Num_Coils = (int)thresh;
  
  cout << "Copy Back" << endl << flush; 
  for(int coil=0; coil< Num_Coils; coil++){
  	int pos = 0;
  	for(int e =0; e< Num_Encodings; e++){
  	for(int slice =0; slice< Num_Slices; slice++){
  	for(int readout =0; readout< Num_Readouts; readout++){
  	for(int i =0; i< Num_Pts; i++){
  		kdata(e,coil)(i,readout,slice) = all_data(pos,coil);
		pos++;
  	}}}}
	cout << "Copied pixels = " << pos << endl;
  }
  
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
	arma::cx_fmat NoiseData = arma::randu<arma::cx_fmat>(noise_samples.length(secondDim),noise_samples.length(firstDim));
  	for(int coil=0; coil < Num_Coils; coil++){
	for(int i=0; i< noise_samples.length(firstDim); i++){
		NoiseData(coil,i) = noise_samples(i,coil);
	}}
	
	
	cout << "Calc Cov" << endl << flush;
	arma::cx_fmat CV =  NoiseData*NoiseData.t();
	CV.save("CovMatrix.dat",arma::raw_binary);
		
	cout << "Whiten" << endl;
	arma::cx_fmat V = chol(CV);
	arma::cx_fmat VT = V.t(); 
	arma::cx_fmat Decorr = VT.i();
			
	// Test Whitening
	arma::cx_fmat W = NoiseData;
	arma::cx_fmat temp = arma::randu<arma::cx_fmat>(Num_Coils);
	for(int i =0; i< noise_samples.length(firstDim); i++){
		
		for(int coil=0; coil < Num_Coils; coil++){
			temp(coil,0)=W(coil,i);
		}
		arma::cx_fmat temp2 = Decorr*temp;
		
		for(int coil=0; coil < Num_Coils; coil++){
			W(coil,i)=temp2(coil,0);
		}
	}
			
	arma::cx_fmat CV_POST = W*W.t();
	CV_POST.save("CovMatrixPost.dat", arma::raw_binary);
	
	
	// Now Whiten Actual Data
	cout << "Whiten all data" << endl;
	for(int e =0; e< Num_Encodings; e++){
  	for(int k =0; k< Num_Slices; k++){
  		#pragma omp parallel for
  		for(int j =0; j< Num_Readouts; j++){
  	
			arma::cx_fmat AA;
  			AA.zeros( Num_Coils,Num_Pts); // Working memory
    	
			// Copy into memory
			for(int coil = 0; coil < Num_Coils; coil++) {
			for(int i =0; i< Num_Pts; i++){
				AA(coil,i) =  kdata(e,coil)(i,j,k);
			}
			}
	
			// Transform to matrix
			arma::cx_fmat temp2 = Decorr*AA; 
		
			// Copy Back
			for(int coil = 0; coil < Num_Coils; coil++) {
 			for(int i =0; i< Num_Pts; i++){
				kdata(e,coil)(i,j,k) = temp2(coil,i);
			}}
  		}
  
  	}}
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



