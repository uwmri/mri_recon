#include "mri_data.h"

using arma::cx_fmat;
using arma::fvec;
using namespace NDarray;

void MRI_DATA::dump_stats(const string name, const Array< Array<float, 3>, 1> & in) {
	cout << name << endl;
	cout << "\tContainer Size = " << in.length(firstDim) << endl;
	cout << "\tElement size = " << in(0).length(firstDim) << " x " << in(0).length(secondDim) << " x " << in(0).length(thirdDim) << endl;
	float high = 0;
	float low = 9999;
	for (Array< Array<float, 3>, 1> ::const_iterator miter = in.begin(); miter != in.end(); miter++) {
		high = max(high, max(*miter));
		low = min(low, min(*miter));
	}
	cout << "\tRange = " << low << " to " << high << endl;
}

void MRI_DATA::dump_stats(const string name, const Array< Array<complex<float>, 3>, 2> & in) {
	cout << name << endl;
	cout << "\tContainer Size = " << in.length(firstDim) << " x " << in.length(secondDim) << endl;
	cout << "\tElement size = " << in(0).length(firstDim) << " x " << in(0).length(secondDim) << " x " << in(0).length(thirdDim) << endl;
	float high = 0;
	float low = 9999;
	for (Array< Array< complex<float>, 3>, 2>::const_iterator miter = in.begin(); miter != in.end(); miter++) {
		float temp = max(abs(*miter));
		high = max(high, temp);

		temp = min(abs(*miter));
		low = min(low, temp);
	}
	cout << "\tRange = " << low << " to " << high << endl;
}

void MRI_DATA::scale_fov(float scale_x, float scale_y, float scale_z) {

	if ((scale_x == 1.0) && (scale_y == 1.0) && (scale_z == 1.0)) {
		return;
	}

	cout << "Scaling the fov by" << scale_x << "," << scale_y << "," << scale_z << endl;
	float scale_kx = ( dft_needed(0) ) ? ( scale_x ) : ( scale_x );
	float scale_ky = ( dft_needed(1) ) ? ( scale_y ) : ( scale_y );
	float scale_kz = ( dft_needed(2) ) ? ( scale_z ) : ( scale_z );
	cout << "Scaling the kspace by" << scale_kx << "," << scale_ky << "," << scale_kz << endl;
		
	// Multiply kspace by inverse
	for (Array< Array<float, 3>, 1>::iterator miter = kx.begin(); miter != kx.end(); miter++) {
		(*miter) *= scale_kx;
	}

	for (Array< Array<float, 3>, 1>::iterator miter = ky.begin(); miter != ky.end(); miter++) {
		(*miter) *= scale_ky;
	}

	for (Array< Array<float, 3>, 1>::iterator miter = kz.begin(); miter != kz.end(); miter++) {
		(*miter) *= scale_kz;
	}

	// Multiple image by scale
	zfov *= scale_z;
	yfov *= scale_y;
	xfov *= scale_x;

	if (sms_type == SMSon) {
		for (Array< Array<float, 3>, 2>::iterator miter = z.begin(); miter != z.end(); miter++) {
			(*miter) /= scale_z;
		}
	}

}

void MRI_DATA::convert_encodes_to_coils( int passes ){
	
	// This function treats encodes as sperate coils
	// Only works if encodes are identical
	int new_Num_Coils = (Num_Encodings/passes)*Num_Coils;
	// int new_Num_Encodings =  passes;
	/*
	this->kx.resizeAndPreserve(1);
	this->ky.resizeAndPreserve(1);
	this->kz.resizeAndPreserve(1);
	this->kw.resizeAndPreserve(1);
	this->kt.resizeAndPreserve(1);
	*/
	
	cout << "Create temp array" << endl;
	Array< Array<complex<float>,3>, 2> kdata_temp;
	kdata_temp.setStorage(ColumnMajorArray<2>());
	kdata_temp.resize(1, new_Num_Coils);
			
	cout << "Copy array" << endl;	
	int count= 0;
	for (int ee = 0; ee < this->Num_Encodings; ee ++) {
		for (int coil = 0; coil < this->Num_Coils; coil++) {
			kdata_temp(0,count).setStorage(ColumnMajorArray<3>());
			kdata_temp(0,count).resize(this->kdata(ee,coil).shape());
			kdata_temp(0,count) = this->kdata(ee,coil);
			count++;
		}
	}

	cout << "Copy Back" << endl;
	this->kdata.resize(1,new_Num_Coils);
	for( int pos=0; pos <  new_Num_Coils; pos++){
			this->kdata(0,pos).setStorage(ColumnMajorArray<3>());
			this->kdata(0,pos).resize(kdata_temp(0,pos).shape());
			this->kdata(0,pos) = kdata_temp(0,pos);
	}
	
	// Reset size
	this->Num_Encodings = 1;
	this->Num_Coils = new_Num_Coils;
	
}

MRI_DATA MRI_DATA::subframe(int eStart, int eStop, int eStride) {

	MRI_DATA data2;
	data2.clone_attributes(*this);

	// Copy info
	data2.Num_Encodings = floor((1 + eStop - eStart) / eStride);
	data2.init_memory();

	int count = 0;
	for (int ee = eStart; ee <= eStop; ee += eStride) {

		data2.kx(count) = this->kx(ee);
		data2.ky(count) = this->ky(ee);
		data2.kz(count) = this->kz(ee);
		data2.kw(count) = this->kw(ee);
		data2.kt(count) = this->kt(ee);

		// SMS need z positions
		if (data2.sms_type == SMSon) {
			for (int sms_pos = 0; sms_pos < sms_factor; sms_pos++) {
				data2.z(sms_pos, count) = this->z(sms_pos, ee);
			}
		}

		// Copy the data	
		for (int coil = 0; coil < data2.Num_Coils; coil++) {
			data2.kdata(count, coil) = this->kdata(ee, coil);
		}
		count = count + 1;
	}

	return(data2);
}

void MRI_DATA::stats(void) {

	if (kx.numElements() == 0) {
		cout << "Kspace does not exist" << endl;
	}
	else {
		dump_stats("Kx", kx);
		dump_stats("Ky", ky);
		dump_stats("Kz", kz);
		dump_stats("Kw", kw);
	}


	if (kdata.numElements() == 0) {
		cout << "Kdata does not exist yet" << endl;
	}
	else {
		dump_stats("Kdata", kdata);
	}

	// gating
	if (ecg.numElements() == 0) {
		cout << "Physiologic data does not exist yet" << endl;
	}
	else {
		cout << "Range ECG = " << Dmin(ecg) << " to " << Dmax(ecg) << endl;
		cout << "Range RESP = " << Dmin(resp) << " to " << Dmax(resp) << endl;
		cout << "Range TIME = " << Dmin(time) << " to " << Dmax(time) << endl;
		cout << "Range PREP = " << Dmin(prep) << " to " << Dmax(prep) << endl;
	}

}

// Data for Whitening
void MRI_DATA::init_noise_samples(int total_samples) {
	noise_samples.setStorage(ColumnMajorArray<2>());
	noise_samples.resize(total_samples, Num_Coils);
	noise_samples = complex<float>(0, 0);
}

void MRI_DATA::clone_attributes(MRI_DATA &data) {

	Num_Encodings = data.Num_Encodings;;
	Num_Coils = data.Num_Coils;

	dft_needed = data.dft_needed;
	trajectory_type = data.trajectory_type;

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
void MRI_DATA::init_memory(void) {

	cout << "Container Size = " << Num_Encodings << " x " << Num_Coils << endl;
	kx.resize(Num_Encodings);
	ky.resize(Num_Encodings);
	kz.resize(Num_Encodings);
	kw.resize(Num_Encodings);
	kt.resize(Num_Encodings);

	if (sms_type == SMSon) {
		z.resize(Num_Encodings, sms_factor);
	}

	kdata.setStorage(ColumnMajorArray<2>());
	kdata.resize(Num_Encodings, Num_Coils);

	// Times 
	time.resize(Num_Encodings);
	ecg.resize(Num_Encodings);
	resp.resize(Num_Encodings);
	prep.resize(Num_Encodings);

}

void MRI_DATA::init_memory(int readouts, int shots, int slices) {
	cout << "Initializing Container for " << readouts << " x " << shots << " x " << slices << endl;
	init_memory();
	for (int e = 0; e < Num_Encodings; e++) {
		init_encode(e, readouts, shots, slices);
	}
}


void MRI_DATA::init_encode(int e, int readouts, int shots, int slices) {

	kx(e).setStorage(ColumnMajorArray<3>());
	kx(e).resize(readouts, shots, slices);

	ky(e).setStorage(ColumnMajorArray<3>());
	ky(e).resize(readouts, shots, slices);

	kz(e).setStorage(ColumnMajorArray<3>());
	kz(e).resize(readouts, shots, slices);

	if (sms_type == SMSon) {
		for (int sms_pos = 0; sms_pos < sms_factor; sms_pos++) {
			z(e, sms_pos).setStorage(ColumnMajorArray<3>());
			z(e, sms_pos).resize(readouts, shots, slices);
		}
	}

	kw(e).setStorage(ColumnMajorArray<3>());
	kw(e).resize(readouts, shots, slices);

	kt(e).setStorage(ColumnMajorArray<3>());
	kt(e).resize(readouts, shots, slices);

	for (int coil = 0; coil < Num_Coils; coil++) {
		kdata(e, coil).setStorage(ColumnMajorArray<3>());
		kdata(e, coil).resize(readouts, shots, slices);
	}

	time(e).setStorage(ColumnMajorArray<2>());
	time(e).resize(shots, slices);

	ecg(e).setStorage(ColumnMajorArray<2>());
	ecg(e).resize(shots, slices);

	resp(e).setStorage(ColumnMajorArray<2>());
	resp(e).resize(shots, slices);

	prep(e).setStorage(ColumnMajorArray<2>());
	prep(e).resize(shots, slices);

}


void MRI_DATA::init_gating_kdata(int gating_samples) {
	
	cout << "Gating samples = " << gating_samples << endl;
	/*kdata_gating.setStorage( ColumnMajorArray<5>());
	kdata_gating.resize(gating_samples,Num_Readouts,Num_Slices,Num_Encodings,Num_Coils);
	*/
}


MRI_DATA::MRI_DATA(void) {

	// Initial Values
	Num_Encodings = -1;
	Num_Coils = -1;
	for(int dir=0; dir < 3; dir++){
		trajectory_type(dir) = NONCARTESIAN;
		dft_needed(dir) = true;
	}
	sms_type = SMSoff;
	sms_factor = 1;
}

//---------------------------------------------------
//    This function allocates and reads all data into memory
//---------------------------------------------------

void MRI_DATA::demod_kdata(float demod) {

	for (int c = 0; c < Num_Coils; c++) {
		for (int e = 0; e < Num_Encodings; e++) {

			// Each Dataset
			for (int slice = 0; slice < kdata(e, c).length(thirdDim); slice++) {

#pragma omp parallel for
				for (int readout = 0; readout < kdata(e, c).length(secondDim); readout++) {
					for (int i = 0; i < kdata(e, c).length(firstDim); i++) {
						kdata(e, c)(i, readout, slice) *= polar<float>(1.0, demod*2.0*arma::datum::pi*kt(e)(i, readout, slice));
					}
				}
			}
		}
	}
}

//---------------------------------------------------
//  Temporary Function to Write Data ( will be replaced by ismrmd ) 
//---------------------------------------------------

void MRI_DATA::write_external_data(const char *fname) {

	HDF5 file = HDF5(fname, "w");

	cout << "Exporting Attributes" << endl;

	// Add dimensions
	file.AddH5Scaler("Kdata", "Num_Encodings", Num_Encodings);
	file.AddH5Scaler("Kdata", "Num_Coils", Num_Coils);

	// 2D/3D Cartesian/Non-Cartesian
	file.AddH5Scaler("Kdata", "trajectory_typeX", (int)trajectory_type(0) );
	file.AddH5Scaler("Kdata", "trajectory_typeY", (int)trajectory_type(1) );
	file.AddH5Scaler("Kdata", "trajectory_typeZ", (int)trajectory_type(2) );

	file.AddH5Scaler("Kdata", "dft_neededX", (int)dft_needed(0) );
	file.AddH5Scaler("Kdata", "dft_neededY", (int)dft_needed(1) );
	file.AddH5Scaler("Kdata", "dft_neededZ", (int)dft_needed(2) );

	for (int encode = 0; encode < kdata.length(firstDim); encode++) {

		cout << "Exporting " << encode << endl;
		{
			try{
				stringstream ss;
				ss << "KT_E" << encode;
				string s = ss.str();
				file.AddH5Array("Kdata", s.c_str(), kt(encode));
			}catch(...){
				cout << "Can't export KT for encode " << encode << endl;
			}
		}


		{
			try{
				stringstream ss;
				ss << "KX_E" << encode;
				string s = ss.str();
				file.AddH5Array("Kdata", s.c_str(), kx(encode));
			}catch(...){
				cout << "Can't export KX for encode " << encode << endl;
			}
		}

		{
			try{
				stringstream ss;
				ss << "KY_E" << encode;
				string s = ss.str();
				file.AddH5Array("Kdata", s.c_str(), ky(encode));
			}catch(...){
				cout << "Can't export KY for encode " << encode << endl;
			}
		}

		{
			try{
				stringstream ss;
				ss << "KZ_E" << encode;
				string s = ss.str();
				file.AddH5Array("Kdata", s.c_str(), kz(encode));
			}catch(...){
				cout << "Can't export KZ for encode " << encode << endl;
			}
			
		}

		if (sms_type == SMSon) {
			cout << "Export SMS Z2" << endl << flush;
			for (int sms_pos = 0; sms_pos < sms_factor; sms_pos++) {
				stringstream ss;
				ss << "Z_E" << encode << "_S" << sms_pos;
				string s = ss.str();
				file.AddH5Array("Kdata", s.c_str(), z(encode, sms_pos));
			}
		}

		{
			try{
				stringstream ss;
				ss << "KW_E" << encode;
				string s = ss.str();
				file.AddH5Array("Kdata", s.c_str(), kw(encode));
			}catch(...){
				cout << "Can't export KW for encode " << encode << endl;
			}
		}
		
		cout << "Exporting data" << endl << flush;
		for (int coil = 0; coil < kdata.length(secondDim); coil++) {
			{
				try{
					stringstream ss;
					ss << "KData_E" << encode << "_C" << coil;
					string s = ss.str();
					file.AddH5Array("Kdata", s.c_str(), kdata(encode, coil));
				}catch(...){
					cout << "Can't export Kdata for encode " << encode << ", coil " << coil << endl;
				}
			}

		}

		// Gating
		if (ecg.numElements() != 0) {
			
			cout << "Exporting Gating " << endl << flush;
			
			{
				try{ 
					stringstream ss;
					ss << "ECG_E" << encode;
					string s = ss.str();
					file.AddH5Array("Gating", s.c_str(), ecg(encode));
				}catch(...){
					cout << "Can't export ECG data" << endl;
				}
			}

			{
				try{
				stringstream ss;
				ss << "RESP_E" << encode;
				string s = ss.str();
				file.AddH5Array("Gating", s.c_str(), resp(encode));
				}catch(...){
					cout << "Can't export Resp data" << endl;
				}

			}


			{
				try{
				stringstream ss;
				ss << "PREP_E" << encode;
				string s = ss.str();
				file.AddH5Array("Gating", s.c_str(), prep(encode));
				}catch(...){
					cout << "Can't export PREP data" << endl;
				}
			}

			{
				try{
				stringstream ss;
				ss << "TIME_E" << encode;
				string s = ss.str();
				file.AddH5Array("Gating", s.c_str(), time(encode));
				}catch(...){
						cout << "Can't export TIME data" << endl;
			
				}
			}

		}
	}

	// Noise Samples
	if (noise_samples.numElements() != 0) {
		file.AddH5Array("Kdata", "Noise", noise_samples);
	}

	if (kdata_gating.numElements() != 0) {
		// TEMP file.AddH5Array( "Gating","kdata_gating",kdata_gating);	
	}
}

void MRI_DATA::read_external_data(const char *fname) {

	HDF5 file = HDF5(fname, "r");

	cout << "Reading External File " << fname << endl;

	// Read atrributes common to all
	file.ReadH5Scaler("Kdata", "Num_Encodings", &Num_Encodings);
	file.ReadH5Scaler("Kdata", "Num_Coils", &Num_Coils);

	init_memory();


	// 2D/3D Cartesian/Non-Cartesian
	int temp1;
	file.ReadH5Scaler("Kdata", "trajectory_typeX", &temp1);
	trajectory_type(0) = static_cast<MRI_DATA::TrajType>(temp1);

	file.ReadH5Scaler("Kdata", "trajectory_typeY", &temp1);
	trajectory_type(1) = static_cast<MRI_DATA::TrajType>(temp1);

	file.ReadH5Scaler("Kdata", "trajectory_typeZ", &temp1);
	trajectory_type(2) = static_cast<MRI_DATA::TrajType>(temp1);


	file.ReadH5Scaler("Kdata", "dft_neededX", &temp1);
	dft_needed(0) = static_cast<bool>(temp1);

	file.ReadH5Scaler("Kdata", "dft_neededY", &temp1);
	dft_needed(1) = static_cast<bool>(temp1);

	file.ReadH5Scaler("Kdata", "dft_neededZ", &temp1);
	dft_needed(2) = static_cast<bool>(temp1);

	for (int encode = 0; encode < kdata.length(firstDim); encode++) {

		cout << "Importing " << encode << endl;

		{
			stringstream ss;
			ss << "KT_E" << encode;
			string s = ss.str();
			file.ReadH5Array("Kdata", s, kt(encode));
		}


		{
			stringstream ss;
			ss << "KX_E" << encode;
			string s = ss.str();
			file.ReadH5Array("Kdata", s, kx(encode));
		}

		{
			stringstream ss;
			ss << "KY_E" << encode;
			string s = ss.str();
			file.ReadH5Array("Kdata", s, ky(encode));
		}

		{
			stringstream ss;
			ss << "KZ_E" << encode;
			string s = ss.str();
			file.ReadH5Array("Kdata", s, kz(encode));
		}

		if (sms_type == SMSon) {
			for (int sms_pos = 0; sms_pos < sms_factor; sms_pos++) {
				stringstream ss;
				ss << "Z_E" << encode << "_S" << sms_pos;
				string s = ss.str();
				file.ReadH5Array("Kdata", s, z(encode, sms_pos));
			}
		}

		{
			stringstream ss;
			ss << "KW_E" << encode;
			string s = ss.str();
			file.ReadH5Array("Kdata", s, kw(encode));
		}

		for (int coil = 0; coil < kdata.length(secondDim); coil++) {
			{
				stringstream ss;
				ss << "KData_E" << encode << "_C" << coil;
				string s = ss.str();
				file.ReadH5Array("Kdata", s, kdata(encode, coil));
			}
		}

		cout << "Read Gating " << endl << flush;
		{
			stringstream ss;
			ss << "ECG_E" << encode;
			string s = ss.str();
			file.ReadH5Array("Gating", s.c_str(), ecg(encode));
		}

		{
			stringstream ss;
			ss << "RESP_E" << encode;
			string s = ss.str();
			file.ReadH5Array("Gating", s.c_str(), resp(encode));
		}


		{
			stringstream ss;
			ss << "PREP_E" << encode;
			string s = ss.str();
			file.ReadH5Array("Gating", s.c_str(), prep(encode));
		}

		{
			stringstream ss;
			ss << "TIME_E" << encode;
			string s = ss.str();
			file.ReadH5Array("Gating", s.c_str(), time(encode));
		}
	}

	// Noise Samples
	cout << "Read Noise Samples" << endl << flush;
	file.ReadH5Array("Kdata", "Noise", noise_samples);

	//file.ReadH5Array( "Gating","kdata_gating",kdata_gating);	

}


/** Coil compress data with a cutoff of thresh*max(SV)
*
*/
void MRI_DATA::coilcompress(float thresh)
{

	int Num_Pixels = 0;
	for (Array< Array<complex<float>, 3>, 2>::const_iterator miter = kdata.begin(); miter != kdata.end(); miter++) {
		Num_Pixels += (*miter).numElements();
	}

	arma::cx_fmat all_data;
	all_data.zeros(Num_Pixels, Num_Coils);
	cout << "Num_pixels = " << Num_Pixels << endl;

	cout << "Collect Data" << endl;
	for (int coil = 0; coil < Num_Coils; coil++) {

		int pos = 0;
		for (int e = 0; e < Num_Encodings; e++) {
			for (Array<complex<float>, 3>::const_iterator miter = kdata(e, coil).begin(); miter != kdata(e, coil).end(); miter++) {
				all_data(pos, coil) = (*miter);
				pos++;
			}
		}
		cout << "Copied pixels = " << pos << endl;
	}

	cout << "SVD " << endl << flush;
	arma::fvec s;
	arma::cx_fmat U;
	arma::cx_fmat V;
	arma::svd_econ(U, s, V, all_data);
	s = s / s(0);

	s.print("S");

	arma::cx_fmat VV = V.cols(0, (int)thresh - 1);
	VV.print("V");
	cout << "Rotate " << endl << flush;
	all_data = all_data*VV;

	cout << "Resize to " << thresh << endl << flush;
	Array< Array<complex<float>, 3>, 2>kdata2(kdata.length(firstDim), (int)thresh, ColumnMajorArray<2>());
	for (int e = 0; e < Num_Encodings; e++) {
		for (int c = 0; c < thresh; c++) {
			kdata2(e, c).reference(kdata(e, c));
		}
	}
	kdata.reference(kdata2);
	Num_Coils = (int)thresh;

	cout << "Copy Back" << endl << flush;
	for (int coil = 0; coil < Num_Coils; coil++) {
		int pos = 0;
		for (int e = 0; e < Num_Encodings; e++) {
			for (Array< complex<float>, 3>::iterator miter = kdata(e, coil).begin(); miter != kdata(e, coil).end(); miter++) {
				(*miter) = all_data(pos, coil);
				pos++;
			}
		}
		cout << "Copied pixels = " << pos << endl;
	}

	cout << "done" << endl;

}

//--------------------------------------------------
//  Whiten Data
//--------------------------------------------------

void MRI_DATA::whiten(void) {

	if (noise_samples.numElements() == 0) {
		cout << "Noise Samples do not exist:: Can't whiten data" << endl;
		return;
	}

	cout << "Noise Samples : " << noise_samples.length(firstDim) << endl;
	cout << "Noise Pre-Whitening" << endl << flush;

	// Copy into matrix
	arma::cx_fmat NoiseData = arma::randu<arma::cx_fmat>(noise_samples.length(secondDim), noise_samples.length(firstDim));
	for (int coil = 0; coil < Num_Coils; coil++) {
		for (int i = 0; i < noise_samples.length(firstDim); i++) {
			NoiseData(coil, i) = noise_samples(i, coil);
		}
	}

	cout << "Calc Cov" << endl << flush;
	arma::cx_fmat CV = NoiseData*NoiseData.t();
	CV.save("CovMatrix.dat", arma::raw_binary);

	cout << "Whiten" << endl;
	arma::cx_fmat V = chol(CV);
	arma::cx_fmat VT = V.t();
	arma::cx_fmat Decorr = VT.i();

	// Test Whitening
	arma::cx_fmat W = NoiseData;
	arma::cx_fmat temp = arma::randu<arma::cx_fmat>(Num_Coils);
	for (int i = 0; i < noise_samples.length(firstDim); i++) {

		for (int coil = 0; coil < Num_Coils; coil++) {
			temp(coil, 0) = W(coil, i);
		}
		arma::cx_fmat temp2 = Decorr*temp;

		for (int coil = 0; coil < Num_Coils; coil++) {
			W(coil, i) = temp2(coil, 0);
		}
	}

	arma::cx_fmat CV_POST = W*W.t();
	CV_POST.save("CovMatrixPost.dat", arma::raw_binary);


	// Now Whiten Actual Data
	cout << "Whiten all data" << endl;
	for (int e = 0; e < Num_Encodings; e++) {

		for (int k = 0; k < kdata(e, 0).length(thirdDim); k++) {
#pragma omp parallel for
			for (int j = 0; j < kdata(e, 0).length(secondDim); j++) {

				arma::cx_fmat AA;
				AA.zeros(Num_Coils, kdata(e, 0).length(firstDim)); // Working memory

				// Copy into memory
				for (int coil = 0; coil < Num_Coils; coil++) {
					for (int i = 0; i < kdata(e, 0).length(firstDim); i++) {
						AA(coil, i) = kdata(e, coil)(i, j, k);
					}
				}

				// Transform to matrix
				arma::cx_fmat temp2 = Decorr*AA;

				// Copy Back
				for (int coil = 0; coil < Num_Coils; coil++) {
					for (int i = 0; i < kdata(e, 0).length(firstDim); i++) {
						kdata(e, coil)(i, j, k) = temp2(coil, i);
					}
				}
			}

		}
	}
}



