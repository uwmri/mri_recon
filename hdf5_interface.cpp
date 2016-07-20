#include "hdf5_interface.h"

using namespace std; 
using namespace NDarray;
using namespace H5;

HDF5::HDF5( string FileName, const string operation){
	
	cout << "HDF5 - Opening file operation = " << operation << endl; 
	
	// OPen the file
	if(operation.compare("r") == 0){
		file =  H5File(FileName,H5F_ACC_RDONLY);	
	}else if( (operation.compare("rw") ==0) || (operation.compare("rw")==0) ){
		file =  H5File(FileName,H5F_ACC_RDWR);	
	}else{
		file =  H5File(FileName,H5F_ACC_TRUNC);	
	}
}

HDF5::HDF5( void  ){

}

H5::Group HDF5::create_group( H5::H5File &act_file, string GroupName ){
	
	// Create Group
	Group group;
	try{
		//cout << "Trying to Open: " << GroupName << endl;
		group = Group(act_file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		//cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(act_file.createGroup( x) );
	}
	
	return(group);
}

H5::Group HDF5::open_group( H5::H5File &act_file, string GroupName ){
	
	// Open Group
	Group group;
	try{
		group = Group(act_file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		//cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		cout << "HDF5 Read : Group " << GroupName << " does not exist" << endl;
	}
	return(group);
}

int HDF5::ReadH5Scaler( string GroupName, string Name, int *out){
	
	Exception::dontPrint();
	
	// Open Group
	Group group = open_group(file,GroupName);
	
	// Now open attribute
	Attribute att;
	try{
		att = group.openAttribute( Name );
	}catch(FileIException error){
		cout << "HDF5 Read :  " << GroupName << " attribute " << Name << " does not exist" << endl;
	}
	
	// Get the data type
	DataType data_type = DataType(att.getDataType());
	DataType int_type(PredType::STD_I32LE); 
	
	if( int_type==data_type){
		att.read( data_type, out ); 
	}else{
		cout << "Expected to read little endian int32 but did not" << endl;
		exit(1);
	}

	return(0);
}

int HDF5::AddH5Scaler( string GroupName, string Name, int A){
	
	Exception::dontPrint();
	
	// Create Group
	Group group = create_group(file,GroupName);
	
	IntType int_type(PredType::STD_I32LE);
	DataSpace att_space(H5S_SCALAR);
	Attribute att = group.createAttribute( Name , int_type, att_space );
	att.write( int_type, &A );
		
	return(0);
}

int HDF5::ReadH5Scaler( string GroupName, string Name, float *out){
	
	Exception::dontPrint();
	
	// Open Group
	Group group = open_group(file,GroupName);
	
	// Now open attribute
	Attribute att;
	try{
		att = group.openAttribute( Name );
	}catch(FileIException error){
		cout << "HDF5 Read :  " << GroupName << " attribute " << Name << " does not exist" << endl;
	}
	
	// Get the data type
	DataType data_type = DataType(att.getDataType());
	DataType int_type(PredType::NATIVE_FLOAT); 
	
	if( int_type==data_type){
		att.read( data_type, out ); 
	}else{
		cout << "Expected to read little endian float32 but did not" << endl;
		exit(1);
	}

	return(0);
}


int HDF5::AddH5Scaler( string GroupName, string Name, float A){
	
	Exception::dontPrint();
	
	// Create Group
	Group group = create_group(file,GroupName);
		
	FloatType float_type(PredType::NATIVE_FLOAT);
	DataSpace att_space(H5S_SCALAR);
	Attribute att = group.createAttribute( Name , float_type, att_space );
	att.write( float_type, &A );
		
	return(0);
}


int HDF5::ReadH5Char( string GroupName, string Name, string & out){
	
	Exception::dontPrint();
	
	// Open Group
	Group group = open_group(file,GroupName);
		
	// Now open attribute
	Attribute att;
	try{
		att = group.openAttribute( Name );
	}catch(FileIException error){
		cout << "HDF5 Read :  " << GroupName << " attribute " << Name << " does not exist" << endl;
	}
	
	// Get the data type
	H5std_string strreadbuf ("");
	StrType str_type(PredType::C_S1,512);
	att.read( str_type, strreadbuf ); 
	
	out = strreadbuf;
	
	return(0);
}

int HDF5::AddH5Char( string GroupName, string Name, string S){
	
	Exception::dontPrint();
	
	// Create Group
	Group group = create_group(file,GroupName);
		
	StrType str_type(PredType::C_S1,512);
	DataSpace att_space(H5S_SCALAR);
	
	Attribute att = group.createAttribute( Name, str_type, att_space );
	att.write( str_type, S );
		
	return(0);
}

// Template to write a float
template< int N>
void H5ArrayWrite(H5::H5File &file, string GroupName, string Name, blitz::Array< float,N>& A){
	
	Exception::dontPrint();
	H5::Group group = HDF5::create_group(file,GroupName);
	
	// Create DataSet
	hsize_t dimsf[N];              // dataset dimensions
	for( int i=0; i < N; i++){
		dimsf[i] = A.length(N-i-1);
	}
	DataSpace dataspace( N, dimsf );

	DataSet dataset( group.createDataSet(Name, PredType::NATIVE_FLOAT,dataspace));
	dataset.write( A.data(),PredType::NATIVE_FLOAT, dataspace);
}

int HDF5::AddH5Array( string GroupName,string Name, Array<float,1> & A){ 	H5ArrayWrite<1>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<float,2> & A){ 	H5ArrayWrite<2>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<float,3> & A){ 	H5ArrayWrite<3>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<float,4> & A){ 	H5ArrayWrite<4>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<float,5> & A){ 	H5ArrayWrite<5>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<float,6> & A){ 	H5ArrayWrite<6>(file,GroupName,Name,A);	return(0); }

// Template to write a double
template< int N>
void H5ArrayWrite(H5::H5File &file, string GroupName, string Name, blitz::Array< double,N>& A){
	
	Exception::dontPrint();
	Group group = HDF5::create_group(file,GroupName);
	
	// Create DataSet
	hsize_t dimsf[N];              // dataset dimensions
	for( int i=0; i < N; i++){
		dimsf[i] = A.length(N-i-1);
	}
	DataSpace dataspace( N, dimsf );

	DataSet dataset( group.createDataSet(Name, PredType::NATIVE_DOUBLE,dataspace));
	dataset.write( A.data(),PredType::NATIVE_DOUBLE, dataspace);
}
int HDF5::AddH5Array( string GroupName,string Name, Array<double,1> & A){ 	H5ArrayWrite<1>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<double,2> & A){ 	H5ArrayWrite<2>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<double,3> & A){ 	H5ArrayWrite<3>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<double,4> & A){ 	H5ArrayWrite<4>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<double,5> & A){ 	H5ArrayWrite<5>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<double,6> & A){ 	H5ArrayWrite<6>(file,GroupName,Name,A);	return(0); }

// Template to write a complex<float>
template< int N>
void H5ArrayWrite(H5::H5File &file, string GroupName, string Name, blitz::Array< complex<float>,N>& A){
	
	Exception::dontPrint();
	Group group = HDF5::create_group(file,GroupName);
	
	// Create DataSet
	hsize_t dimsf[N];              // dataset dimensions
	for( int i=0; i < N; i++){
		dimsf[i] = A.length(N-i-1);
	}
	DataSpace dataspace( N, dimsf );
	
	/*DataType Needed for Complex<float>*/
	CompType datatype(sizeof(complex<float>));
	datatype.insertMember( "real", 0, PredType::NATIVE_FLOAT);
	datatype.insertMember( "imag", sizeof(float), PredType::NATIVE_FLOAT);
		
	DataSet dataset( group.createDataSet(Name,datatype,dataspace));
	dataset.write( A.data(),datatype, dataspace);
}
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<float>,1> & A){ 	H5ArrayWrite<1>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<float>,2> & A){ 	H5ArrayWrite<2>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<float>,3> & A){ 	H5ArrayWrite<3>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<float>,4> & A){ 	H5ArrayWrite<4>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<float>,5> & A){ 	H5ArrayWrite<5>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<float>,6> & A){ 	H5ArrayWrite<6>(file,GroupName,Name,A);	return(0); }

// Template to write a complex<double>
template< int N>
void H5ArrayWrite(H5::H5File &file, string GroupName, string Name, blitz::Array< complex<double>,N>& A){
	
	Exception::dontPrint();
	Group group = HDF5::create_group(file,GroupName);
	
	// Create DataSet
	hsize_t dimsf[N];              // dataset dimensions
	for( int i=0; i < N; i++){
		dimsf[i] = A.length(N-i-1);
	}
	DataSpace dataspace( N, dimsf );
	
	/*DataType Needed for Complex<float>*/
	CompType datatype(sizeof(complex<float>));
	datatype.insertMember( "real", 0, PredType::NATIVE_DOUBLE);
	datatype.insertMember( "imag", sizeof(float), PredType::NATIVE_DOUBLE);
		
	DataSet dataset( group.createDataSet(Name,datatype,dataspace));
	dataset.write( A.data(),datatype, dataspace);
}
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<double>,1> & A){ 	H5ArrayWrite<1>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<double>,2> & A){ 	H5ArrayWrite<2>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<double>,3> & A){ 	H5ArrayWrite<3>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<double>,4> & A){ 	H5ArrayWrite<4>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<double>,5> & A){ 	H5ArrayWrite<5>(file,GroupName,Name,A);	return(0); }
int HDF5::AddH5Array( string GroupName,string Name, Array<complex<double>,6> & A){ 	H5ArrayWrite<6>(file,GroupName,Name,A);	return(0); }

// Template to read a float
template< typename T, int N>
void H5ArrayRead(H5::H5File &file, string GroupName, string Name, blitz::Array<T,N>& A){
	
	// Open the file
	H5::Group group = HDF5::open_group(file,GroupName);
	
	// Get the dataset
	DataSet dataset = group.openDataSet( Name );
	
	// Get the dataspace
	DataSpace dataspace = dataset.getSpace();
	
	// Get the number of dimensions
	int rank = dataspace.getSimpleExtentNdims();
	if( rank != N ){
		cout << "Expected Rank " << N << " but opened a rank " << rank << " file " << endl << flush;
	}
			
	// Create DataSet
	hsize_t dimsf[6];             
	int ndims = dataspace.getSimpleExtentDims( dimsf, NULL);
	
	unsigned int total = dataspace.getSelectNpoints();
	if( total != A.numElements()){
		cout << "Dataspace has " << total << " points while array has " << A.numElements() << endl;
		cout << "Resizing input array" << endl;
		TinyVector<int,N> dims;
		for( int i=0; i < N; i++){
			dims(i) = dimsf[ndims-i-1];
		}
		A.setStorage(ColumnMajorArray<N>());
		A.resize( dims );
	}
	
	DataType data_type = DataType(dataset.getDataType());	
		
	// Read the data
	dataset.read( A.data(),data_type);
}
int HDF5::ReadH5Array( string GroupName,string Name, Array<float,1> & A){ 	H5ArrayRead<float,1>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<float,2> & A){ 	H5ArrayRead<float,2>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<float,3> & A){ 	H5ArrayRead<float,3>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<float,4> & A){ 	H5ArrayRead<float,4>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<float,5> & A){ 	H5ArrayRead<float,5>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<float,6> & A){ 	H5ArrayRead<float,6>(file,GroupName,Name,A);	return(0); }

int HDF5::ReadH5Array( string GroupName,string Name, Array<double,1> & A){ 	H5ArrayRead<double,1>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<double,2> & A){ 	H5ArrayRead<double,2>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<double,3> & A){ 	H5ArrayRead<double,3>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<double,4> & A){ 	H5ArrayRead<double,4>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<double,5> & A){ 	H5ArrayRead<double,5>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<double,6> & A){ 	H5ArrayRead<double,6>(file,GroupName,Name,A);	return(0); }

int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<float>,1> & A){ 	H5ArrayRead<complex<float>,1>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<float>,2> & A){ 	H5ArrayRead<complex<float>,2>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<float>,3> & A){ 	H5ArrayRead<complex<float>,3>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<float>,4> & A){ 	H5ArrayRead<complex<float>,4>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<float>,5> & A){ 	H5ArrayRead<complex<float>,5>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<float>,6> & A){ 	H5ArrayRead<complex<float>,6>(file,GroupName,Name,A);	return(0); }

int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<double>,1> & A){ 	H5ArrayRead<complex<double>,1>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<double>,2> & A){ 	H5ArrayRead<complex<double>,2>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<double>,3> & A){ 	H5ArrayRead<complex<double>,3>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<double>,4> & A){ 	H5ArrayRead<complex<double>,4>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<double>,5> & A){ 	H5ArrayRead<complex<double>,5>(file,GroupName,Name,A);	return(0); }
int HDF5::ReadH5Array( string GroupName,string Name, Array<complex<double>,6> & A){ 	H5ArrayRead<complex<double>,6>(file,GroupName,Name,A);	return(0); }
