#include "hdf5_interface.h"

using namespace std; 
using namespace NDarray;
using namespace H5;

/*
herr_t file_info( hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata){
	hid_t group;
	group = H5Gopen2(loc_id,name,H5P_DEFAULT);
	
	cout << "Name: " << name <<endl;
	H5Gclose(group);
	return(0);
} */

HDF5::HDF5( const char *FileName ){
	
	// OPen the file
	file =  H5File(FileName,H5F_ACC_TRUNC);	

}

HDF5::HDF5( void  ){


}

int HDF5::AddH5Scaler( const char *GroupName, const char *Name, int A){
	
	Exception::dontPrint();
	
	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	
	IntType int_type(PredType::STD_I32LE);
	DataSpace att_space(H5S_SCALAR);
	Attribute att = group.createAttribute( Name , int_type, att_space );
	att.write( int_type, &A );
		
	return(0);
}

int HDF5::AddH5Scaler( const char *GroupName, const char *Name, float A){
	
	Exception::dontPrint();
	
	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	
	FloatType float_type(PredType::NATIVE_FLOAT);
	DataSpace att_space(H5S_SCALAR);
	Attribute att = group.createAttribute( Name , float_type, att_space );
	att.write( float_type, &A );
		
	return(0);
}

int HDF5::AddH5Char( const char *GroupName, const char *Name, char *S){
	
	Exception::dontPrint();
	
	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	StrType str_type(PredType::NATIVE_CHAR);
	DataSpace att_space(H5S_SCALAR);
	Attribute att = group.createAttribute( Name, str_type, att_space );
	att.write( str_type, S );
		
	return(0);
}


int HDF5::AddH5Array( const char *GroupName, const char *Name, Array<float,2> & A){
	
	Exception::dontPrint();

	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	// Create DataSet
	hsize_t dimsf[2];              // dataset dimensions
	dimsf[0] = A.length(secondDim);
	dimsf[1] = A.length(firstDim);
	DataSpace dataspace( 2, dimsf );
	
	/* Write Data*/
	DataSet dataset( group.createDataSet(Name, PredType::NATIVE_FLOAT,dataspace));
	dataset.write( A.data(),PredType::NATIVE_FLOAT, dataspace);

	return(0);
}

	
int HDF5::AddH5Array( const char *GroupName, const char *Name, Array<float,3> & A){
	
	Exception::dontPrint();

	
	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	// Create DataSet
	hsize_t dimsf[3];              // dataset dimensions
	dimsf[0] = A.length(thirdDim);
	dimsf[1] = A.length(secondDim);
	dimsf[2] = A.length(firstDim);
	DataSpace dataspace( 3, dimsf );
	
	/* Write Data*/
	DataSet dataset( group.createDataSet(Name, PredType::NATIVE_FLOAT,dataspace));
	dataset.write( A.data(),PredType::NATIVE_FLOAT, dataspace);

	return(0);
}

int HDF5::AddH5Array( const char *GroupName, const char *Name, Array<float,4> & A){
	
	Exception::dontPrint();

	
	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	// Create DataSet
	hsize_t dimsf[4];              // dataset dimensions
	dimsf[0] = A.length(fourthDim);
	dimsf[1] = A.length(thirdDim);
	dimsf[2] = A.length(secondDim);
	dimsf[3] = A.length(firstDim);
	DataSpace dataspace( 4, dimsf );
	
	/* Write Data*/
	DataSet dataset( group.createDataSet(Name, PredType::NATIVE_FLOAT,dataspace));
	dataset.write( A.data(),PredType::NATIVE_FLOAT, dataspace);

	return(0);
}


int HDF5::AddH5Array( const char *GroupName, const char *Name, Array<float,5> & A){
	
	Exception::dontPrint();

	
	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	// Create DataSet
	hsize_t dimsf[5];              // dataset dimensions
	dimsf[0] = A.length(fifthDim);
	dimsf[1] = A.length(fourthDim);
	dimsf[2] = A.length(thirdDim);
	dimsf[3] = A.length(secondDim);
	dimsf[4] = A.length(firstDim);
	DataSpace dataspace( 5, dimsf );
	
	/* Write Data*/
	DataSet dataset( group.createDataSet(Name, PredType::NATIVE_FLOAT,dataspace));
	dataset.write( A.data(),PredType::NATIVE_FLOAT, dataspace);

	return(0);
}



int HDF5::AddH5Array( const char *GroupName,const char *Name, Array<complex<float>,2> & A){
	
	Exception::dontPrint();
		
	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	// Create DataSet
	hsize_t dimsf[2];              // dataset dimensions
	dimsf[0] = A.length(secondDim);
	dimsf[1] = A.length(firstDim);
	DataSpace dataspace( 2, dimsf );
	
	
	/*DataType Needed for Complex<float>*/
	CompType datatype(sizeof(complex<float>));
	datatype.insertMember( "real", 0, PredType::NATIVE_FLOAT);
	datatype.insertMember( "imag", sizeof(float), PredType::NATIVE_FLOAT);
	
	/* Write Data*/
	DataSet dataset( group.createDataSet(Name, datatype,dataspace));
	dataset.write( A.data(),datatype, dataspace);
	
	return(0);
	
}



int HDF5::AddH5Array( const char *GroupName,const char *Name, Array<complex<float>,3> & A){
	
	Exception::dontPrint();
		
	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	// Create DataSet
	hsize_t dimsf[3];              // dataset dimensions
	dimsf[0] = A.length(thirdDim);
	dimsf[1] = A.length(secondDim);
	dimsf[2] = A.length(firstDim);
	DataSpace dataspace( 3, dimsf );
	
	
	/*DataType Needed for Complex<float>*/
	CompType datatype(sizeof(complex<float>));
	datatype.insertMember( "real", 0, PredType::NATIVE_FLOAT);
	datatype.insertMember( "imag", sizeof(float), PredType::NATIVE_FLOAT);
	
	/* Write Data*/
	DataSet dataset( group.createDataSet(Name, datatype,dataspace));
	dataset.write( A.data(),datatype, dataspace);
	
	return(0);
	
}

int HDF5::AddH5Array( const char *GroupName,const char *Name, Array<complex<float>,4> & A){
	
	Exception::dontPrint();
		
	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	// Create DataSet
	hsize_t dimsf[4];              // dataset dimensions
	dimsf[0] = A.length(fourthDim);
	dimsf[1] = A.length(thirdDim);
	dimsf[2] = A.length(secondDim);
	dimsf[3] = A.length(firstDim);
	DataSpace dataspace( 4, dimsf );
	
	/*DataType Needed for Complex<float>*/
	CompType datatype(sizeof(complex<float>));
	datatype.insertMember( "real", 0, PredType::NATIVE_FLOAT);
	datatype.insertMember( "imag", sizeof(float), PredType::NATIVE_FLOAT);
	
	/* Write Data*/
	DataSet dataset( group.createDataSet(Name, datatype,dataspace));
	dataset.write( A.data(),datatype, dataspace);
	
	return(0);
	
}

int HDF5::AddH5Array( const char *GroupName,const char *Name, Array<complex<float>,5> & A){
	
	Exception::dontPrint();
		
	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	// Create DataSet
	hsize_t dimsf[5];              // dataset dimensions
	dimsf[0] = A.length(fifthDim);
	dimsf[1] = A.length(fourthDim);
	dimsf[2] = A.length(thirdDim);
	dimsf[3] = A.length(secondDim);
	dimsf[4] = A.length(firstDim);
	DataSpace dataspace( 5, dimsf );
	
	/*DataType Needed for Complex<float>*/
	CompType datatype(sizeof(complex<float>));
	datatype.insertMember( "real", 0, PredType::NATIVE_FLOAT);
	datatype.insertMember( "imag", sizeof(float), PredType::NATIVE_FLOAT);
	
	/* Write Data*/
	DataSet dataset( group.createDataSet(Name, datatype,dataspace));
	dataset.write( A.data(),datatype, dataspace);
	
	return(0);
	
}

int HDF5::AddH5Array( const char *GroupName,const char *Name, Array<complex<float>,6> & A){
	
	Exception::dontPrint();
		
	// Create Group
	Group group;
	try{
		cout << "Trying to Open: " << GroupName << endl;
		group = Group(file.openGroup( GroupName ) ); 
	}catch(FileIException error){
		cout << "Group does not exist-Creating" << endl;
		string x(GroupName);
		x.insert(0,"/");
		cout << x << endl;
		group = Group(file.createGroup( x) );
	}
	
	// Create DataSet
	hsize_t dimsf[6];              // dataset dimensions
	dimsf[0] = A.length(sixthDim);
	dimsf[1] = A.length(fifthDim);
	dimsf[2] = A.length(fourthDim);
	dimsf[3] = A.length(thirdDim);
	dimsf[4] = A.length(secondDim);
	dimsf[5] = A.length(firstDim);
	DataSpace dataspace( 6, dimsf );
	
	/*DataType Needed for Complex<float>*/
	CompType datatype(sizeof(complex<float>));
	datatype.insertMember( "real", 0,             PredType::NATIVE_FLOAT);
	datatype.insertMember( "imag", sizeof(float), PredType::NATIVE_FLOAT);
	
	/* Write Data*/
	cout << "Create Dataset " << endl << flush;
	DataSet dataset( group.createDataSet(Name, datatype,dataspace));
	
	cout << "Write Dataset " << endl << flush;
	dataset.write( A.data(),datatype, dataspace);
	
	return(0);
	
}


