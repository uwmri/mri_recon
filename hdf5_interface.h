#include "ArrayTemplates.hpp"

#ifndef HDF_INTERFACE
#define HDF_INTERFACE

#include <H5Cpp.h>
#include <string>

class HDF5{
	
	public:		
	// HDF5 File	
	H5::H5File file;
	
	// Constructor
	HDF5(void); 
	HDF5(std::string name,const std::string); 
	
	// Allow Adding of Arrays
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<float,1> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<float,2> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<float,3> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<float,4> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<float,5> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<float,6> & A);

	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<double,1> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<double,2> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<double,3> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<double,4> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<double,5> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<double,6> & A);

	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,1> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,2> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,3> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,4> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,5> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,6> & A);
	
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,1> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,2> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,3> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,4> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,5> & A);
	int AddH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,6> & A);
	
	// Allow Adding of Arrays
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<float,1> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<float,2> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<float,3> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<float,4> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<float,5> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<float,6> & A);

	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<double,1> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<double,2> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<double,3> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<double,4> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<double,5> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<double,6> & A);

	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,1> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,2> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,3> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,4> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,5> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<float>,6> & A);
	
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,1> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,2> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,3> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,4> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,5> & A);
	int ReadH5Array(std::string GroupName,std::string Name, NDarray::Array<complex<double>,6> & A);
			
	int AddH5Scaler(std::string GroupName,std::string Name, int A);
	int AddH5Scaler(std::string GroupName,std::string Name, float A);
	int AddH5Char(std::string GroupName,std::string Name,std::string A);
	
	int ReadH5Scaler(std::string GroupName,std::string Name, int * A);
	int ReadH5Scaler(std::string GroupName,std::string Name, float * A);
	int ReadH5Char(std::string GroupName,std::string Name, std::string &A);
	
	static H5::Group open_group(  H5::H5File&,std::string GroupName );
	static H5::Group create_group( H5::H5File&,std::string GroupName );
		
	private:
	
};

#endif
