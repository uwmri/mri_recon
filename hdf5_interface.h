#include <H5Cpp.h>
#include "ArrayTemplates.hpp"

class HDF5{
	
	public:		
	// HDF5 File	
	H5::H5File file;
	
	// Constructor
	HDF5(void); 
	HDF5(const char *name); 
	
	// Allow Adding of Arrays
	int AddH5Array( const char *GroupName, const char *Name, NDarray::Array<float,3> & A);
	int AddH5Array( const char *GroupName, const char *Name, NDarray::Array<float,4> & A);
	int AddH5Array( const char *GroupName, const char *Name, NDarray::Array<float,5> & A);
	int AddH5Array( const char *GroupName, const char *Name, NDarray::Array<float,6> & A);
		
	int AddH5Array( const char *GroupName, const char *Name, NDarray::Array<complex<float>,3> & A);
	int AddH5Array( const char *GroupName, const char *Name, NDarray::Array<complex<float>,4> & A);
	int AddH5Array( const char *GroupName, const char *Name, NDarray::Array<complex<float>,5> & A);
	int AddH5Array( const char *GroupName, const char *Name, NDarray::Array<complex<float>,6> & A);
	
	int AddH5Scaler( const char *GroupName, const char *Name, int A);
	int AddH5Scaler( const char *GroupName, const char *Name, float A);
	int AddH5Char( const char *GroupName, const char *Name, char *A);

	private:
	
};

