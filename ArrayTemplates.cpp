#pragma once

// Array Class for MRI Reconstruction
//   Used to allow very specific control over memory
//   and parallelization of key operations 
// Main Classes: 
//		Array2D,Array3D	 - With contigous memory access
//		Array4D,Array5D	 - Non-contigous memory access (3D are contigous)



// Switching to Blitz Based Arrays
#include "blitz/array.h"
using namespace blitz;

/*
 * Class RowMajorArray specializes GeneralArrayStorage to provide column
 * major arrays (column major ordering, base of 0).
 */

template<int N_rank>
class RowMajorArray : public GeneralArrayStorage<N_rank> {
private:
    typedef GeneralArrayStorage<N_rank> T_base;
    typedef _bz_typename T_base::noInitializeFlag noInitializeFlag;
    using T_base::ordering_;
    using T_base::ascendingFlag_;
    using T_base::base_;
public:
    RowMajorArray()
        : GeneralArrayStorage<N_rank>(noInitializeFlag())
    {
	for (int i=0; i < N_rank; ++i)
          ordering_(i) = i;        
	ascendingFlag_ = true;
        base_ = 0;
    }
};


template< typename T, int N>
void ArrayRead( Array< T,N>& temp, char *name){
	 	FILE *fid;
		if( (fid=fopen(name,"r")) == NULL){
			cout << "Array:Can't Open " << name << endl;
			cout << "Exiting" << endl;
			exit(1);
		}else{	
			int j;
			if( (j=fread(temp.data(),sizeof(T),temp.numElements(),fid)) != (int)(temp.numElements())){
				cout << "Array3:Not enough data: only read " << j << "points of" <<  temp.numElements() << endl;
				exit(1);
			}
		}
}



template< typename T, int N>
void ArrayWrite( Array< T,N>& temp, char *name){
	 	FILE *fid;
		if( (fid=fopen(name,"w")) == NULL){
			cout << "ArrayWrite:Can't Open " << name << endl;
			cout << "Exiting" << endl;
			exit(1);
		}else{	
			fwrite(temp.data(),temp.numElements(),sizeof(T),fid);
			fclose(fid);
		}
}
 
template< typename T>
void ArrayWriteMag( Array<complex<T>,2>& temp, char *name){
	 
	 T *buffer = new T[temp.length(0)]; 
	 ofstream ofs(name, ios_base::binary);
	 	for(int j= 0; j<temp.extent(1); j++){
			
			for(int i=0; i<temp.extent(0);i++){
				buffer[i]= abs(temp(i,j));
			}
			ofs.write( (char *)buffer,temp.length(0)*sizeof(T));
     }	
 delete [] buffer;
}
template< typename T>
void ArrayWriteMag( Array<complex<T>,3>& temp, char *name){
	 
	 T *buffer = new T[temp.length(0)]; 
	 ofstream ofs(name, ios_base::binary);
	 for(int k=0; k< temp.extent(2);k++){
	 	for(int j= 0; j<temp.extent(1); j++){
			
			for(int i=0; i<temp.extent(0);i++){
				buffer[i]= abs(temp(i,j,k));
			}
			ofs.write( (char *)buffer,temp.length(0)*sizeof(T));
     }}
	 delete [] buffer;
}

template< typename T>
T ArrayEnergy( Array<complex<T>,5>& temp){
T EE=0;
for(int e=0; e< temp.extent(4);e++){	 
for(int t=0; t< temp.extent(3);t++){
	 for(int k=0; k< temp.extent(2);k++){
	 	for(int j= 0; j<temp.extent(1); j++){
			for(int i=0; i<temp.extent(0);i++){
				EE += norm(temp(i,j,k,t,e));
					
     }}}}}
	 return(EE);
}

template< typename T>
void ArrayWriteMag( Array< complex<T>,4>& temp, char *name){
	 
	 T *buffer = new T[temp.length(0)]; 
	 ofstream ofs(name, ios_base::binary);
	 for(int t=0; t< temp.extent(3);t++){
	 for(int k=0; k< temp.extent(2);k++){
	 	for(int j= 0; j<temp.extent(1); j++){
			
			for(int i=0; i<temp.extent(0);i++){
				buffer[i]= abs(temp(i,j,k,t));
			}
			ofs.write( (char *)buffer,temp.length(0)*sizeof(T));
     }}}
	 delete [] buffer;
}


template< typename T>
void ArrayWriteMag( Array< complex<T>,5>& temp, char *name){
	 
	 T *buffer = new T[temp.length(0)]; 
	 ofstream ofs(name, ios_base::binary);
	 for(int e=0; e< temp.extent(4);e++){
	 for(int t=0; t< temp.extent(3);t++){
	 for(int k=0; k< temp.extent(2);k++){
	 	for(int j= 0; j<temp.extent(1); j++){
			
			for(int i=0; i<temp.extent(0);i++){
				buffer[i]= abs(temp(i,j,k,t,e));
			}
			ofs.write( (char *)buffer,temp.length(0)*sizeof(T));
     }}}}
	 delete [] buffer;
}


