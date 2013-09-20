#pragma once

// Switching to Blitz Based Arrays
#include <blitz/array.h>
#include <fftw3.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string.h>
#include <complex>
#include <omp.h>

// It should be ok at populate namespace with complex
using std::complex;

namespace NDarray {

using namespace blitz;

// FFT Libraries complex Float
void fftshift( Array<complex<float>,3>& temp);
void fft( Array<complex<float>,3>& temp);
void ifft( Array<complex<float>,3>& temp);

void fft( Array<complex<float>,3>& temp,int);
void ifft( Array<complex<float>,3>& temp,int);
void fft3( Array<complex<float>,3>& temp,int,int);

inline void endian_swap( int& x){
    x = ( x<<24 & 0xFF000000) |
        ( x<<8  & 0x00FF0000) |
        ( x>>8  & 0x0000FF00) |
        ( x>>24 & 0x000000FF);
}

template< typename T, int N>
void ArrayRead( blitz::Array< T,N>& temp, const char *name){
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

template< typename T, const int N_rank>
void ArrayWriteMagAppend( Array<complex<T>,N_rank>& temp, const char *name){
	 
	 ofstream ofs(name, ios_base::binary | ios_base::app);
	 for(typename Array<complex<T>,N_rank>::iterator miter=temp.begin(); miter!=temp.end(); miter++){
		T val = abs( *miter);
		ofs.write( (char *)&val,sizeof(T));
     }	
}

template< typename T, const int N_rank>
void ArrayWrite( Array< T , N_rank>& temp, const char *name){
	 ofstream ofs(name, ios_base::binary);
	 for(typename Array<T,N_rank>::iterator miter=temp.begin(); miter!=temp.end(); miter++){
			T val= *miter;
			ofs.write( (char *)&val,sizeof(T));
     }
}

template< typename T, const int N_rank>
void ArrayWriteMag( Array<complex<T>, N_rank>& temp, const char *name){
	 ofstream ofs(name, ios_base::binary);
	 for(typename Array<complex<T>,N_rank>::iterator miter=temp.begin(); miter!=temp.end(); miter++){
			T val=abs( *miter);
			ofs.write( (char *)&val,sizeof(T));
     }
}

template< typename T, const int N_rank >
T ArrayEnergy( Array< complex< T >, N_rank >& temp){
	T EE=0;
	for(typename Array<complex<T>,N_rank>::iterator miter=temp.begin(); miter!=temp.end(); miter++){
		EE+= norm( *miter );
	}
    return(EE);
}

template< typename T, const int N_rank, const int M_rank >
T ArrayEnergy( Array< Array< complex< T >,N_rank>, M_rank >& temp){
	T EE=0;
	for(typename Array< Array< complex<T>,N_rank>,M_rank >::iterator miter=temp.begin(); miter!=temp.end(); miter++){
		EE+= ArrayEnergy( *miter );
	}
    return(EE);
}


template< typename T, const int N_rank>
void ArrayWritePhase( Array<complex<T>,N_rank>& temp, const char *name){
	 
	 ofstream ofs(name, ios_base::binary);
	 for(typename Array<complex<T>,N_rank>::iterator miter=temp.begin(); miter!=temp.end(); miter++){
	 	T val= arg( *miter);
		ofs.write( (char *)&val,sizeof(T));
     }	
}

template< typename T, const int N_rank>
void ArrayWritePhaseAppend( Array<complex<T>,N_rank>& temp, const char *name){
	 
	 ofstream ofs(name, ios_base::binary  | ios_base::app);
	 for(typename Array<complex<T>,N_rank>::iterator miter=temp.begin(); miter!=temp.end(); miter++){
	 	T val= arg( *miter);
		ofs.write( (char *)&val,sizeof(T));
     }	
}

template < typename T > 
Array< Array<T,3>, 1> Alloc4DContainer( int x, int y, int z, int t){
	Array< Array<T,3>,1> temp;
	temp.setStorage(ColumnMajorArray<1>());
	temp.resize( t );
	for(int i1=0; i1< t; i1++){
		temp(i1).setStorage(ColumnMajorArray<3>());
		temp(i1).resize(x,y,z);
		temp(i1)= (T )0;
	}
	return(temp);
}

template < typename T > 
Array< Array<T,3>, 3> Alloc6DContainer( int x, int y, int z, int d1, int d2, int d3){
	Array< Array<T,3>,3> temp;
	temp.setStorage(ColumnMajorArray<3>());
	temp.resize( d1,d2,d3);
	
	for( typename Array<Array<T,3>,3>::iterator miter=temp.begin();   miter !=temp.end(); miter++){
		(*miter).setStorage(ColumnMajorArray<3>());
		(*miter).resize(x,y,z);
		(*miter)= (T )0;
	}
	return(temp);
}

template < typename T > 
Array< Array<T,3>, 2> Alloc5DContainer( int x, int y, int z, int d1, int d2){
	Array< Array<T,3>,2> temp;
	temp.setStorage(ColumnMajorArray<2>());
	temp.resize( d1,d2);
	
	for( typename Array< Array<T,3>,2>::iterator miter=temp.begin();   miter !=temp.end(); miter++){
		(*miter).setStorage(ColumnMajorArray<3>());
		(*miter).resize(x,y,z);
		(*miter)= (T )0;
	}
	return(temp);
}


}// Namespace

