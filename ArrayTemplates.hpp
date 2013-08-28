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

using namespace std;
using namespace blitz;

// FFT Libraries complex Float
void fftshift( Array< complex<float>,3>& temp);
void fft( Array< complex<float>,3>& temp);
void ifft( Array< complex<float>,3>& temp);
void fft( Array< complex<float>,3>& temp,int);
void ifft( Array< complex<float>,3>& temp,int);
void fft3( Array< complex<float>,3>& temp,int,int);

inline void endian_swap( int& x){
    x = ( x<<24 & 0xFF000000) |
        ( x<<8  & 0x00FF0000) |
        ( x>>8  & 0x0000FF00) |
        ( x>>24 & 0x000000FF);
}

template< typename T, int N>
void ArrayRead( Array< T,N>& temp, const char *name){
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
	 
	 fstream filestr;
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


template< typename T, const int N_rank>
void ArrayWritePhase( Array<complex<T>,N_rank>& temp, const char *name){
	 
	 ofstream ofs(name, ios_base::binary);
	 for(typename Array<complex<T>,N_rank>::iterator miter=temp.begin(); miter!=temp.end(); miter++){
	 	T val= arg( *miter);
		ofs.write( (char *)&val,sizeof(T));
     }	
}
