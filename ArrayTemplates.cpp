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
#include <fftw3.h>

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

inline void fftshift( Array< complex<float>,3>& temp){
	for(int k=0; k<temp.extent(thirdDim);k++){
	for(int j=0; j<temp.extent(secondDim);j++){
	for(int i=0; i<temp.extent(firstDim);i++){
		int mod = ((i+j+k)%2) == 0 ? 1 : -1;
        temp(i,j,k) *= mod;
	}}}	
}


inline void ifft( Array< complex<float>,3>& temp){
	fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());
     
	fftwf_complex *ptr = reinterpret_cast<fftwf_complex*>(temp.data());
		
	//cout << " Planning FFT " << endl << flush; 
	fftwf_plan fft_plan = fftwf_plan_dft_3d(temp.length(thirdDim),temp.length(secondDim),temp.length(firstDim),ptr,ptr,FFTW_BACKWARD, FFTW_MEASURE);
	
	fftshift( temp);
	fftwf_execute(fft_plan);	
	fftshift( temp);
}


inline void fft( Array< complex<float>,3>& temp){
	fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());
     
	fftwf_complex *ptr = reinterpret_cast<fftwf_complex*>(temp.data());
		
	//cout << " Planning FFT " << endl << flush; 
	fftwf_plan fft_plan = fftwf_plan_dft_3d(temp.length(thirdDim),temp.length(secondDim),temp.length(firstDim),ptr,ptr,FFTW_FORWARD, FFTW_MEASURE);
	
	fftshift( temp);
	fftwf_execute(fft_plan);	
	fftshift( temp);
}


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


template< typename T>
void ArrayWritePhase( Array<complex<T>,2>& temp, char *name){
	 
	 T *buffer = new T[temp.length(0)]; 
	 ofstream ofs(name, ios_base::binary);
	 	for(int j= 0; j<temp.extent(1); j++){
			
			for(int i=0; i<temp.extent(0);i++){
				buffer[i]= arg(temp(i,j));
			}
			ofs.write( (char *)buffer,temp.length(0)*sizeof(T));
     }	
 delete [] buffer;
}
template< typename T>
void ArrayWritePhase( Array<complex<T>,3>& temp, char *name){
	 
	 T *buffer = new T[temp.length(0)]; 
	 ofstream ofs(name, ios_base::binary);
	 for(int k=0; k< temp.extent(2);k++){
	 	for(int j= 0; j<temp.extent(1); j++){
			
			for(int i=0; i<temp.extent(0);i++){
				buffer[i]= arg(temp(i,j,k));
			}
			ofs.write( (char *)buffer,temp.length(0)*sizeof(T));
     }}
	 delete [] buffer;
}

template< typename T>
void ArrayWritePhase( Array< complex<T>,4>& temp, char *name){
	 
	 T *buffer = new T[temp.length(0)]; 
	 ofstream ofs(name, ios_base::binary);
	 for(int t=0; t< temp.extent(3);t++){
	 for(int k=0; k< temp.extent(2);k++){
	 	for(int j= 0; j<temp.extent(1); j++){
			
			for(int i=0; i<temp.extent(0);i++){
				buffer[i]= arg(temp(i,j,k,t));
			}
			ofs.write( (char *)buffer,temp.length(0)*sizeof(T));
     }}}
	 delete [] buffer;
}


template< typename T>
void ArrayWritePhase( Array< complex<T>,5>& temp, char *name){
	 
	 T *buffer = new T[temp.length(0)]; 
	 ofstream ofs(name, ios_base::binary);
	 for(int e=0; e< temp.extent(4);e++){
	 for(int t=0; t< temp.extent(3);t++){
	 for(int k=0; k< temp.extent(2);k++){
	 	for(int j= 0; j<temp.extent(1); j++){
			
			for(int i=0; i<temp.extent(0);i++){
				buffer[i]= arg(temp(i,j,k,t,e));
			}
			ofs.write( (char *)buffer,temp.length(0)*sizeof(T));
     }}}}
	 delete [] buffer;
}
