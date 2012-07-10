#ifndef hWAVE3DLIB
#define hWAVE3DLIB

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <complex>
using namespace std;
#include <omp.h>
#include "ArrayTemplates.cpp"

#define WAVE_DB2 0
#define WAVE_DB4 1
#define WAVE_DB6 2
#define WAVE_DB8 3
#define WAVE_SYM2 4
#define WAVE_SYM4 5
#define WAVE_SYM6 6
#define WAVE_SYM8 7
#define WAVE_BO33 8

#define THREADS 32

class WAVELET3D{
	public:
		int N[5];
		int L[5];
		int W[5];
		int S[5];
		
		int threads;
		int max_level;
		
		float *lpf[5];
		float *hpf[5];
		float *Slpf[5];
		float *Shpf[5];
		int wN[5];
		
		int transform_in_z;
				
		array5D< complex<float> >Coef;
						
		//WAVELET3D(int,int,int,int,int);
		//WAVELET3D( array3D< complex<float> > *,int,int);
		
		WAVELET3D( array3D< complex<float> > *,int,int);
		WAVELET3D( array3D< complex<float> > *,int *,int *);
		
		WAVELET3D( array4D< complex<float> > *,int,int);
		WAVELET3D( array4D< complex<float> > *,int *,int *);
		
		WAVELET3D( array5D< complex<float> > *,int ,int );
		WAVELET3D( array5D< complex<float> > *,int *,int *);
		~WAVELET3D();
		
		void get_filter_banks();
		void forward(void);		
		void backward(void);
		void random_shift(void);		
		void wave_x(int,int,int,int,int,int);
		void wave_y(int,int,int,int,int,int);
		void wave_z(int,int,int,int,int,int);
		void wave_e(int,int,int,int,int,int);
		void wave_t(int,int,int,int,int,int);
		void wave1D( complex<float> [],complex<float> [],int,int);
		void iwave1D( complex<float> [],complex<float> [],int,int);
		void wave_threshold(float);
		void setup_wavelet_directions( void);
	private:	
		
};


#endif

