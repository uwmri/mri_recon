#ifndef hWAVE3DLIB
#define hWAVE3DLIB

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
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
		int L[3];
		int S[3];
		
		int threads;
		int max_level;
		
		float *lpf;
		float *hpf;
		float *Slpf;
		float *Shpf;
		int wN;
		int wType;
		
		Array< complex<float>, 3 >Coef;
		
		WAVELET3D( Array< complex<float> , 3> &,int *,int);
		
		~WAVELET3D();
		
		void get_filter_banks();
		void forward( void);		
		void backward( void);
		void random_shift(void);		
		void wave_x(int,int,int,int,int,int);
		void wave_y(int,int,int,int,int,int);
		void wave_z(int,int,int,int,int,int);
		void wave1D( complex<float> [],complex<float> [],int,int);
		void iwave1D( complex<float> [],complex<float> [],int,int);
		void wave_threshold(float);
		void setup_wavelet_directions( void);
	private:	
		
};


#endif

