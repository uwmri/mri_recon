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
#include <omp.h>
#include "ArrayTemplates.hpp"

/**
 *  @brief
 *  Class to perform Wavelet Transforms in 3D (spatial dim)
 *
 *  Usage Example:
 *  @code
 *  	include "wavelet3D.h"
 *  	int dirs[3] = {4, 4, 4}; // Maximum 4 levels
 *		WAVELET3D wave3(X,dirs,WAVELET3D::WAVE_DB4);
 *		wave3.forward(); // Transform
 *		wave3.backward(); // UnTransfrom
 *	@endcode
 */
class WAVELET3D{
	public:
		enum WaveType{
			 WAVE_DB2, /*!< Daubauchies 2. */
			 WAVE_DB4, /*!< Daubauchies 4. */
			 WAVE_DB6, /*!< Daubauchies 6. */
			 WAVE_DB8, /*!< Daubauchies 8. */
			 WAVE_SYM2, /*!< Symlet 2. */
			 WAVE_SYM4, /*!< Symlet 4. */
			 WAVE_SYM6, /*!< Symlet 6. */
			 WAVE_SYM8, /*!< Symlet 8. */
			 WAVE_BO33 /*!< Bio-Orthogonal 4. */
		};

		// Attributes Documented Here
		int N[5];	//!<Array Size Place Holder
		int L[3];	//!<Wavelet Levels Place Holder
		int S[3];	//!<Wavelet Shifts Place Holder
		
		int max_level;	//!<Maximum Level Needed.
		
		float *lpf;		//!<Internal Storage for decomposing low pass filter
		float *hpf;		//!<Internal Storage for decomposing high pass filter
		float *Slpf;	//!<Internal Storage for synthesis low pass filter
		float *Shpf;	//!<Internal Storage for synthesis high pass filter
		int wN;			//!<Wavelet Basis Size
		int wType;		//!<Wavelet Basis Type
		
		// Public Member Functions Documented in .cpp
		WAVELET3D();
		WAVELET3D( NDarray::TinyVector<int,3> ,NDarray::TinyVector<int,3> ,int);
		WAVELET3D( NDarray::Array< complex<float> , 3> &,int *,int);
		WAVELET3D( NDarray::Array< complex<float> , 4> &,int *,int);
		WAVELET3D( NDarray::Array< complex<float> , 5> &,int *,int);
		~WAVELET3D();
		
		void forward( NDarray::Array< complex<float> , 3> &);		
		void backward( NDarray::Array< complex<float> , 3> &);
		void forward( NDarray::Array< complex<float> , 4> &);		
		void backward( NDarray::Array< complex<float> , 4> &);
		void forward( NDarray::Array< complex<float> , 5> &);		
		void backward( NDarray::Array< complex<float> , 5> &);

		void random_shift(void);		

	private:
		void forward3D( NDarray::Array< complex<float> , 3> &);
		void backward3D( NDarray::Array< complex<float> , 3> &);

		void setup_wavelet_directions( void);
		void get_filter_banks();
		void setup( NDarray::Array< complex<float> , 3> &,int *,int);

		void wave_x(NDarray::Array< complex<float> , 3> &,int,int,int,int);
		void wave_y(NDarray::Array< complex<float> , 3> &,int,int,int,int);
		void wave_z(NDarray::Array< complex<float> , 3> &,int,int,int,int);

		void wave1D( complex<float> [],complex<float> [],int,int);
		void iwave1D( complex<float> [],complex<float> [],int,int);

};


#endif

