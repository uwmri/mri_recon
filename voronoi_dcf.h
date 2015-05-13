#ifndef hVORDCF
#define hVORDCF

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

#include "voro++/voro++.hh"

/**
 *  @brief
 *  Class to perform Voroni Diagram of MRI Data
 *
 *  Usage Example:
 *  @code
 *  	include "voroni_dcf.h"

 *	@endcode
 */
class VORONOI_DCF{
	public:
		
		static void dcf_3D( NDarray::Array< float,3> &,NDarray::Array< float,3> &,NDarray::Array< float,3> &,NDarray::Array< float,3> &);
		static void dcf_2D( NDarray::Array< float,3> &,NDarray::Array< float,3> &,NDarray::Array< float,3> &);
						
	private:


};


#endif



