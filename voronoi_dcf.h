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

 
#include <armadillo>

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
		
		enum KShape{ SPHERE, CYLINDER,CUBE};
		static void vor_dcf( NDarray::Array< float,3> &,NDarray::Array< float,3> &,NDarray::Array< float,3> &,NDarray::Array< float,3> &,KShape);
		
		
		static void vor_sphere( NDarray::Array< float,1> &,NDarray::Array< float,1> &,NDarray::Array<float,1> &,NDarray::Array< float,1> &);				
	private:


};


#endif



