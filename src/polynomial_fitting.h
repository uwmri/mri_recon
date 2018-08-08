#ifndef hPOLYFIT
#define hPOLYFIT

/**Code to Fit image data to a polynomial*/
 
#include <armadillo>
#include <omp.h>

// For Wrapper (this is the actual recon + NDarray)
#include "ArrayTemplates.hpp"
 
class POLYFIT {
	public:

		/*Code for fitting image to polynomial**/
		double poly3d(float x, float y, float z);
		double poly2d(float x, float y);
				
		/* Fitting */
		void poly_fitting3d(NDarray::Array< float,3> &back_mag /*Binary Matrix*/,blitz::Array<float,3> &image, int);
		
		/* Subtraction */
		void poly_subtract3d(NDarray::Array<float,3> &);
				
	private:
		int number;
		arma::vec alpha;
		arma::vec px;
		arma::vec py;
		arma::vec pz;
};		
#endif
