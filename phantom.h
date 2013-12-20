#pragma once

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
#include <armadillo>
#include "tictoc.hpp"

#ifndef PI
#define PI 3.14159265359
#endif


/**
 *  @brief
 *  Class to generate fractal vascular tree phantom
 */
class FRACTAL3D{
	public:
		enum PerfType {
			ASL,/*!< Use external perfusion file. */
			ELLIPSE/*!< Use simple ellipse model of perfusion 2. */
		};
		class TFRACT_RAND;
		int fractal_pts;//!<Number of endpoints for fractal phantom
		float alpha;//!<Alpha parameter for branch rule
		float beta;//!<Beta parameter for branch rule
		int Nx;	//!<Size of Phantom in X
		int Ny; //!<Size of Phantom in Y
		int Nz; //!<Size of Phantom in Z
		int Nt; //!<Size of Phantom in time

		// Need to maintain these functions for public calls
		void read_commandline(int numarg, char **pstring);
		void calc_image(int,int);
		void write_matlab_truth_script( const char *);
		void calc_image(NDarray::Array<complex<float>,3>&,int,int);
		void build_tree(int Nx, int Ny, int Nz,int Nt);
				
	private:
		NDarray::Array<int,3> synthetic_perfusion(int xs, int ys, int zs, PerfType ptype);
		void update_children(arma::field<TFRACT_RAND>&tree, int pos);
		arma::field<TFRACT_RAND> create_tree( arma::fmat seeds_start, arma::fmat seeds_stop,arma::fmat X);

		bool debug;
		NDarray::Array< float,4>FUZZY;//!<Fuzzy model densities (3D+compartments)
		NDarray::Array< float,4>FUZZYT;//!<Fuzzy model arrival times (3D+compartments)

};


/**
 *  @brief
 *  Class to generate simulated data from known k-space
 *
 *  Usage Example:
 *  @code
 *  	PHANTOM phantom;
 *		phantom.read_commandline(argc,argv);
 *		phantom.init(256,256,256);
 *  	for(int coil =0; coil < data.Num_Coils; coil++){
 *			phantom.update_smap_biotsavart(coil,data.Num_Coils);
 *
 *			for(int t =0; t < recon.rcframes; t++){
 *				// Get Image
 *				phantom.calc_image(t,recon.rcframes);
 *
 *				your_function( phantom.IMAGE );
 *				}
 *			}
 *
 *			// Add Noise
 *			phantom.add_noise( data.kdata );
 *			phantom.write_arma::matlab_truth_script("PhantomData/");
 *	@endcode
 */
class PHANTOM{
	public:
		
		enum PhantomType {
			FRACTAL, /*!< Fractal Vasculature Phantom. */
			SHEPP, /*!< Simple Shepp-Logg Phantom */
			PSF, /*!< Set kdata to ones. */
			EXTERNAL/*!< Use external image domain phantom. */
		};

		// Size of Phantom
		int Nx;	//!<Size of Phantom in X
		int Ny; //!<Size of Phantom in Y
		int Nz; //!<Size of Phantom in Z
		int Nt; //!<Size of Phantom in time
		float over_res; //!<Over-resolution factor for discrete phantoms
		
		// Structures for holding types of phantoms
		FRACTAL3D fractal;

		// Noise Factor
		float phantom_noise; //!<Amount of noise to add in percent
		
		// Type of Phantom
		PhantomType phantom_type; //!<Defines phantom type
		
		// Fractal Input
		
		// External Phantom
		char *external_phantom_name;//!<Name of external phantom
		
		// Arrays for holding Images
		NDarray::Array< complex<float>,3>IMAGE;//!<Final discrete image
		NDarray::Array< complex<float>,3>SMAP;//!<Sensitivity map

		// Constructor,I/O, and Init
		PHANTOM(void);
		void init(int,int,int,int);
		static void help_message(void);
		void read_commandline(int numarg, char **pstring);
		bool phantom_export_smaps;
		
		// Functions
		void add_phase(void);
		void add_noise( NDarray::Array< NDarray::Array<complex<float>,3>,2>&kdata);
		void update_smap_biotsavart(int,int);
		void calc_image(int,int);
		void write_matlab_truth_script( const char *);
		

	private:	
		bool debug;
};










