#pragma once

#include <iostream>
#include <armadillo>
#include <math.h>
#include <fftw3.h>
#include <omp.h>
#include "tictoc.cpp"
#include "ArrayTemplates.hpp"

enum {SP_SQUARE, SP_CIRCLE};
enum {SP_TIK, SP_TSVD};
enum {SP_COIL_PHASE, SP_SMOOTH};


void fftshift3(Array<complex<float>,3>&);
void ifft3(Array< complex<float>,3>&);

class SPIRIT {
  public: 
    Array< complex<float>,5> k;
    Array< complex<float>,5> im;
    
    // enum flags
    int shape;
        
    int mapshrink;
    float mapthresh;
    
    // Size of the kernel (2*kr+1 for squares)
    int krx;
    int kry;
    int krz;
    
    float kr_f;
    float krx_f;
    float kry_f;
    float krz_f;
    
    int cr;
    int crx;
    int cry;
    int crz;
    
    int rcxres;
    int rcyres;
    int rczres;
    int debug;
	
    int ncoils;
    int nV;
    int phase_type;
	
    static void help_message(void);
	
    SPIRIT();
    void read_commandline(int numarg, char **pstring);
    void init(int, int, int, int);
    
    void generateEigenCoils(Array< Array<complex<float>,3>, 1> &);
    void calibrate_ellipsoid(Array< Array<complex<float>,3>, 1>&);
    void prep();
    
    void getcoils( Array< Array<complex<float>,3>, 1>&);
    void phase_correct(Array< Array<complex<float>,3>, 1>&);
    
};

