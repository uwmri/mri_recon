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
enum {SP_LOW_RES_PHASE, SP_COIL_PHASE, SP_SMOOTH,SP_SMOOTH_OLD};


void fftshift3(Array<complex<float>,3>&);
void ifft3(Array< complex<float>,3>&);

class SPIRIT {
  public: 
    Array< complex<float>,5> k;
    Array< complex<float>,5> im;
    
    // enum flags
    int shape;
    int calib_type;
    
    float calib_lam;
    
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
    
    void generateEigenCoils(Array< complex<float>,4> &);
    void calibrate_ellipsoid(Array< complex<float>,4>&);
    void prep();
    
    void getcoils( Array< complex<float>,4>&);
    void phase_correct(Array< complex<float>,4>&maps, Array< complex<float>,4>&ref);
    void interpMaps(Array< complex<float>,4>&, Array< complex<float>,4>&);

};

