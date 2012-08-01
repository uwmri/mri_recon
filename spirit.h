#ifndef hSPIRIT
#define hSPIRIT

#include <iostream>
#include <armadillo>
#include <math.h>
#include <fftw3.h>
#include <omp.h>
#include "tictoc.cpp"
#include "ArrayTemplates.cpp"

#define PI 3.1415926535898

enum {SP_SQUARE, SP_CIRCLE};
enum {SP_TIK, SP_TSVD};

void fftshift3(int,int,int,complex<float>*);
void ifft3(arrayND< complex<float> > &);

class SPIRIT {
  public: 
    arrayND< complex<float> > k;
    arrayND< complex<float> > im;
    
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
    
    int ncoils;
    
    SPIRIT();
    void read_commandline(int numarg, char **pstring);
    void init(int, int, int, int);
    
    complex<float> getImPointDFT(float,float,float,int,int);
    
    void calibrate_ellipsoid(arrayND< complex<float> > &);
    void prep();
    
    void getcoils(array4D< complex<float> > &);
    void rotateCoils(array4D< complex<float> > &maps, array4D< complex<float> > &ref);
    void interpMaps(array4D< complex<float> >&, array4D< complex<float> >&);
};

#endif
