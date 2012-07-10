#ifndef hSPIRIT
#define hSPIRIT

#include <iostream>
#include <armadillo>
#include <math.h>
#include <fftw3.h>
#include <omp.h>
#include "ArrayTemplates.cpp"

//enum spirit_method_t {KSPACE, IM_MEM, IM_FFT};
//enum spirit_type_t {SP_FORWARD, SP_FORWARD_DIFF, SP_BACKWARD_DIFF};

void fftshift3(int,int,int,complex<float>*);
void ifft3(arrayND< complex<float> > &);

class SPIRIT {
  public: 
    arrayND< complex<float> > k;
    arrayND< complex<float> > im;
    
    SPIRIT(int,int,int,int);
    void calibrate(int, int, int, arrayND< complex<float> > &, double);
    void prep(int,int,int,int);
    void getcoils(array4D< complex<float> > &,int,float);
    
    
    
    
    /*
   int kdim[3];
    int imdim[3];
    int ncoils;
    
    int imsize3;
    int imsize4;
    long imsize5;
    
    float gridscale;
    float igridscale;
    
    complexf *k;
    complexf *im;
    
    spirit_method_t spirit_method;
    fftwf_plan fftplan;
    
    float global_lam;
    MAT4_cf global_preobj(MAT4_cf &x);
    float global_obj(MAT4_cf &x, MAT4_cf &dx, float t);
    MAT4_cf global_gradient(MAT4_cf &x);
    
    MAT4_cf pocs(MAT4_cf &x);

    SPIRITREG(int,int,int,int,spirit_method_t,float);
    void apply(MAT4_cf&, spirit_type_t);
    void calibrate(int, int, int, MAT4_cf&, float);
    void calibrate(int, int, int, array3D< complex<float> >*, float);
    void prep(int[4]);
    void saveim(char*);
    
    void getcoils(array3D< complex<float> > *);
    
    void forward(MAT4_cf&, int);
    void backward(MAT4_cf&, int);
    void forward_fft(MAT4_cf&, int);
    void backward_fft(MAT4_cf&, int);
    void calibrate_2D(int, int, int, MAT4_cf&, float);
    void calibrate_3D(int, int, int, MAT4_cf&, float);
    */
    
};

#endif
