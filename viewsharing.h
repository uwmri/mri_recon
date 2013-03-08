#ifndef hVIEWSHARINGLIB
#define hVIEWSHARINGLIB

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <string>
#include <complex>
#include <omp.h>
#include "wavelet3D.h"
#include "ArrayTemplates.hpp"
#include "mri_data.h"
#include "gridFFT.h"
#include "export.h"
#include "tictoc.cpp"
#include "io_templates.cpp"

// View sharing modes
#define VS_NONE 0
#define VS_SLIDING 1
#define VS_TORNADO 2

using namespace std;

#define THREADS 32

class VIEWSHARING {

    public:
        VIEWSHARING(int numarg,char **pstring);
        int a,b,vs_type,vs_shape;                               // width upper, with bottom (shape of the filter)
        int si,sj,sk;
        float fdif;
        int idif;
    char cvs_type[20];

        static void help_message(void);
        void createmask(Array<float,3>&Tw,Array<float,3>&times,int t, int frames);
        void tornado(Array<float,3>&Tw,Array<float,3>&times,int t, int frames);
        void slidingwindow(Array<float,3>&Tw,Array<float,3>&times,int t);

    private:

};

#endif
