# UW-MRI Reconstruction Toolbox
This is reconstruction code developed for MRI reconstructions, primarily non-cartesian data but will work with arbitrary code. It has a focus on iterative reconstructions using coil sensitvities derieved from the data (e.g. ESPIRiT ). 

# Highlighted Publications using this Code
* Rivera-Rivera LA, Johnson KM, Turski PA, Wieben O, Schubert T. Measurement of microvascular cerebral blood volume changes over the cardiac cycle with ferumoxytol-enhanced T2* MRI. Magn Reson Med. 2019 Jun;81(6):3588–3598. PMID: 30756424 [Link](https://www.ncbi.nlm.nih.gov/pubmed/30756424)
* Jimenez JE, Strigel RM, Johnson KM, Henze Bancroft LC, Reeder SB, Block WF. Feasibility of high spatiotemporal resolution for an abbreviated 3D radial breast MRI protocol. Magn Reson Med. 2018;80(4):1452–1466. PMCID: PMC6721961 [Link](https://www.ncbi.nlm.nih.gov/pubmed/29446125)
* Johnson KM. Hybrid radial-cones trajectory for accelerated MRI. Magn Reson Med. 2017;77(3):1068–1081. PMCID: PMC5040626 [Link](https://www.ncbi.nlm.nih.gov/pubmed/27017991)
* Kaiser J, Monawer A, Chaudhary R, Johnson KM, Wieben O, Kijowski R, Thelen DG. Accuracy of model-based tracking of knee kinematics and cartilage contact measured by dynamic volumetric MRI. Med Eng Phys. 2016;38(10):1131–1135. PMCID: PMC5035576 [Link](https://www.ncbi.nlm.nih.gov/pubmed/27387902)
* Bauman G, Johnson KM, Bell LC, Velikina JV, Samsonov AA, Nagle SK, Fain SB. 3D pulmonary perfusion MRI with radial ultra-short echo time and spatial-temporal constrained reconstruction. Magn Reson Med. 2015 Feb;73(2):555–564. PMCID: PMC4156934 [Link](https://www.ncbi.nlm.nih.gov/pubmed/24604452)
* Johnson KM, Block WF, Reeder ScottB, Samsonov A. Improved Least Squares MR Image Reconstruction Using Estimates of k-Space Data Consistency. Magn Reson Med. 2012 Jun;67(6):1600–1608. PMCID: PMC3354012 [Link](https://www.ncbi.nlm.nih.gov/pubmed/22135155) 

# Software dependencies
* Blitz++ a template array library, used for 3D+ arrays [Link](https://github.com/blitzpp/blitz)
* Armadillo a linear algebra package. It is highly recomended you install with the Intel MKL package [Link](http://arma.sourceforge.net/docs.html) 
* HDF5 C++ library [Link](https://www.hdfgroup.org/solutions/hdf5/). Install From source with "--enable-fortran --enable-cxx --enable-static --disable-shared" 
* FFTW3 with OpenMp support [Link](http://www.fftw.org/)  Install From source with "--enable-float --enable-threads --enable-openmp"
* [Optional]Voro++ a 3D voronoi diagram package (used for sampling density estimates). [Link](http://math.lbl.gov/voro++/)

#  GE Orchestra Support
For Orchestra support, you must additionally:
* Install Ochestra 1.8 or higher and set the OX_INSTALL_DIR
* Compile with gcc/g++ 4.x (e.g. 4.8)
## Compiling with GE Orchestra (Example with gcc 4.8)
CC=gcc-4.8 CXX=g++-4.8 FC=gfortran-4.8 cmake ../ -DCMAKE_BUILD_TYPE=Release -DENABLE_ORCHESTRA="yes" CMAKE_INSTALL_PREFIX=${HOME}/local/
CC=gcc-4.8 CXX=g++-4.8 FC=gfortran-4.8 cmake ../ -DCMAKE_BUILD_TYPE=Debug -DENABLE_ORCHESTRA="yes" -DCMAKE_INSTALL_PREFIX=${HOME}/local/

# Coding style
Code is C++99 to support current coding compatibility of Vendor libraries. Style is based on 2 space indentation according to Google's style ( https://google.github.io/styleguide/cppguide.html ) but doesn't have column width limit. 

# Data preparation
This code can be called from external c++ programs to avoid an explicit file. 

# Common Options
Options should be visible by typing "uwmri_recon -h"

### Simple recon using sensitivity maps from the center of k-space
```bash
uwmri_recon -f MRI_Data.h5 -pils
```
### Simple iterative SENSE 
```bash
uwmri_recon -f MRI_Data.h5 -cg -max_iter 20
```
### Simple L1-Wavelet regularized SENSE with FISTA 
```bash
uwmri_recon -f MRI_Data.h5 -fista -max_iter 100 -threshold_type fraction -thresh 0.1
```
### Gated Local Low Rank using Time info in the header to gate the data into 10 frames
```bash
uwmri_recon -f MRI_Data.h5 -fista -max_iter 100 -gating_type time -rcframes 10 -clear_alpha_time 0.001  
```






