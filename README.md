# UW-MRI Reconstruction Toolbox
This is reconstruction code developed for MRI reconstructions, primarily non-cartesian data but will work with arbitrary code. It has a bit of a focus on iterative reconstructions using coil sensitvities derieved from the data (e.g. ESPIRiT ). 

# Highlighted Publications using this Code
* Rivera-Rivera LA, Johnson KM, Turski PA, Wieben O, Schubert T. Measurement of microvascular cerebral blood volume changes over the cardiac cycle with ferumoxytol-enhanced T2* MRI. Magn Reson Med. 2019 Jun;81(6):3588–3598. PMID: 30756424
* Jimenez JE, Strigel RM, Johnson KM, Henze Bancroft LC, Reeder SB, Block WF. Feasibility of high spatiotemporal resolution for an abbreviated 3D radial breast MRI protocol. Magn Reson Med. 2018;80(4):1452–1466. PMCID: PMC6721961
* Johnson KM. Hybrid radial-cones trajectory for accelerated MRI. Magn Reson Med. 2017;77(3):1068–1081. PMCID: PMC5040626
* Kaiser J, Monawer A, Chaudhary R, Johnson KM, Wieben O, Kijowski R, Thelen DG. Accuracy of model-based tracking of knee kinematics and cartilage contact measured by dynamic volumetric MRI. Med Eng Phys. 2016;38(10):1131–1135. PMCID: PMC5035576
* Bauman G, Johnson KM, Bell LC, Velikina JV, Samsonov AA, Nagle SK, Fain SB. 3D pulmonary perfusion MRI with radial ultra-short echo time and spatial-temporal constrained reconstruction. Magn Reson Med. 2015 Feb;73(2):555–564. PMCID: PMC4156934
* Johnson KM, Block WF, Reeder ScottB, Samsonov A. Improved Least Squares MR Image Reconstruction Using Estimates of k-Space Data Consistency. Magn Reson Med. 2012 Jun;67(6):1600–1608. PMCID: PMC3354012

# Software dependencies
* Blitz++ a template array library, used for 3D+ arrays
* Armadillo a linear algebra package. It is highly recomended you install with the Intel MKL package
* HDF5 C++ library
* FFTW with OpenMp support 
* Voro++ a 3D voronoi diagram package (used for sampling density estimates)
* Ochestra 1.8 or higher if using GE library compatibility
* gcc/g++ 4.x for Orchestra code

# Compiling without Orchestra 


# Compiling with GE Orchestra
CC=gcc-4.8 CXX=g++-4.8 FC=gfortran-4.8 cmake ../ -DCMAKE_BUILD_TYPE=Release  -DENABLE_ORCHESTRA="yes" -DOX_INSTALL_DIRECTORY=$SDKTOP -DCMAKE_INSTALL_PREFIX=${HOME}/local/
CC=gcc-4.8 CXX=g++-4.8 FC=gfortran-4.8 cmake ../ -DCMAKE_BUILD_TYPE=Debug -DENABLE_ORCHESTRA="yes" -DOX_INSTALL_DIRECTORY=$SDKTOP -DCMAKE_INSTALL_PREFIX=${HOME}/local/


# Data preparation

# Common Options
