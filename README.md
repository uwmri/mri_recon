# UW-MRI Reconstruction Toolbox
This is reconstruction code developed for MRI reconstructions, primarily non-cartesian data but will work with arbitrary code. It has a bit of a focus on iterative reconstructions using coil sensitvities derieved from the data (e.g. ESPIRiT ). 

# Highlighted Publications using this Code
* Rivera-Rivera LA, Johnson KM, Turski PA, Wieben O, Schubert T. Measurement of microvascular cerebral blood volume changes over the cardiac cycle with ferumoxytol-enhanced T2* MRI. Magn Reson Med. 2019 Jun;81(6):3588-3598. doi: 10.1002/mrm.27670
* Jimenez JE, Strigel RM, Johnson KM, Henze Bancroft LC, Reeder SB, Block WF. Feasibility of high spatiotemporal resolution for an abbreviated 3D radial breast MRI protocol. Magn Reson Med. 2018 Oct;80(4):1452-1466. doi: 10.1002/mrm.27137. 

# Software dependencies
* Blitz++ a template array library, used for 3D+ arrays
* Armadillo a linear algebra package. It is highly recomended you install with the Intel MKL package
* HDF5 C++ library
* FFTW with OpenMp support 
* [Optional]Voro++ a 3D voronoi diagram package (used for sampling density estimates)
* [Optional]Ochestra 1.8 or higher if using GE library compatibility
* [Optional]gcc/g++ 4.x for Orchestra code compatibility

# Compiling without Orchestra 


# Compiling with GE Orchestra (Example with gcc 4.8)
CC=gcc-4.8 CXX=g++-4.8 FC=gfortran-4.8 cmake ../ -DCMAKE_BUILD_TYPE=Release -DENABLE_ORCHESTRA="yes" -DOX_INSTALL_DIRECTORY=$SDKTOP -DCMAKE_INSTALL_PREFIX=${HOME}/local/
CC=gcc-4.8 CXX=g++-4.8 FC=gfortran-4.8 cmake ../ -DCMAKE_BUILD_TYPE=Debug -DENABLE_ORCHESTRA="yes" -DOX_INSTALL_DIRECTORY=$SDKTOP -DCMAKE_INSTALL_PREFIX=${HOME}/local/


# Data preparation

# Common Options
