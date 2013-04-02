
cmake_minimum_required(VERSION 2.8)

project(MRI_RECON)
set(PROGNAME recon_binary)

file(GLOB LOCAL_INCS ${CMAKE_SOURCE_DIR}/*.h)

ADD_DEFINITIONS(-DLINUX -DRECON_64BIT -D_FILE_OFFSET_BITS=64)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-g -Wall -m64 -O2 -fopenmp")
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} "-fopenmp")
set(CMAKE_C_FLAGS "-g -Wno-deprecated")

set(LOCAL_LIBS armadillo acml_mp acml_mv fftw3f_omp fftw3f blitz m)

include_directories(${INCLUDE_DIRECTORIES} ${CMAKE_INCLUDE_PATH})
link_directories(${LINK_DIRECTORIES} ${CMAKE_LIBRARY_PATH})
set(BUILD_SHARED_LIBS OFF)

#--------------------------------------
# Target Library 
#--------------------------------------

add_library(RECONLIB STATIC clear.cpp gating.cpp recon_lib.cpp ArrayTemplates.cpp threshold.cpp temporal_diff.cpp mri_data.cpp spirit.cpp gridFFT.cpp wavelet3D.cpp phantom.cpp)
target_link_libraries(RECONLIB ${LOCAL_LIBS} )

#--------------------------------------
# Support for Command Line Binary
#--------------------------------------
add_executable(recon_binary commandline_recon.cxx ${LOCAL_INCS})
target_link_libraries(recon_binary ${LOCAL_LIBS} RECONLIB)

