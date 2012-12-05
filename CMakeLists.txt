#
# Written by Orhan Unal, November 2012
#

cmake_minimum_required(VERSION 2.8)

project(PCVIPR)
set(PROGNAME recon_binary)

file(GLOB LOCAL_INCS ${CMAKE_SOURCE_DIR}/*.h)

ADD_DEFINITIONS(-DLINUX -DRECON_64BIT -D_FILE_OFFSET_BITS=64)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-g -Wall -m64 -O3 -fopenmp")
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} "-fopenmp")
set(CMAKE_C_FLAGS "-g -Wno-deprecated")

set(LOCAL_LIBS armadillo acml_mp acml_mv fftw3f_omp fftw3f blitzwave blitz pthread m)

# Add Libraries
add_library(PHANTOM_LIB phantom.cpp)
add_library(WAVELET_LIB wavelet3D.cpp)
add_library(ARRAY_LIB ArrayTemplates.cpp)

# Set Soruce
set (SOURCES 
     recon.cxx 
     mri_data.cpp 
     gridFFT.cpp 
	 spirit.cpp 
     temporal_diff.cpp 
     recon_lib.cpp
	 threshold.cpp)

include_directories(${INCLUDE_DIRECTORIES} ${CMAKE_INCLUDE_PATH})
link_directories(${LINK_DIRECTORIES} ${CMAKE_LIBRARY_PATH})

add_executable(${PROGNAME} ${SOURCES} ${LOCAL_INCS})
target_link_libraries(${PROGNAME} WAVELET_LIB PHANTOM_LIB ARRAY_LIB ${LOCAL_LIBS})

