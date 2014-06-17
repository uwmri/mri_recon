
cmake_minimum_required(VERSION 2.8)

project(MRI_RECON)
set(PROGNAME recon_binary)

file(GLOB LOCAL_INCS ${CMAKE_SOURCE_DIR}/*.h)

ADD_DEFINITIONS(-DLINUX -DRECON_64BIT -D_FILE_OFFSET_BITS=64)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-g -Wall -Wextra -m64 -J16 -O2 -fopenmp")
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} )
set(CMAKE_C_FLAGS "-g -Wno-deprecated")

#Temp for HD5 
include_directories( BEFORE "/export/home/kmjohnso/local/include/")

set(LOCAL_LIBS armadillo acml_mp acml_mv fftw3f_omp fftw3f blitz m hdf5_cpp)

set(CMAKE_FIND_STATIC FIRST) 
set(BUILD_SHARED_LIBS ON)
set(STATIC_LIBS ON)

link_directories(${LINK_DIRECTORIES} ${CMAKE_LIBRARY_PATH} )

#--------------------------------------
# Target Library 
#--------------------------------------

#Main Recon
add_library(RECONLIB STATIC polynomial_fitting.cpp hdf5_interface.cpp ArrayTemplates.cpp l2reg.cpp clear.cpp gating.cpp recon_lib.cpp threshold.cpp temporal_diff.cpp mri_data.cpp spirit.cpp gridFFT.cpp wavelet3D.cpp phantom.cpp)
target_link_libraries(RECONLIB ${LOCAL_LIBS})

#--------------------------------------
# Support for Command Line Binary
#--------------------------------------
add_executable(recon_binary commandline_recon.cxx ${LOCAL_INCS})
target_link_libraries(recon_binary ${LOCAL_LIBS} RECONLIB)

#--------------------------------------
# Install 
#--------------------------------------

message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")

install( TARGETS recon_binary DESTINATION bin)
install( TARGETS RECONLIB DESTINATION lib)
install( DIRECTORY ./ DESTINATION include 
	FILES_MATCHING PATTERN "*.h" 
	PATTERN "*.hpp" 
	PATTERN "html" EXCLUDE 
	PATTERN ".git" EXCLUDE
	PATTERN "Cmake*" EXCLUDE
	PATTERN "build" EXCLUDE
	PATTERN ".settings" EXCLUDE
	PATTERN ".extern*" EXCLUDE) 








