
cmake_minimum_required(VERSION 2.8)

project(MRI_RECON)
set(PROGNAME recon_binary)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CmakeModules/")
file(GLOB LOCAL_INCS ${CMAKE_SOURCE_DIR}/*.h)

ADD_DEFINITIONS(-DLINUX -DRECON_64BIT -D_FILE_OFFSET_BITS=64)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-g -Wall -Wextra -m64 -O3 -fopenmp ")
set(CMAKE_CXX_LINK_FLAGS  ${CMAKE_CXX_LINK_FLAGS} )

# Find static libraries for portability
#SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

# Find Required Libraries
message( STATUS "Check for required libraries ")
include(FindBlitz)
include(FindFFTW)
include(FindHDF5)
#include(FindARMADILLO)
include(FindVORO)

find_package(Armadillo REQUIRED)

#Set include directories
include_directories( BEFORE ${ARMADILLO_INCLUDE_DIR}  ${BLITZ_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} ${VORO_INCLUDE_DIR} )
set(LOCAL_LIBS ${VORO_LIBRARIES} ${FFTW3_LIBRARIES} ${ARMADILLO_LIBRARIES} ${BLITZ_LIBRARIES} ${HDF5_LIBRARIES} dl m z )

set( LOCAL_LIBS ${LOCAL_LIBS})

MESSAGE(STATUS "   Local Libraries = ${LOCAL_LIBS}")

set(CMAKE_FIND_STATIC FIRST) 
set(BUILD_SHARED_LIBS ON)
set(STATIC_LIBS ON)

link_directories( ${LINK_DIRECTORIES} ${CMAKE_LIBRARY_PATH} )

#--------------------------------------
# Target Library 
#--------------------------------------

#Main Recon
add_library(RECONLIB STATIC polynomial_fitting.cpp hdf5_interface.cpp ArrayTemplates.cpp voronoi_dcf.cpp l2reg.cpp clear.cpp gating.cpp recon_lib.cpp threshold.cpp temporal_diff.cpp mri_data.cpp spirit.cpp smsEncode.cpp DCFgridFFT.cpp gridFFT.cpp wavelet3D.cpp phantom.cpp)
target_link_libraries(RECONLIB ${LOCAL_LIBS})

#--------------------------------------
# Support for Command Line Binary
#--------------------------------------
add_executable(recon_binary commandline_recon.cxx ${LOCAL_INCS})
target_link_libraries(recon_binary ${LOCAL_LIBS} RECONLIB)

add_executable(ops_test test_array_ops.cpp ${LOCAL_INCS})
target_link_libraries(ops_test ${LOCAL_LIBS} RECONLIB)

add_executable(hdf5_test test_hdf5.cpp ${LOCAL_INCS})
target_link_libraries(hdf5_test ${LOCAL_LIBS} RECONLIB)

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








