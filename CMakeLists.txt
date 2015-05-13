
cmake_minimum_required(VERSION 2.8)

project(MRI_RECON)
set(PROGNAME recon_binary)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CmakeModules/")
file(GLOB LOCAL_INCS ${CMAKE_SOURCE_DIR}/*.h)

ADD_DEFINITIONS(-DLINUX -DRECON_64BIT -D_FILE_OFFSET_BITS=64)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-g -Wall -Wextra -m64 -J16 -O3 -fopenmp")
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} )

# Find Required Libraries
message( STATUS "Check for required libraries ")
include(FindBlitz)
include(FindFFTW)
include(FindHDF5)
include(FindARMADILLO)

#Set include directories
include_directories( BEFORE ${ARMADILLO_INCLUDE_DIR}  ${BLITZ_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} )
set(LOCAL_LIBS voro++ ${FFTW3_LIBRARIES} ${ARMADILLO_LIBRARIES} ${BLITZ_LIBRARIES} ${HDF5_LIBRARIES} dl m z )

set(CMAKE_FIND_STATIC FIRST) 
set(BUILD_SHARED_LIBS ON)
set(STATIC_LIBS ON)

link_directories(${LINK_DIRECTORIES} ${CMAKE_LIBRARY_PATH} )

#--------------------------------------
# Target Library 
#--------------------------------------

#Main Recon
add_library(RECONLIB STATIC polynomial_fitting.cpp hdf5_interface.cpp ArrayTemplates.cpp voronoi_dcf.cpp l2reg.cpp clear.cpp gating.cpp recon_lib.cpp threshold.cpp temporal_diff.cpp mri_data.cpp spirit.cpp DCFgridFFT.cpp gridFFT.cpp wavelet3D.cpp phantom.cpp)
target_link_libraries(RECONLIB ${LOCAL_LIBS})

#--------------------------------------
# Support for Command Line Binary
#--------------------------------------
add_executable(recon_binary commandline_recon.cxx ${LOCAL_INCS})
target_link_libraries(recon_binary ${LOCAL_LIBS} RECONLIB)

add_executable(toy_recon toy_recon.cpp ${LOCAL_INCS})
target_link_libraries(toy_recon ${LOCAL_LIBS} RECONLIB)

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








