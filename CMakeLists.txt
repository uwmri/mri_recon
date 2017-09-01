cmake_minimum_required(VERSION 3.6)

project(MRI_RECON)
set(PROGNAME recon_binary)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CmakeModules/")
file(GLOB LOCAL_INCS ${CMAKE_SOURCE_DIR}/*.h)

set(CMAKE_FIND_STATIC FIRST) 
set(BUILD_SHARED_LIBS ON)
set(STATIC_LIBS ON)

link_directories( ${LINK_DIRECTORIES} ${CMAKE_LIBRARY_PATH} )


if (DEFINED ENABLE_ORCHESTRA )

#--------------------------------------
# Orchestra (GE Recon Library)
#--------------------------------------

# TOPDIR is the root of the source tree
set(OX_INSTALL_DIRECTORY "/usr/local/orchestra16/" CACHE PATH "Install directory of the GE Orchestra SDK")
set(TOPDIR "${OX_INSTALL_DIRECTORY}/include")
set(LARGE_LIBRARY_TYPE STATIC)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CmakeModules/")

message("Configuring GE SDK Recon Build")
set(PRODUCT EXTERNAL_SDK)

# Include recon libraries exported from the SDK CMake build
set(RECON_LIBRARIES_INCLUDE_FILE ${OX_INSTALL_DIRECTORY}/lib/ReconLibraries.cmake)
if(EXISTS ${RECON_LIBRARIES_INCLUDE_FILE})
    include (${RECON_LIBRARIES_INCLUDE_FILE})
else()
    message("Could not find ${RECON_LIBRARIES_INCLUDE_FILE}")
    message(FATAL_ERROR "Verify that the CMake OX_INSTALL_DIRECTORY option is set correctly")
endif()

# Include SDK build configuration
include (${TOPDIR}/recon/SDK/product.cmake)

# Override 3P libs to include only those libs required by SDK executables
# Original set(OX_3P_LIBS ${DCMTK_LIBS} ${HDF5_LIBS} ${BLITZ_LIBS} ${BOOST_LIBS} ${FFTW_LIBS} ${CLAPACK_LIBS} ${Z_LIBS})
set(OX_3P_LIBS ${DCMTK_LIBS} ${HDF5_LIBS} ${BLITZ_LIBS} ${BOOST_LIBS} ${FFTW_LIBS} ${Z_LIBS})

# This code only requires a few libraries
set( LOCAL_LIBS ${LOCAL_LIBS} ${OX_3P_LIBS} ${OX_OS_LIBS} Hdf5 Common)


else ()

#--------------------------------------
# Non-Orchestra
#--------------------------------------

include(FindBlitz)
include(FindFFTW)
include(FindHDF5)

#Set include directories
include_directories( BEFORE  ${BLITZ_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${HDF5_INCLUDE_DIR}  )
set(LOCAL_LIBS ${LOCAL_LIBS} ${FFTW3_LIBRARIES}  ${BLITZ_LIBRARIES} ${HDF5_LIBRARIES} dl m z )

endif()


ADD_DEFINITIONS(-DLINUX -DRECON_64BIT -D_FILE_OFFSET_BITS=64)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-g -Wall -Wextra -m64 -O2 -fopenmp ")
set(CMAKE_CXX_LINK_FLAGS  ${CMAKE_CXX_LINK_FLAGS} )

# Find static libraries for portability
#SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

# Find Required Libraries
message( STATUS "Check for required libraries ")

# Find Voro++
include(FindVORO)
include_directories( BEFORE ${VORO_INCLUDE_DIR} )
set(LOCAL_LIBS ${LOCAL_LIBS} ${VORO_LIBRARIES} )

#Find Armadillo
#find_package(Armadillo REQUIRED)
include(FindARMADILLO)
set( ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARIES} "-lgfortran")
include_directories( BEFORE ${ARMADILLO_INCLUDE_DIRS} )

MESSAGE(STATUS "   Local Libraries = ${LOCAL_LIBS}")



message(STATUS " ARMADILLO_FOUND  = ${ARMADILLO_FOUND}") 
message(STATUS " ARMADILLO_INCLUDE_DIRS  = ${ARMADILLO_INCLUDE_DIRS}") 
message(STATUS " ARMADILLO_LIBRARIES  = ${ARMADILLO_LIBRARIES}") 
message(STATUS " ARMADILLO_VERSION_MAJOR  = ${ARMADILLO_VERSION_MAJOR}") 
message(STATUS " ARMADILLO_VERSION_MINOR  = ${ARMADILLO_VERSION_MINOR}") 
message(STATUS " ARMADILLO_VERSION_PATCH  = ${ARMADILLO_VERSION_PATCH}") 
message(STATUS " ARMADILLO_VERSION_STRING  = ${ARMADILLO_VERSION_STRING}") 
message(STATUS " ARMADILLO_VERSION_NAME  = ${ARMADILLO_VERSION_NAME}") 


#--------------------------------------
# Target Library 
#--------------------------------------

#Main Recon
add_library(RECONLIB STATIC polynomial_fitting.cpp hdf5_interface.cpp ArrayTemplates.cpp voronoi_dcf.cpp l2reg.cpp clear.cpp gating.cpp recon_lib.cpp threshold.cpp temporal_diff.cpp mri_data.cpp spirit.cpp smsEncode.cpp DCFgridFFT.cpp gridFFT.cpp wavelet3D.cpp phantom.cpp)
target_link_libraries(RECONLIB ${LOCAL_LIBS} ${ARMADILLO_LIBRARIES})

#--------------------------------------
# Support for Command Line Binary
#--------------------------------------
add_executable(recon_binary commandline_recon.cxx ${LOCAL_INCS})
target_link_libraries(recon_binary ${LOCAL_LIBS} RECONLIB ${ARMADILLO_LIBRARIES})

add_executable(ops_test test_array_ops.cpp ${LOCAL_INCS})
target_link_libraries(ops_test ${LOCAL_LIBS} RECONLIB ${ARMADILLO_LIBRARIES})

add_executable(hdf5_test test_hdf5.cpp ${LOCAL_INCS})
target_link_libraries(hdf5_test ${LOCAL_LIBS} RECONLIB ${ARMADILLO_LIBRARIES})

add_executable(voronoi_test test_voronoi.cpp ${LOCAL_INCS})
target_link_libraries(voronoi_test ${LOCAL_LIBS} RECONLIB ${ARMADILLO_LIBRARIES})


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








