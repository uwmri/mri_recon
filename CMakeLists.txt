cmake_minimum_required(VERSION 3.13)

project(MRI_RECON)
set(PROGNAME recon_binary)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CmakeModules/")
file(GLOB LOCAL_INCS ${CMAKE_SOURCE_DIR}/*.h)

set(CMAKE_FIND_STATIC FIRST)
set(BUILD_SHARED_LIBS OFF)
set(STATIC_LIBS ON)

link_directories( ${LINK_DIRECTORIES} ${CMAKE_LIBRARY_PATH} )

if ( ENABLE_ORCHESTRA )

	message( STATUS "Using orchestra libraries")
	
	#--------------------------------------
	# Orchestra (GE Recon Library)
	#--------------------------------------

	if(EXISTS $ENV{OX_INSTALL_DIRECTORY})
		set( OX_INSTALL_DIRECTORY $ENV{OX_INSTALL_DIRECTORY})
		message("Install directory of the GE Orchestra SDK = $ENV{OX_INSTALL_DIRECTORY}" )
	else()

		message("Looking for /lib/libFlex.a to find Orchestra")
		find_path(OX_INSTALL_DIRECTORY NAMES /lib/libFlex.a  PATHS /usr/local/orchestra )
		if(EXISTS ${OX_INSTALL_DIRECTORY})
			message("Could not find Orchestra by searching for /lib/libFlex.a")
			message(FATAL_ERROR "Verify orchestra install or update CMakePath")
		endif()
	endif()

	# TOPDIR is the root of the source tree
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
	#include (productUW.cmake)
	include (${TOPDIR}/recon/SDK/product.cmake)

	# Override 3P libs to include only those libs required by SDK executables
	#set(OX_3P_LIBS ${DCMTK_LIBS} ${HDF5_LIBS} ${BLITZ_LIBS} ${BOOST_LIBS} ${FFTW_LIBS} ${CLAPACK_LIBS} ${ITK_LIBS} ${Z_LIBS})
  	set(OX_3P_LIBS ${DCMTK_LIBS} ${HDF5_LIBS} ${BLITZ_LIBS} ${BOOST_LIBS} ${FFTW_LIBS} ${ITK_LIBS} ${Z_LIBS})

	message("Using Orchestra 3p: ${OX_3P_LIBS} ")
	message("Using Orchestra Os: ${OX_OS_LIBS} ")
	message("Using Orchestra Includes: ${OX_3P_INCLUDES} ")
	

	# This code only requires a few libraries
	set( LOCAL_LIBS ${LOCAL_LIBS} ${OX_3P_LIBS} ${OX_OS_LIBS} Common)
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -m64")
else (ENABLE_ORCHESTRA)

	#--------------------------------------
	# Non-Orchestra
	#--------------------------------------

	include(FindBlitz)
	include(FindFFTW)
	include(FindHDF5)

	#Set include directories
	include_directories( BEFORE  ${BLITZ_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${HDF5_INCLUDE_DIR}  )
	set(LOCAL_LIBS ${LOCAL_LIBS} ${FFTW3_LIBRARIES}  ${BLITZ_LIBRARIES} ${HDF5_LIBRARIES} dl m z )

	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
	#set(CMAKE_CXX_LINK_FLAGS  ${CMAKE_CXX_LINK_FLAGS} )

endif(ENABLE_ORCHESTRA)

# Find Required Libraries
message( STATUS "Check for required libraries ")

# Find Voro++
include(FindVORO)
if (VORO_FOUND)
	include_directories( BEFORE ${VORO_INCLUDE_DIR} )
	set(LOCAL_LIBS ${LOCAL_LIBS} ${VORO_LIBRARIES} )
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_VORO")
else (VORO_FOUND)
	message(WARNING " Compiling without Voronoi support - Install Voro++ to enable")
endif (VORO_FOUND)


# Find MKL
set(MKL_INTERFACE lp64)
set(MKL_LINK static)
set(MKL_THREADING gnu_thread)
find_package(MKL CONFIG REQUIRED)
message(STATUS "Found MKL ${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input
message(STATUS "Include MKL ${MKL_INCLUDE_DIRS}") #Provides available list of targets based on input
message(STATUS "MKL Libs ${MKL_LIBRARIES}") #Provides available list of targets based on input
message(STATUS "MKL Link Link $MKL::MKL") #Provides available list of targets based on input
message(STATUS "MKL ${MKL_INCLUDE}")

#Find Armadillo
find_package(Armadillo REQUIRED)
include_directories( BEFORE ${ARMADILLO_INCLUDE_DIRS} )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_DONT_USE_WRAPPER -DARMA_USE_LAPACK -DARMA_USE_BLAS -DARMA_USE_HDF5 -DARMA_USE_OPENMP")

MESSAGE(STATUS "   Local Libraries = ${LOCAL_LIBS}")
message(STATUS " ARMADILLO_FOUND  = ${ARMADILLO_FOUND}")
message(STATUS " ARMADILLO_INCLUDE_DIRS  = ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS " ARMADILLO_LIBRARIES  = ${ARMADILLO_LIBRARIES}")
message(STATUS " ARMADILLO_VERSION_MAJOR  = ${ARMADILLO_VERSION_MAJOR}")
message(STATUS " ARMADILLO_VERSION_MINOR  = ${ARMADILLO_VERSION_MINOR}")
message(STATUS " ARMADILLO_VERSION_PATCH  = ${ARMADILLO_VERSION_PATCH}")
message(STATUS " ARMADILLO_VERSION_STRING  = ${ARMADILLO_VERSION_STRING}")
message(STATUS " ARMADILLO_VERSION_NAME  = ${ARMADILLO_VERSION_NAME}")

find_package(OpenMP)
if (OPENMP_FOUND)
	message(STATUS "Open CXX Flags ${OpenMP_C_FLAGS}" )
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# #TEMP EXPORT ALL
# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()

#--------------------------------------
# Target Library
#--------------------------------------

#Main recon library
set(RECONLIB_SOURCES 
	src/polynomial_fitting.cpp
	src/hdf5_interface.cpp
	src/ArrayTemplates.cpp
	src/l2reg.cpp
	src/clear.cpp
	src/gating.cpp
	src/recon_lib.cpp
	src/threshold.cpp
	src/temporal_diff.cpp
	src/mri_data.cpp
	src/spirit.cpp
	src/smsEncode.cpp
	src/DCFgridFFT.cpp
	src/gridFFT.cpp
	src/gridFFT_CoilThreaded.cpp
	src/wavelet3D.cpp
	src/phantom.cpp)

# Only add if library is found
if (VORO_FOUND)
	set(RECONLIB_SOURCES ${RECONLIB_SOURCES} src/voronoi_dcf.cpp)
endif(VORO_FOUND)

add_library(RECONLIB STATIC ${RECONLIB_SOURCES}	)
target_link_libraries(RECONLIB ${LOCAL_LIBS} )

#--------------------------------------
# Support for Command Line Binary
#--------------------------------------

add_executable(uwmri_recon src/commandline_recon.cxx ${LOCAL_INCS})
target_link_libraries(uwmri_recon PUBLIC ${LOCAL_LIBS} RECONLIB)
target_compile_options(uwmri_recon PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(uwmri_recon PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(uwmri_recon PUBLIC $<LINK_ONLY:MKL::MKL>)

add_executable(fractal_binary src/commandline_fractal.cpp ${LOCAL_INCS})
target_link_libraries(fractal_binary PUBLIC ${LOCAL_LIBS} RECONLIB)
target_compile_options(fractal_binary PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(fractal_binary PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(fractal_binary PUBLIC $<LINK_ONLY:MKL::MKL>)

add_executable(ops_test src/test_array_ops.cpp ${LOCAL_INCS})
target_link_libraries(ops_test PUBLIC ${LOCAL_LIBS} RECONLIB)
target_compile_options(ops_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(ops_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(ops_test PUBLIC $<LINK_ONLY:MKL::MKL>)

add_executable(hdf5_test src/test_hdf5.cpp ${LOCAL_INCS})
target_link_libraries(hdf5_test PUBLIC ${LOCAL_LIBS} RECONLIB)
target_compile_options(hdf5_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(hdf5_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(hdf5_test PUBLIC $<LINK_ONLY:MKL::MKL>)

if(VORO_FOUND)
	add_executable(voronoi_test src/test_voronoi.cpp ${LOCAL_INCS})
	target_link_libraries(voronoi_test PUBLIC ${LOCAL_LIBS} RECONLIB)
	target_compile_options(voronoi_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
	target_include_directories(voronoi_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
	target_link_libraries(voronoi_test PUBLIC $<LINK_ONLY:MKL::MKL>)
endif(VORO_FOUND)

#--------------------------------------
# Install
#--------------------------------------

message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")

install( TARGETS uwmri_recon fractal_binary DESTINATION bin)
install( TARGETS RECONLIB DESTINATION lib)
install( DIRECTORY ./src/ DESTINATION include
	FILES_MATCHING PATTERN "src/*.h"
	PATTERN "src/*.hpp"
	PATTERN "html" EXCLUDE
	PATTERN ".git" EXCLUDE
	PATTERN "Cmake*" EXCLUDE
	PATTERN "build" EXCLUDE
	PATTERN ".settings" EXCLUDE
	PATTERN ".extern*" EXCLUDE)
