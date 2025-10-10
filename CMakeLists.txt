cmake_minimum_required(VERSION 3.17)
project(mri_recon)

project(MRI_RECON)
set(PROGNAME recon_binary)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CmakeModules/")
file(GLOB LOCAL_INCS ${CMAKE_SOURCE_DIR}/*.h)

set(CMAKE_FIND_STATIC FIRST)
set(BUILD_SHARED_LIBS OFF)
set(STATIC_LIBS ON)

link_directories( ${LINK_DIRECTORIES} ${CMAKE_LIBRARY_PATH} )

if ( DISABLE_ORCHESTRA )

	#--------------------------------------
	# Non-Orchestra
	#--------------------------------------

	include(FindBlitz)
	include(FindFFTW)
	include(FindHDF5)

	#Set include directories
	include_directories( BEFORE  ${BLITZ_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${HDF5_INCLUDE_DIR}  )
	set(LOCAL_LIBS ${LOCAL_LIBS} ${FFTW3_LIBRARIES}  ${BLITZ_LIBRARIES} ${HDF5_LIBRARIES} dl m z )

else ( DISABLE_ORCHESTRA )
	message( STATUS "Using orchestra libraries")
	
	#--------------------------------------
	# Orchestra (GE Recon Library)
	#--------------------------------------

	if(NOT ORCHESTRA_SDK_ROOT)
		message(FATAL_ERROR "The ORCHESTRA_SDK_ROOT CMake variable must be set when configuring. Please pass -DORCHESTRA_SDK_ROOT=[/path/to/sdk/root] when configuring.")
	else()
		message(STATUS "ORCHESTRA_SDK_ROOT set to: ${ORCHESTRA_SDK_ROOT}")
	endif()

	option(STATIC_LINK_OPENSSL "Link to the OpenSSL static libraries" ON)

	set(BINARY_OUTPUT_ROOT ${CMAKE_CURRENT_LIST_DIR}/Targets CACHE STRING "Root of output directory")
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BINARY_OUTPUT_ROOT}/${CMAKE_HOST_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}/bin)
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BINARY_OUTPUT_ROOT}/${CMAKE_HOST_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}/lib)
	set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BINARY_OUTPUT_ROOT}/${CMAKE_HOST_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}/lib)

	# Include orchestra platform libraries
	include(${ORCHESTRA_SDK_ROOT}/cmake/orchestra-sdk.cmake)

endif(DISABLE_ORCHESTRA)

set(CMAKE_VERBOSE_MAKEFILE OFF)

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

#Find Armadillo
include(FindARMADILLO)
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
target_compile_features(RECONLIB PUBLIC cxx_std_17)

target_link_libraries(RECONLIB PUBLIC ${ARMADILLO_LIBRARIES} ${LOCAL_LIBS} )
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}Gradwarp)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}Legacy)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}Core)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}CalibrationCommon)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}Foundation)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}Acquisition)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}Control)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}Common)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}ProcessingControl)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}Crucial)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}Dicom)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}Hdf5)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}Math)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}System)
target_link_libraries(RECONLIB PUBLIC ${OX_SDK_IMPORT_PREFIX}FrameDownSampler)


#--------------------------------------
# Support for Command Line Binary
#--------------------------------------
add_executable(uwmri_recon src/commandline_recon.cxx ${LOCAL_INCS})
target_link_libraries(uwmri_recon ${ARMADILLO_LIBRARIES} ${LOCAL_LIBS} RECONLIB)
target_compile_features(uwmri_recon PUBLIC cxx_std_17)


add_executable(fractal_binary src/commandline_fractal.cpp ${LOCAL_INCS})
target_link_libraries(fractal_binary ${ARMADILLO_LIBRARIES} ${LOCAL_LIBS} RECONLIB)
target_compile_features(fractal_binary PUBLIC cxx_std_17)

add_executable(ops_test src/test_array_ops.cpp ${LOCAL_INCS})
target_link_libraries(ops_test ${ARMADILLO_LIBRARIES} ${LOCAL_LIBS} RECONLIB )
target_compile_features(ops_test PUBLIC cxx_std_17)


add_executable(hdf5_test src/test_hdf5.cpp ${LOCAL_INCS})
target_link_libraries(hdf5_test ${ARMADILLO_LIBRARIES} ${LOCAL_LIBS} RECONLIB )
target_compile_features(hdf5_test PUBLIC cxx_std_17)


if(VORO_FOUND)
	add_executable(voronoi_test src/test_voronoi.cpp ${LOCAL_INCS})
	target_link_libraries(voronoi_test ${ARMADILLO_LIBRARIES} ${LOCAL_LIBS} RECONLIB )
	target_compile_features(voronoi_test PUBLIC cxx_std_17)
endif(VORO_FOUND)

#--------------------------------------
# Install
#--------------------------------------

message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")

install( TARGETS uwmri_recon fractal_binary hdf5_test ops_test DESTINATION bin)
if(VORO_FOUND)
	install( TARGETS voronoi_test DESTINATION bin)
endif(VORO_FOUND)

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
