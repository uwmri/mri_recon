# - Try to find Armadillo include dirs and libraries
# Usage of this module as follows:
#
# == Using Armadillo: ==
#
#   find_package( Armadillo RECQUIRED )
#   include_directories(${ARMADILLO_INCLUDE_DIRS})
#   add_executable(foo foo.cc)
#   target_link_libraries(foo ${ARMADILLO_LIBRARIES})
#
#=============================================================================
#
# This module sets the following variables:
#  ARMADILLO_FOUND - set to true if the library is found
#  ARMADILLO_INCLUDE_DIRS - list of required include directories
#  ARMADILLO_LIBRARIES - list of libraries to be linked 
#  ARMADILLO_VERSION_MAJOR - major version number
#  ARMADILLO_VERSION_MINOR - minor version number
#  ARMADILLO_VERSION_PATCH - patch version number
#  ARMADILLO_VERSION_STRING - version number as a string (ex: "1.0.4")
#  ARMADILLO_VERSION_NAME - name of the version (ex: "Antipodean Antileech")


# UNIX paths are standard, no need to write.
find_library(ARMADILLO_LIBRARY  NAMES armadillo
  PATHS "$ENV{ProgramFiles}/Armadillo/lib"  "$ENV{ProgramFiles}/Armadillo/lib64" "$ENV{ProgramFiles}/Armadillo"
  )
  
find_path(ARMADILLO_INCLUDE_DIR
  NAMES armadillo
  PATHS "$ENV{ProgramFiles}/Armadillo/include"
  )


if(ARMADILLO_INCLUDE_DIR)

  # ------------------------------------------------------------------------
  #  Extract version information from <armadillo>
  # ------------------------------------------------------------------------

  # WARNING: Early releases of Armadillo didn't have the arma_version.hpp file.
  # (e.g. v.0.9.8-1 in ubuntu maverick packages (2001-03-15))
  # If the file is missing, set all values to 0  
  set(ARMADILLO_VERSION_MAJOR 0)
  set(ARMADILLO_VERSION_MINOR 0)
  set(ARMADILLO_VERSION_PATCH 0)
  set(ARMADILLO_VERSION_NAME "EARLY RELEASE")

  if(EXISTS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/arma_version.hpp")
	
    # Read and parse armdillo version header file for version number 
    file(READ "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/arma_version.hpp" _armadillo_HEADER_CONTENTS)
    string(REGEX REPLACE ".*#define ARMA_VERSION_MAJOR ([0-9]+).*" "\\1" ARMADILLO_VERSION_MAJOR "${_armadillo_HEADER_CONTENTS}")
    string(REGEX REPLACE ".*#define ARMA_VERSION_MINOR ([0-9]+).*" "\\1" ARMADILLO_VERSION_MINOR "${_armadillo_HEADER_CONTENTS}")
    string(REGEX REPLACE ".*#define ARMA_VERSION_PATCH ([0-9]+).*" "\\1" ARMADILLO_VERSION_PATCH "${_armadillo_HEADER_CONTENTS}")

    # WARNING: The number of spaces before the version name is not one.
    string(REGEX REPLACE ".*#define ARMA_VERSION_NAME\ +\"([0-9a-zA-Z\ _-]+)\".*" "\\1" ARMADILLO_VERSION_NAME "${_armadillo_HEADER_CONTENTS}")
  
  endif(EXISTS "${ARMADILLO_INCLUDE_DIR}/armadillo_bits/arma_version.hpp")

  set(ARMADILLO_VERSION_STRING "${ARMADILLO_VERSION_MAJOR}.${ARMADILLO_VERSION_MINOR}.${ARMADILLO_VERSION_PATCH}")
endif (ARMADILLO_INCLUDE_DIR)





set(ARMA_USE_LAPACK           false)
set(ARMA_USE_BLAS             false)
set(ARMA_USE_ATLAS            false)
set(ARMA_USE_HDF5_ALT         false)
set(ARMA_USE_ARPACK           false)
set(ARMA_USE_EXTERN_CXX11_RNG false)
set(ARMA_USE_SUPERLU          false)


##
## Find LAPACK and BLAS libraries, or their optimised versions
##

if(APPLE)
  
  set(ARMA_OS macos)
  
  set(ARMA_USE_LAPACK true)
  set(ARMA_USE_BLAS   true)
  
  set(ARMA_LIBS ${ARMA_LIBS} "-framework Accelerate")  # or "-framework accelerate" ?
  message(STATUS "MacOS X detected. Added '-framework Accelerate' to compiler flags")
  
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    message(STATUS "Clang compiler on MacOS X detected. Added '-stdlib=libc++' to compiler flags")
  endif()
  
else()

  # MKL outside of Armadillo
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_DONT_USE_WRAPPER ")

#include(FindMKL)
#include_directories( BEFORE ${MKL_INCLUDE_DIRS} )
  
  set(ARMA_OS unix)
  
  include(FindMKL)
  include(ARMA_FindACMLMP)
  include(ARMA_FindACML)
  include(ARMA_FindOpenBLAS)
  include(ARMA_FindATLAS)
  include(ARMA_FindBLAS)
  include(ARMA_FindLAPACK)
  
  message(STATUS "     MKL_FOUND = ${MKL_FOUND}"     )
  message(STATUS "  ACMLMP_FOUND = ${ACMLMP_FOUND}"  )
  message(STATUS "    ACML_FOUND = ${ACML_FOUND}"    )
  message(STATUS "OpenBLAS_FOUND = ${OpenBLAS_FOUND}")
  message(STATUS "   ATLAS_FOUND = ${ATLAS_FOUND}"   )
  message(STATUS "    BLAS_FOUND = ${BLAS_FOUND}"    )
  message(STATUS "  LAPACK_FOUND = ${LAPACK_FOUND}"  )
     
  if(MKL_FOUND OR ACMLMP_FOUND OR ACML_FOUND)
    
    set(ARMA_USE_LAPACK true)
    set(ARMA_USE_BLAS   true)
    
    message(STATUS "")
    message(STATUS "*** If the MKL or ACML libraries are installed in non-standard locations such as")
    message(STATUS "*** /opt/intel/mkl, /opt/intel/composerxe/, /usr/local/intel/mkl")
    message(STATUS "*** make sure the run-time linker can find them.")
    message(STATUS "*** On Linux systems this can be done by editing /etc/ld.so.conf")
    message(STATUS "*** or modifying the LD_LIBRARY_PATH environment variable.")
    message(STATUS "")
    message(STATUS "*** On systems with SELinux enabled (eg. Fedora, RHEL),")
    message(STATUS "*** you may need to change the SELinux type of all MKL/ACML libraries")
    message(STATUS "*** to fix permission problems that may occur during run-time.")
    message(STATUS "*** See README.txt for more information")
    message(STATUS "")
    
    if(MKL_FOUND)
      set(ARMA_LIBS ${ARMA_LIBS} ${MKL_LIBRARIES})
      include_directories( BEFORE ${MKL_INCLUDE_DIRS} )
	  
      if(ACMLMP_FOUND OR ACML_FOUND)
        message(STATUS "*** Intel MKL as well as AMD ACML libraries were found.")
        message(STATUS "*** Using only the MKL library to avoid linking conflicts.")
        message(STATUS "*** If you wish to use ACML instead, please link manually with")
        message(STATUS "*** acml or acml_mp instead of the armadillo wrapper library.")
        message(STATUS "*** Alternatively, remove MKL from your system and rerun")
        message(STATUS "*** Armadillo's configuration using ./configure") 
      endif()
      
    else()
      
      if(ACMLMP_FOUND)
        set(ARMA_LIBS ${ARMA_LIBS} ${ACMLMP_LIBRARIES})
        
        message(STATUS "*** Both single-core and multi-core ACML libraries were found.")
        message(STATUS "*** Using only the multi-core library to avoid linking conflicts.")
      else()
        if(ACML_FOUND)
          set(ARMA_LIBS ${ARMA_LIBS} ${ACML_LIBRARIES})
        endif()
      endif()
      
    endif()
    
  else()
    if(LAPACK_FOUND)
      set(ARMA_USE_LAPACK true)
      set(ARMA_LIBS ${ARMA_LIBS} ${LAPACK_LIBRARIES})
    endif()
      
      
    if(OpenBLAS_FOUND AND ATLAS_FOUND)
      message(STATUS "")
      message(STATUS "*** WARNING: found both OpenBLAS and ATLAS; ATLAS will not be used")
    endif()
    
    if(OpenBLAS_FOUND AND BLAS_FOUND)
      message(STATUS "")
      message(STATUS "*** WARNING: found both OpenBLAS and BLAS; BLAS will not be used")
    endif()
    
    if(OpenBLAS_FOUND)
      
      set(ARMA_USE_BLAS true)
      set(ARMA_LIBS ${ARMA_LIBS} ${OpenBLAS_LIBRARIES})
      
      message(STATUS "")
      message(STATUS "*** If the OpenBLAS library is installed in")
      message(STATUS "*** /usr/local/lib or /usr/local/lib64")
      message(STATUS "*** make sure the run-time linker can find it.")
      message(STATUS "*** On Linux systems this can be done by editing /etc/ld.so.conf")
      message(STATUS "*** or modifying the LD_LIBRARY_PATH environment variable.")
      message(STATUS "")
      
    else()
      
      if(ATLAS_FOUND)
        set(ARMA_USE_ATLAS true)
        set(ARMA_ATLAS_INCLUDE_DIR ${ATLAS_INCLUDE_DIR})
        set(ARMA_LIBS ${ARMA_LIBS} ${ATLAS_LIBRARIES})
        
        message(STATUS "ATLAS_INCLUDE_DIR = ${ATLAS_INCLUDE_DIR}")
      endif()
      
      if(BLAS_FOUND)
        set(ARMA_USE_BLAS true)
        set(ARMA_LIBS ${ARMA_LIBS} ${BLAS_LIBRARIES})
      endif()
      
    endif()
    

  endif()
  
endif()



#Check for header
IF( ARMADILLO_INCLUDE_DIR )
  SET(ARMADILLO_HEADER_FOUND "yes")
  MESSAGE(STATUS "   Found  armadillo header ${ARMADILLO_VERSION_STRING} : ${ARMADILLO_INCLUDE_DIR}")
ELSE( ARMADILLO_INCLUDE_DIR )
  SET(ARMADILLO_HEADER_FOUND "no")
  MESSAGE(FATAL_ERROR "Could not find the armadillo headers")
endif ( ARMADILLO_INCLUDE_DIR)


#Check for library
IF( ARMADILLO_LIBRARY )
#  SET(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARY} ${ARMA_LIBS})
  SET(ARMADILLO_LIBRARIES ${ARMA_LIBS})

  SET(ARMADILLO_FOUND "yes")
  MESSAGE(STATUS "   Found  armadillo ${ARMADILLO_VERSION_STRING} : ${ARMADILLO_LIBRARIES}")
ELSE( ARMADILLO_LIBRARY )
  SET(ARMADILLO_FOUND "no")
  MESSAGE(FATAL_ERROR "Could not find the armadillo library")
endif ( ARMADILLO_LIBRARY)



