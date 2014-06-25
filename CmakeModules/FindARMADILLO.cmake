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

#Find support libraries
include(ARMA_FindACMLMP)
include(ARMA_FindACML)
message(STATUS "ACMLMP_FOUND   = ${ACMLMP_FOUND}")
message(STATUS "  ACML_FOUND   = ${ACML_FOUND}")

if(ACMLMP_FOUND)
      set(ARMA_LIBS ${ARMA_LIBS} ${ACMLMP_LIBRARIES})
      
      message(STATUS "*** Both single-core and multi-core ACML libraries were found.")
      message(STATUS "*** Using only the multi-core library to avoid linking conflicts.")
else()
     if(ACML_FOUND)
        set(ARMA_LIBS ${ARMA_LIBS} ${ACML_LIBRARIES})
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
  SET(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARY} ${ARMA_LIBS})
  SET(ARMADILLO_FOUND "yes")
  MESSAGE(STATUS "   Found  armadillo ${ARMADILLO_VERSION_STRING} : ${ARMADILLO_LIBRARIES}")
ELSE( ARMADILLO_LIBRARY )
  SET(ARMADILLO_FOUND "no")
  MESSAGE(FATAL_ERROR "Could not find the armadillo library")
endif ( ARMADILLO_LIBRARY)



