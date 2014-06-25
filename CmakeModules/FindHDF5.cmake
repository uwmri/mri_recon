# - Find HDF5_cpp library
# Find the native Blitz includes and library
# This module defines
#  HDF5_INCLUDE_DIR
#  HDF5_LIBRARIES
#  HDF5_FOUND

#Look for libraries
SET(HDF5_NAMES ${HDF5_NAMES} hdf5_cpp )
FIND_LIBRARY( HDF5_LIBRARY  NAMES ${HDF5_NAMES} )
FIND_PATH( HDF5_INCLUDE_DIR H5Cpp.h )


#Check for header
IF( HDF5_INCLUDE_DIR )
  SET(HDF5_HEADER_FOUND "yes")
  MESSAGE(STATUS "   Found the HDF5 headers: ${HDF5_INCLUDE_DIR}")
ELSE( HDF5_INCLUDE_DIR )
  SET(HDF5_HEADER_FOUND "no")
  MESSAGE(FATAL_ERROR "Could not find the HDF5 headers")
endif ( HDF5_INCLUDE_DIR)


#Check for library
IF( HDF5_LIBRARY )
  SET(HDF5_LIBRARIES ${HDF5_LIBRARY})
  SET(HDF5_FOUND "yes")
  MESSAGE(STATUS "   Found the HDF5 library: ${HDF5_LIBRARIES}")
ELSE( HDF5_LIBRARY )
  SET(BLITZ_FOUND "no")
  MESSAGE(FATAL_ERROR "Could not find the BLITZ library")
endif ( HDF5_LIBRARY)


