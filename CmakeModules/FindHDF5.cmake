# - Find HDF5_cpp library
# Find the native Blitz includes and library
# This module defines
#  HDF5_INCLUDE_DIR
#  HDF5_LIBRARIES
#  HDF5_FOUND



#Look for libraries
FIND_LIBRARY( HDF5CPP_LIBRARY  NAMES libhdf5_cpp.a hdf5_cpp PATHS /usr/local/ /usr/local/lib64/ /usr/lib/ /usr/lib64/ )
FIND_LIBRARY( HDF5_LIBRARY  NAMES libhdf5.a hdf5 PATHS /usr/local/ /usr/local/lib64/ /usr/lib/ /usr/lib64/ )
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
IF( HDF5_LIBRARY AND HDF5CPP_LIBRARY)
  SET(HDF5_LIBRARIES  ${HDF5CPP_LIBRARY} ${HDF5_LIBRARY})
  SET(HDF5_FOUND "yes")
  MESSAGE(STATUS "   Found the HDF5 library: ${HDF5_LIBRARIES}")
ELSE( HDF5_LIBRARY AND HDF5CPP_LIBRARY )
  SET(HDF5_FOUND "no")
  MESSAGE( STATUS " Could not find HDF5")
ENDIF ( HDF5_LIBRARY AND HDF5CPP_LIBRARY )


