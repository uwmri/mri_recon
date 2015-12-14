# - Find Voro++ library
# Find the native Blitz includes and library
# This module defines
#  VORO_INCLUDE_DIR
#  VORO_LIBRARIES
#  VORO_FOUND


#Look for libraries
FIND_LIBRARY( VORO_LIBRARY  NAMES libvoro++.a )
FIND_PATH( VORO_INCLUDE_DIR voro++/voro++.hh )

#Check for header
IF( VORO_INCLUDE_DIR )
  SET(VORO_HEADER_FOUND "yes")
  MESSAGE(STATUS "   Found the Voro++ headers: ${VORO_INCLUDE_DIR}")
ELSE( VORO_INCLUDE_DIR )
  SET(VORO_HEADER_FOUND "no")
  MESSAGE(FATAL_ERROR "Could not find the Voro++ headers")
endif ( VORO_INCLUDE_DIR)


#Check for library
IF( VORO_LIBRARY )
  SET(VORO_LIBRARIES ${VORO_LIBRARY})
  SET(VORO_FOUND "yes")
  MESSAGE(STATUS "   Found the Voro++ library: ${VORO_LIBRARIES}")
ELSE( VORO_LIBRARY )
  SET(VORO_FOUND "no")
  MESSAGE(FATAL_ERROR "Could not find the Voro++ library")
endif ( VORO_LIBRARY)


