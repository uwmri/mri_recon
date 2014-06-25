# - Find Blitz library
# Find the native Blitz includes and library
# This module defines
#  Blitz_INCLUDE_DIR, where to find tiff.h, etc.
#  Blitz_LIBRARIES, libraries to link against to use Blitz.
#  Blitz_FOUND, If false, do not try to use Blitz.
# also defined, but not for general use are
#  Blitz_LIBRARY, where to find the Blitz library.


#Look for libraries
SET(BLITZ_NAMES ${BLITZ_NAMES} blitz)
FIND_LIBRARY( BLITZ_LIBRARY  NAMES ${BLITZ_NAMES} )
FIND_PATH( BLITZ_INCLUDE_DIR blitz/blitz.h )

#Check for header
IF( BLITZ_INCLUDE_DIR )
  SET(BLITZ_HEADER_FOUND "yes")
  MESSAGE(STATUS "   Found the BLITZ headers: ${BLITZ_INCLUDE_DIR}")
ELSE( BLITZ_INCLUDE_DIR )
  SET(BLITZ_HEADER_FOUND "no")
  MESSAGE(FATAL_ERROR "Could not find the BLITZ headers")
endif ( BLITZ_INCLUDE_DIR)


#Check for library
IF( BLITZ_LIBRARY )
  SET(BLITZ_LIBRARIES ${BLITZ_LIBRARY})
  SET(BLITZ_FOUND "yes")
  MESSAGE(STATUS "   Found the BLITZ library: ${BLITZ_LIBRARIES}")
ELSE( BLITZ_LIBRARY )
  SET(BLITZ_FOUND "no")
  MESSAGE(FATAL_ERROR "Could not find the BLITZ library")
endif ( BLITZ_LIBRARY)


