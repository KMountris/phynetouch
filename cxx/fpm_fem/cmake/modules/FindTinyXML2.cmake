# - Find TinyXML2
#
# Find the native TinyXML2 includes and add_library
#
# TinyXML2_FOUND - True if TinyXML2 found.
# TinyXML2_INCLUDE_DIR - Where to find tinyxml2.h, etc.
# TinyXML2_INCLUDE_DIRS - Set if found tinyxml2.h, etc.
# TinyXML2_LIBRARIES - List of libraries when using TinyXML2.
#

INCLUDE( "FindPackageHandleStandardArgs" )

IF (MSVC)
SET( TinyXML2_DIRS "C:/Program Files (x86)" CACHE PATH "Directory containing TinyXML2 directories" )
ELSE ()
SET( TinyXML2_DIRS "/usr/local" CACHE PATH "Directory containing TinyXML2 directories" )
ENDIF (MSVC)

#Find header files.
FIND_PATH( TinyXML2_INCLUDE_DIR "tinyxml2.h" PATHS ${TinyXML2_DIRS}/include PATH_SUFFIXES "tinyxml2" )

#Find libraries.
FIND_LIBRARY( TinyXML2_LIBRARY NAMES "tinyxml2" PATHS ${TinyXML2_DIRS}/lib PATH_SUFFIXES "tinyxml2" )

# handle the QUIETLY and REQUIRED arguments and set TinyXML2_FOUND to TRUE
# if all listed variables are TRUE.

FIND_PACKAGE_HANDLE_STANDARD_ARGS( "TinyXML2" DEFAULT_MSG TinyXML2_INCLUDE_DIR TinyXML2_LIBRARY )

IF(TinyXML2_FOUND)
    SET( TinyXML2_INCLUDE_DIRS  ${TinyXML2_INCLUDE_DIR} )
    SET( TinyXML2_LIBRARIES  ${TinyXML2_LIBRARY} )
ELSE(TinyXML2_FOUND)
    MESSAGE( FATAL_ERROR "Could not locate TinyXML2..." )
ENDIF()
