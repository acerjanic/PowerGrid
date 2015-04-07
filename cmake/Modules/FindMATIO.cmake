# Find the MATIO headers and library.
#
# MATIO_INCLUDE_DIRS - where to find matio.h, etc.
# MATIO_LIBRARIES - List of libraries.
# MATIO_FOUND - True if matio found.
# Look for the header file.
FIND_PATH(MATIO_INCLUDE_DIR NAMES matio.h
		PATHS /shared/mrfil-data/Software/MRFIL-Modules/PowerGridSupport/include )
MARK_AS_ADVANCED(MATIO_INCLUDE_DIR)
# Look for the library.
FIND_LIBRARY(MATIO_LIBRARY NAMES matio
			   PATH_SUFFIXES lib64 libs lib
		           PATHS  /usr/local /usr /shared/mrfil-data/Software/MRFIL-Modules/PowerGridSupport/lib )
MARK_AS_ADVANCED(MATIO_LIBRARY)
# handle the QUIETLY and REQUIRED arguments and set MATIO_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MATIO DEFAULT_MSG MATIO_LIBRARY MATIO_INCLUDE_DIR)
IF(MATIO_FOUND)
SET(MATIO_LIBRARIES ${MATIO_LIBRARY} ${HDF5_LIBRARIES})
SET(MATIO_INCLUDE_DIRS ${MATIO_INCLUDE_DIR} ${HDF5_INCLUDE_DIR})
ELSE(MATIO_FOUND)
SET(MATIO_LIBRARIES)
SET(MATIO_INCLUDE_DIRS)
ENDIF(MATIO_FOUND)
