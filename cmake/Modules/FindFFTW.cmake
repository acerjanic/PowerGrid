# - Locate FFTW library
# Defines:
#
#  FFTW_FOUND
#  FFTW_INCLUDE_DIR
#  FFTW_INCLUDE_DIRS (not cached)
#  FFTW_LIBRARIES

#find_path(FFTW_INCLUDE_DIR fftw3.h)
#find_library(FFTW_LIBRARIES NAMES fftw3)


if(NOT $ENV{FFTW_DIR})
set(FFTW_INCLUDE_DIR $ENV{FFTW_INC})
set(FFTW_LIBRARIES $ENV{FFTW_DIR}/libfftw3.a)
else()
find_path(FFTW_INCLUDE_DIR fftw3.h)
find_library(FFTW_LIBRARIES NAMES fftw3)
endif(NOT $ENV{FFTW_DIR})

set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW DEFAULT_MSG FFTW_INCLUDE_DIR FFTW_LIBRARIES)

mark_as_advanced(FFTW_FOUND FFTW_INCLUDE_DIR FFTW_LIBRARIES)
