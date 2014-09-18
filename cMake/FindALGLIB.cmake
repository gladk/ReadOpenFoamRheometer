# - Try to find agllib 
# Once done this will define
#  ALGLIB_INCLUDE_DIRS - The octave include directories
#  ALGLIB_LIB - The libraries needed to use octave
#  ALGLIB_FOUND -  True if ALGLIB found.

find_path (ALGLIB_INCLUDE_DIRS ap.h)
find_library (ALGLIB_LIB NAMES alglib)

include (FindPackageHandleStandardArgs)

find_package_handle_standard_args (ALGLIB DEFAULT_MSG ALGLIB_LIB ALGLIB_INCLUDE_DIRS)

mark_as_advanced (ALGLIB_LIB ALGLIB_INCLUDE_DIRS)
