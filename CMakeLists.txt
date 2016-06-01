cmake_minimum_required (VERSION 2.8)
project (YIP)
enable_language (Fortran)

# Default to Release
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose the type of build, options are: None Debug Release." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Require MPI
find_package(MPI REQUIRED)
set( CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

# Setup Compile Flags
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fdefault-real-8")
  set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -fdefault-real-8")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -i4 -r8")
  set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -i4 -r8")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Build Library
add_library (fft STATIC "fft.f")

# Build Executable
add_executable (yip "les.f")
target_link_libraries (yip fft)

install ( TARGETS yip RUNTIME DESTINATION "bin")
