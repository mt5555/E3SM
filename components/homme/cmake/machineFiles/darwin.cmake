# CMake initial cache file for macOS (Intel or Apple Silicon)
# tested with gcc/gfortran, openmpi, netcdf, netcdf-fortran from HOMEBREW
#
#   brew install cmake gcc open-mpi netcdf netcdf-fortran
#
# NetCDF is located with nf-config/nc-config, so any Homebrew prefix
# works (/usr/local on Intel, /opt/homebrew on Apple Silicon).
#
# openmpi by default will use Mac's clang based compiler. to use GCC
# (substitute your gcc version):
#
# setenv OMPI_CXX g++-16
# setenv OMPI_CC gcc-16
# setenv OMPI_FC gfortran-16
#
# also of interest:
# setenv  OMPI_MCA_btl self,tcp    # workaround https://github.com/open-mpi/ompi/issues/6518
#

# openmp is supported by GNU gcc, but not supported by Apple LLVM gcc
SET (ENABLE_OPENMP FALSE CACHE FILEPATH "")

SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicxx CACHE FILEPATH "")

SET (ADD_Fortran_FLAGS "-fbacktrace" CACHE STRING "")
#SET (ADD_Fortran_FLAGS "-fbacktrace -fbounds-check -O0" CACHE STRING "")

SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
EXECUTE_PROCESS(COMMAND nf-config --prefix
  RESULT_VARIABLE NFCONFIG_RESULT
  OUTPUT_VARIABLE NFCONFIG_OUTPUT
  ERROR_VARIABLE  NFCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
IF(NOT NFCONFIG_RESULT EQUAL 0)
  MESSAGE(FATAL_ERROR "nf-config failed: ${NFCONFIG_ERROR}")
ENDIF()
SET (NetCDF_Fortran_PATH "${NFCONFIG_OUTPUT}" CACHE STRING "")

EXECUTE_PROCESS(COMMAND nc-config --prefix
  RESULT_VARIABLE NCCONFIG_RESULT
  OUTPUT_VARIABLE NCCONFIG_OUTPUT
  ERROR_VARIABLE  NCCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
IF(NOT NCCONFIG_RESULT EQUAL 0)
  MESSAGE(FATAL_ERROR "nc-config failed: ${NCCONFIG_ERROR}")
ENDIF()
SET (NetCDF_C_PATH "${NCCONFIG_OUTPUT}" CACHE STRING "")

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

# see ../../README.cmake for options regarding cprnc
# For non-CIME supported systems, assume cprnc is already built:
#SET (CPRNC_DIR /usr/local/bin CACHE FILEPATH "")

SET (USE_MPI_OPTIONS "-oversubscribe" CACHE FILEPATH "")