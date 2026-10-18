# CMake initial cache file for Ubuntu (tested with 24.04)
# Also works in Claude Code cloud sessions (e.g. with a Claude Pro account),
# which run on an Ubuntu 24.04 VM with 4 cores.
#
# stock gcc/gfortran, OpenMPI, NetCDF and LAPACK from apt:
#
#   apt-get install gfortran libopenmpi-dev openmpi-bin \
#                   libnetcdff-dev libnetcdf-dev liblapack-dev libblas-dev
#
# If running as root (e.g. in a container), OpenMPI's mpirun needs:
#   export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicxx CACHE FILEPATH "")

SET (ADD_Fortran_FLAGS "-fbacktrace" CACHE STRING "")

SET (ENABLE_OPENMP FALSE CACHE BOOL "")

SET (WITH_PNETCDF FALSE CACHE BOOL "")
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

# Ubuntu puts libraries in /usr/lib/<multiarch>, but SCORPIO's FindNetCDF
# only searches <prefix>/lib.  Have find_library also try lib/<multiarch>.
EXECUTE_PROCESS(COMMAND gcc -print-multiarch
  OUTPUT_VARIABLE MULTIARCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
IF(MULTIARCH)
  SET (CMAKE_FIND_LIBRARY_CUSTOM_LIB_SUFFIX "/${MULTIARCH}" CACHE STRING "")
ENDIF()

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET (USE_NUM_PROCS 4 CACHE STRING "")
SET (USE_MPIEXEC "mpirun" CACHE STRING "")
SET (USE_MPI_OPTIONS "--allow-run-as-root --oversubscribe" CACHE STRING "")
