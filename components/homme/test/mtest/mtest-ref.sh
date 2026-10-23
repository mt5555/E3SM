#!/bin/bash
#
#SBATCH --job-name mtest-ref
#SBATCH -N 8
#SBATCH -c 2
#SBATCH -C cpu
#SBATCH --time=0:30:00
#SBATCH -q debug
#
# Run the reference HOMME theta-l mountain wave tests (KSP2015) with the
# default settings from the static namelists m1-ne30.nl, m2-ne60.nl, m3-ne60.nl
#
# Submit from this directory (components/homme/test/mtest):
#     sbatch mtest-ref.sh                          # m2 (default)
#     sbatch --export=ALL,testcase=m1 mtest-ref.sh
#
# Output: $SCRATCH/homme/<testcase>-ref/movies/
#

testcase=${testcase:-m2}
case $testcase in
    m1) NAMELIST=m1-ne30.nl; nlev=40 ;;
    m2) NAMELIST=m2-ne60.nl; nlev=60 ;;
    m3) NAMELIST=m3-ne60.nl; nlev=60 ;;
    *) echo "ERROR: unknown testcase '$testcase' (expected m1, m2, or m3)"; exit 1 ;;
esac

export OMP_NUM_THREADS=1
HOMME=$(cd ${SLURM_SUBMIT_DIR:-.}/../.. && pwd)
MACH=$HOMME/cmake/machineFiles/pm-cpu.cmake
NAMELIST=$HOMME/test/mtest/$NAMELIST
if [[ ! -f $NAMELIST ]]; then
    echo "ERROR: namelist not found: $NAMELIST (submit from components/homme/test/mtest)"
    exit 1
fi

# build
bld=$SCRATCH/homme/bld-$nlev
EXEC=$bld/src/theta-l/theta-l
mkdir -p $bld && cd $bld || exit 1
if [[ ! -f CMakeCache.txt ]]; then
    cmake -C $MACH -DQSIZE_D=0 -DPREQX_PLEV=$nlev -DPREQX_NP=4 \
          -DBUILD_HOMME_SWDGX=FALSE                    \
          -DBUILD_HOMME_SWEQX=FALSE                    \
          -DBUILD_HOMME_PRIMDGX=FALSE                  \
          -DPREQX_USE_ENERGY=TRUE $HOMME || exit 1
fi
make -j4 theta-l || exit 1

# run
WDIR=$SCRATCH/homme/$testcase-ref
mkdir -p $WDIR/movies && cd $WDIR || exit 1
cp $NAMELIST input.nl
date
srun -c 2 -N $SLURM_NNODES $EXEC < input.nl
date
