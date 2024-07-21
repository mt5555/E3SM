#!/bin/bash 
#SBATCH --job-name swt
#SBATCH -p acme-medium
#SBATCH --account=condo
#SBATCH -N 23
#SBATCH --time=2:00:00

set -e
HOMME=`pwd`/../..
wdir=~/scratch1/preqx
MACH=$HOMME/cmake/machineFiles/anvil.cmake
NCPU=0


#
#  problem setup
#
tstep=150
#tstep=120
#tstep=30 
#tstep=10

#tstep=0.125
#nu=1.2e-6  #3.5e-8 for hv_scaling=3.2, 1.2e-6 for hv_scaling=4.0
nu=3.4e-8
hvscaling=3
test_case=swtc2
name=${test_case}-tensor

#mesh='none'
#mesh=~/scratch1/mapping/grids/mountain_10_x2.g
#mesh=~/scratch1/mapping/grids/mountain_10_x8.g
mesh=~/scratch1/mapping/grids/TEMPEST_ne30.g
#mesh=~/scratch1/mapping/grids/TEMPEST_ne480.g  (needs -c 2)
#mesh=~/scratch1/mapping/grids/TEMPEST_ne480.g  (needs -c 2)
#mesh=~/scratch1/mapping/grids/CAne32x128_Altamont100m_v2.g

#
# run/build directoryies
#
input=$HOMME/test/sw_conservative
builddir=$wdir/bldsw
rundir=$wdir/$name
mkdir -p $builddir
mkdir -p $rundir


cd $builddir
if ! ([ -e  CMakeCache.txt ]) then
  rm -rf CMakeFiles CMakeCache.txt
  cmake -C $MACH -DSWEQX_PLEV=1  -DSWEQX_NP=4 $HOMME
  exit
fi


#make clean 
make -j4 sweqx
exe=$builddir/src/sweqx/sweqx



mkdir -p $rundir/movies
cd $rundir

let sfreq=300     # every 300 secends
sfreq=`echo "$sfreq / $tstep" | bc`   # convert to timesteps

sed s/tstep.\*/"tstep = $tstep"/  $input/swtc2-tensor-hv.nl |\
sed s/hypervis_scaling.\*/"hypervis_scaling = $hvscaling"/  |\
sed s:mesh_file.\*:mesh_file="'$mesh'":  |\
sed s/nu=.\*/"nu= $nu"/  |\
sed s/nu_s=.\*/"nu_s= $nu"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  \
    > input.nl



if [ -n "$SLURM_NNODES" ]; then
   mpirun="srun -K -c 1 -N $SLURM_NNODES"
else
   mpirun="srun -K -c 1"
fi
echo $mpirun
$mpirun $exe  < input.nl | tee  sweq.out

mv -f sweq.mass $name.mass
mv -f sweq.out $name.out
mv -f movies/swtc21.nc movies/$name.nc



