#!/bin/bash
#
#SBATCH --job-name swirl
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -C cpu
#SBATCH --time=0:30:00
#
# swirl (Nair & Lauritzen 2010 deformational flow) with the split form advection
# levels: 1 cosine bells, 2 Gaussian bells, 3 slotted cylinders, 4 density (=1),
#         5-6 terminator toy chemistry (sub_case 4 only)
# levels 1-3 are mixing ratios (divided by level 4); the flow reverses so the
# exact solution at t=12 days is the initial condition

bld=~/scratch1/homme/swirl
HOMME=~/codes/acme/components/homme
MACH=$HOMME/cmake/machineFiles/darwin.cmake
input=$HOMME/test/sw_conservative
output=$bld                # output will be in $output/subcase#

test_case=swirl

EXEC=$bld/src/sweqx/sweqx

mkdir -p $bld
cd $bld
#rm -Rf CMakeCache.txt CMakeFiles   # remove this file to force re-configure
if [[ ! -f CMakeCache.txt  ]]; then
    echo "running CMAKE to configure the model"
    cmake -C $MACH -DSWEQX_PLEV=6 -DSWEQX_NP=4  $HOMME
    [ $? -ne 0 ] && exit 1
    make -j4 clean
fi
make -j4 sweqx
[ $? -ne 0 ] && exit 1
[ ! -x $EXEC  ] && exit 1

subcase=4          # 1,2: deformational flow, 3,4: 1,2 plus solid body rotation
adv_form=2         # theta_advect_form: 0 = conservation form, 2 = split form

ndays=12           # one period
NE=30
nu=2.0e14

# classic RK4
integration=runge_kutta
LFTfreq=0
rk_stage=-4
tstep=600     # 1800 blows up in the first day (max wind 90 m/s vs 39 m/s for swtc1)

smooth=0                         # only set if using classic LF
limiter=0                        # no limiter (limiters are only used by the RK SSP methods)
sfreq=24  # error output ever 24h
sfreq=`echo "$sfreq * 3600 / $tstep" | bc`
nu_s=$nu




rundir=$output/subcase$subcase
[ $adv_form -eq 2 ] && rundir=$rundir-splitform
[ $rk_stage -gt 0 ] && rundir=$rundir-ssp2s${rk_stage}   # 2nd-order SSP, rk_stage stages
[ $rk_stage -eq -4 ] && rundir=$rundir-rk4
awk -v n="$nu" 'BEGIN{exit !(n+0!=0)}' && rundir=$rundir-nu$nu   # only when nu is nonzero
mkdir -p $rundir
cd $rundir
mkdir -p movies

sed s/^ne.\*/"ne = $NE"/  $input/swtc1.nl |\
sed s/tstep.\*/"tstep = $tstep"/  |\
sed s/limiter_option.\*/"limiter_option = $limiter"/  |\
sed s/smooth.\*/"smooth = $smooth"/  |\
sed s/test_case.\*/"test_case = \'$test_case\', sub_case = $subcase, theta_advect_form = $adv_form"/  |\
sed s/rk_stage_user.\*/"rk_stage_user = $rk_stage"/  |\
sed "s/^integration.*/integration = \"$integration\"/"  |\
sed s/LFTfreq.\*/"LFTfreq = $LFTfreq"/  |\
sed s/ndays.\*/"ndays = $ndays"/  |\
sed s/nu=.\*/"nu= $nu"/  |\
sed s/nu_s=.\*/"nu_s= $nu_s"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  \
> swirl.nl


# get number of performance cores (Apple Silicon), else number of CPUs (Linux)
ncpu=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || getconf _NPROCESSORS_ONLN)
date
if command -v srun >/dev/null 2>&1; then
    srun -c 2 -N ${SLURM_NNODES:-1} $EXEC < swirl.nl | tee sweq.out
else
    echo using mpirun -np $ncpu
    mpirun -np $ncpu  $EXEC < swirl.nl | tee sweq.out
fi
date

# errors at the end of the run, all levels (the .errors files only have level 1)
grep 'l1,l2,linf=' sweq.out | tail -6
