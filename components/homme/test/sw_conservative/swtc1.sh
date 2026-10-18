#!/bin/bash
#
#SBATCH --job-name swtc1
#SBATCH -N 1                                                                                                                             
#SBATCH -c 2                                                                                                                             
#SBATCH -C cpu                                                                                                                           
#SBATCH --time=0:30:00                                                                                                                   

bld=~/scratch1/homme/swtc1
HOMME=~/codes/acme2/components/homme
MACH=$HOMME/cmake/machineFiles/darwin.cmake
input=$HOMME/test/sw_conservative
output=$bld                # output will be in $output/subcase#

test_case=swtc1

EXEC=$bld/src/sweqx/sweqx

mkdir -p $bld
cd $bld
#rm -Rf CMakeCache.txt CMakeFiles   # remove this file to force re-configure
if [[ ! -f CMakeCache.txt  ]]; then
    echo "running CMAKE to configure the model"
    cmake -C $MACH -DSWEQX_PLEV=4 -DSWEQX_NP=4  $HOMME
    [ $? -ne 0 ] && exit 1
    make -j4 clean
fi
make -j4 sweqx
[ $? -ne 0 ] && exit 1
[ ! -x $EXEC  ] && exit 1

# subcase 3, nu=0, 500 days it will go unstable
subcase=4          # 1,2,3,4 for alpha=0,pi/2,pi/2-0.05,pi/4
adv_form=2        # theta_advect_form: 0 = conservation form, 2 = split form

ndays=$((12*30))   # 12 days per revolution
NE=30
#nu=1e15  #  requires dt<=360
nu=2.0e14  # RK4, tstep=1800  2.25e14 RUNS.  2.5e14 starts to blow up 38 days

# HOMME's hacky RK2-4 method
#integration=explicit
#LFTfreq=4       # RK2-4  used in HOMME 3D subcyclign
#rk_stage=0
#tstep=120


# the RK SSP methods:
#integration=runge_kutta
#LFTfreq=0
#rk_stage=3
#tstep=300    # ok: 1200  bad: 1800


# classic RK4
integration=runge_kutta
LFTfreq=0
rk_stage=-4
tstep=1800    # 2160 is a good value.
              # with nu=0, dt=2330 is stable for 3 rev, dt=2340 unstable
              # with nu=1e14, dt=2580 is stable for 3 rev, dt=2590 unstable
#tstep=300

smooth=0                         # only set if using classic LF
limiter=0                        # no limiter
name=${test_case}-sc${subcase}-NE${NE}-t${tstep}-limiter$limiter
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
> swtc1.nl


# get number of performance cores (Apple Silicon), else number of CPUs (Linux)
ncpu=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || getconf _NPROCESSORS_ONLN)
date
if command -v srun >/dev/null 2>&1; then
    srun -c 2 -N ${SLURM_NNODES:-1} $EXEC < swtc1.nl
else
    echo using mpirun -np $ncpu
    mpirun -np $ncpu  $EXEC < swtc1.nl
fi
date

python3 $input/python/plot_swtc1_errors.py .
