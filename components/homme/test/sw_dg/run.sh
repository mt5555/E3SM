#!/bin/bash
#
# Run one sweqx swtc1 (or swirl) transport case with classical RK4, CG or DG.
# Builds sweqx first unless an executable is given.  See README.md.
#
# usage: run.sh [NAME NE SUBCASE TSTEP NDAYS DGFLUX FORM [STATEFREQ] [TESTCASE] [EXEC]]
#   all arguments are optional; "-" (or omitting trailing arguments) means the default
#   NAME      run directory name, created under $output     (default: from the parameters)
#   NE        elements per cube-face edge                    (default 30)
#   SUBCASE   1,2,3,4 for alpha = 0, pi/2, pi/2-0.05, pi/4   (default 4)
#   TSTEP     seconds                                        (default 22500/NE, 750 at NE30)
#   NDAYS     run length in days, 12 days per revolution     (default 36)
#   DGFLUX    -1 = CG (DSS), 0 = DG central, 1 = DG upwind   (default -1)
#             "none" leaves dg_flux out of the namelist (for executables built before the DG option)
#   FORM      theta_advect_form: 0 = conservation form, 2 = split form   (default 2)
#   STATEFREQ steps between error/diagnostic output          (default: about once a day, ending at the last step)
#   TESTCASE  swtc1 or swirl                                 (default swtc1)
#   EXEC      run this sweqx executable instead of building one
#
# environment:
#   NU        hyperviscosity nu = nu_s (default 0; DG requires 0).  For CG, near-optimal is
#             about 1e13*(30/NE)^4; the RK4 limit is roughly nu*dt < 4.5e17*(30/NE)^4.
#   NP        MPI ranks (default: performance cores on Apple Silicon, else all CPUs)
#   MOVIES    1 = keep the template's movie output (daily for 36 days + every 12 days), 0 = none (default 0)
#   PLOT      1 = run sw_conservative/python/plot_swtc1_errors.py at the end (swtc1 only) (default 0)
#   MPIRUN_FLAGS  extra mpirun flags, e.g. "--bind-to none" to run several cases at once with Open MPI
#   bld       build directory (default ~/scratch1/homme/swtc1)
#   output    parent of the run directories (default ~/scratch1/homme/runs)
#   HOMME     components/homme (default: found from this script's location)
#   MACH      cmake machine file (default $HOMME/cmake/machineFiles/darwin.cmake)
#
# The namelist template is sw_conservative/swtc1.nl.
#
# examples:
#   ./run.sh                                     # CG, nu=0, NE30, 3 revolutions
#   ./run.sh - 30 4 - 36 1                       # same, DG upwind
#   NU=1e13 ./run.sh cg_nu_ne30                  # CG with near-optimal hyperviscosity
#   for ne in 15 30 60; do ./run.sh conv_ne$ne $ne 4 - 12 1; done   # convergence, 1 revolution

arg() { if [ -z "$1" ] || [ "$1" = - ]; then echo "$2"; else echo "$1"; fi; }
NE=$(arg "$2" 30)
SC=$(arg "$3" 4)
DT=$(arg "$4" $(python3 -c "print(22500/$NE)"))
ND=$(arg "$5" 36)
DG=$(arg "$6" -1)
FORM=$(arg "$7" 2)
TC=$(arg "$9" swtc1)
NU=${NU:-0}
nsteps=$(python3 -c "print(round($ND*86400/$DT))")
# default: the divisor of nsteps closest to one day, so the last output is at the end of the run
# (or simply one day, if no divisor is within a factor of 2 of it)
SF=$(arg "$8" $(python3 -c "
n=$nsteps; d=max(1.0, 86400/$DT)
k=[k for k in range(1,n+1) if n%k==0 and d/2<=k<=2*d]
print(min(k, key=lambda k: abs(k-d)) if k else min(n, round(d)))"))
name=${TC}-sc${SC}-ne${NE}-dt${DT}-form${FORM}
[ "$DG" != none ] && [ "$DG" -ge 0 ] && name=$name-dg$DG
awk -v n="$NU" 'BEGIN{exit !(n+0!=0)}' && name=$name-nu$NU   # only when nu is nonzero
NAME=$(arg "$1" $name)

HOMME=${HOMME:-$(cd "$(dirname "$0")/../.." && pwd)}
bld=${bld:-~/scratch1/homme/swtc1}
output=${output:-~/scratch1/homme/runs}
MACH=${MACH:-$HOMME/cmake/machineFiles/darwin.cmake}
input=$HOMME/test/sw_conservative

EXEC=${10:-$bld/src/sweqx/sweqx}
if [ -z "${10}" ]; then
    mkdir -p $bld
    cd $bld
    #rm -Rf CMakeCache.txt CMakeFiles   # remove these to force re-configure
    if [[ ! -f CMakeCache.txt ]]; then
        echo "running CMAKE to configure the model"
        cmake -C $MACH -DSWEQX_PLEV=4 -DSWEQX_NP=4 $HOMME
        [ $? -ne 0 ] && exit 1
        make -j4 clean
    fi
    make -j4 sweqx
    [ $? -ne 0 ] && exit 1
fi
[ ! -x $EXEC ] && echo "no executable $EXEC" && exit 1

rundir=$output/$NAME
mkdir -p $rundir/movies
cd $rundir
rm -f swtc1.*.errors sweq.* movies/*.nc     # no stale output from an earlier run

DGSTR=", dg_flux = $DG"; [ "$DG" = none ] && DGSTR=""
MOVIEOPT='-e s/^output_frequency.*/output_frequency=0,0/'; [ "${MOVIES:-0}" = 1 ] && MOVIEOPT=""
sed -e "s/^ne .*/ne = $NE/" \
    -e "s/^tstep.*/tstep = $DT/" \
    -e "s/^smooth.*/smooth = 0/" \
    -e "s/^limiter_option.*/limiter_option = 0/" \
    -e "s/^test_case.*/test_case = '$TC', sub_case = $SC, theta_advect_form = $FORM$DGSTR/" \
    -e "s/^integration.*/integration = \"runge_kutta\"/" \
    -e "s/^rk_stage_user.*/rk_stage_user = -4/" \
    -e "s/^LFTfreq.*/LFTfreq = 0/" \
    -e "s/^ndays.*/nmax = $nsteps/" \
    -e "s/^nu=.*/nu= $NU/" \
    -e "s/^nu_s=.*/nu_s= $NU/" \
    -e "s/^statefreq.*/statefreq = $SF/" \
    $MOVIEOPT $input/swtc1.nl > swtc1.nl

# number of ranks: performance cores (Apple Silicon), else number of CPUs (Linux)
NP=${NP:-$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || getconf _NPROCESSORS_ONLN)}
echo "$rundir: NE=$NE sub_case=$SC dt=$DT days=$ND ($nsteps steps) dg_flux=$DG form=$FORM nu=$NU np=$NP"
if command -v srun >/dev/null 2>&1; then
    srun -c 2 -N ${SLURM_NNODES:-1} $EXEC < swtc1.nl > out.txt 2>&1
else
    mpirun $MPIRUN_FLAGS -np $NP $EXEC < swtc1.nl > out.txt 2>&1
fi
status=$?
echo "$NAME exit=$status $([ "$TC" = swtc1 ] && grep 'k= 2' out.txt | tail -1)"
[ $status -ne 0 ] && tail -5 out.txt && exit $status

[ "$TC" = swtc1 ] && [ "${PLOT:-0}" = 1 ] && python3 $input/python/plot_swtc1_errors.py .
exit 0
