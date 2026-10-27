#!/bin/bash
#
# Extract the matrices of the swtc1 transport operator L (dp/dt = L p) and of
# the hyperviscosity operator B (del^4, as applied in advance_hypervis), for
# eigenvalue analysis with python/opmat_eigs.py.  See src/sweqx/opmat_mod.F90.
#
# usage:  swtc1-opmat.sh NE subcase adv_form velocity [ops] [nranks]
#   subcase   1,2,3,4 for alpha=0,pi/2,pi/2-0.05,pi/4
#   adv_form  theta_advect_form: 0 = conservation form, 2 = split form
#   velocity  analytic (solid body rotation at the GLL nodes) or
#             curlpsi  (DSS(curl(psi)), discretely divergence free)
#   ops       L, B or LB (default LB).  B depends only on NE.
#   nranks    MPI ranks (default 1).  Any number works; files are written per rank.
#
# output: $outdir/ne${NE}-sc${subcase}-f${adv_form}-${velocity}.*
#
# with the default coloring the number of operator applications is independent
# of NE (384 for L, 864 for B with np=4, nlev=4).  SWEQX_OPMAT_COLOR=0 probes one
# column per level per application (n/4 applications), for checking.
#
bld=${bld:-~/scratch/opmat-build}           # built with -DSWEQX_PLEV=4 -DSWEQX_NP=4
HOMME=${HOMME:-~/codes/acme/components/homme}
outdir=${outdir:-~/scratch/opmat}
input=$HOMME/test/sw_conservative
EXEC=$bld/src/sweqx/sweqx

NE=${1:?NE}
subcase=${2:?subcase}
adv_form=${3:?adv_form}
velocity=${4:?velocity}
ops=${5:-LB}
nranks=${6:-1}

[ ! -x $EXEC ] && { echo "missing $EXEC"; exit 1; }
mkdir -p $outdir
outdir=$(cd $outdir && pwd)
name=ne${NE}-sc${subcase}-f${adv_form}-${velocity}

rundir=$outdir/run-$name
mkdir -p $rundir/movies
cd $rundir

# nu=0: the extraction does not use nu or nu_s.  hypervis_order=2 is required for B.
sed s/^ne.\*/"ne = $NE"/  $input/swtc1.nl |\
sed s/tstep.\*/"tstep = 600"/  |\
sed s/limiter_option.\*/"limiter_option = 0"/  |\
sed s/smooth.\*/"smooth = 0"/  |\
sed s/test_case.\*/"test_case = \'swtc1\', sub_case = $subcase, theta_advect_form = $adv_form"/  |\
sed s/rk_stage_user.\*/"rk_stage_user = -4"/  |\
sed "s/^integration.*/integration = \"runge_kutta\"/"  |\
sed s/ndays.\*/"ndays = 1"/  |\
sed s/nu=.\*/"nu= 0"/  |\
sed s/nu_s=.\*/"nu_s= 0"/  |\
sed s/hypervis_order.\*/"hypervis_order = 2"/  |\
sed s/interp_nlon.\*/"interp_nlon = 64"/  |\
sed s/interp_nlat.\*/"interp_nlat = 32"/  \
> swtc1.nl

export SWEQX_TC1_VELOCITY=$velocity
export SWEQX_OPMAT=$outdir/$name
export SWEQX_OPMAT_OPS=$ops
if command -v srun >/dev/null 2>&1; then
    srun -n $nranks $EXEC < swtc1.nl > sweqx.log 2>&1
else
    mpirun -np $nranks $EXEC < swtc1.nl > sweqx.log 2>&1
fi
status=$?
grep -E "SWEQX_OPMAT|DSS\(curl psi\)|prescribed velocity" sweqx.log
exit $status
