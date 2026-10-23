#!/bin/bash                                                                                                                              
#                                                                                                                                        
#SBATCH --job-name d21-theta                                                                                                             
#SBATCH -N 8                                                                                                                             
#SBATCH -c 2                                                                                                                             
#SBATCH -C cpu                                                                                                                           
#SBATCH --time=0:45:00
#SBATCH -q regular
# regular or debug (default)                                                                                                             
#
# Output: moviesVL-L<nlev>-nu_top<nu_top>-<sponge_start>-rf<dcmip_rayleigh_depth>[-n<nmax>]
#
# option:
#     testcase=dcmip,m1,m2,m3
#     dcmip_rayleigh_depth=-1 (default)  use test case Rayleigh friction sponge
#     dcmip_rayleigh_depth=0             no Rayleigh friction
#     dcmip_rayleigh_depth=<value>       Rayleigh friction sponge depth (m) below model top
#     nu_top=<value>        EAM/EAMxx sponge layer coefficient.  default 0 (disabled)
#     sponge_start=<value>  reference pressure (hPa) where the top sponge starts (default 2)
#     nlev:                 not set, use default test analytic levels (dcmip=60L, m1=40L, m2=60L)
#     nlev=144, 152, 160, 176, 184, 192, 256    use levels from vcoord file
#     nmax:                 not set, use template value (2h for all tests)
#     nmax=<steps>          number of timesteps (output dirs get a -n<steps> suffix)
#
# Examples:
#     sbatch --export=ALL,nu_top=4000,sponge_start=100,dcmip_rayleigh_depth=0 mtest1.sh
#     sbatch --export=ALL,dcmip_rayleigh_depth=20000 mtest1.sh
#



export OMP_NUM_THREADS=1
SCRATCH=/pscratch/sd/t/taylorm/homme
HOMME=/global/u2/t/taylorm/codes/acme/components/homme
MACH=$HOMME/cmake/machineFiles/pm-cpu.cmake

testcase=${testcase:-m2}
# vlevels=analytic (vanalytic=1, hardcoded levels) unless the user
# explicitly sets nlev, in which case we switch to file-based vcoords
# read from vcoord-<nlev>{m,i}.ascii
if [[ -z ${nlev+x} ]]; then
    vlevels=analytic
else
    vlevels=file
fi
case $testcase in
    dcmip) TEMPLATE=$HOMME/test/mtest/dcmip2012-nh-template.nl; TESTDIR=dcmip2012-2.1; OUTFILE=dcmip2012_test2_11.nc; nlev=${nlev:-60}; nu_top=${nu_top:-0} ;;
    m1)    TEMPLATE=$HOMME/test/mtest/m1-ne30-template.nl;      TESTDIR=m1;            OUTFILE=mtest11.nc;            nlev=${nlev:-40}; nu_top=${nu_top:-0} ;;
    m2)    TEMPLATE=$HOMME/test/mtest/m2-ne60-template.nl;      TESTDIR=m2;            OUTFILE=mtest21.nc;            nlev=${nlev:-60}; nu_top=${nu_top:-0} ;;
    *) echo "ERROR: unknown testcase '$testcase' (expected dcmip, m1, or m2)"; exit 1 ;;
esac
dcmip_rayleigh_depth=${dcmip_rayleigh_depth:--1}
sponge_start=${sponge_start:-2}
qsize=0

if [[ $vlevels == analytic ]]; then
    VANALYTIC=1
    VFILE_MID=unused
    VFILE_INT=unused
else
    VANALYTIC=0
    VFILE_MID=$HOMME/test/mtest/vcoord-${nlev}m.ascii
    VFILE_INT=$HOMME/test/mtest/vcoord-${nlev}i.ascii
fi
for f in $TEMPLATE; do
    if [[ ! -f $f ]]; then
        echo "ERROR: required file not found: $f"
        exit 1
    fi
done
if [[ $vlevels == file ]]; then
    for f in $VFILE_MID $VFILE_INT; do
        if [[ ! -f $f ]]; then
            echo "ERROR: required file not found: $f"
            exit 1
        fi
    done
fi

# build case for nlev, lat/lon output
bld=$SCRATCH/bld-$nlev
EXEC=$bld/src/theta-l/theta-l
mkdir -p $bld
cd $bld
#rm -Rf CMakeCache.txt CMakeFiles   # remove this file to force re-configure
if [[ ! -f CMakeCache.txt  ]]; then
    echo "running CMAKE to configure the model"
    cmake -C $MACH -DQSIZE_D=$qsize -DPREQX_PLEV=$nlev -DPREQX_NP=4  \
          -DBUILD_HOMME_SWDGX=FALSE                     \
          -DBUILD_HOMME_SWEQX=FALSE                     \
          -DBUILD_HOMME_PRIMDGX=FALSE                   \
          -DPREQX_USE_ENERGY=TRUE  $HOMME
    [ $? -ne 0 ] && exit 1
    make -j4 clean
fi
make -j4 theta-l
[ $? -ne 0 ] && exit 1
[ ! -x $EXEC  ] && exit 1



WDIR=$SCRATCH/$TESTDIR
TAG=L${nlev}-nu_top${nu_top}-${sponge_start}-rf${dcmip_rayleigh_depth}${nmax:+-n$nmax}
MOVIES_DIR=moviesVL-$TAG
mkdir -p $WDIR/$MOVIES_DIR
cd $WDIR
# nonhydrostatic theta: generate input.nl from the shared template,
# substituting in nu_top, the per-nu_top/sponge_start/dcmip_rayleigh_depth movies
# output subdirectory, and (if vlevels=file) the vertical coordinate
# ascii files for this nlev. The namelist filename is tagged with nlev,
# nu_top, sponge_start, and dcmip_rayleigh_depth so that concurrent jobs (differing
# in any of those) running in this same WDIR don't clobber each other's
# input file.
INPUT_NL=input-$TAG.nl
NMAX_SED="s#^\s*nmax\s*=.*#  nmax              = ${nmax}                     ! set by mtest1.sh#"
[[ -z $nmax ]] && NMAX_SED=""
sed -e "s#@NU_TOP@#$nu_top#" \
    -e "$NMAX_SED" \
    -e "s#@SPONGE_START@#$sponge_start#" \
    -e "s#@MOVIES_DIR@#$MOVIES_DIR#" \
    -e "s#@VANALYTIC@#$VANALYTIC#" \
    -e "s#@VFILE_MID@#$VFILE_MID#" \
    -e "s#@VFILE_INT@#$VFILE_INT#" \
    -e "s#@DCMIP_RAYLEIGH_DEPTH@#$dcmip_rayleigh_depth#" \
    $TEMPLATE > $INPUT_NL

echo WDIR = $WDIR
echo NAMELIST = $INPUT_NL
echo OUTPUT = $MOVIES_DIR

date
if [ $SLURM_NNODES  ]; then
    srun -c 2 -N $SLURM_NNODES $EXEC < $INPUT_NL
    NCFILE=$MOVIES_DIR/$OUTFILE
    PLOTXZ=$HOME/codes/nclscript/contour/plotxz.py
    if [[ $testcase == dcmip ]]; then
        python3 $PLOTXZ -i $NCFILE -c -1.2,1.2 Tprime
        python3 $PLOTXZ -i $NCFILE -c -3.5,3.5 w
    elif [[ $testcase == m1 || $testcase == m2 ]]; then
        python3 $PLOTXZ -i $NCFILE -c -3.0,3.0 thetaprime
        python3 $PLOTXZ -i $NCFILE -c -2.0,2.0 w
    else
        echo "unknown testcase '$testcase'; no plotting configured"
    fi
else
    echo SLRUM_NNODES not set, skipping srun
fi
date


