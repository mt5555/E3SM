#!/bin/bash
#
# Extract the swtc1 transport operator L for the velocity DSS(curl psi) with
#   psi = psi_solidbody + eps*a*u0*Y_l^m,   eps = r / sqrt(1.5*l*(l+1)),
# i.e. r = rms|u_perturbation| / rms|u_solidbody| (before the rescaling of the
# total velocity to max|u| = u0).  Wraps swtc1-opmat.sh; analyze with
# python/opmat_flows.py.
#
# usage:  swtc1-opmat-ylm.sh NE subcase l m r [forms] [nranks]
#   forms   theta_advect_form values to extract (default "0 2")
#   nranks  MPI ranks (default 4)
#
# output: $outdir/ne${NE}-sc${subcase}-f${form}-curlpsi-ylm${l}_${m}_${eps}.*
#
NE=${1:?NE}
subcase=${2:?subcase}
l=${3:?l}
m=${4:?m}
r=${5:?r}
forms=${6:-"0 2"}
nranks=${7:-4}

here=$(cd $(dirname $0) && pwd)
eps=$(awk -v r=$r -v l=$l 'BEGIN{printf "%.6g", r/sqrt(1.5*l*(l+1))}')
status=0
for f in $forms; do
    SWEQX_TC1_PSI=ylm SWEQX_TC1_PSI_PARAMS=$l,$m,$eps \
        $here/swtc1-opmat.sh $NE $subcase $f curlpsi L $nranks | grep -E "max\|DSS div u\|" | head -1
    status=$(( status + ${PIPESTATUS[0]} ))
done
exit $status
