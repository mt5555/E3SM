#!/bin/bash
#
# Reproduce the runs behind the swtc1-DG and CG vs DG galleries, with the run
# directory names that conv_dg.py, plots_dg.py and plots_cvd.py expect.
# All runs: swtc1, classical RK4, split form unless noted.  See README.md.
#
# usage: gallery_runs.sh [SECTION ...]      (default: all sections, in this order)
#   forms   swtc1-DG gallery: NE30, dt=900, 30 revolutions, DG central/upwind x split/conservation form
#   conv    1 revolution, dt = 22500/NE: CG nu=0, CG near-optimal nu, DG central, DG upwind
#           (sub_case 4, NE8-120), the same for sub_case 1 (NE8-60), and conservation form DG (NE8-60)
#   cvd     CG vs DG gallery: NE30, dt=750, 3 revolutions, the four methods, with movies
#   nuscan  CG hyperviscosity scan, NE30 and NE60, 1 revolution
#   stab    RK4 stability limit of DG central and upwind by bisection, NE30, 1 revolution
#
# The NE120 runs take about 1-1.5 hours each on 4 cores; everything else is minutes.
# Builds sweqx on the first run (see run.sh for bld, output, NP, MACH, ...).

RUN=$(cd "$(dirname "$0")" && pwd)/run.sh
output=${output:-~/scratch1/homme/runs}; export output
sections=${*:-forms conv cvd nuscan stab}
nuc=1e13     # CG hyperviscosity nu = nuc*(30/NE)^4; within 5% of optimal at NE30 and NE60
nu_of() { python3 -c "print('%e' % ($nuc*(30/$1)**4))"; }

for s in $sections; do case $s in
forms)
    for f in 2 0; do for dg in 0 1; do
        MOVIES=1 $RUN gal_dg${dg}_f$f 30 4 900 360 $dg $f
    done; done ;;
conv)
    for ne in 8 15 30 60 120; do
        for dg in -1 0 1; do $RUN conv_s4_ne${ne}_dg$dg $ne 4 - 12 $dg 2; done
        NU=$(nu_of $ne) $RUN conv_s4_ne${ne}_cgnu $ne 4 - 12 -1 2
    done
    for ne in 8 15 30 60; do
        for dg in -1 0 1; do $RUN conv_s1_ne${ne}_dg$dg $ne 1 - 12 $dg 2; done
        for dg in 0 1; do $RUN convf0_s4_ne${ne}_dg$dg $ne 4 - 12 $dg 0; done
    done ;;
cvd)
    for dg in -1 0 1; do MOVIES=1 $RUN cvd_var_dg$dg 30 4 750 36 $dg 2 24; done
    NU=$nuc MOVIES=1 $RUN cvd_var_cgnu 30 4 750 36 -1 2 24 ;;
nuscan)
    for ne in 30 60; do for c in 3e12 6e12 2e13 4e13; do   # c = 1e13 is conv_s4_ne*_cgnu
        NU=$(python3 -c "print('%e' % ($c*(30/$ne)**4))") $RUN nuscan_ne${ne}_c$c $ne 4 - 12 -1 2
    done; done ;;
stab)
    # stable = finite Gaussian l2 error below 1e-2 after 12 days
    for dg in 1 0; do lo=600; hi=1800
        while [ $((hi - lo)) -gt 20 ]; do
            mid=$(( (lo + hi) / 2 ))
            $RUN stab_dg${dg}_$mid 30 4 $mid 12 $dg 2 > /dev/null
            l2=$(tail -1 $output/stab_dg${dg}_$mid/swtc1.gaus.l2.errors 2>/dev/null | awk '{print $2}')
            if python3 -c "import sys; v=float('$l2' or 'nan'); sys.exit(0 if v==v and v<1e-2 else 1)"; then lo=$mid; else hi=$mid; fi
            echo "stab dg_flux=$dg dt=$mid l2=$l2 -> stable range [$lo, $hi)"
        done
    done ;;
*) echo "unknown section $s"; exit 1 ;;
esac; done
