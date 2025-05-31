#!/bin/bash
#
#XXSBATCH --account=FY150001
#XXSBATCH -p ec
#SBATCH --account=condo
#SBATCH -p acme-medium
#SBATCH -N 6
#SBATCH --time=2:00:00
#
# Anvil: 6 nodes, 1h for all runs
#
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=16M     #  Cori has 96GB per node. had to lower to 8M on 3K nodes
export MV2_ENABLE_AFFINITY=0
set | grep NODE


# hydrostatic
EXEC=../../../test_execs/theta-l-nlev30/theta-l-nlev30
#EXEC=../../../test_execs/thetada-nlev30/thetada-nlev30
\rm -f input.nl
mkdir restart


function run { 
local tstep=$1
#rsplit1=`echo "print( round((2*1200)/$tstep))"  | python3`
local rsplit1=$2


logfile1=log1-${tstep}.out
\rm -f $logfile1 

# change statefreq to every 6h
sfreq=`echo "print( round((6*3*1200)/$tstep))"  | python3`


# run 15 days, write a restart file
namelist=checkE-run1.nl
sed s/tstep.\*1200/"tstep=$tstep"/  $namelist |
sed s/statefreq.\*/"statefreq=$sfreq"/  |
sed s/rsplit.\*/"rsplit=$rsplit1"/ > input.nl

echo "tstep=$tstep rsplit=$rsplit1,  run to 15 days and write restart file..."
\rm -f restart/*
mpirun $EXEC < input.nl  > $logfile1

grep SHA $logfile1




echo "diagnostics time=6d dt=$tstep"           | tee -a day6.log
grep "KE,d/dt" $logfile1 | tail -37 | head -1  | tee -a day6.log
grep "IE,d/dt" $logfile1 | tail -37 | head -1  | tee -a day6.log
grep "PE,d/dt" $logfile1 | tail -37 | head -1  | tee -a day6.log
grep " E,d/dt" $logfile1 | tail -37 | head -1  | tee -a day6.log
grep "E-E0" $logfile1    | tail -37 | head -1  | tee -a day6.log

echo "diagnostics time=7d dt=$tstep"           | tee -a day7.log
grep "KE,d/dt" $logfile1 | tail -33 | head -1  | tee -a day7.log
grep "IE,d/dt" $logfile1 | tail -33 | head -1  | tee -a day7.log
grep "PE,d/dt" $logfile1 | tail -33 | head -1  | tee -a day7.log
grep " E,d/dt" $logfile1 | tail -33 | head -1  | tee -a day7.log
grep "E-E0" $logfile1    | tail -33 | head -1  | tee -a day7.log

echo "diagnostics time=9d dt=$tstep"           | tee -a day9.log
grep "KE,d/dt" $logfile1 | tail -25 | head -1  | tee -a day9.log
grep "IE,d/dt" $logfile1 | tail -25 | head -1  | tee -a day9.log
grep "PE,d/dt" $logfile1 | tail -25 | head -1  | tee -a day9.log
grep " E,d/dt" $logfile1 | tail -25 | head -1  | tee -a day9.log
grep "E-E0" $logfile1    | tail -25 | head -1  | tee -a day9.log

echo "CAAR diagnostics time=15d dt=$tstep"     | tee -a day15.log
grep "KE,d/dt" $logfile1 | tail -1             | tee -a day15.log
grep "IE,d/dt" $logfile1 | tail -1             | tee -a day15.log
grep "PE,d/dt" $logfile1 | tail -1             | tee -a day15.log
grep " E,d/dt" $logfile1 | tail -1             | tee -a day15.log
grep "E-E0" $logfile1    | tail -1             | tee -a day15.log


# debug:
#  grep "E-E0" $logfile2 | tail -1                   time=15d
#                        ! tail -25 | head -1        time=9d
#                        ! tail -37 | head -1        time=6d
#                        ! tail -41 | head -1        time=5d
# 

}

rm -f day*.log

#rsplit=0   # Euerlian
#rsplit=2   # Lagrangian
rsplit=1

date
#run 1200  $rsplit
run 600   $rsplit
run 300   $rsplit
run 150   $rsplit
run 75    $rsplit
run 37.5  $rsplit
run 18.75 $rsplit
date

