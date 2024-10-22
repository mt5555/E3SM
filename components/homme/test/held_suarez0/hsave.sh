#!/bin/bash
#
#  Postprocess native grid HS data
#  run zonal mean plots
#  make omega500 map
#
#SBATCH --job-name ncclimo
#SBATCH -p acme-small  # anvil
#XXXSBATCH -p acme-centos6      #chyrsalis
#SBATCH --account=condo
#SBATCH -N 1
#SBATCH --time=12:00:00


if [[ ! -f hszonal.nc ]]; then
  python ~/codes/acme2/components/homme/test/held_suarez0/hsave.py    # time average
  remap pyave.nc                                                      # map to latlon
  ncwa -O  -a lon pyave.latlon.nc hszonal.nc                          # zonal means
fi
contour.py -i pyave.nc -y ngl -o 1 -c -.15,.15 -r 300x600 -p 500 -m andes omega
ncl ~/codes/acme2/components/homme/test/held_suarez0/hsave2.ncl








exit 0
# very slow:  2h for ne30L26 1000 snapshots, 10G each
time ncap2 -v -O -s 'u2=u^2' $orig hsu2.nc  &
time ncap2 -v -O -s 'v2=v^2' $orig hsv2.nc  &
time ncap2 -v -O -s 'T2=T^2' $orig hsT2.nc  &
wait

# 2min for ne30L26 1000 snapshots
time ncra -O -d time,0,2000  $orig hsave.nc  &
time ncra -O -d time,0,2000  hsu2.nc hsaveu2.nc         &
time ncra -O -d time,0,2000  hsv2.nc hsavev2.nc         &
time ncra -O -d time,0,2000  hsT2.nc hsaveT2.nc         &
wait

time ncks -A -v u2 hsaveu2.nc hsave.nc
time ncks -A -v v2 hsavev2.nc hsave.nc
time ncks -A -v T2 hsaveT2.nc hsave.nc


