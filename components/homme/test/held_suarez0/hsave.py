#!/usr/bin/env python3
import xarray as xr
import numpy as np
import time
#
# read in a held_suarez01.nc file of snapshots and compute all needed time averages
# not sure why, but much faster to do it this way than to use ncap2 and ncra
#
# input:  held_suarez01.nc
# ouput:  pyave.nc
#
# to plot zonal means:
#  remap pyave.nc                             # remap to native grid
#  ncwa -a lon pyave.latlon.nc hszonal.nc     # compute zonal averages
#  ncl hsave2.ncl                             # make plots
#
# to plot omega500:
#  contour.py -i pyave.nc -y ngl -o 1 -c -.15,.15 -r 300x600 -p 500 -m andes omega
#

# Define input and output file names
input_file = 'held_suarez01.nc'  # Replace with your input NetCDF file name
#input_file = 'hstest.nc'  # Replace with your input NetCDF file name
output_file = 'pyave.nc'  # Replace with your desired output NetCDF file name

print(f"Opening {input_file}")
time1=time.perf_counter()
ds = xr.open_dataset(input_file)
time2=time.perf_counter()
print(f"  time={time2-time1:.4f} seconds")

u = ds['u']
v = ds['v']
T = ds['T']
omega = ds['omega']
ps = ds['ps']
u2 = np.zeros((u[0,]).shape).squeeze()
v2 = np.zeros((v[0,]).shape).squeeze()
T2 = np.zeros((T[0,]).shape).squeeze()
uave = np.zeros((u[0,]).shape).squeeze()
vave = np.zeros((v[0,]).shape).squeeze()
Tave = np.zeros((T[0,]).shape).squeeze()
omegaave = np.zeros((omega[0,]).shape).squeeze()
psave = np.zeros((ps[0,]).shape).squeeze()
# compute u^2 inside averaging   loop, 30% faster than "u2=u*u", less memory:
print("computing u^2, v^2 and T^2 averages")
time3=time.perf_counter()
nt=u.shape[0]
for i in range(nt):
    u2=u2+u[i,]*u[i,]/nt
    v2=v2+v[i,]*v[i,]/nt
    T2=T2+T[i,]*T[i,]/nt
    uave=uave+u[i,]/nt    # for these variables, this is faster than using u.mean(dim='time')
    vave=vave+v[i,]/nt
    Tave=Tave+T[i,]/nt
    omegaave=omegaave+omega[i,]/nt
    psave=psave+ps[i,]/nt
    if (i % 100 == 0 ): print(f"i={i}/{nt}")
time4=time.perf_counter()
print(f"  time={time4-time3:.4f} seconds")


print(f"writing output")
output_ds = xr.Dataset({'u2': u2,'v2': v2,'T2': T2,
          'u':uave,'v':vave,'T':Tave,'omega':omegaave,'ps':psave,
          'hyam':ds['hyam'],'hybm':ds['hybm'],
          'hyai':ds['hyai'],'hybi':ds['hybi'],
})
output_ds.to_netcdf(output_file, encoding={var: {"_FillValue": None} for var in output_ds.variables})
#output_ds.to_netcdf(output_file)

