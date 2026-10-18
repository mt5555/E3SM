#!/usr/bin/env python3
# Animate one swtc1 shape (geop level) from movies/swtc11.nc as a GIF.
# usage: animate_swtc1.py rundir out.gif [level] [symlog]
#   level: 1 cosine bell (default), 2 gaussian, 3 cylinder, 4 slotted cylinder
#   symlog: diverging color map on a symmetric log scale, to show small +/- noise
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import SymLogNorm
from scipy.io import netcdf_file

rundir, out = sys.argv[1], sys.argv[2]
lev = int(sys.argv[3]) if len(sys.argv) > 3 else 1
name = ['cosine bell', 'gaussian', 'cylinder', 'slotted cylinder'][lev - 1]
hmax = 1.0 if lev == 4 else 1000.0
units = '' if lev == 4 else ' m'
symlog = len(sys.argv) > 4 and sys.argv[4] == 'symlog'
f = netcdf_file(os.path.join(rundir, 'movies', 'swtc11.nc'), 'r', mmap=False)
bell = f.variables['geop'][:, lev - 1]   # (time, lat, lon)
lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
time = f.variables['time'][:]
tag = os.path.basename(os.path.normpath(rundir))

fig, ax = plt.subplots(figsize=(10, 5.4), constrained_layout=True)
if symlog:
    lt = 1e-4 * hmax                    # linear below 0.1 m (bell) / 1e-4 (slotted cylinder)
    norm = SymLogNorm(linthresh=lt, vmin=-hmax, vmax=hmax, base=10)
    pc = ax.pcolormesh(lon, lat, bell[0], cmap='RdBu_r', norm=norm, shading='auto')
    ticks = [s * 10.0**e for s in (-1, 1) for e in range(int(np.log10(lt)), int(np.log10(hmax)) + 1)] + [0]
    cb = fig.colorbar(pc, ax=ax, shrink=0.85, extend='both', ticks=sorted(ticks))
    cb.ax.set_yticklabels(['%g' % t for t in sorted(ticks)])
    cb.set_label('h' + (' [m]' if units else '') + ' (symmetric log scale)')
else:
    pc = ax.pcolormesh(lon, lat, bell[0], cmap='Blues', vmin=0, vmax=hmax, shading='auto')
    fig.colorbar(pc, ax=ax, shrink=0.85, extend='both', label='h' + (' [m]' if units else ''))
ax.set_xlim(0, 360); ax.set_ylim(-90, 90); ax.set_aspect('equal')
ax.set_xticks(range(0, 361, 60)); ax.set_yticks(range(-90, 91, 30))
ax.set_xlabel('longitude'); ax.set_ylabel('latitude')
ax.grid(alpha=0.3)
title = ax.set_title('')

def draw(it):
    pc.set_array(bell[it].ravel())
    title.set_text('swtc1 %s, %s   day %2d  (rev %.2f)   h: min %.3g, max %.3g%s'
                   % (name, tag, time[it], time[it] / 12, bell[it].min(), bell[it].max(), units))
    return pc, title

anim = FuncAnimation(fig, draw, frames=len(time), blit=False)
anim.save(out, writer=PillowWriter(fps=6), dpi=80)
print('wrote', out, len(time), 'frames')
