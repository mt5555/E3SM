#!/usr/bin/env python3
# Animate one swirl tracer (geop level) from movies/swirl1.nc as a GIF.
# usage: animate_swirl.py rundir out.gif [level] [symlog]
#   level: 1 cosine bells (default), 2 gaussian bells, 3 slotted cylinders
#   symlog: diverging color map on a symmetric log scale of q - background,
#           to show small +/- noise
# geop holds the mixing ratio q (tracer / density level 4)
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
name = ['cosine bells', 'gaussian bells', 'slotted cylinders'][lev - 1]
qbg = 0.0 if lev == 2 else 0.1          # background value
qmax = 0.95 if lev == 2 else 1.0        # peak value
symlog = len(sys.argv) > 4 and sys.argv[4] == 'symlog'
f = netcdf_file(os.path.join(rundir, 'movies', 'swirl1.nc'), 'r', mmap=False)
q = f.variables['geop'][:, lev - 1]     # (time, lat, lon)
lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
time = f.variables['time'][:]
tag = os.path.basename(os.path.normpath(rundir))

fig, ax = plt.subplots(figsize=(10, 5.4), constrained_layout=True)
if symlog:
    hmax = qmax - qbg
    lt = 1e-4 * hmax                    # linear within 1e-4 of the shape height
    norm = SymLogNorm(linthresh=lt, vmin=-hmax, vmax=hmax, base=10)
    pc = ax.pcolormesh(lon, lat, q[0] - qbg, cmap='RdBu_r', norm=norm, shading='auto')
    ticks = [s * 10.0**e for s in (-1, 1) for e in range(int(np.floor(np.log10(lt))) + 1, 0)] + [0]
    cb = fig.colorbar(pc, ax=ax, shrink=0.85, extend='both', ticks=sorted(ticks))
    cb.ax.set_yticklabels(['%g' % t for t in sorted(ticks)])
    cb.set_label('q - %g (symmetric log scale)' % qbg)
else:
    pc = ax.pcolormesh(lon, lat, q[0], cmap='Blues', vmin=qbg, vmax=qmax, shading='auto')
    fig.colorbar(pc, ax=ax, shrink=0.85, extend='both', label='q')
ax.set_xlim(0, 360); ax.set_ylim(-90, 90); ax.set_aspect('equal')
ax.set_xticks(range(0, 361, 60)); ax.set_yticks(range(-90, 91, 30))
ax.set_xlabel('longitude'); ax.set_ylabel('latitude')
ax.grid(alpha=0.3)
title = ax.set_title('')

def draw(it):
    pc.set_array((q[it] - qbg if symlog else q[it]).ravel())
    title.set_text('swirl %s, %s   day %4.1f   q: min %.3g, max %.3g'
                   % (name, tag, time[it], q[it].min(), q[it].max()))
    return pc, title

anim = FuncAnimation(fig, draw, frames=len(time), blit=False)
anim.save(out, writer=PillowWriter(fps=4), dpi=80)
print('wrote', out, len(time), 'frames')
