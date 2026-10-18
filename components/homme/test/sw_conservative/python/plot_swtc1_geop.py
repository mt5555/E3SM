#!/usr/bin/env python3
# Plot swtc1 geop movie output (movies/swtc11.nc): shapes at day 0 vs day 36, and the bell's path.
# usage: plot_swtc1_geop.py rundir outdir
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.io import netcdf_file

rundir, outdir = sys.argv[1], sys.argv[2]
f = netcdf_file(os.path.join(rundir, 'movies', 'swtc11.nc'), 'r', mmap=False)
geop = f.variables['geop'][:]          # (time, lev, lat, lon)
lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
time = f.variables['time'][:]
names = ['cosine bell', 'gaussian', 'cylinder', 'slotted cylinder']
heights = [1000, 1000, 1000, 1]
tag = os.path.basename(os.path.normpath(rundir))

# ---- figure 1: day 0, day 36, difference, zoomed on each shape --------------
i0 = 0
i1 = int(np.argmin(abs(time - 36)))
# window centered on the initial bell
j, i = np.unravel_index(np.argmax(geop[0, 1]), geop[0, 1].shape)
latc, lonc = lat[j], lon[i]
jw = abs(lat - latc) <= 35
dl = (lon - lonc + 180) % 360 - 180
iw = abs(dl) <= 50
order = np.argsort(dl[iw])
x = dl[iw][order] + lonc
y = lat[jw]

fig, axes = plt.subplots(4, 3, figsize=(13, 13), constrained_layout=True)
for k in range(4):
    a = geop[i0, k][np.ix_(jw, iw)][:, order]
    b = geop[i1, k][np.ix_(jw, iw)][:, order]
    h = heights[k]
    for col, (fld, title) in enumerate([(a, 'day %g' % time[i0]), (b, 'day %g' % time[i1])]):
        ax = axes[k, col]
        pc = ax.pcolormesh(x, y, fld, cmap='Blues', vmin=0, vmax=h, shading='auto')
        ax.set_title('%s, %s   (min %.3g, max %.3g)' % (names[k], title, fld.min(), fld.max()),
                     fontsize=10)
    fig.colorbar(pc, ax=axes[k, :2], shrink=0.9, label='h' + (' [m]' if h > 1 else ''))
    d = b - a
    lim = max(abs(d).max(), 1e-12)
    ax = axes[k, 2]
    pc = ax.pcolormesh(x, y, d, cmap='RdBu_r', norm=TwoSlopeNorm(0, -lim, lim), shading='auto')
    ax.set_title('%s, day %g minus day %g' % (names[k], time[i1], time[i0]), fontsize=10)
    fig.colorbar(pc, ax=ax, shrink=0.9)
for ax in axes.flat:
    ax.set_aspect('equal')
    ax.tick_params(labelsize=8)
for ax in axes[-1]:
    ax.set_xlabel('longitude')
for ax in axes[:, 0]:
    ax.set_ylabel('latitude')
fig.suptitle('swtc1 %s: after 3 revolutions (36 days)' % tag)
out1 = os.path.join(outdir, 'swtc1_geop_day0_vs_day36.png')
fig.savefig(out1, dpi=110)

# ---- figure 2: cosine bell path during the first revolution ------------------
fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
days = [d for d in range(0, 13, 2)]
cmap = plt.get_cmap('Blues')
for n, d in enumerate(days):
    it = int(np.argmin(abs(time - d)))
    c = cmap(0.35 + 0.65 * n / (len(days) - 1))
    ax.contour(lon, lat, geop[it, 0], levels=[100, 500, 900], colors=[c], linewidths=1.5)
    if d == 12:
        continue   # back at the start: labeled together with day 0
    jj, ii = np.unravel_index(np.argmax(geop[it, 0]), geop[it, 0].shape)
    label = 'day 0 & 12' if d == 0 else 'day %g' % time[it]
    ax.annotate(label, (lon[ii], lat[jj]), fontsize=9, ha='center',
                va='bottom', xytext=(0, 6), textcoords='offset points')
ax.set_xlim(0, 360); ax.set_ylim(-90, 90); ax.set_aspect('equal')
ax.set_xlabel('longitude'); ax.set_ylabel('latitude')
ax.grid(alpha=0.3)
ax.set_title('swtc1 %s: cosine bell every 2 days over the first revolution '
             '(contours at h = 100, 500, 900 m)' % tag, fontsize=11)
out2 = os.path.join(outdir, 'swtc1_bell_path.png')
fig.savefig(out2, dpi=110)
print('wrote', out1, out2)
