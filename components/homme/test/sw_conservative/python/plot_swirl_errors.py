#!/usr/bin/env python3
# Plot swirl relative errors (l2 and max) for levels 1-4 vs time.
# usage: plot_swirl_errors.py [rundir]   (default: current directory)
# The swirl.*.errors files only have level 1, so all levels are read from the
# 'k= ... l1,l2,linf=' lines in sweq.out (swirl-splitform.sh tees the log there).
# Levels 5-6 (toy chemistry) are left out: they do not return to their initial state.
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rundir = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else '.')
days_per_period = 12.0
names = {1: 'cosine bells', 2: 'gaussian bells', 3: 'slotted cylinders', 4: 'density'}

with open(os.path.join(rundir, 'swirl.nl')) as fnl:
    nl = fnl.read()
sub_case = int(re.search(r'sub_case\s*=\s*(\d+)', nl).group(1))
m = re.search(r'theta_advect_form\s*=\s*(\d+)', nl)
adv_form = int(m.group(1)) if m else 0
m = re.search(r'^\s*nu\s*=\s*([^\s,!]+)', nl, re.MULTILINE)
nu = m.group(1) if m else None
m = re.search(r'^\s*tstep\s*=\s*([^\s,!]+)', nl, re.MULTILINE)
tstep = m.group(1) if m else None

# written with format (a,i2,f6.2,a,3e15.7): level and day can run together, e.g. 'k= 3263.00'
err = {}
with open(os.path.join(rundir, 'sweq.out')) as f:
    for line in f:
        if line.startswith('k=') and 'l1,l2,linf=' in line:
            k = int(line[2:4])
            day = float(line[4:10])
            l1, l2, li = (float(x) for x in line.split('=')[-1].split())
            err.setdefault(k, []).append((day, l1, l2, li))

# tracer variance rho*q^2 per level, if this run wrote it: columns are time, level 1..nlev
varfile = os.path.join(rundir, 'sweq.var')
var = np.loadtxt(varfile) if os.path.exists(varfile) else None

# with hyperviscosity on rho the density error can be tiny; leave it out so the
# tracer curves are readable, and say so on the plot
skip_density = False
if 4 in err:
    e = np.array(err[4])
    e = e[(e[:, 0] > 0) & np.isclose(np.mod(e[:, 0], days_per_period), 0)]
    skip_density = len(e) > 0 and e[:, 2].max() < 1e-6

if var is None:
    fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5))
else:
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(19, 5))
for k in sorted(names):
    if k not in err or (k == 4 and skip_density):
        continue
    e = np.array(err[k])
    e = e[e[:, 0] > 0]                   # errors are 0 at day 0
    # the exact solution is the initial state only at whole periods; in between the
    # 'error' is mostly the shapes being somewhere else, so only plot whole periods
    p = e[np.isclose(np.mod(e[:, 0], days_per_period), 0)]
    for ax, col in ((a, 2), (b, 3)):
        ax.plot(p[:, 0]/days_per_period, p[:, col], 'o-', ms=3,
                label='%d %s' % (k, names[k]))
for ax, norm in ((a, 'l2'), (b, 'max (linf)')):
    ax.set_yscale('log')
    ax.axhline(1, color='k', lw=0.8, alpha=0.5)
    ax.grid(True, which='both', alpha=.3)
    ax.set_xlabel('time (periods of 12 days)')
    ax.set_ylabel('relative error')
    ax.set_title(norm + ' at whole periods')
    ax.set_xlim(left=0)
    if skip_density:
        ax.text(0.02, 0.02, '4 density: l2 error below 1e-6, not plotted',
                transform=ax.transAxes, fontsize=9, color='0.3')
a.legend()

if var is not None:
    for k in (1, 2, 3):                  # tracers; level 4 is the density
        c.plot(var[:, 0]/days_per_period, var[:, k]/var[0, k] - 1, label='%d %s' % (k, names[k]))
    c.set_yscale('symlog', linthresh=1e-8)
    c.axhline(0, color='k', lw=0.8, alpha=0.5)
    c.grid(True, which='both', alpha=.3)
    c.set_xlabel('time (periods of 12 days)')
    c.set_ylabel('(V - V0)/V0')
    c.set_title('tracer variance V = integral of rho q^2')
    c.set_xlim(left=0)
    c.legend()

title = 'swirl errors, sub_case %d' % sub_case
title += ', splitform' if adv_form == 2 else ', conservation form'
if nu is not None:
    title += ', nu=%s' % nu
if tstep is not None:
    title += ', dt=%s' % tstep
fig.suptitle(title)
fig.tight_layout()
out = os.path.join(rundir, 'swirl_errors.png')
fig.savefig(out, dpi=120)
print('wrote', out)
