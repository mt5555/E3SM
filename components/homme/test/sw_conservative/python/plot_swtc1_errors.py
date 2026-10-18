#!/usr/bin/env python3
# Plot swtc1 cosine bell relative errors (l2 and max) from the run directory.
# usage: plot_swtc1_errors.py [rundir [max_revolutions]]   (default: current directory, all output)
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rundir = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else '.')
days_per_rev = 12.0
ytop_3rev = 0.04     # y-axis top of the linear plot of the first 3 revolutions
max_rev = float(sys.argv[2]) if len(sys.argv) > 2 else None

# alpha for each sub_case; must match tc1_init_state in shallow_water_mod.F90
alpha_label = {1: '0', 2: 'pi/2', 3: 'pi/2 - 0.05', 4: 'pi/4'}
alpha_tag = {1: '0', 2: 'pi2', 3: 'pi-05', 4: 'pi4'}   # used in the file name
with open(os.path.join(rundir, 'swtc1.nl')) as fnl:
    nl = fnl.read()
sub_case = int(re.search(r'sub_case\s*=\s*(\d+)', nl).group(1))
m = re.search(r'theta_advect_form\s*=\s*(\d+)', nl)
adv_form = int(m.group(1)) if m else 0
m = re.search(r'^\s*nu\s*=\s*([^\s,!]+)', nl, re.MULTILINE)
nu = m.group(1) if m else None

l2 = np.loadtxt(os.path.join(rundir, 'swtc1.cosi.l2.errors'))
li = np.loadtxt(os.path.join(rundir, 'swtc1.cosi.linf.errors'))
if max_rev is not None:
    l2 = l2[l2[:, 0] <= max_rev*days_per_rev]
    li = li[li[:, 0] <= max_rev*days_per_rev]

fig, (a, b) = plt.subplots(1, 2, figsize=(13, 5))
for ax in (a, b):
    ax.plot(l2[:, 0]/days_per_rev, l2[:, 1], label='l2')
    ax.plot(li[:, 0]/days_per_rev, li[:, 1], label='max (linf)')
    ax.grid(True, which='both', alpha=.3)
    ax.set_xlabel('time (revolutions)')
    ax.set_ylabel('relative error')

nrev = 3
m = l2[:, 0] <= nrev*days_per_rev
# fixed y-axis so different runs can be compared; only expanded if the errors exceed it
ytop = max(ytop_3rev, 1.1*max(l2[m, 1].max(), li[m, 1].max()))
a.set_xlim(0, nrev)
a.set_ylim(0, ytop)
a.set_title('first %d revolutions (linear)' % nrev)
a.legend()

b.set_yscale('log')
b.set_ylim(1e-3, 10)
b.set_title('full run (log)' if max_rev is None else 'first %g revolutions (log)' % max_rev)

title = 'swtc1 cosine bell errors, alpha = %s' % alpha_label[sub_case]
if adv_form == 2:
    title += ', splitform'
if nu is not None:
    title += ', nu=%s' % nu
fig.suptitle(title)
fig.tight_layout()
out = os.path.join(rundir, 'swtc1_cosi_errors_%s.png' % alpha_tag[sub_case])
fig.savefig(out, dpi=120)
print('wrote', out)
