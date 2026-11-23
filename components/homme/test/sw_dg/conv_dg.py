#!/usr/bin/env python3
# Convergence figures for the swtc1-DG Run Gallery: Gaussian l2 (solid) and max (dashed) at day 12.
# usage: conv_dg.py RUNS_DIR OUT_DIR
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R, OUT = sys.argv[1], sys.argv[2]
NES = [8, 15, 30, 60, 120]

def last(run, norm):
    f = os.path.join(R, run, 'swtc1.gaus.%s.errors' % norm)
    if not os.path.exists(f):
        return None
    d = np.loadtxt(f, ndmin=2)
    if d.size == 0 or d[-1, 0] < 11.9:
        return None
    return d[-1, 1]

def series(fmt):
    ne, l2, mx = [], [], []
    for n in NES:
        a, b = last(fmt % n, 'l2'), last(fmt % n, 'linf')
        if a is not None:
            ne.append(n); l2.append(a); mx.append(b)
    return np.array(ne), np.array(l2), np.array(mx)

def order(ne, e):
    return np.log(e[-2] / e[-1]) / np.log(ne[-1] / ne[-2]) if len(ne) > 1 else np.nan

def plot(curves, title, fname):
    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
    for fmt, lab, c, mk in curves:
        ne, l2, mx = series(fmt)
        if len(ne) == 0:
            continue
        ax.loglog(ne, l2, '-', color=c, marker=mk, lw=2.4, ms=8,
                  label='%s   l2 %.2f, max %.2f  (NE%d-%d)' % (lab, order(ne, l2), order(ne, mx), ne[-2], ne[-1]))
        ax.loglog(ne, mx, '--', color=c, marker=mk, lw=1.2, ms=4, alpha=.8)
    x = np.array([15, 120.])
    for p, ls in ((3, ':'), (4, '-.')):
        ax.loglog(x, 1e-2 * (x / 15) ** -p, ls, color='0.45', lw=1)
        ax.text(x[-1] * 0.93, 1e-2 * (x[-1] / 15) ** -p, 'slope %d' % p, color='0.45', ha='right', va='bottom')
    ax.set_xticks(NES); ax.set_xticklabels([str(n) for n in NES]); ax.minorticks_off()
    ax.set_xlabel('NE  (elements per cube-face edge; dt = 22500 s / NE, RK4)')
    ax.set_ylabel('relative error at day 12  (solid: l2, dashed: max)')
    ax.grid(alpha=.3, which='both'); ax.set_title(title)
    ax.legend(title='orders between the two finest grids', loc='lower left', fontsize=9)
    fig.savefig(os.path.join(OUT, fname), dpi=110); plt.close(fig)

plot([('conv_s1_ne%d_dg1', 'upwind, α = 0 (equator)', '#2c7a4b', 'o'),
      ('conv_s4_ne%d_dg1', 'upwind, α = π/4 (cube corners)', '#1f5f8b', 'D'),
      ('conv_s1_ne%d_dg0', 'central, α = 0', '#e08a2e', 's'),
      ('conv_s4_ne%d_dg0', 'central, α = π/4', '#8e3b8e', 'v')],
     'swtc1-DG Gaussian, split form: upwind vs central flux', 'conv_flux.png')
plot([('conv_s4_ne%d_dg1', 'upwind, split form', '#1f5f8b', 'D'),
      ('convf0_s4_ne%d_dg1', 'upwind, conservation form', '#6fa8d6', 'o'),
      ('conv_s4_ne%d_dg0', 'central, split form', '#8e3b8e', 'v'),
      ('convf0_s4_ne%d_dg0', 'central, conservation form', '#d29ad2', 's')],
     'swtc1-DG Gaussian, α = π/4: split vs conservation form', 'conv_forms.png')
for fmt in ('conv_s4_ne%d_dg1', 'conv_s4_ne%d_dg0', 'conv_s1_ne%d_dg1', 'conv_s1_ne%d_dg0', 'convf0_s4_ne%d_dg1', 'convf0_s4_ne%d_dg0'):
    ne, l2, mx = series(fmt)
    print(fmt, list(ne), ['%.3g' % v for v in l2], ['%.3g' % v for v in mx])
