#!/usr/bin/env python3
"""
Figures from opmat_shnorm.py results: the worst-case error per spherical
harmonic degree after one revolution, in two panels (left: linear in the
normalized wavenumber k dx/pi over wavelengths >= 4.5 dx; right: log-log over
the resolved range).  Curves from different NE line up by wavelength in units
of dx, the mean GLL spacing.

usage: opmat_shnorm_figs.py --dir DIR [--cons DIR] --out OUT
writes OUT/fig_forms.png (split vs conservation form, nu=0 and nu>0),
OUT/fig_alpha.png (alpha=0 vs pi/4, NE10 and NE20) and OUT/fig_hypervis.png
(split form, nu=0 vs nu>0) for the cases whose files exist.
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter, FixedLocator, FuncFormatter

from opmat_eigs import style_ax, INK, MUTED, GRID

ACOL = {1: '#2a78d6', 4: '#eb6834'}          # alpha=0, alpha=pi/4
CCOL = '#6250d6'                              # conservation form
ALAB = {1: r'$\alpha=0$', 4: r'$\alpha=\pi/4$'}
PPW_MIN = 4.5                                 # beyond this, interpolated harmonics alias
SAT = 0.5                                     # worst-case error above this: saturated

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': 0.6, 'legend.frameon': False,
                     'figure.dpi': 130})


def case(path, label, color, ls='-', lw=1.6):
    """one curve: worst case E(N) against k dx/pi, or None if the file is missing"""
    if not os.path.exists(path):
        print('missing', path)
        return None
    d = np.load(path)
    dx = np.sqrt(4 * np.pi / int(d['n']))
    nmax = min(len(d['E']) - 1, int(2 * np.pi / (PPW_MIN * dx)))
    N = np.arange(1, nmax + 1)
    return dict(x=N * dx / np.pi, E=d['E'][N], label=label, color=color, ls=ls, lw=lw)


def twopanel(curves, title, out, note=None):
    curves = [c for c in curves if c is not None]
    if not curves:
        return
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.9))
    xsat = min(c['x'][np.argmax(c['E'] > SAT)] for c in curves if (c['E'] > SAT).any())
    for c in curves:
        for ax in axs:
            ax.plot(c['x'], c['E'], c['ls'], color=c['color'], lw=c['lw'], label=c['label'])
    for ax in axs:
        ax.set_yscale('log')
    axs[1].set_xscale('log')
    xmax = 2 / PPW_MIN
    axs[0].set_xlim(0, xmax)
    axs[0].axvspan(xsat, xmax, color=GRID, alpha=0.5, lw=0)
    axs[0].text(xsat + 0.005, 0.55, 'saturated\n(phase random)', transform=axs[0].get_xaxis_transform(),
                color=MUTED, fontsize=8, va='bottom')
    xlo = 0.9 * min(c['x'][0] for c in curves)
    axs[1].set_xlim(xlo, 1.05 * xsat)
    emin = min(c['E'].min() for c in curves)
    for ax in axs:
        ax.set_ylim(10 ** np.floor(np.log10(emin / 2)), 3)
    style_ax(axs[0], 'full resolved range: linear k')
    style_ax(axs[1], 'low wavenumbers: log-log')
    wl = FuncFormatter(lambda v, _: f'{v:g}' + r'$\Delta x$')
    inv = lambda x: 2 / np.maximum(x, 1e-6)
    for ax, ticks in ((axs[0], [5, 6, 8, 10, 15, 20, 40]), (axs[1], [8, 10, 15, 20, 30, 50, 100, 200])):
        ax.set_xlabel(r'$k\,\Delta x/\pi$   ($\Delta x$ = mean GLL spacing; 1 = node Nyquist)')
        ax.set_ylabel('worst-case relative l2 error after one revolution')
        top = ax.secondary_xaxis('top', functions=(inv, inv))
        top.xaxis.set_major_locator(FixedLocator(ticks))
        top.xaxis.set_minor_locator(FixedLocator([]))
        top.xaxis.set_major_formatter(wl)
        top.set_xlabel('wavelength', color=MUTED)
    axs[1].xaxis.set_major_formatter(ScalarFormatter())
    axs[1].xaxis.set_minor_formatter(NullFormatter())
    axs[1].set_xticks([t for t in (0.01, 0.02, 0.05, 0.1, 0.2) if xlo <= t <= 1.05 * xsat])
    axs[0].legend(fontsize=7.5, loc='lower right')
    axs[1].legend(fontsize=7.5, loc='upper left')
    fig.suptitle(title, x=0.01, ha='left', color=INK)
    if note:
        fig.text(0.01, 0.005, note, color=MUTED, fontsize=7.5, ha='left', va='bottom')
    fig.tight_layout(rect=(0, 0.03 if note else 0, 1, 1))
    fig.savefig(out)
    plt.close(fig)
    print('wrote', out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--dir', required=True, help='opmat_shnorm.py output (split form, conservation nu=0)')
    p.add_argument('--cons', help='opmat_shnorm.py output for the conservation form with nu > 0')
    p.add_argument('--out', required=True)
    p.add_argument('--ne', type=int, default=10)
    p.add_argument('--ne2', type=int, default=20, help='second resolution for the alpha figure')
    p.add_argument('--hv', default='nu30_1e+14_dx3', help='hyperviscosity tag of the nu > 0 runs')
    p.add_argument('--hvlabel', default=r'$\nu=2.7\times10^{15}$ ($10^{14}$ at NE30, $\propto\Delta x^3$)')
    p.add_argument('--sc', type=int, default=4, help='alpha (sub_case) for the forms figure')
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    ne, D, C = args.ne, args.dir, args.cons or args.dir

    def f(d, ne_, sc, form, tag):
        return os.path.join(d, f'ne{ne_}-sc{sc}-f{form}-analytic-ne{ne_}-sc{sc}-{tag}-sd-shnorm.npz')

    sc = args.sc
    twopanel([case(f(D, ne, sc, 2, 'nu0'), r'split form, $\nu=0$', ACOL[sc]),
              case(f(D, ne, sc, 0, 'nu0'), r'conservation form, $\nu=0$', CCOL),
              case(f(D, ne, sc, 2, args.hv), 'split form, ' + args.hvlabel, ACOL[sc], '--'),
              case(f(C, ne, sc, 0, args.hv), 'conservation form, ' + args.hvlabel, CCOL, '--')],
             f'NE{ne}, {ALAB[sc]}, semi-discrete: split vs conservation form, with and without hyperviscosity',
             os.path.join(args.out, 'fig_forms.png'))
    twopanel([case(f(D, ne, 1, 2, 'nu0'), f'{ALAB[1]}, NE{ne}', ACOL[1], '--', 1.2),
              case(f(D, ne, 4, 2, 'nu0'), f'{ALAB[4]}, NE{ne}', ACOL[4], '--', 1.2),
              case(f(D, args.ne2, 1, 2, 'nu0'), f'{ALAB[1]}, NE{args.ne2}', ACOL[1]),
              case(f(D, args.ne2, 4, 2, 'nu0'), f'{ALAB[4]}, NE{args.ne2}', ACOL[4])],
             r'split form, $\nu=0$, semi-discrete: flow along the grid vs across the cube edges',
             os.path.join(args.out, 'fig_alpha.png'),
             note='Curves from both resolutions line up by wavelength in units of the mean GLL spacing.')
    twopanel([case(f(D, ne, 1, 2, 'nu0'), f'{ALAB[1]}, ' + r'$\nu=0$', ACOL[1]),
              case(f(D, ne, 4, 2, 'nu0'), f'{ALAB[4]}, ' + r'$\nu=0$', ACOL[4]),
              case(f(D, ne, 1, 2, args.hv), f'{ALAB[1]}, ' + args.hvlabel, ACOL[1], '--'),
              case(f(D, ne, 4, 2, args.hv), f'{ALAB[4]}, ' + args.hvlabel, ACOL[4], '--')],
             f'NE{ne}, split form, semi-discrete: effect of hyperviscosity',
             os.path.join(args.out, 'fig_hypervis.png'))


if __name__ == '__main__':
    main()
