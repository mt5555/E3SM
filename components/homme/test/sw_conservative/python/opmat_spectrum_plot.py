#!/usr/bin/env python3
"""
Complex-plane plots of whole spectra from opmat_dense.py: the convex hull (or a scatter) of
the eigenvalues of each file, with the RK4 stability boundary for given time steps.

usage:
  opmat_spectrum_plot.py --files F1.npz [F2.npz ...] --out OUT.png [--labels L1 L2 ...]
                         [--dt 1800 ...] [--scatter] [--zoom XMIN] [--title TEXT]

F*.npz are opmat_dense.py outputs (lam already divided by omega = u0/a, plus omega, nu30,
case).  Hull mode (the default) fills each spectrum's convex hull, so several spectra fit on
one axis; a spectrum on the imaginary axis (|Re| below 1e-6 of its extent, e.g. the split
form at nu = 0) is drawn as a segment of the axis.  --scatter plots every eigenvalue instead
(one panel per file).  --dt draws the boundary |R(dt lam)| = 1 of RK4 (R from
opmat_eigs.rk4_amp) for each dt in seconds.  --zoom XMIN adds a panel limited to
Re lam/omega >= XMIN (and wide enough on the right for any Re lam > 0).  The legend gives each spectrum's range of Re and its max |Im|.

Example (NE30 solid body, alpha = pi/4, split form, nu = 0, 1e14, 1e15, from
`opmat_dense.py --sym --nu30 0 1e14 1e15`):
  opmat_spectrum_plot.py --files D/ne30-sc4-f2-analytic-dense-nu{1e+14,1e+15,0}.npz \\
      --dt 1800 --zoom -70 --out hull.png
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

from opmat_eigs import INK, MUTED, GRID, rk4_amp

SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']   # fixed order


def load(f):
    d = np.load(f, allow_pickle=True)
    return d['lam'], float(d['omega']), d


def default_label(f, d):
    case = str(d['case']) if 'case' in d else os.path.basename(f)
    nu = rf', $\nu_{{30}}$ = {float(d["nu30"]):.3g}' if 'nu30' in d else ''
    return case + nu


def hull_outline(lam):
    """closed outline (x, y) of the convex hull, or None if the spectrum lies on a vertical line"""
    P = np.c_[lam.real, lam.imag]
    width, height = np.ptp(P[:, 0]), np.ptp(P[:, 1])
    if width <= 1e-6 * max(height, 1e-300):
        return None
    h = P[ConvexHull(P).vertices]
    h = np.vstack([h, h[:1]])
    return h[:, 0], h[:, 1]


def style(ax):
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.grid(True, color=GRID, lw=0.6)
    for s in ax.spines.values():
        s.set_color(MUTED)
    ax.set_xlabel('Re λ / ω')


def rk4_boundary(axs, dts, om, xlim, ylim):
    """contour of |R(dt lam)| = 1 on each axis, over the region the spectra occupy"""
    X, Y = np.meshgrid(np.linspace(*xlim, 1200), np.linspace(*ylim, 1200))
    for j, dt in enumerate(dts):
        R = np.abs(rk4_amp((X + 1j * Y) * dt * om))
        ls = ['--', ':', '-.'][j % 3]
        for ax in axs:
            ax.contour(X, Y, R, levels=[1], colors=INK, linewidths=1.5, linestyles=ls)
            ax.plot([], [], color=INK, lw=1.5, ls=ls, label=f'RK4 stability boundary, dt = {dt:g} s')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--files', nargs='+', required=True)
    p.add_argument('--labels', nargs='+')
    p.add_argument('--out', required=True)
    p.add_argument('--dt', type=float, nargs='*', default=[])
    p.add_argument('--scatter', action='store_true')
    p.add_argument('--zoom', type=float, help='extra panel with Re lam/omega >= ZOOM')
    p.add_argument('--title')
    a = p.parse_args()
    if a.labels and len(a.labels) != len(a.files):
        raise SystemExit('--labels needs one label per file')
    if len(a.files) > len(SERIES):
        raise SystemExit(f'at most {len(SERIES)} files per plot')

    data = []
    for i, f in enumerate(a.files):
        lam, om, d = load(f)
        lab = a.labels[i] if a.labels else default_label(f, d)
        lab += f'  (Re λ/ω in [{lam.real.min():.3g}, {lam.real.max():.2g}], |Im λ|/ω ≤ {np.abs(lam.imag).max():.0f})'
        data.append((lam, om, lab, SERIES[i]))
    oms = [om for _, om, _, _ in data]
    if max(oms) - min(oms) > 1e-12 * max(oms):
        raise SystemExit('the files have different omega')
    om = oms[0]

    npan = (len(data) if a.scatter else 1) + (a.zoom is not None)
    fig, axs = plt.subplots(1, npan, figsize=(5.5 * npan + (1 if npan == 1 else 0), 4.8), sharey=True,
                            constrained_layout=True, squeeze=False)
    axs = list(axs[0])
    main_axs = axs[:-1] if a.zoom is not None else axs
    zoom_ax = axs[-1] if a.zoom is not None else None

    for i, (lam, _, lab, c) in enumerate(data):
        targets = [main_axs[i] if a.scatter else main_axs[0]] + ([zoom_ax] if zoom_ax else [])
        for ax in targets:
            if a.scatter:
                ax.scatter(lam.real, lam.imag, s=1.5, lw=0, color=c, rasterized=True, label=lab)
                continue
            o = hull_outline(lam)
            if o is None:                       # on the imaginary axis: a segment
                x = lam.real.mean()
                ax.plot([x, x], [lam.imag.min(), lam.imag.max()], color=c, lw=3, solid_capstyle='butt',
                        label=lab, zorder=5)
                ax.plot([x, x], [lam.imag.min(), lam.imag.max()], 'o', color=c, ms=8, mec='white', mew=2,
                        zorder=6)
            else:
                ax.fill(o[0], o[1], color=c, alpha=0.15, lw=0)
                ax.plot(o[0], o[1], color=c, lw=2, label=lab)

    allre = np.concatenate([lam.real for lam, *_ in data])
    allim = np.concatenate([lam.imag for lam, *_ in data])
    ymax = np.abs(allim).max()
    xmin = allre.min()
    if a.dt:
        reach = max(2.83 / (dt * om) for dt in a.dt)    # RK4 reaches 2.83 (imag) / 2.79 (real)
        ymax = max(ymax, 1.04 * reach)
        xmin = min(xmin, -1.0 * reach)
        span = max(-xmin, ymax)
        rk4_boundary(axs, a.dt, om, (xmin - 0.05 * span, 0.15 * span), (-1.1 * ymax, 1.1 * ymax))
    for ax in axs:
        style(ax)
        ax.set_ylim(-1.1 * ymax, 1.1 * ymax)
    axs[0].set_ylabel('Im λ / ω')
    if zoom_ax:
        zoom_ax.set_xlim(a.zoom, max(-0.07 * a.zoom, 1.15 * allre.max()))
        zoom_ax.set_title(f'zoom: Re λ/ω ≥ {a.zoom:g}', loc='left', color=INK, fontsize=9)
    for ax in main_axs:
        ax.legend(loc='upper left', fontsize=7.5, frameon=False, markerscale=6 if a.scatter else 1)
    if not a.scatter:
        main_axs[0].set_title('convex hull of each spectrum', loc='left', color=INK, fontsize=9)
    fig.suptitle(a.title or f'eigenvalues, ω = u0/a = {om:.4g} 1/s', color=INK, fontsize=10, x=0.01, ha='left')
    fig.savefig(a.out, dpi=150)
    print(a.out)


if __name__ == '__main__':
    main()
