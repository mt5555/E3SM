"""
Figures from opmat_timeint.py runs: the variance history of every run, and maps of the tracer
at the snapshot times.

    python opmat_timeint_figs.py --runs OUT --center 30 -35 --out FIGDIR \
        [--maps f2:0 f0:0 f0:1e15] [--times 0 0.25 1 6] [--rates 0:40.75 2.5e15:0.388]

--rates NU30:RE draws the growth exp(2 Re lambda t) of the variance predicted by the rightmost
eigenvalue (Re lambda in units of omega) for that conservation-form run.  Each map panel has its
own symmetric colour scale, normalized by its max|q| (written in the panel), so the shape is
visible while the amplitude grows by hundreds of orders of magnitude.
"""
import argparse
import glob
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from opmat_sym import xyz
from opmat_noise_map import project, element_lines

INK, MUTED, GRID = '#22252a', '#6b7079', '#d9dce1'
SPLIT = '#2a78d6'            # split form: blue
CONS = plt.get_cmap('Oranges')   # conservation form: orange, darker = smaller nu (more unstable)


def load_runs(d):
    runs = {}
    for f in sorted(glob.glob(os.path.join(d, 'timeint-f*-nu*.npz'))):
        z = np.load(f)
        runs[(int(z['form']), float(z['nu30']))] = z
    return runs


def label(form, n30):
    nu = '0' if n30 == 0 else rf'{n30 / 1e15:g}\times10^{{15}}'
    return rf'{"split" if form == 2 else "conservation"}, $\nu_{{30}} = {nu}$'


def style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)


def variance_fig(runs, rates, out):
    om = 38.6396503844300128 / 6.376e6
    cons = sorted(n for f, n in runs if f == 0)
    colors = {n: CONS(0.95 - 0.5 * i / max(1, len(cons) - 1)) for i, n in enumerate(cons)}
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw=dict(width_ratios=[1, 1.25]))
    for ax, zoom in zip(axs, (False, True)):
        for (form, n30), z in sorted(runs.items(), key=lambda kv: (-kv[0][0], kv[0][1])):
            c = SPLIT if form == 2 else colors[n30]
            ls = '-' if (form == 0 or n30 == 0) else (0, (5, 2))
            ax.plot(z['t'], z['log10var'], color=c, lw=2, ls=ls, label=label(form, n30))
        for n30, re in rates.items():
            if (0, n30) in runs:
                t = runs[(0, n30)]['t']
                yl = 2 * re * om * t * 86400 / np.log(10)
                ax.plot(t, yl, color=INK, lw=0.9, ls=(0, (2, 2)), zorder=5)
                if zoom != (yl[-1] < 5):       # label each line in the panel where it fits
                    continue
                ylim = 5 if zoom else yl[-1]
                k = min(np.searchsorted(yl, 0.55 * ylim), len(t) - 1)
                ax.annotate(rf'Re$\,\lambda$ = {re:g}$\,\omega$', (t[k], yl[k]), xytext=(-8, 6),
                            textcoords='offset points', ha='right', fontsize=8.5, color=MUTED)
        style(ax)
        ax.set_xlabel('time (days)', color=INK)
        ax.set_xlim(0, max(z['t'][-1] for z in runs.values()))
    axs[0].set_ylabel('log₁₀ variance / initial variance', color=INK)
    axs[0].set_title('all runs', color=INK, fontsize=11, loc='left')
    axs[1].set_ylim(-0.6, 5)
    axs[1].set_title('near the threshold (zoom)', color=INK, fontsize=11, loc='left')
    # direct labels at the right end of the zoomed panel, for the runs that end inside it
    ends = sorted((float(z['log10var'][-1]), key) for key, z in runs.items() if -0.6 < z['log10var'][-1] < 5)
    ys, gap = [], 0.28                       # spread labels that would collide
    for y, _ in ends:
        ys.append(max(y, ys[-1] + gap) if ys else y)
    for (y, key), yl in zip(ends, ys):
        z = runs[key]
        axs[1].annotate(label(*key), (z['t'][-1], y), xytext=(z['t'][-1] + 0.4, yl), textcoords='data',
                        va='center', fontsize=8.5, color=INK)
    axs[1].set_xlim(0, runs[next(iter(runs))]['t'][-1] * 1.32)
    axs[0].legend(frameon=False, fontsize=8.5, loc='upper left', labelcolor=INK)
    axs[0].text(0.03, 0.52, 'dotted: exp(2 Re$\\,\\lambda\\,t$), $\\lambda$ the rightmost\neigenvalue of $L - \\nu B$', transform=axs[0].transAxes,
                ha='left', va='top', fontsize=8.5, color=MUTED)
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close(fig)
    print(out)


def maps_fig(runs, keys, times, center, radius, out):
    z0 = next(iter(runs.values()))
    X = xyz(z0['lat'], z0['lon'])
    c = xyz(np.radians([center[0]]), np.radians([center[1]]))[0]
    x, y = project(X, c)
    sel = (np.abs(x) < radius * 1.1) & (np.abs(y) < radius * 1.1) & (X @ c > 0)
    tri = mtri.Triangulation(x[sel], y[sel])
    lines = element_lines(int(z0['ne']), c, radius)
    fig, axs = plt.subplots(len(keys), len(times), figsize=(3.1 * len(times), 3.0 * len(keys)), squeeze=False)
    for i, key in enumerate(keys):
        z = runs[key]
        st = list(np.round(z['snap_t'], 6))
        for j, t in enumerate(times):
            ax = axs[i, j]
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_color(GRID)
            if round(t, 6) not in st:
                ax.text(0.5, 0.5, 'not saved', ha='center', va='center', color=MUTED, transform=ax.transAxes)
                continue
            k = st.index(round(t, 6))
            q = z['snap_q'][k]
            amp = np.abs(q).max()
            l10 = np.log10(amp) + float(z['snap_log10scale'][k])
            ax.tripcolor(tri, q[sel] / amp, cmap='RdBu_r', vmin=-1, vmax=1, shading='gouraud', rasterized=True)
            for P, cube in lines:
                px, py = project(P, c)
                ok = P @ c > 0
                ax.plot(np.where(ok, px, np.nan), np.where(ok, py, np.nan), color='k' if cube else '0.6',
                        lw=1.2 if cube else 0.3)
            ax.set_xlim(-radius, radius); ax.set_ylim(-radius, radius); ax.set_aspect('equal')
            ax.text(0.03, 0.03, rf'max$|q|$ = 10$^{{{l10:.1f}}}$' if abs(l10) >= 2 else rf'max$|q|$ = {10 ** l10:.3g}',
                    transform=ax.transAxes, fontsize=8, color=INK,
                    bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85))
            if i == 0:
                tl = f'{t * 24:g} h' if t < 1 else f'day {t:g}'
                ax.set_title('t = 0' if t == 0 else tl, color=INK, fontsize=10)
            if j == 0:
                ax.set_ylabel(label(*key), color=INK, fontsize=9.5)
    plt.tight_layout()
    plt.savefig(out, dpi=130)
    plt.close(fig)
    print(out)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--runs', required=True)
    p.add_argument('--center', type=float, nargs=2, default=[30.0, -35.0])
    p.add_argument('--radius', type=float, default=12.0)
    p.add_argument('--maps', nargs='+', default=['f2:0', 'f0:0', 'f0:1e15'])
    p.add_argument('--times', type=float, nargs='+', default=[0, 0.25, 1, 6])
    p.add_argument('--rates', nargs='*', default=['0:40.75', '2.5e15:0.388'])
    p.add_argument('--out', required=True)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    runs = load_runs(a.runs)
    rates = {float(k): float(v) for k, v in (r.split(':') for r in a.rates)}
    variance_fig(runs, rates, os.path.join(a.out, 'variance.png'))
    keys = [(int(f.lstrip('f')), float(n)) for f, n in (m.split(':') for m in a.maps)]
    maps_fig(runs, keys, a.times, a.center, a.radius, os.path.join(a.out, 'tracer-maps.png'))


if __name__ == '__main__':
    main()
