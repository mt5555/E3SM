#!/usr/bin/env python3
"""
Analysis and plots for opmat_proj.py: the swtc1 Gaussian's day-12 error from
the extracted operator, its eigen-decomposition, and where it lives on the
cubed sphere.

usage:
  opmat_proj_plot.py --proj DIR [DIR ...] --out OUT [--conv CONVDIR]

DIR holds ne{NE}-sc{sc}-f2-analytic-proj.npz from opmat_proj.py (the files with
eigenpairs are used for the modal plots; for each NE the file with eigenpairs
wins).  CONVDIR holds HOMME runs ne{NE}-sc{sc}-f2-dt{dt}/swtc1.gaus.{l2,linf}.errors.
Writes results_proj.txt and figures to OUT.
"""
import argparse
import glob
import os
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from opmat_eigs import cube_distances, style_ax, INK, MUTED, GRID, leff

ALPHA = {1: '0', 4: 'pi/4'}
SCCOL = {1: '#2a78d6', 4: '#eb6834'}
GOOD = 0.1          # |g| = |exp(lambda T) - 1| below this: the mode's phase error after one revolution < 0.016 rev

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': 0.6, 'legend.frameon': False,
                     'figure.dpi': 130})


def collect(dirs):
    files = {}
    for dr in dirs:
        for f in glob.glob(os.path.join(dr, 'ne*-sc*-f2-analytic-proj.npz')):
            ne, sc = map(int, re.match(r'ne(\d+)-sc(\d+)', os.path.basename(f)).groups())
            d = np.load(f)
            if (ne, sc) not in files or 'lam' in d:
                files[(ne, sc)] = f
    return files


def homme_errors(conv):
    out = {}
    if not conv:
        return out
    for dr in glob.glob(os.path.join(conv, 'ne*-sc*-f2-dt*')):
        ne, sc = map(int, re.match(r'ne(\d+)-sc(\d+)', os.path.basename(dr)).groups())
        try:
            l2 = float(open(os.path.join(dr, 'swtc1.gaus.l2.errors')).read().split()[-1])
            li = float(open(os.path.join(dr, 'swtc1.gaus.linf.errors')).read().split()[-1])
        except (OSError, IndexError, ValueError):
            continue
        out[(ne, sc)] = (l2, li)
    return out


def norms(e, p0, mass):
    return (np.sqrt((mass * np.abs(e) ** 2).sum() / (mass * p0 ** 2).sum()),
            np.abs(e).max() / np.abs(p0).max())


def orders(nes, errs):
    nes, errs = np.asarray(nes, float), np.asarray(errs, float)
    return np.log(errs[:-1] / errs[1:]) / np.log(nes[1:] / nes[:-1])


def edge_split(d, ne):
    """semi-discrete error split into the GLL nodes on cube edges and the rest"""
    e, p0, mass = d['e_sd_direct'].real, d['p0'], d['mass']
    dc, de = cube_distances(d['lat'], d['lon'])
    on = de < 1e-9
    nrm2, pmax = (mass * p0 ** 2).sum(), np.abs(p0).max()
    return dict(l2_on=np.sqrt((mass[on] * e[on] ** 2).sum() / nrm2),
                l2_off=np.sqrt((mass[~on] * e[~on] ** 2).sum() / nrm2),
                max_on=np.abs(e[on]).max() / pmax, max_off=np.abs(e[~on]).max() / pmax,
                area_on=mass[on].sum() / mass.sum(), de=de, dc=dc)


def profile(d, ne, nbin=40):
    """rms error (relative to max|p0|) vs distance to the nearest cube edge, in element widths"""
    e, mass = d['e_sd_direct'].real, d['mass']
    _, de = cube_distances(d['lat'], d['lon'])
    x = de / ((np.pi / 2) / ne)
    on = de < 1e-9
    bins = np.linspace(0, 4, nbin + 1)
    idx = np.digitize(x[~on], bins) - 1
    w, ee = mass[~on], e[~on] ** 2
    rms = np.array([np.sqrt((w[idx == i] * ee[idx == i]).sum() / w[idx == i].sum()) if (idx == i).any() else np.nan
                    for i in range(nbin)])
    rms_on = np.sqrt((mass[on] * e[on] ** 2).sum() / mass[on].sum())
    return 0.5 * (bins[1:] + bins[:-1]), rms / np.abs(d['p0']).max(), rms_on / np.abs(d['p0']).max()


def modal_summary(d):
    om, mass, p0 = float(d['omega']), d['mass'], d['p0']
    nrm = np.sqrt((mass * p0 ** 2).sum())
    a = np.abs(d['c']) / nrm
    g = np.abs(d['g_sd'])
    E2 = (g * a) ** 2
    tot = E2.sum()
    lam = d['lam'] / om
    good = g < GOOD
    o = np.argsort(-E2)
    cum = np.cumsum(E2[o]) / tot
    split = None
    if 'e_res' in d:
        er = d['e_res']
        eu = d['e_sd'] - er
        split = norms(er, p0, mass) + norms(eu, p0, mass)
    return dict(a=a, g=g, E2=E2, tot=tot, lam=lam, good=good, split=split,
                nmodes=[int(np.searchsorted(cum, f) + 1) for f in (0.5, 0.9, 0.99)],
                frac_bad=E2[~good].sum() / tot, a2_bad=(a[~good] ** 2).sum(),
                edge_w=(E2 * d['edge']).sum() / tot, corner_w=(E2 * d['corner']).sum() / tot,
                edge_a=(a ** 2 * d['edge']).sum() / (a ** 2).sum())


def cube_edge_lines(npt=60):
    """lon, lat (degrees) of the 12 cube edges, HOMME orientation (corners at (+-1,+-1,+-1))"""
    lines = []
    for ax in range(3):
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                t = np.linspace(-1, 1, npt)
                P = np.zeros((npt, 3))
                P[:, ax] = t
                P[:, (ax + 1) % 3] = s1
                P[:, (ax + 2) % 3] = s2
                P /= np.linalg.norm(P, axis=1)[:, None]
                lon = np.degrees(np.arctan2(P[:, 1], P[:, 0])) % 360
                lat = np.degrees(np.arcsin(P[:, 2]))
                brk = np.where(np.abs(np.diff(lon)) > 180)[0]
                for seg in np.split(np.arange(npt), brk + 1):
                    lines.append((lon[seg], lat[seg]))
    return lines


def plot_maps(files, nes, out):
    """log10 |semi-discrete error| / max|p0| at day 12 on the GLL nodes"""
    sel = [ne for ne in nes if all((ne, sc) in files for sc in (1, 4))]
    if not sel:
        return
    fig, axs = plt.subplots(len(sel), 2, figsize=(13, 3.2 * len(sel)), squeeze=False)
    edges = cube_edge_lines()
    for row, ne in enumerate(sel):
        for col, sc in enumerate((1, 4)):
            d = np.load(files[(ne, sc)])
            e = np.abs(d['e_sd_direct'].real) / np.abs(d['p0']).max()
            emax = e.max()
            ax = axs[row][col]
            o = np.argsort(e)
            sz = 0.5 if ne >= 60 else (1.5 if ne >= 30 else 4)
            s = ax.scatter(np.degrees(d['lon'][o]), np.degrees(d['lat'][o]), c=np.log10(np.maximum(e[o], 1e-30)),
                           s=sz, cmap='magma_r', vmin=np.log10(emax) - 3, vmax=np.log10(emax), linewidths=0)
            for lo, la in edges:
                ax.plot(lo, la, color='#2a78d6', lw=0.6)
            ax.set_xlim(0, 360); ax.set_ylim(-90, 90); ax.grid(False)
            style_ax(ax, f'NE{ne} alpha={ALPHA[sc]}: log10 |e|/max|p0| (max {emax:.2e}); blue: cube edges')
            fig.colorbar(s, ax=ax, shrink=0.9)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'proj_error_maps.png'))
    plt.close(fig)


def plot_top_modes(d, ne, sc, out):
    """|v|^2 of the eigenmodes with the largest error contributions (one of each conjugate pair)"""
    om = float(d['omega'])
    E, lam, V = d['top_E'], d['top_lam'] / om, d['top_v']
    keep = [0]
    for j in range(1, len(E)):
        if not any(abs(lam[j] - np.conj(lam[i])) < 1e-9 * abs(lam[j]) + 1e-12 and abs(E[j] - E[i]) < 1e-9 * E[i] for i in keep):
            keep.append(j)
    keep = keep[:4]
    fig, axs = plt.subplots(1, len(keep), figsize=(4.2 * len(keep), 2.8), squeeze=False)
    edges = cube_edge_lines()
    tot = np.sqrt((d['mass'] * d['p0'] ** 2).sum())
    enorm = np.sqrt((d['mass'] * np.abs(d['e_sd']) ** 2).sum())
    for ax, j in zip(axs[0], keep):
        w = np.abs(V[:, j]) ** 2
        o = np.argsort(w)
        ax.scatter(np.degrees(d['lon'][o]), np.degrees(d['lat'][o]), c=w[o] / w.max(), s=1 + 30 * w[o] / w.max(), cmap='Blues', linewidths=0)
        for lo, la in edges:
            ax.plot(lo, la, color='#eb6834', lw=0.6)
        ax.set_xlim(0, 360); ax.set_ylim(-90, 90); ax.grid(False)
        style_ax(ax, f'lambda/omega={lam[j].imag:+.3f}i, E_j/||e||={E[j] * tot / enorm:.3f}')
    fig.suptitle(f'NE{ne} alpha={ALPHA[sc]}: |v|^2 of the largest error contributors (orange: cube edges)',
                 x=0.01, ha='left', color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'proj_top_modes_ne{ne}_sc{sc}.png'))
    plt.close(fig)


def source_analysis(src, out, rep):
    """orders of the error driven by the truncation error of each region (opmat_source.py),
    and the error history scaled by NE^4"""
    rep.append('\n== semi-discrete error split by where the truncation error is generated (opmat_source.py) ==')
    for sc in (1, 4):
        nes = sorted(int(re.match(r'ne(\d+)', os.path.basename(f)).group(1))
                     for f in glob.glob(os.path.join(src, f'ne*-sc{sc}-f2-analytic-source.npz')))
        if not nes:
            continue
        rows = {}
        for ne in nes:
            d = np.load(os.path.join(src, f'ne{ne}-sc{sc}-f2-analytic-source.npz'))
            regs = [str(r) for r in d['regions']]
            rows[ne] = [norms(d['E'][:, j], d['p0'], d['mass']) for j in range(len(regs))] + \
                [norms(d['E'].sum(axis=1), d['p0'], d['mass'])]
        regs = regs + ['total']
        rep.append(f'alpha={ALPHA[sc]}: relative l2 / max per region; orders between successive NE')
        for ne in nes:
            rep.append(f'  NE{ne:3d}: ' + '  '.join(f'{r} {x[0]:.3e}/{x[1]:.3e}' for r, x in zip(regs, rows[ne])))
        ev = [ne for ne in nes if ne % 2 == 0]
        for a, b in zip(ev[:-1], ev[1:]):
            f = np.log(b / a)
            rep.append(f'  orders {a}->{b}: ' + '  '.join(
                f'{r} {np.log(x[0] / y[0]) / f:.2f}/{np.log(x[1] / y[1]) / f:.2f}' for r, x, y in zip(regs, rows[a], rows[b])))
    hist = {}
    for f in glob.glob(os.path.join(src, 'ne*-sc*-f2-analytic-history.npz')):
        ne, sc = map(int, re.match(r'ne(\d+)-sc(\d+)', os.path.basename(f)).groups())
        hist[(ne, sc)] = np.load(f)
    if not hist:
        return
    fig, axs = plt.subplots(3, 2, figsize=(13, 8.5), sharex=True)
    for j, sc in enumerate((1, 4)):
        nes = sorted(ne for ne, s in hist if s == sc and ne % 2 == 0)
        cols = plt.cm.viridis(np.linspace(0, 0.85, len(nes)))
        for ne, c in zip(nes, cols):
            d = hist[(ne, sc)]
            t = d['t'] / 86400
            axs[0][j].plot(t, d['l2'] * ne ** 4, color=c, lw=1.2, label=f'NE{ne}')
            axs[1][j].plot(t, d['max'] * ne ** 4, color=c, lw=1.2, label=f'NE{ne}')
        axs[2][j].plot(t, np.degrees(d['de']), color='#2a78d6', label='to the nearest cube edge')
        axs[2][j].plot(t, np.degrees(d['dc']), color='#eb6834', label='to the nearest cube corner')
        style_ax(axs[0][j], f'alpha={ALPHA[sc]}: relative l2 error x NE^4 (curves collapse if 4th order)')
        style_ax(axs[1][j], f'alpha={ALPHA[sc]}: relative max error x NE^4')
        style_ax(axs[2][j], 'distance of the Gaussian center (degrees); Gaussian e-folding radius 19 deg')
        axs[2][j].set_xlabel('days')
        for ax in axs[:, j]:
            ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'proj_error_history.png'))
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--proj', nargs='+', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--conv')
    p.add_argument('--source', help='directory with opmat_source.py output (*-source.npz, *-history.npz)')
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    files = collect(args.proj)
    homme = homme_errors(args.conv)
    rep = []
    res = {}

    # ------------------------------------------------ convergence and the cube-edge split
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    for sc in (1, 4):
        nes = sorted(ne for (ne, s) in files if s == sc)
        rows = []
        for ne in nes:
            d = np.load(files[(ne, sc)])
            fd = norms(d['e_fd_direct'], d['p0'], d['mass'])
            sd = norms(d['e_sd_direct'], d['p0'], d['mass'])
            es = edge_split(d, ne)
            md = modal_summary(d) if 'lam' in d else None
            res[(ne, sc)] = dict(fd=fd, sd=sd, es=es, md=md)
            rows.append((ne, fd, sd, es))
        rep.append(f'\n== alpha={ALPHA[sc]}: Gaussian day-12 error (relative l2 / max) ==')
        rep.append('NE    HOMME l2     HOMME max  | matrix RK4 l2  max        | semi-discrete l2  max      '
                   '| cube-edge nodes: area  l2         max   | other nodes: l2   max')
        for ne, fd, sd, es in rows:
            h = homme.get((ne, sc), (np.nan, np.nan))
            rep.append(f'{ne:4d}  {h[0]:.5e} {h[1]:.5e} | {fd[0]:.5e} {fd[1]:.5e} | {sd[0]:.5e} {sd[1]:.5e} | '
                       f'{es["area_on"]:.4f} {es["l2_on"]:.4e} {es["max_on"]:.4e} | {es["l2_off"]:.4e} {es["max_off"]:.4e}')
        nh = sorted(ne for (ne, s) in homme if s == sc)
        if len(nh) > 1:
            o2 = orders(nh, [homme[(ne, sc)][0] for ne in nh])
            oi = orders(nh, [homme[(ne, sc)][1] for ne in nh])
            rep.append('HOMME orders l2/max: ' + ', '.join(f'{a}->{b}: {x:.2f}/{y:.2f}'
                                                          for a, b, x, y in zip(nh[:-1], nh[1:], o2, oi)))
        ev = [r for r in rows if r[0] % 2 == 0 and r[0] >= 8]
        # even NE only (odd NE puts the equator, the alpha=0 path, through element centers)
        for lab, k in (('semi-discrete l2', lambda r: r[2][0]), ('semi-discrete max', lambda r: r[2][1]),
                       ('edge-node l2', lambda r: r[3]['l2_on']), ('edge-node max', lambda r: r[3]['max_on']),
                       ('other-node l2', lambda r: r[3]['l2_off']), ('other-node max', lambda r: r[3]['max_off'])):
            nn = [r[0] for r in ev]
            o = orders(nn, [k(r) for r in ev])
            rep.append(f'orders, {lab} (even NE): ' + ', '.join(f'{a}->{b}: {x:.2f}' for a, b, x in zip(nn[:-1], nn[1:], o)))

        col = SCCOL[sc]
        nn = np.array([r[0] for r in rows])
        axs[0].loglog(nn, [r[2][0] for r in rows], 'o-', color=col, ms=4, label=f'alpha={ALPHA[sc]} l2')
        axs[0].loglog(nn, [r[2][1] for r in rows], 's--', color=col, ms=4, mfc='none', label=f'alpha={ALPHA[sc]} max')
        if nh:
            axs[0].loglog(nh, [homme[(ne, sc)][0] for ne in nh], 'x', color=INK, ms=6)
            axs[0].loglog(nh, [homme[(ne, sc)][1] for ne in nh], '+', color=INK, ms=7)
        axs[1 if sc == 4 else 2].loglog(nn, [r[3]['max_on'] for r in rows], 'o-', color='#c2410c', ms=4, label='max, cube-edge nodes')
        axs[1 if sc == 4 else 2].loglog(nn, [r[3]['max_off'] for r in rows], 's-', color='#1d4ed8', ms=4, label='max, other nodes')
        axs[1 if sc == 4 else 2].loglog(nn, [r[3]['l2_on'] for r in rows], 'o--', color='#c2410c', ms=4, mfc='none', label='l2, cube-edge nodes')
        axs[1 if sc == 4 else 2].loglog(nn, [r[3]['l2_off'] for r in rows], 's--', color='#1d4ed8', ms=4, mfc='none', label='l2, other nodes')
    for ax in axs:
        x = np.array([8, 120.0])
        for k, ls in ((3, ':'), (4, '-.')):
            ax.loglog(x, 3e-1 * (x / 8) ** -k, ls, color=MUTED, lw=1, label=f'NE^-{k}')
        ax.set_xlabel('NE')
    axs[0].set_ylabel('relative error at day 12')
    style_ax(axs[0], 'Gaussian, semi-discrete exp(LT) (lines) and HOMME RK4 runs (x l2, + max)')
    style_ax(axs[1], 'alpha=pi/4: semi-discrete error on cube-edge GLL nodes vs the rest')
    style_ax(axs[2], 'alpha=0: same split')
    for ax in axs:
        ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'proj_convergence.png'))
    plt.close(fig)

    # ------------------------------------------------ error profile across the cube edges
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for ax, sc in zip(axs, (1, 4)):
        nes = sorted(ne for (ne, s) in files if s == sc and ne >= 15)
        cols = plt.cm.viridis(np.linspace(0, 0.85, len(nes)))
        for ne, c in zip(nes, cols):
            x, rms, rms_on = profile(np.load(files[(ne, sc)]), ne)
            ax.semilogy(x, rms * ne ** 4, '-', color=c, lw=1.4, label=f'NE{ne}')
            ax.semilogy([0], [rms_on * ne ** 4], 'o', color=c, ms=5)
        ax.set_xlabel('distance to the nearest cube edge (element widths); dot: nodes on the edge')
        style_ax(ax, f'alpha={ALPHA[sc]}: rms semi-discrete error x NE^4 (collapses if O(dx^4))')
        ax.legend(fontsize=7)
    axs[0].set_ylabel('rms error / max|p0|  x NE^4')
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'proj_edge_profile.png'))
    plt.close(fig)

    plot_maps(files, [20, 30, 60, 120], args.out)
    if args.source:
        source_analysis(args.source, args.out, rep)

    # ------------------------------------------------ modal budget
    mod = sorted(k for k in res if res[k]['md'] is not None)
    if mod:
        rep.append('\n== modal decomposition (eigenvectors of the skew part K of M^1/2 L M^-1/2) ==')
        rep.append(f'a_j = |c_j|/||p0||, E_j = |exp(lambda_j T)-1| a_j, "resolved" = |exp(lambda_j T)-1| < {GOOD}')
        rep.append('NE sc  n      #modes 50/90/99% of E^2   E^2 from unresolved  a^2 unresolved  '
                   'E^2-weighted edge / corner enrichment  (a^2-weighted edge)')
        for ne, sc in mod:
            m = res[(ne, sc)]['md']
            rep.append(f'{ne:2d} {sc}  {len(m["a"]):6d} {m["nmodes"]}  {m["frac_bad"]:.3f}  {m["a2_bad"]:.3e}  '
                       f'{m["edge_w"]:.2f} / {m["corner_w"]:.2f}  ({m["edge_a"]:.2f})')
        rep.append('error field from resolved / unresolved modes (relative l2, max):')
        for ne, sc in mod:
            sp_ = res[(ne, sc)]['md']['split']
            if sp_ is not None:
                rep.append(f'  NE{ne} alpha={ALPHA[sc]}: resolved l2={sp_[0]:.4e} max={sp_[1]:.4e}; '
                           f'unresolved l2={sp_[2]:.4e} max={sp_[3]:.4e}')
        for sc in (1, 4):
            nn = [ne for ne, s in mod if s == sc and ne % 2 == 0]
            if len(nn) > 1:
                r = [res[(ne, sc)]['md'] for ne in nn]
                ob = orders(nn, [np.sqrt(x['E2'][~x['good']].sum()) for x in r])
                og = orders(nn, [np.sqrt(x['E2'][x['good']].sum()) for x in r])
                oa = orders(nn, [np.sqrt(x['a2_bad']) for x in r])
                rep.append(f'alpha={ALPHA[sc]} orders (even NE): error from unresolved modes '
                           + ', '.join(f'{a}->{b}: {x:.2f}' for a, b, x in zip(nn[:-1], nn[1:], ob))
                           + '; from resolved modes ' + ', '.join(f'{x:.2f}' for x in og)
                           + '; amplitude in unresolved modes ' + ', '.join(f'{x:.2f}' for x in oa))

        # scatter: amplitude vs |Im lambda|, colored by |g|
        nes = sorted({ne for ne, _ in mod})
        sel = [ne for ne in nes if ne in (8, 12, 16, 20)] or nes
        fig, axs = plt.subplots(2, len(sel), figsize=(4 * len(sel), 7), sharey=True, squeeze=False)
        for col, ne in enumerate(sel):
            for row, sc in enumerate((1, 4)):
                ax = axs[row][col]
                if (ne, sc) not in res or res[(ne, sc)]['md'] is None:
                    ax.axis('off'); continue
                m = res[(ne, sc)]['md']
                o = np.argsort(m['g'])
                s = ax.scatter(np.abs(m['lam'].imag[o]), np.maximum(m['a'][o], 1e-16), c=np.log10(np.maximum(m['g'][o], 1e-12)),
                               s=3, cmap='plasma', vmin=-8, vmax=0.3, linewidths=0)
                ax.set_yscale('log'); ax.set_ylim(1e-14, 1)
                style_ax(ax, f'NE{ne} alpha={ALPHA[sc]}: a_j vs |Im lambda|/omega')
                if row == 1:
                    ax.set_xlabel('|Im lambda|/omega')
        for ax in axs[:, 0]:
            ax.set_ylabel('Gaussian amplitude a_j = |c_j|/||p0||')
        fig.colorbar(s, ax=axs, label='log10 |exp(lambda T) - 1|', shrink=0.8)
        fig.savefig(os.path.join(args.out, 'proj_amplitudes.png'))
        plt.close(fig)

        # error spectrum: E^2 binned by |Im lambda|/omega, resolved vs unresolved
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
        for ax, sc in zip(axs, (1, 4)):
            nn = [ne for ne, s in mod if s == sc]
            cols = plt.cm.viridis(np.linspace(0, 0.85, len(nn)))
            for ne, c in zip(nn, cols):
                m = res[(ne, sc)]['md']
                x = np.abs(m['lam'].imag)
                bins = np.arange(0, x.max() + 2)
                for good, ls in ((True, '-'), (False, ':')):
                    hsum = np.histogram(x[m['good'] == good], bins, weights=m['E2'][m['good'] == good])[0]
                    ax.semilogy(bins[:-1] + 0.5, np.maximum(hsum, 1e-30), ls, color=c, lw=1.2,
                                label=f'NE{ne} ' + ('resolved' if good else 'unresolved'))
            ax.set_ylim(1e-16, 1e-3)
            ax.set_xlim(0, 60)
            ax.set_xlabel('|Im lambda|/omega (unit bins)')
            style_ax(ax, f'alpha={ALPHA[sc]}: E_j^2 summed per bin')
            ax.legend(fontsize=6, ncol=2)
        axs[0].set_ylabel('error energy, relative to ||p0||^2')
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, 'proj_error_spectrum.png'))
        plt.close(fig)

    for sc in (1, 4):
        nn = [ne for ne, s in mod if s == sc and 'top_v' in np.load(files[(ne, sc)])]
        if nn:
            plot_top_modes(np.load(files[(max(nn), sc)]), max(nn), sc, args.out)

    txt = '\n'.join(rep)
    print(txt)
    with open(os.path.join(args.out, 'results_proj.txt'), 'w') as f:
        f.write(txt + '\n')


if __name__ == '__main__':
    main()
