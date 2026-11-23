#!/usr/bin/env python3
"""
Initial-condition-independent accuracy of the swtc1 transport operator, per
spherical harmonic degree N:

    E(N) = || (exp(A T) - I) Q_N ||_2      (T = one revolution, A = L - nu B)

in the mass-weighted l2 norm, where the columns of Q_N are the 2N+1 degree-N
spherical harmonics at the GLL nodes, orthonormalized in the mass inner
product.  For every initial condition of degree N, ||e(T)|| <= E(N) ||p0||,
with equality for the worst one (the leading right singular vector).  The
degree-N space is closed under rotation, so E(N) covers every degree-N
pattern in every position on the grid.

All harmonics up to --nmax are propagated as one block with
scipy.sparse.linalg.expm_multiply (semi-discrete), or with RK4 at
dt = dt30*30/NE/sub (--rk4; sub=1 is HOMME's time step).  Also reported per
degree: the in-degree part ||Q_N^H X|| (phase error within degree N) and the
scattered part ||(I - Q_N Q_N^H) X|| (energy sent to other degrees), the RMS
over the degree (||X||_F/sqrt(2N+1)), and, as a check, the bound
sum_N E(N) ||p0_N|| for the swtc1 Gaussian against its actual error.

usage: opmat_shnorm.py --mats DIR --out OUT --ne 10 [--sc 1] [--forms 2 0]
                       [--nmax 60] [--nu30 0] [--rk4 [--sub 4]] [--dense]
--dense also computes exp(A T) as a dense matrix (small NE only) to check the
propagation and to give ||exp(A T)||_2.
"""
import argparse
import os
import time
import numpy as np
import scipy.linalg as la
import scipy.sparse.linalg as spla
from scipy.special import sph_harm_y
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from opmat_io import load
from opmat_sym import xyz
from opmat_eigs import name, omega, rotation_axis, style_ax, INK, MUTED, GRID
from opmat_proj import gaussian, timestep, rk4_run, RR, LATC, LONC

FORM = {0: 'conservation form', 2: 'split form'}
FCOL = {2: '#eb6834', 0: '#6250d6'}
ALPHA = {1: '0', 4: 'pi/4'}
SCCOL = {1: '#2a78d6', 4: '#eb6834'}

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': 0.6, 'legend.frameon': False,
                     'figure.dpi': 130})


def real_harmonics(lat, lon, alpha, nmax):
    """real spherical harmonics of degree 0..nmax at the nodes, pole on the swtc1 rotation
    axis; returns Y (n x (nmax+1)^2) and the degree of each column"""
    X = xyz(lat, lon)
    a, e1, e2 = rotation_axis(alpha)
    colat = np.arccos(np.clip(X @ a, -1, 1))
    az = np.arctan2(X @ e2, X @ e1)
    cols, deg = [], []
    for N in range(nmax + 1):
        for m in range(-N, N + 1):
            y = sph_harm_y(N, abs(m), colat, az)
            cols.append(y.real if m == 0 else np.sqrt(2) * (y.real if m > 0 else y.imag))
            deg.append(N)
    return np.stack(cols, axis=1), np.array(deg)


def cosine_bell(lat, lon):
    """swtc1 level 1 (tc1_geopotential, k=1): h0/2 (1 + cos(pi r/R)) for r < R, R = 1/3 radian,
    same center as the Gaussian"""
    A = np.sin(LATC) * np.sin(lat) + np.cos(LATC) * np.cos(lat) * np.cos(lon - LONC)
    r = np.arccos(np.clip(A, -1, 1))
    return np.where(r < RR, 500.0 * (1.0 + np.cos(np.pi * r / RR)), 0.0)


def degree_spectrum(p, Qs, sq):
    """||p_N|| / ||p|| per degree N in the mass norm, with Qs the mass-orthonormal bases per degree"""
    q = sq * p
    return np.array([np.linalg.norm(Q.T @ q) for Q in Qs]) / np.linalg.norm(q)


def degree_bases(Y, deg, sq):
    """mass-orthonormal basis Q_N of each degree (W^1/2 Y_N = Q R), as in degree_norms"""
    return [np.linalg.qr(sq[:, None] * Y[:, deg == N])[0] for N in np.unique(deg)]


def propagate(L, Y, T, dt=None, nstep=None):
    if dt is None:
        return spla.expm_multiply(L * T, Y)
    return rk4_run(L, Y, dt, nstep)


def degree_norms(P, Y, deg, sq):
    """per degree: largest singular value, RMS, in-degree and scattered parts of
    X = W^1/2 P R^-1 - Q, with W^1/2 Y_N = Q R"""
    out = {k: [] for k in ('E', 'rms', 'in', 'out')}
    Qs = []
    for N in np.unique(deg):
        c = deg == N
        Q, R = np.linalg.qr(sq[:, None] * Y[:, c])
        X = la.solve_triangular(R, (sq[:, None] * P[:, c]).T, trans='T').T - Q
        QX = Q.T @ X
        out['E'].append(la.svdvals(X)[0])
        out['rms'].append(np.linalg.norm(X) / np.sqrt(c.sum()))
        out['in'].append(la.svdvals(QX)[0])
        out['out'].append(la.svdvals(X - Q @ QX)[0])
        Qs.append(Q)
    return {k: np.array(v) for k, v in out.items()}, Qs


def plot(curves, title, n, nmax, ppw_min, outfile):
    """curves: list of (label, color, result dict with E, rms, in, out, cN)"""
    dx = np.sqrt(4 * np.pi / n)               # mean node spacing on the unit sphere
    xmax = min(nmax, 2 * np.pi / (ppw_min * dx))   # degree with wavelength ppw_min * dx
    nplot = min(nmax, int(np.ceil(xmax)))           # one point past xmax so curves reach the edge
    N = np.arange(nplot + 1)
    fig, axs = plt.subplots(1, 3, figsize=(19, 4.6))
    for lab, col, r in curves:
        r = {k: v[:nplot + 1] for k, v in r.items()}
        # degree 0 is plotted too; the y-limits below come from N >= 1, so a roundoff-level
        # N = 0 value (e.g. the in-degree part, or E(0) for a divergence-free velocity)
        # falls below the axis instead of stretching it
        axs[0].semilogy(N, r['E'], '-', color=col, lw=1.6, label=f'{lab}: worst case E(N)')
        axs[1].semilogy(N, r['in'], '-', color=col, lw=1.6, label=f'{lab}: in-degree')
        axs[1].semilogy(N, r['out'], '--', color=col, lw=1.4, label=f'{lab}: scattered')
    # third panel: degree spectra of the swtc1 initial conditions (independent of form and nu)
    axs[2].semilogy(N, curves[0][2]['cN'][:nplot + 1], color=INK, lw=1.4, label='Gaussian')
    if 'cN_bell' in curves[0][2]:
        axs[2].semilogy(N, curves[0][2]['cN_bell'][:nplot + 1], color=INK, lw=1.4, ls='--',
                        label='cosine bell')
    fig.suptitle(title, color=INK, fontsize=10)    # the case; panel titles say what is plotted
    style_ax(axs[0], 'worst-case error after one revolution: E(N) = ||(exp(AT) - I) Q_N||')
    style_ax(axs[1], 'worst-case error split: within degree N vs scattered to other degrees')
    style_ax(axs[2], 'degree spectrum of the swtc1 initial conditions')
    # left: fit the error curves; middle: the in-degree part reaches much lower values
    emin = min(r['E'][1:nplot + 1].min() for _, _, r in curves)
    axs[0].set_ylim(10 ** np.floor(np.log10(emin / 2)), 3)
    imin = min(r['in'][1:nplot + 1].min() for _, _, r in curves)
    axs[1].set_ylim(10 ** np.floor(np.log10(imin / 2)), 3)
    axs[2].set_ylim(1e-8, 1)
    for i, ax in enumerate(axs):
        degree_axis(ax, dx, xmax, bound=i < 2)
    for ax in axs:
        ax.legend(fontsize=7, loc='lower right')
    axs[0].set_ylabel('relative l2 error at day 12')
    axs[2].set_ylabel('||p0_N|| / ||p0||')
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)


def degree_axis(ax, dx, xmax, bound=True):
    """x axis in degree N with a top axis in wavelength / mean dx; bound: the line E = 2"""
    ax.set_xlabel('spherical harmonic degree N')
    ax.set_xlim(0, xmax)
    if bound:
        ax.axhline(2, color=MUTED, lw=0.8, ls='--')
    top = ax.secondary_xaxis('top', functions=(lambda x: 2 * np.pi / (np.maximum(x, 1e-9) * dx),
                                                lambda y: 2 * np.pi / (np.maximum(y, 1e-9) * dx)))
    top.set_xlabel(r'wavelength $2\pi a/N$ in units of mean $\Delta x$', color=MUTED)
    top.set_xticks([5, 6, 8, 10, 15, 20, 40])


SCAT_COL = {'in': '#eb6834', 'out': '#2a9d8f'}    # in-degree (split-form orange), scattered


def plot_compare(curves, title, n, nmax, ppw_min, outfile, split=None):
    """E(N) of several cases on one axis, plus the initial-condition spectra.
    curves: list of (label, color, linestyle, result dict with E, cN, cN_bell);
    split: optional list of (label, linestyle, result dict) whose in-degree and scattered
    parts get a middle panel (color = part, line style = case)"""
    dx = np.sqrt(4 * np.pi / n)
    xmax = min(nmax, 2 * np.pi / (ppw_min * dx))
    nplot = min(nmax, int(np.ceil(xmax)))
    N = np.arange(nplot + 1)
    if split:
        fig, axs = plt.subplots(1, 3, figsize=(19, 4.6), gridspec_kw={'width_ratios': [3, 3, 2]})
    else:
        fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={'width_ratios': [3, 2]})
    ax_e, ax_s = axs[0], axs[-1]
    for lab, col, ls, r in curves:
        ax_e.semilogy(N, r['E'][:nplot + 1], ls, color=col, lw=1.6, label=lab)
    r = curves[0][3]
    ax_s.semilogy(N, r['cN'][:nplot + 1], color=INK, lw=1.4, label='Gaussian')
    ax_s.semilogy(N, r['cN_bell'][:nplot + 1], color=INK, lw=1.4, ls='--', label='cosine bell')
    fig.suptitle(title, color=INK, fontsize=10)
    style_ax(ax_e, 'worst-case error after one revolution: E(N) = ||(exp(AT) - I) Q_N||')
    style_ax(ax_s, 'degree spectrum of the swtc1 initial conditions')
    emin = min(c[3]['E'][1:nplot + 1].min() for c in curves)
    ax_e.set_ylim(10 ** np.floor(np.log10(emin / 2)), 3)
    ax_s.set_ylim(1e-8, 1)
    degree_axis(ax_e, dx, xmax)
    degree_axis(ax_s, dx, xmax, bound=False)
    if split:
        ax_m = axs[1]
        for lab, ls, r in split:
            ax_m.semilogy(N, r['in'][:nplot + 1], ls, color=SCAT_COL['in'], lw=1.6, label=f'in-degree, {lab}')
            ax_m.semilogy(N, r['out'][:nplot + 1], ls, color=SCAT_COL['out'], lw=1.6, label=f'scattered, {lab}')
        style_ax(ax_m, 'split form: worst-case error within degree N vs scattered to other degrees')
        imin = min(r['in'][1:nplot + 1].min() for _, _, r in split)
        ax_m.set_ylim(10 ** np.floor(np.log10(imin / 2)), 3)
        degree_axis(ax_m, dx, xmax)
    for ax in axs:
        ax.legend(fontsize=7, loc='lower right')
    ax_e.set_ylabel('relative l2 error at day 12')
    ax_s.set_ylabel('||p0_N|| / ||p0||')
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)


def load_curve(args, sc, f):
    """saved results of one case (nu from args), with the cosine bell spectrum computed for older files"""
    t = f'ne{args.ne}-sc{sc}-{nu_tag(args)}' + (f'-rk4sub{args.sub}' if args.rk4 else '-sd')
    r = np.load(os.path.join(args.out, f'{name(args.ne, sc, f, args.velocity)}-{t}-shnorm.npz'))
    cv = {k: r[k] for k in ('E', 'rms', 'in', 'out', 'cN')}
    if 'cN_bell' in r:
        cv['cN_bell'] = r['cN_bell']
    else:   # older files: compute the cosine bell spectrum from the node coordinates
        d = load(os.path.join(args.mats, name(args.ne, sc, f, args.velocity)), ops=())
        Y, deg = real_harmonics(d['lat'], d['lon'], d['meta']['alpha'], len(r['E']) - 1)
        sq = np.sqrt(d['mass'])
        cv['cN_bell'] = degree_spectrum(cosine_bell(d['lat'], d['lon']), degree_bases(Y, deg, sq), sq)
    return cv, int(r['n'])


def nu_tag(args):
    if args.nu30 == 0:
        return 'nu0'
    return f'nu30_{args.nu30:g}' + ('' if args.nu_exp == 4 else f'_dx{args.nu_exp:g}')


def nu_label(args):
    if args.nu30 == 0:
        return 'nu=0'
    return f'nu={args.nu30 * (30 / args.ne) ** args.nu_exp:.3g} (nu30={args.nu30:g}, ~dx^{args.nu_exp:g})'


def sub_label(args):
    return 'semi-discrete exp(AT)' if not args.rk4 else f'RK4, dt/{args.sub}'


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mats', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--ne', type=int, default=10)
    p.add_argument('--sc', type=int, default=1)
    p.add_argument('--forms', type=int, nargs='+', default=[2, 0])
    p.add_argument('--velocity', default='analytic')
    p.add_argument('--nmax', type=int, default=60)
    p.add_argument('--nu30', type=float, default=0.0,
                   help='hyperviscosity at NE30; nu = nu30 (30/NE)^nu_exp')
    p.add_argument('--nu-exp', type=float, default=4.0,
                   help='resolution scaling exponent of nu (4: nu ~ dx^4, 3: nu ~ dx^3)')
    p.add_argument('--rk4', action='store_true')
    p.add_argument('--sub', type=int, default=1)
    p.add_argument('--dt30', type=float, default=1800.0)
    p.add_argument('--dense', action='store_true')
    p.add_argument('--ppw-min', type=float, default=6.0,
                   help='plot degrees with at least this many mean GLL points per wavelength '
                        '(beyond it interpolated harmonics alias onto lower degrees)')
    p.add_argument('--replot', action='store_true', help='only redo the plot from the saved npz files')
    p.add_argument('--compare-alpha', action='store_true',
                   help='with --replot: overlay alpha=0 (sc1) and alpha=pi/4 (sc4) for each of --forms')
    p.add_argument('--compare-nu', type=float, nargs='+', metavar='NU30',
                   help='with --replot: E(N) for these nu30 values (0 allowed, scaled with --nu-exp) and '
                        'each of --forms on one axis; color = form, line style = nu')
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rep, res = [], {}
    tag = f'ne{args.ne}-sc{args.sc}-{nu_tag(args)}' + (f'-rk4sub{args.sub}' if args.rk4 else '-sd')
    if args.replot and args.compare_nu:
        curves, split = [], []
        for f in args.forms:
            for nu30, ls in zip(args.compare_nu, ('-', '--', ':', '-.')):
                a = argparse.Namespace(**vars(args))
                a.nu30 = nu30
                try:
                    cv, n = load_curve(a, args.sc, f)
                except FileNotFoundError as e:     # run not finished yet: plot the rest
                    print(f'skipping {FORM[f]}, {nu_label(a)}: {e.filename} not found')
                    continue
                curves.append((f'{FORM[f]}, {nu_label(a)}', FCOL[f], ls, cv))
                if f == 2:      # the in-degree / scattered panel shows the split form only
                    split.append((nu_label(a).split(' (')[0], ls, cv))
        forms = ''.join(map(str, args.forms))
        nus = '_'.join(f'{v:g}' for v in args.compare_nu)
        dxs = '' if args.nu_exp == 4 else f'_dx{args.nu_exp:g}'
        out = os.path.join(args.out, f'shnorm_ne{args.ne}-sc{args.sc}-f{forms}-nu30_{nus}{dxs}'
                                     + (f'-rk4sub{args.sub}' if args.rk4 else '-sd') + '.png')
        plot_compare(curves, f'NE{args.ne}, alpha={ALPHA[args.sc]}, {sub_label(args)}',
                     n, args.nmax, args.ppw_min, out, split=split)
        print(out)
        return
    if args.replot:
        scs = (1, 4) if args.compare_alpha else (args.sc,)
        curves = []
        for f in args.forms:
            for sc in scs:
                if args.compare_alpha:
                    lab = f'alpha={ALPHA[sc]}' + (f', {FORM[f]}' if len(args.forms) > 1 else '')
                    col = SCCOL[sc] if len(args.forms) == 1 else (SCCOL[sc] if f == 2 else FCOL[0])
                else:
                    lab, col = FORM[f], FCOL[f]
                cv, n = load_curve(args, sc, f)
                curves.append((lab, col, cv))
        if args.compare_alpha:
            forms = ', '.join(FORM[f] for f in args.forms)
            title = f'NE{args.ne}, {forms}, {nu_label(args)}, {sub_label(args)}'
            out = os.path.join(args.out, f'shnorm_ne{args.ne}-alpha-f{"".join(map(str, args.forms))}'
                                         f'-{nu_tag(args)}' + (f'-rk4sub{args.sub}' if args.rk4 else '-sd') + '.png')
        else:
            title = f'NE{args.ne}, alpha={ALPHA[args.sc]}, {nu_label(args)}, {sub_label(args)}'
            out = os.path.join(args.out, f'shnorm_{tag}.png')
        plot(curves, title, n, args.nmax, args.ppw_min, out)
        return
    for f in args.forms:
        nm = name(args.ne, args.sc, f, args.velocity)
        d = load(os.path.join(args.mats, nm), ops=('L',))
        A = d['L']
        if args.nu30 > 0:
            B = load(os.path.join(args.mats, name(args.ne, 1, 2, 'analytic')), ops=('B',))['B']
            A = (A - args.nu30 * (30 / args.ne) ** args.nu_exp * B).tocsr()
        om = omega(d['meta'])
        T = 2 * np.pi / om
        mass, sq = d['mass'], np.sqrt(d['mass'])
        dt = nstep = None
        if args.rk4:
            dt, nstep = timestep(args.ne, args.dt30, T)
            dt, nstep = dt / args.sub, nstep * args.sub
        Y, deg = real_harmonics(d['lat'], d['lon'], d['meta']['alpha'], args.nmax)
        p0 = gaussian(d['lat'], d['lon'])
        t0 = time.time()
        PP = propagate(A, np.concatenate([Y, p0[:, None]], axis=1), T, dt, nstep)
        P, pT = PP[:, :-1], PP[:, -1]
        tprop = time.time() - t0
        r, Qs = degree_norms(P, Y, deg, sq)
        # Gaussian: degree spectrum, bound and actual error
        q0 = sq * p0
        cN = np.array([np.linalg.norm(Q.T @ q0) for Q in Qs])
        resid = np.linalg.norm(q0 - sum(Q @ (Q.T @ q0) for Q in Qs)) / np.linalg.norm(q0)
        bound = (r['E'] * cN).sum() / np.linalg.norm(q0)
        actual = np.linalg.norm(sq * (pT - p0)) / np.linalg.norm(q0)
        line = (f'{nm} {nu_label(args)} {"RK4 dt=%.0f" % dt if dt else "semi-discrete"}: n={d["n"]}, '
                f'degrees 0..{args.nmax} ({len(deg)} columns, {tprop:.0f} s); Gaussian: actual l2 error '
                f'{actual:.4e}, bound sum E(N)||p0_N|| = {bound:.4e} (Gaussian content beyond degree '
                f'{args.nmax}: {resid:.1e})')
        rep.append(line)
        print(line, flush=True)
        if args.dense:
            t0 = time.time()
            S = (np.diag(sq) @ A.toarray() @ np.diag(1 / sq)) * T
            Ed = la.expm(S)
            Pd = (Ed @ (sq[:, None] * Y)) / sq[:, None]
            rd, _ = degree_norms(Pd, Y, deg, sq)
            lam = np.linalg.eigvals(S / T) / om
            rep.append(f'    dense check ({time.time() - t0:.0f} s): max |E(N) - E_dense(N)| / E_dense(N) = '
                       f'{np.max(np.abs(r["E"] - rd["E"]) / rd["E"]):.1e}; ||exp(A T)||_2 = {la.norm(Ed, 2):.4f}, '
                       f'||exp(A T) - I||_2 = {la.norm(Ed - np.eye(len(sq)), 2):.4f}, max Re(lambda)/omega = '
                       f'{lam.real.max():.3e}')
            print(rep[-1], flush=True)
        cN_bell = degree_spectrum(cosine_bell(d['lat'], d['lon']), Qs, sq)
        res[f] = dict(r, cN=cN / np.linalg.norm(q0), cN_bell=cN_bell)
        np.savez(os.path.join(args.out, f'{nm}-{tag}-shnorm.npz'), N=np.arange(args.nmax + 1), n=d['n'], **r,
                 cN=cN / np.linalg.norm(q0), cN_bell=cN_bell, actual=actual, bound=bound)
        rep.append('    N  E(N)        RMS         in-degree   scattered   Gaussian ||p0_N||/||p0||')
        for N in range(0, args.nmax + 1, 5 if args.nmax > 20 else 1):
            rep.append(f'   {N:2d}  {r["E"][N]:.3e}   {r["rms"][N]:.3e}   {r["in"][N]:.3e}   {r["out"][N]:.3e}   '
                       f'{cN[N] / np.linalg.norm(q0):.3e}')

    plot([(FORM[f], FCOL[f], res[f]) for f in args.forms],
         f'NE{args.ne}, alpha={ALPHA[args.sc]}, {nu_label(args)}, {sub_label(args)}', d['n'], args.nmax, args.ppw_min,
         os.path.join(args.out, f'shnorm_{tag}.png'))
    txt = '\n'.join(rep)
    print(txt)
    with open(os.path.join(args.out, f'shnorm_{tag}.txt'), 'w') as fo:
        fo.write(txt + '\n')


if __name__ == '__main__':
    main()
