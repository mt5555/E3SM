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
from opmat_proj import gaussian, timestep, rk4_run

FORM = {0: 'conservation form', 2: 'split form'}
FCOL = {2: '#eb6834', 0: '#6250d6'}
ALPHA = {1: '0', 4: 'pi/4'}

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


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mats', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--ne', type=int, default=10)
    p.add_argument('--sc', type=int, default=1)
    p.add_argument('--forms', type=int, nargs='+', default=[2, 0])
    p.add_argument('--velocity', default='analytic')
    p.add_argument('--nmax', type=int, default=60)
    p.add_argument('--nu30', type=float, default=0.0)
    p.add_argument('--rk4', action='store_true')
    p.add_argument('--sub', type=int, default=1)
    p.add_argument('--dt30', type=float, default=1800.0)
    p.add_argument('--dense', action='store_true')
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rep, res = [], {}
    tag = f'ne{args.ne}-sc{args.sc}-nu{args.nu30:g}' + (f'-rk4sub{args.sub}' if args.rk4 else '-sd')
    for f in args.forms:
        nm = name(args.ne, args.sc, f, args.velocity)
        d = load(os.path.join(args.mats, nm), ops=('L',))
        A = d['L']
        if args.nu30 > 0:
            B = load(os.path.join(args.mats, name(args.ne, 1, 2, 'analytic')), ops=('B',))['B']
            A = (A - args.nu30 * (30 / args.ne) ** 4 * B).tocsr()
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
        line = (f'{nm} nu30={args.nu30:g} {"RK4 dt=%.0f" % dt if dt else "semi-discrete"}: n={d["n"]}, '
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
        res[f] = dict(r, cN=cN / np.linalg.norm(q0))
        np.savez(os.path.join(args.out, f'{nm}-{tag}-shnorm.npz'), N=np.arange(args.nmax + 1), n=d['n'], **r,
                 cN=cN / np.linalg.norm(q0), actual=actual, bound=bound)
        rep.append('    N  E(N)        RMS         in-degree   scattered   Gaussian ||p0_N||/||p0||')
        for N in range(0, args.nmax + 1, 5 if args.nmax > 20 else 1):
            rep.append(f'   {N:2d}  {r["E"][N]:.3e}   {r["rms"][N]:.3e}   {r["in"][N]:.3e}   {r["out"][N]:.3e}   '
                       f'{cN[N] / np.linalg.norm(q0):.3e}')

    # plot
    n = d['n']
    dx = np.sqrt(4 * np.pi / n)               # mean node spacing on the unit sphere
    N = np.arange(args.nmax + 1)
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6))
    for f in args.forms:
        r = res[f]
        axs[0].semilogy(N[1:], r['E'][1:], '-', color=FCOL[f], lw=1.6, label=f'{FORM[f]}: worst case E(N)')
        axs[0].semilogy(N[1:], r['rms'][1:], ':', color=FCOL[f], lw=1.4, label=f'{FORM[f]}: RMS over degree N')
        axs[1].semilogy(N[1:], r['in'][1:], '-', color=FCOL[f], lw=1.6, label=f'{FORM[f]}: in-degree')
        axs[1].semilogy(N[1:], r['out'][1:], '--', color=FCOL[f], lw=1.4, label=f'{FORM[f]}: scattered')
    axs[0].semilogy(N[1:], res[args.forms[0]]['cN'][1:], color=MUTED, lw=1, label='Gaussian ||p0_N||/||p0||')
    sub = 'semi-discrete exp(AT)' if not args.rk4 else f'RK4, dt/{args.sub}'
    style_ax(axs[0], f'NE{args.ne}, alpha={ALPHA[args.sc]}, nu30={args.nu30:g}, {sub}: '
                     f'||(exp(AT)-I) Q_N|| per degree')
    style_ax(axs[1], 'worst-case error split: within degree N vs scattered to other degrees')
    for ax in axs:
        ax.set_xlabel('spherical harmonic degree N')
        ax.set_ylim(1e-12, 3)
        ax.axhline(2, color=MUTED, lw=0.8, ls='--')
        top = ax.secondary_xaxis('top', functions=(lambda x: 2 * np.pi / (np.maximum(x, 1e-9) * dx),
                                                    lambda y: 2 * np.pi / (np.maximum(y, 1e-9) * dx)))
        top.set_xlabel('mean GLL points per wavelength', color=MUTED)
        top.set_xticks([4, 5, 6, 8, 10, 15, 20, 40])
        ax.legend(fontsize=7, loc='lower right')
    axs[0].set_ylabel('relative l2 error at day 12')
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, f'shnorm_{tag}.png'))
    plt.close(fig)
    txt = '\n'.join(rep)
    print(txt)
    with open(os.path.join(args.out, f'shnorm_{tag}.txt'), 'w') as fo:
        fo.write(txt + '\n')


if __name__ == '__main__':
    main()
