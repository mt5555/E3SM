#!/usr/bin/env python3
"""
Project the swtc1 Gaussian onto the eigenvectors of the extracted transport
operator L (dp/dt = L p, see opmat_eigs.py) and predict its day-12 error.

usage:
  opmat_proj.py --mats DIR --out OUT --ne 8 15 --sc 1 4 [--form 2] [--velocity analytic]
                [--dt30 1800] [--eig 8 15]

For each case the one-revolution (T = 2 pi/omega) error of the Gaussian is
computed three ways:

  fd     RK4 with dt = dt30*30/NE (7200 s at NE8, as in the convergence study):
         R4(dt L)^N p0 - p0 with sparse products.  Reproduces the HOMME run.
  sd     semi-discrete, exp(L T) p0 - p0 (scipy expm_multiply)
  modal  with the eigenpairs L v_j = lambda_j v_j (NE in --eig):
         p0 = sum_j c_j v_j, e = sum_j c_j g_j v_j, with
         g_j = exp(lambda_j T) - 1 (sd) or R4(dt lambda_j)^N - 1 (fd).
         The exact solution is p0 again at T, so exp(-i m omega T) = 1 for every m
         and no wavenumber assignment is needed.

Errors are those of HOMME's tc1_errors: relative l2 with the assembled GLL mass
(= global_integral) and relative max.  Per mode, a_j = |c_j| ||v_j||_M is the
Gaussian's amplitude and E_j = |g_j| a_j its error contribution; sum E_j^2 equals
||e||_M^2 when the eigenvectors are M-orthogonal (split form: to ~1e-6).

Dense eigenvectors are computed one symmetry block at a time (opmat_sym.py);
NE20 alpha=pi/4 (two blocks of 10801) needs about 6 GB.  Use
OPENBLAS_NUM_THREADS=1 when the cores are shared.
"""
import argparse
import os
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from opmat_io import load
from opmat_sym import symmetry_blocks, xyz
from opmat_eigs import name, omega, rk4_amp, rotation_axis, cube_distances

RR, LATC, LONC = 1.0 / 3.0, 0.0, 1.5 * np.pi     # tc1_geopotential, level 2 (Gaussian)


def gaussian(lat, lon):
    A = np.sin(LATC) * np.sin(lat) + np.cos(LATC) * np.cos(lat) * np.cos(lon - LONC)
    return 1000.0 * np.exp(-(np.arccos(np.clip(A, -1, 1)) / RR) ** 2)


def norms(e, p0, mass):
    return (np.sqrt((mass * np.abs(e) ** 2).sum() / (mass * p0 ** 2).sum()),
            np.abs(e).max() / np.abs(p0).max())


def rk4_run(L, p, dt, nstep):
    for _ in range(nstep):
        k1 = L @ p
        k2 = L @ (p + 0.5 * dt * k1)
        k3 = L @ (p + 0.5 * dt * k2)
        k4 = L @ (p + dt * k3)
        p = p + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return p


def timestep(ne, dt30, T):
    """(dt, nstep): dt30*30/NE rounded to a whole number of steps per revolution
    (NE8: 7200 s, as in the convergence study)"""
    nstep = 144 if ne == 8 else int(round(T / (dt30 * 30 / ne)))
    return T / nstep, nstep


def mode_diagnostics(Vc, d, alpha, ne, mass, B):
    """per column of Vc: M-norm, dominant wavenumber m* about the rotation axis and its
    power fraction, corner/edge energy enrichment, B Rayleigh quotient"""
    X = xyz(d['lat'], d['lon'])
    _, e1, e2 = rotation_axis(alpha)
    phi = np.arctan2(X @ e2, X @ e1)
    ms = np.arange(-12 * ne, 12 * ne + 1)
    E = np.exp(-1j * np.outer(ms, phi)) * mass
    C = np.abs(E @ Vc) ** 2
    imax = C.argmax(axis=0)
    w = mass[:, None] * np.abs(Vc) ** 2
    s2 = w.sum(axis=0)
    w /= s2
    h = (np.pi / 2) / ne
    dc, de = cube_distances(d['lat'], d['lon'])
    area = mass / mass.sum()
    cm, em = dc < 1.5 * h, de < 0.5 * h
    r = dict(snorm=np.sqrt(s2), mstar=ms[imax], purity=C[imax, np.arange(Vc.shape[1])] / C.sum(axis=0),
             corner=w[cm].sum(axis=0) / area[cm].sum(), edge=w[em].sum(axis=0) / area[em].sum())
    if B is not None:
        r['q'] = (np.conj(Vc) * (mass[:, None] * (B @ Vc))).sum(axis=0).real / s2
    return r


def skew_part(d):
    """K = (S - S^T)/2, S = M^1/2 L M^-1/2"""
    sq = np.sqrt(d['mass'])
    S = (sp.diags(sq) @ d['L'] @ sp.diags(1 / sq)).tocsr()
    return ((S - S.T) / 2).tocsr(), sq


GOOD = 0.1      # |exp(lambda T) - 1| < GOOD: "resolved" mode (phase error < 0.016 revolutions)


def modal(d, p0, T, dt, nstep, B, skew, chunk=1000, ntop=8):
    """eigen-decomposition of p0 and of the day-T error, block by block, for
    S = M^1/2 L M^-1/2 (nearly normal for the split form, so its eigenvectors are
    nearly orthonormal and the Euclidean norm is the M-norm)"""
    mass, alpha, ne = d['mass'], d['meta']['alpha'], d['meta']['ne']
    sq = np.sqrt(mass)
    S = (sp.diags(sq) @ d['L'] @ sp.diags(1 / sq)).tocsr()
    q0 = sq * p0
    _, blocks = symmetry_blocks(d['lat'], d['lon'], alpha)
    m = len(blocks)
    out = {k: [] for k in ('lam', 'c', 'k', 'snorm', 'mstar', 'purity', 'corner', 'edge', 'q')}
    e_sd = np.zeros(d['n'], complex)
    e_fd = np.zeros(d['n'], complex)
    e_res = np.zeros(d['n'], complex)        # semi-discrete error from resolved modes
    top = []                                 # (E_j, lambda_j, v_j) of the largest error contributions
    gram = 0.0
    prev = None
    for k, Q in blocks:
        t0 = time.time()
        if m == 4 and k == 3:                       # block 3 = conj(block 1) for real L, p0
            lam, U = np.conj(prev[0]), np.conj(prev[1])
        else:
            Ak = (Q.conj().T @ (S @ Q)).toarray()
            if skew:
                # skew-Hermitian part: i K_k is Hermitian, orthonormal eigenvectors
                w, U = np.linalg.eigh(0.5j * (Ak - Ak.conj().T))
                lam = -1j * w
            else:
                lam, U = np.linalg.eig(Ak)
            del Ak
        # departure from orthonormality of the eigenvectors (unit 2-norm columns)
        G = U.conj().T @ U
        gram = max(gram, np.abs(G - np.eye(len(lam))).max())
        del G
        c = np.linalg.solve(U, Q.conj().T @ q0)
        g_sd = np.expm1(lam * T)
        g_fd = rk4_amp(dt * lam) ** nstep - 1
        for j0 in range(0, len(lam), chunk):
            sl = slice(j0, j0 + chunk)
            Vc = Q @ U[:, sl]
            e_sd += Vc @ (c[sl] * g_sd[sl])
            e_fd += Vc @ (c[sl] * g_fd[sl])
            good = np.abs(g_sd[sl]) < GOOD
            e_res += Vc[:, good] @ (c[sl][good] * g_sd[sl][good])
            Ej = np.abs(c[sl] * g_sd[sl])
            for j in np.argsort(-Ej)[:ntop]:
                top.append((Ej[j], lam[sl][j], Vc[:, j] / sq))
            top = sorted(top, key=lambda x: -x[0])[:ntop]
            r = mode_diagnostics(Vc / sq[:, None], d, alpha, ne, mass, B)
            for key in r:
                out[key].append(r[key])
        out['lam'].append(lam); out['c'].append(c); out['k'].append(np.full(len(lam), k))
        if m == 4 and k == 1:
            prev = (lam, U)
        del U
        print(f'    block {k}: size {Q.shape[1]}, {time.time() - t0:.0f} s', flush=True)
    r = {key: np.concatenate(v) for key, v in out.items() if v}
    r['e_sd'], r['e_fd'], r['e_res'] = e_sd / sq, e_fd / sq, e_res / sq
    r['top_E'] = np.array([x[0] for x in top]) / np.sqrt((q0 ** 2).sum())
    r['top_lam'] = np.array([x[1] for x in top])
    r['top_v'] = np.stack([x[2] for x in top], axis=1)
    r['g_sd'] = np.expm1(r['lam'] * T)
    r['g_fd'] = rk4_amp(dt * r['lam']) ** nstep - 1
    r['gram'] = gram
    return r


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mats', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--ne', type=int, nargs='+', default=[8, 15])
    p.add_argument('--sc', type=int, nargs='+', default=[1, 4])
    p.add_argument('--form', type=int, default=2)
    p.add_argument('--velocity', default='analytic')
    p.add_argument('--dt30', type=float, default=1800.0)
    p.add_argument('--eig', type=int, nargs='*', default=[8, 15])
    p.add_argument('--full', action='store_true',
                   help='eigenvectors of S itself (default: of its skew part K, orthonormal)')
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rep = []
    for ne in args.ne:
        dB = None
        bfile = os.path.join(args.mats, name(ne, 1, 2, 'analytic'))
        if os.path.exists(bfile + '.meta'):
            dB = load(bfile, ops=('B',)).get('B')
        for sc in args.sc:
            nm = name(ne, sc, args.form, args.velocity)
            d = load(os.path.join(args.mats, nm), ops=('L',))
            om = omega(d['meta'])
            T = 2 * np.pi / om
            dt, nstep = timestep(ne, args.dt30, T)
            p0 = gaussian(d['lat'], d['lon'])
            t0 = time.time()
            e_fd = rk4_run(d['L'], p0, dt, nstep) - p0
            e_sd = spla.expm_multiply(d['L'] * T, p0) - p0
            l2fd, lifd = norms(e_fd, p0, d['mass'])
            l2sd, lisd = norms(e_sd, p0, d['mass'])
            rep.append(f'{nm} dt={dt:.0f} nstep={nstep}: RK4 l2={l2fd:.6e} max={lifd:.6e}; '
                       f'semi-discrete l2={l2sd:.6e} max={lisd:.6e}; '
                       f'||e_fd-e_sd||/||e_fd||={norms(e_fd - e_sd, e_fd, d["mass"])[0]:.3e} '
                       f'({time.time() - t0:.0f} s)')
            print(rep[-1], flush=True)
            K, sq = skew_part(d)
            e_sdK = spla.expm_multiply(K * T, sq * p0) / sq - p0
            e_fdK = rk4_run(K, sq * p0, dt, nstep) / sq - p0
            rep.append(f'    skew part K only: RK4 l2={norms(e_fdK, p0, d["mass"])[0]:.6e} max={norms(e_fdK, p0, d["mass"])[1]:.6e}; '
                       f'semi-discrete l2={norms(e_sdK, p0, d["mass"])[0]:.6e}; ||e_fd(K)-e_fd(L)||/||e_fd(L)||='
                       f'{norms(e_fdK - e_fd, e_fd, d["mass"])[0]:.3e}')
            print(rep[-1], flush=True)
            res = dict(e_fd_direct=e_fd, e_sd_direct=e_sd, p0=p0, mass=d['mass'], lat=d['lat'], lon=d['lon'],
                       T=T, dt=dt, nstep=nstep, omega=om)
            if ne in args.eig:
                t0 = time.time()
                r = modal(d, p0, T, dt, nstep, dB, skew=not args.full)
                res.update(r)
                es, ef = (e_sd, e_fd) if args.full else (e_sdK, e_fdK)
                res['skew'] = not args.full
                a = np.abs(r['c']) * r['snorm']
                Esd, Efd = np.abs(r['g_sd']) * a, np.abs(r['g_fd']) * a
                nrm = np.sqrt((d['mass'] * p0 ** 2).sum())
                rep.append(f'    modal ({time.time() - t0:.0f} s): max|U^H U - I|={r["gram"]:.1e}, reconstruction err sd {norms(r["e_sd"] - es, es, d["mass"])[0]:.2e}, '
                           f'fd {norms(r["e_fd"] - ef, ef, d["mass"])[0]:.2e}; sqrt(sum a^2)/||p0||='
                           f'{np.sqrt((a**2).sum()) / nrm:.6f}; sqrt(sum E_sd^2)/||p0||={np.sqrt((Esd**2).sum()) / nrm:.6e}, '
                           f'sum E_fd^2: {np.sqrt((Efd**2).sum()) / nrm:.6e}')
                print(rep[-1], flush=True)
            np.savez(os.path.join(args.out, nm + '-proj.npz'), **res)
    with open(os.path.join(args.out, 'proj.txt'), 'a') as f:
        f.write('\n'.join(rep) + '\n')


if __name__ == '__main__':
    main()
