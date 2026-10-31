#!/usr/bin/env python3
"""
Eigenvalue analysis of the sweqx swtc1 transport operator L (dp/dt = L p) and
the hyperviscosity operator B (approximately del^4), extracted with
SWEQX_OPMAT (src/sweqx/opmat_mod.F90, ../swtc1-opmat.sh).

usage:
  opmat_eigs.py --mats DIR [--mats30 DIR30] --out OUT [--ne 8 15]

DIR holds ne{NE}-sc{1,4}-f{0,2}-{analytic,curlpsi}.* for each NE (B in the
ne{NE}-sc1-f2-analytic files), DIR30 the NE30 files (only extremal eigenvalues
are computed there, with ARPACK).  Results are cached in OUT/cache, plots and
results.txt are written to OUT.

The exact operator -u.grad for solid body rotation with angular speed
omega = u0/rearth has eigenvalues -i m omega, m integer (e^{i m phi'}, phi' the
longitude about the rotation axis), for any alpha.  Eigenvalues are reported in
units of omega.  Dense eigenvalue problems are reduced with the symmetries
of the grid and the flow (a rotation and the antipodal map, opmat_sym.py) and
solved one block at a time.

NE15 takes ~10 minutes on 4 cores (all cases, B, and L - nu B for 4 nu and
2 cases).  If the cores are shared with other jobs, set OPENBLAS_NUM_THREADS=1:
an oversubscribed multithreaded OpenBLAS made the dense eigensolver ~15x slower.
"""
import argparse
import os
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from opmat_io import load
from opmat_sym import block_eigs, eigvals_by_symmetry, symmetry_blocks, xyz

CASES = [(sc, f, v) for sc in (1, 4) for f in (2, 0) for v in ('analytic', 'curlpsi')]
FORM = {0: 'conservation', 2: 'split'}
ALPHA = {1: '0', 4: 'pi/4'}
# fixed categorical order (one color per case family)
COLORS = {('f2', 'curlpsi'): '#2a78d6', ('f2', 'analytic'): '#eb6834',
          ('f0', 'curlpsi'): '#1baf7a', ('f0', 'analytic'): '#6250d6'}
INK, MUTED, GRID = '#0b0b0b', '#898781', '#e1e0d9'
NU30 = [1e13, 5.4e13, 1e14]
SQRT8 = 2 * np.sqrt(2)

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': 0.6, 'legend.frameon': False,
                     'figure.dpi': 130})


def name(ne, sc, f, v):
    return f'ne{ne}-sc{sc}-f{f}-{v}'


def label(sc, f, v):
    return f'alpha={ALPHA[sc]}, {FORM[f]}, {v}'


def omega(meta):
    return meta['u0'] / meta['rearth']


class Cache:
    def __init__(self, out):
        self.dir = os.path.join(out, 'cache')
        os.makedirs(self.dir, exist_ok=True)

    def get(self, key, fn):
        f = os.path.join(self.dir, key + '.npz')
        if os.path.exists(f):
            return dict(np.load(f, allow_pickle=True))
        t = time.time()
        r = fn()
        np.savez(f, **r)
        print(f'  computed {key} in {time.time() - t:.0f} s', flush=True)
        return r


# ---------------------------------------------------------------- diagnostics
def mass_sym(A, mass):
    """S = M^1/2 A M^-1/2 (sparse)"""
    return sp.diags(np.sqrt(mass)) @ A @ sp.diags(1 / np.sqrt(mass))


def skew_sym_err(A, mass):
    S = mass_sym(A, mass)
    nS = sp.linalg.norm(S)
    return sp.linalg.norm(S + S.T) / nS, sp.linalg.norm(S - S.T) / nS


def rk4_amp(z):
    return 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24


def rk4_dtmax(lam, gtol):
    """largest dt such that the RK4 growth rate max log|R(dt lam)|/dt does not exceed
    max(max Re lam, 0) + gtol, scanning dt upwards (lam: eigenvalues, 1/s)"""
    lam = lam[np.abs(lam) > 0]
    base = max(lam.real.max(), 0.0)
    d0 = SQRT8 / np.abs(lam).max()
    dts = np.linspace(0.3 * d0, 3.0 * d0, 2701)
    prev = None
    for dt in dts:
        g = np.log(np.abs(rk4_amp(dt * lam))).max() / dt
        if g > base + gtol:
            return prev
        prev = dt
    return prev


def rotation_axis(alpha):
    a = np.array([-np.sin(alpha), 0, np.cos(alpha)])
    e2 = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(e2, a)
    return a, e1, e2


def cube_distances(lat, lon):
    """angular distance (rad) of each node to the nearest cube corner and cube edge
    (HOMME orientation: corners at (+-1,+-1,+-1)/sqrt(3))"""
    X = xyz(lat, lon)
    corners = np.array([[i, j, k] for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)]) / np.sqrt(3)
    dc = np.arccos(np.clip(X @ corners.T, -1, 1)).min(axis=1)
    A = np.sort(np.abs(X), axis=1)        # A[:,2] = |largest|, the face direction
    de = np.arcsin(np.clip((A[:, 2] - A[:, 1]) / np.sqrt(2), 0, 1))
    return dc, de


def vector_diagnostics(V, d, alpha, ne, B=None):
    """per eigenvector: dominant rotated wavenumber m*, its power fraction,
    participation fraction PN/n, energy enrichment near cube corners and edges,
    and (with B) the M-Rayleigh quotient q = <v,Bv>/<v,v> (first-order damping rate nu*q)"""
    mass = d['mass']
    X = xyz(d['lat'], d['lon'])
    _, e1, e2 = rotation_axis(alpha)
    phi = np.arctan2(X @ e2, X @ e1)
    mmax = 4 * ne * 3
    ms = np.arange(-mmax, mmax + 1)
    E = np.exp(-1j * np.outer(ms, phi)) * mass
    C = np.abs(E @ V) ** 2
    imax = C.argmax(axis=0)
    mstar = ms[imax]
    purity = C[imax, np.arange(V.shape[1])] / C.sum(axis=0)
    w = mass[:, None] * np.abs(V) ** 2
    w /= w.sum(axis=0)
    pn = 1 / (w**2).sum(axis=0) / len(mass)
    h = (np.pi / 2) / ne
    dc, de = cube_distances(d['lat'], d['lon'])
    area = mass / mass.sum()
    cm, em = dc < 1.5 * h, de < 0.5 * h
    corner = w[cm].sum(axis=0) / area[cm].sum()
    edge = w[em].sum(axis=0) / area[em].sum()
    r = dict(mstar=mstar, purity=purity, pn=pn, corner=corner, edge=edge)
    if B is not None:
        r['q'] = (np.conj(V) * (mass[:, None] * (B @ V))).sum(axis=0).real / (mass[:, None] * np.abs(V)**2).sum(axis=0)
    return r


# ---------------------------------------------------------------- per-case solves
def solve_L(cache, d, key, vectors, B=None, chunk=1000):
    """eigenvalues of L block by block; with vectors, the per-mode diagnostics are
    computed in chunks of eigenvectors, so the n x n eigenvector matrix is never formed"""
    def fn():
        ne = d['meta']['ne']
        alpha = d['meta']['alpha']
        lams, ks, diag, cand = [], [], {}, []
        off = 0.0
        for lab, lam, Q, U, info in block_eigs(d['L'], d['lat'], d['lon'], alpha, vectors=vectors):
            lams.append(lam); ks.append(np.full(len(lam), lab))
            off = max(off, info['offblock'])
            if not vectors:
                continue
            for j0 in range(0, len(lam), chunk):
                V = Q @ U[:, j0:j0 + chunk]
                for kk, v in vector_diagnostics(V, d, alpha, ne, B).items():
                    diag.setdefault(kk, []).append(v)
                # candidates for the saved maps (largest Re, largest |lam| in this chunk)
                lc = lam[j0:j0 + chunk]
                for j in np.unique(np.r_[np.argsort(-lc.real)[:6], np.argsort(-np.abs(lc))[:2]]):
                    cand.append((lc[j], np.abs(V[:, j]) ** 2))
        lam = np.concatenate(lams)
        r = dict(lam=lam, ks=np.concatenate(ks), commute=info['commute'], offblock=off)
        if vectors:
            r.update({kk: np.concatenate(v) for kk, v in diag.items()})
            # save |v|^2 maps of a few modes: largest Re, largest |lam|
            sel = np.unique(np.r_[np.argsort(-lam.real)[:6], np.argsort(-np.abs(lam))[:2]])
            cl = np.array([c[0] for c in cand])
            r['sel'] = sel
            r['selmaps'] = np.stack([cand[np.argmin(np.abs(cl - lam[s]))][1] for s in sel], axis=1)
        return r
    return cache.get(key, fn)


def solve_sym_B(cache, d, key):
    """eigenvalues of B (M-self-adjoint) with the symmetry blocks of S = M^1/2 B M^-1/2"""
    def fn():
        S = mass_sym(d['B'], d['mass'])
        Ss = ((S + S.T) / 2).tocsr()
        _, blocks = symmetry_blocks(d['lat'], d['lon'], 0.0)
        lam = {}
        for lab, Q, src in blocks:
            if src is not None:          # conjugate block: same (real) eigenvalues
                lam[lab] = lam[src]
                continue
            Bk = (Q.conj().T @ Ss @ Q).toarray()
            lam[lab] = np.linalg.eigvalsh((Bk + Bk.conj().T) / 2)
        return dict(lam=np.sort(np.concatenate(list(lam.values()))))
    return cache.get(key, fn)


def solve_A(cache, d, B, nu, key):
    def fn():
        A = (d['L'] - nu * B).tocsr()
        lam, _, info, ks = eigvals_by_symmetry(A, d['lat'], d['lon'], d['meta']['alpha'])
        return dict(lam=lam, commute=info['commute'], offblock=info['offblock'])
    return cache.get(key, fn)


# ---------------------------------------------------------------- plots
def style_ax(ax, title=None):
    ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    if title:
        ax.set_title(title, color=INK, fontsize=9, loc='left')


def plot_spectra(res, ne, om, out):
    fig, axs = plt.subplots(2, 4, figsize=(13, 6.5), sharey=True)
    for ax, (sc, f, v) in zip(axs.flat, CASES):
        lam = res[(ne, sc, f, v)]['lam'] / om
        ax.scatter(lam.real, lam.imag, s=3, color=COLORS[(f'f{f}', v)], linewidths=0)
        style_ax(ax, label(sc, f, v))
        ax.set_xlabel('Re(lambda)/omega')
    for ax in axs[:, 0]:
        ax.set_ylabel('Im(lambda)/omega')
    fig.suptitle(f'NE{ne}: spectrum of L (exact: -i m, m integer, on the imaginary axis)',
                 x=0.01, ha='left', color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'spectrum_L_ne{ne}.png'))
    plt.close(fig)


def plot_accuracy_nearest_integer(res, ne, om, out):
    """distance of lambda/omega to the nearest exact eigenvalue -i m (a lower bound for
    the error), vs |Im|/omega.  rows: split / conservation form, columns: alpha"""
    fig, axs = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True, sharey=True)
    for row, f in enumerate((2, 0)):
        for ax, sc in zip(axs[row], (1, 4)):
            for v in ('analytic', 'curlpsi'):
                lam = res[(ne, sc, f, v)]['lam'] / om
                dist = np.abs(lam + 1j * np.round(-lam.imag))   # m = nearest integer to -Im
                ax.scatter(np.abs(lam.imag), np.maximum(dist, 1e-15), s=2 if v == 'analytic' else 4,
                           color=COLORS[(f'f{f}', v)], linewidths=0, label=v,
                           alpha=1.0 if v == 'analytic' else 0.6)
            ax.set_yscale('log')
            ax.set_ylim(1e-15, 1)
            style_ax(ax, f'NE{ne}, alpha={ALPHA[sc]}, {FORM[f]} form')
            ax.legend(markerscale=4, loc='lower right')
    for ax in axs[1]:
        ax.set_xlabel('|Im(lambda)|/omega')
    for ax in axs[:, 0]:
        ax.set_ylabel('|lambda/omega + i m|, m nearest integer')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'accuracy_nearest_integer_ne{ne}.png'))
    plt.close(fig)


def plot_accuracy_mstar(res, ne, om, out):
    """eigenvalue error vs the dominant rotated wavenumber m* of the eigenvector"""
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for ax, sc in zip(axs, (1, 4)):
        for f, v in [(2, 'curlpsi'), (2, 'analytic'), (0, 'curlpsi'), (0, 'analytic')]:
            r = res[(ne, sc, f, v)]
            lam = r['lam'] / om
            err = np.abs(lam + 1j * r['mstar'])
            good = r['purity'] > 0.5
            ax.scatter(np.abs(r['mstar'][good]), err[good] + 1e-16, s=3,
                       color=COLORS[(f'f{f}', v)], linewidths=0, label=f'{FORM[f]}, {v}')
        ax.set_yscale('log')
        ax.set_ylim(1e-14, 30)
        ax.set_xlabel('|m*| (dominant wavenumber about the rotation axis)')
        style_ax(ax, f'NE{ne}, alpha={ALPHA[sc]}  (modes with >50% power in m*)')
    axs[0].set_ylabel('|lambda/omega + i m*|')
    axs[1].legend(markerscale=4, loc='lower right')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'accuracy_mstar_ne{ne}.png'))
    plt.close(fig)


def plot_sorted_imag(res, ne, om, out):
    """sorted Im(lambda)/omega >= 0 compared to a spherical harmonic truncation with
    the same number of degrees of freedom (m has multiplicity N+1-|m|)"""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    n = None
    for sc, ls in ((1, '-'), (4, '--')):
        for f, v in [(2, 'curlpsi'), (0, 'analytic')]:
            lam = res[(ne, sc, f, v)]['lam'] / om
            n = len(lam)
            y = np.sort(np.abs(lam.imag))
            ax.plot(np.arange(n) / n, y, ls, lw=1.5, color=COLORS[(f'f{f}', v)],
                    label=f'alpha={ALPHA[sc]}, {FORM[f]}, {v}')
    N = int(np.sqrt(n)) - 1
    ref = np.sort(np.abs(np.concatenate([np.full(N + 1 - abs(m), m) for m in range(-N, N + 1)])))
    ax.plot(np.arange(len(ref)) / len(ref), ref, color=MUTED, lw=1, label=f'SH truncation N={N}')
    ax.set_xlabel('fraction of eigenvalues')
    ax.set_ylabel('|Im(lambda)|/omega, sorted')
    style_ax(ax, f'NE{ne}: sorted |Im(lambda)|/omega')
    ax.legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'sorted_imag_ne{ne}.png'))
    plt.close(fig)


def plot_localization(res, ne, om, out):
    """corner enrichment of each eigenvector; growing modes colored by their growth rate"""
    fig, axs = plt.subplots(2, 4, figsize=(13, 6.5), sharex=True, sharey=True)
    for ax, (sc, f, v) in zip(axs.flat, CASES):
        r = res[(ne, sc, f, v)]
        lam = r['lam'] / om
        ce = np.clip(r['corner'], 1e-2, 1e2)
        g = lam.real > 1e-9
        ax.scatter(np.abs(lam.imag[~g]), ce[~g], s=3, color='#c3c2b7', linewidths=0, label='Re <= 1e-9 omega')
        if g.any():
            o = np.argsort(lam.real[g])
            scat = ax.scatter(np.abs(lam.imag[g][o]), ce[g][o], c=np.log10(lam.real[g][o]), s=6,
                              cmap='viridis', vmin=-9, vmax=0, linewidths=0)
        style_ax(ax, label(sc, f, v) + f'\n{g.sum()} modes with Re > 1e-9 omega')
        ax.set_yscale('log')
        ax.set_ylim(1e-2, 1e2)
    fig.colorbar(scat, ax=axs, label='log10(Re(lambda)/omega) of growing modes', shrink=0.8)
    for ax in axs[1]:
        ax.set_xlabel('|Im(lambda)|/omega')
    for ax in axs[:, 0]:
        ax.set_ylabel('corner enrichment (clipped)')
    fig.suptitle(f'NE{ne}: eigenvector energy within 1.5 elements of the cube corners / area fraction '
                 f'(1 = uniform)', x=0.01, ha='left', color=INK)
    fig.savefig(os.path.join(out, f'localization_corner_ne{ne}.png'))
    plt.close(fig)


def plot_mode_maps(res, d, ne, sc, f, v, om, out):
    r = res[(ne, sc, f, v)]
    lam = r['lam'] / om
    sel, maps = r['sel'], r['selmaps']
    k = min(len(sel), 8)
    fig, axs = plt.subplots(2, 4, figsize=(14, 5.5))
    lon = np.degrees(d['lon']); lon = np.where(lon > 180, lon - 360, lon)
    lat = np.degrees(d['lat'])
    for ax, i in zip(axs.flat, range(k)):
        ax.scatter(lon, lat, c=maps[:, i] / maps[:, i].max(), s=4, cmap='Blues', linewidths=0)
        for clat in (35.26, -35.26):
            ax.scatter([45, 135, -45, -135], [clat] * 4, marker='x', color='#e34948', s=20)
        style_ax(ax, f'lambda/omega={lam[sel[i]].real:.3g}{lam[sel[i]].imag:+.4g}i')
        ax.set_xlim(-180, 180); ax.set_ylim(-90, 90)
        ax.grid(False)
    for ax in axs.flat[k:]:
        ax.axis('off')
    fig.suptitle(f'NE{ne} {label(sc, f, v)}: |v|^2 of the modes with largest Re(lambda) and '
                 f'largest |lambda| (x: cube corners)', x=0.01, ha='left', color=INK)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'modes_{name(ne, sc, f, v)}.png'))
    plt.close(fig)


def plot_B(resB, ne_list, rearth, out):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for ne, c in zip(ne_list, ('#2a78d6', '#eb6834', '#1baf7a')):
        lam = resB[ne]['lam']
        n = len(lam)
        ls = np.arange(int(np.sqrt(n)) + 2)
        ex = np.repeat((ls * (ls + 1)) ** 2 / rearth**4, 2 * ls + 1)[:n]
        ax.plot(np.sqrt(np.arange(1, n)), lam[1:] / ex[1:], lw=1.5, color=c, label=f'NE{ne}')
    ax.set_xlabel('sqrt(index) ~ spherical harmonic degree l')
    ax.set_ylabel('lambda(B) / (l(l+1)/r^2)^2')
    ax.set_yscale('log')
    style_ax(ax, 'Hyperviscosity B: sorted eigenvalues / exact del^4 eigenvalues')
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'spectrum_B.png'))
    plt.close(fig)


def plot_hypervis(resA, ne, cases, nus, om, out):
    """spectra of A = L - nu B in the complex plane (symlog real axis)"""
    fig, axs = plt.subplots(1, len(cases), figsize=(6.5 * len(cases), 4.8), sharey=True, squeeze=False)
    cols = ['#898781', '#2a78d6', '#eb6834', '#6250d6']
    for ax, case in zip(axs[0], cases):
        for nu, c in zip(nus[case], cols):
            lam = resA[(ne,) + case + (nu,)]['lam'] / om
            ax.scatter(lam.real, lam.imag, s=3, color=c, linewidths=0,
                       label='nu=0' if nu == 0 else f'nu={nu:.3g} (nu30={nu / (30 / ne)**4:.3g})')
        ax.set_xscale('symlog', linthresh=1e-2)
        ax.set_xlim(right=1)
        ax.set_xlabel('Re(lambda)/omega  (symlog; -1/(2 pi) = 1/e damping per revolution)')
        style_ax(ax, f'NE{ne} {label(*case)}: spectrum of A = L - nu B')
    axs[0][0].set_ylabel('Im(lambda)/omega')
    axs[0][-1].legend(markerscale=4, loc='center right')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'hypervis_ne{ne}.png'))
    plt.close(fig)


def leff(q, rearth):
    """effective spherical harmonic degree of a mode from its B Rayleigh quotient"""
    return (-1 + np.sqrt(1 + 4 * rearth**2 * np.sqrt(np.maximum(q, 0)))) / 2


def plot_error_vs_leff(res, ne, om, rearth, out):
    """eigenvalue error of the L eigenvectors vs their effective degree l_eff from <v,Bv>"""
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for ax, sc in zip(axs, (1, 4)):
        for f, v in [(2, 'curlpsi'), (2, 'analytic'), (0, 'curlpsi'), (0, 'analytic')]:
            r = res[(ne, sc, f, v)]
            err = np.abs(r['lam'] / om + 1j * r['mstar'])
            good = r['purity'] > 0.5
            ax.scatter(leff(r['q'][good], rearth), err[good] + 1e-16, s=3, color=COLORS[(f'f{f}', v)],
                       linewidths=0, label=f'{FORM[f]}, {v}')
        ax.set_yscale('log')
        ax.set_ylim(1e-14, 30)
        ax.set_xlabel('l_eff: l(l+1)/r^2 = sqrt(<v,Bv>/<v,v>)  (hyperviscosity damping rate nu*(l(l+1)/r^2)^2)')
        style_ax(ax, f'NE{ne}, alpha={ALPHA[sc]}: eigenvalue error vs l_eff (modes with >50% power in m*)')
    axs[0].set_ylabel('|lambda/omega + i m*|')
    axs[1].legend(markerscale=4, loc='lower right')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'error_vs_leff_ne{ne}.png'))
    plt.close(fig)


# ---------------------------------------------------------------- NE30 (sparse)
def ne30_analysis(cache, mats30, rep):
    files = sorted(f[:-5] for f in os.listdir(mats30) if f.startswith('ne30') and f.endswith('.meta'))
    if not files:
        return
    rep.append('\n== NE30 (ARPACK, extremal eigenvalues only) ==')
    dB = None
    for nm in files:
        d = load(os.path.join(mats30, nm))
        om = omega(d['meta'])
        if 'L' in d:
            r = cache.get(nm + '-LM', lambda: dict(lam=spla.eigs(d['L'], k=16, which='LM', tol=1e-10,
                                                                 return_eigenvectors=False)))
            lam = r['lam']
            sk, _ = skew_sym_err(d['L'], d['mass'])
            dtm = rk4_dtmax(lam, 1e-3 * om)
            rep.append(f'{nm}: n={d["n"]} nnz(L)={d["L"].nnz} ||S+S^T||/||S||={sk:.2e} '
                       f'max|lambda|={np.abs(lam).max():.6e}/s ={np.abs(lam).max() / om:.1f} omega, '
                       f'2.83/max|lambda|={SQRT8 / np.abs(lam).max():.1f} s, RK4 dt_max(16 LM eigs)={dtm:.1f} s, '
                       f'Re of LM eigs in [{lam.real.min() / om:.3g},{lam.real.max() / om:.3g}] omega')
        if 'B' in d and dB is None:
            dB = d
    if dB is not None:
        S = mass_sym(dB['B'], dB['mass'])
        Ss = (S + S.T) / 2
        r = cache.get('ne30-B-LA', lambda: dict(lam=spla.eigsh(Ss, k=4, which='LA', tol=1e-10,
                                                               return_eigenvectors=False)))
        lmax = r['lam'].max()
        dB['lmax'] = lmax
        rep.append(f'NE30 B: nnz={dB["B"].nnz} ||S-S^T||/||S||={sp.linalg.norm(S - S.T) / sp.linalg.norm(S):.2e} '
                   f'lambda_max(B)={lmax:.6e} m^-4; forward-Euler hypervis limit dt*nu*lambda_max<=2: '
                   + ', '.join(f'nu={nu:.3g}: dt<={2 / (nu * lmax):.0f} s' for nu in (1e14, 2e14, 1e15)))
        return dB, lmax


def ne30_rk4_hv_growth(cache, mats30, dB, rep):
    """per-step amplification of the HOMME step G = (I - dt nu B) R4(dt L) (RK4, then
    forward-Euler hyperviscosity) for the largest-|lambda| eigenpairs of L at NE30,
    sc4 split analytic (the swtc1.sh configuration), to first order in the
    hyperviscosity:  |R4(dt lambda)| * (1 - dt nu <v,Bv>/<v,v>)"""
    nm = 'ne30-sc4-f2-analytic'
    if not os.path.exists(os.path.join(mats30, nm + '.meta')):
        return
    d = load(os.path.join(mats30, nm), ops=('L',))
    L, B, mass = d['L'], dB['B'], d['mass']

    def fn():
        lam, V = spla.eigs(L, k=150, which='LM', tol=1e-10)
        q = (np.conj(V) * (mass[:, None] * (B @ V))).sum(axis=0).real / (mass[:, None] * np.abs(V)**2).sum(axis=0)
        return dict(lam=lam, q=q)
    r = cache.get(nm + '-LMvec150', fn)
    lam, q = r['lam'], r['q']
    om = omega(d['meta'])
    rep.append(f'{nm}: 150 largest-|lambda| modes (|lambda|/omega in [{np.abs(lam).min() / om:.1f},{np.abs(lam).max() / om:.1f}]): <v,Bv>/<v,v> in [{q.min():.3e},{q.max():.3e}] m^-4 '
               f'(lambda_max(B) fraction {q.max() / dB["lmax"]:.2f})')
    rep.append('  HOMME step amplification max over these modes, |R4(dt lam)|(1-dt nu q), and growth over 36 days:')
    for nu, dt in [(0, 2330), (0, 2340), (1e14, 2330), (1e14, 2580), (1e14, 2590), (1e14, 2700), (2e14, 1800)]:
        g = (np.abs(rk4_amp(dt * lam)) * (1 - dt * nu * q)).max()
        rep.append(f'    nu={nu:g} dt={dt}: max amp-1={g - 1:+.3e}, amp^(36 d/dt)={g ** (36 * 86400 / dt):.3g}')
    for nu in (0, 1e14, 2e14):
        dts = np.arange(1500, 3400, 5.0)
        amp = np.array([(np.abs(rk4_amp(dt * lam)) * (1 - dt * nu * q)).max() for dt in dts])
        with np.errstate(over='ignore'):
            bad = dts[amp ** (36 * 86400 / dts) > 100]
        rep.append(f'    nu={nu:g}: smallest dt with growth >100 over 36 days (these modes, first order): '
                   f'{bad.min() if len(bad) else np.nan:.0f} s')


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mats', required=True)
    p.add_argument('--mats30')
    p.add_argument('--out', required=True)
    p.add_argument('--ne', type=int, nargs='+', default=[8, 15])
    p.add_argument('--vectors-ne', type=int, nargs='*', default=[8],
                   help='NE values for which all eigenvectors are computed')
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    cache = Cache(args.out)
    rep = []
    res, resB, resA = {}, {}, {}
    om = rearth = None

    for ne in args.ne:
        rep.append(f'\n== NE{ne} ==')
        rep.append('case: n, test residual ||L x - y||/||y|| (max over 4 fields), ||L 1||_inf/||L||_inf, '
                   '||S+S^T||/||S|| (S=M^1/2 L M^-1/2), symmetry check, max Re/omega, max|lambda|/omega, '
                   '2.83/max|lambda| [s], RK4 dt_max [s]')
        dref = None
        for sc, f, v in CASES:
            nm = name(ne, sc, f, v)
            path = os.path.join(args.mats, nm)
            if not os.path.exists(path + '.meta'):
                rep.append(f'{nm}: missing')
                continue
            d = load(path)
            om = omega(d['meta']); rearth = d['meta']['rearth']
            if 'B' in d and dref is None:
                dref = d
            x, y = d['Ltest']
            tres = max(np.linalg.norm(d['L'] @ x[:, k] - y[:, k]) / np.linalg.norm(y[:, k]) for k in range(x.shape[1]))
            l1 = np.abs(d['L'] @ np.ones(d['n'])).max() / abs(d['L']).max()
            sk, _ = skew_sym_err(d['L'], d['mass'])
            r = solve_L(cache, d, nm, vectors=ne in args.vectors_ne,
                        B=dref['B'] if dref is not None else None)
            res[(ne, sc, f, v)] = r
            lam = r['lam']
            dtm = rk4_dtmax(lam, 1e-3 * om)
            rep.append(f'{nm}: n={d["n"]} test={tres:.1e} L1={l1:.1e} skew={sk:.2e} '
                       f'sym(comm={float(r["commute"]):.0e},off={float(r["offblock"]):.0e}) '
                       f'maxRe={lam.real.max() / om:.3e} max|lam|={np.abs(lam).max() / om:.2f} '
                       f'2.83/max|lam|={SQRT8 / np.abs(lam).max():.0f} RK4dt={dtm:.0f}')
            if ne in args.vectors_ne:
                g = lam.real > 1e-6 * om
                rep.append(f'    eigvecs: median PN/n={np.median(r["pn"]):.3f}; growing modes '
                           f'(Re>1e-6 omega): {g.sum()}, median PN/n={np.median(r["pn"][g]) if g.any() else np.nan:.3f}, '
                           f'median corner enrichment={np.median(r["corner"][g]) if g.any() else np.nan:.2f} '
                           f'(all modes {np.median(r["corner"]):.2f}), median edge enrichment='
                           f'{np.median(r["edge"][g]) if g.any() else np.nan:.2f} (all {np.median(r["edge"]):.2f})')
                err0 = np.abs(lam / om + 1j * r['mstar'])
                top = np.argsort(-lam.real)[:3]
                rep.append('    top growing: ' + '; '.join(
                    f'{lam[i].real / om:.3g}{lam[i].imag / om:+.4g}i m*={r["mstar"][i]} PN/n={r["pn"][i]:.3f} '
                    f'corner={r["corner"][i]:.2f} edge={r["edge"][i]:.2f}' for i in top))
                if 'q' in r:
                    T = 2 * np.pi / om
                    good = r['purity'] > 0.5
                    for nu30 in NU30:
                        nu = nu30 * (30 / ne) ** 4
                        fac = np.exp(-nu * r['q'] * T)
                        acc, bad = good & (err0 < 1e-2), good & (err0 > 1e-1)
                        rep.append(f'    first-order hypervis nu30={nu30:.2g} (nu={nu:.3g}): per-revolution factor '
                                   f'exp(-nu<v,Bv>T), median over modes with error<1e-2: {np.median(fac[acc]):.3f} '
                                   f'({acc.sum()} modes), error>0.1: {np.median(fac[bad]):.2e} ({bad.sum()} modes)')
                # accuracy summary by |m*|
                err = np.abs(lam / om + 1j * r['mstar'])
                good = r['purity'] > 0.5
                for mm in (4, 8, 16, 32):
                    s = good & (np.abs(r['mstar']) == mm)
                    if s.any():
                        rep.append(f'    |m*|={mm}: {s.sum()} modes, median |lam/omega+i m*|={np.median(err[s]):.2e}, '
                                   f'max {err[s].max():.2e}')
        # B
        if dref is not None:
            S = mass_sym(dref['B'], dref['mass'])
            symerr = sp.linalg.norm(S - S.T) / sp.linalg.norm(S)
            xb, yb = dref['Btest']
            tres = max(np.linalg.norm(dref['B'] @ xb[:, k] - yb[:, k]) / (sp.linalg.norm(dref['B']) * np.linalg.norm(xb[:, k]) / np.sqrt(dref['n'])) for k in range(xb.shape[1]))
            rB = solve_sym_B(cache, dref, f'ne{ne}-B')
            resB[ne] = rB
            lb = rB['lam']
            ls = np.arange(int(np.sqrt(len(lb))) + 2)
            ex = np.repeat((ls * (ls + 1)) ** 2 / rearth**4, 2 * ls + 1)[:len(lb)]
            rep.append(f'B NE{ne}: nnz={dref["B"].nnz} ||S-S^T||/||S||={symerr:.2e} test(rel to ||B||)={tres:.1e} '
                       f'lambda_min={lb[0]:.3e} (lambda_max*eps={lb[-1] * 2.2e-16:.1e}), next={lb[1]:.4e} '
                       f'(exact l=1: {4 / rearth**4:.4e}), #negative={np.sum(lb < 0)}, lambda_max={lb[-1]:.4e}, '
                       f'lambda_max/exact(l_max={int(np.sqrt(len(lb)))-1})={lb[-1] / ex[-1]:.2f}')
            rep.append(f'    B/exact at l=10: {lb[120] / ex[120]:.4f}, l=20: {lb[440] / ex[440]:.4f}, '
                       f'l=40: {lb[1680] / ex[1680]:.4f}' if len(lb) > 1700 else
                       f'    B/exact at l=10: {lb[120] / ex[120]:.4f}, l=20: {lb[440] / ex[440]:.4f}')
            # hyperviscosity spectra
            if ne == max(args.ne):
                nus = [0.0] + [nu * (30 / ne) ** 4 for nu in NU30]
                hv_cases = [(1, 2, 'analytic'), (4, 2, 'analytic')]
                hv_nus = {c: nus for c in hv_cases}
                rep.append(f'A = L - nu B, NE{ne}: max Re(lambda)/omega and damping per revolution')
                for case in hv_cases:
                    dc = load(os.path.join(args.mats, name(ne, *case)), ops=('L',))
                    for nu in hv_nus[case]:
                        rA = solve_A(cache, dc, dref['B'], nu, f'{name(ne, *case)}-A-nu{nu:.3g}')
                        resA[(ne,) + case + (nu,)] = rA
                        lam = rA['lam']
                        T = 2 * np.pi / om
                        fac = np.exp(lam.real * T)
                        low = np.abs(lam.imag) / om < 10
                        rep.append(f'  {label(*case)} nu={nu:.3g}: maxRe={lam.real.max() / om:.3e}, '
                                   f'per-rev factor: |Im|<10 omega median {np.median(fac[low]):.4f}, '
                                   f'|Im|>=10 median {np.median(fac[~low]):.3g}, min {fac.min():.2e}')
                # does the smallest hyperviscosity remove the conservation-form growth?
                for case in [(4, 0, 'analytic'), (1, 0, 'analytic')]:
                    dc = load(os.path.join(args.mats, name(ne, *case)), ops=('L',))
                    rA = solve_A(cache, dc, dref['B'], nus[1], f'{name(ne, *case)}-A-nu{nus[1]:.3g}')
                    rep.append(f'  {label(*case)} nu={nus[1]:.3g}: maxRe={rA["lam"].real.max() / om:.3e} '
                               f'(nu=0: {res[(ne,) + case]["lam"].real.max() / om:.3e})')
                plot_hypervis(resA, ne, hv_cases, hv_nus, om, args.out)

        if all((ne,) + c in res for c in CASES):
            plot_spectra(res, ne, om, args.out)
            plot_accuracy_nearest_integer(res, ne, om, args.out)
            plot_sorted_imag(res, ne, om, args.out)
            if ne in args.vectors_ne:
                plot_accuracy_mstar(res, ne, om, args.out)
                plot_localization(res, ne, om, args.out)
                if 'q' in res[(ne,) + CASES[0]]:
                    plot_error_vs_leff(res, ne, om, rearth, args.out)
                for case in [(4, 0, 'analytic'), (4, 2, 'analytic'), (1, 0, 'analytic')]:
                    plot_mode_maps(res, load(os.path.join(args.mats, name(ne, *case)), ops=()), ne, *case, om, args.out)
    if resB:
        plot_B(resB, sorted(resB), rearth, args.out)
    if args.mats30:
        r = ne30_analysis(cache, args.mats30, rep)
        if r:
            ne30_rk4_hv_growth(cache, args.mats30, r[0], rep)

    txt = '\n'.join(rep)
    print(txt)
    with open(os.path.join(args.out, 'results.txt'), 'w') as f:
        f.write(f'omega = u0/rearth = {om:.6e} 1/s\n' + txt + '\n')


if __name__ == '__main__':
    main()
