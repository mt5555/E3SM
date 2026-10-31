#!/usr/bin/env python3
"""
Checks of the symmetry-block eigensolver (opmat_sym.py) on extracted operators.

usage: opmat_sym_test.py --mats DIR --test blocks|full|old|vectors|proj --ne 8 10 [--old OLD.py]
                         [--ref DIR --new DIR]
                         [--sc 1 4] [--forms 2 0]

  blocks   Q orthonormal, block sizes adding up to n, P_R Q = w^k Q, P_J Q = s Q,
           commutation and off-block residuals (L of every case, and B)
  full     block eigenvalues against a dense eigensolve of the whole matrix, matched
           one to one (linear assignment), for L and for the symmetrized B
  old      block eigenvalues against another symmetry module (--old, e.g. the
           rotation-only opmat_sym.py from git), matched one to one, with timings
  vectors  eigenpair residuals ||L v - lam v|| / (||L|| ||v||) of every block, and
           orthonormality of the eigenvectors of the skew part over all blocks
  proj     opmat_proj.py results (--new DIR) against earlier ones (--ref DIR): error
           fields, modal budget, and the per-mode diagnostics of isolated eigenvalues
"""
import argparse
import importlib.util
import os
import time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linear_sum_assignment

from opmat_io import load
from opmat_eigs import name, omega, mass_sym
from opmat_sym import symmetry_blocks, block_eigs, eigvals_by_symmetry, check_commutes


def match(a, b):
    """max |a_i - b_p(i)| over the one-to-one assignment p minimizing the total distance"""
    i, j = linear_sum_assignment(np.abs(a[:, None] - b[None, :]))
    return np.abs(a[i] - b[j]).max()


def cases(args, ne):
    for sc in args.sc:
        for f in args.forms:
            path = os.path.join(args.mats, name(ne, sc, f, 'analytic'))
            if os.path.exists(path + '.meta'):
                yield f'ne{ne}-sc{sc}-f{f}', load(path, ops=('L',))


def symB(args, ne):
    path = os.path.join(args.mats, name(ne, 1, 2, 'analytic'))
    if not os.path.exists(path + '.B.r0000.bin'):
        return None, None
    d = load(path, ops=('B',))
    S = mass_sym(d['B'], d['mass'])
    return d, ((S + S.T) / 2).tocsr()


def test_blocks(args, ne):
    for nm, d in cases(args, ne):
        alpha = d['meta']['alpha']
        gens, blocks = symmetry_blocks(d['lat'], d['lon'], alpha)
        m = gens[0][1]
        n = d['n']
        Ps = [sp.csr_matrix((np.ones(n), (perm, np.arange(n))), shape=(n, n)) for perm, _, _ in gens]
        orth = eig = off = 0.0
        for lab, Q, src in blocks:
            k, b = lab % m, lab // m
            orth = max(orth, abs(Q.conj().T @ Q - sp.eye(Q.shape[1])).max())
            eig = max(eig, abs(Ps[0] @ Q - np.exp(2j * np.pi * k / m) * Q).max(),
                      abs(Ps[1] @ Q - (-1.0) ** b * Q).max())
            off = max(off, sp.linalg.norm(d['L'] @ Q - Q @ (Q.conj().T @ d['L'] @ Q)) / sp.linalg.norm(d['L'] @ Q))
        comm = [check_commutes(d['L'], perm) for perm, _, _ in gens]
        sizes = [Q.shape[1] for _, Q, _ in blocks]
        print(f'{nm}: n={n} blocks={sizes} (sum {sum(sizes)}), real={[not np.iscomplexobj(Q.data) for _, Q, _ in blocks]}\n'
              f'    max|Q^H Q - I|={orth:.1e}  max|P Q - chi Q|={eig:.1e}  commute R {comm[0]:.1e} J {comm[1]:.1e}'
              f'  off-block {off:.1e}', flush=True)
    d, Ss = symB(args, ne)
    if d is not None:
        gens, _ = symmetry_blocks(d['lat'], d['lon'], 0.0)
        print(f'ne{ne} B: commute R {check_commutes(Ss, gens[0][0]):.1e} J {check_commutes(Ss, gens[1][0]):.1e}', flush=True)


def test_full(args, ne):
    for nm, d in cases(args, ne):
        om = omega(d['meta'])
        t = time.time()
        lam, _, info, _ = eigvals_by_symmetry(d['L'], d['lat'], d['lon'], d['meta']['alpha'])
        tb = time.time() - t
        t = time.time()
        ref = np.linalg.eigvals(d['L'].toarray())
        tf = time.time() - t
        print(f'{nm}: blocks {info["sizes"]}, max|lam_block - lam_full| = {match(lam, ref) / om:.2e} omega '
              f'(max|lam| {np.abs(ref).max() / om:.0f} omega); {tb:.1f} s vs {tf:.1f} s full', flush=True)
    d, Ss = symB(args, ne)
    if d is not None:
        _, blocks = symmetry_blocks(d['lat'], d['lon'], 0.0)
        lam = []
        for lab, Q, src in blocks:
            Bk = (Q.conj().T @ Ss @ Q).toarray()
            lam.append(np.linalg.eigvalsh((Bk + Bk.conj().T) / 2))
        lam = np.sort(np.concatenate(lam))
        ref = np.linalg.eigvalsh(Ss.toarray())
        print(f'ne{ne} B: max|lam_block - lam_full| / max lam = {np.abs(lam - ref).max() / ref.max():.2e}', flush=True)


def test_old(args, ne):
    spec = importlib.util.spec_from_file_location('opmat_sym_old', args.old)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    for nm, d in cases(args, ne):
        om = omega(d['meta'])
        t = time.time()
        lam, _, info, _ = eigvals_by_symmetry(d['L'], d['lat'], d['lon'], d['meta']['alpha'])
        tn = time.time() - t
        t = time.time()
        ref, _, info0, _ = old.eigvals_by_symmetry(d['L'], d['lat'], d['lon'], d['meta']['alpha'])
        to = time.time() - t
        print(f'{nm}: new blocks {info["sizes"]} {tn:.0f} s, old blocks {info0["sizes"]} {to:.0f} s; '
              f'max|lam_new - lam_old| = {match(lam, ref) / om:.2e} omega', flush=True)


def test_vectors(args, ne):
    for nm, d in cases(args, ne):
        alpha = d['meta']['alpha']
        L = d['L']
        nL = np.abs(L).sum(axis=1).max()
        res = 0.0
        for lab, lam, Q, U, info in block_eigs(L, d['lat'], d['lon'], alpha, vectors=True):
            V = Q @ U
            res = max(res, (np.linalg.norm(L @ V - V * lam, axis=0) / np.linalg.norm(V, axis=0)).max() / nL)
        # skew part: orthonormal over all blocks (Euclidean in M^1/2 coordinates)
        S = mass_sym(L, d['mass']).tocsr()
        herm = lambda Ak: 0.5j * (Ak - Ak.conj().T)
        Vs = [Q @ U for _, _, Q, U, _ in block_eigs(S, d['lat'], d['lon'], alpha, vectors=True, hermitian=herm)]
        V = np.concatenate(Vs, axis=1)
        orth = np.abs(V.conj().T @ V - np.eye(V.shape[1])).max()
        print(f'{nm}: max eigenpair residual (L) {res:.1e} (relative to ||L||_inf); '
              f'skew part: max|V^H V - I| = {orth:.1e} over {V.shape[1]} vectors', flush=True)


def test_proj(args, ne):
    for sc in args.sc:
        nm = name(ne, sc, 2, 'analytic') + '-proj.npz'
        a, b = np.load(os.path.join(args.ref, nm)), np.load(os.path.join(args.new, nm))
        om, mass = float(a['omega']), a['mass']
        nrm = lambda x: np.sqrt((mass * np.abs(x) ** 2).sum())
        rel = lambda x, y: nrm(x - y) / nrm(x)
        ia, ib = np.argsort(a['lam'].imag), np.argsort(b['lam'].imag)
        la, lb = a['lam'][ia], b['lam'][ib]
        # per-mode quantities only for isolated eigenvalues (basis-independent)
        gap = np.minimum(np.r_[np.inf, np.diff(la.imag)], np.r_[np.diff(la.imag), np.inf])
        iso = gap > 1e-6 * om
        amp = lambda r, i: (np.abs(r['c']) * r['snorm'])[i]
        Ea, Eb = np.abs(a['g_sd'][ia]) * amp(a, ia), np.abs(b['g_sd'][ib]) * amp(b, ib)
        q = (np.abs(a['q'][ia] - b['q'][ib]) / np.abs(a['q'][ia]).max())[iso].max() if 'q' in a.files and 'q' in b.files else np.nan
        print(f'{nm}: n={len(la)}, max|dlam|={np.abs(la - lb).max() / om:.1e} omega; e_sd {rel(a["e_sd"], b["e_sd"]):.1e}, '
              f'e_fd {rel(a["e_fd"], b["e_fd"]):.1e}, e_res {rel(a["e_res"], b["e_res"]):.1e} (relative differences)\n'
              f'    sqrt(sum E_sd^2): {np.sqrt((Ea ** 2).sum()):.6e} vs {np.sqrt((Eb ** 2).sum()):.6e}; '
              f'top_E {np.abs(a["top_E"] - b["top_E"]).max() / a["top_E"].max():.1e}, '
              f'top_lam {np.abs(np.sort_complex(a["top_lam"]) - np.sort_complex(b["top_lam"])).max() / om:.1e} omega\n'
              f'    isolated modes ({iso.sum()} of {len(la)}): max|dE_j|/max E {np.abs(Ea - Eb)[iso].max() / Ea.max():.1e}, '
              f'm* equal {np.mean(a["mstar"][ia][iso] == b["mstar"][ib][iso]):.4f} (m* is roundoff where the '
              f'sphere-integrated Fourier power vanishes), '
              f'max|dcorner| {np.abs(a["corner"][ia] - b["corner"][ib])[iso].max():.1e}, '
              f'max|dedge| {np.abs(a["edge"][ia] - b["edge"][ib])[iso].max():.1e}, max|dq|/max q {q:.1e}', flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mats', required=True)
    p.add_argument('--test', required=True, choices=['blocks', 'full', 'old', 'vectors', 'proj'])
    p.add_argument('--ne', type=int, nargs='+', default=[8])
    p.add_argument('--sc', type=int, nargs='+', default=[1, 4])
    p.add_argument('--forms', type=int, nargs='+', default=[2, 0])
    p.add_argument('--old', help='symmetry module to compare with (--test old)')
    p.add_argument('--ref', help='earlier opmat_proj.py output (--test proj)')
    p.add_argument('--new', help='new opmat_proj.py output (--test proj)')
    args = p.parse_args()
    for ne in args.ne:
        globals()['test_' + args.test](args, ne)


if __name__ == '__main__':
    main()
