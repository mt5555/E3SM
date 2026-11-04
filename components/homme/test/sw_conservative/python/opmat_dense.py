"""
All eigenvalues of A = L - nu B by dense LAPACK (dgeev, eigenvalues only), for a big machine:
the NE30 matrices (n = 48,602) are 18.9 GB dense.

    OPENBLAS_NUM_THREADS=64 python opmat_dense.py --prefix MATS/ne30-sc1-f0-curlpsi-<case> \
        --mats MATS --nu30 0 2.5e15 3e15 --out OUT [--vectors 4] [--dry]
    python opmat_dense.py --prefix MATS/ne30-sc1-f2-curlpsi-<case> --out OUT --skew

For each nu30 (nu = nu30*(30/NE)^4, B from MATS/ne<NE>-sc1-f2-analytic) it writes
OUT/<case>-dense-nu<nu30>.npz with all eigenvalues lam/omega (omega = u0/a) and prints the
rightmost ones (null eigenvalues, |lam| < 1e-8 omega, such as the constant mode, are
counted separately).  Results already on disk
are skipped, so a run can be restarted.  Memory: one dense n x n real array (8 n^2 bytes,
overwritten in place by dgeev) plus O(n) workspace, about 20 GB at NE30.  Cost: about 15 n^3
flops (Hessenberg reduction 10/3 n^3, then the Hessenberg QR), 1.7e15 at NE30; the QR phase
does not scale well with threads, so expect hours.  --dry prints these estimates and stops.

--vectors K: the eigenvectors of the K rightmost eigenvalues (constant mode excluded) by
inverse iteration with a sparse LU of A - lam I, with their residuals |A v - lam v|/|lam|
(columns of `vectors`, unit 2-norm, in the dof order of the opmat files).

--skew (split form, nu = 0): S = M^1/2 L M^-1/2 is skew-symmetric for theta_advect_form=2 with
a DSS(curl psi) velocity, so its eigenvalues are +-i sigma with sigma the singular values of S;
dgesdd values-only is symmetric-type work and parallelizes far better than dgeev.  The skew
error |S + S^T|/|S| is reported, and the Re lam = 0 claim rests on it.

--sym: for flows that keep the symmetries of the cubed-sphere grid (the solid-body rotation
with alpha = 0 or pi/4, opmat_sym.py), solve the dense eigenvalue problem block by block:
8 blocks of about n/8 (alpha = 0, two of them complex conjugates of others) or 4 real blocks
of about n/4 (alpha = pi/4).  NE30 alpha = pi/4 then takes 15-20 min and 2.6 GB on 4 threads,
so it runs on a laptop-sized machine.  The rotation and the antipodal map must commute with
A = L - nu B (checked first, to 1e-12 relative); P1 and other localized flows break them and
are refused.  Each block is checkpointed in OUT/<case>-dense-nu<nu30>-ckpt/ and the run
resumes from there; the block label of each eigenvalue is saved as `block`.
"""
import argparse
import os
import shutil
import time
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from opmat_io import load     # numpy and scipy only, nothing else needed on the big machine


def nu_from30(n30, ne):
    return n30 * (30.0 / ne) ** 4


def load_B(mats, ne, lat, lon):
    """B from MATS/ne<NE>-sc1-f2-analytic, checked to use the dof order of L (as opmat_flows_hv.load_B)"""
    dB = load(os.path.join(mats, f'ne{ne}-sc1-f2-analytic'), ops=('B',))
    if 'B' not in dB:
        raise SystemExit(f'{mats}/ne{ne}-sc1-f2-analytic: no B')
    if np.abs(dB['lat'] - lat).max() > 1e-12 or np.abs(np.angle(np.exp(1j * (dB['lon'] - lon)))).max() > 1e-12:
        raise SystemExit('B: dof ordering differs from the L files')
    return dB['B']


def report(lam, om_tol=1e-6, null=1e-8, top=12):
    """lam already divided by omega; |lam| < null (x omega) counts as a null eigenvalue"""
    isnull = np.abs(lam) < null
    r = lam[~isnull]
    idx = np.argsort(-r.real)[:top]
    print(f'  {isnull.sum()} null eigenvalue(s); max Re = {r.real.max():.6g} omega; '
          f'{(r.real > om_tol).sum()} with Re > {om_tol:g} omega; '
          f'max |Im| = {np.abs(r.imag).max():.4g} omega', flush=True)
    for j in idx:
        print(f'    {r[j].real:+.6f} {r[j].imag:+.6f}i', flush=True)
    return r[idx]


def sym_eigvals(A, d, f, om, tol=1e-12):
    """all eigenvalues of A (divided by om) via the symmetry blocks of opmat_sym, with a
    checkpoint per block in f[:-4] + '-ckpt'; returns (lam, block label)"""
    from opmat_sym import generators, check_commutes, block_eigs
    try:
        gens = generators(d['lat'], d['lon'], d['meta']['alpha'])
    except (ValueError, KeyError) as e:
        raise SystemExit(f'--sym: no grid symmetry for this case ({e})')
    comm = max(check_commutes(A, perm) for perm, _, _ in gens)
    if comm > tol:
        raise SystemExit(f'--sym: the grid symmetries do not commute with A (relative error {comm:.1e}); '
                         'run without --sym')
    ck = f[:-4] + '-ckpt'
    os.makedirs(ck, exist_ok=True)
    def ckf(k):
        return os.path.join(ck, f'block{k}.npz')
    key = np.array([A.shape[0], spla.norm(A)])
    def load_ck(k):
        if not os.path.exists(ckf(k)):
            return None
        c = np.load(ckf(k))
        return c['lam'] if np.allclose(c['key'], key, rtol=1e-12, atol=0) else None
    skip = [k for k in range(8) if load_ck(k) is not None]
    if skip:
        print(f'  resuming: blocks {skip} from {ck}', flush=True)
    t0 = time.time()
    lams, labs = [], []
    for lab, lam, Q, U, info in block_eigs(A, d['lat'], d['lon'], d['meta']['alpha'], skip=skip):
        if lam is None:
            lam, how = load_ck(lab), 'from checkpoint'
        else:
            lam = lam / om
            tmp = ckf(lab) + '.tmp.npz'
            np.savez(tmp, lam=lam, key=key)
            os.replace(tmp, ckf(lab))
            how = f'{time.time() - t0:.0f} s, off-block residual {info["offblock"]:.1e}'
        print(f'  block {lab} of {info["m"]}: size {info["size"]}, {how}, commute {info["commute"]:.1e}, '
              f'max Re {lam.real.max():.4g} omega', flush=True)
        lams.append(lam); labs.append(np.full(len(lam), lab))
    return np.concatenate(lams), np.concatenate(labs), ck


def inverse_vectors(A, lams, om, iters=12):
    n = A.shape[0]
    Ac = A.astype(complex).tocsc()
    I = sp.identity(n, format='csc', dtype=complex)
    V, res = [], []
    rng = np.random.default_rng(0)
    for lam in lams:
        sg = lam * om * (1 + 1e-9) + 1e-9 * om
        lu = spla.splu((Ac - sg * I).tocsc(), permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0.1)
        v = rng.standard_normal(n) + 0j
        for _ in range(iters):
            v = lu.solve(v)
            v /= np.linalg.norm(v)
        lv = (v.conj() @ (Ac @ v))                 # Rayleigh quotient, close to lam*om
        r = np.linalg.norm(Ac @ v - lv * v) / abs(lv)
        print(f'    vector at {lam:.6f}: Rayleigh {lv / om:.6f}, residual {r:.2e}', flush=True)
        V.append(v); res.append(r)
    return np.array(V).T, np.array(res)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--prefix', required=True, help='opmat prefix of L, e.g. MATS/ne30-sc1-f0-curlpsi-<case>')
    p.add_argument('--mats', help='directory with ne<NE>-sc1-f2-analytic (B), needed for nu30 > 0')
    p.add_argument('--nu30', type=float, nargs='+', default=[0.0])
    p.add_argument('--out', required=True)
    p.add_argument('--vectors', type=int, default=0)
    p.add_argument('--skew', action='store_true', help='split form at nu=0: singular values of S')
    p.add_argument('--sym', action='store_true',
                   help='solve by grid-symmetry blocks (solid-body rotation, alpha = 0 or pi/4)')
    p.add_argument('--dry', action='store_true')
    a = p.parse_args()

    d = load(a.prefix, ops=('L',))
    n, ne, meta = d['n'], d['meta']['ne'], d['meta']
    om = meta['u0'] / meta['rearth']
    case = os.path.basename(a.prefix)
    print(f'{case}: n = {n}, dense array {8 * n * n / 1e9:.1f} GB, about {15 * n ** 3:.2e} flops per dgeev '
          f'({4 * n ** 3:.2e} for --skew)', flush=True)
    if a.sym and a.skew:
        raise SystemExit('--sym and --skew together are not supported')
    if a.sym:
        from opmat_sym import symmetry_blocks
        nk = [Q.shape[1] for _, Q, src in symmetry_blocks(d['lat'], d['lon'], meta['alpha'])[1] if src is None]
        print(f'  --sym: {len(nk)} blocks to solve, sizes {nk}, about {sum(15 * k ** 3 for k in nk):.2e} flops, '
              f'largest dense block {8 * max(nk) ** 2 / 1e9:.2f} GB', flush=True)
    if a.dry:
        return
    os.makedirs(a.out, exist_ok=True)

    if a.skew:
        f = os.path.join(a.out, f'{case}-dense-skew.npz')
        if os.path.exists(f):
            print(f'{f} exists'); return
        sq = np.sqrt(d['mass'])
        S = (sp.diags(sq) @ d['L'] @ sp.diags(1 / sq)).tocsr()
        err = spla.norm(S + S.T) / spla.norm(S)
        print(f'  skew error |S+S^T|/|S| = {err:.2e}', flush=True)
        t0 = time.time()
        sig = la.svdvals(S.toarray(order='F'), overwrite_a=True, check_finite=False) / om
        print(f'  svdvals: {time.time() - t0:.0f} s; max sigma = {sig.max():.6g} omega, '
              f'{(sig < 1e-8 * sig.max()).sum()} zero', flush=True)
        sig_s = np.sort(sig)                      # eigenvalues +-i sigma, each pair gives sigma twice
        lam = 1j * sig_s * np.where(np.arange(n) % 2 == 0, 1, -1)
        np.savez(f, sigma=sig, lam=lam, skew_err=err, n=n, ne=ne, nu30=0.0, omega=om, case=case,
                 time=time.time() - t0)
        print(f)
        return

    B = load_B(a.mats, ne, d['lat'], d['lon']) if any(x > 0 for x in a.nu30) else None
    for n30 in a.nu30:
        f = os.path.join(a.out, f'{case}-dense-nu{n30:.4g}.npz')
        if os.path.exists(f):
            print(f'{f} exists: skipped'); continue
        A = d['L'] if n30 == 0 else (d['L'] - nu_from30(n30, ne) * B).tocsr()
        t0 = time.time()
        ck = None
        if a.sym:
            print(f'nu30 = {n30:.4g} (nu = {nu_from30(n30, ne):.4g}): by symmetry blocks', flush=True)
            lam, block, ck = sym_eigvals(A, d, f, om)
        else:
            lam = la.eigvals(A.toarray(order='F'), overwrite_a=True, check_finite=False) / om
        dt = time.time() - t0
        print(f'nu30 = {n30:.4g} (nu = {nu_from30(n30, ne):.4g}): {"blocks" if a.sym else "dgeev"} {dt:.0f} s',
              flush=True)
        top = report(lam)
        out = dict(lam=lam, nu30=n30, nu=nu_from30(n30, ne), n=n, ne=ne, omega=om, time=dt, case=case)
        if a.sym:
            out.update(block=block)
        if a.vectors:
            V, res = inverse_vectors(A, top[:a.vectors], om)
            out.update(vectors=V, vector_lam=top[:a.vectors], vector_res=res,
                       lat=d['lat'], lon=d['lon'], mass=d['mass'])
        np.savez(f, **out)
        if ck is not None:
            shutil.rmtree(ck)
        print(f, flush=True)


if __name__ == '__main__':
    main()
