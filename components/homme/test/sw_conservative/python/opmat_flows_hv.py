"""
Hyperviscosity needed to stabilize the conservation-form transport operator L for the
general flows of opmat_flows.py (u = DSS(curl psi), psi = psi_sb + eps*a*u0*Y_l^m).
B is HOMME's hyperviscosity operator (hypervis_order=2, B p = rspheremp*DSS(biharmonic_wk p)),
read from the ne{NE}-sc1-f2-analytic files (B depends only on NE).

    python opmat_flows_hv.py --mats DIR --out OUT [--ne 8] [--sc 1] [--select 'ylm*_0_*' ...]
                             [--spectral 12] [--discrete 4] [--dt DT]

Three thresholds, all reported as nu30 = nu*(NE/30)^4 (m^4/s):
  * energy:     nu_E(tol), the smallest nu with sym(S) - nu*Bs <= tol*omega (S = M^1/2 L M^-1/2,
                Bs = M^1/2 B M^-1/2): no growth of the mass-weighted l2 norm, even transient.
                One generalized symmetric eigenproblem per tol (constants, the common null
                vector of sym(S) and Bs, are shifted out).  Upper bound for nu_spec.
  * first order: for the top --nfo growing eigenvalues of L, nu_j = Re lam_j / Re(y^T B v / y^T v)
                (v, y right/left eigenvectors, shift-invert).  A screen only: L is far from normal.
  * spectral:   nu_spec, the smallest nu above which max Re lambda(L - nu B) <= tol*omega,
                from a scan over nu and bisection (dense eigenvalues), for the --spectral
                fields with the largest screening value.
  * discrete:   HOMME's step p <- (I - dt nu B) R4(dt L) p (advance_nonstag_rk4, then
                advance_hypervis with hypervis_subcycle=1): the window [nu_lo, nu_hi] of nu
                with growth rate log(rho)/dt <= tol*omega, for the --discrete worst fields.
Cached in OUT/cache/hv-*.npz; table in OUT/results_hv-<tag>.txt; figures in OUT.
"""
import argparse
import glob
import os
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from opmat_io import load, read_meta
from opmat_eigs import Cache
from opmat_proj import timestep
from opmat_flows import analyze, r_of, regions, enrichment, style, cube_edge_lines, GROW

TOL = 1e-6            # stable: growth rate <= TOL * omega
NU30_SCAN = np.concatenate([[0.0], np.logspace(10, 16, 25)])


def nu30(nu, ne):
    return nu * (ne / 30.0) ** 4


def nu_from30(n30, ne):
    return n30 * (30.0 / ne) ** 4


def load_B(mats, ne, lat, lon):
    pB = os.path.join(mats, f'ne{ne}-sc1-f2-analytic')
    dB = load(pB, ops=('B',))
    if 'B' not in dB:
        raise SystemExit(f'{pB}: no B (extract with swtc1-opmat.sh {ne} 1 2 analytic B)')
    if np.abs(dB['lat'] - lat).max() > 1e-12 or np.abs(np.angle(np.exp(1j * (dB['lon'] - lon)))).max() > 1e-12:
        raise SystemExit(f'{pB}: dof ordering differs from the L files')
    return dB['B']


# ---------------------------------------------------------------- thresholds
def sym_parts(L, B, mass):
    sq = np.sqrt(mass)
    S = (sp.diags(sq) @ L @ sp.diags(1 / sq)).toarray()
    Bs = (sp.diags(sq) @ B @ sp.diags(1 / sq)).toarray()
    Ss = 0.5 * (S + S.T)
    Bs = 0.5 * (Bs + Bs.T)
    q = sq / np.linalg.norm(sq)          # M^1/2 * constants
    return Ss, Bs, q


def energy_nu(Ss, Bs, q, om, tol):
    """smallest nu with lambda_max(Ss - nu Bs) <= tol*om, and the maximizing vector (M^1/2-scaled)"""
    n = len(q)
    c = np.abs(np.diag(Bs)).max()
    w, V = sla.eigh(Ss - tol * om * np.eye(n), Bs + c * np.outer(q, q), subset_by_index=[n - 1, n - 1])
    return max(w[0], 0.0), V[:, 0]


def first_order(L, B, lam, om, k):
    """first-order nu_j for the top k growing eigenvalues (one per conjugate pair)"""
    o = np.argsort(-lam.real)
    sel = [j for j in o if lam[j].imag >= 0 and lam[j].real > GROW * om][:k]
    Lc, LTc = L.astype(complex), L.T.astype(complex).tocsr()
    out = []
    for j in sel:
        s = lam[j] + 1e-9 * om * (1 + 1j)
        ev, v = spla.eigs(Lc, k=1, sigma=s)
        _, y = spla.eigs(LTc, k=1, sigma=s)
        v, y = v[:, 0], y[:, 0]
        rate = ((y @ (B @ v)) / (y @ v)).real
        out.append((ev[0], rate, ev[0].real / rate if rate > 0 else np.inf, v))
    return out


def max_re(L, B, nu):
    return np.linalg.eigvals((L - nu * B).toarray()).real.max()


def spectral_nu(L, B, om, ne, tol, nbis=8):
    """nu_spec from a scan over NU30_SCAN and bisection (in log nu) after the last unstable point"""
    g = np.array([max_re(L, B, nu_from30(n30, ne)) for n30 in NU30_SCAN])
    bad = np.where(g > tol * om)[0]
    if len(bad) == 0:
        return 0.0, g
    i = bad[-1]
    if i == len(NU30_SCAN) - 1:
        return np.inf, g
    lo, hi = max(NU30_SCAN[i], 1e9), NU30_SCAN[i + 1]
    for _ in range(nbis):
        mid = np.sqrt(lo * hi)
        if max_re(L, B, nu_from30(mid, ne)) > tol * om:
            lo = mid
        else:
            hi = mid
    return hi, g


def rk4_matrix(L, dt):
    Z = dt * L.toarray()
    n = Z.shape[0]
    R = np.eye(n) + Z / 4
    R = np.eye(n) + Z @ R / 3
    R = np.eye(n) + Z @ R / 2
    return np.eye(n) + Z @ R


def discrete_window(L, B, om, ne, dt, tol, nbis=7):
    """growth rate of (I - dt nu B) R4(dt L) over NU30_SCAN, and the stable window"""
    R = rk4_matrix(L, dt)
    Bd = B.toarray()

    def rate(n30):
        A = R - dt * nu_from30(n30, ne) * (Bd @ R)
        return np.log(np.abs(np.linalg.eigvals(A)).max()) / dt

    g = np.array([rate(n30) for n30 in NU30_SCAN])
    ok = g <= tol * om
    if not ok.any():
        return np.nan, np.nan, g
    i0 = np.argmax(ok)                       # first stable point
    i1 = i0 + np.argmin(ok[i0:]) - 1 if not ok[i0:].all() else len(ok) - 1
    lo, hi = (max(NU30_SCAN[i0 - 1], 1e9), NU30_SCAN[i0]) if i0 > 0 else (0.0, 0.0)
    for _ in range(nbis if i0 > 0 else 0):
        mid = np.sqrt(lo * hi)
        lo, hi = (lo, mid) if rate(mid) <= tol * om else (mid, hi)
    nlo = hi
    if i1 == len(ok) - 1:
        nhi = np.inf
    else:
        lo, hi = NU30_SCAN[i1], NU30_SCAN[i1 + 1]
        for _ in range(nbis):
            mid = np.sqrt(lo * hi)
            lo, hi = (mid, hi) if rate(mid) <= tol * om else (lo, mid)
        nhi = lo
    return nlo, nhi, g


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--mats', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--ne', type=int, default=8)
    p.add_argument('--sc', type=int, default=1)
    p.add_argument('--select', default=['*'], nargs='+')
    p.add_argument('--baseline', action='store_true')
    p.add_argument('--nfo', type=int, default=20, help='growing modes in the first-order screen')
    p.add_argument('--spectral', type=int, default=12, help='fields for the spectral bisection')
    p.add_argument('--discrete', type=int, default=4, help='fields for the RK4+hyperviscosity window')
    p.add_argument('--dt', type=float, default=None,
                   help='time step (default: as opmat_proj.timestep, 1800 s*30/NE per 12-day revolution, 7200 s at NE8)')
    p.add_argument('--tol', type=float, default=TOL)
    p.add_argument('--tag', default='hv')
    args = p.parse_args()
    ne, tol = args.ne, args.tol
    dt = args.dt or timestep(ne, 1800.0, 12 * 86400.0)[0]
    os.makedirs(args.out, exist_ok=True)
    cache = Cache(args.out)

    files = sorted(set(f[:-5] for sel in args.select for f in
                       glob.glob(os.path.join(args.mats, f'ne{ne}-sc{args.sc}-f0-curlpsi-{sel}.meta'))))
    if args.baseline:
        files = [os.path.join(args.mats, f'ne{ne}-sc{args.sc}-f0-curlpsi')] + files
    B = None
    rows = []
    for path in files:
        nm = os.path.basename(path)
        meta = read_meta(path)
        r = analyze(cache, path, 3)
        om = float(r['om'])

        def screen():
            d = load(path, ops=('L',))
            Bl = load_B(args.mats, ne, d['lat'], d['lon'])
            Ss, Bs, q = sym_parts(d['L'], Bl, d['mass'])
            nE0, _ = energy_nu(Ss, Bs, q, om, 0.0)
            nE, vE = energy_nu(Ss, Bs, q, om, tol)
            fo = first_order(d['L'], Bl, r['lam'], om, args.nfo)
            j = int(np.argmax([f[2] for f in fo])) if fo else -1
            return dict(nuE0=nE0, nuE=nE, vE=vE / np.sqrt(d['mass']),
                        fo_lam=np.array([f[0] for f in fo]), fo_rate=np.array([f[1] for f in fo]),
                        fo_nu=np.array([f[2] for f in fo]),
                        fo_v=fo[j][3] if fo else np.zeros(len(d['mass']), complex))
        h = cache.get('hv-screen-' + nm, screen)
        reg = regions(r, ne, r['speed'])
        row = dict(name=nm, path=path, l=meta.get('psi_l', 1), m=meta.get('psi_m', 0), r=r_of(meta), om=om,
                   maxre=r['lam'].real.max() / om, ngrow=int((r['lam'].real > GROW * om).sum()),
                   nuE0=nu30(float(h['nuE0']), ne), nuE=nu30(float(h['nuE']), ne),
                   nufo=nu30(h['fo_nu'].max(), ne) if len(h['fo_nu']) else 0.0,
                   enrE=enrichment(h['vE'], r['mass'], reg), enrfo=enrichment(h['fo_v'], r['mass'], reg),
                   h=h, res=r)
        rows.append(row)
        print(f'{nm}: maxRe={row["maxre"]:.3f} nu30 first-order={row["nufo"]:.2e} '
              f'energy(tol)={row["nuE"]:.2e} energy(0)={row["nuE0"]:.2e}', flush=True)

    # spectral and discrete thresholds for the fields with the largest screening value
    order = sorted(rows, key=lambda row: -max(row['nufo'], 0))
    for row in order[:args.spectral]:
        nm = row['name']

        def spec():
            d = load(row['path'], ops=('L',))
            Bl = load_B(args.mats, ne, d['lat'], d['lon'])
            ns, g = spectral_nu(d['L'], Bl, row['om'], ne, tol)
            return dict(nus=ns, g=g)
        s = cache.get(f'hv-spec-tol{tol:g}-' + nm, spec)
        row['nus'], row['gs'] = float(s['nus']), s['g'] / row['om']
        print(f'{nm}: nu30 spectral={row["nus"]:.2e}', flush=True)
    order = sorted([row for row in rows if 'nus' in row], key=lambda row: -row['nus'])
    for row in order[:args.discrete]:
        nm = row['name']

        def disc():
            d = load(row['path'], ops=('L',))
            Bl = load_B(args.mats, ne, d['lat'], d['lon'])
            lo, hi, g = discrete_window(d['L'], Bl, row['om'], ne, dt, tol)
            return dict(lo=lo, hi=hi, g=g)
        s = cache.get(f'hv-disc-dt{dt:g}-tol{tol:g}-' + nm, disc)
        row['dlo'], row['dhi'], row['gd'] = float(s['lo']), float(s['hi']), s['g'] / row['om']
        print(f'{nm}: discrete window nu30 [{row["dlo"]:.2e}, {row["dhi"]:.2e}] at dt={dt:g}', flush=True)

    # ---------------------------------------------------------------- report
    rep = [f'NE{ne} sc{args.sc}: hyperviscosity thresholds as nu30 = nu*(NE/30)^4 [m^4/s]; stable = growth rate '
           f'<= {tol:g} omega; dt = {dt:g} s for the discrete window',
           '  first-order: max_j Re lam_j / Re(y^T B v / y^T v) over the top growing modes; energy(tol), energy(0): '
           'no growth of the mass-weighted l2 norm; spectral: max Re lambda(L - nu B) <= tol omega; '
           'enrichment [corner cubeedge elemedge slow] of the hardest first-order mode',
           f'{"case":38s} {"l":>3s} {"r":>6s} {"maxRe/om":>9s} {"ngrow":>5s} {"first-ord":>9s} {"spectral":>9s} '
           f'{"energy":>9s} {"energy0":>9s} {"window lo":>9s} {"window hi":>9s}  enrichment']
    for row in sorted(rows, key=lambda row: (row['l'], row['r'])):
        e = row['enrfo']
        rep.append(f'{row["name"]:38s} {row["l"]:3d} {row["r"]:6.2f} {row["maxre"]:9.3e} {row["ngrow"]:5d} '
                   f'{row["nufo"]:9.2e} {row.get("nus", np.nan):9.2e} {row["nuE"]:9.2e} {row["nuE0"]:9.2e} '
                   f'{row.get("dlo", np.nan):9.2e} {row.get("dhi", np.nan):9.2e}  '
                   f'[{e["corner"]:.1f} {e["cubeedge"]:.1f} {e["elemedge"]:.1f} {e["slow"]:.1f}]')
    txt = '\n'.join(rep)
    print(txt)
    with open(os.path.join(args.out, f'results_hv-{args.tag}.txt'), 'w') as f:
        f.write(txt + '\n')
    plots(rows, ne, args, dt)


def plots(rows, ne, args, dt):
    tag = f'ne{ne}-sc{args.sc}-{args.tag}'
    # thresholds against max Re lambda
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
    ax = axs[0]
    x = np.array([row['maxre'] for row in rows])
    ax.scatter(x, [row['nufo'] for row in rows], s=12, color='#1baf7a', label='first order')
    sel = [row for row in rows if 'nus' in row]
    ax.scatter([row['maxre'] for row in sel], [row['nus'] for row in sel], s=30, marker='D', facecolor='none',
               edgecolor='#eb6834', label='spectral (bisection)')
    ax.scatter(x, [row['nuE'] for row in rows], s=12, marker='s', color='#6250d6', label=f'energy, tol')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('max Re λ(L) / ω')
    ax.set_ylabel('ν₃₀ needed [m⁴/s]')
    ax.legend(fontsize=7)
    style(ax)
    ax = axs[1]
    for row in sorted(sel, key=lambda row: -row['nus'])[:6]:
        ax.plot(np.maximum(NU30_SCAN[1:], 1), np.maximum(row['gs'][1:], 1e-12), 'o-', ms=3,
                label=row['name'].split('curlpsi-')[-1])
    ax.axhline(args.tol, color='#898781', lw=0.8, ls=':')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('ν₃₀ [m⁴/s]')
    ax.set_ylabel('max Re λ(L − νB) / ω')
    ax.legend(fontsize=7)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, f'{tag}-nu.png'))
    plt.close(fig)
    # hardest-to-damp modes (first order) of the worst fields
    worst = sorted(sel, key=lambda row: -row['nus'])[:4]
    if not worst:
        return
    edges = cube_edge_lines()
    fig, axs = plt.subplots(len(worst), 2, figsize=(9, 2.3 * len(worst)), squeeze=False)
    for i, row in enumerate(worst):
        r = row['res']
        lon, lat = np.degrees(r['lon']) % 360, np.degrees(r['lat'])
        for k, (v, ttl) in enumerate(((row['h']['fo_v'], 'hardest first-order mode'),
                                      (row['h']['vE'], 'energy-maximizing vector'))):
            ax = axs[i][k]
            for lo, la in edges:
                ax.plot(lo, la, color='#2a78d6', lw=0.5)
            e = np.abs(v) ** 2
            e /= e.max()
            o = np.argsort(e)
            ax.scatter(lon[o], lat[o], c=e[o], s=2 + 14 * e[o], cmap='magma_r', vmin=0, vmax=1, linewidths=0)
            ax.set_xlim(0, 360)
            ax.set_ylim(-90, 90)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f'{row["name"].split("curlpsi-")[-1]}: {ttl}', fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, f'{tag}-modes.png'))
    plt.close(fig)


if __name__ == '__main__':
    main()
