"""
Eigenvalues of the swtc1 transport operator L for general steady divergence-free
velocities u = DSS(curl psi), psi = psi_solidbody + eps*a*u0*Y_l^m (SWEQX_TC1_PSI=ylm,
see ../swtc1-opmat-ylm.sh), conservation form (f0) against split form (f2).

    python opmat_flows.py --mats DIR --out OUT [--ne 8] [--sc 1] [--select 'ylm*_0_*']
                          [--modes 3] [--maps 4] [--scan-plot l]

For every ne{NE}-sc{SC}-f0-curlpsi[-ylm...] case in DIR (and its f2 partner):
  * test residual ||L x - y||/||y|| against sweqx's own L x (Ltest)
  * ||S+S^T||/||S||, S = M^1/2 L M^-1/2 (zero for a skew operator)
  * dense eigenvalues: max Re, min Re (= -max Re of the reversed flow -u), number with
    Re > 1e-8 omega, max|lambda| and the RK4 limit 2.83/max|lambda|;
    omega = u_max/a = u0/a (every ylm field is rescaled to max|u| = u0 at the GLL nodes)
  * velocity-gradient scale of the field: G = max ||grad u|| and Gs = max ||sym grad u||
    (2-norms of the tangential gradient, from the analytic psi; G = omega, Gs = 0 for
    solid body)
  * for the top --modes growing modes of f0 (shift-invert eigenvectors): energy
    enrichment (energy fraction / area fraction) near cube corners, along cube edges,
    on element-edge GLL nodes, and where |u| < 0.2 u_max (near stagnation points)
Results are cached in OUT/cache; table in OUT/results_flows.txt; figures in OUT.
"""
import argparse
import glob
import os
import re
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from opmat_io import load, read_meta
from opmat_eigs import Cache, mass_sym, skew_sym_err, cube_distances, SQRT8, INK, MUTED, GRID
from opmat_sym import xyz
from opmat_proj_plot import cube_edge_lines

GROW = 1e-8           # growing: Re lambda > GROW * omega
SLOW = 0.2            # low-speed region: |u| < SLOW * u_max


# ---------------------------------------------------------------- the velocity field
def ylm(l, m, lat, lon):
    """real Y_l^m with unit mean square on the sphere (as tc1_ylm in shallow_water_mod.F90)"""
    ma = abs(m)
    x, c = np.sin(lat), np.cos(lat)
    pmm = np.ones_like(x)
    if ma >= 1:
        pmm = np.sqrt(3.0) * c
    for k in range(2, ma + 1):
        pmm = np.sqrt((2 * k + 1) / (2.0 * k)) * c * pmm
    if l == ma:
        pl = pmm
    else:
        p2, p1 = pmm, np.sqrt(2.0 * ma + 3) * x * pmm
        for k in range(ma + 2, l + 1):
            a = np.sqrt((2.0 * k - 1) * (2 * k + 1) / ((k - ma) * (k + ma)))
            b = np.sqrt((2.0 * k + 1) * (k + ma - 1) * (k - ma - 1) / ((k - ma) * (k + ma) * (2 * k - 3.0)))
            p2, p1 = p1, a * x * p1 - b * p2
        pl = p1
    return pl * (np.cos(m * lon) if m > 0 else np.sin(ma * lon) if m < 0 else 1.0)


def psi_func(meta):
    """psi(lat, lon) / (a u0), before the rescaling to max|u| = u0"""
    al = meta['alpha']
    sb = lambda la, lo: np.sin(la) * np.cos(al) - np.cos(lo) * np.cos(la) * np.sin(al)
    if meta.get('psi', 'solidbody') != 'ylm':
        return sb
    l, m, eps = meta['psi_l'], meta['psi_m'], meta['psi_eps']
    return lambda la, lo: sb(la, lo) + eps * ylm(l, m, la, lo)


def velocity_3d(meta, X):
    """u (m/s) as 3-D vectors at unit vectors X: u = scale * u0 * grad(psi/(a u0)) x rhat
    (HOMME curl_sphere convention), gradient by centered differences on the sphere"""
    f = psi_func(meta)
    scale = meta.get('psi_uscale', 1.0) if meta.get('psi', 'solidbody') == 'ylm' else 1.0
    h = 1e-5
    lat = np.arcsin(np.clip(X[:, 2], -1, 1))
    lon = np.arctan2(X[:, 1], X[:, 0])
    # tangent basis (east, north); finite differences along great circles
    east = np.stack([-np.sin(lon), np.cos(lon), 0 * lon], axis=1)
    north = np.stack([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)], axis=1)

    def fx(Y):
        return f(np.arcsin(np.clip(Y[:, 2], -1, 1)), np.arctan2(Y[:, 1], Y[:, 0]))

    def walk(t, s):
        return np.cos(s) * X + np.sin(s) * t

    de = (fx(walk(east, h)) - fx(walk(east, -h))) / (2 * h)
    dn = (fx(walk(north, h)) - fx(walk(north, -h))) / (2 * h)
    grad = de[:, None] * east + dn[:, None] * north          # grad psi / u0 (unit sphere)
    return scale * meta['u0'] * np.cross(grad, X)


def gradient_scale(meta, X):
    """max over X of ||grad u|| and ||sym grad u|| (tangential 2x2 gradient, 1/s)"""
    a = meta['rearth']
    lat = np.arcsin(np.clip(X[:, 2], -1, 1))
    lon = np.arctan2(X[:, 1], X[:, 0])
    east = np.stack([-np.sin(lon), np.cos(lon), 0 * lon], axis=1)
    north = np.stack([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)], axis=1)
    h = 1e-4
    J = np.zeros((len(X), 2, 2))
    for j, t in enumerate((east, north)):
        dU = (velocity_3d(meta, np.cos(h) * X + np.sin(h) * t) -
              velocity_3d(meta, np.cos(h) * X - np.sin(h) * t)) / (2 * h * a)
        J[:, 0, j] = (dU * east).sum(1)
        J[:, 1, j] = (dU * north).sum(1)
    G = np.linalg.norm(J, ord=2, axis=(1, 2))
    Gs = np.linalg.norm(0.5 * (J + J.transpose(0, 2, 1)), ord=2, axis=(1, 2))
    return G.max(), Gs.max()


# ---------------------------------------------------------------- grid regions
def element_edge_nodes(lat, lon, ne):
    """True for GLL nodes on element boundaries (equiangular coordinate on a grid line)"""
    X = xyz(lat, lon)
    k = np.argmax(np.abs(X), axis=1)
    xf = np.abs(X[np.arange(len(X)), k])
    o1, o2 = X[np.arange(len(X)), (k + 1) % 3], X[np.arange(len(X)), (k + 2) % 3]
    on = np.zeros(len(X), bool)
    for c in (np.arctan(o1 / xf), np.arctan(o2 / xf)):
        s = (c + np.pi / 4) / (np.pi / 2 / ne)
        on |= np.abs(s - np.round(s)) < 1e-8
    return on


def regions(d, ne, speed):
    dc, de = cube_distances(d['lat'], d['lon'])
    he = np.pi / 2 / ne
    corner = dc < he
    return {'corner': corner, 'cubeedge': (de < 0.5 * he) & ~corner,
            'elemedge': element_edge_nodes(d['lat'], d['lon'], ne), 'slow': speed < SLOW * speed.max()}


def enrichment(v, mass, reg):
    w = mass * np.abs(v) ** 2
    w /= w.sum()
    return {k: w[r].sum() / (mass[r].sum() / mass.sum()) if r.any() else np.nan for k, r in reg.items()}


# ---------------------------------------------------------------- one case
def parse(nm):
    m = re.match(r'ne(\d+)-sc(\d+)-f(\d)-(\w+?)(?:-ylm(-?\d+)_(-?\d+)_(\S+))?$', nm)
    return m


def r_of(meta):
    if meta.get('psi', 'solidbody') != 'ylm':
        return 0.0
    l = meta['psi_l']
    return meta['psi_eps'] * np.sqrt(1.5 * l * (l + 1))


def analyze(cache, path, nmodes):
    nm = os.path.basename(path)

    def compute():
        d = load(path, ops=('L',))
        meta, L, mass = d['meta'], d['L'], d['mass']
        om = meta['u0'] / meta['rearth']
        x, y = d['Ltest']
        tres = max(np.linalg.norm(L @ x[:, k] - y[:, k]) / np.linalg.norm(y[:, k]) for k in range(x.shape[1]))
        sk, _ = skew_sym_err(L, mass)
        lam = np.linalg.eigvals(L.toarray())
        X = xyz(d['lat'], d['lon'])
        U = velocity_3d(meta, X)
        speed = np.linalg.norm(U, axis=1)
        G, Gs = gradient_scale(meta, X)
        r = dict(meta=np.array(meta, dtype=object), lam=lam, tres=tres, skew=sk, om=om, G=G, Gs=Gs, umax=speed.max(), rearth=meta['rearth'],
                 lat=d['lat'], lon=d['lon'], mass=mass, speed=speed)
        # top growing modes (one of each conjugate pair), eigenvectors by shift-invert
        o = np.argsort(-lam.real)
        sel = [j for j in o if lam[j].imag >= 0][:nmodes]
        vecs, lams = [], []
        for j in sel:
            if lam[j].real <= GROW * om:
                break
            ev, V = spla.eigs(L.astype(complex), k=1, sigma=lam[j] + 1e-9 * om * (1 + 1j))
            vecs.append(V[:, 0])
            lams.append(ev[0])
        r['vlam'] = np.array(lams)
        r['V'] = np.array(vecs).T if vecs else np.zeros((len(mass), 0), complex)
        return r

    return cache.get(nm, compute)


def summarize(r, ne):
    lam, om = r['lam'], float(r['om'])
    out = dict(maxre=lam.real.max() / om, minre=lam.real.min() / om,
               ngrow=int((lam.real > GROW * om).sum()), maxabs=np.abs(lam).max() / om,
               dt=SQRT8 / np.abs(lam).max(), skew=float(r['skew']), tres=float(r['tres']),
               G=float(r['G']) / om, Gs=float(r['Gs']) / om, umax=float(r['umax']) / (om * r['rearth']))
    reg = regions(r, ne, r['speed'])
    out['enr'] = [enrichment(r['V'][:, k], r['mass'], reg) for k in range(r['V'].shape[1])]
    out['vlam'] = r['vlam'] / om
    return out


# ---------------------------------------------------------------- figures
def style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def plot_scan(rows, key, out, fname):
    """max Re lambda/omega (f0, f2), and normalized by the velocity gradient, vs `key`"""
    xs = np.array([row[key] for row in rows])
    o = np.argsort(xs)
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4))
    f0 = np.array([row['f0']['maxre'] for row in rows])[o]
    f0m = np.array([-row['f0']['minre'] for row in rows])[o]
    f2 = np.array([row['f2']['maxre'] if row.get('f2') else np.nan for row in rows])[o]
    G = np.array([row['f0']['G'] for row in rows])[o]
    Gs = np.array([row['f0']['Gs'] for row in rows])[o]
    ax = axs[0]
    ax.plot(xs[o], f0, 'o-', color='#1baf7a', label='conservation, u')
    ax.plot(xs[o], f0m, 's--', color='#6250d6', mfc='none', label='conservation, -u')
    ax.plot(xs[o], np.abs(f2), '^-', color='#2a78d6', label='split |max Re|')
    ax.set_yscale('log')
    ax.set_xlabel(key)
    ax.set_ylabel('max Re λ / (u_max/a)')
    ax.legend(fontsize=8)
    style(ax)
    ax = axs[1]
    ax.plot(xs[o], G, 'o-', color=INK, label='max ||∇u|| / (u_max/a)')
    ax.plot(xs[o], Gs, 's-', color=MUTED, label='max ||sym ∇u|| / (u_max/a)')
    ax.set_xlabel(key)
    ax.legend(fontsize=8)
    style(ax)
    ax = axs[2]
    ax.plot(xs[o], f0 / G, 'o-', color='#1baf7a', label='max Re λ / max||∇u||')
    ax.plot(xs[o], f0 / np.where(Gs > 1e-3, Gs, np.nan), 's-', color='#eb6834', label='max Re λ / max||sym ∇u||')
    ax.set_xlabel(key)
    ax.legend(fontsize=8)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out, fname))
    plt.close(fig)


def plot_spectra(cases, out, fname):
    n = len(cases)
    fig, axs = plt.subplots(1, n, figsize=(3.2 * n, 3.2), squeeze=False)
    for ax, (title, r) in zip(axs[0], cases):
        om = float(r['om'])
        lam = r['lam'] / om
        ax.plot(lam.real, lam.imag, '.', ms=2, color='#1baf7a')
        ax.set_title(title, fontsize=8)
        ax.set_xlabel('Re λ / ω')
        style(ax)
    axs[0][0].set_ylabel('Im λ / ω')
    fig.tight_layout()
    fig.savefig(os.path.join(out, fname))
    plt.close(fig)


def plot_modes(cases, out, fname, nmodes):
    """|v|^2 of the top growing modes, lon-lat, with the cube edges; contours of |u|"""
    edges = cube_edge_lines()
    fig, axs = plt.subplots(len(cases), nmodes, figsize=(4.2 * nmodes, 2.3 * len(cases)), squeeze=False)
    for row, (title, r) in enumerate(cases):
        om = float(r['om'])
        lon, lat = np.degrees(r['lon']) % 360, np.degrees(r['lat'])
        # |u|/u_max of the analytic field on a lon-lat grid, for the low-speed contour
        glo, gla = np.meshgrid(np.linspace(0, 360, 361), np.linspace(-89.5, 89.5, 180))
        sp_ = np.linalg.norm(velocity_3d(r['meta'][()], xyz(np.radians(gla.ravel()), np.radians(glo.ravel()))), axis=1)
        sg = (glo, gla, sp_.reshape(glo.shape) / r['speed'].max())
        for k in range(nmodes):
            ax = axs[row][k]
            for lo, la in edges:
                ax.plot(lo, la, color='#2a78d6', lw=0.5)
            ax.set_xlim(0, 360)
            ax.set_ylim(-90, 90)
            ax.set_xticks([])
            ax.set_yticks([])
            if k >= r['V'].shape[1]:
                continue
            e = np.abs(r['V'][:, k]) ** 2
            e /= e.max()
            o = np.argsort(e)
            ax.scatter(lon[o], lat[o], c=e[o], s=2 + 14 * e[o], cmap='magma_r', vmin=0, vmax=1, linewidths=0)
            ax.contour(sg[0], sg[1], sg[2], levels=[SLOW], colors=MUTED, linewidths=0.6)
            lam = r['vlam'][k] / om
            ax.set_title(f'{title}  λ/ω = {lam.real:.3f}{lam.imag:+.2f}i', fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, fname))
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--mats', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--ne', type=int, nargs='+', default=[8])
    p.add_argument('--sc', type=int, default=1)
    p.add_argument('--select', default='*', help="glob on the psi tag, e.g. 'ylm*_0_*'")
    p.add_argument('--baseline', action='store_true', help='include the solid-body curlpsi case')
    p.add_argument('--modes', type=int, default=3)
    p.add_argument('--maps', type=int, default=4, help='mode maps for this many cases (largest max Re)')
    p.add_argument('--scan-plot', default='l', choices=['l', 'm', 'r'])
    p.add_argument('--tag', default='scan')
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    cache = Cache(args.out)
    rep = []
    for ne in args.ne:
        pats = sorted(glob.glob(os.path.join(args.mats, f'ne{ne}-sc{args.sc}-f0-curlpsi-{args.select}.meta')))
        if args.baseline and os.path.exists(os.path.join(args.mats, f'ne{ne}-sc{args.sc}-f0-curlpsi.meta')):
            pats = [os.path.join(args.mats, f'ne{ne}-sc{args.sc}-f0-curlpsi.meta')] + pats
        rows, res = [], {}
        for mf in pats:
            p0 = mf[:-5]
            nm = os.path.basename(p0)
            meta = read_meta(p0)
            row = dict(name=nm, l=meta.get('psi_l', 1), m=meta.get('psi_m', 0), r=r_of(meta),
                       eps=meta.get('psi_eps', 0.0))
            print(nm, flush=True)
            r0 = analyze(cache, p0, args.modes)
            row['f0'] = summarize(r0, ne)
            res[nm] = r0
            p2 = p0.replace('-f0-', '-f2-')
            if os.path.exists(p2 + '.meta'):
                row['f2'] = summarize(analyze(cache, p2, 0), ne)
            rows.append(row)
        if not rows:
            continue
        rep.append(f'NE{ne} sc{args.sc}: max Re, min Re, max|lambda| in units of omega = u_max/a; '
                   f'G, Gs = max||grad u||, max||sym grad u|| / omega; ngrow = #(Re > {GROW:g} omega); '
                   f'skew = ||S+S^T||/||S||; test = ||L x - y||/||y||; umax = max|u| of the analytic psi at the nodes / u0')
        rep.append(f'{"case":38s} {"l":>3s} {"m":>3s} {"r":>5s} {"G":>6s} {"Gs":>6s} | {"f0 maxRe":>9s} '
                   f'{"minRe":>9s} {"ngrow":>5s} {"max|l|":>7s} {"skew":>8s} {"test":>7s} {"umax":>6s} | {"f2 maxRe":>9s} '
                   f'{"skew":>8s} | top f0 modes: lam/omega [enrichment corner cubeedge elemedge slow]')
        for row in rows:
            a, b = row['f0'], row.get('f2')
            s = (f'{row["name"]:38s} {row["l"]:3d} {row["m"]:3d} {row["r"]:5.2f} {a["G"]:6.2f} {a["Gs"]:6.2f} | '
                 f'{a["maxre"]:9.2e} {a["minre"]:9.2e} {a["ngrow"]:5d} {a["maxabs"]:7.2f} {a["skew"]:8.2e} '
                 f'{a["tres"]:7.1e} {a["umax"]:6.4f} | ')
            s += f'{b["maxre"]:9.2e} {b["skew"]:8.2e} |' if b else f'{"-":>9s} {"-":>8s} |'
            for lam, e in zip(a['vlam'], a['enr']):
                s += (f' {lam.real:.3f}{lam.imag:+.2f}i [{e["corner"]:.1f} {e["cubeedge"]:.1f} '
                      f'{e["elemedge"]:.1f} {e["slow"]:.1f}]')
            rep.append(s)
        rep.append('')
        tag = f'ne{ne}-sc{args.sc}-{args.tag}'
        plot_scan(rows, args.scan_plot, args.out, f'{tag}-maxre.png')
        worst = sorted(rows, key=lambda row: -row['f0']['maxre'])[:args.maps]
        sel = sorted(worst, key=lambda row: row[args.scan_plot])
        plot_modes([(f'{row["name"].split("curlpsi-")[-1]}', res[row['name']]) for row in sel],
                   args.out, f'{tag}-modes.png', args.modes)
        spec = sorted(rows, key=lambda row: row[args.scan_plot])
        spec = [spec[i] for i in np.unique(np.linspace(0, len(spec) - 1, min(4, len(spec))).astype(int))]
        plot_spectra([(row['name'].split('curlpsi-')[-1], res[row['name']]) for row in spec],
                     args.out, f'{tag}-spectra.png')
    txt = '\n'.join(rep)
    print(txt)
    with open(os.path.join(args.out, f'results_flows-{args.tag}.txt'), 'w') as f:
        f.write(txt + '\n')


if __name__ == '__main__':
    main()
