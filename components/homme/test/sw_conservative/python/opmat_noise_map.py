"""
Regional maps of a curlpsi velocity field from an opmat extraction: stream function,
speed with velocity arrows, DSS'd vorticity, and (optionally) the growing eigenmode of
L - nu B that sets the spectral hyperviscosity threshold.  Element and cube edges are drawn.

    python opmat_noise_map.py --prefix DIR/ne30-sc1-f0-curlpsi-noise-ne30-lcinf-s2-w30_-35_6_3_bg0 \
        --center 30 -35 --radius 14 [--mode 2.6e15] [--mats DIR] [--out map.png]

The stream function is read from the meta's psi_file, in physical units
psi_uscale*rearth*u0*psi_file (m^2/s), the same scaling sweqx applied.  --mode NU30
computes the rightmost eigenvector of L - nu B at that nu30 (just below the threshold;
needs B in --mats, see opmat_flows_hv.load_B) by shift-invert near --shift (x omega).
Projection: local azimuthal equidistant about --center (distances in degrees of arc).
"""
import argparse
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm

from opmat_io import load
from opmat_sym import xyz
from opmat_flows_hv import load_B, nu_from30


def project(X, c):
    """azimuthal equidistant about unit vector c: (east, north) in degrees of arc"""
    z = np.array([0, 0, 1.0])
    e = np.cross(z, c); e /= np.linalg.norm(e)
    n = np.cross(c, e)
    d = np.arccos(np.clip(X @ c, -1, 1))
    ang = np.arctan2(X @ n, X @ e)
    return np.degrees(d) * np.cos(ang), np.degrees(d) * np.sin(ang)


def element_lines(ne, c, rmax, npts=60):
    """great-circle segments of the element edges (equiangular cube) within rmax of c"""
    t = np.linspace(-np.pi / 4, np.pi / 4, npts * ne // 10 + 2)
    ks = -np.pi / 4 + np.arange(ne + 1) * np.pi / (2 * ne)
    lines = []
    for face in range(6):
        ax, sg = face % 3, (1 if face < 3 else -1)
        o1, o2 = [i for i in range(3) if i != ax]
        for k in ks:
            for which in (0, 1):
                P = np.zeros((len(t), 3))
                P[:, ax] = sg
                P[:, o1] = np.tan(k) if which == 0 else np.tan(t)
                P[:, o2] = np.tan(t) if which == 0 else np.tan(k)
                P /= np.linalg.norm(P, axis=1)[:, None]
                if np.degrees(np.arccos(np.clip(P @ c, -1, 1))).min() > rmax * 1.5:
                    continue
                lines.append((P, abs(abs(k) - np.pi / 4) < 1e-12))
    return lines


def threshold_mode(d, mats, n30, shift):
    L, ne, om = d['L'], d['meta']['ne'], d['meta']['u0'] / d['meta']['rearth']
    B = load_B(mats, ne, d['lat'], d['lon'])
    A = (L - nu_from30(n30, ne) * B).tocsc()
    n, sg = A.shape[0], shift * om
    lu = spla.splu((A - sg * sp.identity(n, format='csc')).tocsc(), permc_spec='MMD_AT_PLUS_A',
                   diag_pivot_thresh=0.1)
    op = spla.LinearOperator((n, n), matvec=lu.solve, dtype=float)
    lam, V = spla.eigs(A, k=16, sigma=sg, OPinv=op, which='LM', tol=1e-8, maxiter=300)
    re = np.where(np.abs(lam) < 1e-8 * om, -np.inf, lam.real)   # skip the constant (null) mode
    j = re.argmax()
    return lam[j] / om, V[:, j]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--prefix', required=True)
    p.add_argument('--center', type=float, nargs=2, required=True, metavar=('LAT', 'LON'))
    p.add_argument('--radius', type=float, default=14.0, help='map half-width (degrees of arc)')
    p.add_argument('--zoom', type=float, default=5.0, help='half-width of the arrow inset (degrees)')
    p.add_argument('--mode', type=float, help='nu30 for the threshold eigenmode (omit to skip)')
    p.add_argument('--shift', type=float, default=0.05, help='shift-invert target for --mode, x omega')
    p.add_argument('--mats', help='directory with ne<NE>-sc1-f2-analytic B files (for --mode)')
    p.add_argument('--out', default='flow-map.png')
    a = p.parse_args()

    d = load(a.prefix, ops=('L',) if a.mode else ())
    meta, ne = d['meta'], d['meta']['ne']
    om = meta['u0'] / meta['rearth']
    rec = np.fromfile(meta['psi_file'], dtype=[('g', '<i8'), ('p', '<f8')])
    g2i = {g: i for i, g in enumerate(rec['g'])}
    psi = rec['p'][[g2i[g] for g in d['gdof']]] * meta['psi_uscale'] * meta['rearth'] * meta['u0']
    X = xyz(d['lat'], d['lon'])
    c = xyz(np.radians([a.center[0]]), np.radians([a.center[1]]))[0]
    x, y = project(X, c)
    sel = (np.abs(x) < a.radius * 1.1) & (np.abs(y) < a.radius * 1.1) & (X @ c > 0)
    tri = mtri.Triangulation(x[sel], y[sel])
    # unit east/north at each node, rotated into the map frame (exact at the centre)
    lat, lon = d['lat'], d['lon']
    east = np.stack([-np.sin(lon), np.cos(lon), 0 * lon], 1)
    north = np.stack([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)], 1)
    V3 = d['vel'][:, :1] * east + d['vel'][:, 1:] * north
    h = 1e-3
    xp, yp = project(X + h * V3 / np.maximum(np.linalg.norm(V3, axis=1), 1e-30)[:, None], c)
    sp_ = np.linalg.norm(d['vel'], axis=1)
    ux, uy = (xp - x) / h * sp_ / np.maximum(np.hypot((xp - x) / h, (yp - y) / h), 1e-30), \
             (yp - y) / h * sp_ / np.maximum(np.hypot((xp - x) / h, (yp - y) / h), 1e-30)

    panels = 4 if a.mode else 3
    fig, axs = plt.subplots(1, panels, figsize=(5.2 * panels, 5.0))
    lines = element_lines(ne, c, a.radius)

    def frame(ax, r, title):
        for P, cube in lines:
            px, py = project(P, c)
            ok = P @ c > 0
            ax.plot(np.where(ok, px, np.nan), np.where(ok, py, np.nan),
                    color='k' if cube else '0.55', lw=1.4 if cube else 0.4, zorder=3)
        ax.set_xlim(-r, r); ax.set_ylim(-r, r); ax.set_aspect('equal')
        ax.set_xlabel('east (deg of arc)'); ax.set_ylabel('north (deg of arc)'); ax.set_title(title, fontsize=10)

    lim = np.abs(psi[sel]).max()
    im = axs[0].tripcolor(tri, psi[sel] / 1e6, cmap='RdBu_r', vmin=-lim / 1e6, vmax=lim / 1e6, shading='gouraud')
    plt.colorbar(im, ax=axs[0], shrink=0.8, label=r'$\psi$ (10$^6$ m$^2$/s)')
    frame(axs[0], a.radius, 'stream function (nodal values, linear between nodes)')

    im = axs[1].tripcolor(tri, sp_[sel], cmap='viridis', shading='gouraud')
    plt.colorbar(im, ax=axs[1], shrink=0.8, label='|u| (m/s)')
    zsel = sel & (np.abs(x) < a.zoom) & (np.abs(y) < a.zoom)
    frame(axs[1], a.zoom, f'speed and velocity at the GLL nodes (max {sp_.max():.1f} m/s)')
    axs[1].quiver(x[zsel], y[zsel], ux[zsel], uy[zsel], color='w', scale=25 * sp_.max() / 2, width=0.004, zorder=4)

    z = d['zeta'] / om
    lz = np.percentile(np.abs(z[sel]), 99.5)
    im = axs[2].tripcolor(tri, z[sel], cmap='PuOr_r', vmin=-lz, vmax=lz, shading='gouraud')
    plt.colorbar(im, ax=axs[2], shrink=0.8, label=r'$\zeta/\omega$')
    frame(axs[2], a.radius, rf'DSS vorticity (max $|\zeta|$ = {np.abs(z).max():.0f} $\omega$)')

    if a.mode:
        lam, v = threshold_mode(d, a.mats, a.mode, a.shift)
        e = d['mass'] * np.abs(v) ** 2
        e /= e.max()
        im = axs[3].tripcolor(tri, np.maximum(e[sel], 1e-4), cmap='magma_r', norm=LogNorm(1e-4, 1), shading='gouraud')
        plt.colorbar(im, ax=axs[3], shrink=0.8, label='mode energy (rel. to max)')
        frame(axs[3], a.zoom, rf'growing mode of L$-\nu$B, $\nu_{{30}}$={a.mode:.2e}: Re$\lambda$={lam.real:.3f}$\omega$')
        print(f'mode: lambda/omega = {lam:.4f}')
    for ax in axs[:3:2]:
        for r in (6, 9):
            t = np.linspace(0, 2 * np.pi, 200)
            ax.plot(r * np.cos(t), r * np.sin(t), ':', color='0.2', lw=0.8, zorder=4)
    fig.suptitle(f'{a.prefix.split("/")[-1]}   (centre {a.center[0]:g}N {a.center[1]:g}E; dotted: '
                 f'window flat-top and taper radii, for patches)', fontsize=10)
    plt.tight_layout()
    plt.savefig(a.out, dpi=120)
    print(a.out)


if __name__ == '__main__':
    main()
