#!/usr/bin/env python3
"""
Where is the swtc1 Gaussian's error generated?  The semi-discrete error
e = p - p_ex (p_ex the exact solution at the GLL nodes) obeys

    de/dt = L e + tau(t),   tau = L p_ex - d p_ex/dt,   e(0) = 0,

so e(T) = sum_r e_r(T), where e_r is driven by tau restricted to region r.
Regions (h = element width pi/(2 NE)):
  corner : within 1.5 h of a cube corner
  edge   : GLL nodes on a cube edge (outside the corner regions)
  layer  : within one element width of a cube edge, not on it
  inner  : the rest
Each e_r is integrated with RK4 at dt = dt30*30/NE/sub.

usage: opmat_source.py --mats DIR --out OUT --ne 16 20 30 60 [--sc 1 4] [--sub 4]
"""
import argparse
import os
import time
import numpy as np

from opmat_io import load
from opmat_sym import xyz
from opmat_eigs import name, omega, cube_distances
from opmat_proj import RR, LATC, LONC, norms

REGIONS = ('corner', 'edge', 'layer', 'inner')


class Exact:
    """p_ex(x, t) = F(c(t).x), the Gaussian centered at c(t) = R(omega t) c0 about k"""
    def __init__(self, X, alpha, om):
        self.X, self.om = X, om
        self.k = np.array([-np.sin(alpha), 0, np.cos(alpha)])
        self.c0 = np.array([np.cos(LATC) * np.cos(LONC), np.cos(LATC) * np.sin(LONC), np.sin(LATC)])

    def center(self, t):
        th = self.om * t
        k, c = self.k, self.c0
        return c * np.cos(th) + np.cross(k, c) * np.sin(th) + k * (k @ c) * (1 - np.cos(th))

    def __call__(self, t):
        """(p_ex, d p_ex/dt) at the nodes"""
        c = self.center(t)
        s = np.clip(self.X @ c, -1, 1)
        r = np.arccos(s)
        F = 1000.0 * np.exp(-(r / RR) ** 2)
        sinr = np.sqrt(np.maximum(1 - s * s, 0))
        # dF/ds = F 2 r/(RR^2 sin r); r/sin r -> 1 at the center, F ~ 1e-36 at the antipode
        ratio = np.where(sinr > 1e-12, r / np.maximum(sinr, 1e-12), np.where(r < 1, 1.0, 0.0))
        Fp = F * 2 / RR ** 2 * ratio
        return F, self.om * Fp * (self.X @ np.cross(self.k, c))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mats', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--ne', type=int, nargs='+', required=True)
    p.add_argument('--sc', type=int, nargs='+', default=[1, 4])
    p.add_argument('--dt30', type=float, default=1800.0)
    p.add_argument('--sub', type=int, default=4)
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rep = []
    for ne in args.ne:
        for sc in args.sc:
            nm = name(ne, sc, 2, 'analytic')
            d = load(os.path.join(args.mats, nm), ops=('L',))
            L, mass = d['L'], d['mass']
            om = omega(d['meta'])
            T = 2 * np.pi / om
            nstep = int(round(T / (args.dt30 * 30 / ne))) * args.sub
            dt = T / nstep
            X = xyz(d['lat'], d['lon'])
            ex = Exact(X, d['meta']['alpha'], om)
            h = (np.pi / 2) / ne
            dc, de = cube_distances(d['lat'], d['lon'])
            corner = dc < 1.5 * h
            masks = np.stack([corner, ~corner & (de < 1e-9), ~corner & (de >= 1e-9) & (de < h),
                              ~corner & (de >= h)], axis=1).astype(float)
            assert np.all(masks.sum(axis=1) == 1)

            def rhs(E, t):
                pe, dpe = ex(t)
                tau = L @ pe - dpe
                return L @ E + masks * tau[:, None]

            t0 = time.time()
            E = np.zeros((d['n'], len(REGIONS)))
            for i in range(nstep):
                t = i * dt
                k1 = rhs(E, t)
                k2 = rhs(E + 0.5 * dt * k1, t + 0.5 * dt)
                k3 = rhs(E + 0.5 * dt * k2, t + 0.5 * dt)
                k4 = rhs(E + dt * k3, t + dt)
                E += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            p0 = ex(0.0)[0]
            etot = E.sum(axis=1)
            line = (f'{nm} dt={dt:.0f}: total l2={norms(etot, p0, mass)[0]:.5e} max={norms(etot, p0, mass)[1]:.5e} | '
                    + ' | '.join(f'{r} (area {masks[:, j] @ mass / mass.sum():.4f}) l2={norms(E[:, j], p0, mass)[0]:.4e} '
                                 f'max={norms(E[:, j], p0, mass)[1]:.4e}' for j, r in enumerate(REGIONS))
                    + f'  ({time.time() - t0:.0f} s)')
            rep.append(line)
            print(line, flush=True)
            np.savez(os.path.join(args.out, nm + '-source.npz'), E=E, masks=masks, p0=p0, mass=mass,
                     lat=d['lat'], lon=d['lon'], regions=np.array(REGIONS))
    with open(os.path.join(args.out, 'source.txt'), 'a') as f:
        f.write('\n'.join(rep) + '\n')


if __name__ == '__main__':
    main()
