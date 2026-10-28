#!/usr/bin/env python3
"""
Error-minimizing hyperviscosity for the swtc1 Gaussian, from the extracted
operators: HOMME's classic-RK4 step with hyperviscosity (advance_nonstag_rk4,
then advance_hypervis with hypervis_subcycle=1) is

    p <- (I - dt nu B) R4(dt L) p

with B = rspheremp DSS(biharmonic_wk).  Many nu are run at once as the columns
of one block.  A coarse pass over nu30 = nu (NE/30)^4 is followed by a fine
pass around the minima of the relative l2 and max errors at day 12.

usage: opmat_nuopt.py --mats DIR --out OUT --ne 16 20 30 60 [--sc 4] [--dt30 1800]
B is read from ne{NE}-sc1-f2-analytic (it does not depend on the flow).
"""
import argparse
import os
import time
import numpy as np

from opmat_io import load
from opmat_eigs import name, omega
from opmat_proj import gaussian, norms, timestep


def run(L, B, p0, dt, nstep, nus):
    """day-T fields for each nu in nus (columns)"""
    X = np.repeat(p0[:, None], len(nus), axis=1)
    hv = dt * np.asarray(nus)[None, :]
    for _ in range(nstep):
        k1 = L @ X
        k2 = L @ (X + 0.5 * dt * k1)
        k3 = L @ (X + 0.5 * dt * k2)
        k4 = L @ (X + dt * k3)
        X += dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        X -= hv * (B @ X)
    return X


def errors(X, p0, mass):
    e = X - p0[:, None]
    l2 = np.sqrt((mass[:, None] * e ** 2).sum(axis=0) / (mass * p0 ** 2).sum())
    return l2, np.abs(e).max(axis=0) / np.abs(p0).max()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mats', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--ne', type=int, nargs='+', required=True)
    p.add_argument('--sc', type=int, nargs='+', default=[4])
    p.add_argument('--dt30', type=float, default=1800.0)
    p.add_argument('--coarse', type=float, nargs=3, default=[11.0, 14.5, 15],
                   help='log10 nu30 range and count of the coarse pass (0 is added)')
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for ne in args.ne:
        B = load(os.path.join(args.mats, name(ne, 1, 2, 'analytic')), ops=('B',))['B']
        for sc in args.sc:
            nm = name(ne, sc, 2, 'analytic')
            d = load(os.path.join(args.mats, nm), ops=('L',))
            om = omega(d['meta'])
            T = 2 * np.pi / om
            dt, nstep = timestep(ne, args.dt30, T)
            p0 = gaussian(d['lat'], d['lon'])
            scale = (30.0 / ne) ** 4
            t0 = time.time()
            nu30 = np.r_[0.0, np.logspace(*args.coarse[:2], int(args.coarse[2]))]
            l2, mx = errors(run(d['L'], B, p0, dt, nstep, nu30 * scale), p0, d['mass'])
            # fine pass: 8 values in log nu30 around each coarse minimum (skip nu=0)
            fine = []
            for err in (l2, mx):
                i = 1 + np.argmin(err[1:])
                lo, hi = np.log10(nu30[max(i - 1, 1)]), np.log10(nu30[min(i + 1, len(nu30) - 1)])
                fine.append(np.logspace(lo, hi, 9)[1:-1])
            nu30f = np.unique(np.concatenate(fine))
            l2f, mxf = errors(run(d['L'], B, p0, dt, nstep, nu30f * scale), p0, d['mass'])
            nu30 = np.r_[nu30, nu30f]
            l2, mx = np.r_[l2, l2f], np.r_[mx, mxf]
            o = np.argsort(nu30)
            nu30, l2, mx = nu30[o], l2[o], mx[o]
            np.savez(os.path.join(args.out, nm + '-nuopt.npz'), nu30=nu30, nu=nu30 * scale, l2=l2, max=mx,
                     ne=ne, sc=sc, dt=dt)
            i2, im = np.argmin(l2), np.argmin(mx)
            line = (f'{nm} dt={dt:.0f}: nu=0 l2={l2[0]:.4e} max={mx[0]:.4e} | l2-opt nu30={nu30[i2]:.3e} '
                    f'(nu={nu30[i2] * scale:.3e}) l2={l2[i2]:.4e} | max-opt nu30={nu30[im]:.3e} '
                    f'(nu={nu30[im] * scale:.3e}) max={mx[im]:.4e}  ({time.time() - t0:.0f} s)')
            print(line, flush=True)
            with open(os.path.join(args.out, 'nuopt.txt'), 'a') as f:
                f.write(line + '\n')


if __name__ == '__main__':
    main()
