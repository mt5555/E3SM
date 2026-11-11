#!/usr/bin/env python3
#
#   Input:
#   transport matrix   A
#   nodal coordinates: lat,lon
#   weights:           M
#
#   step 1:  compute spherical harmonics, orthonormalize with QR decomp
#            at nodal values given by lat,lon.
#            verify Q^T M Q = I  (orthonormal)
#   step 2:   P = expm_multiply(A * T, Q)
#   step 3:   compute max svd of   M^1/2 (P - Q)
#
#
#   for NE15:
#     wave length     degree N
#       2 dx            90            should be horrible
#       4 dx            45            modes above this are horrible
#       8 dx            22            phase error ~1%, so transport
#                                     error will be large after 12 days
#      16 dx            11            very well resolved, expect small errors
#
#
"""
Worst-case error per spherical harmonic degree N after time T, in nodal variables:

    E(N) = max over p0 in span(Y_N) of ||(exp(AT) - I) p0||_M / ||p0||_M
         = sigma_max( M^1/2 (exp(AT) - I) Q_N ),   with Q_N^T M Q_N = I

A minimal, one-degree-at-a-time version of opmat_shnorm.py (which propagates all degrees
as one block and adds plots); both give the same E(N) to roundoff.  In the table above
dx = 2 pi a/(12 NE) (the GLL spacing along the equator) and N = 2 pi a / wavelength.

Inputs: A sparse n x n (dp/dt = A p, 1/s), mass the assembled diagonal mass (GLL weights,
any units), lat/lon in radians, T in seconds (one swtc1 revolution: 12 days).

usage: opmat_shnorm_simple.py --prefix MATS/ne10-sc1-f2-analytic [--nmax 30] [--ndays 12]
                              [--ref OUT/<case>-...-shnorm.npz]
"""
import argparse
import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import expm_multiply
from scipy.special import sph_harm_y        # scipy >= 1.15; older: sph_harm(m, n, azimuth, colatitude)


def real_harmonics(lat, lon, N):
    """the 2N+1 real spherical harmonics of degree N at the nodes (n x (2N+1)), orthonormal
    on the unit sphere.  The pole can be anywhere: span(Y_N) is the same for every pole,
    so E(N) does not depend on it."""
    colat, az = np.pi / 2 - lat, lon
    cols = []
    for m in range(-N, N + 1):
        y = sph_harm_y(N, abs(m), colat, az)            # complex Y_N^|m|
        if m == 0:
            cols.append(y.real)
        else:
            cols.append(np.sqrt(2) * (y.real if m > 0 else y.imag))
    return np.stack(cols, axis=1)


def orthonormalize(Y, mass):
    """step 1: Q = Y R^-1 from the QR factorization M^1/2 Y = Qh R, so that Q^T M Q = I.
    Returns Q, R and max |Q^T M Q - I|."""
    sq = np.sqrt(mass)
    Qh, R = np.linalg.qr(sq[:, None] * Y)               # reduced QR: Qh is n x k, R is k x k
    Q = Qh / sq[:, None]                                # back to nodal values: Q = M^-1/2 Qh
    G = Q.T @ (mass[:, None] * Q)                       # Q^T M Q, should be the identity
    return Q, R, np.abs(G - np.eye(G.shape[0])).max()


def degree_error(A, Q, mass, T):
    """steps 2 and 3: P = exp(AT) Q, then E(N) = sigma_max(M^1/2 (P - Q)).  Also returns the
    largest singular values of the part of the error that stays in degree N and of the part
    scattered to other degrees."""
    sq = np.sqrt(mass)
    P = expm_multiply(A * T, Q)                         # step 2: action of exp(AT), never formed
    X = sq[:, None] * (P - Q)                           # step 3: error in mass-weighted coordinates
    E = la.svdvals(X)[0]
    Qh = sq[:, None] * Q                                # Euclidean-orthonormal version of Q
    Xin = Qh @ (Qh.T @ X)                               # projection onto degree N
    return E, la.svdvals(Xin)[0], la.svdvals(X - Xin)[0]


def shnorm(A, mass, lat, lon, T, Nmax):
    """[(N, E(N), in-degree, scattered, max|Q^T M Q - I|, cond(R))] for N = 1..Nmax"""
    out = []
    for N in range(1, Nmax + 1):
        Q, R, orth_err = orthonormalize(real_harmonics(lat, lon, N), mass)
        E, Ein, Eout = degree_error(A, Q, mass, T)
        out.append((N, E, Ein, Eout, orth_err, np.linalg.cond(R)))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--prefix', required=True, help='opmat prefix of L, e.g. MATS/ne10-sc1-f2-analytic')
    p.add_argument('--nmax', type=int, default=30)
    p.add_argument('--ndays', type=float, default=12.0)
    p.add_argument('--ref', help='opmat_shnorm.py npz to compare E(N) with')
    a = p.parse_args()
    from opmat_io import load
    d = load(a.prefix, ops=('L',))
    T = a.ndays * 86400.0
    ref = np.load(a.ref)['E'] if a.ref else None
    print(f'{a.prefix}: n = {d["n"]}, T = {a.ndays:g} days')
    print(' N   E(N)        in-degree   scattered   |Q^T M Q - I|  cond(R)' + ('   rel diff to --ref' if a.ref else ''))
    for N, E, Ein, Eout, oe, cR in shnorm(d['L'], d['mass'], d['lat'], d['lon'], T, a.nmax):
        extra = f'   {abs(E - ref[N]) / ref[N]:.1e}' if ref is not None else ''
        print(f'{N:2d}   {E:.4e}  {Ein:.3e}   {Eout:.3e}   {oe:.1e}        {cR:.3f}{extra}', flush=True)


if __name__ == '__main__':
    main()
