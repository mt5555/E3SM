"""
Time integration of a tracer with the extracted operators, HOMME's way:

    q <- (I - (dt/s) nu B)^s R4(dt L) q

(RK4 for the transport, then s forward-Euler hyperviscosity subcycles, as advance_hypervis with
hypervis_subcycle=s; s is the smallest with (dt/s) nu lambda_max(B) <= 1.5).  Tracks the
mass-weighted variance  Var = sum m (q - qbar)^2 / sum m,  the mass  sum m q,  and max|q|;
q is rescaled whenever it grows past 1e100 (the log of the factor is kept), so a run of the
unstable conservation form does not overflow.

    python opmat_timeint.py --mats MATS --case noise-ne30-lcinf-s2-w30_-35_6_3_bg0 --ne 30 \
        --runs f0:0 f0:1e15 f0:2.5e15 f0:3e15 f2:0 f2:1e15 --days 24 --dt 1800 \
        --gauss 30 -35 6 --snap 0 0.25 0.5 1 2 6 12 --out OUT

Run <form>:<nu30> uses MATS/ne<NE>-sc1-f<form>-curlpsi-<case> and nu = nu30 (30/NE)^4, B from
MATS/ne<NE>-sc1-f2-analytic.  Initial tracer (--gauss LAT LON R): q0 = exp(-(d/R)^2), d the
great-circle distance in degrees.  Output OUT/timeint-<form>-nu<nu30>.npz: t (days), log10
of Var/Var0, mass drift, log10 max|q|, and q at the --snap times (with their log10 scale).
"""
import argparse
import os
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from opmat_io import load
from opmat_dense import load_B, nu_from30
from opmat_sym import xyz


def var(q, m):
    qb = (m * q).sum() / m.sum()
    return (m * (q - qb) ** 2).sum() / m.sum()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--mats', required=True)
    p.add_argument('--case', required=True)
    p.add_argument('--ne', type=int, default=30)
    p.add_argument('--runs', nargs='+', default=['f0:0', 'f2:0'])
    p.add_argument('--days', type=float, default=12.0)
    p.add_argument('--dt', type=float, default=1800.0)
    p.add_argument('--gauss', type=float, nargs=3, default=[30.0, -35.0, 6.0], metavar=('LAT', 'LON', 'R'))
    p.add_argument('--snap', type=float, nargs='*', default=[0, 1, 12], help='snapshot times (days)')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)

    B = lmax = None
    for run in a.runs:
        form, n30 = run.split(':')
        form = form.lstrip('f')
        n30 = float(n30)
        f = os.path.join(a.out, f'timeint-f{form}-nu{n30:.4g}.npz')
        if os.path.exists(f):
            print(f'{f} exists: skipped'); continue
        d = load(os.path.join(a.mats, f'ne{a.ne}-sc1-f{form}-curlpsi-{a.case}'), ops=('L',))
        L, m = d['L'].tocsr(), d['mass']
        nu = nu_from30(n30, a.ne)
        s = 0
        if nu > 0:
            if B is None:
                B = load_B(a.mats, a.ne, d['lat'], d['lon']).tocsr()
                lmax = np.abs(spla.eigs(B, k=1, which='LM', return_eigenvectors=False, maxiter=5000, tol=1e-4)).max()
                print(f'lambda_max(B) = {lmax:.4e}', flush=True)
            s = max(1, int(np.ceil(a.dt * nu * lmax / 1.5)))
        X = xyz(d['lat'], d['lon'])
        la, lo, R = a.gauss
        c = xyz(np.radians([la]), np.radians([lo]))[0]
        q = np.exp(-(np.degrees(np.arccos(np.clip(X @ c, -1, 1))) / R) ** 2)
        nstep = int(round(a.days * 86400 / a.dt))
        v0, mass0 = var(q, m), (m * q).sum()
        snaps = {int(round(t * 86400 / a.dt)): t for t in a.snap}
        logs = 0.0                                 # q_true = q * 10^logs
        T, LV, MD, LMX, SQ, ST, SS = [], [], [], [], [], [], []
        t0 = time.time()
        for k in range(nstep + 1):
            if k in snaps:
                SQ.append(q.copy()); ST.append(snaps[k]); SS.append(logs)
            T.append(k * a.dt / 86400)
            LV.append(np.log10(var(q, m) / v0) + 2 * logs)
            MD.append(((m * q).sum() * 10 ** min(logs, 300) - mass0) / mass0 if logs < 300 else np.nan)
            LMX.append(np.log10(np.abs(q).max()) + logs)
            if k == nstep or not np.isfinite(LV[-1]):
                break
            h = a.dt
            k1 = L @ q
            k2 = L @ (q + 0.5 * h * k1)
            k3 = L @ (q + 0.5 * h * k2)
            k4 = L @ (q + h * k3)
            q = q + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            for _ in range(s):
                q = q - (h / s) * nu * (B @ q)
            mx = np.abs(q).max()
            if mx > 1e100:
                q /= mx; logs += np.log10(mx)
        print(f'f{form} nu30={n30:.4g}: {nstep} steps of {a.dt:g} s, {s} hypervis subcycles, '
              f'log10 Var/Var0 at end {LV[-1]:.2f}, max mass drift {np.nanmax(np.abs(MD)):.1e} '
              f'({time.time() - t0:.0f} s)', flush=True)
        np.savez(f, t=np.array(T), log10var=np.array(LV), massdrift=np.array(MD), log10max=np.array(LMX),
                 snap_t=np.array(ST), snap_q=np.array(SQ), snap_log10scale=np.array(SS),
                 form=int(form), nu30=n30, nu=nu, dt=a.dt, subcycles=s, gauss=np.array(a.gauss),
                 lat=d['lat'], lon=d['lon'], mass=m, case=a.case, ne=a.ne)
        print(f, flush=True)


if __name__ == '__main__':
    main()
