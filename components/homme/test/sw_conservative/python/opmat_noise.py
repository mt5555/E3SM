"""
Random discretely divergence-free velocities for the swtc1 transport operator L:
u = DSS(curl psi) with psi drawn at the unique GLL nodes (nothing is sampled from a
continuous field, so there is no aliasing), optionally smoothed with HOMME's own
hyperviscosity operator B and windowed to a region, then passed to sweqx with
SWEQX_TC1_PSI=file (see shallow_water_mod:tc1_init_velocity).  With r given, the
velocity is u_sb + c*u_noise with rms|c*u_noise| = r*rms|u_sb|; without r it is the
noise alone.  Every field is rescaled to max|u| = u0 at the GLL nodes.

    python opmat_noise.py gen     --ref PREFIX --psi DIR --lc 4 8 16 32 inf --seeds 10
                                  [--window LAT LON R]
    python opmat_noise.py patch   --ref PREFIX --psi DIR --mats DIR --src noise-ne30-lcinf-s2
                                  --window LAT LON FLAT TAPER [--bg 0 1]
    python opmat_noise.py extract --psi DIR --mats DIR [--select 'noise*'] [--r none 0.3 1]
                                  [--forms 0 2] [--ne 8]
    python opmat_noise.py export  --prefix MATS/ne30-sc1-f0-curlpsi-<case> --out FILE.nc [--attr k=v ...]
    python opmat_noise.py analyze --mats DIR --out OUT [--select 'noise*'] [--spectral 6]
                                  [--maps 4] [--tag noise]

gen: white noise w ~ N(0,1) at the nodes (seeded by --seed0+k only, so a seed gives
the same noise at every lc), psi = (I + (a^4/(lc(lc+1))^2) B)^-1 w (B ~ (l(l+1)/a^2)^2 on
Y_l^m, so a smooth low-pass at degree lc; lc=inf: no smoothing, grid-scale noise),
times exp(-(d/R)^2) for --window (d = great-circle distance, degrees), mass-weighted
mean removed, scaled to max|psi| = 1.  Files: noise-ne<NE>-lc<lc>-s<seed>[-w<lat>_<lon>_<R>].bin
(stream records int64 gdofP, real64 psi), gdofP from the --ref .dofs files.

analyze: for every ne<NE>-sc<SC>-f0-curlpsi-noise* case: the opmat_flows.analyze
diagnostics (dense eigenvalues, skew, top growing modes and where they sit), the shear
scale zeta_max = max|DSS vorticity|, umax/urms, the split-form skew when the f2 partner
exists, and the spectral hyperviscosity threshold nu30 for the --spectral worst fields.
"""
import argparse
import glob
import os
import re
import subprocess
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from opmat_io import load, read_meta, psi_file
from opmat_eigs import Cache, skew_sym_err, MUTED
from opmat_sym import xyz
from opmat_flows import analyze, summarize, regions, style, cube_edge_lines, GROW
from opmat_flows_hv import load_B, spectral_nu, nu30, nu_from30, TOL

HERE = os.path.dirname(os.path.abspath(__file__))
NAME = re.compile(r'ne(\d+)-sc(\d+)-f(\d)-curlpsi-noise-ne\d+-lc(\w+?)-s(\d+)(?:-w(\S+?))?(?:_r(\S+))?$')


def lc_str(lc):
    return 'inf' if np.isinf(lc) else f'{lc:g}'


# ---------------------------------------------------------------- gen
def gen(args):
    d = load(args.ref, ops=('B',))
    if 'B' not in d:
        raise SystemExit(f'{args.ref}: no B')
    n, a, B, mass = d['n'], d['meta']['rearth'], d['B'], d['mass']
    ne = d['meta']['ne']
    os.makedirs(args.psi, exist_ok=True)
    X = xyz(d['lat'], d['lon'])
    win, wtag = None, ''
    if args.window:
        la, lo, R = args.window
        c = xyz(np.radians([la]), np.radians([lo]))[0]
        dist = np.degrees(np.arccos(np.clip(X @ c, -1, 1)))
        win = np.exp(-(dist / R) ** 2)
        wtag = f'-w{la:g}_{lo:g}_{R:g}'
    I = sp.identity(n, format='csc')
    for lc in [float(x) for x in args.lc]:
        if not np.isinf(lc):
            lu = spla.splu((I + (a ** 4 / (lc * (lc + 1)) ** 2) * B).tocsc())
        for s in range(args.seed0, args.seed0 + args.seeds):
            w = np.random.default_rng(s).standard_normal(n)
            psi = w if np.isinf(lc) else lu.solve(w)
            if win is not None:
                psi = psi * win
            psi -= (mass * psi).sum() / mass.sum()
            psi /= np.abs(psi).max()
            rec = np.zeros(n, dtype=[('g', '<i8'), ('p', '<f8')])
            rec['g'], rec['p'] = d['gdof'], psi
            f = os.path.join(args.psi, f'noise-ne{ne}-lc{lc_str(lc)}-s{s}{wtag}.bin')
            rec.tofile(f)
            print(f)


# ---------------------------------------------------------------- patch
def patch(args):
    """psi = bg*sin(lat) + uscale*w*psi_src: the noise of an existing file, at its physical
    amplitude (uscale = psi_uscale of the extracted pure-noise case, so that inside the patch the
    velocity is the one that case was analyzed with), in a flat-top window w = 1 within FLAT deg
    of the centre and exp(-((d-FLAT)/TAPER)^2) outside, on a solid-body background bg*u0
    (bg=0: patch alone).  sweqx still rescales the total to max|u| = u0."""
    d = load(args.ref, ops=())
    src = np.fromfile(os.path.join(args.psi, args.src + '.bin'), dtype=[('g', '<i8'), ('p', '<f8')])
    g2i = {g: i for i, g in enumerate(src['g'])}
    psi_n = src['p'][[g2i[g] for g in d['gdof']]]
    # with no background the amplitude is irrelevant (sweqx rescales to max|u| = u0)
    uscale = 1.0 if not any(args.bg) else \
        read_meta(os.path.join(args.mats, f'ne{d["meta"]["ne"]}-sc1-f0-curlpsi-{args.src}'))['psi_uscale']
    X = xyz(d['lat'], d['lon'])
    la, lo, flat, taper = args.window
    c = xyz(np.radians([la]), np.radians([lo]))[0]
    dist = np.degrees(np.arccos(np.clip(X @ c, -1, 1)))
    w = np.exp(-(np.maximum(dist - flat, 0) / taper) ** 2)
    for bg in args.bg:
        psi = bg * np.sin(d['lat']) + uscale * w * psi_n
        rec = np.zeros(len(psi), dtype=[('g', '<i8'), ('p', '<f8')])
        rec['g'], rec['p'] = d['gdof'], psi
        f = os.path.join(args.psi, f'{args.src}-w{la:g}_{lo:g}_{flat:g}_{taper:g}_bg{bg:g}.bin')
        rec.tofile(f)
        print(f)


# ---------------------------------------------------------------- extract
def extract(args):
    script = os.path.join(HERE, '..', 'swtc1-opmat.sh')
    env0 = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1',
                HOMME=os.environ.get('HOMME', os.path.abspath(os.path.join(HERE, '..', '..', '..'))),
                outdir=args.mats, SWEQX_TC1_PSI='file')
    files = sorted(set(f for sel in args.select for f in glob.glob(os.path.join(args.psi, sel + '.bin'))))
    for f in files:
        for r in args.r:
            for form in args.forms:
                tag = os.path.basename(f)[:-4] + ('' if r == 'none' else f'_r{r}')
                pre = os.path.join(args.mats, f'ne{args.ne}-sc{args.sc}-f{form}-curlpsi-{tag}')
                if os.path.exists(pre + '.meta'):
                    continue
                env = dict(env0, SWEQX_TC1_PSI_FILE=f)
                env.pop('SWEQX_TC1_PSI_PARAMS', None)
                if r != 'none':
                    env['SWEQX_TC1_PSI_PARAMS'] = r
                p = subprocess.run([script, str(args.ne), str(args.sc), str(form), 'curlpsi', 'L',
                                    str(args.nranks)], env=env, capture_output=True, text=True)
                if p.returncode != 0 or not os.path.exists(pre + '.meta'):
                    raise SystemExit(f'{pre}: extraction failed\n{p.stdout}\n{p.stderr}')
                print(pre, flush=True)


# ---------------------------------------------------------------- analyze
def describe(row):
    w = f' w{row["win"]}' if row['win'] else ''
    r = 'noise' if np.isinf(row['r']) else f'r={row["r"]:g}'
    return f'{r} lc={row["lc"]} s{row["seed"]}{w}'


def plot_scan(rows, out, fname):
    """max Re lambda / omega vs lc for each (window, r): seeds as dots, median as a line;
    and max Re lambda against the vorticity scale zeta_max"""
    keys = sorted(set((row['win'], row['r']) for row in rows), key=lambda k: (k[0], k[1]))
    lcs = sorted(set(row['lc'] for row in rows), key=lambda s: float(s))
    xpos = {lc: i for i, lc in enumerate(lcs)}
    pal = ['#1baf7a', '#2a78d6', '#eb6834', '#6250d6', '#c2409a', '#898781', '#0b0b0b', '#d6a320']
    fig, axs = plt.subplots(1, 3, figsize=(14, 3.8))
    for ik, k in enumerate(keys):
        c = pal[ik % len(pal)]
        lab = ('noise' if np.isinf(k[1]) else f'r={k[1]:g}') + (f', w{k[0]}' if k[0] else '')
        rs = [row for row in rows if (row['win'], row['r']) == k]
        med, mx = [], []
        for lc in lcs:
            v = [row['maxre'] for row in rs if row['lc'] == lc]
            if not v:
                med.append(np.nan); mx.append(np.nan)
                continue
            jit = (np.arange(len(v)) - len(v) / 2) * 0.02 + ik * 0.06
            axs[0].plot(xpos[lc] + jit, np.maximum(v, 1e-16), '.', color=c, ms=4, alpha=0.6)
            med.append(np.median(v)); mx.append(np.max(v))
        axs[0].plot(np.arange(len(lcs)) + ik * 0.06, np.maximum(med, 1e-16), '-', color=c, lw=1.5, label=lab)
        axs[1].plot([row['zmax'] for row in rs], [row['maxre'] for row in rs], 'o', color=c, ms=3, label=lab)
        axs[2].plot([xpos[row['lc']] + ik * 0.06 for row in rs], [row['maxre'] / row['zmax'] for row in rs],
                    'o', color=c, ms=3)
    axs[0].set_yscale('log')
    axs[0].set_ylabel('max Re λ / ω   (line: median over seeds)')
    axs[1].set_xscale('log'); axs[1].set_yscale('log')
    axs[1].set_xlabel('max |ζ| / ω'); axs[1].set_ylabel('max Re λ / ω')
    axs[2].set_ylabel('max Re λ / max |ζ|')
    for ax in (axs[0], axs[2]):
        ax.set_xticks(range(len(lcs)))
        ax.set_xticklabels(lcs)
        ax.set_xlabel('smoothing degree lc (inf = grid-scale)')
    for ax in axs:
        style(ax)
    axs[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, fname))
    plt.close(fig)


def plot_modes(rows, out, fname, nmodes=2):
    """per case: speed |u|/u_max, then |v|^2 of the top growing modes"""
    edges = cube_edge_lines()
    fig, axs = plt.subplots(len(rows), nmodes + 1, figsize=(4.2 * (nmodes + 1), 2.3 * len(rows)), squeeze=False)
    for i, row in enumerate(rows):
        r = row['res']
        om = float(r['om'])
        lon, lat = np.degrees(r['lon']) % 360, np.degrees(r['lat'])
        for k in range(nmodes + 1):
            ax = axs[i][k]
            for lo, la in edges:
                ax.plot(lo, la, color='#2a78d6', lw=0.5)
            ax.set_xlim(0, 360); ax.set_ylim(-90, 90)
            ax.set_xticks([]); ax.set_yticks([])
            if k == 0:
                sp_ = r['speed'] / r['speed'].max()
                ax.scatter(lon, lat, c=sp_, s=3, cmap='viridis', vmin=0, vmax=1, linewidths=0)
                ax.set_title(f'{describe(row)}: |u|/u_max', fontsize=7)
                continue
            if k - 1 >= r['V'].shape[1]:
                continue
            e = np.abs(r['V'][:, k - 1]) ** 2
            e /= e.max()
            o = np.argsort(e)
            ax.scatter(lon[o], lat[o], c=e[o], s=2 + 14 * e[o], cmap='magma_r', vmin=0, vmax=1, linewidths=0)
            lam = r['vlam'][k - 1] / om
            ax.set_title(f'mode {k}: λ/ω = {lam.real:.3f}{lam.imag:+.2f}i', fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, fname))
    plt.close(fig)


def analyze_all(args):
    os.makedirs(args.out, exist_ok=True)
    cache = Cache(args.out)
    pat = [os.path.join(args.mats, f'ne{args.ne}-sc{args.sc}-f0-curlpsi-{sel}.meta') for sel in args.select]
    files = sorted(set(f[:-5] for p in pat for f in glob.glob(p)))
    rows = []
    for path in files:
        nm = os.path.basename(path)
        m = NAME.match(nm)
        if not m:
            continue
        res = analyze(cache, path, args.modes)
        s = summarize(res, args.ne)
        meta = res['meta'][()]
        om = float(res['om'])
        urms = np.sqrt((res['mass'] * res['speed'] ** 2).sum() / res['mass'].sum())
        p2 = path.replace(f'-f0-', '-f2-')
        sk2 = np.nan
        if os.path.exists(p2 + '.meta'):
            sk2 = float(cache.get('skew-' + os.path.basename(p2),
                                  lambda: dict(skew=skew_sym_err(load(p2, ops=('L',))['L'], res['mass'])[0]))['skew'])
        row = dict(name=nm, path=path, lc=m.group(4), seed=int(m.group(5)), win=m.group(6) or '',
                   r=meta['psi_r'] if meta['psi_r'] >= 0 else np.inf, maxre=s['maxre'], ngrow=s['ngrow'],
                   zmax=float(res['zmax']) / om, umax_rms=res['speed'].max() / urms, skew=s['skew'], skew2=sk2,
                   slow5=s['slow5'], enr=list(s['enr'][0].values()) if s['enr'] else [np.nan] * 4, maxabs=s['maxabs'], om=om,
                   res=res)
        rows.append(row)
        print(f'{nm}: maxRe={row["maxre"]:.3e} ngrow={row["ngrow"]} zmax={row["zmax"]:.1f} '
              f'skew f0={row["skew"]:.3f} f2={sk2:.1e}', flush=True)

    # spectral hyperviscosity threshold for the worst fields
    worst = sorted(rows, key=lambda row: -row['maxre'])
    for row in worst[:args.spectral]:
        def spec():
            d = load(row['path'], ops=('L',))
            Bl = load_B(args.mats, args.ne, d['lat'], d['lon'])
            ns, g = spectral_nu(d['L'], Bl, row['om'], args.ne, TOL)
            return dict(nus=ns, g=g)
        sres = cache.get(f'hv-spec-tol{TOL:g}-' + row['name'], spec)
        row['nus'] = float(sres['nus'])
        print(f'{row["name"]}: nu30 spectral={row["nus"]:.2e}', flush=True)

    out = os.path.join(args.out, f'results_noise-{args.tag}.txt')
    with open(out, 'w') as f:
        f.write(f'NE{args.ne} sc{args.sc}: u = DSS(curl psi), psi random at the GLL nodes (see opmat_noise.py); '
                f'omega = u0/a, max|u| = u0; conservation form (f0)\n')
        f.write('  zeta = max|DSS vorticity|/omega; skew = ||S+S^T||/||S|| (f0, f2); nu30 = spectral threshold '
                '[m^4/s]; enrichment [corner cubeedge elemedge slow] of the top mode\n')
        f.write(f'{"case":40s} {"r":>5s} {"lc":>4s} {"seed":>4s} {"maxRe/om":>9s} {"ngrow":>5s} {"zeta":>7s} '
                f'{"maxRe/z":>7s} {"umax/rms":>8s} {"skew f0":>7s} {"skew f2":>7s} {"slow5":>5s} {"nu30":>8s}  enrichment\n')
        for row in sorted(rows, key=lambda row: (row['win'], row['r'], float(row['lc']), row['seed'])):
            tag = row['name'].split('-curlpsi-')[1]
            f.write(f'{tag:40s} {row["r"]:5g} {row["lc"]:>4s} {row["seed"]:4d} {row["maxre"]:9.3e} {row["ngrow"]:5d} '
                    f'{row["zmax"]:7.1f} {row["maxre"] / row["zmax"]:7.4f} {row["umax_rms"]:8.2f} {row["skew"]:7.3f} '
                    f'{row["skew2"]:7.1e} {row["slow5"]:5.3f} {row.get("nus", np.nan):8.2e}  '
                    f'[{" ".join(f"{e:.1f}" for e in row["enr"])}]\n')
        f.write('\nsummary over seeds: median / max of max Re lambda / omega\n')
        for k in sorted(set((row['win'], row['r'], row['lc']) for row in rows), key=lambda k: (k[0], k[1], float(k[2]))):
            v = [row['maxre'] for row in rows if (row['win'], row['r'], row['lc']) == k]
            f.write(f'  w={k[0] or "-":12s} r={k[1]:<5g} lc={k[2]:>4s}: n={len(v):2d}  median {np.median(v):.3e}  '
                    f'max {np.max(v):.3e}  growing in {sum(x > GROW for x in v)}/{len(v)}\n')
    print(open(out).read())
    plot_scan(rows, args.out, f'ne{args.ne}-sc{args.sc}-{args.tag}-scan.png')
    if args.maps:
        plot_modes(worst[:args.maps], args.out, f'ne{args.ne}-sc{args.sc}-{args.tag}-modes.png')


# ---------------------------------------------------------------- sparse (large NE)
def sym_bound(A, mass):
    """largest eigenvalue of sym(M^1/2 A M^-1/2): an upper bound for max Re lambda(A)"""
    sq = np.sqrt(mass)
    S = sp.diags(sq) @ A @ sp.diags(1 / sq)
    return spla.eigsh((0.5 * (S + S.T)).tocsc(), k=1, which='LA', return_eigenvectors=False)[0]


def rightmost(A, om, bound, k=16, floor=1e-3, fracs=(0.5, 0.2, 0.05)):
    """max Re lambda(A) by shift-invert Arnoldi at real shifts from the symmetric-part bound
    down to floor*omega (the growing modes found so far have Im lambda ~ 0); also the
    eigenvalues found"""
    shifts = [bound * 1.02 + floor * om] + [bound * f for f in fracs] + [floor * om]
    lam = []
    A = A.tocsc()
    for sg in sorted(set(x for x in shifts if x >= floor * om), reverse=True):
        try:
            lam.append(spla.eigs(A, k=k, sigma=sg, which='LM', return_eigenvectors=False, tol=1e-10))
        except spla.ArpackNoConvergence as e:
            lam.append(e.eigenvalues)
    lam = np.concatenate(lam)
    return lam.real.max(), lam


def nu_sparse(L, B, mass, om, ne, tol, bound, grid=(1e12, 3e12, 1e13, 3e13, 1e14, 3e14, 1e15, 3e15, 1e16),
              nbis=7, fracs=(0.5, 0.2, 0.05)):
    """spectral nu30 threshold, max Re lambda(L - nu B) <= tol*omega, with rightmost():
    coarse grid, then bisection in log nu30 above the last unstable grid point.  bound =
    lambda_max(sym S) of L, also a bound for L - nu B (M^1/2 B M^-1/2 is positive semidefinite)"""
    def unstable(n30):
        return rightmost(L - nu_from30(n30, ne) * B, om, bound, fracs=fracs)[0] > tol * om
    lo = None
    for n30 in grid:
        if not unstable(n30):
            hi = n30
            break
        lo = n30
    else:
        return np.inf
    if lo is None:
        return grid[0]
    for _ in range(nbis):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if unstable(mid) else (lo, mid)
    return hi


def eigs_near(A, sg, k=8, maxiter=40, tol=1e-8):
    """converged eigenvalues of A nearest the real shift sg, by shift-invert with a fill-reducing
    LU (minimum degree on A^T+A: about 3x less fill than the default COLAMD for L - nu B, whose B
    stencil is wide) and capped restarts.  Far from any eigenvalue ARPACK converges slowly, and
    with the default restart limit (10 n) a single call can run for hours at NE30; with the cap,
    the eigenvalues that did converge are returned (possibly none)"""
    n = A.shape[0]
    lu = spla.splu((A - sg * sp.identity(n, format='csc')).tocsc(), permc_spec='MMD_AT_PLUS_A',
                   diag_pivot_thresh=0.1)
    op = spla.LinearOperator((n, n), matvec=lu.solve, dtype=float)
    try:
        return spla.eigs(A, k=k, sigma=sg, OPinv=op, which='LM', return_eigenvectors=False, tol=tol,
                         maxiter=maxiter)
    except spla.ArpackNoConvergence as e:
        return e.eigenvalues


def nu_track(L, B, om, ne, tol, lam0, grid=(3e14, 1e15, 2e15, 4e15, 8e15, 1.6e16), nbis=5,
             fracs=(1.05, 0.4, 0.15), low=0.03, log=lambda s: print(s, flush=True)):
    """spectral nu30 threshold for large NE (NE30), where nu_sparse() stalls.  A test of L - nu B
    places real shifts at fracs x ref, ref = the rightmost Re lambda found at the largest unstable
    nu so far (lam0 = max Re lambda(L) to start; the growing modes have Im lambda ~ 0 and Re lambda
    falls as nu grows), then one at low*omega with k=16.  Unstable as soon as a converged eigenvalue
    has Re > tol*omega.  Coarse grid, then bisection in log nu30 (resolution (2^(1/2^nbis)) for a
    factor-2 bracket)"""
    ref = [lam0]

    def unstable(n30):
        t0 = time.time()
        A = (L - nu_from30(n30, ne) * B).tocsc()
        shifts = [(f * ref[0], 8) for f in fracs if f * ref[0] > low * om] + [(low * om, 16)]
        for sg, k in shifts:
            lam = eigs_near(A, sg, k=k)
            mr = lam.real.max() if len(lam) else -np.inf
            if mr > tol * om:
                ref[0] = mr
                log(f'    nu30={n30:.3e}: unstable, Re/om={mr / om:.4f} (shift {sg / om:.3f}, {time.time() - t0:.0f} s)')
                return True
        log(f'    nu30={n30:.3e}: stable ({len(shifts)} shifts, {time.time() - t0:.0f} s)')
        return False
    lo = None
    for n30 in grid:
        if not unstable(n30):
            hi = n30
            break
        lo = n30
    else:
        return np.inf
    if lo is None:
        return grid[0]
    for _ in range(nbis):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if unstable(mid) else (lo, mid)
    return hi


def sparse_all(args):
    """max Re lambda, symmetric-part bound, vorticity scale and the spectral nu30 threshold
    without dense matrices (large NE); --dense also computes the dense max Re for checking"""
    os.makedirs(args.out, exist_ok=True)
    cache = Cache(args.out)
    rows = []
    for ne in args.ne:
        pat = [os.path.join(args.mats, f'ne{ne}-sc{args.sc}-f0-curlpsi-{sel}.meta') for sel in args.select]
        for path in sorted(set(f[:-5] for p in pat for f in glob.glob(p))):
            nm = os.path.basename(path)
            m = NAME.match(nm)
            if not m:
                continue

            def comp():
                t0 = time.time()
                d = load(path, ops=('L',))
                L, mass, meta = d['L'], d['mass'], d['meta']
                om = meta['u0'] / meta['rearth']
                bnd = sym_bound(L, mass)
                mr, lam = rightmost(L, om, bnd)
                print(f'  {nm}: bound {bnd / om:.2f}, max Re {mr / om:.3f} ({time.time() - t0:.0f} s)', flush=True)
                dmr = np.linalg.eigvals(L.toarray()).real.max() if args.dense else np.nan
                speed = np.linalg.norm(d['vel'], axis=1)
                out = dict(maxre=mr / om, bound=bnd / om, dense=dmr / om, zmax=np.abs(d['zeta']).max() / om,
                           umax_rms=speed.max() / np.sqrt((mass * speed ** 2).sum() / mass.sum()),
                           sk=skew_sym_err(L, mass)[0], nus=np.nan)
                if args.nu:
                    Bl = load_B(args.mats, ne, d['lat'], d['lon'])
                    if args.track:
                        out['nus'] = nu_track(L, Bl, om, ne, TOL, mr, **({'grid': args.grid} if args.grid else {}))
                    else:
                        kw = dict(grid=(1e14, 3e14, 1e15, 3e15, 1e16), nbis=5, fracs=(0.3,)) if args.fast else {}
                        out['nus'] = nu_sparse(L, Bl, mass, om, ne, TOL, bnd, **kw)
                out['time'] = time.time() - t0
                return out
            mode = '-track' if args.track else '-fast' if args.fast else ''
            if args.track and args.grid:
                mode += f'-g{args.grid[0]:.0e}'
            key = f'sparse{"-nu" if args.nu else ""}{mode}{"-dense" if args.dense else ""}-' + nm
            r = cache.get(key, comp)
            row = dict(name=nm, ne=ne, lc=m.group(4), seed=int(m.group(5)), win=m.group(6) or '',
                       r=float(m.group(7)) if m.group(7) else np.inf, **{k: float(v) for k, v in r.items()})
            rows.append(row)
            print(f'{nm}: maxRe={row["maxre"]:.3f} (dense {row["dense"]:.3f}) bound={row["bound"]:.2f} '
                  f'zeta={row["zmax"]:.1f} nu30={row["nus"]:.2e}  {row["time"]:.0f} s', flush=True)
    out = os.path.join(args.out, f'results_noise-sparse-{args.tag}.txt')
    with open(out, 'w') as f:
        f.write('random nodal psi, conservation form: max Re lambda by shift-invert (dense check where given); '
                'bound = lambda_max(sym S); zeta = max|vorticity|/omega; nu30 = nu*(NE/30)^4 spectral threshold; '
                'nu_actual = nu30*(30/NE)^4; K = nu_actual/(max|zeta| dx^4), dx = 2 pi a/(12 NE) the mean GLL spacing\n')
        f.write(f'{"case":44s} {"NE":>3s} {"maxRe/om":>9s} {"dense":>9s} {"bound":>7s} {"zeta":>7s} {"maxRe/z":>7s} '
                f'{"umax/rms":>8s} {"skew":>6s} {"nu30":>9s} {"nu_actual":>9s} {"K":>7s}\n')
        for row in sorted(rows, key=lambda row: (row['lc'], row['r'], row['seed'], row['ne'])):
            a, om = 6.37122e6, 2 * np.pi / (12 * 86400.0)
            dx = 2 * np.pi * a / (12 * row['ne'])
            nua = row['nus'] * (30 / row['ne']) ** 4
            K = nua / (row['zmax'] * om * dx ** 4)
            f.write(f'{row["name"].split("-curlpsi-")[1]:44s} {row["ne"]:3d} {row["maxre"]:9.3f} {row["dense"]:9.3f} '
                    f'{row["bound"]:7.2f} {row["zmax"]:7.1f} {row["maxre"] / row["zmax"]:7.4f} {row["umax_rms"]:8.2f} '
                    f'{row["sk"]:6.3f} {row["nus"]:9.2e} {nua:9.2e} {K:7.4f}\n')
    print(open(out).read())


# ---------------------------------------------------------------- export
def export(args):
    """the velocity field of an extracted curlpsi case at the GLL nodes, as netCDF-3 for other codes"""
    from scipy.io import netcdf_file
    d = load(args.prefix, ops=())
    meta = d['meta']
    a, u0 = meta['rearth'], meta['u0']
    rec = np.fromfile(psi_file(args.prefix, meta), dtype=[('g', '<i8'), ('p', '<f8')])
    g2i = {g: i for i, g in enumerate(rec['g'])}
    psi_f = rec['p'][[g2i[g] for g in d['gdof']]]
    f = netcdf_file(args.out, 'w', version=2)
    f.createDimension('node', d['n'])
    def var(name, data, typ, units, long):
        v = f.createVariable(name, typ, ('node',))
        v[:] = data
        v.units, v.long_name = units, long
    var('lat', np.degrees(d['lat']), 'd', 'degrees_north', 'latitude of the unique GLL node')
    var('lon', np.degrees(d['lon']), 'd', 'degrees_east', 'longitude of the unique GLL node')
    var('gdofP', d['gdof'].astype(np.int32), 'i', '1', 'HOMME elem%gdofP id of the node (the key of the psi file)')
    var('mass', d['mass'] * a ** 2, 'd', 'm2', 'assembled GLL quadrature weight (1/rspheremp times a^2)')
    var('psi_file', psi_f, 'd', '1', 'nodal values read by sweqx from SWEQX_TC1_PSI_FILE (psi/(a u0) up to the '
        'rescaling)')
    var('psi', psi_f * meta['psi_uscale'] * a * u0, 'd', 'm2 s-1',
        'stream function: psi_uscale*a*u0*psi_file, so that (u,v) = DSS(curl psi)')
    var('u', d['vel'][:, 0], 'd', 'm s-1', 'eastward velocity, DSS(curl psi) at the GLL node')
    var('v', d['vel'][:, 1], 'd', 'm s-1', 'northward velocity, DSS(curl psi) at the GLL node')
    var('zeta', d['zeta'], 'd', 's-1', 'DSS relative vorticity of (u,v)')
    f.title = os.path.basename(args.prefix)
    f.grid = (f'HOMME cubed sphere, equiangular gnomonic, ne={meta["ne"]}, np={meta["np"]} (GLL), '
              f'{d["n"]} unique nodes; nodes are identified by lat/lon (or gdofP within HOMME)')
    for k, v in (('rearth', a), ('u0', u0), ('omega', u0 / a), ('psi_uscale', meta['psi_uscale'])):
        setattr(f, k, np.array([v], dtype='>f8'))   # a bare float is written as float32
    f.velocity = ('u = DSS(curl psi) with the DSS (assembled, mass-weighted average) at shared nodes, '
                  'then scaled by psi_uscale to max|u| = u0 at the nodes; discretely divergence-free')
    f.source = f'sweqx SWEQX_TC1_VELOCITY=curlpsi SWEQX_TC1_PSI=file SWEQX_TC1_PSI_FILE={os.path.basename(meta["psi_file"])}'
    for kv in args.attr:
        k, v = kv.split('=', 1)
        setattr(f, k, v)
    f.close()
    print(args.out)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    g = sub.add_parser('gen')
    g.add_argument('--ref', required=True, help='opmat prefix with B and the dofs (e.g. ne8-sc1-f2-analytic)')
    g.add_argument('--psi', required=True)
    g.add_argument('--lc', nargs='+', default=['inf'])
    g.add_argument('--seeds', type=int, default=10)
    g.add_argument('--seed0', type=int, default=0)
    g.add_argument('--window', type=float, nargs=3, metavar=('LAT', 'LON', 'R'))
    pt = sub.add_parser('patch')
    pt.add_argument('--ref', required=True, help='opmat prefix with the dofs (e.g. ne30-sc1-f2-analytic)')
    pt.add_argument('--psi', required=True)
    pt.add_argument('--mats', help='where the pure-noise case of --src was extracted (psi_uscale; '
                    'needed only with a nonzero --bg)')
    pt.add_argument('--src', required=True, help='noise file name without .bin, e.g. noise-ne30-lcinf-s2')
    pt.add_argument('--window', type=float, nargs=4, required=True, metavar=('LAT', 'LON', 'FLAT', 'TAPER'))
    pt.add_argument('--bg', type=float, nargs='+', default=[0.0], help='solid-body background amplitudes (x u0)')
    e = sub.add_parser('extract')
    e.add_argument('--psi', required=True)
    e.add_argument('--mats', required=True)
    e.add_argument('--select', nargs='+', default=['noise*'])
    e.add_argument('--r', nargs='+', default=['none'], help="'none' (noise alone) or r values")
    e.add_argument('--forms', nargs='+', default=['0', '2'])
    e.add_argument('--ne', type=int, default=8)
    e.add_argument('--sc', type=int, default=1)
    e.add_argument('--nranks', type=int, default=4)
    a = sub.add_parser('analyze')
    a.add_argument('--mats', required=True)
    a.add_argument('--out', required=True)
    a.add_argument('--select', nargs='+', default=['noise*'])
    a.add_argument('--ne', type=int, default=8)
    a.add_argument('--sc', type=int, default=1)
    a.add_argument('--modes', type=int, default=2)
    a.add_argument('--spectral', type=int, default=6)
    a.add_argument('--maps', type=int, default=4)
    a.add_argument('--tag', default='noise')
    q = sub.add_parser('sparse')
    q.add_argument('--mats', required=True)
    q.add_argument('--out', required=True)
    q.add_argument('--select', nargs='+', default=['noise*'])
    q.add_argument('--ne', type=int, nargs='+', default=[8])
    q.add_argument('--sc', type=int, default=1)
    q.add_argument('--nu', action='store_true', help='also the spectral nu30 threshold')
    q.add_argument('--dense', action='store_true', help='also the dense max Re lambda (small NE only)')
    q.add_argument('--fast', action='store_true', help='nu search: 3 shifts, grid from nu30=1e14, 5 bisections '
                   '(resolution about 1.04x); assumes nu30 > 1e14')
    q.add_argument('--track', action='store_true', help='nu search with nu_track(): shifts placed from the '
                   'rightmost eigenvalue at the last unstable nu, capped ARPACK restarts (needed at NE30); grid '
                   'from nu30=3e14')
    q.add_argument('--grid', type=float, nargs='+', help='--track: coarse nu30 grid (default 3e14 1e15 2e15 '
                   '4e15 8e15 1.6e16); the result is grid[0] when that is already stable')
    q.add_argument('--tag', default='sparse')
    x = sub.add_parser('export')
    x.add_argument('--prefix', required=True, help='opmat prefix of an extracted curlpsi case')
    x.add_argument('--out', required=True, help='netCDF file')
    x.add_argument('--attr', nargs='*', default=[], help='extra global attributes, key=value')
    args = p.parse_args()
    {'gen': gen, 'patch': patch, 'extract': extract, 'analyze': analyze_all, 'sparse': sparse_all, 'export': export}[args.cmd](args)


if __name__ == '__main__':
    main()
