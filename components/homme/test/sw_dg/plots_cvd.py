#!/usr/bin/env python3
# Figures for the CG vs DG gallery (swtc1, alpha = pi/4, split form, RK4).
# usage: plots_cvd.py RUNS_DIR OUT_DIR
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.colors import SymLogNorm

R, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
REV = 12.0
NES = [8, 15, 30, 60, 120]
SH = [('cosi', 'cosine bell'), ('gaus', 'Gaussian'), ('cyli', 'cylinder'), ('slcy', 'slotted cylinder')]
NUC = os.environ.get('NUC', '1e13')
# key, label, convergence run pattern, 3-revolution run, colour, marker, line style
M = [('cg0', r'CG, $\nu$ = 0', 'conv_s4_ne%d_dg-1', 'cvd_var_dg-1', '#e08a2e', 'o', '-'),
     ('dg0', 'DG, central flux', 'conv_s4_ne%d_dg0', 'cvd_var_dg0', '#8e3b8e', 'D', '--'),
     ('cgn', r'CG, $\nu$ = %s$\,(30/NE)^4$' % NUC, 'conv_s4_ne%d_cgnu', 'cvd_var_cgnu', '#1f5f8b', 'o', '-'),
     ('dg1', 'DG, upwind flux', 'conv_s4_ne%d_dg1', 'cvd_var_dg1', '#2c7a4b', 'D', '--')]
plt.rcParams.update({'font.size': 10})

def load(run, shape, norm):
    f = os.path.join(R, run, 'swtc1.%s.%s.errors' % (shape, norm))
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        return None
    d = np.loadtxt(f, ndmin=2)
    return d[:, 0], d[:, 1]

def at_day(run, shape, norm, day):
    e = load(run, shape, norm)
    if e is None or e[0][-1] < day - 0.2:
        return None
    return e[1][np.argmin(abs(e[0] - day))]

def series(pat, shape, norm):
    ne, v = [], []
    for n in NES:
        x = at_day(pat % n, shape, norm, 12)
        if x is not None:
            ne.append(n); v.append(x)
    return np.array(ne), np.array(v)

def order(ne, e):
    return np.log(e[-2] / e[-1]) / np.log(ne[-1] / ne[-2]) if len(ne) > 1 else np.nan

def var(run):
    f = os.path.join(R, run, 'sweq.var')
    if not os.path.exists(f):
        return None
    d = np.loadtxt(f, ndmin=2)
    _, i = np.unique(d[:, 0], return_index=True)
    d = d[np.sort(i)]
    return d[:, 0] / REV, (d[:, 1:5] - d[0, 1:5]) / d[0, 1:5]

# ---- Gaussian convergence ---------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
for k, lab, pat, _, c, mk, _ls in M:
    ne, l2 = series(pat, 'gaus', 'l2'); _, mx = series(pat, 'gaus', 'linf')
    if len(ne) == 0:
        continue
    ax.loglog(ne, l2, '-', color=c, marker=mk, lw=2.4, ms=8,
              label='%s   l2 %.2f, max %.2f  (NE%d–%d)' % (lab, order(ne, l2), order(ne, mx), ne[-2], ne[-1]))
    ax.loglog(ne, mx, '--', color=c, marker=mk, lw=1.1, ms=4, alpha=.75)
x = np.array([15, 120.])
for p, ls in ((3, ':'), (4, '-.')):
    ax.loglog(x, 1e-2 * (x / 15) ** -p, ls, color='0.45', lw=1)
    ax.text(x[-1] * 0.93, 1e-2 * (x[-1] / 15) ** -p, 'slope %d' % p, color='0.45', ha='right', va='bottom')
ax.set_xticks(NES); ax.set_xticklabels([str(n) for n in NES]); ax.minorticks_off()
ax.set_xlabel('NE  (elements per cube-face edge; dt = 22500 s / NE, RK4)')
ax.set_ylabel('relative error at day 12  (solid: l2, dashed: max)')
ax.grid(alpha=.3, which='both'); ax.set_title('swtc1 Gaussian, α = π/4, split form: CG vs DG')
ax.legend(title='orders between the two finest grids', loc='lower left', fontsize=9)
fig.savefig(os.path.join(OUT, 'conv_gauss.png'), dpi=110); plt.close(fig)

# ---- all shapes, l2 vs NE ---------------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
for a, (s, name) in zip(axs.flat, SH):
    for k, lab, pat, _, c, mk, _ls in M:
        ne, l2 = series(pat, s, 'l2')
        if len(ne):
            a.loglog(ne, l2, '-', color=c, marker=mk, lw=2, ms=6, label='%s (%.2f)' % (lab, order(ne, l2)))
    a.set_xticks(NES); a.set_xticklabels([str(n) for n in NES]); a.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    a.yaxis.set_minor_formatter(matplotlib.ticker.LogFormatterSciNotation(labelOnlyBase=False, minor_thresholds=(2, 0.5)))
    a.set(title=name, xlabel='NE', ylabel='relative l2 error at day 12'); a.grid(alpha=.3, which='both')
    a.legend(fontsize=8, title='order, two finest grids', title_fontsize=8)
fig.suptitle('swtc1, α = π/4, split form: l2 error after one revolution, all four shapes')
fig.savefig(os.path.join(OUT, 'conv_shapes.png'), dpi=110); plt.close(fig)

# ---- l2 error over 3 revolutions --------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(13, 8.6), constrained_layout=True)
for a, (s, name) in zip(axs.flat, SH):
    lo = []
    for k, lab, _, run, c, mk, ls in M:
        e = load(run, s, 'l2')
        if e is not None:
            a.semilogy(e[0] / REV, e[1], ls, color=c, lw=2, label=lab); lo.append(e[1][e[0] > 0.6].min())
    if lo:
        a.set_ylim(bottom=0.6 * min(lo))
    a.set(title=name, xlabel='time (revolutions)', ylabel='relative l2 error', xlim=(0, 3)); a.grid(alpha=.3, which='both')
axs[0, 0].legend()
fig.suptitle('swtc1 NE30, dt = 750 s, α = π/4, split form: l2 error over 3 revolutions')
fig.savefig(os.path.join(OUT, 'errors_3rev.png'), dpi=110); plt.close(fig)

# ---- variance over 3 revolutions --------------------------------------------
V = {k: var(run) for k, _, _, run, *_ in M}
fig, axs = plt.subplots(2, 2, figsize=(13, 8.6), constrained_layout=True)
for j, (a, (s, name)) in enumerate(zip(axs.flat, SH)):
    vals = []
    for k, lab, _, run, c, mk, ls in M:
        v = V[k]
        if v is None:
            continue
        a.plot(v[0], v[1][:, j], ls, color=c, lw=2, label=lab); vals.append(np.abs(v[1][v[0] > 0.5, j]))
    med = [np.median(x) for x in vals if len(x) and np.median(x) > 0]
    lt = 10.0 ** (np.floor(np.log10(min(med))) - 1) if med else 1e-12
    a.set_yscale('symlog', linthresh=lt); a.axhline(0, color='0.4', lw=.8)
    a.set(title=name, xlabel='time (revolutions)', ylabel='(V − V0) / V0', xlim=(0, 3)); a.grid(alpha=.3, which='both')
axs[0, 0].legend()
fig.suptitle('Relative change in tracer variance V = ∫ (g h)², NE30, dt = 750 s, 3 revolutions')
fig.savefig(os.path.join(OUT, 'variance_3rev.png'), dpi=110); plt.close(fig)

# ---- day-36 maps, 4 methods ----------------------------------------------------
def movie(run):
    from scipy.io import netcdf_file
    f = os.path.join(R, run, 'movies', 'swtc11.nc')
    if not os.path.exists(f):
        return None
    nc = netcdf_file(f, 'r', mmap=False)
    return nc.variables['lon'][:], nc.variables['lat'][:], nc.variables['time'][:], nc.variables['geop'][:]

MV = {k: movie(run) for k, _, _, run, *_ in M}
if all(v is not None for v in MV.values()):
    for lev, tag, mode in ((1, 'bell', 'symlog'), (4, 'slotted', 'symlog'), (4, 'slotted', 'linear')):
        hmax = 1.0 if lev == 4 else 1000.0
        fig, axs = plt.subplots(2, 2, figsize=(13, 7.4), constrained_layout=True)
        for a, (k, lab, *_r) in zip(axs.flat, M):
            lon, lat, t, g = MV[k]
            it = np.argmin(abs(t - 36)); fld = g[it, lev - 1]
            if mode == 'symlog':
                lt = 1e-4 * hmax
                pc = a.pcolormesh(lon, lat, fld, cmap='RdBu_r', shading='auto',
                                  norm=SymLogNorm(linthresh=lt, vmin=-hmax, vmax=hmax, base=10))
            else:
                pc = a.pcolormesh(lon, lat, fld, cmap='Blues', vmin=0, vmax=hmax, shading='auto')
            a.set_aspect('equal'); a.set_xlim(0, 360); a.set_ylim(-90, 90)
            a.set_title('%s   day %.0f   min %.3g  max %.3g' % (lab, t[it], fld.min(), fld.max()), fontsize=10)
        fig.colorbar(pc, ax=axs, shrink=0.8, extend='both')
        fig.suptitle('%s at day 36 (3 revolutions), %s color scale' % ('cosine bell [m]' if lev == 1 else 'slotted cylinder', mode))
        fig.savefig(os.path.join(OUT, 'map36_%s_%s.png' % (tag, mode)), dpi=100); plt.close(fig)

# ---- tables -------------------------------------------------------------------
with open(os.path.join(OUT, 'table.txt'), 'w') as f:
    f.write('day 12, NE30 l2/max; day 36 (NE30 dt750) l2/max; variance change at end\n')
    for s, name in SH:
        j = [x[0] for x in SH].index(s)
        row = [name]
        for k, lab, pat, run, *_ in M:
            a = at_day(pat % 30, s, 'l2', 12); b = at_day(pat % 30, s, 'linf', 12)
            c = at_day(run, s, 'l2', 36); d = at_day(run, s, 'linf', 36)
            v = V[k][1][-1, j] if V[k] is not None else np.nan
            row.append('%s: %s/%s | %s/%s | dV %.3g' % (k, *['%.3g' % x if x is not None else '-' for x in (a, b, c, d)], v))
        f.write('  '.join(row) + '\n')
    for k, lab, pat, *_ in M:
        ne, l2 = series(pat, 'gaus', 'l2'); _, mx = series(pat, 'gaus', 'linf')
        f.write('%s gaus NE %s l2 %s max %s\n' % (k, list(ne), ['%.3g' % x for x in l2], ['%.3g' % x for x in mx]))
print(open(os.path.join(OUT, 'table.txt')).read())
