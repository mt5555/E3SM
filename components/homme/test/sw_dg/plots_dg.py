#!/usr/bin/env python3
# Figures for the swtc1-DG Run Gallery (styles follow the swtc1 Run Gallery).
# usage: plots_dg.py RUNS_DIR OUT_DIR
import os, sys, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
REV = 12.0
SH = [('cosi', 'cosine bell'), ('gaus', 'gaussian'), ('cyli', 'cylinder'), ('slcy', 'slotted cylinder')]
RUNS = {('central', 'split'): 'gal_dg0_f2', ('upwind', 'split'): 'gal_dg1_f2',
        ('central', 'cons'): 'gal_dg0_f0', ('upwind', 'cons'): 'gal_dg1_f0'}
FORM = {'split': 'split form', 'cons': 'conservation form'}
C_SPLIT, C_CONS = '#1f5f8b', '#e08a2e'
plt.rcParams.update({'font.size': 10})

def err(run, shape, norm):
    f = os.path.join(R, run, 'swtc1.%s.%s.errors' % (shape, norm))
    if not os.path.exists(f):
        return None
    d = np.loadtxt(f, ndmin=2)
    if d.size == 0 or d.shape[1] < 2:
        return None
    return d[:, 0] / REV, d[:, 1]

def var(run):
    f = os.path.join(R, run, 'sweq.var')
    if not os.path.exists(f):
        return None
    d = np.loadtxt(f, ndmin=2)
    if d.size == 0:
        return None
    _, i = np.unique(d[:, 0], return_index=True)     # first line is written twice
    d = d[np.sort(i)]
    return d[:, 0] / REV, (d[:, 1:5] - d[0, 1:5]) / d[0, 1:5]

# ---- cosine bell error norms, one figure per flux (split form) --------------
for flux in ('central', 'upwind'):
    run = RUNS[(flux, 'split')]
    if err(run, 'cosi', 'l2') is None:
        continue
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    for norm, lab in (('l2', 'l2'), ('linf', 'max (linf)')):
        t, e = err(run, 'cosi', norm)
        ax[0].plot(t, e, label=lab); ax[1].semilogy(t, e, label=lab)
    ax[0].set(xlim=(0, 3), ylim=(0, 0.04), title='first 3 revolutions (linear)')
    ax[1].set(ylim=(1e-3, 10), title='full run (log)')
    for a in ax:
        a.set_xlabel('time (revolutions)'); a.set_ylabel('relative error'); a.grid(alpha=.3, which='both')
    ax[0].legend()
    fig.suptitle('swtc1-DG cosine bell errors, alpha = pi/4, split form, %s flux' % flux)
    fig.savefig(os.path.join(OUT, 'errors_%s.png' % flux), dpi=110); plt.close(fig)

# ---- split vs conservation form, per flux -----------------------------------
for flux in ('central', 'upwind'):
    rs, rc = RUNS[(flux, 'split')], RUNS[(flux, 'cons')]
    if err(rs, 'cosi', 'l2') is None or err(rc, 'cosi', 'l2') is None:
        continue
    for nrev, logy in ((3, False), (30, True)):
        fig, axs = plt.subplots(2, 2, figsize=(13, 8.6), constrained_layout=True)
        for a, (s, name) in zip(axs.flat, SH):
            for run, lab, c, ls in ((rs, 'split form', C_SPLIT, '-'), (rc, 'conservation form', C_CONS, '--')):
                t, e = err(run, s, 'l2')
                m = t <= nrev + 1e-9
                (a.semilogy if logy else a.plot)(t[m], e[m], ls, color=c, lw=2, label='%s, %s' % (lab, flux))
            a.set_xlim(0, nrev)
            if not logy:
                a.set_ylim(bottom=0)
            a.set(title=name, xlabel='time (revolutions)', ylabel='relative l2 error'); a.grid(alpha=.3, which='both')
        axs[0, 0].legend()
        fig.suptitle('swtc1-DG l2 errors, %s flux: split form vs conservation form (%s)'
                     % (flux, 'first 3 revolutions' if nrev == 3 else '30 revolutions, log axes'))
        fig.savefig(os.path.join(OUT, 'forms_%s_%drev.png' % (flux, nrev)), dpi=110); plt.close(fig)

# ---- tracer variance -------------------------------------------------------
def symlog_axis(a, curves):
    # linear threshold one decade below the smallest typical |value| of any curve (after 1 revolution)
    med = [np.median(np.abs(c[1][c[0] > 1])) for c in curves if np.any(c[0] > 1)]
    med = [m for m in med if m > 0]
    lt = 10.0 ** (np.floor(np.log10(min(med))) - 1) if med else 1e-12
    a.set_yscale('symlog', linthresh=lt); a.axhline(0, color='0.4', lw=.8)

for flux in ('central', 'upwind'):
    vs, vc = var(RUNS[(flux, 'split')]), var(RUNS[(flux, 'cons')])
    if vs is None:
        continue
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.6), constrained_layout=True)
    for k, (a, (s, name)) in enumerate(zip(axs.flat, SH)):
        cur = [(vs[0], vs[1][:, k])]
        a.plot(vs[0], vs[1][:, k], '-', color=C_SPLIT, lw=2, label='split form')
        if vc is not None:
            a.plot(vc[0], vc[1][:, k], '--', color=C_CONS, lw=2, label='conservation form'); cur.append((vc[0], vc[1][:, k]))
        symlog_axis(a, cur)
        a.set(title=name, xlabel='time (revolutions)', ylabel='(V - V0) / V0', xlim=(0, 30)); a.grid(alpha=.3, which='both')
    axs[0, 0].legend()
    fig.suptitle('DG %s flux: relative change in tracer variance V = integral of (g h)^2' % flux)
    fig.savefig(os.path.join(OUT, 'variance_%s.png' % flux), dpi=110); plt.close(fig)

runs4 = [(k, var(v)) for k, v in RUNS.items()]
runs4 = [(k, v) for k, v in runs4 if v is not None]
if runs4:
    for lin in (False, True):
        fig, a = plt.subplots(figsize=(8, 5), constrained_layout=True)
        sty = {('central', 'split'): ('#1f5f8b', '-'), ('central', 'cons'): ('#1f5f8b', '--'),
               ('upwind', 'split'): ('#c0392b', '-'), ('upwind', 'cons'): ('#c0392b', '--')}
        alld = []
        for (flux, form), v in runs4:
            m = v[0] <= 10
            c, ls = sty[(flux, form)]
            a.plot(v[0][m], v[1][m, 0], ls, color=c, lw=2, label='%s, %s' % (flux, FORM[form])); alld.append((v[0][m], v[1][m, 0]))
        if not lin:
            symlog_axis(a, alld)
        a.set(xlim=(0, 10), xlabel='time (revolutions)', ylabel='(V - V0) / V0',
              title='DG cosine bell variance change, first 10 revolutions' + (' (linear axis)' if lin else ''))
        a.grid(alpha=.3, which='both'); a.legend()
        fig.savefig(os.path.join(OUT, 'variance_bell_4runs%s.png' % ('_linear' if lin else '')), dpi=110); plt.close(fig)

# ---- error table at the last common time -----------------------------------
rows = []
for s, name in SH:
    row = [name]
    for key in (('central', 'split'), ('upwind', 'split'), ('central', 'cons'), ('upwind', 'cons')):
        for norm in ('l2', 'linf'):
            e = err(RUNS[key], s, norm)
            row.append('%.3g@%.0f' % (e[1][-1], e[0][-1] * REV) if e is not None else '-')
    rows.append(row)
with open(os.path.join(OUT, 'table.txt'), 'w') as f:
    for r in rows:
        f.write('\t'.join(r) + '\n')
print(open(os.path.join(OUT, 'table.txt')).read())
