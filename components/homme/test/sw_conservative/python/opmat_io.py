"""
Read the operator matrices written by sweqx with SWEQX_OPMAT=<prefix>
(see components/homme/src/sweqx/opmat_mod.F90 and ../swtc1-opmat.sh).

    from opmat_io import load
    d = load(prefix)          # d['L'], d['B'] scipy.sparse CSR (n x n), may be absent
                              # d['lat'], d['lon'] (radians), d['mass'] (assembled GLL mass)
                              # d['meta'] dict, d['Ltest'] = (x, y) with y = L x from sweqx
"""
import glob
import numpy as np
import scipy.sparse as sp

DOF = np.dtype([('dof', '<i4'), ('gdof', '<i8'), ('lat', '<f8'), ('lon', '<f8'), ('mass', '<f8')])
TRIP = np.dtype([('row', '<i4'), ('col', '<i4'), ('val', '<f8')])


def read_meta(prefix):
    meta = {}
    with open(prefix + '.meta') as f:
        for line in f:
            if '=' in line:
                k, v = line.split('=', 1)
                v = v.strip()
                for typ in (int, float):
                    try:
                        v = typ(v)
                        break
                    except ValueError:
                        pass
                meta[k.strip()] = v
    return meta


def _files(prefix, kind):
    return sorted(glob.glob(f'{prefix}.{kind}.r[0-9][0-9][0-9][0-9].bin'))


def load(prefix, ops=('L', 'B')):
    meta = read_meta(prefix)
    n, nlev = meta['n'], meta['nlev']
    d = {'meta': meta, 'n': n}
    dofs = np.concatenate([np.fromfile(f, dtype=DOF) for f in _files(prefix, 'dofs')])
    if len(dofs) != n or not np.array_equal(np.sort(dofs['dof']), np.arange(n)):
        raise ValueError(f'{prefix}: dof files do not cover 0..{n-1} exactly once')
    order = np.argsort(dofs['dof'])
    dofs = dofs[order]
    for k in ('gdof', 'lat', 'lon', 'mass'):
        d[k] = dofs[k].copy()
    for op in ops:
        files = _files(prefix, op)
        if not files:
            continue
        t = np.concatenate([np.fromfile(f, dtype=TRIP) for f in files])
        A = sp.coo_matrix((t['val'], (t['row'], t['col'])), shape=(n, n)).tocsr()
        if A.nnz != len(t):
            raise ValueError(f'{prefix}.{op}: duplicate entries')
        d[op] = A
        tdt = np.dtype([('dof', '<i4'), ('x', '<f8', (nlev,)), ('y', '<f8', (nlev,))])
        tt = np.concatenate([np.fromfile(f, dtype=tdt) for f in _files(prefix, op + 'test')])
        tt = tt[np.argsort(tt['dof'])]
        d[op + 'test'] = (tt['x'], tt['y'])
    return d
