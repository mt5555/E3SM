"""
Block diagonalize swtc1 operators with a rotational symmetry of the cubed
sphere grid that also preserves the solid-body rotation:

  alpha = 0    : rotation by 90 deg about the z axis      (group Z4)
  alpha = pi/4 : rotation by 180 deg about the rotation axis a = (-1,0,1)/sqrt(2),
                 which passes through the midpoints of two cube edges (group Z2)

If P (the permutation of the GLL nodes induced by the rotation) commutes with L,
then L maps each eigenspace of P (eigenvalue w^k, w = exp(2 pi i/m)) into
itself, so the spectrum of L is the union of the spectra of the m blocks
L_k = Q_k^H L Q_k, where the orthonormal columns of Q_k span that eigenspace.
Dense eigenvalues then cost ~m^2 times less.
"""
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree


def xyz(lat, lon):
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=1)


def rotation(alpha):
    """(R, m): rotation matrix and group order for the swtc1 rotation axis alpha"""
    if abs(alpha) < 1e-12:
        c, s = 0.0, 1.0      # 90 deg about z
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]]), 4
    if abs(alpha - np.pi / 4) < 1e-12:
        a = np.array([-np.sin(alpha), 0, np.cos(alpha)])
        return 2 * np.outer(a, a) - np.eye(3), 2    # 180 deg about a
    raise ValueError('no symmetry implemented for this alpha')


def node_permutation(lat, lon, R, tol=1e-9):
    """perm[j] = index of the node at R x_j"""
    X = xyz(lat, lon)
    dist, perm = cKDTree(X).query(X @ R.T)
    if dist.max() > tol or len(np.unique(perm)) != len(perm):
        raise ValueError(f'grid is not invariant under the rotation (max dist {dist.max():.2e})')
    return perm


def symmetry_blocks(lat, lon, alpha):
    """list of (k, Q_k) with Q_k sparse n x n_k, orthonormal columns, P Q_k = w^k Q_k"""
    R, m = rotation(alpha)
    perm = node_permutation(lat, lon, R)
    n = len(lat)
    w = np.exp(2j * np.pi / m)
    # orbits of the nodes
    seen = np.zeros(n, bool)
    orbits = []
    for j in range(n):
        if seen[j]:
            continue
        orb = [j]
        while perm[orb[-1]] != j:
            orb.append(perm[orb[-1]])
        seen[orb] = True
        orbits.append(orb)
    blocks = []
    for k in range(m):
        rows, cols, vals = [], [], []
        ncol = 0
        for orb in orbits:
            o = len(orb)
            # P^s e_r = e_{orb[s mod o]};  sum_s w^{-ks} P^s e_r is nonzero iff w^{k o} = 1
            if (k * o) % m:
                continue
            for s in range(o):
                rows.append(orb[s]); cols.append(ncol); vals.append(w ** (-k * s) / np.sqrt(o))
            ncol += 1
        Q = sp.csc_matrix((vals, (rows, cols)), shape=(n, ncol))
        if m == 2 or k % 2 == 0:
            Q = Q.real.tocsc()
        blocks.append((k, Q))
    assert sum(Q.shape[1] for _, Q in blocks) == n
    return perm, blocks


def check_commutes(A, perm):
    """||P A - A P||_F / ||A||_F with (P x)[perm[j]] = x[j]"""
    n = len(perm)
    P = sp.csr_matrix((np.ones(n), (perm, np.arange(n))), shape=(n, n))
    return sp.linalg.norm(P @ A - A @ P) / sp.linalg.norm(A)


def block_matrices(A, blocks):
    """dense blocks Q_k^H A Q_k and the relative norm of the neglected off-block part"""
    out = []
    off = 0.0
    for k, Q in blocks:
        AQ = A @ Q
        Ak = (Q.conj().T @ AQ).toarray()
        # A Q_k should lie in span(Q_k): residual of the projection
        res = AQ - Q @ sp.csc_matrix(Ak)
        off = max(off, sp.linalg.norm(res) / max(sp.linalg.norm(AQ), 1e-300))
        out.append((k, Ak))
    return out, off


def eigvals_by_symmetry(A, lat, lon, alpha, vectors=False):
    """eigenvalues (and optionally eigenvectors, n x n) of A using the symmetry blocks.
    returns lam, V (or None), info dict, symmetry class k of each eigenvalue"""
    perm, blocks = symmetry_blocks(lat, lon, alpha)
    comm = check_commutes(A, perm)
    mats, off = block_matrices(A, blocks)
    m = len(blocks)
    lams, Vs, ks = [], [], []
    for (k, Ak), (_, Q) in zip(mats, blocks):
        if m == 4 and k == 3:
            # conjugate of block 1 for real A
            l1 = lams[1]
            lams.append(np.conj(l1)); ks.append(np.full(len(l1), 3))
            if vectors:
                Vs.append(np.conj(Vs[1]))
            continue
        if vectors:
            lam, U = np.linalg.eig(Ak)
            Vs.append(Q @ U)
        else:
            lam = np.linalg.eigvals(Ak)
        lams.append(lam); ks.append(np.full(len(lam), k))
    lam = np.concatenate(lams)
    V = np.concatenate(Vs, axis=1) if vectors else None
    info = {'m': m, 'commute': comm, 'offblock': off, 'sizes': [Q.shape[1] for _, Q in blocks]}
    return lam, V, info, np.concatenate(ks)
