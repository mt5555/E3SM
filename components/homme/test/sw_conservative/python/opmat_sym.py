"""
Block diagonalize swtc1 operators with the symmetries of the cubed sphere grid
that also preserve the solid-body rotation:

  R, a rotation about a grid axis:
    alpha = 0    : 90 deg about the z axis                              (order m = 4)
    alpha = pi/4 : 180 deg about the rotation axis a = (-1,0,1)/sqrt(2),
                   which passes through the midpoints of two cube edges  (order m = 2)
  J, the antipodal map x -> -x (order 2).  It maps the grid onto itself and
    commutes with any solid-body rotation, since u(-x) = -u(x).

R and J commute, so they generate the abelian group Z_m x Z_2 (8 elements for
alpha=0, 4 for alpha=pi/4).  If the node permutations P_R and P_J commute with
L, then L maps each joint eigenspace of (P_R, P_J), eigenvalues (w^k, s) with
w = exp(2 pi i/m) and s = +-1, into itself.  The spectrum of L is then the union
of the spectra of the blocks L_ks = Q_ks^H L Q_ks, where the orthonormal columns
of Q_ks span that eigenspace.  Each block is about n/(2m) in size, so dense
eigenvalues cost about (2m)^2 times less than for the whole matrix.  Q_ks is
real when w^k is real; for real L the block (m-k, s) is the complex conjugate of
(k, s), so for alpha=0 the blocks k=3 are taken from k=1 without a solve.

inversion=False keeps only R (the original Z_m blocks).
"""
import numpy as np
import scipy.linalg
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
        raise ValueError(f'grid is not invariant under the map (max dist {dist.max():.2e})')
    return perm


def generators(lat, lon, alpha, inversion=True):
    """[(perm, order, name)] of the group generators"""
    R, m = rotation(alpha)
    gens = [(node_permutation(lat, lon, R), m, 'R')]
    if inversion:
        gens.append((node_permutation(lat, lon, -np.eye(3)), 2, 'J'))
    return gens


def symmetry_blocks(lat, lon, alpha, inversion=True):
    """(gens, blocks): blocks is a list of (label, Q, src) with Q sparse n x n_k,
    orthonormal columns, P_R Q = w^k Q and P_J Q = s Q, label = k + m*(s == -1),
    and src the label of the block whose complex conjugate this one is (None if
    it must be solved).  (P x)[perm[j]] = x[j]."""
    gens = generators(lat, lon, alpha, inversion)
    n = len(lat)
    orders = [o for _, o, _ in gens]
    # image of every node under every group element g = prod_i G_i^{a_i}
    expo = np.array(np.meshgrid(*[np.arange(o) for o in orders], indexing='ij')).reshape(len(gens), -1).T
    img = np.empty((n, len(expo)), dtype=np.int64)
    for c, a in enumerate(expo):
        j = np.arange(n)
        for (perm, _, _), ai in zip(gens, a):
            for _ in range(ai):
                j = perm[j]
        img[:, c] = j
    reps = np.flatnonzero(img.min(axis=1) == np.arange(n))     # one node per orbit
    m = orders[0]
    blocks = []
    for lab in range(int(np.prod(orders))):
        k, b = lab % m, lab // m                                # b = 1: s = -1
        if m == 4 and k == 3:
            blocks.append((lab, None, 1 + m * b))
            continue
        # basis vector of orbit r: sum_g conj(chi(g)) e_{g r}, chi(R^a J^c) = w^(k a) s^c
        chi = np.exp(2j * np.pi * k * expo[:, 0] / m)
        if len(gens) > 1:
            chi = chi * (-1.0) ** (b * expo[:, 1])
        rows = img[reps].ravel()
        cols = np.repeat(np.arange(len(reps)), len(expo))
        vals = np.tile(np.conj(chi), len(reps))
        Q = sp.csc_matrix((vals, (rows, cols)), shape=(n, len(reps)))
        Q.sum_duplicates()
        nrm = np.sqrt(np.asarray(abs(Q).power(2).sum(axis=0))).ravel()
        keep = nrm > 1e-8
        Q = Q[:, keep] @ sp.diags(1 / nrm[keep])
        if (2 * k) % m == 0:
            Q = Q.real
        Q = sp.csc_matrix(Q)
        Q.eliminate_zeros()
        blocks.append((lab, Q, None))
    # fill in the conjugate blocks
    byl = {lab: Q for lab, Q, _ in blocks}
    blocks = [(lab, Q if src is None else byl[src].conj().tocsc(), src) for lab, Q, src in blocks]
    assert sum(Q.shape[1] for _, Q, _ in blocks) == n
    return gens, blocks


def check_commutes(A, perm):
    """||P A - A P||_F / ||A||_F with (P x)[perm[j]] = x[j]"""
    n = len(perm)
    P = sp.csr_matrix((np.ones(n), (perm, np.arange(n))), shape=(n, n))
    return sp.linalg.norm(P @ A - A @ P) / sp.linalg.norm(A)


def block_matrix(A, Q):
    """dense block Q^H A Q and the relative norm of the neglected off-block part"""
    AQ = A @ Q
    Ak = (Q.conj().T @ AQ).toarray()
    # A Q should lie in span(Q): residual of the projection
    res = AQ - Q @ sp.csc_matrix(Ak)
    return Ak, sp.linalg.norm(res) / max(sp.linalg.norm(AQ), 1e-300)


def eig_unbalanced(A, maxscale=1e3, seed=0):
    """np.linalg.eig, except when LAPACK's balancing would rescale A strongly (a
    range above maxscale, e.g. 7e7 in some conservation-form blocks): its eigenvectors
    then have residuals ||A v - lam v|| up to ~1e-9 ||A||.  A random orthogonal
    similarity Z^T A Z defeats the balancing and restores residuals ~1e-15."""
    s, _ = scipy.linalg.matrix_balance(A, permute=False, separate=True)[1]
    if s.max() / s.min() <= maxscale:
        return np.linalg.eig(A)
    Z, _ = np.linalg.qr(np.random.default_rng(seed).standard_normal(A.shape))
    lam, U = np.linalg.eig(Z.T @ A @ Z)
    return lam, Z @ U


def block_eigs(A, lat, lon, alpha, vectors=False, inversion=True, hermitian=None, driver='evd'):
    """generator over the symmetry blocks of A, one at a time:
    yields (label, lam, Q, U, info), with A (Q U) = (Q U) diag(lam) (U None without
    vectors).  hermitian: optional function Ak -> Hermitian matrix H, then lam = -i
    eig(H) (for i K with K skew: the eigenpairs of the skew part).  info holds the
    commutation errors of the generators, the off-block residual and the block size."""
    gens, blocks = symmetry_blocks(lat, lon, alpha, inversion)
    comm = max(check_commutes(A, perm) for perm, _, _ in gens)
    done = {}
    for lab, Q, src in blocks:
        if src is not None:
            lam, U, off = done.pop(src)
            lam, U = np.conj(lam), (None if U is None else np.conj(U))
        else:
            Ak, off = block_matrix(A, Q)
            if hermitian is not None:
                H = hermitian(Ak)
                del Ak
                # in place; driver 'evd' (divide and conquer, orthogonal to 1e-15) needs ~2 n^2
                # complex of workspace, 'evr' (MRRR) O(n) but orthogonal only to ~1e-11
                r = scipy.linalg.eigh(H, eigvals_only=not vectors, overwrite_a=True, check_finite=False,
                                      driver=driver)
                del H
                w, U = r if vectors else (r, None)
                lam = -1j * w
            else:
                if vectors:
                    lam, U = eig_unbalanced(Ak)
                else:
                    lam, U = np.linalg.eigvals(Ak), None
                del Ak
            if any(s == lab for _, _, s in blocks):
                done[lab] = (lam, U, off)
        yield lab, lam, Q, U, dict(commute=comm, offblock=off, size=Q.shape[1], m=len(blocks))


def eigvals_by_symmetry(A, lat, lon, alpha, vectors=False, inversion=True):
    """eigenvalues (and optionally eigenvectors, n x n) of A using the symmetry blocks.
    returns lam, V (or None), info dict, symmetry class (block label) of each eigenvalue"""
    lams, Vs, ks = [], [], []
    off, sizes = 0.0, []
    for lab, lam, Q, U, inf in block_eigs(A, lat, lon, alpha, vectors, inversion):
        lams.append(lam); ks.append(np.full(len(lam), lab))
        if vectors:
            Vs.append(Q @ U)
        off = max(off, inf['offblock']); sizes.append(inf['size'])
    lam = np.concatenate(lams)
    V = np.concatenate(Vs, axis=1) if vectors else None
    info = {'m': len(sizes), 'commute': inf['commute'], 'offblock': off, 'sizes': sizes}
    return lam, V, info, np.concatenate(ks)
