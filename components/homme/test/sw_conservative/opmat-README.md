# sweqx operator matrices: extraction, eigenvalue analysis and status

This note documents the operator-matrix tools on branch `mt5555/homme/sweqx-operator-matrix`
and serves as the handoff for continuing the eigenvalue work in a new session. It covers
what the tools do, how to run them, what has been found so far, and the open questions.

## Background: the convergence question

The tools were written during a convergence study of swtc1 (solid-body advection of four
shapes on the cubed sphere, np=4 / p=3, classical RK4, dt ∝ dx, one 12-day revolution):

- With no hyperviscosity, the Gaussian converges at 4th order only when the flow follows the
  grid (α=0: l2 4.02, max 3.93). For oblique flow it loses order: α=0.6 gives 3.60 / 3.39,
  and α=π/4 (whose path runs along two cube edges and through four cube corners) gives
  3.51 / 3.04. p=6 at α=π/4 gives 6.0 (= p). The split and conservation forms have the same
  orders, and a discretely divergence-free velocity (DSS of curl ψ) does not change them.
- Hyperviscosity ν = ν₃₀ (30/NE)^k: k=3 gives exactly 3; k=4 with ν₃₀=1e13 (close to the
  error-minimizing ν) gives 3.71; ν₃₀=1e14 gives a clean 4.00 because its own dx⁴ error
  dominates. A fitted error model predicts that the error-minimizing ν scales like dx^3.75 and
  gives order ≈ 3.75.
- Independent Python SEM codes (private repo `mt5555/sem-advection-tests`) show that p=3 SEM
  keeps 4th order on uniform and smoothly mapped meshes (curved, non-orthogonal, oblique flow:
  4.01 at N=256→512), and loses a few tenths where the mesh mapping is not smooth (kinked grid
  lines, skewed patches); a single valence-3 vertex costs nothing on its own. The cubed
  sphere is smooth inside faces and non-smooth at the cube edges, which fits the swtc1 pattern.
  HOMME's α=π/4 max order (3.0) is still lower than any comparable planar case (≥3.6).

The question the eigenvalue analysis addresses: is HOMME's discrete operator behaving like a
correct SEM operator, and which modes are responsible for the loss of order?

## What the tools do

`src/sweqx/opmat_mod.F90` extracts, as sparse matrices:

- **L**: the swtc1 transport tendency, dp/dt = L p, exactly as `compute_and_apply_rhs` computes
  it (called with dt2=0, pmean=0, so it returns the DSS'd, mass-divided tendency in `ptens`;
  `advance_mod.F90` is not modified for this), for the prescribed velocity and the namelist
  `theta_advect_form`.
- **B**: the hyperviscosity operator of `advance_hypervis` with `hypervis_order=2`,
  B p = rspheremp·DSS(biharmonic_wk(p)), a discrete ∇⁴ (positive semidefinite).

Columns are probed with continuous delta functions at each unique GLL node. The 4 levels are
independent, and with graph coloring many far-apart columns are probed per call: 384 operator
applications for L and 864 for B, independent of NE. Rows and columns use HOMME's contiguous
unique-point ids (0-based); n = 6·NE²·(np−1)² + 2.

Related option (commit 68d165ee): `SWEQX_TC1_VELOCITY=curlpsi` replaces the analytic swtc1
velocity by u = DSS(curl ψ) for the rigid-body stream function, which has zero DSS'd discrete
divergence (about 1e-18). Default runs are unchanged (bit-for-bit checked at NE8).

### Switches (environment variables, so the shared namelist code is untouched)

| Variable | Meaning |
|---|---|
| `SWEQX_OPMAT=<prefix>` | extract, write `<prefix>.*`, then stop (only in `sweq_rk`, i.e. `integration="runge_kutta"`) |
| `SWEQX_OPMAT_OPS=L\|B\|LB` | operators to extract (default LB) |
| `SWEQX_OPMAT_COLOR=0` | one column per level per call (no coloring), for checking |
| `SWEQX_TC1_VELOCITY=analytic\|curlpsi` | swtc1 prescribed velocity |

Requirements checked at run time: 1 thread, `kmass=-1`, `test_case='swtc1'` for L,
`hypervis_order=2` for B.

### Running

Build sweqx in its own directory with `-DSWEQX_PLEV=4 -DSWEQX_NP=4` (see `swtc1.sh`). On a
fresh clone this needs the submodules cime, externals/scorpio and externals/ekat and their
nested submodules (cime/CIME/non_py/cprnc, externals/ekat/extern/kokkos):
`git -c url.https://github.com/.insteadOf=git@github.com: submodule update --init --depth=1 <path>`.
A configure run before cprnc is fetched leaves a stub `Macros.cmake` in its directory; delete
it before fetching. Then:

```bash
# test/sw_conservative/swtc1-opmat.sh NE subcase adv_form velocity [ops] [nranks]
#   subcase 1 (alpha=0) or 4 (alpha=pi/4); adv_form 0 (conservation) or 2 (split)
bld=~/scratch/opmat-build HOMME=~/codes/acme/components/homme outdir=~/scratch/opmat \
  ./swtc1-opmat.sh 15 4 2 analytic LB 4
```

Output per MPI rank: `<prefix>.L.r####.bin` and `.B.r####.bin` (binary int32 row, int32 col,
float64 value), `.dofs.r####.bin` (int32 dof, int64 gdofP, float64 lat, lon, assembled mass),
`.Ltest.r####.bin` (y = L x for test fields, to validate), and a `.meta` file. NE30 L takes
about 25 s and B about 2 min on 4 ranks. File sizes at NE30: L about 6.7 MB, B about 90 MB.
L at NE60 takes about 70 s and NE120 about 5 min on 2 ranks (NE120 L: 3.5 M entries per rank).
B takes 2 min at NE60 and 9 min at NE120 on 4 ranks (NE120 B: 94.5 M entries, 1.5 GB).

### Python (test/sw_conservative/python/)

- `opmat_io.py`: `load(prefix)` returns `L`, `B` (scipy CSR), `lat`, `lon`, `mass`, `meta`, `Ltest`.
- `opmat_sym.py`: block-diagonalizes L with the symmetries of grid and flow: a rotation
  (90° about z for α=0, order 4; 180° about the α=π/4 rotation axis, which passes through
  two cube-edge midpoints, order 2) and the antipodal map x → −x, which maps the grid onto
  itself and commutes with any solid-body rotation. That gives 8 blocks for α=0 and 4 for
  α=π/4, each about n/8 or n/4. `block_eigs()` builds and solves one block at a time (conjugate
  blocks come from their partner without a solve), optionally for the Hermitian form of the
  skew part. Both generators commute with L, B and every case to roundoff (1e-16 at NE8–20,
  1e-15 at NE30). Block eigenvalues match full dense ones to 3e-12 ω (NE8, NE10). LAPACK's
  balancing rescales some conservation-form blocks by up to 7e7, which left eigenvector
  residuals near 1e-9 ‖L‖; `eig_unbalanced()` sidesteps it with a random orthogonal similarity.
- `opmat_sym_test.py --mats DIR --test blocks|full|old|vectors|proj`: checks of the blocks
  (orthonormality, symmetry eigenspaces, off-block residual), eigenvalues against full dense
  solves and another symmetry module, eigenpair residuals, and `opmat_proj.py` results against
  earlier runs. With the antipodal map added, NE16 and NE20 reproduce the earlier modal budgets
  (error fields to 3e-9, per-mode error and diagnostics of isolated modes to 1e-8).
- `opmat_eigs.py --mats DIR [--mats30 DIR30] --out OUT [--ne 8 15]`: the full analysis
  (checks, spectra, accuracy against −i m ω, localization, hyperviscosity), cached in
  `OUT/cache`, plots and `results.txt` in `OUT`. NE15 takes about 10 minutes on 4 cores.
  Set `OPENBLAS_NUM_THREADS=1` if the cores are shared: oversubscribed OpenBLAS was about
  15× slower.

- `opmat_proj.py --mats DIR --out OUT --ne 8 16 20 --eig 8 16 20`: the Gaussian's day-12
  error from L: RK4 with the convergence-study time step (dt = 1800 s·30/NE, 7200 s at NE8),
  which reproduces the HOMME runs; the semi-discrete error exp(LT)p0 − p0 (any NE, sparse); and
  for the NE in `--eig` its decomposition over the eigenvectors of the skew part K of
  M^½ L M^-½ (orthonormal, so the per-mode budget adds up exactly). Saves
  `OUT/<case>-proj.npz` with the error fields, per-mode amplitudes, eigenvalues and diagnostics,
  the error from resolved modes, and the 8 largest error-carrying eigenvectors. Dense cost with
  `OPENBLAS_NUM_THREADS=1`: NE16 2.3 min (α=0) and 7.5 min (α=π/4), NE20 8 min and 30 min
  (four blocks of about 5,460, peak 3.3 GB). With 4 threads, one NE20 α=π/4 block takes 2 min.
  Peak memory is about 4 × 16 n_k² bytes per block (n_k the block size), so NE30 α=π/4
  (blocks of about 12,150) needs about 10 GB; expected time about 1.5 h on 4 threads.
- **Caveat on m\***: `mstar` and `purity` (in `opmat_eigs.py` and `opmat_proj.py`) project a
  mode onto e^{imφ} summed over the whole sphere. For about half the modes (those odd about the
  plane perpendicular to the rotation axis) that sum vanishes, and m\* is roundoff noise; the
  purity can still come out above 0.5. Results that select modes by m\* include such modes.
- `opmat_source.py --mats DIR --out OUT --ne 30 60`: integrates de/dt = L e + τ(t) with the
  truncation error τ = L p_ex − ∂ₜp_ex restricted to cube corners, cube-edge nodes, the
  one-element layer beside the edges and the interior, to show where the error is generated
  (the parts add up to the semi-discrete error). `--history` writes ‖p(t) − p_ex(t)‖ through
  the revolution and the distance of the Gaussian's center to the cube edges and corners.
  NE60 takes 10 min, NE120 75 min (history 13 min).
- `opmat_nuopt.py --mats DIR --out OUT --ne 30 60 [--sc 4] [--coarse 11 14.5 15]`: scans the
  hyperviscosity coefficient with HOMME's step (I − dt ν B) R4(dt L), many ν per pass as the
  columns of one block (a coarse pass in log ν₃₀, then a fine pass around the l2 and max
  minima), and writes the day-12 Gaussian errors per ν to `OUT/<case>-nuopt.npz` and
  `OUT/nuopt.txt`. B is read from the `sc1-f2-analytic` files.
- `opmat_proj_plot.py --proj DIR [DIR ...] --out OUT [--conv RUNS] [--source SRC]`: tables
  (`results_proj.txt`) and figures: convergence, error maps with the cube edges, error vs
  distance to the edges, error history scaled by NE⁴, Gaussian amplitude per mode, error per
  frequency band, and maps of the top error-carrying modes.

### HOMME details found along the way

- `elem%gdofP` is not contiguous (it ranges up to nelem·np²); use the unique-point ids.
- `elem%spheremp` is not the assembled mass; the assembled mass is `1/elem%rspheremp`.

## Spectrum results (NE8, NE15 dense; NE30 extremal eigenvalues with ARPACK)

ω = u0/rearth = 6.060e-6 s⁻¹ (one revolution per 12 days). Exact spectrum of the solid-body
rotation operator: −i m ω, m integer, the same for every α.

**Validation against HOMME runs (NE30, α=π/4)**

| Quantity | From the matrices | Observed |
|---|---|---|
| RK4 dt limit, ν=0 | 2329.4 s | stable 2330, unstable 2340 |
| hyperviscosity limit at dt=1800 | ν ≤ 2.3e14 | 2.25e14 runs, 2.5e14 blows up |
| dt limit at ν=1e15 | 417 s | 360 s |
| RK4 onset with ν=1e14 (first-order estimate) | about 2610 s | 2580–2590 s |

**Structure**

- B: M^½ B M^-½ symmetric to 1.8e-16, smallest eigenvalue ≈ 0 (constants), matches
  (l(l+1)/a²)² to 1e-4 up to l=20 (NE8) / l=40 (NE15), then 0.85–1.3, rising to 3.4× at the
  grid scale.
- Split form, curl ψ velocity: skew in the mass inner product to roundoff
  (‖S+Sᵀ‖/‖S‖ = 1–6e-15 at NE8/15/30), max Re λ ≈ 1e-14 ω.
- Split form, analytic velocity: ‖S+Sᵀ‖/‖S‖ = 9e-6 / 7.5e-7 / 4.7e-8 (NE8/15/30), consistent
  with a symmetric part ½·DSS(m·div_h u) from the discrete divergence of the analytic velocity.
  α=0: spectrum imaginary to roundoff. α=π/4: max Re λ = 2.7e-5 ω (NE8), 4.4e-6 ω (NE15);
  these are the largest-|λ| modes, localized just inside the cube corners.
- Conservation form: ‖S+Sᵀ‖/‖S‖ = 4.5e-2 / 2.5e-2 / 1.2e-2, even with the curl ψ velocity.
  α=π/4: max Re λ = 0.13 ω (NE8), 0.18 ω (NE15): grid-scale growing modes along element
  edges, not corner-localized. α=0: about 0 at NE8, 0.074 ω at NE15. ν₃₀=1e13 removes all
  growth at NE15. This explains the conservation-form blow-ups in long runs.

**Stability**: max|λ| is about 33 / 67 / 142 ω (α=0) and 47 / 95 / 200 ω (α=π/4) at
NE8 / 15 / 30, a ratio of about √2 set by the grid orientation, nearly independent of form and
velocity. RK4 limit 2.83/max|λ|: 3292 s (α=0) and 2329 s (α=π/4) at NE30.

**Accuracy (NE8)**: |m| ≤ 3 errors 1e-9 to 1e-5 ω; about 1e-2 ω at |m|=8, 0.1 at |m|=16, O(1)
beyond |m|≈25. At α=0 many high-l, small-m (near-zonal) modes stay accurate; at α=π/4 the
error grows smoothly with |Im λ| and the spectrum reaches √2 further.

**Hyperviscosity**: at NE15, ν₃₀=1e13 already makes max Re λ of L − νB roundoff for both α.
At NE8, first-order damping per revolution for modes with error > 0.1 is about 5e-3 at
ν₃₀=1e13; accurate modes keep 0.98 / 0.92 / 0.86 per revolution at ν₃₀ = 1e13 / 5.4e13 / 1e14
(α=π/4, split, analytic). At α=0 the accurate near-zonal high-l modes are damped strongly
(median about 0.2 per revolution at ν₃₀=1e13).

Caveats: m* comes from a mass-weighted Fourier fit, reliable only for modes with more than 50%
of their power in one m; NE15 has eigenvalues only (eigenvectors at NE8); the NE30
hyperviscosity numbers are first-order Rayleigh-quotient estimates from the 150 largest-|λ|
modes.

## Gaussian error budget (open question 1, October 2026)

Setup as in the convergence study: Gaussian (level 2 of swtc1), split form, analytic velocity,
ν=0, RK4 with dt = 1800 s·30/NE; error at day 12 (one revolution), relative l2 with the
assembled mass and relative max, as in `tc1_errors`.

**The matrices reproduce the runs.** Stepping the extracted L with RK4 gives HOMME's day-12
errors to 6 digits at NE 8, 15, 20, 30 and 60 (both α). Matrix RK4 NE60→120 orders: α=0
4.02 / 3.93 (l2 / max), α=π/4 3.51 / 3.04, the measured values.

**The order loss is spatial.** The semi-discrete error exp(LT)p0 − p0 (no time stepping) has
the same orders: α=π/4 3.53 / 3.02 (NE30→60) and 3.51 / 3.05 (NE60→120); α=0 4.06 / 4.13 and
4.03 / 3.84. RK4 changes the error field by about 25% but not its norm. Orders at NE ≤ 20 are
pre-asymptotic (4–7), so dense eigenvectors (NE ≤ 20 here) cannot show the asymptotic order by
themselves. Odd NE is special at α=0: the equator (the path) runs through element centers
instead of along element edges, and NE15 is 2× more accurate than NE16; compare even NE only.

**Modal form of the prediction.** With e^{−i m ω T} = 1 for every m at T = one revolution, the
semi-discrete error is Σ_j c_j (e^{λ_j T} − 1) v_j: no wavenumber assignment is needed. The
eigenvectors of L itself are far from orthogonal in near-degenerate clusters (Σ|c_j|² ≈ 3‖p0‖²),
so the decomposition uses the skew part K of S = M^½ L M^-½ instead (Hermitian iK, `eigh`):
orthonormal to 1e-14, and K alone reproduces L's error to 2e-4 (α=0) and 3e-3 (α=π/4).
Then ‖e‖² = Σ_j a_j² |e^{λ_j T} − 1|² exactly, with a_j = |c_j|/‖p0‖. "Resolved" below means
|e^{λ_j T} − 1| < 0.1 (phase error < 0.016 revolutions).

| NE | α | share of ‖e‖² from unresolved modes | modes for 90% of ‖e‖² | E²-weighted edge enrichment |
|---|---|---|---|---|
| 8 | 0 / π/4 | 0.98 / 0.96 | 39 / 309 | 0.67 / 1.11 |
| 12 | 0 / π/4 | 0.89 / 0.91 | 81 / 824 | 0.70 / 1.27 |
| 16 | 0 / π/4 | 0.81 / 0.96 | 63 / 888 | 0.70 / 1.61 |
| 20 | 0 / π/4 | 0.57 / 0.98 | 53 / 844 | 0.69 / 2.42 |

- α=0: the error moves to resolved modes as NE grows. Its largest contributors at NE20 are
  smooth equatorial waveguide modes with eigenvalues close to integers (−8.965i ω, +5.990i ω):
  ordinary phase error, which converges at 4th order.
- α=π/4: unresolved modes keep 96–98% of the error, and their error converges more slowly
  with each step (orders 6.7, 6.1, 4.4, 3.9 for NE 8→10→12→16→20), while the resolved part
  converges at about 5.6. The error-carrying modes become more concentrated on the cube edges.
- **Edge-trapped modes.** At NE16 and NE20 the largest single contributors at α=π/4 are a
  conjugate pair of modes trapped on exactly the two cube edges the path follows, corner to
  corner (lon 315→45 at the top, 135→225 at the bottom), grid-scale across the edge. Their
  eigenvalue does not scale with 1/h: ±3.654i, 3.719i, 3.738i ω at NE10, 16, 20 (3.615 at NE8;
  odd NE15 has a similar pair at 5.58i). The 0.5% most edge-concentrated modes carry 0.1%,
  3% and 16% of ‖e‖² at NE8, 16 and 20. A spurious mode with an O(ω) frequency is not averaged
  out by the moving Gaussian, so the truncation error that drives it accumulates.
- Removing the unresolved modes from the initial condition does not help: the removed part is
  itself an error at day 12 of the same size (NE12 α=π/4: 2.3e-3, against a day-12 error of
  1.2e-3).

**Where and when the error appears (sparse, up to NE120).**

- The truncation error τ is O(h³) pointwise everywhere for both α (max order 2.97–3.06 on the
  edges, near them and in the face interiors, with the Gaussian centered on an edge), as for
  p=3 CG SEM on any mesh (a 1-D periodic check gives τ order 3.0 and solution order 4.3–4.4).
  4th order at α=0 comes from cancellation, not from a small τ. On and next to the edges, τ's
  l2 converges at 3.5 (O(h³) amplitude on an O(h)-wide strip).
- History (‖e(t)‖·NE⁴, NE16–120): at α=π/4 the curves collapse until the Gaussian reaches the
  first cube corner (day 1.8), separate only during the two transits along cube edges and
  corners (days 1.8–4.2 and 7.8–10.2), and stay parallel while the Gaussian crosses face
  interiors. At α=0 (edges crossed at right angles) they collapse for the whole revolution.
- Maps: at α=π/4 the max error sits on the two transit edges (on a corner at NE30, on an edge
  at NE60 and NE120), and thin streaks run from them across the faces. More than 10° from any
  edge (62% of the sphere) the error still converges at l2 3.52 / max 2.85 (NE60→120), but the
  area holding half of that error's energy shrinks 2.4–3.3× per doubling: an O(h³) amplitude on
  O(h)-wide rays, whose l2 the skew operator conserves after they are generated. This is
  consistent with l2 order 3.5 and max order 3.
- Source split (where τ is generated), α=π/4, l2 / max orders NE60→120: corners 4.55 / 4.29,
  cube-edge nodes 3.53 / 2.94, near-edge layer 3.43 / 2.72, interior 3.41 / 2.62 (α=0, NE30→60:
  all regions 4.1–5.4). The interior part is slow too, so the failure of cancellation is a
  property of the operator's modes near the edges (the edge-trapped modes above), not of the
  nodes on which τ sits.

**Interpretation.** Cancellation of the O(h³) truncation error fails while the Gaussian moves
along a cube edge, because the operator has edge-trapped grid-scale modes with O(ω)
frequencies that the moving truncation error drives. That leaves an O(h³) error in an
O(h)-wide band that is then carried along rays; energy conservation of the split-form operator
keeps its l2 at O(h^3.5). Velocity changes at truncation-error level cannot remove it (τ
involves higher derivatives of p; the curl ψ velocity already left the orders unchanged).

## Error-minimizing hyperviscosity (October 2026)

**Heuristic.** With hyperviscosity ν ∝ h^q, damping limits the resonant edge error to
O(h^{3+4−q}) pointwise while the hyperviscosity's own error on the smooth field is O(h^q);
balancing the l2 terms (O(h^{7.5−q}) against O(h^q)) gives q = 3.75 and order 3.75, the values
of the convergence study's fitted error model. The same balance for the max norm
(O(h^{7−q}) against O(h^q)) gives 3.5.

**Measured** with `opmat_nuopt.py`: HOMME's classic-RK4 step with hyperviscosity,
p ← (I − dt ν B) R4(dt L) p (`advance_nonstag_rk4`, then `advance_hypervis` with
`hypervis_subcycle=1`, ν_s = ν), applied to many ν at once. It reproduces a sweqx run with ν to
6 digits (NE20, α=π/4, ν=2.85e13: l2 7.18577e-5, max 1.99942e-4). Gaussian, α=π/4, split form,
analytic velocity, dt = 1800 s·30/NE, error at day 12; the optimum is a parabola fit of
log(error) against log ν around the minimum of the scan. ν₃₀ = ν·(NE/30)⁴.

| NE | l2-optimal ν₃₀ | minimum l2 | l2 at ν=0 | gain | minimum max | max at ν=0 |
|---|---|---|---|---|---|---|
| 16 | 5.5e12 | 1.841e-4 | 3.315e-4 | 1.8 | 4.26e-4 | 6.47e-4 |
| 20 | 5.7e12 | 7.189e-5 | 1.362e-4 | 1.9 | 1.67e-4 | 3.00e-4 |
| 30 | 6.6e12 | 1.480e-5 | 3.086e-5 | 2.1 | 4.24e-5 | 1.04e-4 |
| 60 | 8.8e12 | 1.096e-6 | 2.662e-6 | 2.4 | 3.94e-6 | 1.29e-5 |
| 120 | 1.15e13 | 8.323e-8 | 2.336e-7 | 2.8 | 3.47e-7 | 1.57e-6 |

| NE step | l2-optimal ν ∝ h^q | order of the minimum l2 | ν=0 l2 order | order of the minimum max |
|---|---|---|---|---|
| 16→20 | 3.88 | 4.21 | 3.99 | 4.20 |
| 20→30 | 3.61 | 3.90 | 3.66 | 3.38 |
| 30→60 | 3.59 | 3.76 | 3.53 | 3.43 |
| 60→120 | 3.61 | 3.72 | 3.51 | 3.50 |

- The heuristic gets the orders about right: 3.72–3.76 in l2 (predicted 3.75) and 3.43–3.50
  in max (predicted 3.5), against 3.51 and 3.04 at ν=0. The gain over ν=0 therefore grows with
  resolution, from 2.1× at NE30 to 2.8× at NE120.
- It gets the optimal-ν exponent wrong: ν_opt ∝ h^3.6, steady over three steps, not h^3.75.
  A sum of two power laws forces order = q; the data have the order 0.1–0.15 above q at every
  step, so the model is too crude (part of the hyperviscosity's error probably cancels part of
  the scheme's error instead of adding to it), and the agreement on the order is partly
  coincidental.
- For this test the l2-optimal ν is about 6.6e12·(30/NE)^3.6 m⁴/s. The convergence study's
  ν₃₀ = 1e13 ∝ h⁴ is close to optimal at NE120 and about 1.5× too large at NE30, consistent
  with its measured order 3.71. The E3SM default ν₃₀ = 1e15 is 100–150× the optimum for this
  Gaussian.
- The max-optimal ν is erratic (the max error is flat in ν, with several local minima); the
  minimum max error itself is well defined.
- At α=0 (NE16, 20 only) the best ν also halves the l2 error, with the same optimal ν₃₀
  (6.2e12, 5.8e12): at these resolutions most of the error is in unresolved modes, which
  hyperviscosity damps.

Cost on 4 cores: B extraction 30 s (NE30), 2 min (NE60), 9 min (NE120, 94.5 M entries, 1.5 GB);
scans with about 30 ν values 1.5 min (NE30), 10 min (NE60) and 2 h (NE120, 20 ν values).

## Open questions and suggested next steps

1. **Link the spectrum to the convergence orders: done** (section above). Remaining: a proof
   of the 3.5 / 3 orders. Suggested route: a planar model of one cube edge (two uniform SEM
   grids glued along a line with the local kink angle and stretching jump, HOMME's split-form
   weak form), Fourier-transformed along the edge into small 1-D eigenproblems per along-edge
   wavenumber. It is scale-free (eigenvalues (u/h)·f(kh, θ)), so it gives the edge-trapped
   modes, the interface truncation error and the resonant response for every h, with controls
   (θ = 0, flow across the edge). Corners (valence 3) need a separate 2-D local model.
   `mt5555/sem-advection-tests` already has planar kinked-mesh SEM codes.
2. **NE30 dense spectra and eigenvectors on a bigger machine** (e.g. a Perlmutter CPU node, 512 GB). n=48,602:
   dense L is 19 GB, `dgeev` eigenvalues-only needs about 40–60 GB and ~1e15 flops (hours).
   For the split form with the curl ψ velocity, S = M^½ L M^-½ is skew, so use a symmetric
   solver on −S² (or Hermitian iS): much faster and parallel (dsyevd, ELPA, cuSOLVER).
   Regenerate the NE30 matrices with `swtc1-opmat.sh` (minutes); they are not in git.
   NE30 eigenvectors would give the edge-trapped modes and the modal budget in the asymptotic
   regime (`opmat_proj.py --eig 30`). With the antipodal symmetry (October 2026) the skew-part
   eigenvectors fit on a 15 GB machine: four Hermitian blocks of about 12,150 for α=π/4.
3. **Why α=π/4's max norm drops to p=3: largely answered** (edge-trapped modes, section
   above). Worth checking whether the planar kinked-mesh tests had flow running parallel to a
   kink for a long time, which is what the α=π/4 path does. Possible ways to recover 4th order, to test:
   a grid that is smooth across cube edges (a conformal cubed sphere is singular only at the
   corners, whose source converged at 4.3–4.6); or dissipation applied only in a one-element
   strip along the cube edges (heuristic: rate ∝ h^{-1/2} gives l2 4, max 3.5). The second can be
   tested with the extracted matrices: exp((L − ν χBχ)T)p0 with χ the strip mask, which stays
   dissipative in the mass inner product (needs B at NE30 and NE60).
4. **Optional extensions**: np=7 (p=6) operators; the element-local map
   (`cubed_sphere_map=2`, the theta-model default, settable in the namelist); B scaled as
   dx³ vs dx⁴ combined with L to study the error-minimizing ν analytically.

## Related material

- Gallery with the convergence figures: https://claude.ai/artifact/Utx6BesN44TNqMYDdTiQbp
- Independent planar and hexagon tests: https://claude.ai/artifact/JzjnkXpN77xLTbgyRpGrsG and
  the repo `mt5555/sem-advection-tests`
- Convergence-study notes: https://claude.ai/artifact/2RPeSzKDFKivusjdLeanFD
