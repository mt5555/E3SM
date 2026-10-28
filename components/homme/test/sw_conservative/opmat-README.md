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

Build sweqx in its own directory with `-DSWEQX_PLEV=4 -DSWEQX_NP=4` (see `swtc1.sh`), then:

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

### Python (test/sw_conservative/python/)

- `opmat_io.py`: `load(prefix)` returns `L`, `B` (scipy CSR), `lat`, `lon`, `mass`, `meta`, `Ltest`.
- `opmat_sym.py`: block-diagonalizes L with the rotational symmetry of grid and flow
  (90° about z for α=0, group Z4; 180° about the α=π/4 rotation axis, which passes through
  two cube-edge midpoints, group Z2). Block eigenvalues match full dense ones to 1e-12 ω.
- `opmat_eigs.py --mats DIR [--mats30 DIR30] --out OUT [--ne 8 15]`: the full analysis
  (checks, spectra, accuracy against −i m ω, localization, hyperviscosity), cached in
  `OUT/cache`, plots and `results.txt` in `OUT`. NE15 takes about 10 minutes on 4 cores.
  Set `OPENBLAS_NUM_THREADS=1` if the cores are shared: oversubscribed OpenBLAS was about
  15× slower.

### HOMME details found along the way

- `elem%gdofP` is not contiguous (it ranges up to nelem·np²); use the unique-point ids.
- `elem%spheremp` is not the assembled mass; the assembled mass is `1/elem%rspheremp`.

## Results so far (NE8, NE15 dense; NE30 extremal eigenvalues with ARPACK)

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

## Open questions and suggested next steps

1. **Link the spectrum to the convergence orders.** The convergence loss at α=π/4 is a
   property of how well the operator represents the Gaussian's modes. Project the initial
   Gaussian onto the eigenvectors (NE8, NE15; NE20 is feasible with 15 GB) and compute the
   predicted day-12 error Σ |c_j|²·|e^{λ_j T} − e^{−i m_j ω T}|² for α=0 and α=π/4. Does it
   reproduce the measured errors and orders? Which modes carry the error, and are they
   corner/edge-localized?
2. **NE30 dense spectra on a bigger machine** (e.g. a Perlmutter CPU node, 512 GB). n=48,602:
   dense L is 19 GB, `dgeev` eigenvalues-only needs about 40–60 GB and ~1e15 flops (hours).
   For the split form with the curl ψ velocity, S = M^½ L M^-½ is skew, so use a symmetric
   solver on −S² (or Hermitian iS): much faster and parallel (dsyevd, ELPA, cuSOLVER).
   Regenerate the NE30 matrices with `swtc1-opmat.sh` (minutes); they are not in git.
3. **Why α=π/4's max norm drops to p=3** when planar non-smooth meshes stay ≥3.6: look at the
   eigenvectors and the error budget near the cube edges the α=π/4 path follows.
4. **Optional extensions**: np=7 (p=6) operators; the element-local map
   (`cubed_sphere_map=2`, the theta-model default, settable in the namelist); B scaled as
   dx³ vs dx⁴ combined with L to study the error-minimizing ν analytically.

## Related material

- Gallery with the convergence figures: https://claude.ai/artifact/Utx6BesN44TNqMYDdTiQbp
- Independent planar and hexagon tests: https://claude.ai/artifact/JzjnkXpN77xLTbgyRpGrsG and
  the repo `mt5555/sem-advection-tests`
- Convergence-study notes: https://claude.ai/artifact/2RPeSzKDFKivusjdLeanFD
