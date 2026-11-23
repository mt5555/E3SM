# DG option in sweqx (`dg_flux`)

sweqx can couple its spectral elements with discontinuous Galerkin (DG) fluxes
instead of the DSS.  The volume terms are HOMME's usual strong-form operators;
only the element coupling changes.  This directory has the run and plot
scripts for comparing CG and DG on the swtc1 transport test.

| `dg_flux` | method |
| --- | --- |
| -1 (default) | CG: DSS, bit-for-bit as before |
| 0 | DG, central fluxes (conserves q^2 / energy semi-discretely) |
| 1 | DG, upwind tracer flux (transport) or Rusanov Bernoulli flux (shallow water) |

Requirements: `integration = "runge_kutta"`, `rk_stage_user = -4` (RK4),
`nu = nu_s = 0`.  For shallow water (any test case without a prescribed wind)
the flux form is required (`theta_advect_form /= 2`).  The run aborts with a
message otherwise.

## Scripts

All runs use classical RK4 and the namelist template
`../sw_conservative/swtc1.nl`.

| file | what it does |
| --- | --- |
| `run.sh` | Builds sweqx (unless given an executable) and runs one swtc1 or swirl case. The parameters are optional arguments: `NAME NE SUBCASE TSTEP NDAYS DGFLUX FORM [STATEFREQ] [TESTCASE] [EXEC]`. `NU`, `NP`, `MOVIES`, `PLOT` and the build and output directories come from the environment. See the header for details and defaults. |
| `gallery_runs.sh` | Every run behind the two galleries, by section (`forms conv cvd nuscan stab`), with the run directory names the plot scripts expect. |
| `conv_dg.py` | swtc1-DG gallery: Gaussian convergence figures, upwind vs central flux and split vs conservation form. |
| `plots_dg.py` | swtc1-DG gallery: error and variance histories over 3 and 30 revolutions. |
| `plots_cvd.py` | CG vs DG gallery: convergence (Gaussian and all shapes), errors and variance over 3 revolutions, day-36 maps, and a table of the numbers. |

### Quick start

```bash
./run.sh                          # CG, nu = 0, NE30, dt = 750 s, 3 revolutions
./run.sh - 30 4 - 36 1            # DG upwind, same case
NU=1e13 ./run.sh cg_nu_ne30       # CG with near-optimal hyperviscosity
./run.sh - 20 4 240 12 1 2 - swirl   # swirl, DG upwind
```

Runs go to `$output/NAME` (default `~/scratch1/homme/runs`).  Each writes the
`swtc1.<shape>.<norm>.errors` files, `sweq.var` (tracer variance) and, with
`MOVIES=1`, `movies/swtc11.nc`.

### Reproducing the galleries

```bash
./gallery_runs.sh                         # all sections; the NE120 runs take ~1-1.5 h each
./gallery_runs.sh cvd conv nuscan         # only what the CG vs DG gallery needs
python3 conv_dg.py ~/scratch1/homme/runs img
python3 plots_dg.py ~/scratch1/homme/runs img
python3 plots_cvd.py ~/scratch1/homme/runs img
```

The animated GIFs come from `../sw_conservative/python/animate_swtc1.py`, e.g.

```bash
python3 ../sw_conservative/python/animate_swtc1.py ~/scratch1/homme/runs/cvd_var_dg1 bell_dg1.gif 1 symlog
```

where the level is 1 = cosine bell, 2 = Gaussian, 3 = cylinder, 4 = slotted cylinder.

Settings used by the galleries:
- **Time step:** dt = 22500 s / NE for the convergence runs, 750 s for the
  3-revolution NE30 runs (900 s in the swtc1-DG gallery). This is below every
  method's RK4 limit.
- **CG hyperviscosity:** nu = 1e13 (30/NE)^4, with `nu_s = nu`,
  `hypervis_order = 2` and `hypervis_subcycle = 1`. The `nuscan` section shows
  this is within 5% of the best value at NE30 and NE60.
- **DG:** requires nu = 0.

### Output

Use the interpolated lat/lon movie output (`PIO_INTERP` builds, the default).
Each lat/lon point is evaluated from one element's own values, so it is correct
for DG.  The native output (`shal_movie_mod`) writes each shared GLL point once
(`UniquePoints`), which drops the other elements' edge values and is not
correct for DG.  The error norms and invariants integrate each element with its
local mass `spheremp`, so they are correct for both.

## Method

Code: `src/sweqx/advance_mod.F90` (`compute_and_apply_rhs`,
`dg_edge_correction`, `dg_swe_correction`); boundary operator
`element_boundary_integral_jcontra` in `src/share/derivative_mod.F90`; halo
unpack `edgeDGVunpack_nlyr` in `src/share/edge_mod.F90`.

- Each element exchanges the edge values of `p` and of the lat-lon velocity.
  Normal flux densities a = metdet u^n are computed from both sides' velocity
  with the element's own metric, so no orientation logic is needed.
- DG tendency: `ptens = (-spheremp*div(F) - B(Fhat - F^-))/spheremp`, with B the
  element boundary integral.  `test_edge_flux` checks
  `spheremp*div - div_wk = B` per element.
- Fluxes follow Ricardo, Lee & Duru (JCP 2024; SISC 2024):
  - transport: mass (level `kmass`, or rho = 1) central `{{rho a}}`; tracers
    `qhat*{{rho a}}` with `qhat` central or upwind.  Works with the flux and the
    split form.
  - shallow water (vector invariant): central mass flux `{{Phi a}}`; Bernoulli
    edge term `Ghat = {{G}} + alpha [[F]].n`, alpha = 0 or max(c/Phi)/2 with
    c = |u| + sqrt(Phi); absolute vorticity with the tangential `{{u}}` edge term.

`SWEQX_TC1_VELOCITY=curlpsi` gives swtc1 a discretely divergence-free wind
(curl of the stream function), used for exact variance checks.

## Time step

DG needs a smaller RK4 time step than CG.  Largest stable dt, swtc1 sub_case 4,
NE30, np=4, nu=0 (DG: 12-day runs, bisection, `gallery_runs.sh stab`; CG: 3
revolutions):

| method | dt limit (s) | ratio to CG | 1D linear analysis, ratio to CG |
| --- | --- | --- | --- |
| CG | 2329 | 1 | 1 |
| DG upwind | ~1600 (1650 grows) | 0.70 | 0.56 |
| DG central | ~1370 | 0.59 | 0.60 |

The 1D analysis is periodic constant-velocity advection, np=4: the DG
spectral radius is larger because the edge correction is divided by the end
point GLL weight (1/w_end = 6).  In 2D the upwind flux gains ~20% over central,
which linear 1D theory does not show; most likely it damps the slightly
non-skew modes of the variable-wind central operator.

Other runs: swirl at NE20 with DG at 240 s; shallow water TC2 at NE15 at
300 s (dt*NE = 4500).

## Results (np = 4, RK4, nu = 0)

Gaussian relative l2 error at day 12 (swtc1, split form, dt*NE = 22500):

| case | NE | CG | DG upwind | DG central |
| --- | --- | --- | --- | --- |
| sub_case 4 | 15 | 4.24e-4 | 1.52e-4 | 6.36e-4 |
| sub_case 4 | 30 | 3.07e-5 | 8.95e-6 | 7.86e-5 |
| sub_case 4 | 60 | 2.66e-6 | 5.58e-7 | 9.81e-6 |
| sub_case 4 | 120 | 2.33e-7 | 3.49e-8 | 1.23e-6 |
| sub_case 1 | 30 | 2.54e-5 | 1.56e-5 | 7.89e-5 |

Orders NE60-120 (sub_case 4): DG upwind 4.0 (l2 and max), DG central 3.0
(l2; max 2.6), CG 3.5 (l2; max 3.0).  CG with near-optimal hyperviscosity
(nu = 1e13 (30/NE)^4) has l2 8.12e-8 at NE120, order 3.7.

Shallow water TC2, height l2 at day 5 (dt*NE = 4500):

| NE | CG | DG Rusanov | DG central |
| --- | --- | --- | --- |
| 8 | 9.70e-6 | 3.84e-5 | 2.29e-4 |
| 15 | 6.53e-7 | 4.50e-6 | 3.35e-5 |
| 30 | 4.76e-8 | 4.29e-7 | 4.23e-6 |

Conservation checks:

- Mass: conserved to ~1e-13 (transport and shallow water).
- q^2, swtc1 with curl-psi wind, central, split form: change 1.4e-7, 4.5e-9,
  1.4e-10, 4.9e-12 at dt = 1800, 900, 450, 225 (NE15): proportional to dt^5,
  the RK4 truncation, so the semi-discrete scheme conserves it.
- rho q^2, swirl with density level, central, split form: change falls 31x when
  dt is halved (dt^5).
- q = const (swirl, Cly with the chemistry forcing off): preserved to 1e-14 by
  both DG fluxes.  With the toy chemistry on, its clipping of negative values
  adds ~1e-9 for schemes that oscillate (CG, central DG).
- Energy, TC6, central: change -4.4e-12 at dt = 150 and -3.1e-13 at dt = 75
  (NE15); the Rusanov flux dissipates.
- TC5 and TC6 run 15 days stably with the Rusanov flux and no hyperviscosity.

Not yet checked: absolute vorticity (the PV diagnostic uses the element-local
vorticity without the DG edge term), discrete geostrophic balance.
