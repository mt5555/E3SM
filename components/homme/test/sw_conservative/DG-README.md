# DG option in sweqx (`dg_flux`)

sweqx can couple its spectral elements with discontinuous Galerkin (DG) fluxes
instead of the DSS.  The volume terms are HOMME's usual strong-form operators;
only the element coupling changes.

| `dg_flux` | method |
| --- | --- |
| -1 (default) | CG: DSS, bit-for-bit as before |
| 0 | DG, central fluxes (conserves q^2 / energy semi-discretely) |
| 1 | DG, upwind tracer flux (transport) or Rusanov Bernoulli flux (shallow water) |

Requirements: `integration = "runge_kutta"`, `rk_stage_user = -4` (RK4),
`nu = nu_s = 0`.  For shallow water (any test case without a prescribed wind)
the flux form is required (`theta_advect_form /= 2`).  The run aborts with a
message otherwise.  The run scripts take `dg_flux` from the environment, e.g.
`dg_flux=1 ./swtc1.sh` (also set `nu=0` and a smaller `tstep`).

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

DG needs a smaller RK4 time step than CG.  swtc1 at NE30 (np=4): CG is stable
at 1800 s; DG blows up at 1800 s and runs at 900 s.  Swirl at NE20: DG at 240 s.
Shallow water TC2 at NE15: 300 s (dt*NE = 4500).

## Results (np = 4, RK4, nu = 0)

Gaussian relative l2 error at day 12 (swtc1, split form, dt*NE = 22500):

| case | NE | CG | DG upwind | DG central |
| --- | --- | --- | --- | --- |
| sub_case 4 | 15 | see conv runs | | |
| sub_case 4 | 30 | 3.07e-5 | 8.95e-6 | 7.86e-5 |
| sub_case 4 | 60 | 2.66e-6 | | |
| sub_case 1 | 30 | 2.54e-5 | 1.56e-5 | 7.89e-5 |

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
