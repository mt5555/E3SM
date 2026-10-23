# mtest1.sh: HOMME theta-l mountain-wave tests

`mtest1.sh` builds the `theta-l` HOMME dycore, runs one of the
non-hydrostatic mountain-wave test cases on a reduced-radius planet at a
chosen vertical resolution, and generates lon/height cross-section plots
of the wave field using the matplotlib-based `plotxz.py` script.

## Running

The script must be submitted with `sbatch` (it will not run the model if
invoked directly, since it needs `$SLURM_NNODES` to launch with `srun`):

```bash
cd components/homme/test/mtest
sbatch mtest1.sh
```

Since `sbatch` jobs do not inherit shell variables unless explicitly
exported, set options via `--export`, e.g.:

```bash
sbatch --export=ALL,testcase=m2,nlev=144,dcmip_rayleigh_depth=10000 mtest1.sh
sbatch --export=ALL,nu_top=4000,sponge_start=100,dcmip_rayleigh_depth=0 mtest1.sh
```

### GPU dycore (`theta-l_kokkos`)

These tests run only on CPUs. The GPU dycore (`theta-l_kokkos`) works,
but the test-case Rayleigh friction forcing is still computed in Fortran
on the CPU. With `test_with_forcing=.true.`, the full state is copied from
GPU to host and back every timestep, which makes GPU runs horribly slow.
For example, L144 m2 on 1 GPU node was about 13× slower than on 8 CPU nodes.

## Options

### `testcase` (default `m2`)

| testcase | template | test | grid | model top (analytic levels) |
|---|---|---|---|---|
| `dcmip` | `dcmip2012-nh-template.nl` | DCMIP2012 test 2-1, mountain waves without shear | ne20, rearth = a/500 | 30 km, 60 levels |
| `m1` | `m1-ne30-template.nl` | KSP2015 M1, quasi-2D mountain ridge | ne30, rearth = a/166.7 | 20 km, 40 levels |
| `m2` | `m2-ne60-template.nl` | KSP2015 M2, circular mountain, no shear | ne60, rearth = a/166.7 | 30 km, 60 levels |

### `nlev` (default: unset)

- **Unset:** use analytic, evenly spaced height levels (`vanalytic=1`)
  with the per-test level count and model top listed above.
- **Set:** use `vanalytic=0` and read vertical coordinates from
  `vcoord-<nlev>m.ascii` / `vcoord-<nlev>i.ascii` in this directory
  (available for 144, 152, 160, 176, 184, 192, 256). The model top is then
  set by the file's top interface (z = H log(1/eta)), which is much higher
  than the analytic-level top.

### `dcmip_rayleigh_depth` (default `-1`)

Controls the Rayleigh friction sponge applied by the test case forcing
(namelist variable `dcmip_rayleigh_depth` in `&ctl_nl`):

- `-1`: use the test case default sponge (dcmip and m2: starts at
  20 km; m1: starts at 10 km), ramping up to the model top.
- `0`: Rayleigh friction off.
- `>0`: sponge depth in meters below the model top (starts at
  `ztop - dcmip_rayleigh_depth`). The model aborts if the depth exceeds
  the model top.

The damping coefficient ramps as sin^2 from zero at the sponge base to
`1/tau` at the model top. `tau` is hardcoded (25 s for dcmip, 75 s for
the m tests). The model log reports the sponge base, top, depth, and the
levels it covers.

### `nmax` (default: unset)

- **Unset:** use the template's `nmax`. Every template is set to the
  test case default of 2 h: dcmip 72000 × 0.1 s, m1 4000 × 1.8 s, and
  m2 8000 × 0.9 s.
- **Set:** override the number of timesteps. The output directory and
  namelist names get a `-n<nmax>` suffix. Example: 12 h of m2 is
  `nmax=48000`. Increase the job time limit to match (e.g.
  `sbatch -N 16 -t 1:30:00 --export=ALL,nmax=48000,...`).

### `nu_top` (default `0`, disabled)

Coefficient of HOMME's top-of-model Laplacian (second-order) diffusion
sponge, the same sponge used by EAM/EAMxx. It is applied only when
`nu_top > 0` and `hypervis_subcycle_tom > 0` (set in each template).

### `sponge_start` (default `2`)

Reference pressure (hPa) where the `nu_top` Laplacian sponge starts
(namelist `tom_sponge_start`). It has no effect when `nu_top=0`.

## How output is organized

- **Build directory:** `$SCRATCH/bld-$nlev`. There is one shared build
  per `nlev` (60 or 40 when `nlev` is unset), reused across test cases and
  reruns. The job runs `make theta-l` at startup, so concurrent jobs with
  the same `nlev` can collide at link time. Build first, or stagger
  submissions. It is only reconfigured if `CMakeCache.txt` is missing.
- **Working directory:** `$SCRATCH/<testdir>` (`WDIR`), where `<testdir>`
  is `dcmip2012-2.1`, `m1`, or `m2`.
- **History output:**
  `$WDIR/moviesVL-L<nlev>-nu_top<nu_top>-<sponge_start>-rf<dcmip_rayleigh_depth>[-n<nmax>]/<outfile>`,
  where `<outfile>` is `dcmip2012_test2_11.nc`, `mtest11.nc`, or
  `mtest21.nc`. Each option combination gets its own subdirectory.
- **Namelist:** generated per run at
  `$WDIR/input-L<nlev>-nu_top<nu_top>-<sponge_start>-rf<dcmip_rayleigh_depth>[-n<nmax>].nl`
  from the test's template. The template placeholders are `@NU_TOP@`,
  `@SPONGE_START@`, `@MOVIES_DIR@`, `@VANALYTIC@`, `@VFILE_MID@`,
  `@VFILE_INT@`, and `@DCMIP_RAYLEIGH_DEPTH@`. If `nmax` is set, the
  template's `nmax` line is replaced.

`$SCRATCH` is currently `/pscratch/sd/t/taylorm/homme`.

## Plotting

After the model run completes, `mtest1.sh` calls `plotxz.py` (from
`~/codes/nclscript/contour/`) on the history file:

- `dcmip`: `Tprime` (`-c -1.2,1.2`) and `w` (`-c -3.5,3.5`)
- `m1`, `m2`: `thetaprime` (`-c -3.0,3.0`) and `w` (`-c -2.0,2.0`)

The symmetric `-c` ranges select `plotxz.py`'s diverging colormap. The
plots are written next to the history file as
`<outfile-base>.<var>_t<NNN>.png` (by default, the last output time).
See `plotxz.py --help` for options, including `-t` to select a different
time index.
