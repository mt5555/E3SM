# Mixing Ratios and Dry Mass

This page describes how EAMxx defines the mass of a layer, which mixing ratios
the tracers hold, and which functions to use to convert between them and to
compute $R^*$, virtual temperature and layer thickness.

## Layer masses

- `pseudo_density` ($\delta p$, in Pa) is the total mass of a layer per unit
  area, times $g$. It is the same as HOMME's `dp3d`. It includes dry air and
  water vapor and, when [water loading](#water-loading) is on, the condensates.
- `pseudo_density_dry` ($\delta p_{dry}$) is the dry-air mass of the layer. It
  is computed in HommeDynamics (`Homme::ElementOps::compute_dp_dry`) as

$$
\delta p_{dry} = \delta p \, (1 - q_v - q_{cond}),
$$

where $q_{cond}$ is the sum of the condensate mixing ratios, and is zero when
water loading is off.

## Wet and dry mixing ratios

The tracers in the field manager (`qv`, `qc`, `qr`, `qi`, ...) are **wet**
mixing ratios: mass per unit total mass $\delta p$. This matches HOMME, where
`Q = Qdp/dp3d`.

Some processes, such as P3 and MAM4xx, work with **dry** mixing ratios (mass per
unit dry-air mass $\delta p_{dry}$) internally. They convert the tracers on
entry and convert back before returning them. The Confluence page
[Wet or Dry Mixing Ratios in EAMxx Processes](https://e3sm.atlassian.net/wiki/spaces/DOC/pages/3020980232/Wet+or+Dry+Mixing+Ratios+in+EAMxx+Processes)
lists the convention each process uses.

Use the `PhysicsFunctions` conversions based on the layer masses:

| Function | Formula |
|---|---|
| `calculate_drymmr_from_wetmmr_dp_based(wetmmr, pseudo_density, pseudo_density_dry)` | $q^{dry} = q^{wet} \, \delta p / \delta p_{dry}$ |
| `calculate_wetmmr_from_drymmr_dp_based(drymmr, pseudo_density, pseudo_density_dry)` | $q^{wet} = q^{dry} \, \delta p_{dry} / \delta p$ |

These are exact with or without condensates. The older pair
`calculate_drymmr_from_wetmmr(wetmmr, qv_wet)` and
`calculate_wetmmr_from_drymmr(drymmr, qv_dry)` divides or multiplies by
$1 - q_v$ or $1 + q_v^{dry}$. That assumes the layer holds only dry air and
vapor, so it is exact only when water loading is off.

## Dry-air mass fraction and $R^*$

The dry-air mass fraction of a layer is

$$
q_{dry} = \frac{\delta p_{dry}}{\delta p} = 1 - q_v - q_{cond}.
$$

The gas constant of moist air, $R^*$, involves only dry air and vapor:

$$
R^* = R_d \, q_{dry} + R_v \, q_v ,
$$

where $q_{dry}$ and $q_v$ are both per unit total mass, so $q_v$ is the **wet**
vapor mixing ratio. The condensates have no partial pressure and enter only
through $q_{dry}$ (the zero volume condensates approximation). The virtual
temperature and the layer thickness follow from $R^*$:

$$
T_v = T \, \frac{R^*}{R_d} = T \left( q_{dry} + \frac{R_v}{R_d} \, q_v \right),
\qquad
\delta z = \frac{R_d}{g} \, \frac{\delta p \, T_v}{p}.
$$

These functions take $(q_{dry}, q_v)$:

| Where | Function |
|---|---|
| HOMME (C++) | `ElementOps::compute_q_dry(qv, qcond)`, `ElementOps::get_q_dry(...)`, `ElementOps::get_R_star(..., q_dry, qv, R)` |
| HOMME (Fortran) | `element_ops`: `get_q_dry(q_dry, Q)`, `get_R_star(R_star, q_dry, qv)` |
| EAMxx `PhysicsFunctions` | `calculate_virtual_temperature(T, q_dry, qv)`, `calculate_temperature_from_virtual_temperature(T_virtual, q_dry, qv)`, `calculate_dz(pseudo_density, p_mid, T_mid, q_dry, qv, water_loading)` |

In physics, compute $q_{dry}$ as `pseudo_density_dry/pseudo_density`.

!!! warning "Pass the wet $q_v$"
    $q_{dry}$ and $q_v$ must be per unit total mass. Passing a dry vapor mixing
    ratio, or a dry mass as `pseudo_density`, gives the wrong $R^*$ and layer
    thickness. See
    [E3SM-Project/E3SM#8939](https://github.com/E3SM-Project/E3SM/issues/8939).

The two-argument `calculate_virtual_temperature(T, qv)` and the four-argument
`calculate_dz(pseudo_density, p_mid, T_mid, qv)` compute
$R^* = R_d + (R_v - R_d) \, q_v$. That equals the formula above only when
there are no condensates ($q_{dry} = 1 - q_v$). With water loading off, the
$(q_{dry}, q_v)$ functions use this formula, so results are bit-for-bit with
earlier versions. Those branches are marked `BFB-only` in the code.

The `VirtualTemperature`, `AtmosphereDensity`, `VerticalLayer` (z, geopotential,
height, dz), `VertContract` (dz weighting) and `AeroComCld` diagnostics always
use $q_{dry}$ from `pseudo_density_dry`.

## Water loading

Water loading makes the condensates part of the layer mass. It is controlled
by two `driver_options` parameters:

- `water_loading` (default `false`)
- `water_loading_condensates` (default `[qc, qr, qi]`): the tracers counted as
  condensates. `qm` is not in the list because it is the rime part of `qi`.

When water loading is on:

- $\delta p$ is total mass, including the condensates, and `pseudo_density_dry`
  and $q_{dry}$ exclude them.
- When HOMME applies the physics tendencies, the dry-air mass stays fixed, so
  changes in vapor and condensates (such as precipitation) change $\delta p$
  and the surface pressure.
- Physics processes never sum the condensates themselves. They get the effect
  of water loading through `pseudo_density_dry`.

When water loading is off, results are bit-for-bit with earlier versions.
