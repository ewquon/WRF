# Mesoscale-to-Microscale Coupling with WRF

This document is intended to provide some background for various mesoscale-to-microscale coupling (MMC) approaches.

## Variable definitions

These variable names and descriptions were obtained from [`Registry/Registry.EM_COMMON`](https://github.com/a2e-mmc/WRF/blob/master/Registry/Registry.EM_COMMON) whenever available.

### Pressure/geopotential variables

| variable name | symbol            | description               | units  |
|:-------------:|:-----------------:|---------------------------|--------|
| `p`           | $p'$              | perturbation pressure     | Pa     |
| `pb`          | $\bar{p}$         | base state pressure       | Pa     |
|               | $p_d$             | dry hydrostatic pressure  |        |
|               | $p_0$             | reference sea-level pressure | Pa  |
| `ph`          | $\phi'$           | perturbation geopotential | m2 s-2 |
| `phb`         | $\overline{\phi}$ | base-state geopotential   | m2 s-2 |
| `php`         | $\phi$            | geopotential, $\phi=gz$   | m2 s-2 |

Note that the vertical coordinate (*geopotential* height) is recovered by (`PH`+`PHB`)/9.81, i.e., $z = \frac{1}{g}(\overline{\phi} + \phi')$.

From [`dyn_em/start_em.F`](https://github.com/a2e-mmc/WRF/blob/master/dyn_em/start_em.F):

> Base state pressure is a function of eta level and terrain, only, plus the hand full of constants: `p00` (sea level pressure, Pa), `t00` (sea level temperature, K), `A` (temperature difference, from 1000 mb to 300 mb, K), `tiso` (isothermal temperature at tropopause/lower stratosphere), `p_strat` (pressure at top of isothermal layer), `A_strat` (lapse rate in stratosphere above isothermal layer)


### Other common variables

| variable name | symbol | description | units |
|:-------------:|:------:|-------------|-------|
|               | $\eta$   | terrain-following hybrid sigma-pressure vertical coordinate | - |
|               | $\mu_d$  | vertical coordinate metric
| `mu`          |          | perturbation dry air mass in column | Pa |
| `t`           | $\theta$ or $\theta_m$ | potential temperature, depending on `use_theta_m` | K |
| `al`          | $\alpha'$ | inverse perturbation density of air | m3 kg-1 |
| `alt`         | $\alpha$ | inverse density of air | m3 kg-1 |
| `alb`         | $\bar\alpha$ | inverse base density of air | m3 kg-1 |

> With the WRF v4.0 release, the default behavior is to include the moist potential temperature option and to include the hybrid vertical coordinate. 
> *(ARW Version 4 Modeling System User's Guide, January 2019)*

Recall: $\theta_m = \theta (1 + \frac{461.6}{287}q_v)$

### Variables for Eulerian mass coordinate dynamics

| variable name | description | units |
|:-------------:|-------------|-------|
| `ru_tend`     |  | ? |
| `rv_tend`     |  | ? |
| `mut`         | total column dry air mass | ? |
| `muu`         | `mut` at u points | ? |
| `muv`         | `mut` at v points | ? |


### Tendency variables

The tendency variables described below are used for "internal coupling" MMC, based on code developed by Javier Sanz Rodrigo and adapted by Dries Allaerts for [WRF-to-SOWFA coupling](https://github.com/a2e-mmc/assessment/blob/study/coupling_comparison/studies/coupling_comparison/preprocessing/internal/wrf_to_sowfa.ipynb).

| field name     | definition | description | units |
|:--------------:|------------|------------------------------|---|
| `ru_tend_adv`  |            | u-advection tendency         | ? |
| `ru_tend_pgf`  | $-\mu\dfrac{m_x}{m_y}\dfrac{1}{\rho}\dfrac{\partial p}{\partial x}$ | u-pressure gradient tendency | ? |
| `ru_tend_cor`  |            | u-coriolis tendency          | ? |
| `ru_tend_phys` |            | total u physics tendency     | ? |
| `rv_tend_adv`  |            | v-advection tendency         | ? |
| `rv_tend_pgf`  | $-\mu\dfrac{m_y}{m_x}\dfrac{1}{\rho}\dfrac{\partial p}{\partial y}$ | v-pressure gradient tendency | ? |
| `rv_tend_cor`  |            | v-coriolis tendency          | ? |
| `rv_tend_phys` |            | total v physics tendency     | ? |
| `t_tend_adv`   |            | t-advection tendency         | ? |


### What is provided by auxiliary output


### What is provided by `tslist` output

These are useful as virtual met towers. Time series output of surface/near-surface quantities are also provided, detailed in [README.tslist](https://github.com/a2e-mmc/WRF/blob/master/run/README.tslist). TS and profile data are written out in [wrf_timeseries.F](https://github.com/a2e-mmc/WRF/blob/master/share/wrf_timeseries.F):

| variable name    | symbol            | description               | units |
|:----------------:|:-----------------:|---------------------------|-------|
| `ts_u_profile`   | $u$               | u wind component, earth-relative | m s-1 |
| `ts_v_profile`   | $v$               | v wind component, earth-relative | m s-1 |
| `ts_w_profile`   | $w$               | w wind component          | m s-1 |
| `ts_th_profile`  | $\theta$          | dry potential temperature, regardless of `use_theta_m` | K |
| `ts_gph_profile` | $z$ (*not* $\phi'$) | geopotential height ($=\phi/g$) | m |
| `ts_qv_profile`  | $q_v$             | water vapor mixing ratio  | kg kg-1 |
| `ts_p_profile`   | $p$               | full pressure ($\bar{p} + p'$) | Pa |


## Code calculations

### Diagnostic quantities

In [`dyn_em/module_big_step_utilities_em.F`](https://github.com/a2e-mmc/WRF/blob/master/dyn_em/module_big_step_utilities_em.F):

- subroutine `calculate_full` calculates literally `mut` = `mub` + `mu`, where `mub` and `mu` are the the perturbation and base-state dry air mass in column, respectively
- subrouinte `calc_mu_uv` calculates `muu` and `muv` by summing `mub` and `mu` (previously calculated `mut` is not used) and then destaggering

### Tendency terms

The calculations described here are based on code in [`dyn_em/module_em.F`](https://github.com/a2e-mmc/WRF/blob/master/dyn_em/module_em.F), starting from the subroutine `rk_tendency`. Note that `r*_tend` is the cumulative sum of all tendencies:

- subroutine `advect_*` or `advect_weno_*` gives `r*_tend_adv` (defined in [`dyn_em/module_advect_em.F`](https://github.com/a2e-mmc/WRF/blob/master/dyn_em/module_advect_em.F))

- subroutine `advect_scalar` or `advect_scalar_weno` gives `t_tend_adv` (defined in [`dyn_em/module_advect_em.F`](https://github.com/a2e-mmc/WRF/blob/master/dyn_em/module_advect_em.F)))

- subroutine `horizontal_pressure_gradient` gives `r*_tend_pgf` (defined in [`dyn_em/module_big_step_utilities_em.F`](https://github.com/a2e-mmc/WRF/blob/master/dyn_em/module_big_step_utilities_em.F))

	> calculates the horizontal pressure gradient terms for the large-timestep tendency in the horizontal momentum equations (u,v)
	
- subroutine `coriolis` gives `r*_tend_cor` (defined in [`dyn_em/module_big_step_utilities_em.F`](https://github.com/a2e-mmc/WRF/blob/master/dyn_em/module_big_step_utilities_em.F))

	> calculates the large timestep tendency terms in the u, v, and w momentum equations arise from the coriolis force

- subroutine `curvature` gives `r*_tend_cor` (defined in [`dyn_em/module_big_step_utilities_em.F`](https://github.com/a2e-mmc/WRF/blob/master/dyn_em/module_big_step_utilities_em.F))

	> calculates the large timestep tendency terms in the u, v, and w momentum equations arise from the curvature terms

### Horizontal pressure gradient

From [`dyn_em/module_big_step_utilities_em.F`](https://github.com/a2e-mmc/WRF/blob/master/dyn_em/module_big_step_utilities_em.F):
>  mu dp/dx = mu alpha partial dp'/dx + (nu mu partial dmubar/dx) alpha'
>           + mu partial dphi'/dx + (partial dphi/dx)*(partial dp'/dnu - mu')

>  mu dp/dy = mu alpha partial dp'/dy + (nu mu partial dmubar/dy) alpha'
>           + mu partial dphi'/dy + (partial dphi/dy)*(partial dp'/dnu - mu')
	
E.g., the west-east gradient terms:

1. $\mu \alpha \dfrac{\partial p'}{\partial x}$
2. $\nu \mu \dfrac{\partial\bar{\mu}}{\partial x} \alpha'$
3. $\mu \dfrac{\partial\phi'}{\partial x}$
4. $\dfrac{\partial\phi}{\partial y} \left(\dfrac{\partial p'}{\partial \nu} - \mu'\right)$

From the code:

-  the first three terms are essentially $\dfrac{m_x}{m_y}\mu$ times:
	- average finite difference of $\phi'$, averaged between `k` and `k+1` (term 3)
	- average $\alpha$ (between `i` and `i-1`) times the finite difference of $p'$ (term 1)
	- average $\alpha'$ (between `i` and `i-1`) times the finite difference of $\bar{p}$ (term 2?)
- the fourth term (for non-hydrostatic solver) is $\dfrac{m_x}{m_y}\mu$ times the finite difference of $\phi$ times ...

Notes:

- `cq_`: moisture coefficients, `cqv` is for `v`; calculated with `call_cq` from `rk_step_prep` in [`dyn_em/module_em.F`](https://github.com/a2e-mmc/WRF/blob/master/dyn_em/module_em.F)
- `msf__`: `msfvy` is the "map scale factor" for `v` in the `y` direction, etc
- `rdy` is $1/\Delta y$
- `muu`, `muv` correspond to $\mu$ in `x` and `y`
- **unclear why the actual $\mu$ is `c1h(k)*muu(i,j) + c2h(k)` and `c1h(k)*muv(i,j) + c2h(k)`**
