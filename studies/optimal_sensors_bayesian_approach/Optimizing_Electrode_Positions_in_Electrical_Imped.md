# Handoff: free (theta, z) sensor optimization at fixed radius (MDEIT)

Session summary written 2026-07-06 (Claude Sonnet 5). Goal: hand this to a new session to run
the next experiment in the Bayesian sensor-position-optimization study. Everything below is in
`studies/optimal_sensors_bayesian_approach/` unless stated otherwise.

## Context (read this first, don't re-derive it)

A previous session built and tuned `example_anomaly_circle.m`: it optimizes MDEIT magnetometer
**azimuths** `theta_m` on a circle of **fixed radius** `r_s = 1.5*tank_radius` and **fixed height**
`z_s = tank_height/2`, for a conductive spherical anomaly off-center at half height. Full
background, theory, notation, code map and settled parameters are in:
- `HANDOFF_sensor_optimization_example.md` — narrative log of that session (tuning sweep, final
  numbers, bug fixes applied to `main.m`). Read the "UPDATE (2026-07-06 session...)" block at the
  bottom for the final settled configuration and results.
- `sensor_position_optimization_theory.tex` / `.pdf` — the theory companion. **New**:
  Section "Case study: sensors confined to a circle at fixed height" documents the exact setup,
  results and diagnosis from that session. Read this section before starting; it explains *why*
  the A-optimal gain was modest (4.6% trace reduction) while D-optimal gain was substantial (17.6
  nats), and ends with the open question this handoff is about.

### Settled results to date (theta free only, r and z fixed)
With `prior_std_background = 0.005*background_conductivity`, `prior_std_roi =
1.00*background_conductivity`, `roi_radius_factor = 1.0`, `d_target = 100`, `n_sensors = 8`:
- **A-optimality** (`trace(H^-1)`): even-spacing phi = 0.28791, optimized phi = 0.27476 →
  **4.6% trace reduction**.
- **D-optimality** (`-logdet(H)`): even-spacing phi = -55552.7, optimized phi = -55570.3 →
  **17.57 nats (25.35 bits)** extra information.
- Diagnosis (see theory .tex): the A-opt gain is small because A-optimality (an *averaged*
  posterior variance) is comparatively insensitive to how sensors are distributed on a rotationally
  near-symmetric ring, whereas D-optimality (log-det, which penalizes redundancy/correlation
  between sensors) benefits much more from clustering. This was confirmed by showing that (a)
  raising the SNR budget (`d_target` 60→100→150) does not change the ~4-6% A-opt figure, and (b)
  moving the anomaly closer to the tank wall barely moves A-opt while roughly doubling D-opt.

## The new task

**Question**: does letting each sensor's height `z_m` vary *in addition to* its azimuth
`theta_m` (radius still fixed at `r_s = 1.5*tank_radius`) produce a *more considerable* A-optimal
(and/or D-optimal) cost reduction than the theta-only design above?

This directly tests the closing conjecture of the theory .tex case-study section: that the modest
A-opt gain is a symmetry artifact of the overly restrictive fixed-radius-fixed-height search domain,
and that giving the optimizer one more degree of freedom per sensor (height) might break enough of
that symmetry to matter. If it does not help much either, that is itself an interesting (and
reportable) result — the task explicitly says to understand why if the reduction is still not
"considerable," the same standard as before.

## What needs to change in `example_anomaly_circle.m`

Copy `example_anomaly_circle.m` to a new script, e.g. `example_anomaly_circle_thetaz.m`, and make
the following changes. Keep everything else identical (anomaly, prior, noise calibration,
`d_target = 100`, `prior_std_background = 0.005*background_conductivity`, 3-axis MDEIT, `fminunc`
quasi-Newton with analytic gradient) so the results are directly comparable to the theta-only
baseline above.

### 1. Parametrization and chain rule

Currently (`theta_to_locations`, `costgrad_theta_3axis`): sensor `m` sits at
`(r_s cos(theta_m), r_s sin(theta_m), z_s)`, one free scalar `theta_m` per sensor, and the chain
rule collapses to
```
dphi/dtheta_m = -r_s*sin(theta_m)*dphidp(1,m) + r_s*cos(theta_m)*dphidp(2,m)
```
where `dphidp(p,m)`, `p=1,2,3`, are the full Cartesian partials already computed inside
`costgrad_theta_3axis` (from `compute_dJxyz_xyz`) — **the z-partial `dphidp(3,m)` is already
computed today and simply discarded** because z was fixed. That is the key simplification: most of
the machinery needed for the new task already exists in the script; only the parametrization and
what is done with `dphidp(3,:)` needs to change.

Two free coordinates per sensor, `(theta_m, z_m)`, both at fixed `r_s`:
- `x_m = r_s*cos(theta_m)`, `y_m = r_s*sin(theta_m)`, `z_m = z_m` (unchanged, i.e. z **is** already
  a Cartesian coordinate, so no transform is needed for it):
  ```
  dphi/dtheta_m = -r_s*sin(theta_m)*dphidp(1,m) + r_s*cos(theta_m)*dphidp(2,m)   % unchanged
  dphi/dz_m     = dphidp(3,m)                                                    % new, direct
  ```
- Design vector doubles from `M` to `2M`: `q = [theta(1..M); z(1..M)]`.

### 2. Bounding z (open decision — pick one)

`theta_m` is inherently unconstrained (periodic, `cos`/`sin` absorb any real value), which is why
`fminunc` works directly on it. `z_m` is not periodic: an unconstrained optimizer could push
sensors far above/below the tank, which is unlikely to be physically meaningful and may create
flat/badly-scaled directions. Two options:

- **(a) Sigmoid reparametrization (recommended — keeps `fminunc`, one analytic gradient, no
  solver switch).** `main.m` already has exactly this pattern for a "cylinder shell" domain
  (fixed radius, free theta, bounded z) at lines ~279-309:
  `map_x_to_q_cyl_shell` / `map_q_to_x_cyl_shell`, reproduced here for convenience:
  ```matlab
  function x = map_q_to_x_cyl_shell(q,r0,zmin,zmax)
  sigmoid = @(y) 1./(1+exp(-y));
  n_sensors = numel(q)/2;
  thetam = q(1:n_sensors);
  etam   = q(n_sensors+1:2*n_sensors);
  zm = zmin + (zmax-zmin).*sigmoid(etam);
  x = [r0.*cos(thetam(:)); r0.*sin(thetam(:)); zm(:)];
  end
  ```
  Do **not** just import this function as-is from `main.m` without re-deriving/re-checking it —
  several unrelated bugs were already found and fixed in `main.m` this session (see
  `HANDOFF_sensor_optimization_example.md` section 3), so treat any `main.m` code as a *reference
  for the idea*, not as pre-verified. Re-implement it as a local function in the new script and
  validate with the finite-difference gradient check (see below) before trusting it.
  With this map, the extra chain-rule factor for the free variable `eta_m` is
  `dz_m/deta_m = (zmax-zmin)*sigmoid(eta_m)*(1-sigmoid(eta_m))`, so
  `dphi/deta_m = dphidp(3,m) * dz_m/deta_m`.
  Suggested bounds: `zmin = 0`, `zmax = tank_height`, so that `eta_m = 0` maps to
  `z_m = tank_height/2` — i.e. the theta-only baseline's `z_s` — making `q0 = [theta_even;
  zeros(M,1)]` a natural, identical-to-baseline starting point.
- **(b) Switch to `fmincon` with bound constraints on z (theta unconstrained).** Simpler
  conceptually (no sigmoid derivative), but changes the optimizer/options and loses direct
  comparability with the `fminunc` convergence behavior used throughout the rest of the study.
  Only do this if (a) turns out to behave badly (e.g. sigmoid saturation stalling the optimizer).

### 3. Everything else that needs updating

- `assign_sensor_locations` / `theta_to_locations` → generalize to build
  `[r_s*cos(theta), r_s*sin(theta), z]` from the unpacked `q`.
- `costgrad_theta_3axis` → `costgrad_thetaz_3axis(img, q, rs, inv_Gamma_prior, inv_Gamma_noise, A,
  opt_mode, alpha_rep)`: unpack `theta = q(1:M)`, `eta = q(M+1:2*M)` (or `z` directly if not using
  the sigmoid map), assemble sensor locations, compute cost/gradient exactly as today, then return
  `dphi = [dtheta(:); deta(:)]` per the chain rule above.
- `compute_posterior_variance_diag` → also needs the 2-vector `q` instead of `theta` alone.
- **Finite-difference gradient check**: reuse `check_gradient_fd` unchanged (it is agnostic to the
  parametrization) but call it on the new `2*M`-length `q`, checking a few components from *both*
  halves (some `theta` indices and some `z`/`eta` indices) — do not skip the `z`/`eta` half, since
  that is exactly the newly-derived part of the gradient.
- Figures: the existing polar plot (azimuths only) no longer shows the full picture; add a
  side view or 3D scatter of `(theta_m, z_m)` for even/a-opt/d-opt so height migration is visible.
- Save results to a new file, e.g. `data/example_anomaly_circle_thetaz_results.mat`, so the
  theta-only baseline results are not overwritten.

## Suggested plan for the next session

1. Implement the changes above in a new script; keep all prior/noise/anomaly settings identical to
   the settled theta-only configuration (`bg_std=0.005`, `d_target=100`, anomaly at `R/2`).
2. `quick_test`-style smoke test first: 3 `fminunc` iterations + finite-difference gradient check
   on a handful of `theta` **and** `z`/`eta` components (target: rel. err. 1e-6 to 1e-8, matching
   the standard already established for this study). Do not proceed to a full run until this
   passes — every gradient bug found so far in this study (both in `main.m` and transiently while
   writing the theta-only example) was caught exactly this way.
3. Full run (40 iterations, both criteria, `n_starts=1` first). Compare the converged A-opt trace
   reduction and D-opt information gain against the theta-only baseline (4.6% / 17.57 nats).
4. If the A-opt reduction is now considerable (say >=15%): report it, and check qualitatively
   whether sensors migrate toward the anomaly in *height* as well as azimuth (expected: sensors
   should also drift toward `z = H/2`... but the anomaly is already at `z=H/2`, so the more
   interesting question is whether they spread in z to reduce redundancy, or converge tightly in z
   near the anomaly's height — look at the new diagnostic figure).
5. If the A-opt reduction is still modest: this would argue *against* the "restrictive domain"
   hypothesis from the theory .tex and *for* the "A-optimality genuinely tolerates replication"
   explanation being closer to the whole story (i.e., freeing z doesn't help because the
   redundancy A-optimality tolerates isn't really about the 1-D circular symmetry, it's more
   fundamental to the criterion itself). Either outcome should be written back into the theory
   .tex case-study discussion as a follow-up remark/subsection, using the same notation.
6. Given the design space doubles in dimension (`2M = 16` free variables), run a multistart check
   (`n_starts = 3` to 5) — more degrees of freedom means more opportunity for poor local minima
   than the 8-dimensional theta-only problem had (which turned out to have none).
7. Update `HANDOFF_sensor_optimization_example.md` (or a new dated section) and the theory `.tex`
   with the final results, following the same style as the existing "UPDATE" block and case-study
   section — same notation, same standard of reporting a tuning sweep, converged numbers, and an
   honest diagnosis if the result isn't clearly positive.

## Practical notes
- Always run with **R2024a** (`"C:\Program Files\MATLAB\R2024a\bin\matlab.exe"`) — Optimization
  Toolbox (`fminunc`) is only installed there; R2025a/b lack it.
- The mesh/model is cached under the repo `models/` folder; changing `model_parameters` re-meshes
  via Netgen automatically (fast, a few seconds, for this model size).
- `pdflatex` is available at
  `C:\Users\RuiMartins\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe` (2 passes needed
  for cross-references); clean up `.aux`/`.log`/`.out` after compiling, only the `.tex`/`.pdf` are
  kept in the folder.
- Model size/timing reference (theta-only case): 3425 elements, 817 nodes, 16 stims (1 ring x 16
  electrodes at half height), 8 sensors, ~2-4 s per `fminunc` iteration with 6 parallel workers.
  Expect roughly similar or somewhat slower per-iteration cost for the 2M-dimensional problem
  (the Jacobian/adjoint cost per iterate is dominated by PCG solves over sensors and is largely
  unchanged; only the chain-rule assembly and `fminunc`'s internal linear algebra scale with the
  number of free variables).
