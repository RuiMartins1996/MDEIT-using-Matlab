# Handoff: Tier 3 completion, sensor-optimization payoff analysis, brute-force validation, 2D domain

Written 2026-07-24. Continues from `HANDOFF_tier3_data_space.md` (Tier 3 is
now DONE, see below) in the same folder,
`studies/optimal_sensors_bayesian_approach/`. This document summarizes
everything done in the session after that handoff, for a new Claude session
picking this up cold.

## 1. Tier 3 (data-space Woodbury reformulation) -- DONE, validated

`HANDOFF_tier3_data_space.md`'s plan was implemented in
`example_anomaly_circle.m`:
- `costgrad_theta_z`'s cost block and gradient-weight block now use the
  data-space objects `P = J*Gamma_prior`, `S`, `Ls = chol(S,'lower')`,
  `Y = S\P` instead of forming the dense `n x n` `Gamma_post`.
- `compute_posterior_variance_diag` signature changed to take
  `(Gamma_prior, Gamma_noise)` directly (no more `inv_Gamma_prior`/
  `inv_Gamma_noise` args) and uses the same data-space diag formula.
- **The d-opt cost formula (Tier 1 item "1.1") was included** (user chose
  this explicitly when asked) -- `phi_data` for `'d-opt'` is
  `sum(log(diag(Gamma_noise))) - sum(log(d)) - logdet_S`. `opt_modes` at
  the top of the script is still overridden to `{'a-opt'}` only
  (`opt_modes = {'a-opt','d-opt'}; opt_modes = {'a-opt'};`) -- the d-opt
  math works and was smoke-tested (gradient check rel err ~1e-7, cost
  decreases monotonically under fminunc) but is not the default run mode.

Validated three ways (regression against a reconstructed pre-Tier-3
baseline, FD gradient check ~1e-8, posterior-variance-sum match) -- see
git history / prior conversation for the exact numbers if needed; not
re-derived here since the code has since evolved (see section 2 below).

## 2. `example_anomaly_circle.m` has been heavily retuned since Tier 3

The user iteratively found which parameter changes make sensor
optimization actually matter (starting point: ROI posterior-variance
ratio optimized/even was ~0.997, i.e. optimization barely helped). Key
findings, in case future tuning needs the reasoning:

- **ROI size vs. `d_target` mismatch is the dominant ceiling.** If the ROI
  has far more elements than `d_target` (data-dominated modes), no sensor
  arrangement can reduce A-opt cost much -- this is a hard information
  ceiling, independent of geometry. (Original config: ROI=1608 elements,
  d_target=75 -- ~20x mismatch, ratio~0.997.)
- **Moving sensors closer to the domain helps a lot.** The Biot-Savart
  kernel falls off steeply with distance; sensors far away (`rs` large)
  see an angle-smoothed field where azimuthal position barely matters.
  Closer sensors sit in the kernel's steep region, giving the optimizer
  real leverage.
- **Fewer electrodes (stimulation patterns) helps.** With many current
  patterns, the excitation basis alone is already rich, making sensor
  placement a secondary readout choice. Fewer patterns forces more of the
  information burden onto sensor placement.
- Current committed parameters in `example_anomaly_circle.m` (as of this
  writing): `numOfElectrodesPerRing=4`, `n_sensors=4`,
  `rs=1.05*model_parameters.radius` (close to domain), `height=0.6`
  (only used for `zs=height/2`), anomaly `material.radius=radius/5` at
  `position=[3*radius/4,0,height/2]`, `d_target=10`,
  `roi_radius_factor=1.0`, `alpha_rep=1e-5`, `maxsz=radius/20`, recon mesh
  `maxsz=radius/8`, `max_iterations=30`. With this config,
  `phi_even=3.6307e-01`, ROI posterior-variance ratio (fminunc/even)
  ~0.974 (or better, see section 3 below for the true global optimum).
- **On quantifying "was optimizing worth it": use the ROI posterior-
  variance ratio** (already printed at the end of the script), not raw
  SNR -- SNR doesn't capture the Fisher-information/geometry argument the
  OED criteria are built on. Secondary metric: A-opt trace-reduction
  percentage (also printed); tertiary: Monte-Carlo-averaged reconstruction
  error in the ROI (not currently implemented, would need averaging over
  many noise draws).

## 3. `brute_force_4x4.m` -- NEW, validates fminunc against exhaustive search (3D)

Created to check whether `fminunc` (gradient-based, local) actually finds
the global A-opt minimum for the 4-electrode/4-magnetometer config in
`example_anomaly_circle.m`. Mirrors that script's exact model setup, then:

1. Runs `fminunc` from the even-spacing start.
2. Exhaustive grid search over sensor angles: since the 4 sensors are
   unordered, uses `nchoosek(1:grid_n,4)` combinations (not `grid_n^4`
   permutations). Default `grid_n=24` (15 deg steps, 10626 combos);
   `quick_test=true` uses `grid_n=10` for a fast smoke test.
3. Locally polishes the best grid point with `fminunc` (removes grid
   discretization error) -- this is the fairest "true global minimum"
   estimate.
4. Reports `phi_even`/`phi_fminunc`/`phi_brute_force`/
   `phi_brute_force_polished`, resulting angles, and ROI
   posterior-variance ratios. Saves to
   `data/brute_force_4x4_results.mat` (includes the full `phi_grid`, so a
   future session can plot the cost landscape without rerunning).

**Full-run result (grid_n=24, max_iterations=30):**
`phi_even=3.6307e-01`, `phi_fminunc=3.5523e-01`,
`phi_brute_force_polished=3.5513e-01` -- **gap only 0.027% of phi_even**:
with enough iterations, fminunc essentially finds the global optimum. Both
converge to the same physical configuration (all four sensors clustered
near azimuth 0 deg, where the anomaly is).

**Important caveat found via a coarse/few-iteration smoke test:** with
`max_iterations=3` (i.e. `quick_test=true`), fminunc gets stuck very near
the even-spacing start (a near-stationary/symmetric point) and looks like
it barely helps -- this is NOT evidence the optimizer is broken, it just
needs enough iterations to escape the symmetric saddle region. Don't judge
optimizer quality from a `quick_test` run alone.

Local functions in this script are a **self-contained copy** of
`example_anomaly_circle.m`'s local functions (tet-based Biot-Savart
quadrature, 35-point rule) -- it does not import from that file (MATLAB
scripts' local functions aren't exported).

## 4. `example_anomaly_circle_2d.m` -- NEW, 2D (flat disk) domain version

Same study as `example_anomaly_circle.m` but the conductivity mesh is a
flat 2D triangular mesh (`fmdl.nodes` is `Nx2`, `fmdl.elems` is `Nx3`), not
a 3D tetrahedral cylinder. Sensors remain 3D points at height `zs` above
the plane; Biot-Savart is integrated over the flat current sheet at z=0 --
this is EIDORS's own documented 2D MDEIT convention (see
`functions/fwd_model/compute_r_matrices.m`,
`compute_r_matrices_quadrature_2d`), not an approximation invented for
this script.

**Key technical facts established while building this (read before
touching the 2D pipeline again):**
- `create_kai_2d_model_parameters` (not `_3d`) + the repo's existing
  `is2D=true` branch in `mk_mdeit_model.m` builds a genuine flat mesh via
  `ng_mk_cyl_models(cyl_shape=[0,radius,maxsz],...)`.
- The anomaly must use `material.type='cylindrical'` (a disk) --
  `'spherical'` is **not implemented** for 2D in
  `create_model_with_material_2d` (errors out).
- **`calc_Jz`, `grad_z`, `build_ctx`, `costgrad_theta_z`,
  `compute_posterior_variance_diag` are all dimension-agnostic and were
  reused UNCHANGED from the 3D script.** This works because
  `compute_gradient_matrix.m` sets `G.Gz` to an exact zero sparse matrix
  for a 2D mesh (confirmed by reading that function) -- the z-terms in the
  shared code vanish on their own, no branching needed.
- **Only the Biot-Savart quadrature needed a 2D rewrite**: tet/35-point
  local functions (`compute_R`, `compute_dR`, `precompute_quad_geometry`,
  `compute_r_matrices_local`, `quad_rule_35`) were replaced with
  triangle/37-point-TOMS analogs (`compute_R_2d`, `compute_dR_2d`,
  `precompute_quad_geometry_2d`, `compute_r_matrices_local` (2D version),
  `quad_rule_tri37`), modeled directly on the production
  `functions/fwd_model/compute_r_matrices.m` implementation
  (`compute_r_matrices_quadrature_2d`). The key trick: quadrature points
  `xi` are `numQuadPts x 2 x numElements` (in-plane only); the sensor's
  z-height is added back explicitly as a constant third coordinate
  (`dz = sensor_z`, since the current sheet is at z=0) before computing
  the same closed-form Biot-Savart kernel and its derivative -- the
  derivative formula is literally identical to the 3D case since it only
  depends on the 3-component vector from quadrature point to sensor, not
  on the source mesh's dimensionality.
- ROI/ centroid distance calc uses `anomaly_center(1:2)` only (2D
  centroids have no z column).
- Figures: 3D ROI-sphere `surf`+`view(3)` replaced with a 2D disk outline
  (`fill`) and plain 2D axes; the posterior-std figure no longer needs a
  z-slice filter since the whole domain already IS the one 2D slice.
- Current parameters (after a user/linter edit post-creation):
  `maxsz=radius/30`, recon mesh `maxsz=radius/15` (finer than what was
  first written, `radius/20`/`radius/8`) -- otherwise same 4-electrode/
  4-sensor/`d_target=10`/`roi_radius_factor=1.0` config as the 3D script.

**Validated:** gradient check rel err ~1e-7 (same order as the 3D script),
full pipeline (optimize -> reconstruct -> posterior-variance diagnostics
-> figures) runs to completion with no errors. Smoke-tested and one
full-parameter run done (`n_elem=1686`, `phi_even=1.6076e-01` in the
`brute_force_4x4_2d.m` run below -- the standalone
`example_anomaly_circle_2d.m` run used a slightly different quick_test
config so its own printed numbers differ slightly; both are internally
consistent, don't cross-compare their raw numbers without re-running).

## 5. `brute_force_4x4_2d.m` -- NEW, 2D counterpart of section 3

Identical structure to `brute_force_4x4.m` but reuses the 2D quadrature
local functions from section 4 instead of the tet ones (again a
self-contained copy, not an import). Same model config as
`example_anomaly_circle_2d.m`'s current 4-electrode/4-sensor setup.

**Full-run result (grid_n=24, max_iterations=30):**
`phi_even=1.6076e-01`, `phi_fminunc=1.5545e-01` (fully converged --
"Local minimum found", gradient norm ~8e-6, NOT iteration-limited),
`phi_brute_force_polished=1.5350e-01` -- **gap = 1.2% of phi_even**,
notably larger than the 3D case's 0.027%.

**Important finding, worth investigating further if anyone picks up 2D
work:** unlike the 3D case, fminunc here converges (fully, not just
iteration-limited) to a **different, meaningfully worse local minimum**
than the true global one:
- `theta_fminunc` ~ `[5.9, 244.6, 345.6, 19.4]` deg -- three sensors
  clustered near 0-20 deg, one isolated near 245 deg (roughly opposite
  side).
- `theta_brute_polished` ~ `[7.2, 20.1, 340.1, 354.5]` deg -- all four
  sensors clustered together near the anomaly's azimuth (0 deg).

This means the 2D A-opt cost landscape (at least for this config) has
multiple distinct local minima that a single gradient-based run from even
spacing can fall into. **Recommendation for future work on the 2D script:
use multistart (`n_starts > 1` in `example_anomaly_circle_2d.m`, which
already has the multistart loop wired up but set to `n_starts=1`) rather
than trusting a single fminunc run.** It's not yet established whether
this is a generic 2D-vs-3D difference or specific to this mesh/parameter
combination -- worth checking multistart results and/or plotting the
saved `phi_grid` (in `data/brute_force_4x4_2d_results.mat`) before drawing
a general conclusion.

## Files in this folder after this session

- `example_anomaly_circle.m` -- modified (Tier 3 + retuned params), see
   1,  2.
- `brute_force_4x4.m` -- new, see  3.
- `example_anomaly_circle_2d.m` -- new, see  4.
- `brute_force_4x4_2d.m` -- new, see  5.
- `HANDOFF_tier3_data_space.md` -- prior handoff, now fully actioned
  (superseded by this document for anything beyond Tier 3 itself).
- Saved results: `data/example_anomaly_circle_results.mat`,
  `data/brute_force_4x4_results.mat`,
  `data/example_anomaly_circle_2d_results.mat`,
  `data/brute_force_4x4_2d_results.mat`.
- No stray log/temp files were left behind (cleaned up after each run in
  this session).

## Suggested next steps (not started)

1. Run `example_anomaly_circle_2d.m` with `n_starts=3` or more (uncomment/
   set at the top) to check whether multistart on the ORIGINAL script
   (not just brute force) reliably finds the better minimum found by
   brute force in  5.
2. If the multi-local-minima behavior in 2D turns out to be robust (not
   just this one config), consider adding multistart as the *default* in
   `example_anomaly_circle_2d.m` rather than leaving `n_starts=1`.
3. The d-opt cost path (section 1) is implemented but not exercised as the
   default `opt_modes`; if d-opt results are wanted for the paper/study,
   flip `opt_modes = {'a-opt','d-opt'}` (already present, just
   commented-out-by-override) and re-validate the printed "nats" numbers.
4. Tier 3.2 (fusing `grad_z` to avoid materializing intermediate
   `[num_stim x n_elem]` blocks) mentioned as optional/deferred in the
   original Tier 3 handoff -- still not done, likely low priority given
   current mesh sizes are small (hundreds to a few thousand elements) and
   FEM/Biot-Savart assembly, not linear algebra, dominates runtime now.
