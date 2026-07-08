# HANDOFF: Bayesian sensor-position optimization on the human thorax (MDEIT)

Plan written 2026-07-06 (Claude Fable 5, planning session). This file is the working
handoff for the HumanThorax phase of the sensor-optimization study. **Protocol: each
agent session works on exactly ONE task below, then appends an UPDATE section to this
file (same style as `../HANDOFF_sensor_optimization_example.md`) so the next agent can
continue.** The next agent to pick this up starts with **Task A**.

All paths are relative to `studies/optimal_sensors_bayesian_approach/` unless stated.

---

## 0. Required reading before touching anything

1. `../HANDOFF_sensor_optimization_example.md` — the finished cylinder-geometry study.
   Contains the methodology standards (noise calibration, diagnostics, verification
   protocol, result tables) that THIS study must follow, plus the list of bugs that
   were found and fixed. Read it fully.
2. `../sensor_position_optimization_theory.tex` — the theory write-up. Especially
   §"Case study: sensors confined to a circle at fixed height" (`\label{sec:case-study}`)
   and its subsections (ROI prior, noise calibration by target data-dominated rank,
   discussion of why A-opt gains stay modest on a symmetric ring). **You will extend
   this document** (see per-task instructions). Compile with MiKTeX `pdflatex`, 2 passes.
3. `../example_anomaly_circle.m` and `../example_anomaly_circle_thetaz.m` — the
   validated cylinder scripts (analytic gradients FD-checked to rel. err ~1e-7).
4. `HumanThorax/human_thorax_example.m` — the current thorax script (this folder).
5. `functions/optimize_sensor_configuration.m` (repo-level `functions/`) — the shared
   optimizer used by the thorax script (cost + gradient for a-opt/d-opt, 1/3-axis,
   generalized-coordinate maps).

## 1. Goal of this phase

The cylinder study is done. We now repeat the Bayesian OED sensor-placement study on a
**more complicated geometry**: the adult-male thorax from EIDORS `shape_library`,
extruded to 3D. The standing standard (keep it): **electrodes, sensors, and the
anomaly/ROI all live at half height `z = height/2`** — the model is 3D but the design
plane is the mid-height slice.

Current state of this folder: `human_thorax_example.m` optimizes **only the sensor
radii** (`xi` = logit-transformed radius, free; `theta` and `z` fixed), with fixed
white noise. The result is that all sensors slide to the minimum allowed radius
`rmin = 1.01`. **This is the mathematically guaranteed outcome, not a bug**: with a
position-independent noise covariance, moving a sensor closer to the conductor
monotonically increases the signal (Biot–Savart kernel ~ 1/r²) and therefore the
information, for both A- and D-optimality. The radius-only design problem is trivial.

We want **non-trivial optimal positions**. Two tasks:

- **Task A — angle free, radius fixed.** Non-triviality comes for free from the
  geometry: the thorax cross-section is strongly non-circular, so a fixed-radius
  circle has an azimuth-dependent standoff from the body, and the ROI prior breaks
  the remaining symmetry. Expect a genuinely structured optimal azimuth distribution.
- **Task B — angle AND radius free.** With fixed noise this must collapse to `rmin`
  (verify and document that as the control). To make the radius non-trivial, extend
  the noise model to be **sensor-position dependent** (physically: a magnetometer very
  close to the body sits in the huge primary field of the injected currents —
  dynamic-range/gain-error noise grows with the local field magnitude). This creates
  a real signal-vs-noise tradeoff in radius and is a genuine methodological extension
  (new gradient term, new theory subsection).

## 2. Facts you need (verified 2026-07-06, planning session)

### 2.1 The shared optimizer is trustworthy
`functions/optimize_sensor_configuration.m` was inspected against the bug list from the
cylinder handoff:
- Cholesky solves use the CORRECT order everywhere: `Hinv = L'\(L\eye(n_elem))`
  (the wrong-order variant is present only as commented-out code).
- The 3-axis `dJ` assembly uses the corrected `dlambda{dim,p}` 2-D cell indexing
  (`compute_dJxyz_xyz_optimized`, comment "THIS WAS WRONG BEFORE" marks the fix).
- The coordinate-transform Jacobian is re-evaluated at the CURRENT iterate:
  `g()` calls `jac_coor_transf(sensor_locations)` every gradient evaluation. Your
  `jac_coord_transf` handle must therefore compute everything (r, theta) from the
  `sensor_locations` argument it receives, never from captured initial values.

So: build both tasks on `optimize_sensor_configuration` + generalized-coordinate maps
(the pattern already used in `human_thorax_example.m`), NOT on the self-contained
cylinder scripts. Still run finite-difference gradient checks (see §5) — they are the
arbiter, every time the parametrization or the noise model changes.

### 2.2 Generalized-coordinate conventions in the shared optimizer
- `q_to_x(q)` must return the cartesian vector `x = [x_1..x_M; y_1..y_M; z_1..z_M]`
  (3M×1, coordinate-major). `x_to_q(x)` is its inverse. See
  `sensor_locations_to_vector_cartesian` in the shared function.
- `jac_coor_transf(sensor_locations)` returns a `qmax × 3 × n_sensors` array:
  entry `(a, dim, m)` = ∂p_dim(m)/∂q_a(m) for sensor m (per-sensor chain rule).
- **Ordering trap:** the gradient is flattened as `reshape(dphidq.',1,[])` where
  `dphidq` is `qmax × n_sensors`. That yields
  `[dphi/dq1_sensor1..M, dphi/dq2_sensor1..M, ...]` — i.e. the q-vector is
  **coordinate-major, all sensors' first coordinate, then all sensors' second**.
  Your `q_to_x` / `x_to_q` maps must use the SAME ordering (for Task B:
  `q = [xi_1..xi_M; theta_1..theta_M]`), or the optimizer will silently scramble
  gradients across sensors. The FD check will catch this — run it.

### 2.3 Traps in the current `human_thorax_example.m` (do NOT inherit them)
1. **Prior-switching guard (worst trap):** the smooth lung prior is built only inside
   `if ~exist('data/sensor_positions_opt.mat','file')`. If that file exists (it does —
   from the radius-only runs), a re-run silently uses `Gamma_prior = variance_prior*speye`
   (identity prior) instead. Results across runs are then NOT comparable. In your new
   scripts, decouple "which prior" from "does a results file exist".
2. **Ad-hoc noise level:** `noise_std = max_B * 1e-1`. Replace with the validated
   **`d_target` calibration** from the cylinder study (noise std = d-th singular value
   of the whitened Jacobian `J0*sqrt(Gpr)` at the even/initial configuration; see
   `.tex` §"Noise calibration by target data-dominated rank" and
   `../example_anomaly_circle.m`). Also port the `roi_energy` diagnostic (mean ROI
   energy of the d leading right singular vectors) — it tells you whether the
   data-dominated subspace actually looks at the ROI.
3. `fmdl_inhom = fmdl_hom` — the "inhomogeneous" model is the same mesh; the
   inhomogeneity is only a material-property change on the lung regions
   (`mat_idx{2}`,`mat_idx{3}`). The CSG spherical anomaly is commented out because
   Netgen/EIDORS misidentified it as an electrode. **You do not need a meshed anomaly**:
   the OED cost is linearized at the prior mean and only needs `J`, `Gpr`, `Gn` — the
   ROI lives in the PRIOR (see §3). Don't fight the Netgen anomaly bug.
4. **Model cache:** meshes are cached in `HumanThorax/model/*.mat`. If you change any
   geometry parameter (height, mesh sizes, electrode count/radius/positions), delete
   the corresponding `.mat` or the change silently does nothing. Note
   `human_thorax_example.m` uses `height = 0.5`, `inj/meas = [0 3]`, 12 electrodes,
   16 sensors, coarse mesh (`maxsz_mesh = 1`, `maxsz_electrode = 0.15`); the cached
   models match those. `human_thorax_validation.m` used `height = 0.25`, `[0 1]`,
   2 sensors — a different (brute-force 2-sensor contour validation) setup, useful as
   a pattern for sanity-checking a new cost landscape, not as a base.
5. **Existing results files** `data/sensor_positions_init.mat` / `sensor_positions_opt.mat`
   are from the radius-only (trivial) run. Archive or rename them (e.g. into
   `data/radius_only_legacy/`); do not overwrite silently, and make your new scripts
   write task-specific filenames (see per-task specs).

### 2.4 Practical notes (carried over from the cylinder study — still apply)
- Run with **MATLAB R2024a** only: `"C:\Program Files\MATLAB\R2024a\bin\matlab.exe"`
  (Optimization Toolbox / `fminunc` is not installed for R2025a/b).
  Invocation pattern:
  `matlab -batch "cd('<repo>/studies/optimal_sensors_bayesian_approach/HumanThorax'); quick_test=true; <script_name>"`
- The scripts start a parpool and run `pctRunOnAll setupEidors(...)` — keep that.
- Support workspace-overridable knobs (`if ~exist('quick_test','var')...`,
  same for `max_iterations`, `n_starts`) as in `example_anomaly_circle_thetaz.m`,
  so runs can be driven from `-batch` without editing the script.
- Expected runtime: the thorax coarse mesh (~4–6k elements) with 16 sensors and 12
  electrodes is comparable to the cylinder runs (~7 s/iteration with 3425 elem,
  8 sensors, 16 stims). Budget minutes for quick tests, ~0.5–1.5 h for full runs.
  If iterations are much slower, drop to 8 sensors first, finer knobs later.
- Long background runs can survive a Windows sleep/wake (parpool reconnects), but
  check the log for a stall before trusting it.

## 3. Shared design choices for both tasks (decide once, in Task A, then reuse)

- **ROI/prior (primary):** Hyvönen "Case 1" style diagonal prior, as in the cylinder
  study — elevated std `prior_std_roi` on elements whose centroid lies within a ball
  of radius `roi_radius_factor * anomaly.radius` around
  `anomaly.position = (0.1, 0.1, height/2)`, and small `prior_std_background`
  elsewhere. This needs no meshed anomaly, matches the validated cylinder methodology,
  and makes results directly comparable across geometries. Start from the cylinder's
  final regime (`bg_std = 0.005*sigma_bg`, `roi_std = 1.0*sigma_bg`, `d_target = 100`)
  and re-tune with the `roi_energy` + "even-config ROI variance reduction vs prior"
  diagnostics exactly as documented in the cylinder handoff (§"Tuning sweep performed").
- **ROI/prior (optional secondary, nice-to-have):** the existing smooth lung prior
  (lungs as the ROI — a physiological target). Only after the primary works; report
  it as a comparison, don't let it block the task.
- **Criteria:** both `a-opt` and `d-opt`, 3-axis MDEIT (`'mdeit3'`), like everything
  else in this study.
- **Baselines to compare against:** (i) evenly-spaced sensors at the same fixed radius
  (this is also the optimizer's start), and (ii) the "boundary configuration" from
  `find_sensor_positions_boundary` (sensors at constant standoff `delta` from the body
  along equiangular rays) — the thorax analogue of "as close as allowed", i.e. the
  trivial answer. A non-trivial claim = the optimizer beats or meaningfully differs
  from BOTH.
- **Report per configuration:** phi(even), phi(opt), % trace reduction (a-opt),
  information gain in nats/bits (d-opt), ROI posterior-variance ratio opt/even,
  number of data-dominated modes, `roi_energy`. Same table format as the cylinder
  handoff updates.

---

## 4. TASK A (next agent — do this one): azimuth free, radius fixed

**Question:** with 16 sensors on a circle of fixed radius `R0 = 1.5` at `z = height/2`
around the thorax, what is the optimal azimuth distribution, and how much does it gain
over even spacing? Is it non-trivial (structured by the thorax shape + ROI), unlike
the cylinder case where A-opt was nearly rotation-invariant?

**Hypothesis to test (state the outcome explicitly in your update):** on the cylinder,
A-opt gains were tiny (~4.6%) because a fixed circle around a rotationally symmetric
body conserves the average sensitivity budget under azimuth moves. The thorax breaks
that symmetry twice (non-circular boundary → azimuth-dependent standoff; off-center
ROI). Expect sensors to migrate toward (a) the ROI side and (b) azimuths where the
body is close to the ring; expect a LARGER a-opt gain than the cylinder's 4.6%.

Steps:

1. **New script** `human_thorax_theta.m` (this folder), based on
   `human_thorax_example.m` but with the traps of §2.3 removed (standalone prior
   construction, `d_target` noise calibration, task-specific output filenames).
2. **Parametrization:** `q = theta (M×1)`, radius fixed at `R0 = 1.5`, `z = height/2`.
   - `map_q_to_x_theta(q)`: `x = [R0*cos(q); R0*sin(q); z0]`.
   - `map_x_to_q_theta(x)`: `q = atan2(y, x)` (unconstrained optimization over
     periodic theta is fine — no box needed).
   - `jac_coord_transf_theta(sensor_locations)`: `qmax = 1`;
     `J(1,1,m) = -r_m*sin(theta_m)`, `J(1,2,m) = r_m*cos(theta_m)`, `J(1,3,m) = 0`,
     with `r_m, theta_m` computed from the CURRENT `sensor_locations` argument
     (this is `\eqref{eq:chainrule1d}` in the `.tex`).
   - Round-trip assert `sensor_positions_0 == q_to_x(x_to_q(...))` like the existing
     script does.
3. **FD gradient check** (mandatory, see §5) at the even start AND at one random
   perturbed configuration, for both a-opt and d-opt.
4. **Runs:** quick test (3 it.) → tune prior/noise per §3 diagnostics → full runs
   (both criteria, `max_iterations = 40`, converged or justify) → multistart check
   (≥3 random-azimuth starts; draw random theta uniform on [0,2π)).
   Watch for sensor collisions/pile-ups; if two sensors coalesce, note it — the
   `.tex` has a repulsion-penalty section (§"Sensor repulsion penalty") you can
   enable, but first report the un-penalized optimum.
5. **Figures** (save `.fig` + `.png` here): FEM + sensors for even / a-opt / d-opt /
   boundary-config overlaid on the thorax slice (top view, like the existing
   `human_thorax_opt_config.fig`); polar/azimuth plot of optimal angles vs ROI azimuth;
   mid-height posterior-std slice for even vs optimized.
   Data → `data/thorax_theta_results.mat`.
6. **Interpret:** where do sensors go — toward the ROI, toward the close-standoff
   azimuths (left/right flanks), or a mixture? Does d-opt cluster tighter than a-opt
   (as on the cylinder)? Quantify vs both baselines of §3.
7. **Update `../sensor_position_optimization_theory.tex`:** add a new section
   "Case study: human thorax, sensors on a fixed circle at half height" (suggest
   `\label{sec:case-study-thorax}`) after the existing case-study sections. Reuse the
   established notation (`\Phi_A`, `\Phi_D`, `\eqref{eq:chainrule1d}`,
   `\eqref{eq:roiprior}`, `\eqref{eq:noisecal}`, `\eqref{eq:roienergy}`). Describe the
   geometry, the two baselines, results table, and the symmetry-breaking discussion.
   Compile cleanly (2 passes), clean up `.aux/.log/.out`.
8. **Update THIS file:** append an UPDATE section with parameters, results table,
   diagnosis, open items — and adjust the Task B section below if anything you
   learned changes it (e.g. final prior/noise regime, mesh/sensor-count choices,
   runtimes). Task B's agent starts from your even/theta-opt numbers as baselines.

## 5. Verification protocol (both tasks, non-negotiable)

- Central finite-difference gradient check on ≥3 components covering EVERY block of
  the design vector (theta block; and for Task B the xi block too — pass an explicit
  index pool as `example_anomaly_circle_thetaz.m` does with `check_gradient_fd`).
  Target rel. err ≤ 1e-6 (the study's standard is 1e-7–1e-8). If it fails, stop and
  fix — do not run optimizations on an unverified gradient.
- Confirm `phi(even)` is identical (all printed digits) between scripts that claim
  the same starting configuration.
- Every reported optimum: state iterations, stopping reason (converged vs MaxIter),
  and multistart agreement, like the cylinder handoff tables do.

## 6. TASK B (following agent): azimuth AND radius free

**Question:** free `q = [xi_1..xi_M; theta_1..theta_M]` (coordinate-major, §2.2), with
`r_m = rmin + (rmax-rmin)*sigmoid(xi_m)`, `rmin = 1.01`, `rmax = 2.0`, `z = height/2`.
Where do sensors go, and can we make the radial optimum non-trivial?

**Part B1 — control (fixed noise).** New script `human_thorax_rtheta.m`. Combine the
existing xi-map (already in `human_thorax_example.m`) with Task A's theta-map into a
`qmax = 2` transform (`J(1,:,m)` = xi row: `cos(theta_m)*drm/dxi_m`,
`sin(theta_m)*drm/dxi_m`, 0; `J(2,:,m)` = theta row as in Task A — all from current
positions). FD-check both blocks. Run both criteria. **Expected result: radial collapse
to the boundary of the feasible set (all `r → rmin`), azimuths roughly as in Task A.**
Document this as the control: it confirms the collapse is a property of the fixed-noise
model, not of the radius-only parametrization. Compare the collapsed configuration
against Task A's fixed-`R0` optimum and the boundary config — with fixed noise, closer
should win. Watch for sigmoid saturation (|xi| large → vanishing gradient): if sensors
pin at `rmin` with huge |xi|, report the achieved radii, and cap displayed convergence
claims accordingly.

**Part B2 — position-dependent noise (the real contribution).** Make the per-channel
noise std grow with the local primary (homogeneous-model) field:

    sigma_i(p_m) = sigma_floor + c * ||B_hom(p_m)||        (channel i of sensor m)

where `B_hom(p_m)` is the 3-axis homogeneous forward field at sensor m's position
(available from `fwd_solve_mdeit` on the homogeneous image; per-stim or aggregated —
start with the per-sensor RMS over stims/axes for smoothness, and say what you chose).
Physics: near the body the primary field of the injected currents is orders of
magnitude larger than the anomaly signal; finite dynamic range / gain error makes the
effective noise scale with it. `B_hom ~ (standoff)^-2` while the ROI signal decays
with distance-to-ROI — so there is a genuine interior optimum in radius.

- **Cost:** `Gn(p)` is diagonal → cheap to rebuild each iterate. Choose
  `sigma_floor` via the `d_target` calibration at the initial config, and `c` such
  that at `r = rmin` the field-dependent term dominates (say 3–10× floor) while at
  `r = rmax` it is subdominant — sweep `c` over a decade and report the sensitivity
  of the optimal radii to it. This knob IS the result; be explicit about it.
- **Gradient:** two routes, in order of preference:
  1. **Analytic (target):** the criteria gain an extra term from `dGn^{-1}/dp`.
     For `Phi_A = tr(H^{-1})`, `H = J' Gn^{-1} J + Gpr^{-1}`:
     `dPhi_A = -tr(H^{-2} dH)` with
     `dH = dJ' Gn^{-1} J + J' Gn^{-1} dJ + J' d(Gn^{-1}) J`; the first two pieces are
     what the code already computes, the new piece is
     `d(Gn^{-1})_ii = -2*(dsigma_i/dp_m)/sigma_i^3` on the channels of sensor m only.
     `dB_hom/dp_m` comes from the existing `dRmkj_xyz_optimized` kernels (B is linear
     in the R-matrices for fixed u_hom). Analogous term for `Phi_D` via
     `dPhi_D = -tr(H^{-1} dH)`. Derive carefully, write it in the `.tex` (new
     methodology subsection, e.g. "Position-dependent noise and the extended
     gradient"), and FD-verify — the FD check is the arbiter, as always.
  2. **Fallback if the analytic term stalls you:** freeze `Gn` at each outer
     iteration (recompute noise from current positions, optimize a few inner
     iterations with fixed `Gn`, repeat — a lagged/majorized scheme). Cheaper to
     implement, no new gradient term, but document that the FD check then only holds
     for the inner (fixed-`Gn`) problem, and verify the outer loop decreases the TRUE
     (position-dependent-noise) cost each round. Do NOT silently drop the `dGn` term
     inside a single `fminunc` run — that produces a wrong gradient that the optimizer
     will chase into nonsense.
- **Runs:** FD check → quick test → full runs both criteria → multistart (random
  theta uniform; random xi ~ 2*randn like the thetaz study did for eta) → sweep `c`.
- **Deliverables:** `data/thorax_rtheta_results.mat` (+ control results), figures
  (sensor maps for control vs noise-aware optimum; radius histogram / r-vs-theta
  scatter; posterior-std slices), `.tex` update (methodology subsection + thorax
  case-study results), and an UPDATE section in this file.

**Success criterion for Task B:** the noise-aware optimum places sensors at radii
strictly inside `(rmin, rmax)` (not pinned at either end), reproducibly across
multistarts, with a stated `c`-sensitivity — i.e. a non-trivial radial design driven
by a physically motivated tradeoff, plus the control showing why the naive model
cannot produce one.

---

## 7. Updates

(Task agents: append `## UPDATE (date, model): <task>` sections here. Include: final
parameters, FD-check numbers, result tables, stopping reasons, diagnosis, `.tex`
sections touched, files written, and anything the next agent must know.)

## UPDATE (2026-07-07, Claude Opus 4.8): TASK A done (azimuth free, radius fixed)

### Files written
- `human_thorax_theta.m` (this folder) — the Task A script. Drives the shared
  `optimize_sensor_configuration.m` (Sec 2.1) with theta coordinate maps
  (`map_q_to_x_theta` / `map_x_to_q_theta` / `compute_jac_coord_transf_theta`), plus
  self-contained FD-validated cost/grad helpers (copied from `example_anomaly_circle.m`)
  used ONLY for the FD check, `phi(even/opt/boundary)`, and the posterior/roi_energy
  diagnostics. All three §2.3 traps removed: prior built unconditionally, `d_target`
  noise calibration, task-specific filenames. Workspace knobs: `quick_test`,
  `max_iterations`, `n_starts`, `d_target`, `prior_std_background_factor`,
  `prior_std_roi_factor`, `roi_radius_factor`, and `tune_sweep` (early-exit diagnostics table).
  Writes a `data/thorax_theta_checkpoint.mat` after every (mode,start) — a shutdown mid-run
  cannot waste completed optimizations (one already happened; see note at end).
- `data/thorax_theta_results.mat` — full results struct (theta_even, theta_opt{a,d},
  theta_starts, phi_even/opt/boundary, posterior variances, roi, noise/roi_energy, params).
- Figures (this folder): `thorax_theta_sensor_configs.{fig,png}` (even/a-opt/d-opt/boundary
  over the thorax slice, ROI star marked), `thorax_theta_azimuths.{fig,png}` (polar, ROI at
  45°), `thorax_theta_posterior.{fig,png}` (mid-slice std: prior/even/a-opt).
- Legacy radius-only files moved to `data/radius_only_legacy/` (Sec 2.3.5).
- `.tex`: new `\section{Case study: human thorax...}` `\label{sec:case-study-thorax}` added
  before "Symbol--code correspondence", reusing `\eqref{eq:chainrule1d,roiprior,noisecal,roienergy}`.
  Compiled clean (2 pdflatex passes, no undefined refs; only the pre-existing OMS/cmtt font
  warning). `.aux/.log/.out` cleaned; PDF regenerated.

### Final parameters (the chosen regime)
`R0=1.5`, `zs=height/2`, `height=0.5`, 16 sensors, 12 electrodes, inj/meas `[0 3]`, coarse
cached mesh (`n_el=5955`, `n_nodes=1525`, `n_stim=12`, `n_data=576`). Linearization point =
homogeneous-lung image (`bg=1.0`, `lungs=0.3`). Prior: `sigma_bg=0.005`, `sigma_roi=1.0`,
**`roi_radius_factor=2.0`** (see tuning note), `c0=(0.1,0.1,H/2)`. Noise: `d_target=100`
→ `noise_std=1.384e-6`. `max_iterations=40`.

### Tuning (Sec 3 diagnostics — sweep is in the script, `tune_sweep=true`)
On this coarse mesh a factor-1 ROI ball holds only **14 elements** → roi_energy 0.11
(pathological, ceiling ≈ n_roi/d_modes). Enlarged to **factor-2 (129 elements**, comparable to
the cylinder's ~100) → **roi_energy 0.33, even-config ROI variance reduction 26%** (cylinder's
chosen point was 0.54 / 30%). Kept `bg=0.005`, `d_target=100` as on the cylinder. Sweep table
(bg=0.005): rf1/d100 → E=0.11,red=76%(!); rf2/d100 → E=0.33,red=26%; rf3/d100 → E=0.60,red=11%;
rf4/d20 → E=1.0,red=2%. Higher roi_energy trades off against having any headroom — rf2/d100 is
the honest analogue of the cylinder regime.

### FD gradient check (Sec 5, the arbiter) — PASS
Central FD on sensors {1,7,8,11,16}, both criteria, even AND perturbed configs:
**rel. err 1.6e-8 … 3.0e-6** (one d-opt perturbed point at 2.8e-6; all others ≤ 5e-7). The
self-contained `phi(even)` matches the shared optimizer's printed `f0` to all digits
(a-opt 9.55589e1, d-opt −6.204432e4), cross-validating the two gradient paths.

### Results (40 it; stopping reasons noted)
| Criterion | phi(even) | phi(opt) | phi(boundary) | gain vs even | vs boundary |
|---|---|---|---|---|---|
| a-opt tr(H⁻¹) | 95.559 | 95.415 | **71.080** | **0.15%** | −25% (opt worse) |
| d-opt −logdetH | −62044.3 | −62048.2 | **−63124.2** | **3.9 nats** | −1076 nats (opt worse) |
- a-opt even-start: "Local minimum found". d-opt: "Local minimum possible" (grad-limited).
- Multistart: every random-azimuth start (2 in this run + 3 more in the earlier interrupted
  run, all uniform on [0,2π)) converges to the **identical** value — a-opt 95.41485 (even-start
  lands at 95.41765, i.e. random starts find a hair better), d-opt −62048.2287 to 7 digits. The
  even-spacing start is essentially at the optimum; additional starts add no information.
- Posterior ROI variance: opt/even = 0.9985, even/prior = 0.7397.

### Diagnosis (state the hypothesis outcome explicitly)
**The hypothesis is REFUTED: the thorax a-opt azimuth gain (0.15%) is much SMALLER than the
cylinder's (4.6%), not larger.** Cause is geometric and quantified: the specified ROI centre
`c0=(0.1,0.1)` sits at radius **0.14** — nearly central relative to the `R0=1.5` ring — so the
azimuthal ROI→sensor distance contrast is only **±9.4%**, versus the cylinder anomaly (0.5R vs
1.5R ring) at **±32.5%**. The ROI is effectively *more* central here, so the fixed-radius circle
conserves the averaged sensitivity budget under azimuth moves even more nearly than on the
cylinder. The body's strongly non-circular boundary (mid-height boundary radius ranges ~0.5–1.0,
node extent x∈[−0.94,1.00], y∈[−0.65,0.68]) does NOT rescue this: with sensors far out at R0=1.5,
δB is set by sensor-to-ROI distance, not sensor-to-surface standoff, so boundary shape barely
enters ROI sensitivity. The optimum is still weakly *structured* in the expected direction
(sensors within ±45° of the ROI azimuth: 3→5 for both criteria; a-opt puts 9/16 sensors in the
0–144° arc, max single move ~36°; d-opt circ-mean biased toward ROI), just with mild magnitude
(resultant |Σe^{iθ}|/M ≈ 0.17 a-opt, 0.07 d-opt). **The decisive lever is radial standoff, not
azimuth** — the boundary ("as close as allowed") config beats every fixed-R0 circle by 25% (a-opt)
/ 1076 nats (d-opt). This is exactly the setup for Task B.

### IMPORTANT adjustments for Task B (read before starting)
1. **Reuse this exact regime** for comparability: `bg=0.005`, `roi=1.0`, **`roi_radius_factor=2.0`**
   (NOT 1.0 — the factor-1 ROI is under-resolved on the coarse mesh), `d_target=100`,
   `c0=(0.1,0.1,H/2)`, 16 sensors, cached models. Your baselines are the numbers above:
   a-opt even 95.559 / theta-opt 95.415; d-opt even −62044.3 / theta-opt −62048.2.
2. **Feasibility-box caveat (matters for B1's collapse claim):** the boundary/proximity optimum
   lives at radii **below** Task B's `rmin=1.01` (the body only reaches radius ~1.0 on its widest
   axis; boundary sensors sit at body+0.05 ≈ r 0.5–1.0). So the huge proximity advantage seen in
   Task A is *unreachable* inside `[rmin,rmax]=[1.01,2.0]` — the whole feasible annulus is already
   outside the body. B1 (fixed noise) will still collapse to `r→rmin=1.01` (closer = more info),
   but note rmin is only marginally inside the body's bounding circle and the standoff at rmin is
   very azimuth-dependent (tight on the ±x flanks, loose on ±y). Expect θ to again do little
   (near-central ROI) so the (r,θ) optimum ≈ {all r=rmin, θ ≈ Task A's mild clustering}.
3. **B2 is where the physics is:** since proximity dominates so strongly, the position-dependent
   noise term (σ grows with the local primary field ~ standoff⁻²) is exactly what's needed to
   create the interior radial optimum; the `c`-sweep is the real result, as the plan says.
4. The `tune_sweep`/checkpoint/`cost_at_positions`/`compute_jac_coord_transf_theta` machinery in
   `human_thorax_theta.m` is directly reusable; for B combine the xi-map (already in
   `human_thorax_example.m`) with the theta-map into the qmax=2 transform (coordinate-major
   `q=[xi_1..xi_M; theta_1..theta_M]`, Sec 2.2) and FD-check BOTH blocks.

### Runtime / operational notes
- Full run (both criteria, 3 starts, 40 it): ~50 min on this machine; a-opt ~500 s/start,
  d-opt ~300 s/start (direct solves, `n_el=5955`). Quick test (`quick_test=true`) ~4 min.
- **A background run was interrupted by a Windows shutdown mid-optimization** (the parpool does
  NOT always survive shutdown, unlike sleep/wake). The per-(mode,start) checkpoint added to the
  script recovers completed starts; `data/thorax_theta_results.mat` is written only at the very
  end, so trust the checkpoint if a run dies. The reported numbers are from a clean completed
  re-run (deterministic, `rng(1)`); they reproduce the interrupted run's partial values exactly.
- Run with MATLAB R2024a only (Optimization Toolbox). The pre-existing `sparse`-shadowing warning
  from `setupEidors` is harmless.
