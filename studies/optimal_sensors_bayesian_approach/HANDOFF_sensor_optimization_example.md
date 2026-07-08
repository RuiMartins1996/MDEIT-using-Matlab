# Handoff: Bayesian sensor-position optimization example (MDEIT)

Session summary written 2026-07-06. Goal: hand this to a new session (Claude Sonnet)
to finish the work. Everything below is in
`studies/optimal_sensors_bayesian_approach/` unless stated otherwise.

## Goal of the task

Create a non-trivial, self-contained example (similar to `main.m`) that optimizes
magnetometer positions for MDEIT following the Bayesian OED approach of Hyvönen,
Seppänen & Staboulis (arXiv:1404.7300, PDF in this folder):

- Sensors constrained to a **circle at half height** of the 3D cylindrical FEM model,
  at radius `rs = 1.5*tank_radius` (only the azimuth angle θ of each sensor is free).
- A **conductive spherical anomaly** off-center at half height
  (`position = [R/2, 0, H/2]`, radius `R/3`, conductivity 5× background).
- Both **A-optimality** (`trace(H^-1)`) and **D-optimality** (`-logdet(H)`),
  `H = J'*inv(Gn)*J + inv(Gpr)`, linearized at the **prior mean** (homogeneous image).
- Requirement: the optimized positions must give a **considerable cost reduction**
  vs. evenly spaced sensors; if not possible, understand why.

## What was done

### 1. Theory document (finished)
`sensor_position_optimization_theory.tex` (+ compiled PDF): full write-up of the
theory (CEM, Biot–Savart discretization with R/G/Γ matrices, linearized Gaussian
posterior, A/D-optimality, gradient derivation with W = Γn⁻¹JH⁻² / Γn⁻¹JH⁻¹,
dJ/dp split into adjoint-derivative and kernel-derivative parts, cylindrical chain
rule, symbol↔code table). Compiles cleanly with MiKTeX pdflatex (2 passes).

### 2. Example script (working, needs final tuning)
`example_anomaly_circle.m` — self-contained script (local functions at the end,
adapted from `main.m`). Key design:

- Only free variables: sensor azimuths `theta` (M×1). Chain rule
  `dphi/dtheta = -rs*sin(th)*dphi/dx + rs*cos(th)*dphi/dy` evaluated at the
  **current** iterate.
- Cost and gradient share ONE Jacobian evaluation per iterate
  (`costgrad_theta_3axis`), used with `fminunc` quasi-Newton,
  `SpecifyObjectiveGradient=true`.
- 3-axis MDEIT (J = [Jx;Jy;Jz]); sensor axes are identity (set via
  `model_parameters.sensorPositions`, see `mk_mdeit_model`/`assign_magnetometers`).
- Prior: diagonal, `prior_std_roi` inside a ball around the (suspected) anomaly
  location, `prior_std_background` elsewhere (Hyvönen "Case 1" style).
- Noise: white; **auto-calibrated** from the spectrum of the whitened Jacobian
  `J0*sqrt(Gpr)` at the even configuration so that exactly `d_target` modes are
  data-dominated (`noise_std = s(d_target)`).
- Diagnostics printed: #data-dominated modes, mean ROI energy fraction of the
  data-dominated right singular vectors (added last; NOT yet exercised in a run),
  posterior variance sums over ROI/background for prior/even/optimized.
- `quick_test = true` in base workspace → 3 iterations + finite-difference
  gradient check (3 random components, central differences).
- Figures (saved to `figures/`): FEM + sensors (even/a-opt/d-opt), polar plot of
  azimuths vs anomaly azimuth 0, mid-height slice scatter of posterior std.
- Results saved to `data/example_anomaly_circle_results.mat`.

### 3. Bugs found in `main.m` (NOT yet fixed there; fixed in the example)
1. **`compute_dJxyz_xyz` (line ~1627)**: uses `dlambda{p}` to index a 3×3 cell
   `dlambda{dim,p}` — linear indexing gives `dlambda{p,1}`, i.e. the wrong
   adjoint derivative for most (dim,p) combinations in the 3-axis gradient.
2. **Cholesky solve order for H⁻¹** (e.g. line ~2152): `Y = L'\eye; invH = L\Y`
   computes `inv(L)*inv(L')` = (LᵀL)⁻¹, but H = L*Lᵀ so H⁻¹ = L⁻ᵀL⁻¹.
   Correct order: `X = L\eye; invH = L'\X`. Affects ALL gradient functions in
   main.m (a-opt and d-opt, 1-axis and 3-axis). Costs are unaffected.
3. **Stale coordinate transform** (line 157): the cylindrical chain-rule Jacobian
   is computed once at the initial sensor locations and captured in the gradient
   handles; it depends on the current (r,θ) and must be re-evaluated per iterate.
4. Performance: `compute_r_matrices_local` in main.m is an interpreted
   double loop (sensors × elements × 35 quad pts); the example has a vectorized
   version (~100× faster). Worth back-porting.

With fixes 1+2 applied, the example's analytic gradients match central finite
differences to **rel. err 1e-7…1e-8** (verified in MATLAB runs, see below).

### 4. Verification runs (MATLAB R2024a required!)
- Optimization Toolbox is only installed for **R2024a**
  (`C:\Program Files\MATLAB\R2024a\bin\matlab.exe`); R2025a/b lack `fminunc`.
- Run command used:
  `matlab -batch "cd('<this folder>'); quick_test=true; example_anomaly_circle"`
- Model: 3425 elements, 817 nodes, 16 stims (1 ring × 16 electrodes at H/2),
  8 sensors. ~7 s per optimizer iteration incl. gradient (parpool, 6 workers).
- Gradient checks pass every run (1e-6…1e-8).

### 5. Results so far and DIAGNOSIS (the important part)
Three quick-test configurations tried (3 optimizer iterations only!):

| config | d.o.m. modes | a-opt reduction | ROI post var opt/even | d-opt gain |
|---|---|---|---|---|
| noise=maxB/20, roi_factor=1.5, bg_std=0.05σ | 1 | — | — | — |
| d_target=30, roi_factor=1.5, bg_std=0.05σ | 29 | 0.8% | 0.992 | 14.0 nats |
| d_target=100, roi_factor=1.0, bg_std=0.05σ | 99 | 1.2% | 0.987 | 3.7 nats |

Why the a-opt reduction is small so far:
- The data can constrain at most ~#(data-dominated modes) directions of the prior.
  With ROI ≫ that many elements, no sensor arrangement helps much (run 2).
- Even with 99 modes (run 3), the ROI posterior variance only dropped ~10% below
  the prior → the **data-dominated modes are not aligned with the ROI**. Raw MDEIT
  sensitivity is much higher near the electrodes/boundary than at the interior
  anomaly, so with weak prior contrast (bg std = 0.05σ) the data budget is spent
  "learning" the (already known) background.
- Conclusion: increase the prior contrast (bg std → 0.01σ or lower) so the whitened
  Jacobian's leading singular vectors point into the ROI. The `roi_energy`
  diagnostic (mean ROI energy of the data-dominated right singular vectors,
  printed by the current script) directly measures this: want it close to 1.

### Current parameter state in `example_anomaly_circle.m`
`d_target=100`, `prior_std_background=0.05*background`, `prior_std_roi=1.0*background`,
`roi_radius_factor=1.0`, `n_sensors=8`, `max_iterations=40`, `n_starts=1`,
`alpha_rep=0`. (A pending change to `bg_std=0.01`, `d_target=60` was NOT applied.)

## Plan for the next session (suitable for Claude Sonnet)

1. **Tune the prior/noise balance for ROI alignment.** In
   `example_anomaly_circle.m` set `prior_std_background = 0.01*background_conductivity`
   and try `d_target` in {60, 100}. Run the quick test
   (`matlab -batch "cd(...); quick_test=true; example_anomaly_circle"`, R2024a!)
   and check the printed `mean ROI energy of the data-dominated modes` — want
   ≳0.7. If still low, lower bg std to 0.005σ, or as a fallback implement
   **ROI-weighted A-optimality** (cost = trace of the ROI block of H⁻¹; gradient
   W = Γn⁻¹·J·H⁻¹·diag(roi)·H⁻¹ — a 5-line change in `costgrad_theta_3axis`).
2. **Confirm the regime**, then check posterior sums: `even` should now reduce the
   ROI variance substantially below prior (say 30–60%), leaving headroom for the
   optimizer.
3. **Full run** (`quick_test=false`; 40 iterations, both criteria; ~15–25 min):
   `matlab -batch "cd('<folder>'); example_anomaly_circle"`. Verify:
   - a-opt: trace reduction vs even, target "considerable" (≥15–30%);
     also report ROI posterior variance ratio optimized/even.
   - d-opt: information gain in nats/bits (tens of nats is already substantial).
   - Sensors should cluster toward azimuth 0 (anomaly side) — check the polar figure.
4. **If reduction still modest, iterate on physics** (in this order): more sensors
   is NOT the point — instead (a) raise `d_target` (deeper SNR), (b) move ROI
   closer to the wall or enlarge the anomaly, (c) ROI-weighted A-opt from step 1.
   Each knob's effect is explained by the printed diagnostics.
5. **Optionally increase `n_starts`** (multistart with random azimuths) to check
   the even-start optimum isn't a poor local minimum; report the best.
6. **Report bugs back into `main.m`** (user decision whether to apply): the
   `dlambda{dim,p}` indexing, the Cholesky solve order (both described above),
   the stale cylindrical chain-rule Jacobian, and optionally back-port the
   vectorized `compute_r_matrices_local`.
7. Keep `sensor_position_optimization_theory.tex` in sync if anything conceptual
   changes (e.g. add a short section on ROI-weighted A-optimality if implemented).

## Practical notes
- Always run with **R2024a** (`"C:\Program Files\MATLAB\R2024a\bin\matlab.exe"`).
- The mesh/model is cached under the repo `models/` folder; changing
  `model_parameters` (incl. sensorPositions) re-meshes via Netgen automatically.
- Quick-test logs from this session are in the session scratchpad
  (`example_quicktest*.log`); the run prints everything needed to the console.
- `fminunc` display is `iter`; each iteration ≈ 7 s on this machine (6 workers).

## UPDATE (2026-07-06 session, Claude Sonnet 5): tuning finished, final results

### Final parameters in `example_anomaly_circle.m`
`prior_std_background = 0.005*background_conductivity`, `prior_std_roi =
1.00*background_conductivity`, `roi_radius_factor = 1.0`, `d_target = 100`,
`n_sensors = 8`, `max_iterations = 40`, `n_starts = 1`. Anomaly stays at the
originally specified `[R/2, 0, H/2]`, radius `R/3`.

### Tuning sweep performed
| bg_std | d_target | roi_energy | even ROI var reduction vs prior |
|---|---|---|---|
| 0.05σ | 100 | 0.54 | 15% |
| 0.01σ | 60  | 0.54 | 18% |
| 0.005σ | 60 | 0.65 | 22% |
| 0.005σ | 100 | 0.54 | **30%** ← chosen |
| 0.005σ | 150 | 0.46 | 38% |

`roi_energy` and "even variance reduction" trade off against each other as
`d_target` grows (more, less ROI-aligned modes cross the noise threshold).
`bg_std=0.005`, `d_target=100` was kept as a representative point already
past the "considerable data effect" threshold from the plan (≥30% reduction);
pushing `d_target` further (150) did not meaningfully change the a-opt
outcome (see below), so it wasn't worth the reduced `roi_energy`.

### Final full-run (40 it., converged to `OptimalityTolerance=1e-8`) results
- **a-opt**: phi(even) = 0.28791, phi(opt) = 0.27476 → **trace reduction = 4.6%**.
- **d-opt**: phi(even) = -55552.7, phi(opt) = -55570.3 →
  **extra information = 17.57 nats (25.35 bits)**.
- ROI posterior variance ratio optimized/even = 0.954 (a-opt).
- Multistart check (`n_starts=3`): all 3 a-opt starts converge to the exact
  same minimum (0.274765); the 3 d-opt starts agree to within 0.03% (best:
  18.00 nats). The even-spacing start is essentially already at the global
  optimum for this design space — not a poor local minimum.
- Sensor azimuth figure confirms both configurations migrate toward azimuth 0
  (anomaly side), with d-opt clustering visibly tighter than a-opt.

### Diagnosis: why a-opt reduction stays modest (~4-5%) but d-opt is strong
This was chased down as the "important part" of the task:
- Raising `d_target` from 60→100→150 did **not** meaningfully change the
  converged a-opt reduction (stayed ~4-6% across all three), ruling out a
  data-budget/SNR explanation.
- **Confirming experiment**: moved the anomaly from `R/2` to `0.6R` (closer to
  the tank wall, same radius, same prior/noise settings) as a one-off
  diagnostic (reverted afterward — not part of the delivered example). Result
  (3-iteration quick test): a-opt reduction only ticked up from 2.7%→4.4%,
  while d-opt's extra information nearly doubled (13.81→24.79 nats).
- **Conclusion**: A-optimality (`trace(H^-1)`, an *averaged* posterior
  variance) is close to rotation-invariant when sensors are constrained to a
  fixed-radius, fixed-height circle: moving all sensors around the ring
  changes each sensor's *individual* distance to the anomaly, but the
  *average* sensitivity budget over the ring is nearly conserved by symmetry,
  regardless of anomaly depth. D-optimality (`-logdet(H)`) instead rewards
  reducing *redundancy/correlation* between sensor measurements — clustering
  sensors on the anomaly side reduces how much they duplicate each other's
  information, which is a genuinely different (and more available) resource
  than raw averaged sensitivity. This is why d-opt shows a "considerable"
  gain (tens of nats, as anticipated by the plan) while a-opt does not, and
  why enlarging the SNR budget alone can't fix it — it's a symmetry property
  of the fixed-circle parametrization, not an SNR-starved regime.
- This explanation, and the confirming diagnostic, should be folded into the
  final write-up/report to the user rather than into
  `sensor_position_optimization_theory.tex` (no new formula or methodology
  was introduced — ROI-weighted A-optimality was considered but not
  implemented, since the ROI already dominates the (unweighted) trace by
  construction: prior ROI variance sum is 4.094e-1 vs background 1.8e-4, i.e.
  >99.9%, so an ROI-weighted trace would differ negligibly from the current
  unweighted one).

### Remaining open items
- `main.m` bugs (dlambda indexing, Cholesky order in all 4 gradient
  functions, stale chain-rule Jacobian, vectorized `compute_r_matrices_local`)
  were **applied to `main.m` directly** in this session, per the user's
  request. Verified with `checkcode('main.m')` (no syntax errors, only
  pre-existing style lint warnings). NOT verified with a full end-to-end
  `main.m` run (it is a large multi-branch driver script with `fmincon` and
  random initial sensor positions; the fixes are mathematically identical to
  the ones already validated via finite-difference gradient checks in
  `example_anomaly_circle.m`, rel. err. 1e-7-1e-8, so correctness carries
  over, but a smoke test of `main.m` itself would still be worthwhile before
  relying on it for new results).
- `sensor_position_optimization_theory.tex` was reviewed but not modified
  (no conceptual change this session).

## UPDATE (2026-07-06 session, Claude Sonnet 5): freeing z, follow-up task done

This continues directly from the handoff in
`Optimizing_Electrode_Positions_in_Electrical_Imped.md` (same file this section lives next to).
That handoff asked: does letting each sensor's height `z_m` vary in addition to its azimuth
`theta_m` (radius still fixed at `rs = 1.5R`) produce a *more considerable* A-optimal gain than the
theta-only baseline above (4.6% / 17.57 nats)?

### What was built
New sibling script `example_anomaly_circle_thetaz.m`, a copy of `example_anomaly_circle.m` with:
- Design vector `q = [theta(1:M); eta(1:M)]` (2M = 16 free variables for M=8 sensors).
- `z_m = zmin + (zmax-zmin)*sigmoid(eta_m)`, `zmin=0`, `zmax=tank_height`, so `eta_m=0` reproduces
  the old fixed `z_s = tank_height/2` exactly — `q0 = [theta_even; zeros(M,1)]` is bit-for-bit the
  same starting configuration as the theta-only baseline (confirmed: `phi_even` matches to all
  reported digits in both scripts).
- Chain rule: `dphi/dtheta_m` unchanged from the theta-only script; new term
  `dphi/deta_m = dphidp(3,m) * (zmax-zmin)*sigmoid(eta_m)*(1-sigmoid(eta_m))`, reusing the
  z-partial `dphidp(3,:)` that the theta-only script already computed and discarded.
- `check_gradient_fd` extended with an optional `idx_pool` argument so the smoke test explicitly
  exercises **both** the theta half and the eta/z half of `q` (previously random sampling over the
  whole vector could miss one half by chance).
- New diagnostic figure: unrolled `(theta_m, z_m)` scatter (replaces the polar-only view, which
  can no longer show height).
- Results saved to `data/example_anomaly_circle_thetaz_results.mat` (theta-only results file
  untouched).

### Gradient check (quick_test = true)
All rel. errors 3e-8 to 9e-7 for both halves of `q` (a-opt and d-opt) — same standard as the rest
of this study. See "Practical notes" above for how to invoke MATLAB R2024a.

### Results
| Run | a-opt: phi(even)&rarr;phi(opt) | trace reduction | d-opt: phi(even)&rarr;phi(opt) | info gain |
|---|---|---|---|---|
| theta-only baseline (for reference) | 0.28791 &rarr; 0.27476 | 4.6% | -55552.7 &rarr; -55570.3 | 17.57 nats |
| theta+z, 40 it., n_starts=1 | 0.287906 &rarr; 0.260086 | **9.7%** | -55552.7 &rarr; -55581.76 | **29.01 nats** |
| theta+z, 150 it., n_starts=1 (fully converged) | 0.287906 &rarr; 0.259857 | 9.7% | -55552.7 &rarr; -55581.8 | 29.03 nats |
| theta+z, multistart (4 starts, 60 it.), best | 0.287906 &rarr; 0.259509 | **9.9%** | -55552.7 &rarr; -55583.00 | **30.26 nats** |

The 150-iteration run confirms the 40-iteration numbers are not an artifact of premature
stopping: a-opt reaches "Local minimum found" (gradient norm 7.6e-9, below
`OptimalityTolerance=1e-8`); d-opt halts because it cannot decrease further along the search
direction (gradient norm 2.4e-5), at essentially the same value as at iteration 40. The
4-start multistart check (60 it. each) lands within 0.3% of the same a-opt value and within
0.003% of the same d-opt value on every start — same conclusion as the theta-only case (the
even-spacing start is already close to the global optimum of this design space), just with
slightly more spread than the theta-only case's exact repeat convergence, consistent with the
larger (2M=16) search space having some near-degenerate configurations rather than one exact
minimum.

### Answering the question: "more considerable"?
**Partially yes.** Freeing `z_m` roughly **doubles** both gains (A-opt: 4.6%&rarr;~9.7-9.9%;
D-opt: 17.57&rarr;~29-30 nats) — a real, reproducible effect, confirming that a meaningful part of
the theta-only A-opt shortfall *was* a symmetry artifact of the overly restrictive
fixed-radius-fixed-height domain, as conjectured. But ~10% still falls short of the "considerable"
bar (>=15%) used elsewhere in this study, so the restrictive-domain hypothesis is only **partially**
confirmed, not fully.

**Diagnostic (see the new unrolled `(theta,z)` figure / `inspect_thetaz_results.m`-style
inspection):** sensors still cluster in azimuth toward the anomaly's side, exactly as in the
theta-only case. In height, they do **not** converge onto the anomaly's own height (`z=H/2`);
instead they **spread out** across roughly a fifth to a fifth-and-a-half of the full shell height
(`std(z_m)` ~= 0.36-0.40 m for a-opt, ~= 0.21 m for d-opt, out of a 1.75 m range). This is the same
decorrelation/de-clustering mechanism already identified for azimuth in the theta-only case-study
discussion, now also operating in the height direction: spreading in `z` reduces redundancy between
sensor rows (helps both criteria, more so D-optimality) rather than concentrating averaged
sensitivity near the anomaly (which is what would move the A-optimal trace substantially further).

**Conclusion or the whole story:** both hypotheses from the theory `.tex` closing remark hold in
part. The fixed-height restriction was indeed too restrictive and explains roughly half of the gap
between the criteria's gains; the remainder is consistent with the more fundamental
"A-optimality tolerates redundancy" explanation being real too, independent of how many geometric
degrees of freedom are on offer. This has been written up as a new subsection
"Follow-up: freeing the height $z_m$" (`\S`\ref{sec:case-study-thetaz}) in
`sensor_position_optimization_theory.tex`, immediately after the existing case-study discussion,
using the same notation (`\eqref{eq:chainrule1d}`, `Phi_A`/`Phi_D`, etc.). Compiled cleanly with 2
`pdflatex` passes (only pre-existing font/hyperref warnings, no undefined references); `.aux`/
`.log`/`.out` cleaned up afterward.

### Practical notes specific to this run
- Multistart random `eta` starts were drawn as `2*randn(n_sensors,1)` (not uniform — `eta` is an
  unconstrained real, not periodic like `theta`), giving a reasonably diverse spread of initial
  heights across the shell without extreme sigmoid saturation.
- `max_iterations` and `n_starts` in `example_anomaly_circle_thetaz.m` are now overridable from the
  workspace before calling the script (`if ~exist('max_iterations','var') ... end`, same pattern as
  `quick_test`), which is how the 40/150/multistart runs above were all driven from the one script
  without editing it between runs.
- A background MATLAB run survived a Windows sleep/wake cycle mid-optimization (the parallel pool
  briefly reported "shutting down" / "Job Queued" around the wake, then reconnected with 6 workers
  and resumed iterating normally) — worth knowing if a long run gets interrupted by a sleep, it is
  not necessarily lost, but it is not guaranteed either; check the log for a stall before assuming
  it is still working.

## NEXT PHASE (2026-07-06): human thorax geometry

The cylinder-geometry study is complete. The study continues on the human-thorax
geometry; the active handoff/plan for that phase is
`HumanThorax/HANDOFF_human_thorax.md` (Task A: azimuth-only optimization at fixed
radius; Task B: azimuth+radius with a position-dependent noise model to make the
radial optimum non-trivial). New agents should go there.
