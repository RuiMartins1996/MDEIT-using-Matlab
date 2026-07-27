# Handoff: the A-optimality ceiling, and the EIT (Hyvönen 2014) comparison

Written 2026-07-25. For a new Claude session that will **discuss and extend**
the finding that MDEIT sensor optimization yields a ~1% A-optimality
reduction while the EIT electrode-placement study yields ~70%.

Continues the handoff chain in this folder (`HANDOFF_tier3_data_space.md` →
`HANDOFF_2d_and_brute_force.md` → this one). Read those two first only if you
need the MDEIT code's history; **this document is self-contained for the
A-optimality question.**

---

## 0. TL;DR

1. A reference implementation of Hyvönen, Seppänen & Staboulis (2014),
   *Optimizing electrode positions in EIT*, was built from scratch in
   `../hyvonen2014_electrode_optimization/`. Fully validated (8 checks,
   including brute-force global-optimum agreement to −0.02%).
2. The MDEIT study's small A-optimality gain was diagnosed. **It is caused by
   the diagonal (white) prior, not by the physics, the design space, or the
   optimizer.** A physics-free bound shows no sensor arrangement could ever
   exceed a few percent with that prior; the optimizer already gets about
   half of that ceiling.
3. Swapping in a Gaussian smoothness prior on the *unchanged* MDEIT model
   takes the reduction from **1.1% → 55.1%**.
4. Two bugs were found in `example_anomaly_circle_2d.m` (one latent and
   silent — see §5).

Full write-up with all numbers: **`WHY_AOPT_REDUCTION_IS_SMALL.md`** (this
folder). Read it before the next discussion; this handoff assumes it.

---

## 1. The bound (the core result — verify it yourself before building on it)

Write `Γpr = LL'`, whitened Jacobian `J̃ = Γn^(-1/2) J L`, `J̃'J̃ = V S² V'`.
Then `Γ_* = L (I + J̃'J̃)⁻¹ L'` and

```
ρ := 1 − trace(Γ_*)/trace(Γpr) = [ Σᵢ (sᵢ²/(1+sᵢ²)) wᵢ ] / [ Σᵢ wᵢ ] ,   wᵢ = ‖L vᵢ‖²
```

`rank(J̃) ≤ n_data`, only modes with `sᵢ² ≳ 1` contribute, and each factor is
`< 1`. Maximizing `Σ wᵢ = Σ vᵢ'Γpr vᵢ` over orthonormal `k`-sets gives (Ky Fan)
the sum of the `k` largest eigenvalues of `Γpr`:

```
ρ  ≤  ( Σ_{i=1..d_eff} μᵢ ) / ( Σᵢ μᵢ )        μᵢ = eig(Γ_prior), descending
```

No forward model, no sensor positions, no physics enter. Computed in
`../hyvonen2014_electrode_optimization/analysis_aopt_reduction_ceiling.m`:

| | prior | effective rank¹ | ceiling | measured |
|---|---|---|---|---|
| MDEIT (this study) | **diagonal** | **485.6** | 2.3–3.6% | **1.1%** |
| Hyvönen EIT | dense smoothness λ=0.5 | **1.4** | 99.2% | 69.6% |

¹ participation ratio `(Σμ)²/Σμ²`.

**Loose end worth tightening:** the MDEIT ceiling above was computed at
`maxsz = radius/30` (7822 elements) while the measured 1.1% run used
`radius/20`. Both are single-digit-% so the conclusion is robust, but a
next session should recompute the ceiling at the *exact* run configuration
for a clean apples-to-apples number.

---

## 2. Why D-optimality escaped the same fate (important, and it resolves the
old "asymmetry" hypothesis)

`logdet(Γ_*) = logdet(Γpr) − Σᵢ log(1+sᵢ²)`, so the D-optimal information
gain is `½ Σᵢ log(1+sᵢ²)` — **a function of the `sᵢ` alone**. It does not
care whether the directions being contracted carry much prior *trace*.
A-optimality does: it weights each contracted direction by `wᵢ`.

That is the whole explanation for "A-opt gains little, D-opt gains a lot"
in the earlier MDEIT results. The previous hypothesis on record —
*"A-optimality tolerates the redundancy of a near-symmetric sensor ring"*
(see `HANDOFF_sensor_optimization_example.md` and the theory `.tex`
case-study section) — **should be considered superseded.** If you touch the
theory `.tex`, that section needs revising.

---

## 3. What reaching 55.1% required, and what it does *not* prove

`example_anomaly_circle_2d_smoothprior.m` + `sweep_smoothprior_configs.m`,
15 configs over 3 rounds. Best: 4 magnetometers, 2 electrodes, anomaly at
`0.75 R`, smoothness prior `λ = 1.6 R`, `d_target = 1` →
`φ_even = 10.473 → φ_opt = 4.702` (**55.1%**), ROI posterior-variance ratio
0.264, FD gradient check 1e-8…1e-10.

Two levers, and **be honest about the second**:

1. **Low effective prior rank.** Sets the ceiling. λ = 1× → 4× ROI radius
   took 20.4% → 39.0%.
2. **Noise level such that the even array is marginal** (`s₁² ≈ 1`). In the
   rank-1 limit `ρ = 1 − (1+s₁²_even)/(1+s₁²_opt)`, maximized near
   `s₁²_even ≈ 1`; at high SNR both arrays saturate and the gap closes. The
   winning config has SNR(even) = 0.24 dB, SNR(opt) = 9 dB.

So the honest claim is: **sensor placement matters most near the detection
threshold.** A 55% number quoted without its SNR context would be
misleading. If the real instrument sits comfortably above threshold, a small
A-opt gain is a true physical statement about the experiment.

**Negative results (all tested, don't re-run them blindly):**
- More sensors *hurts* this metric: 4→8→12 gave 32.7→24.0→22.6% despite the
  ceiling rising to 97.7%. The metric is optimized-vs-**even**, and a dense
  even ring already puts a sensor near the anomaly. Fewer also hurts
  (3 sensors: 28.1%). 4 is near a sweet spot for this geometry.
- Lowering `prior_std_background`: hurts (39.0 → 27.4%).
- Moving the anomaly wallward alone: no help (20.4 → 18.8%).

**Idea already checked and rejected — don't propose it again:** using
Hyvönen's weight matrix `A` in eq (3.3) (they set `A = I`; a ROI-selector `A`
would give "posterior variance in the ROI only") does **not** rescue a
diagonal prior. With `A` = ROI indicator the weighted prior is `0.16·I₃₉₈`,
all eigenvalues equal, so the ceiling is still `d_eff/398 ≈ 2.5%`. Reducing
the prior's effective rank (or massively increasing `d_eff`) is the only
lever.

---

## 4. Open questions worth discussing next session

1. **What prior does the actual MDEIT application warrant?** This is the
   crux and it is a *modelling* question, not a numerical one. The diagonal
   prior is not "wrong" so much as it encodes "every element independent",
   which is not what anyone believes about a conductivity field. What
   correlation length is defensible for the target application (phantom
   tank? thorax? brain?), and can it be justified from physiology/geometry
   rather than tuned to make the optimization look good?
2. **How large can `d_eff` actually get for MDEIT?** The other route to a
   high ceiling is many noise-resolvable data modes. Biot–Savart is a
   smoothing operator, so the singular values of `J` decay fast and `d_eff`
   may saturate well below `n_data`. Concrete experiment: plot the whitened
   spectrum of `J` vs. `n_sensors` (4, 8, 16, 32) and find where `#{sᵢ²>1}`
   stops growing. This tells you whether hardware (more magnetometers) can
   ever substitute for a better prior. **Nobody has measured this yet.**
3. **Is the mesh-dependence of the white-prior ceiling real?** Prediction:
   refining the mesh drives ρ → 0 for a diagonal prior and leaves it
   ~constant for a smoothness prior. Cheap to test (no FEM needed, just
   `eig(Γpr)` at several `maxsz`). If confirmed it is a clean argument that
   the white prior is not a discretization of any continuum prior, and worth
   a paragraph in the theory `.tex`.
4. **Should the study's headline metric change?** Optimized-vs-even is
   sensitive to how good "even" happens to be (hence "more sensors hurts").
   Candidates: reduction vs. the *worst* configuration; reduction vs. the
   brute-force global optimum; or the ROI posterior-variance ratio (already
   printed). Worth deciding deliberately rather than by inheritance.
5. **2D MDEIT has multiple local minima; EIT Case 1 does not.** `fminunc` on
   `brute_force_4x4_2d.m` converged to a meaningfully worse minimum than
   brute force (1.2% gap), whereas EIT Case 1 matched the polished global
   optimum to −0.02%. Is that a generic difference or an artefact of the two
   parameter choices? Would need matched-difficulty priors on both to say.
6. **EIT side, not started:** Case 3 (ellipse/peanut/complicated shapes,
   paper Fig. 6) and the Monte-Carlo MAP evaluation (Fig. 5, `N_draw=500`,
   paper gets a 0.75 MSE ratio). The geometry hooks exist
   (`boundary_curve.m` has only a `'disk'` branch; the general arc-length
   machinery is written but unvalidated).

---

## 5. Bugs found (one is silent — do not lose this)

1. **Latent and silent.** `costgrad_theta_z` in `example_anomaly_circle_2d.m`
   computes `Yd = Y .* d.'`, the diagonal shortcut for `Y*Γpr`. **Only valid
   for a diagonal prior.** The file already contains a `smooth_prior()`
   function and commented-out lines that would enable a dense prior — enable
   them and the *cost* stays correct while the *gradient is wrong*, with no
   error raised; `fminunc` converges somewhere meaningless. Fix:
   `Yd = Y*Gamma_prior` (done in the `_smoothprior` variant).
   Same assumption in the noise calibration:
   `svd(J0 .* sqrt(prior_variance)')` must become
   `svd(J0*chol(Gamma_prior,'lower'))`, and `d_target` is decisive (§3), so
   this is not cosmetic.
2. `imgr.elem_data = J_even(8,:)` hardcodes row 8; errors whenever
   `n_stim*n_sensors < 8` (2 or 3 sensors).

These are fixed **only in the `_smoothprior` variant**. The originals are
untouched — decide deliberately whether to port the fixes over.

---

## 6. File map

**This folder (MDEIT):**
- `WHY_AOPT_REDUCTION_IS_SMALL.md` — the full write-up. Start here.
- `example_anomaly_circle_2d.m` — **unmodified baseline.** Measured
  `φ_even = 1.706398e+01 → φ_opt = 1.688250e+01`, 1.1%, ROI ratio 0.987.
- `example_anomaly_circle_2d_smoothprior.m` — variant: smoothness prior +
  the three fixes. Sweep hooks `sw_*` (see its header).
- `sweep_smoothprior_configs.m` — sweep driver.
- `data/sweep_*.mat` — per-config summaries.

**EIT implementation (`../hyvonen2014_electrode_optimization/`):**
- `PLAN_implementation.md` — the spec that was implemented, incl. the two
  deliberate deviations from the paper (elementwise conductivity; fixed mesh
  with continuous-arc electrodes instead of per-iterate remeshing).
- `RESULTS.md` — all validation results and the cross-study comparison.
- `analysis_aopt_reduction_ceiling.m` — the bound in §1.
- `lib/` — 25 function files (CEM forward, adjoint Jacobian, exact discrete
  shape derivative, Woodbury OED criterion, Algorithm 1).
- `validate_forward_and_gradient.m` — V1–V4. **Run this first if you change
  anything in `lib/`.**

---

## 7. Practical gotchas

- **MATLAB R2024a only** (`"C:\Program Files\MATLAB\R2024a\bin\matlab.exe"`);
  R2025a/b lack the Optimization Toolbox. From Git Bash pass **native
  Windows** paths — `$(pwd)` gives an MSYS `/c/...` path MATLAB mangles.
- **Never loop over configs calling the MDEIT script in one workspace.** It
  shares the caller's workspace and reuses common names (it assigns `c` at
  the reconstruction step), silently corrupting the driver's loop state.
  `sweep_smoothprior_configs.m` spawns one MATLAB process per config for
  exactly this reason.
- **`data/` is gitignored** — nothing in it is recoverable from git. During
  this session the variant script initially inherited the original's `save`
  path and clobbered `data/example_anomaly_circle_2d_results.mat`; it was
  regenerated by re-running the unmodified script, and the variant now writes
  to its own filename. Check `save` paths when copying scripts.
- The dense smoothness prior is `n_elem × n_elem`. Keep `n_elem` in the low
  thousands (`maxsz = radius/20` → ~3000) or `eig`/`chol` get slow.
- Gradient checks are the repo's standard safety net and have caught every
  gradient bug so far. Target rel. err. 1e-6…1e-8. Run them before trusting
  any optimization result.
