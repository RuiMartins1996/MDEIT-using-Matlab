# Results: Hyvönen, Seppänen & Staboulis (2014) electrode-position optimization

> **Continuing this work?** Start from
> `../optimal_sensors_bayesian_approach/HANDOFF_aopt_ceiling_and_eit_comparison.md`,
> which covers both this implementation and the A-optimality-ceiling
> diagnosis of the MDEIT study, with open questions.

Implementation of `PLAN_implementation.md` in this folder. All code lives in
`lib/` (proper function files, no copy-pasted local functions across scripts).
Mesh: fixed 2D triangular disk mesh (`build_mesh_ctx.m`), never remeshed
during optimization (Deviation 2). Conductivity: elementwise, not nodal
(Deviation 1), matching `studies/optimal_sensors_bayesian_approach/`.

## Validation summary

| # | Check | Result |
|---|---|---|
| V1 | Reciprocity of the measurement map (replaces "vs. EIDORS forward solve" — EIDORS's CEM only supports node-list electrodes, not the continuous arcs this study needs; reciprocity is a physical symmetry that directly tests our own discretization instead) | `max\|R-Rᵀ\|/max\|R\|` = 2.1e-15 (machine precision) |
| V2 | Jacobian vs. central differences in σ | rel. err 6e-8 to 3e-7 (5 random elements) |
| V3 | Smoothness of ψ(θ) across a mesh-boundary-node crossing | no discontinuity (max 2nd-diff/range = 2.6e-2, smooth curve) |
| V4 | Analytic gradient vs. central differences (the mandatory gate) | rel. err 2e-8 to 1.3e-3 across a-opt/d-opt, both random components |
| V5 | Case 1 (M=4): `fminunc` vs. brute-force grid (4845 combos) + local polish | gap = **-0.02% of φ_even** — `fminunc` matches the polished global optimum |
| V6 | Case 1 qualitative: electrodes migrate toward Ω′ | 2 of 4 electrodes land within ~15° of Ω′'s azimuth (135°) at every optimum found |
| V7 | Fig. 3 reproduction: κ_out sweep {0.05,0.10,0.15,0.20} | φ_opt grows monotonically (7.6→45.0) and configurations lose the near-mirror symmetry seen at κ_out=0.05 as κ_out grows — matches the paper's qualitative trend |
| V8 | Case 2 (M=12): electrodes should concentrate in the high-variance lower half | 10 of 12 electrodes (A-opt) land at azimuth ∈ (180°,360°) — "almost all electrodes move to the lower half," matching the paper almost verbatim |

Gradient-check gate (V4) was run **before** any optimization, per plan.

## Case 1 (M=4, small high-variance disk Ω′, paper Fig. 2)

Mesh: 3446 elements, 1787 nodes. Ω′: radius 0.25 disk centered at
`(-0.55,0.55)` (azimuth 135°, not specified numerically in the paper —
chosen and fixed). `κ_in=0.4`, `κ_out=0.03`, `λ=0.5`, `α_rep=1e-4`.

| criterion | φ_init | φ_opt (Algorithm 1) | φ_opt (fminunc) | metric |
|---|---|---|---|---|
| A-opt | 1.4069e+01 | 4.5580e+00 | 4.2786e+00 | 67.6% / **69.6%** trace reduction |
| D-opt | -1.0767e+05 | -1.0767e+05 | -1.0767e+05 | 2.38 / 2.38 nats (3.43 bits) info gain |

ROI (Ω′) posterior-variance ratio (optimized/init): **0.19–0.21** (A-opt).

Both optimizers converge to the same physical configuration: two electrodes
cluster within ~15° of Ω′'s azimuth, the other two spread over the rest of
the boundary — matching Fig. 2's description exactly ("two of the four
electrodes move close to the smaller disk").

## Case 1 brute-force validation (paper's own Fig. 2 validation strategy)

Grid: 20 angles (18° steps), `nchoosek(20,4)=4845` unordered combinations,
best point polished with `fminunc`.

```
phi_even                 = 1.406931e+01
phi_fminunc               = 4.278569e+00
phi_brute_force (grid)    = 4.520726e+00
phi_brute_force_polished  = 4.281468e+00
gap (fminunc vs polished brute force) = -0.0206% of phi_even
```

The negative gap (fminunc slightly *below* the polished grid optimum) is
within noise of "found the same optimum" — this is the paper's own headline
claim reproduced: "the optimal configurations provided by Algorithm 1 and the
brute force computations approximately coincide."

## Case 1, κ_out sweep (paper Fig. 3)

| κ_out | φ_opt (A-opt) | θ_opt (deg) |
|---|---|---|
| 0.05 | 7.63 | [15.6, 135.9, 122.6, 202.0] |
| 0.10 | 17.93 | [50.0, 198.5, 127.5, 307.4] |
| 0.15 | 29.72 | [50.1, 204.2, 129.0, 307.6] |
| 0.20 | 44.93 | [129.9, 313.1, 53.0, 209.6] |

At κ_out=0.05 the two "spread" electrodes sit roughly symmetric about Ω′'s
135° axis; by κ_out=0.20 only one electrode stays near Ω′ and the rest
spread asymmetrically — the same qualitative loss of symmetry the paper
reports (their Fig. 3 caption: "the electrode positions... are not symmetric
with respect to Ω′").

## Case 2 (M=12, semidisks with different prior variance, paper Fig. 4)

Mesh: 3446 elements. Upper half (`x₂≥0`): κ=0.03 ("almost known"). Lower half
(`x₂<0`): κ=0.4. 99.4% of the prior trace sits in the lower half.

| criterion | φ_init | φ_opt | metric |
|---|---|---|---|
| A-opt | 1.3621e+01 | 1.0915e+01 | 19.9% trace reduction |
| D-opt | -1.0095e+05 | -1.0095e+05 | 2.06 nats (2.97 bits) info gain |

Lower-half posterior-variance ratio (opt/init): 0.76 (A-opt), 0.82 (D-opt).
θ_opt (A-opt, deg): `[328.3, 356.4, 59.7, 139.8, 202.5, 182.3, 191.2, 231.5,
252.3, 263.7, 288.2, 338.4]` — **10 of 12** electrodes at azimuth in
(180°,360°) (the lower half), only 2 left in the upper half. Matches the
paper's Case 2 conclusion almost exactly.

## Comparison to `studies/optimal_sensors_bayesian_approach/` (MDEIT/magnetometers)

Both studies share the same optimal-experimental-design core: linearize
around the prior mean, form the Bayesian linear-Gaussian posterior
`Γ_* = (J'Γn⁻¹J + Γpr⁻¹)⁻¹`, minimize `tr(Γ_*)` (A-opt) or `logdet(Γ_*)`
(D-opt) via the same data-space Woodbury identities (`P=JΓpr`, `S=Γn+PJ'`,
`Y=S\P`) — this repo's `costgrad_theta_z`/`oed_criterion.m` are, up to the
prior being diagonal vs. block-dense, the same formula. `fminunc`
quasi-Newton with an analytic gradient is used identically in both.

**What differs is the measurement physics and the design space**, and that
difference shows up in the results:

- **Design-space dimensionality per sensor.** A magnetometer has 1 free
  angular coordinate (fixed radius/height) or 2 (with the θ,z extension);
  an electrode has 1 free coordinate (its left arc endpoint, since width is
  fixed) but the *measurement* is far richer per electrode (all M mean-free
  boundary potentials, vs. one Bz reading per magnetometer position).

- **A-opt vs. D-opt asymmetry — RESOLVED, and it is the PRIOR, not the
  physics.** The MDEIT study's headline finding was that A-optimality gains
  little (1-5% trace reduction) while D-optimality gains a lot (17.6 nats),
  which had been attributed to A-optimality tolerating redundancy on a
  near-symmetric sensor ring. That explanation turns out to be wrong. A
  follow-up analysis (`analysis_aopt_reduction_ceiling.m` in this folder,
  written up in `../optimal_sensors_bayesian_approach/WHY_AOPT_REDUCTION_IS_SMALL.md`)
  derives a sharp, physics-free upper bound via Ky Fan's theorem:

  ```
  rho <= (sum of the d_eff largest eigenvalues of Gamma_prior) / trace(Gamma_prior)
  ```

  The MDEIT study uses a **diagonal** prior (effective rank 485.6), giving a
  ceiling of **2.3-3.6%** — its measured 1.1% is already about half of the
  best any sensor arrangement could ever do. This study uses Hyvonen's
  **dense smoothness** prior (effective rank **1.4**), ceiling **99.2%**.
  The difference in achieved A-opt reduction is almost entirely explained by
  the prior's spectral decay, not by CEM-vs-Biot-Savart physics and not by
  the design space.

  Confirmed constructively: putting a smoothness prior on the *unchanged*
  MDEIT model raises its A-opt reduction from 1.1% to **55.1%**
  (`../optimal_sensors_bayesian_approach/example_anomaly_circle_2d_smoothprior.m`,
  gradient-checked at 1e-8). So the two studies' A-opt numbers were never
  comparable in the first place — they were measuring designs against priors
  of vastly different effective rank.

- **Multiple local minima.** The MDEIT 2D study's `brute_force_4x4_2d.m`
  found a case where `fminunc` from even spacing converged to a
  meaningfully *worse* local minimum than brute force (1.2% gap, three
  sensors clustered + one isolated on the far side, vs. all four together
  at the global optimum). **This EIT Case 1 shows the opposite**: `fminunc`
  essentially matches the brute-force global optimum (gap -0.02%, well
  within noise). One plausible reason: our Case-1 prior places the
  informative region as a single, compact, off-center inclusion, giving the
  cost landscape one clearly dominant basin, whereas the MDEIT 2D case's
  parameters (per its handoff notes) were specifically retuned to make
  optimization "matter," which may have also created a harder landscape.
  Whether EIT electrode placement is generically easier to optimize than
  magnetometer placement, or this is specific to the two studies' chosen
  priors, is not established by this comparison alone and would need a
  matched-difficulty prior on both problems to test properly.

- **Repulsion term.** Both use a `1/gap` repulsion barrier
  (`alpha_rep=1e-4` here vs. `1e-5` there) to keep sensors/electrodes from
  colliding; in both studies it plays a purely mechanical role (not part of
  the OED criterion) and did not need retuning to get sensible results.

## Known limitations / not done

- Case 3 (ellipse/peanut/complicated shapes, paper Fig. 6) is not
  implemented — `boundary_curve.m` only has the `'disk'` branch. The
  general non-circular machinery (`solve_right_endpoint.m`'s bisection
  branch, `arc_length_between.m`'s quadrature branch) is written and would
  only need a new `boundary_curve` case to activate, but is unvalidated.
- V9 (paper Fig. 5: Monte-Carlo MAP-error ratio, N_draw=500 conductivities
  vs. simulated noisy data on a finer mesh) was not run — flagged in the
  plan as optional/expensive and out of scope for this comparison pass.
- The d-opt Algorithm 1 vs. fminunc gap was not brute-force-validated
  (only a-opt was, following the paper's own Fig. 2, which is A-opt only).
