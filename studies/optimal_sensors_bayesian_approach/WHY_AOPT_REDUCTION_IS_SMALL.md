# Why the MDEIT A-optimality reduction is small, and how to make it ~50%

Written 2026-07-25. Answers two questions raised when comparing this study
against `studies/hyvonen2014_electrode_optimization/` (EIT electrode
placement, Hyvönen–Seppänen–Staboulis 2014), which achieves a ~70%
A-optimality cost reduction where this study achieves ~2–5%.

**Short answer: the diagonal prior, not the physics and not the optimizer.**
With the prior currently used in `example_anomaly_circle_2d.m`, *no* sensor
arrangement whatsoever can do better than ~3%. The optimizer is already at
the information ceiling. Swapping in a smoothness prior lifts the ceiling to
~99% and a real run reaches **55.1%**.

---

## 1. The bound

Write `Γpr = LL'` and define the prior-whitened Jacobian `J̃ = Γn^(-1/2) J L`,
with `J̃'J̃ = V S² V'`. Then

```
Γ_*      = (J'Γn⁻¹J + Γpr⁻¹)⁻¹ = L (I + J̃'J̃)⁻¹ L'
trace(Γ_*)  = Σᵢ (1+sᵢ²)⁻¹ wᵢ ,     wᵢ := ‖L vᵢ‖²
trace(Γpr)  = Σᵢ wᵢ
```

so the fractional reduction is *exactly*

```
ρ = 1 − trace(Γ_*)/trace(Γpr) = [ Σᵢ (sᵢ²/(1+sᵢ²)) wᵢ ] / [ Σᵢ wᵢ ]        (*)
```

Two facts bound (*):

1. `rank(J̃) ≤ n_data`, so at most `n_data` of the `sᵢ` are nonzero — and
   only those with `sᵢ² ≳ 1` (signal above the noise floor) contribute
   meaningfully. Call that count `d_eff`; it is what the script's noise
   calibration sets via `d_target`.
2. Each surviving factor `sᵢ²/(1+sᵢ²) < 1`.

Maximizing `Σ wᵢ = Σ vᵢ'Γpr vᵢ` over orthonormal sets of size `k` gives, by
**Ky Fan's theorem**, the sum of the `k` largest eigenvalues of `Γpr`. Hence

```
ρ  ≤  ( Σ_{i=1..d_eff} μᵢ ) / ( Σᵢ μᵢ ) ,     μᵢ = eigenvalues of Γ_prior   (**)
```

**(\*\*) contains no forward model, no sensor positions, no physics.** It is a
hard information ceiling set by the prior alone.

## 2. Evaluating the ceiling (`../hyvonen2014_electrode_optimization/analysis_aopt_reduction_ceiling.m`)

| | prior | effective rank¹ | ceiling | achieved |
|---|---|---|---|---|
| **MDEIT** (this study) | **diagonal**, 7822 elements | **485.6** | **2.3%** (`d_eff=10`) / 3.6% (`n_data=16`) | ~2–5% |
| **Hyvönen EIT** | dense smoothness, λ=0.5 | **1.4** | **99.2%** | 69.6% |

¹ participation ratio `(Σμ)²/Σμ²`.

A diagonal prior says *every element is independently uncertain*. Its
eigenvalues are 398 identical copies of `prior_std_roi²` inside the ROI. With
only ~10 noise-resolvable data modes you can remove at most 10 of those 398
directions — about 2% of the prior trace — **regardless of where the
magnetometers go**. Hyvönen's smoothness prior with λ ≈ domain size collapses
to effective rank 1.4, so a *single* data mode already captures 83% of the
prior trace.

**Corollary worth noting:** with a diagonal prior the ceiling *shrinks as the
mesh is refined* (more elements ⇒ more prior variance in directions no
measurement can see). That mesh-dependence is a symptom that the white prior
is not the discretization of a well-posed continuum prior. A smoothness prior
is mesh-stable.

Same MDEIT mesh and ROI, prior swapped for a smoothness prior:

| λ (× tank radius) | ceiling (`d_eff=10`) |
|---|---|
| 0.10 | 87.1% |
| 0.25 | 93.8% |
| 0.50 | 98.8% |
| 1.00 | 99.9% |

## 3. Reaching ~50% for real

`example_anomaly_circle_2d_smoothprior.m` is a copy of
`example_anomaly_circle_2d.m` with the prior replaced by the Hyvönen
smoothness prior (eq 5.1–5.2). `sweep_smoothprior_configs.m` drives it over
configurations. All runs pass the finite-difference gradient check at 1e-8.

| round | config | eff. rank | ceiling | φ(even) | φ(opt) | **reduction** |
|---|---|---|---|---|---|---|
| — | *original, diagonal prior*² | 485.6 | 2.3–3.6% | 17.064 | 16.883 | **1.1%** |
| 1 | λ=1×ROI radius | 2.93 | 81.2% | 9.941 | 7.912 | 20.4% |
| 1 | 8 sensors | 1.74 | 95.1% | 4.143 | 3.149 | 24.0% |
| 1 | 12 sensors | 1.75 | 97.7% | 3.008 | 2.327 | 22.6% |
| 1 | λ=2×ROI radius | 1.74 | 87.8% | 6.616 | 4.456 | 32.7% |
| 2 | λ=4×ROI radius | 1.46 | 95.6% | 4.214 | 2.569 | 39.0% |
| 2 | λ=4, lower background std | 1.11 | 99.1% | 0.868 | 0.630 | 27.4% |
| 2 | 3 sensors | — | — | 0.989 | 0.711 | 28.1% |
| 3 | λ=8, `d_target`=4 | 1.38 | 99.1% | 2.591 | 1.824 | 29.6% |
| 3 | λ=16, `d_target`=4 | 1.36 | 99.8% | 2.333 | 1.755 | 24.8% |
| **3** | **λ=8×ROI radius, `d_target`=1** | **1.38** | **84.1%** | **10.473** | **4.702** | **55.1%** |

² Measured by re-running `example_anomaly_circle_2d.m` unmodified with its
currently committed parameters (2026-07-25): `phi_even = 1.706398e+01`,
`phi_opt = 1.688250e+01`, ROI posterior-variance ratio 0.987. This sits at
roughly half the ceiling — i.e. the optimizer is doing its job; there is
simply almost nothing available to win.

Best configuration (`r3_lam8_d1`): 2 electrodes, 4 magnetometers at
`rs = 1.05 R`, anomaly at `0.75 R`, smoothness prior with
`λ = 1.6 R`, noise calibrated at `d_target = 1`.
ROI posterior-variance ratio (optimized/even) = **0.264**.

### What actually drives it

Two things, and only two:

1. **Low effective prior rank** (the smoothness prior). Necessary — it sets
   the ceiling. Going from λ=1× to 4× ROI radius took 20.4% → 39.0%.
2. **Noise level such that the *even* array is marginal**, `s₁² ≈ 1`. In the
   rank-1 limit `φ ≈ μ₁/(1+s₁²) + bg`, so
   `ρ_opt-vs-even = 1 − (1+s₁²_even)/(1+s₁²_opt)`. If the SNR is *too high*,
   both arrays saturate (`φ → bg`) and the gap closes; if too low, neither
   sees anything. The maximum is near `s₁²_even ≈ 1`. That is exactly what
   `d_target=1` produces: SNR(even) = 0.24 dB, SNR(opt) = 9 dB — an 8.8 dB
   (≈7.5×) power gain purely from moving the sensors, which the rank-1
   formula turns into `10.47 → 4.70`.

**This is a real regime, not a trick, but it must be reported with its SNR
context:** sensor placement matters most when the measurement is near the
detection threshold. If your instrument is comfortably above threshold, the
A-optimality gain from placement is genuinely small — and that is a physical
statement about the experiment, not a defect of the optimizer.

### Things that did *not* help (all tested)

- **More sensors hurts** this metric: 4 → 8 → 12 sensors gave 32.7 → 24.0 →
  22.6%, even though the ceiling *rose* to 97.7%. The metric is
  optimized-vs-**even**, and a dense even ring already places a sensor near
  the anomaly, leaving little for optimization to win. (Same effect this
  study's own handoff noted for electrodes: "fewer electrodes helps.")
  Fewer than 4 also hurts (3 sensors: 28.1%) — 4 is near a sweet spot here.
- **Lowering the background prior std** hurts (39.0 → 27.4%). The background
  trace is not the binding constraint.
- **Moving the anomaly nearer the wall**, on its own, barely moved it
  (20.4 → 18.8%).

## 4. Latent bug this uncovered

`example_anomaly_circle_2d.m`'s `costgrad_theta_z` computes the A-optimality
gradient weight as

```matlab
Yd = Y .* d.';        % d = diag(Gamma_prior)
W  = Yd - (Yd*J.')*Y;
```

which is the diagonal shortcut for `Y*Γpr`. **It is only valid for a diagonal
prior.** The file already contains a `smooth_prior()` function and
commented-out lines that would enable a dense prior — if anyone uncomments
them, the *cost* stays correct but the *gradient is silently wrong*, and
`fminunc` would converge to the wrong place with no error raised. The fix is
`Yd = Y*Gamma_prior` (done in the `_smoothprior` variant).

The noise calibration has the same assumption:
`svd(J0 .* sqrt(prior_variance)')` must become `svd(J0*chol(Gamma_prior,'lower'))`,
otherwise `d_target` — which the section above shows is decisive — is
miscalibrated.

Also fixed in the variant: `imgr.elem_data = J_even(8,:)` (a figure line)
hardcodes row 8 and errors out whenever `n_stim*n_sensors < 8`, e.g. with 2
or 3 sensors.

## 5. Files

- `../hyvonen2014_electrode_optimization/analysis_aopt_reduction_ceiling.m` —
  derivation + ceiling evaluation for both studies (no FEM needed).
- `example_anomaly_circle_2d_smoothprior.m` — variant with the smoothness
  prior and the three fixes above. Saves to its own
  `data/example_anomaly_circle_2d_smoothprior_results.mat`.
- `sweep_smoothprior_configs.m` — sweep driver (one MATLAB process per
  config; the MDEIT script reuses caller variable names such as `c`, so a
  single-workspace loop corrupts itself).
- `data/sweep_*.mat` — per-config summaries.
