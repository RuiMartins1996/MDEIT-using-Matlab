# Handoff: Tier 3 — route the OED criterion through the data space (MDEIT)

Written 2026-07-24. For a new Claude session to implement **Tier 3** of the
performance work on
`studies/optimal_sensors_bayesian_approach/example_anomaly_circle.m`.

Tiers 1 (correctness), the z-channel-only reformulation, the objective-scaled
sensor repulsion, and **Tier 2** (precompute-once context) are already done and
validated (finite-difference gradient check, rel. err. ~1e-8; full pipeline runs
to completion). This document is *only* about Tier 3.

## Problem being fixed

`costgrad_theta_z` (the per-iterate cost/gradient) and
`compute_posterior_variance_diag` currently form the **dense `n_elem × n_elem`
posterior covariance** `Gamma_post = H^-1`:

```matlab
J_times_Gamma_prior = J * Gamma_prior;              % [nd x n]
S = Gamma_noise + J_times_Gamma_prior * J.';        % [nd x nd]  (small, SPD)
X = S \ J_times_Gamma_prior;                         % [nd x n]
Gamma_post = Gamma_prior - Gamma_prior * J.' * X;    % [n  x n]   <-- DENSE n x n
```

with `n = n_elem` (e.g. **17572** at `maxsz = radius/20`) and
`nd = n_data = n_stim*n_sensors` (**128**, z-only). That `Gamma_post` is
`17572^2 * 8 bytes ≈ 2.5 GB` and is allocated **every** function evaluation
(`trace(Gamma_post)` for the cost, and again inside `W`). It is now the dominant
cost — it dwarfs the FEM work Tier 2 already removed.

Everything the criterion needs factors through the **128-dimensional data space**.
Tier 3 never forms an `n × n` matrix: memory drops from `O(n^2)` to `O(nd·n)`,
flops from `O(nd·n^2)` to `O(nd^2·n)` (~130× fewer here), and the 2.5 GB
allocation disappears.

## The math (self-contained; verify before coding)

Notation: `J` is `[nd x n]`, `Gpr = Gamma_prior` (diagonal, prior variances
`d = diag(Gpr)`), `Gn = Gamma_noise` (here `σ²·I`, but keep it general).
`H = J' Gn^-1 J + Gpr^-1`, posterior covariance `H^-1`.

Woodbury:
```
H^-1 = Gpr - Gpr J' S^-1 J Gpr,   S = Gn + J Gpr J'   (nd x nd, SPD)
```

Define the two reusable data-space objects (both `[nd x n]`, computed once via a
single Cholesky of `S`):
```
P = J Gpr
Y = S^-1 P            (= S \ P)
```

Then everything is a contraction of `P` and `Y`:

| quantity | formula (no n×n) | shape |
|---|---|---|
| A-opt cost `φ_A = trace(H^-1)` | `sum(d) - sum(sum(P .* Y))` | scalar |
| posterior variance `diag(H^-1)` | `d - sum(P .* Y, 1).'` | `[n x 1]` |
| D-opt weight `W_D = Gn^-1 J H^-1` | `Y` | `[nd x n]` |
| A-opt weight `W_A = Gn^-1 J H^-2` | `Yd - (Yd J') Y`, with `Yd = Y .* d.'` | `[nd x n]` |

Derivation of the two weights (so the next session can check, not just trust):

- `Q := J H^-1 = J(Gpr - Gpr J' S^-1 J Gpr) = P - (J Gpr J') S^-1 P`.
  Since `J Gpr J' = S - Gn`, `Q = P - (S-Gn) S^-1 P = Gn S^-1 P = Gn·Y`.
  Hence `W_D = Gn^-1 J H^-1 = Gn^-1 Q = Y`.
- `W_A = Gn^-1 J H^-2 = (Gn^-1 Q) H^-1 = Y H^-1 = Y(Gpr - Gpr J' S^-1 J Gpr)`
  `= Y Gpr - (Y Gpr J') S^-1 P = Yd - (Yd J') Y`, using `P' = Gpr J'`,
  `Y Gpr = Yd`, and `S^-1 P = Y`.

These are algebraically identical to the current
`W = Gamma_noise\(J*Gamma_post*Gamma_post)` (A) and
`Gamma_noise\(J*Gamma_post)` (D) — so **`grad_z`, `calc_Jz`, the chain rule and
the repulsion terms do not change at all**. Only the block that produces the
scalar cost and the `[nd x n]` weight `W` changes.

Factor `S` once and reuse: `Ls = chol(S,'lower'); Y = Ls.'\(Ls\P);`.

## Exact edits

### 1. `costgrad_theta_z` — replace the cost block

Find:
```matlab
% Posterior covariance via the Woodbury/data-space form (S is small):
%   Gamma_post = H^-1,  H = J'*inv(Gn)*J + inv(Gpr)
J_times_Gamma_prior = J * Gamma_prior;
S = Gamma_noise + J_times_Gamma_prior * J.';
X = S \ J_times_Gamma_prior;
Gamma_post = Gamma_prior - Gamma_prior * J.' * X;

switch opt_mode
    case 'a-opt'
        phi_data = trace(Gamma_post);        % trace(H^-1)
    case 'd-opt'
        phi_data = nan;                      % (d-opt cost still open, see 1.1)
end
```
Replace with:
```matlab
% Data-space (Woodbury) form: never forms the n x n posterior covariance.
%   H^-1 = Gpr - Gpr J' S^-1 J Gpr,   S = Gn + J Gpr J'   (nd x nd, SPD)
P  = J * Gamma_prior;                  % [nd x n]
S  = Gamma_noise + P * J.';            % [nd x nd]
Ls = chol(S,'lower');                  % factor once (reused for Y and logdet)
Y  = Ls.' \ (Ls \ P);                  % S^-1 P  [nd x n]
d  = diag(Gamma_prior);                % prior variances [n x 1]

switch opt_mode
    case 'a-opt'
        phi_data = sum(d) - sum(sum(P .* Y));          % trace(H^-1)
    case 'd-opt'
        % -logdet(H) = logdet(Gn) + logdet(Gpr) - logdet(S).  The first two
        % terms are theta-independent constants; keep them so the reported
        % "nats" stay comparable across configs (and with the old numbers).
        logdet_S = 2*sum(log(diag(Ls)));
        phi_data = sum(log(diag(Gamma_noise))) ...
                 - sum(log(d)) - logdet_S;
end
```
Notes:
- `sum(sum(P.*Y))` = `trace(P' Y)` = `trace(Gpr J' S^-1 J Gpr)` — the correction
  term subtracted from `trace(Gpr)`.
- The d-opt cost is the **optional** 1.1 fix, essentially free here because `Ls`
  is already computed. If you want to keep this PR strictly Tier 3, set
  `phi_data = nan` for d-opt instead and leave 1.1 for later — **confirm with the
  user** which they want, since re-enabling d-opt changes `opt_modes` at the top
  of the script.

### 2. `costgrad_theta_z` — replace the gradient weight block

Find:
```matlab
% ---- Gradient ----
% dphi = -2*<W, dJ/dp>_F  with  W_A = Gn^-1*J*H^-2,  W_D = Gn^-1*J*H^-1
switch opt_mode
    case 'a-opt'
        W = Gamma_noise\(J*Gamma_post*Gamma_post);   % Gn^-1*J*H^-2
    case 'd-opt'
        W = Gamma_noise\(J*Gamma_post);              % Gn^-1*J*H^-1
end
```
Replace with:
```matlab
% ---- Gradient ----
% dphi = -2*<W, dJ/dp>_F  with  W_A = Gn^-1*J*H^-2,  W_D = Gn^-1*J*H^-1.
% Both come from the same data-space objects P, Y (no n x n matrix).
switch opt_mode
    case 'a-opt'
        Yd = Y .* d.';               % Y*Gpr   [nd x n]
        W  = Yd - (Yd*J.')*Y;        % Gn^-1*J*H^-2
    case 'd-opt'
        W  = Y;                      % Gn^-1*J*H^-1
end
```
`P`, `Y`, `d` from the cost block are in scope (same function). Nothing else in
`costgrad_theta_z` changes — the repulsion scaling (`phi_data*(1+alpha_rep*Grep)`),
`grad_z`, and the cylindrical chain rule are untouched.

### 3. `compute_posterior_variance_diag` — drop the n×n Cholesky

Change the signature from `(...,inv_Gamma_prior,inv_Gamma_noise)` to
`(...,Gamma_prior,Gamma_noise)` and replace the body. Find:
```matlab
function post_var = compute_posterior_variance_diag(ctx,theta,rs,zs,...
    inv_Gamma_prior,inv_Gamma_noise)

n_elem = ctx.n_elem;

J = calc_Jz(ctx,theta_to_locations(theta,rs,zs));

H = full(J.'*inv_Gamma_noise*J + inv_Gamma_prior);
L = chol(H,'lower');
X = L\eye(n_elem);                    % inv(L)
post_var = sum(X.^2,1)';              % diag(H^-1) = diag(inv(L)'*inv(L))
end
```
Replace with:
```matlab
function post_var = compute_posterior_variance_diag(ctx,theta,rs,zs,...
    Gamma_prior,Gamma_noise)

J = calc_Jz(ctx,theta_to_locations(theta,rs,zs));

P = J*Gamma_prior;                          % [nd x n]
S = Gamma_noise + P*J.';                     % [nd x nd]
Y = S \ P;                                   % S^-1 P
post_var = diag(Gamma_prior) - sum(P.*Y,1)'; % diag(H^-1)
end
```

### 4. Update the two call sites (body, "Posterior variance diagnostics")

Find:
```matlab
post_var_even = compute_posterior_variance_diag(ctx,theta_even,rs,zs,...
    inv_Gamma_prior,inv_Gamma_noise);
post_var_aopt = compute_posterior_variance_diag(ctx,theta_opt{1},rs,zs,...
    inv_Gamma_prior,inv_Gamma_noise);
```
Replace with:
```matlab
post_var_even = compute_posterior_variance_diag(ctx,theta_even,rs,zs,...
    Gamma_prior,Gamma_noise);
post_var_aopt = compute_posterior_variance_diag(ctx,theta_opt{1},rs,zs,...
    Gamma_prior,Gamma_noise);
```
After this, `inv_Gamma_prior` / `inv_Gamma_noise` may become unused in the body.
`inv_Gamma_prior` is still built at "Prior covariance"; leave it or remove it
(cosmetic). Do **not** remove `Gamma_prior` / `Gamma_noise`.

## Validation (do all three; R2024a required — see below)

1. **Regression against the current (Tier 2) code — the key check.** Before
   editing, run the current script and record `phi(even spacing)` printed under
   `a-opt` (with the committed params: `roi_radius_factor=2.0`,
   `prior_std_background=0.05`, `d_target=75`, `alpha_rep=1e-5`,
   `maxsz=radius/20`). After Tier 3, `phi_even` must match to ~1e-10 relative —
   it is the *same* math, only reorganized. A mismatch beyond roundoff means a
   sign/transpose slip in `W_A` or the cost.
2. **Finite-difference gradient check.** `quick_test = true` runs
   `check_gradient_fd` (3 random components, central differences). Expect
   rel. err. ~1e-8, same as now. This catches any error in `W_A` (the gradient
   uses `W`, the cost does not, so #1 and #2 together pin both down).
3. **Posterior-variance diagnostic** printed at the end ("Posterior variance
   sums (ROI | background)") must match the Tier 2 numbers to ~1e-8.

Run:
```
"C:\Program Files\MATLAB\R2024a\bin\matlab.exe" -batch ^
  "cd('C:\Users\RuiMartins\Desktop\MDEIT-using-Matlab-main\studies\optimal_sensors_bayesian_approach'); quick_test=true; example_anomaly_circle"
```
(From Git Bash, pass the path as a **native Windows** string — `$(pwd)` yields a
`/c/...` MSYS path that MATLAB mangles into `C:\c\...`.)

Optional: `checkcode('example_anomaly_circle.m')` for lint. `parpool` is no
longer used (Tier 2 removed all `parfor`), so no worker startup.

## Practical notes / expected outcome

- **Numbers must not change.** Tier 3 is a pure reorganization of the linear
  algebra; the optimizer trajectory, `phi_even`, `phi_opt`, and all figures
  should be identical to the Tier 2 run within roundoff. If the *results* change,
  something is wrong — this is not a modeling change.
- **Speed/memory:** the 2.5 GB `Gamma_post` allocation per evaluation is gone;
  per-eval linear-algebra flops drop ~130×. At `maxsz=radius/20` the FEM solves
  (Tier 2 already made them one factorization + block back-substitution) and the
  Biot–Savart `R`/`dR` assembly become the remaining costs. If more speed is
  wanted afterwards, the next lever is **Tier 3.2**: fuse `grad_z` so it never
  materializes intermediate `[num_stim x n_elem]` blocks (accumulate the
  `sum(W.*dJ)` inner products directly) — but profile first; it may already be
  dominated by `R`/`dR`.
- **Numerical robustness:** `S` is SPD (`Gn` positive-definite + `J Gpr J'` PSD),
  so `chol(S)` is safe. If a pathological config ever makes `chol` fail, add a
  tiny jitter (`S = S + 1e-12*trace(S)/nd*eye(nd)`) or fall back to `Y = S\P`
  with `decomposition(S,'ldl')`; not expected to be needed here.
- The `d-opt` cost formula (step 1) is the one still-open Tier 1 item (1.1). It is
  cheap to include here but re-enables a second optimization mode; **ask the user**
  whether to bundle it or ship Tier 3 a-opt-only first.
- Keep `sensor_position_optimization_theory.tex` unchanged — Tier 3 introduces no
  new methodology, only a cheaper evaluation of the existing criterion.
