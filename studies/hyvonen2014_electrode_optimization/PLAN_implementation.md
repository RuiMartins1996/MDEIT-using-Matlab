# Implementation plan: Hyvönen, Seppänen & Staboulis (2014), *Optimizing electrode positions in EIT*

Target: a **Claude Sonnet session** implementing the paper's Algorithm 1 in this repo, in
MATLAB, so the result can be compared against the magnetometer-configuration
optimization already built in `studies/optimal_sensors_bayesian_approach/`.

Paper: `hyvonen_2014_optimizing_electrode_positions_in_electrical_impedance_tomography.pdf`
(repo root). SIAM J. Appl. Math. 74(6), 1831–1851. Read §3.1, §4 and §5 before coding;
§2 (Fréchet-derivative theory) is background you do **not** need to reproduce — see
"Deviation 2" below for why.

Deliverable folder: `studies/hyvonen2014_electrode_optimization/` (this file's folder).
Do **not** modify anything in `studies/optimal_sensors_bayesian_approach/` — that study is
the comparison baseline and must stay as-is.

---

## 0. What the paper actually does (condensed, with equation numbers)

EIT with the **complete electrode model** (CEM), M finite-width electrodes on the boundary
of a 2D star-shaped domain Ω. The electrodes' **angular positions** are the design variables.

1. **Bayesian linearized posterior.** Prior `σ ~ N(σ_*, Γ_pr)`, additive noise
   `ε ~ N(0, Γ_noise)`. Linearize the current-to-voltage map about the prior mean σ_*
   (3.8). The posterior covariance is then explicit and **independent of the data** (3.10):

   ```
   Γ_*(θ) = ( J(θ)' Γ_noise^-1 J(θ) + Γ_pr^-1 )^-1 ,    J = ∂_σ R[σ_*]
   ```

2. **Design criteria** (4.3), including a repulsion term that stops electrodes colliding:

   ```
   ψ(θ) = α · Σ_{m=1..M} 1/g_m(θ)  +  { tr(Γ_*(θ))      (A-optimality, from 3.11)
                                       { log det(Γ_*(θ)) (D-optimality, from 3.12)
   ```
   with `g_m` = arc-length gap between electrode m and m+1 (mod M), and `α = 1e-4`.
   Minimizing `log det Γ_*` == maximizing the prior→posterior information gain.

3. **Algorithm 1** — normalized steepest descent with a line search:
   ```
   θ_new = θ - t_min · ∇ψ(θ)/|∇ψ(θ)| ,   t_min = argmin_{t ∈ △} ψ(θ - t ∇ψ/|∇ψ|)
   ```
   where `△ ⊂ [0,∞)` is chosen so all gaps stay positive. Repeat to convergence.

4. **Gradient** `∂ψ/∂θ` requires `∂J/∂θ` — a *second* derivative of the forward map (once
   w.r.t. conductivity, once w.r.t. electrode position). The paper computes it from the
   continuum Fréchet formula (2.21)/(4.1). We compute the exact **discrete** analogue
   instead (§4 below) — same quantity, cheaper, and it passes finite-difference checks to
   machine precision.

5. **Numerics** (§4, §5): fixed electrode widths, so `θ⁺_m = θ⁺_m(θ⁻_m)` via the
   arc-length relation (4.2); `log det` evaluated through a Cholesky factor to avoid
   overflow; current patterns `I^(j) = e_1 − e_{j+1}`, j = 1..M−1 (one fixed feeding
   electrode); prior mean `σ_* ≡ 1`; contact impedances `z_m = 1`;
   `Γ_noise = (1e-3 · max_{k,l}|U_k(σ_*,θ_init) − U_l(σ_*,θ_init)|)² · 𝟙`, fixed at the
   *initial* configuration; Gaussian smoothness prior (5.1).

---

## 1. Two deliberate deviations from the paper (read this before you start)

### Deviation 1 — element-wise conductivity, no separate background mesh
The paper carries the conductivity on a *fixed* uniform triangulation of a background
domain `D ⊃ Ω`, with a projection matrix `P` onto the forward mesh, **because they remesh
Ω for every electrode configuration** (§4, "we are forced to carry out several
technicalities"). We do not remesh (Deviation 2), so:

- `D = Ω`, `P = I`, and there is **no** projection machinery.
- Conductivity is **piecewise constant per element** (`n_elem` unknowns), not nodal.
  This matches `studies/optimal_sensors_bayesian_approach/` (which is elementwise), which
  is what makes the two studies comparable. Evaluate the prior (5.1) at **element
  centroids** rather than node coordinates.

### Deviation 2 — fixed mesh + continuous-arc electrodes, instead of remeshing
The paper rebuilds the FE mesh at every evaluation of ψ so that electrode edges land on
mesh edges. That is slow and, worse, makes ψ(θ) **non-smooth** (mesh noise), which would
wreck both the line search and any finite-difference validation.

Instead: **one fixed triangular mesh of Ω**, and electrodes defined by *continuous* arc
endpoints `[θ⁻_m, θ⁺_m]` that are allowed to cut boundary edges partway. The CEM electrode
boundary integrals

```
(1/z_m) ∫_{E_m} (u − U_m)(v − V_m) ds
```

are then evaluated by 1-D Gauss quadrature over the **overlap** of each boundary edge with
the arc. This is a legitimate discretization of the same continuum CEM (2.2)/(2.4); it
makes ψ(θ) smooth (C¹ in θ) and makes the exact discrete gradient available. It also means
you never call Netgen inside the optimization loop.

State both deviations in the header comment of the main script.

---

## 2. Repository context you need

- **EIDORS** lives at `eidors-v3.12-ng/eidors` and is on the path via the study scripts'
  `prepare_workspace` / `addpath(genpath(<repo>/functions))` preamble. Copy the preamble
  from `studies/optimal_sensors_bayesian_approach/example_anomaly_circle_2d.m`
  (lines ~40–55).
- **Use EIDORS only for meshing and cross-validation**, not inside the optimization loop:
  - `ng_mk_cyl_models([0, 1, maxsz], 0, [])` → a 2D unit-disk triangular mesh (`fmdl.nodes`
    is Nx2, `fmdl.elems` is Nx3). Electrode-free; we attach our own arc electrodes.
  - `system_mat_1st_order` / `fwd_solve` / `calc_jacobian` → reference implementations to
    validate our CEM assembler and Jacobian against (validation step V1/V2 in §7).
- **MATLAB R2024a is mandatory**: only that install has the Optimization Toolbox.
  ```
  "C:\Program Files\MATLAB\R2024a\bin\matlab.exe" -batch "cd('C:\Users\RuiMartins\Desktop\MDEIT-using-Matlab-main\studies\hyvonen2014_electrode_optimization'); quick_test=true; hyvonen_case1"
  ```
  (From Git Bash pass a **native Windows** path — `$(pwd)` gives an MSYS `/c/...` path that
  MATLAB mangles.)
- **Code organization — do this differently from the existing study.** The MDEIT scripts
  duplicate their local functions across `example_anomaly_circle_2d.m` and
  `brute_force_4x4_2d.m` (MATLAB scripts can't export local functions), which the previous
  handoff flags as a maintenance problem. Here, put every reusable routine in
  `studies/hyvonen2014_electrode_optimization/lib/` as **proper function files**, and have
  each script `addpath(fullfile(script_folder,'lib'))`. No copy-paste duplication.

---

## 3. The forward model to implement (`lib/`)

### 3.1 Geometry: star-shaped boundary `γ: [0,2π) → R²`
```matlab
function [g, gdot] = boundary_curve(theta, shape)
```
- `shape.type = 'disk'` → `γ(θ) = [cos θ, sin θ]`, `|γ̇| = 1` (Cases 1 and 2).
- `shape.type = 'ellipse' | 'peanut'` → Case 3, optional (Phase 6).
Return both the point and the derivative; `|γ̇|` is needed for arc lengths and for the
gradient chain rule.

### 3.2 Electrodes
Design vector `θ⁻ = [θ⁻_1 … θ⁻_M]' ∈ R^M` (left endpoints, kept in increasing order).
Widths are **fixed in arc length** `w`; the right endpoint solves

```
∫_{θ⁻_m}^{θ⁺_m} |γ̇(t)| dt = w        ⇒   dθ⁺_m/dθ⁻_m = |γ̇(θ⁻_m)| / |γ̇(θ⁺_m)|      (4.2)
```

For the disk this is just `θ⁺ = θ⁻ + w` and `dθ⁺/dθ⁻ = 1`; implement the general form
anyway (`solve_right_endpoint.m`) so Case 3 works later.

Gaps: `g_m = arclength(θ⁺_m → θ⁻_{m+1})`, cyclically, with `g_M` wrapping through 2π.

### 3.3 CEM system matrix with partial-edge electrode integrals
```matlab
function [A, dA] = cem_system_matrix(mesh, sigma, thetaMinus, thetaPlus, z_contact, shape)
```
Block structure (standard CEM FEM, P1 triangles), unknowns `x = [u_nodes; U_electrodes]`:

```
A = [ K(σ) + Σ_m (1/z_m) B_m        −Σ_m (1/z_m) c_m
      −Σ_m (1/z_m) c_m'              diag(|E_m|/z_m) ]
```
- `K(σ)` = usual P1 stiffness, `K = Σ_k σ_k K_k`, with `K_k` the per-element stiffness
  contribution (store these — you need them for the Jacobian).
- `B_m[i,j] = ∫_{E_m} φ_i φ_j ds`, `c_m[i] = ∫_{E_m} φ_i ds`, `|E_m| = ∫_{E_m} ds`.

**Partial-edge quadrature (the crux).** For each boundary edge `e` with endpoint
parameters `[t_a, t_b]`, intersect `[t_a,t_b]` with `[θ⁻_m, θ⁺_m]` (handle 2π wraparound).
On the overlap `[s_0, s_1]` use a Gauss–Legendre rule (5–7 points is plenty for P1) in the
parameter `t`, with `ds = |γ̇(t)| dt` and the P1 shape functions evaluated at the local
edge coordinate of `γ(t)` (linear in arc position along the straight edge; on a polygonal
mesh approximate `γ(t)` by its projection onto the edge). Accumulate into `B_m`, `c_m`,
`|E_m|`.

**Derivatives w.r.t. the endpoints are free.** Because the only θ-dependence is the
integration limit, by the fundamental theorem of calculus

```
∂B_m/∂θ⁺_m = +|γ̇(θ⁺_m)| · φ(θ⁺_m) φ(θ⁺_m)' ,   ∂B_m/∂θ⁻_m = −|γ̇(θ⁻_m)| · φ(θ⁻_m) φ(θ⁻_m)'
```
and analogously for `c_m` (`φ(·)` instead of the outer product) and `|E_m|` (`±|γ̇|`).
These are **rank-1/point evaluations at the electrode ends** — exactly the discrete
counterpart of the paper's (2.20)/(4.1) ("the integrals over the boundaries of the
electrodes in (2.21) are reduced to point evaluations at the electrode edges"). Return
them as `dA{m}` (sparse, a handful of nonzeros each), already combined via (4.2):

```
dA_m := ∂A/∂θ⁻_m|_direct  +  (dθ⁺_m/dθ⁻_m) · ∂A/∂θ⁺_m
```

Note `A = K(σ) + C(θ)`: the stiffness block depends on σ only, the electrode blocks on θ
only. **The mixed second derivative `∂²A/∂σ∂θ` is therefore zero**, which is what makes §4
clean.

### 3.4 Grounding
The CEM potential is defined up to an additive constant (the paper works in `C^M/C`).
Enforce `Σ_m U_m = 0` — either by a Lagrange-multiplier row/column or by projecting the
measurement operator onto the mean-free subspace. Pick one and use it consistently in the
forward solve, the Jacobian and the derivative; a sloppy ground is the #1 source of
"gradient check fails by exactly a constant" bugs here.

### 3.5 Currents and measurements
- Patterns (§5): `I^(j) = e_1 − e_{j+1}`, `j = 1..N` with `N = M−1`.
- Measurement: **all M electrode potentials per pattern** (mean-free), so
  `data dim n_data = M·N = M(M−1)`. For Case 1 (M=4) that is 12; for Case 2 (M=12), 132.
  Tiny — exploit this (§5).

### 3.6 Jacobian `J = ∂U/∂σ` (n_data × n_elem)
Standard adjoint: with `A x^{(j)} = b^{(j)}` and measurement selector `m_i`
(`U_i` extraction, mean-free), and `y_i = A^{-1} m_i`,

```
J_{(j,i),k} = −y_i' [K_k ⊕ 0] x^{(j)}
```
One factorization of `A` (`decomposition(A,'ldl')` or `chol` — A is symmetric) serves all
forward and adjoint solves.

---

## 4. The gradient of ψ — the mathematical core

### 4.1 Reduce to `⟨W, ∂J/∂θ⟩` (identical in form to the MDEIT study)
Differentiating (3.10) and using `∂Γ_*/∂θ = −Γ_* (∂J' Γn^-1 J + J' Γn^-1 ∂J) Γ_*`:

```
∂/∂θ tr(Γ_*)     = −2 ⟨ Γn^-1 J Γ_*² , ∂J/∂θ ⟩_F      →  W_A = Γn^-1 J Γ_*²
∂/∂θ logdet(Γ_*) = −2 ⟨ Γn^-1 J Γ_*  , ∂J/∂θ ⟩_F      →  W_D = Γn^-1 J Γ_*
```

This is **exactly** the `dphi = -2*<W, dJ/dp>_F` structure already used in
`example_anomaly_circle_2d.m` (`costgrad_theta_z`). Reuse the data-space Woodbury forms
derived there (`HANDOFF_tier3_data_space.md` has the full derivation — read it, it applies
verbatim with `J` = EIT Jacobian):

```
P = J Γpr          [n_data × n_elem]
S = Γn + P J'      [n_data × n_data]  SPD, small
Ls = chol(S,'lower');  Y = Ls' \ (Ls \ P)

tr(Γ_*)      = sum(diag(Γpr)) − sum(sum(P .* Y))
logdet(Γ_*)  = logdet(Γpr) + logdet(Γn) − 2*sum(log(diag(Ls)))
W_D = Y
W_A = Yd − (Yd J') Y ,  Yd = Y * Γpr        (diagonal Γpr: Yd = Y .* d')
```
**Caveat vs. the MDEIT code:** there `Γpr` is diagonal, so `Yd = Y .* d'`. Here the paper's
prior (5.1) is **dense** (a Gaussian smoothness kernel), so use the matrix product
`Yd = Y * Gamma_prior` and `sum(d)` → `trace(Gamma_prior)`. Write `lib/oed_criterion.m`
to accept a full (or block-diagonal) `Gamma_prior`.

This also removes any need for `Γ_pr^-1` and for the `n_elem × n_elem` posterior — the
paper's Cholesky-of-`Γ_*` overflow workaround (§4) becomes unnecessary because `logdet` now
goes through the `n_data × n_data` Cholesky of `S`.

### 4.2 `∂J/∂θ` contracted with W — the exact discrete shape derivative
With `A(σ,θ) x^{(j)} = b^{(j)}`, `A = K(σ) + C(θ)`, `A_{σ_k} = K_k`, `A_{θ_m} = dA_m`, and
`∂²A/∂σ∂θ = 0`:

```
x_k        = −A^-1 K_k x                       (σ-derivative of the state)
x_θ        = −A^-1 dA_m x                      (θ-derivative of the state)
x_{kθ}     = −A^-1 ( dA_m x_k + K_k x_θ )      (mixed)
∂J_{(j,i),k}/∂θ_m = m_i' x^{(j)}_{kθ} = −y_i'( dA_m x^{(j)}_k + K_k x^{(j)}_θ )
```

Now contract with `W_{(j,i),k}` **before** ever forming `∂J/∂θ` (this is what keeps the
cost independent of `n_elem` and of the number of design variables):

```
Term 1:  Σ_k W_{(j,i),k} x^{(j)}_k  =  z_{ji} := −A^-1 K(w_{ji}) x^{(j)} ,
         where w_{ji} := W_{(j,i),:}' is treated as an "element-data vector" and
         K(w) = Σ_k w_k K_k is assembled exactly like a stiffness matrix.
         contribution:  −Σ_{j,i} y_i' dA_m z_{ji}

Term 2:  Σ_{j,i} y_i' K(w_{ji}) x^{(j)}_θ
       = Σ_j q_j' x^{(j)}_θ  with  q_j := Σ_i K(w_{ji}) y_i
       = −Σ_j (A^-1 q_j)' dA_m x^{(j)}
         contribution:  +Σ_j (A^-1 q_j)' dA_m x^{(j)}   (mind the two minus signs)

∂ψ_data/∂θ_m = −2 · [ Term1 + Term2 ]
```

**Solve count per gradient evaluation** (all with the *same* factorization of `A`):
`N` forward + `M` adjoints `y_i` + `N·M` solves for `z_{ji}` + `N` solves for `A^-1 q_j`.
For M=12 that is ~150 back-substitutions of a small sparse system — cheap. Crucially, the
`dA_m` contractions at the end are O(1) each (a few nonzeros), so adding design variables
costs almost nothing.

### 4.3 Repulsion term
`ψ_rep = α Σ_m 1/g_m`. With `g_m = arclength(θ⁺_m → θ⁻_{m+1})`,
`∂g_m/∂θ⁻_{m+1} = |γ̇(θ⁻_{m+1})|` and `∂g_m/∂θ⁻_m = −|γ̇(θ⁺_m)|·dθ⁺_m/dθ⁻_m`, so

```
∂ψ_rep/∂θ⁻_m = α [ |γ̇(θ⁺_m)| (dθ⁺_m/dθ⁻_m) / g_m²  −  |γ̇(θ⁻_m)| / g_{m-1}² ]
```
Derive this yourself and confirm against FD — don't trust the sign above blindly.

---

## 5. Model configuration (paper §5) — put these in one `lib/hyvonen_config.m`

| quantity | value |
|---|---|
| prior mean | `σ_* ≡ 1` |
| contact impedances | `z_m = 1`, all m |
| currents | `I^(j) = e_1 − e_{j+1}`, `j=1..M−1` |
| electrode width | `π/16` (paper's motivating example) unless stated otherwise; keep configurable |
| repulsion | `α = 1e-4` |
| prior (5.1) | `Γ_ij = κ_ij² exp(−|x_i − x_j|² / (2λ²))`, `x` = element centroids |
| noise | `Γ_noise = (1e-3 · max_{k,l}|U_k(σ_*,θ_init) − U_l(σ_*,θ_init)|)² · 𝟙`, computed once at θ_init and **held fixed** |

**Case 1** (validation case, M = 4): Ω = unit disk; `Ω' ⊂ Ω` a smaller disk (see Fig. 2 —
radius ≈ 0.25 centred near `(−0.55, 0.55)`; exact placement is not given, so pick one, note
it in the header, and keep it fixed across all Case-1 runs). `Γ_pr = Γ(λ=0.5, κ)` with (5.2):
```
κ_ij = 0.4   if both centroids ∈ Ω'
       0.03  if both ∉ Ω'
       0     if one in, one out          → Γ_pr is BLOCK DIAGONAL (2 blocks)
```
**Case 2** (M = 12): unit disk, `Γ_pr = Γ(0.5, κ)` with (5.3): `κ = 0.03` if both
centroids have `x₂ ≥ 0`, `0.4` if both `x₂ < 0`, `0` across — again block diagonal.

**Case 3** (optional): ellipse / peanut / the complicated shape of Fig. 6, homogeneous
`Γ(0.5, 0.4)`, M = 12, `logdet` criterion.

Practical notes on `Γ_pr`:
- It is **dense within each block**. Keep `n_elem` in the low thousands (`maxsz ≈ 1/20`
  on the unit disk gives ~1500 triangles) — memory is `n_elem²·8` bytes.
- The block-zero structure is a *feature*: build it as `blkdiag` of two dense blocks.
- `logdet(Γ_pr)` is θ-independent but keep it in the reported number so D-opt values stay
  comparable across configurations (same convention as the MDEIT study).
- Add a tiny jitter (`+1e-10 · trace/n · I`) if `chol` complains; a squared-exponential
  kernel at `λ=0.5` on a fine mesh is numerically rank-deficient.

---

## 6. Optimizer

Implement **Algorithm 1 as written** (`lib/steepest_descent_algorithm1.m`) — this is the
paper's contribution and the thing being compared:
1. `d = −∇ψ(θ)/|∇ψ(θ)|` (normalized steepest descent).
2. `t_max` = largest step keeping every gap `> g_min` (use `g_min = 0.05·w`, say); the
   paper says only "△ chosen so that the gaps remain positive".
3. `t_min = fminbnd(@(t) psi(θ + t·d), 0, t_max)`.
4. Iterate until `|∇ψ|` or the relative decrease drops below tolerance, or `max_iters`.

Also expose `optimizer = 'fminunc'` (quasi-Newton, `SpecifyObjectiveGradient`, exactly the
options block used in `example_anomaly_circle_2d.m`) as an alternative. **This matters for
the comparison the user wants**: the MDEIT study uses `fminunc`, so running both here
separates "the criterion behaves differently" from "the optimizer behaves differently".
Report both.

**Brute force (Case 1 only, M=4).** Mirror the structure of
`studies/optimal_sensors_bayesian_approach/brute_force_4x4_2d.m`: grid over `θ⁻ ∈ [0,2π)^4`
using `nchoosek(1:grid_n, 4)` (electrodes are unordered ⇒ combinations, not permutations;
reject combinations whose gaps go negative), then **locally polish the best grid point**
with the optimizer to remove grid discretization error. `grid_n = 36` (10° steps) for the
full run, `grid_n = 10` for `quick_test`. Save the whole `psi_grid` so the landscape can be
plotted later without rerunning. This reproduces the paper's own Case-1 validation
("verified by brute force simulations", Fig. 2).

---

## 7. Validation ladder — do these **in order**, do not skip

| # | check | target |
|---|---|---|
| V1 | Our CEM forward solve vs. EIDORS `fwd_solve` on a configuration whose electrode arcs coincide with mesh boundary nodes (so both discretizations agree) | rel. diff ≤ 1e-6 |
| V2 | Our `J` vs. (a) EIDORS `calc_jacobian` on the same aligned config, and (b) central differences in σ | 1e-6 / 1e-8 |
| V3 | `ψ(θ)` smoothness: plot ψ along a 1-D slice in one `θ⁻_m` at fine resolution | visibly smooth, no mesh-induced sawtooth |
| V4 | **`∇ψ` vs. central differences** (reuse `check_gradient_fd` from the MDEIT study; copy it into `lib/`), several components, both criteria, at a *perturbed* (non-symmetric) θ | rel. err. 1e-7 … 1e-8 |
| V5 | Case 1: Algorithm 1 output vs. brute-force-polished global minimum | ψ agrees to ≲0.1%; configurations coincide (paper Fig. 2) |
| V6 | Case 1 qualitative: two of the four electrodes migrate close to `Ω'` (the high-prior-variance disk) | matches Fig. 2 |
| V7 | Fig. 3 trend: sweep `κ_out ∈ {0.05, 0.1, 0.15, 0.2}` with A-optimality; optimal positions should spread out as the background uncertainty grows, and become non-symmetric w.r.t. Ω' at the two largest values | matches Fig. 3 |
| V8 | Case 2: with 12 electrodes, almost all electrodes move to the lower half (the high-variance half); A- and D-optimal results qualitatively similar | matches Fig. 4 |

V4 is the single most important check — every gradient bug in this repo's prior work was
caught exactly this way (see `Optimizing_Electrode_Positions_in_Electrical_Imped.md`, which
despite its filename is a handoff note, not the paper). **Do not run any full optimization
before V4 passes.**

Optional V9 (paper Fig. 5, the "was it worth it" experiment): draw `N_draw = 500`
conductivities from `N(σ_*, Γ_pr)`, simulate data on a **finer** mesh (avoid the inverse
crime) with noise `N(0, Γ_noise)`, compute a Gauss–Newton MAP estimate for both `θ_init`
and `θ_opt`, and report the ratio of mean squared errors. Paper gets ≈ 0.75. This is
expensive; make it a separate script `hyvonen_mc_evaluation.m` guarded by a flag.

---

## 8. Suggested phasing (each phase ends with something runnable)

1. **Mesh + CEM forward.** `lib/`: `make_disk_mesh.m`, `boundary_curve.m`,
   `electrode_arc_integrals.m`, `cem_system_matrix.m`, `cem_fwd_solve.m`. → V1, V3.
2. **Jacobian.** `lib/cem_jacobian.m` (adjoint, reusing the factorization). → V2.
3. **Criterion.** `lib/oed_criterion.m` (Woodbury `tr`/`logdet` + `W_A`/`W_D`),
   `lib/prior_covariance.m` (5.1 with the κ blocks), `lib/noise_covariance.m`. Cost-only
   evaluation of ψ at even spacing works at this point.
4. **Gradient.** `lib/cem_dJ_contract.m` (§4.2) + `lib/repulsion_term.m` +
   `lib/costgrad_electrodes.m` tying it together. → **V4** (gate).
5. **Optimization + Case 1.** `hyvonen_case1.m`, `hyvonen_case1_brute_force.m`. → V5, V6.
6. **Case 1 prior sweep + Case 2.** `hyvonen_case1_kappa_sweep.m`, `hyvonen_case2.m`.
   → V7, V8.
7. *(optional)* **Case 3 shapes** (`boundary_curve` gains 'ellipse'/'peanut') and
   **V9 Monte-Carlo MAP evaluation**.

Each script: `quick_test` flag at the top (coarse mesh, few iterations, gradient check on)
following the convention of the existing study; save to
`studies/hyvonen2014_electrode_optimization/data/<script_name>_results.mat`; figures to
`figures/`.

---

## 9. Reporting — make it directly comparable to the MDEIT study

For every case print and save the same quantities the magnetometer study prints, so the two
can be put side by side:

- `psi_init`, `psi_opt`, and **% reduction** for A-optimality (the MDEIT baseline reports
  "A-opt trace reduction %"; theirs was 4.6% for θ-only, ~2% for the retuned 4×4 config).
- D-optimality **information gain in nats and bits** (`(psi_init − psi_opt)/2` nats for
  `logdet`; the MDEIT study reported 17.57 nats / 25.35 bits).
- **Posterior-variance ratio in the region of interest** (here: `Ω'` for Case 1, the lower
  half-disk for Case 2), optimized/initial. The previous handoff argues this — not SNR — is
  the right "was optimizing worth it" metric; use the same one.
- Electrode angles (degrees) at init and optimum; iteration count; wall-clock time.
- Figures matching the paper: posterior point-variance map (diagonal of `Γ_*`) with
  electrodes drawn on the boundary, for init / A-opt / D-opt (paper Figs. 2, 4).

Then write a short `RESULTS.md` in this folder with a **comparison section** against
`studies/optimal_sensors_bayesian_approach/`: same criteria, same Woodbury evaluation, same
optimizer options — the substantive difference is the measurement physics (boundary
voltages under CEM vs. Biot–Savart magnetic-field sensors) and the design space (electrode
arcs constrained to ∂Ω vs. magnetometers free on a circle outside the domain). Flag in
particular whether the EIT case shows the same "A-optimality gains little, D-optimality
gains a lot" asymmetry the MDEIT study found, and whether the multiple-local-minima
behaviour seen in the 2D MDEIT case (`brute_force_4x4_2d.m`: fminunc converged to a
meaningfully worse minimum than brute force) also appears here.

---

## 10. Known traps

- **Grounding.** See §3.4. A gradient that is right up to an additive constant per
  electrode means the mean-free projection is inconsistent between `J` and `∂J/∂θ`.
- **2π wraparound** in the arc/edge intersection and in the gap computation. Write a unit
  test for `electrode_arc_integrals` with an electrode straddling θ = 0.
- **Electrode ordering.** Algorithm 1 can push electrodes past each other if `t_max` is not
  enforced; the repulsion term blows up first, but the line search must still be bounded.
- **`Γ_pr` conditioning.** Squared-exponential at λ=0.5 on a fine mesh is near-singular.
  The Woodbury form never inverts it, but `logdet(Γ_pr)` still needs `chol` → jitter.
- **`κ_ij = 0` cross-blocks** make `Γ_pr` block diagonal, *not* singular. Build it blockwise
  rather than forming the full dense outer product.
- **Sanity of the noise level.** `Γ_noise` is defined from the *voltage spread* at the
  initial configuration and then frozen (paper explicitly notes this simplification). Do not
  recompute it inside the objective — that would change what is being optimized.
- **Don't call Netgen in the loop.** If you ever find yourself remeshing per iteration,
  Deviation 2 has been lost.
