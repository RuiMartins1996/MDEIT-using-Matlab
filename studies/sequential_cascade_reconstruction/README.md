# Sequential (cascaded) EIT + MDEIT reconstruction

Implementation of **Approach 3** from `joint_eit_mdeit_strategies.pdf`:
*sequential inversion with cross-modality priors*. Instead of stacking the EIT
and MDEIT data residuals into one augmented least-squares problem — which the
rank analysis shows adds no numerical row space and which re-imports EIT's
vertical elongation — one modality is used to shape the **regularizer** of the
other. The information transferred is *structural* (where the anomaly is), so it
is immune to the weighting pathology of the stacked cost.

## Pipeline (Algorithm 1)

1. **Whiten** both datasets by their noise level so `lambda` is comparable
   across modalities.
2. **Stage 1** — single-modality Tikhonov reconstructions (GCV `lambda`):
   `dsE` (EIT) and `dsM` (3-axis MDEIT).
3. **Consensus weight map** `W`:
   - the EIT image contributes only its **z-max-projected transverse support**
     (directional hygiene — its vertical smear is discarded before it can enter
     the prior);
   - a **geometric-mean / min consensus** with MDEIT demands agreement, so any
     support only one modality lights up is vetoed;
   - soft weights `w_k = (epsilon + s_k)^(-p)`, normalized to median 1. Small
     `w_k` = "anomaly likely, let the data act freely"; large `w_k` = "background,
     suppress." Flat `s` ⇒ flat `w` ⇒ Stage 2 reduces to plain Tikhonov (graceful
     degradation to MDEIT-only).
4. **Stage 2 (IRLS)** — re-solve MDEIT with penalty `lambda * ||W dsigma||^2`,
   recomputing `W` from `|dsigma^(t)|` each step while holding the EIT transverse
   factor fixed. With `p = 1/2` this is IRLS for an ℓ1-type / edge-preserving
   prior — the class of priors that can actually defeat minimum-norm vertical
   smearing. The Stage-2 data term is MDEIT only; **EIT acts purely as a
   veto/confirmation signal on support.**

## Files

| File | Role |
|------|------|
| `sequential_cascade_reconstruction.m` | Main driver: builds the 3-D tank model, generates data, sweeps SNR, runs the cascade, plots, reports metrics, saves results. |
| `run_consensus_cascade.m` | The whole cascade (Algorithm 1). Reusable; takes Jacobians + data + centroids. |
| `build_consensus_weight_map.m` | Consensus soft-weight map, including the z-max projection. |
| `solve_weighted_tikhonov.m` | Weighted zeroth-order Tikhonov solve `min ||Jx-dy||^2 + mu||Wx||^2` via the `z = Wx` substitution + SVD, with GCV selection of the ridge (`mu = m*lambda`, matching `functions/generalized_cross_validation.m`). |
| `cascade_metrics.m` | Localization error and a moment-based `rmsZ/rmsXY` elongation proxy (robust on unstructured tetrahedral meshes). |

## Running

```matlab
% from MATLAB, in this folder:
sequential_cascade_reconstruction
```

Key knobs at the top of the driver: `SNR_list` (default `[100 20 1]`),
`anomaly_case` (`'side'` or `'center'`), `contrast`, the two mesh sizes
`maxsz_data` / `maxsz_reconstruction`, and the cascade hyperparameters
(`cascade_T`, `cascade_epsilon`, `cascade_p`, `cascade_rule`). Coarser meshes run
faster; refine for accuracy. Results are saved to `data/cascade_<case>_anomaly.mat`.

## What to look for

- **Elongation** (`rmsZ/rmsXY`) reduced by the cascade relative to MDEIT-only,
  especially for the **side** anomaly (where EIT's transverse localization is
  good).
- **Null test** returns near-flat images (no support hallucinated from noise).
- **Held-out EIT misfit** stays near its noise level — the weight map is not
  overconfident.
- Graceful degradation to MDEIT-only at **SNR = 1**, where EIT is uninformative.

## Notes / caveats

- Whitening here is a **scalar** per modality (matches how the studies generate
  noise). The structure is ready for a full per-channel `Sigma^{-1/2}`: whiten
  `J` and `dy` before calling `run_consensus_cascade` and set `sigma_* = 1`.
- The z-max projection buckets elements by `(x,y)` centroid into cells of size
  `proj_cell`; it extrudes the transverse support along `z`.
- The elongation metric is a second-moment proxy for `FWHMz/FWHMxy`, chosen to
  avoid gridding the unstructured mesh.
