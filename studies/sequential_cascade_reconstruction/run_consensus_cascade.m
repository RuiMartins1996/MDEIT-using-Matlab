function out = run_consensus_cascade(J_M, dy_M, J_E, dy_E, centroids, opts)
%RUN_CONSENSUS_CASCADE  Sequential (cascaded) EIT+MDEIT inversion, Approach 3.
%
%   out = run_consensus_cascade(J_M, dy_M, J_E, dy_E, centroids, opts)
%
%   Implements the consensus-weighted sequential cascade (Algorithm 1 of
%   joint_eit_mdeit_strategies.pdf).  Instead of stacking the two data
%   residuals into one least-squares problem -- which the rank analysis shows
%   cannot add information and which re-imports EIT's vertical elongation --
%   this uses one modality to shape the *regularizer* of the other.  The
%   information transferred is structural (where the anomaly is), so it is
%   immune to the weighting pathology of the stacked cost.
%
%   Pipeline:
%     1. Whiten both datasets by their (scalar) noise level so lambda is
%        comparable across modalities.
%     2. Stage 1: single-modality Tikhonov reconstructions (GCV lambda)
%        of EIT (dsE) and 3-axis MDEIT (dsM).
%     3. Build a consensus soft-weight map W: the EIT image contributes only
%        its z-max-projected transverse support, the MDEIT image caps the
%        vertical extent; a geometric-mean/min consensus vetoes any support
%        only one modality sees.
%     4. Stage 2 (IRLS, T iterations): re-solve MDEIT with penalty
%        lambda||W dsigma||^2, recomputing W from |dsigma^(t)| each step
%        while keeping the EIT transverse factor fixed.  With the soft
%        weights and p = 1/2 this is exactly IRLS for an l1-type prior, a
%        cheap approximation to edge-preserving joint reconstruction.
%
%   The data term in Stage 2 is MDEIT only, so EIT acts purely as a
%   veto/confirmation signal on support.  When EIT is uninformative the
%   consensus map flattens and the method degrades gracefully to MDEIT-only.
%
%   Inputs
%     J_M        m_M x n MDEIT Jacobian (3-axis stack [Jx;Jy;Jz])
%     dy_M       m_M x 1 MDEIT difference data (inhomogeneous - homogeneous)
%     J_E        m_E x n EIT Jacobian
%     dy_E       m_E x 1 EIT difference data
%     centroids  n x 3 element centroid coordinates
%     opts       struct with fields (all optional except where noted):
%                  .sigma_M       MDEIT noise std for whitening (default 1)
%                  .sigma_E       EIT   noise std for whitening (default 1)
%                  .lambda_vector candidate lambdas for GCV (default
%                                 logspace(-15,3,30); applied ridge is m*lambda)
%                  .T             IRLS iterations (default 3)
%                  .epsilon       soft-weight floor (default 0.05)
%                  .p             soft-weight exponent (default 0.5)
%                  .consensus     'geomean' (default) | 'min'
%                  .use_projection z-max-project EIT map (default true)
%                  .proj_cell     transverse projection cell (REQUIRED if
%                                 use_projection is true)
%                  .verbose       print progress (default true)
%
%   Output struct `out`
%     .dsigma        final cascaded reconstruction (n x 1)
%     .dsE, .dsM     Stage-1 EIT-only / MDEIT-only reconstructions
%     .lamE, .lamM   Stage-1 GCV ridge parameters (mu = m*lambda)
%     .w, .s         final weight map and consensus map
%     .sEp           fixed z-max-projected EIT factor
%     .history       1 x T struct array: .dsigma, .lambda, .w per iteration
%
%   See also SOLVE_WEIGHTED_TIKHONOV, BUILD_CONSENSUS_WEIGHT_MAP.

if nargin < 6, opts = struct(); end
opts = set_default(opts, 'sigma_M',        1);
opts = set_default(opts, 'sigma_E',        1);
opts = set_default(opts, 'lambda_vector',  logspace(-15, 3, 30));
opts = set_default(opts, 'T',              3);
opts = set_default(opts, 'epsilon',        0.05);
opts = set_default(opts, 'p',              0.5);
opts = set_default(opts, 'consensus',      'geomean');
opts = set_default(opts, 'use_projection', true);
opts = set_default(opts, 'proj_cell',      []);
opts = set_default(opts, 'verbose',        true);

n = size(J_M, 2);
assert(size(J_E, 2) == n, 'J_E and J_M must have the same number of columns');
assert(size(centroids, 1) == n, 'centroids must have one row per element');

vprint = @(varargin) opts.verbose && fprintf(varargin{:});

% --- 1. Whitening (scalar per modality) ---
JM = J_M ./ opts.sigma_M;   dM = dy_M(:) ./ opts.sigma_M;
JE = J_E ./ opts.sigma_E;   dE = dy_E(:) ./ opts.sigma_E;

% --- 2. Stage 1: single-modality reconstructions ---
vprint('[cascade] Stage 1a: EIT reconstruction\n');
[dsE, lamE] = solve_weighted_tikhonov(JE, dE, ones(n,1), opts.lambda_vector);
vprint('[cascade] Stage 1b: 3-axis MDEIT reconstruction\n');
[dsM, lamM] = solve_weighted_tikhonov(JM, dM, ones(n,1), opts.lambda_vector);

% --- 3. Consensus weight map (fix the EIT transverse factor once) ---
wopts = struct('consensus', opts.consensus, 'epsilon', opts.epsilon, ...
               'p', opts.p, 'use_projection', opts.use_projection, ...
               'proj_cell', opts.proj_cell);
[w, s, winfo] = build_consensus_weight_map(abs(dsE), abs(dsM), centroids, wopts);
sEp = winfo.sEp;                       % held fixed across the IRLS iterations

% --- 4. Stage 2: IRLS-reweighted MDEIT solves ---
history(opts.T) = struct('dsigma', [], 'lambda', [], 'w', []);
ds_cur = dsM;
for t = 1:opts.T
    vprint('[cascade] Stage 2, IRLS iteration %d/%d\n', t, opts.T);

    % Reweight from the current MDEIT estimate; EIT factor stays fixed.
    wopts_t = wopts;
    wopts_t.precomputed_sEp = sEp;
    [w, s] = build_consensus_weight_map([], abs(ds_cur), centroids, wopts_t);

    [ds_cur, lam_t] = solve_weighted_tikhonov(JM, dM, w, opts.lambda_vector);

    history(t).dsigma = ds_cur;
    history(t).lambda = lam_t;
    history(t).w      = w;
end

out.dsigma  = ds_cur;
out.dsE     = dsE;
out.dsM     = dsM;
out.lamE    = lamE;
out.lamM    = lamM;
out.w       = w;
out.s       = s;
out.sEp     = sEp;
out.history = history;
end


% -------------------------------------------------------------------------
function opts = set_default(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end
