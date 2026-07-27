function [w, s, info] = build_consensus_weight_map(sE, sM, centroids, opts)
%BUILD_CONSENSUS_WEIGHT_MAP  Consensus soft-weight map for the EIT+MDEIT cascade.
%
%   [w, s, info] = build_consensus_weight_map(sE, sM, centroids, opts)
%
%   Implements the weight-map construction of Approach 3 (Section 5.3-5.4 and
%   Algorithm 1 of joint_eit_mdeit_strategies):
%
%     1. Normalize each single-modality magnitude map to [0,1].
%     2. Directional hygiene (safeguard 5.6.3): replace the EIT map by its
%        z-max-projection, so only its transverse (in-plane) localization is
%        used and its vertical elongation is discarded before it can enter
%        the prior.
%     3. Consensus map  s_k = sqrt(sE_proj_k * sM_k)  (geometric mean) or
%        s_k = min(sE_proj_k, sM_k).  Both demand agreement: a support that
%        only one modality lights up (e.g. EIT's vertical smear) is vetoed.
%     4. Soft weights  w_k = (epsilon + s_k)^(-p),  normalized to median 1
%        (keeps lambda comparable to the single-modality value).  Small w_k
%        where both modalities agree there is an anomaly; large w_k in the
%        background.  Flat s => flat w => Stage 2 reduces to plain Tikhonov,
%        so the cascade degrades gracefully when a stage is uninformative.
%
%   Inputs
%     sE         n x 1 nonnegative EIT stage-1 magnitude map (|dsigma_EIT|)
%     sM         n x 1 nonnegative MDEIT magnitude map (|dsigma_MDEIT|)
%     centroids  n x 3 element centroid coordinates (columns x, y, z)
%     opts       struct with fields
%                  .consensus       'geomean' (default) | 'min'
%                  .epsilon         floor in the soft weight (default 0.05)
%                  .p               exponent in the soft weight (default 0.5)
%                  .use_projection  z-max-project the EIT map (default true)
%                  .proj_cell       transverse grid cell size for the
%                                   projection (same units as centroids)
%                  .precomputed_sEp (optional) n x 1 fixed, already-projected
%                                   EIT factor.  When supplied, sE and the
%                                   projection are ignored and this factor is
%                                   reused -- used across IRLS iterations to
%                                   keep the EIT transverse factor fixed.
%
%   Outputs
%     w      n x 1 soft weights (median(w) == 1)
%     s      n x 1 normalized consensus map in [0,1]
%     info   struct: sEp (EIT factor actually used), sMn (normalized MDEIT)
%
%   See also RUN_CONSENSUS_CASCADE, SOLVE_WEIGHTED_TIKHONOV.

if nargin < 4, opts = struct(); end
opts = set_default(opts, 'consensus',      'geomean');
opts = set_default(opts, 'epsilon',        0.05);
opts = set_default(opts, 'p',              0.5);
opts = set_default(opts, 'use_projection', true);
opts = set_default(opts, 'proj_cell',      []);
opts = set_default(opts, 'precomputed_sEp', []);

sE = sE(:); sM = sM(:);
normalize = @(v) v ./ max([max(v), eps]);

sMn = normalize(sM);

if ~isempty(opts.precomputed_sEp)
    sEp = opts.precomputed_sEp(:);
else
    sEn = normalize(sE);
    if opts.use_projection
        assert(~isempty(opts.proj_cell) && opts.proj_cell > 0, ...
            'opts.proj_cell must be a positive scalar when use_projection is true');
        sEp = zmax_project(sEn, centroids, opts.proj_cell);
    else
        sEp = sEn;
    end
end

% --- Consensus (agreement) map ---
switch lower(opts.consensus)
    case 'geomean'
        s = sqrt(max(sEp, 0) .* max(sMn, 0));
    case 'min'
        s = min(sEp, sMn);
    otherwise
        error('Unknown consensus rule "%s" (use ''geomean'' or ''min'').', opts.consensus);
end
s = normalize(s);

% --- Soft weights, normalized to unit median ---
w = (opts.epsilon + s) .^ (-opts.p);
w = w ./ median(w);

info.sEp = sEp;
info.sMn = sMn;
info.s   = s;
end


% -------------------------------------------------------------------------
function sp = zmax_project(s, centroids, cell)
%ZMAX_PROJECT  Maximum projection along z onto a transverse grid.
%   Buckets elements by their (x,y) centroid into square cells of side
%   `cell`, assigns every element in a cell the maximum of s over that whole
%   column of z, and returns the result.  This extrudes the transverse
%   support along z: it keeps in-plane localization and deliberately erases
%   the vertical profile of the input map.

xy = centroids(:, 1:2);
ix = floor((xy(:,1) - min(xy(:,1))) / cell);
iy = floor((xy(:,2) - min(xy(:,2))) / cell);

key = ix + (max(ix) + 1) * iy;         % unique transverse-cell id
[~, ~, g] = unique(key);

% Column-wise max via accumarray, then scatter back to elements
cell_max = accumarray(g, s, [], @max);
sp = cell_max(g);
end


% -------------------------------------------------------------------------
function opts = set_default(opts, field, value)
if ~isfield(opts, field) || isempty(opts.(field))
    opts.(field) = value;
end
end
