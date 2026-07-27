function Gamma_pr = prior_covariance(ctx, lambda, category, kappa_lookup)
%PRIOR_COVARIANCE  Gaussian smoothness prior, paper eq (5.1):
%
%   Gamma_ij = kappa_ij^2 * exp( -|x_i - x_j|^2 / (2*lambda^2) )
%
% evaluated at ELEMENT CENTROIDS x_i (Deviation 1 in PLAN_implementation.md:
% elementwise conductivity, not nodal, so this differs from the paper's
% nodal/background-mesh treatment but matches
% studies/optimal_sensors_bayesian_approach's convention).
%
%   Gamma_pr = prior_covariance(ctx, lambda, category, kappa_lookup)
%
% category    : [n_elem x 1] integer label per element (e.g. 1 = inside a
%               region of interest, 2 = background)
% kappa_lookup: [K x K] symmetric matrix, kappa_lookup(a,b) = kappa_ij for
%               any i in category a, j in category b. A zero entry makes
%               Gamma_pr exactly block-diagonal between those categories
%               (paper Cases 1/2 both use this to decouple ROI from
%               background) -- built blockwise below rather than as one
%               dense n_elem^2 exp(), so the zero blocks cost nothing.
%
% Gamma_pr is returned dense (n_elem in the low thousands for the mesh
% sizes used here; see PLAN sec 5 for the memory note).

n_elem = ctx.n_elem;
centroids = element_centroids(ctx);
K = max(category);

Gamma_pr = zeros(n_elem,n_elem);
for a = 1:K
    ia = find(category==a);
    for b = a:K
        kab = kappa_lookup(a,b);
        if kab == 0, continue; end
        ib = find(category==b);
        xa = centroids(ia,:); xb = centroids(ib,:);
        D2 = sum(xa.^2,2) + sum(xb.^2,2).' - 2*(xa*xb.');
        D2 = max(D2,0);
        block = kab^2 * exp(-D2/(2*lambda^2));
        Gamma_pr(ia,ib) = block;
        if b ~= a
            Gamma_pr(ib,ia) = block.';
        end
    end
end
end
