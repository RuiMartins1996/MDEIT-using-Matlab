function [dsigma, lambda_opt, info] = solve_weighted_tikhonov(J, dy, w, lambda_vector)
%SOLVE_WEIGHTED_TIKHONOV  Diagonally-weighted Tikhonov solve, GCV lambda.
%
%   [dsigma, lambda_opt, info] = solve_weighted_tikhonov(J, dy, w, lambda_vector)
%
%   Solves the linearized, spatially-weighted, zeroth-order Tikhonov problem
%
%       min_x  || J x - dy ||^2  +  mu * || W x ||^2 ,        W = diag(w),
%
%   which is the Stage-2 problem of the sequential (cascaded) EIT+MDEIT
%   reconstruction (Approach 3, eq. (8) of joint_eit_mdeit_strategies).
%   Small w_k lets the data act freely at element k ("anomaly likely here");
%   large w_k suppresses element k towards the background ("probably
%   background").
%
%   The weighting is handled by the standard substitution
%
%       z = W x ,   Jt = J * W^{-1} ,   min_z ||Jt z - dy||^2 + mu||z||^2 ,
%
%   which is ordinary Tikhonov on z, solved through the (economy) SVD of Jt.
%   The regularization parameter is chosen by Generalized Cross-Validation
%   (GCV), using the same cost and convention as functions/
%   generalized_cross_validation.m: the value scanned in lambda_vector is
%   related to the ridge actually applied to the normal equations by
%   mu = m*lambda, where m is the number of measurements.  Passing
%   w = ones(n,1) recovers the plain single-modality reconstruction used in
%   Stage 1.
%
%   Inputs
%     J             m x n Jacobian (noise-whitened upstream if desired)
%     dy            m x 1 difference data (inhomogeneous - homogeneous)
%     w             n x 1 strictly-positive weights
%     lambda_vector vector of candidate lambdas scanned by GCV (the applied
%                   ridge is mu = m*lambda)
%
%   Outputs
%     dsigma        n x 1 reconstructed conductivity update x
%     lambda_opt    selected ridge mu = m*lambda actually used in the solve
%     info          struct: V_lambda (GCV cost over lambda_vector), svals
%                   (singular values of Jt), lambda_vector, optimal_id
%
%   See also RUN_CONSENSUS_CASCADE, BUILD_CONSENSUS_WEIGHT_MAP,
%   GENERALIZED_CROSS_VALIDATION.

w  = w(:);
dy = dy(:);

assert(all(isfinite(w)) && all(w > 0), 'weights w must be finite and > 0');
assert(size(J,1) == numel(dy), 'J and dy dimensions are inconsistent');
assert(size(J,2) == numel(w),  'J and w dimensions are inconsistent');

% --- Substitution J -> Jt = J * diag(1/w) (scales the columns of J) ---
Jt = J .* (1./w).';

% --- Economy SVD of the weighted operator ---
[U, S, V] = svd(Jt, 'econ');
s  = diag(S);
Uy = U.' * dy;                 % data expressed in the left-singular basis
m  = size(Jt, 1);              % number of measurements

% --- Generalized Cross-Validation over the lambda grid ---
% Mirrors generalized_cross_validation.m: the shrinkage factors use m*lambda,
% and the reported optimum mu = m*lambda is the ridge applied to (Jt'Jt + mu I).
nl       = numel(lambda_vector);
V_lambda = zeros(nl, 1);
for i = 1:nl
    lam   = lambda_vector(i);
    gamma = s.^2 ./ (s.^2 + m*lam);            % data-space shrinkage factors
    omg   = 1 - gamma;
    den   = (1/m) * sum(omg);
    if den < eps                               % degenerate (mu ~ 0): reject
        V_lambda(i) = Inf;
    else
        V_lambda(i) = ((1/m) * sum((omg .* Uy).^2)) / den^2;
    end
end

[~, optimal_id] = min(V_lambda);               % min() ignores NaN by default
mu_opt = m * lambda_vector(optimal_id);        % effective ridge parameter

% --- Reconstruct at mu_opt and undo the weighting substitution ---
z      = V * ((s ./ (s.^2 + mu_opt)) .* Uy);
dsigma = z ./ w;

lambda_opt = mu_opt;

info.V_lambda      = V_lambda;
info.svals         = s;
info.lambda_vector = lambda_vector(:);
info.optimal_id    = optimal_id;
end
