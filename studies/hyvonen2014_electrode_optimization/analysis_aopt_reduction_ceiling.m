%% Why is the A-optimality cost reduction large for Hyvonen-EIT but tiny for MDEIT?
%
% This script computes a SHARP, MODEL-INDEPENDENT UPPER BOUND on the
% fractional A-optimality reduction that any sensor/electrode configuration
% can possibly achieve, and evaluates it for
%   (a) studies/optimal_sensors_bayesian_approach/example_anomaly_circle_2d.m
%       (MDEIT magnetometers, DIAGONAL prior)
%   (b) studies/hyvonen2014_electrode_optimization/hyvonen_case1.m
%       (EIT electrodes, DENSE Gaussian smoothness prior)
%
% ---------------------------------------------------------------------
% THE BOUND
%
% Write Gamma_pr = L*L'. Whitened Jacobian Jt = Gn^(-1/2) * J * L, with
% eigendecomposition Jt'*Jt = V*S^2*V' (V orthonormal, s_i >= 0). Then
%
%   Gamma_* = (J'Gn^-1 J + Gpr^-1)^-1 = L (I + Jt'Jt)^-1 L'
%   trace(Gamma_*)  = sum_i (1+s_i^2)^-1 * w_i ,   w_i := ||L v_i||^2
%   trace(Gamma_pr) = sum_i w_i
%
% so the fractional reduction is exactly
%
%   rho = 1 - trace(Gamma_*)/trace(Gamma_pr)
%       = [ sum_i (s_i^2/(1+s_i^2)) * w_i ] / [ sum_i w_i ]                (*)
%
% Two facts bound (*):
%   1. rank(Jt) <= n_data, so at most n_data of the s_i are nonzero.
%      More sharply, only modes with s_i^2 >> 1 (i.e. signal above the
%      noise floor) contribute meaningfully -- that count is what the
%      MDEIT script calls d_target / d_modes.
%   2. Each surviving factor s_i^2/(1+s_i^2) < 1.
%
% Maximizing sum of w_i = v_i'*Gamma_pr*v_i over orthonormal sets of size
% k gives, by Ky Fan's theorem, the sum of the k LARGEST EIGENVALUES of
% Gamma_pr. Hence
%
%   rho  <=  ( sum_{i=1}^{k} mu_i ) / ( sum_i mu_i ),     k = d_eff       (**)
%
% where mu_1 >= mu_2 >= ... are the eigenvalues of the PRIOR COVARIANCE.
%
% (**) involves NO forward model, NO sensor positions, NO physics -- only
% the prior and the number of noise-resolvable data modes. It is a hard
% information ceiling: no sensor arrangement whatsoever can beat it.
% ---------------------------------------------------------------------

fullpath = mfilename('fullpath');
script_folder = fileparts(fullpath);
cd(script_folder);
addpath(fullfile(script_folder,'lib'));
grandparent_folder = fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));
prepare_workspace(script_folder);

rng(1);

ceiling = @(mu,k) sum(mu(1:min(k,numel(mu))))/sum(mu);

%% =====================================================================
%% (a) MDEIT: rebuild example_anomaly_circle_2d.m's model + prior exactly
%% =====================================================================
z0 = 0.0058; l0 = 40e-3; I0 = 2.4e-3;
sigma0 = l0/z0;
mp = create_kai_2d_model_parameters(l0, z0, sigma0, I0);
mp.maxsz  = mp.radius/30;
mp.height = 0.6;
mp.numOfElectrodesPerRing = 4;
mp.numOfRings = 1;
mp.material.name = 'disk';
mp.material.type = 'cylindrical';
mp.material.radius = mp.radius/5;
mp.material.position = [+sqrt((mp.radius/2)^2/2), +sqrt((mp.radius/2)^2/2), 0];

n_sensors_mdeit = 4;
mp.numOfSensors = n_sensors_mdeit;
mp.sensorRadius = 1.05*mp.radius;
th_even_m = 2*pi*(0:n_sensors_mdeit-1)'/n_sensors_mdeit;   % theta_to_locations, inlined
mp.sensorPositions = [1.05*mp.radius*cos(th_even_m), 1.05*mp.radius*sin(th_even_m), ...
    (mp.height/2)*ones(n_sensors_mdeit,1)];

model_folder = fullfile(grandparent_folder,'models');
[~,fmdls] = mk_mdeit_model(mp, model_folder);
fmdl_m = fmdls{1};

n_elem_m = size(fmdl_m.elems,1);
n_stim_m = numel(fmdl_m.stimulation);
n_data_m = n_stim_m*n_sensors_mdeit;         % z-channel only

prior_std_background_m = 0.03;
prior_std_roi_m        = 0.4;
roi_radius_factor_m    = 1.0;

anomaly_center_m = mp.material.position;
roi_radius_m = roi_radius_factor_m*mp.material.radius;
cd_m = sqrt(sum((fmdl_m.elem_centroids - anomaly_center_m(1:2)).^2,2));
roi_m = cd_m <= roi_radius_m;

prior_var_m = prior_std_background_m^2*ones(n_elem_m,1);
prior_var_m(roi_m) = prior_std_roi_m^2;

% DIAGONAL prior -> its eigenvalues ARE the prior variances
mu_diag = sort(prior_var_m,'descend');

d_target_m = 10;   % the script's noise calibration makes ~d_target modes data-dominated

fprintf('\n===================== (a) MDEIT (diagonal prior) =====================\n');
fprintf('n_elem = %i, n_stim = %i, n_sensors = %i -> n_data = %i\n', ...
    n_elem_m, n_stim_m, n_sensors_mdeit, n_data_m);
fprintf('ROI elements: %i of %i (%.2f%% of mesh)\n', nnz(roi_m), n_elem_m, 100*nnz(roi_m)/n_elem_m);
fprintf('trace(Gamma_pr) = %.4f  (ROI share %.1f%%)\n', sum(prior_var_m), ...
    100*sum(prior_var_m(roi_m))/sum(prior_var_m));
fprintf('CEILING on A-opt reduction (k = d_target = %i): %.2f%%\n', d_target_m, 100*ceiling(mu_diag,d_target_m));
fprintf('CEILING on A-opt reduction (k = n_data   = %i): %.2f%%\n', n_data_m, 100*ceiling(mu_diag,n_data_m));
fprintf('  (reported achieved reduction in that study: ~2-5%%)\n');

%% What if the SAME MDEIT model used a Hyvonen-style smoothness prior?
fprintf('\n--- same MDEIT mesh/ROI, but a DENSE Gaussian smoothness prior (paper eq 5.1) ---\n');
centroids_m = fmdl_m.elem_centroids;
cat_m = 2*ones(n_elem_m,1); cat_m(roi_m) = 1;

lambda_list = [0.1 0.25 0.5 1.0]*mp.radius;
for lam = lambda_list
    G = build_smoothness_prior(centroids_m, cat_m, ...
        [prior_std_roi_m, 0; 0, prior_std_background_m], lam);
    mu = sort(eig((G+G')/2),'descend');
    mu = max(mu,0);
    fprintf('lambda = %.3f (= %.2f x tank radius): ceiling(k=%2i) = %6.2f%% | ceiling(k=%2i) = %6.2f%%\n', ...
        lam, lam/mp.radius, d_target_m, 100*ceiling(mu,d_target_m), n_data_m, 100*ceiling(mu,n_data_m));
end

%% =====================================================================
%% (b) Hyvonen Case 1: load the prior actually used
%% =====================================================================
S = load(fullfile(script_folder,'data','hyvonen_case1_results.mat'));
Gpr_h = S.cfg.Gamma_prior;
n_elem_h = size(Gpr_h,1);
M_h = 4;
n_data_h = M_h*(M_h-1);

mu_h = sort(eig((Gpr_h+Gpr_h')/2),'descend');
mu_h = max(mu_h,0);

fprintf('\n============ (b) Hyvonen Case 1 (dense smoothness prior) ============\n');
fprintf('n_elem = %i, M = %i -> n_data = %i\n', n_elem_h, M_h, n_data_h);
fprintf('trace(Gamma_pr) = %.4f\n', sum(mu_h));
fprintf('top-1 eigenvalue captures %.2f%% of the prior trace\n', 100*mu_h(1)/sum(mu_h));
fprintf('top-3                     %.2f%%\n', 100*sum(mu_h(1:3))/sum(mu_h));
fprintf('CEILING on A-opt reduction (k = n_data = %i): %.2f%%\n', n_data_h, 100*ceiling(mu_h,n_data_h));
fprintf('  (achieved in hyvonen_case1.m: 69.6%%)\n');

%% =====================================================================
%% Effective rank comparison
%% =====================================================================
eff_rank = @(mu) sum(mu)^2/sum(mu.^2);   % participation ratio
fprintf('\n===================== effective rank of the prior =====================\n');
fprintf('MDEIT diagonal prior      : effective rank = %8.1f  (n_elem = %i)\n', eff_rank(mu_diag), n_elem_m);
fprintf('Hyvonen smoothness prior  : effective rank = %8.1f  (n_elem = %i)\n', eff_rank(mu_h), n_elem_h);
fprintf('\nA design can only ever remove ~min(d_eff, effective rank) worth of prior trace.\n');

save(fullfile(script_folder,'data','aopt_ceiling_analysis.mat'), ...
    'mu_diag','mu_h','n_data_m','n_data_h','d_target_m','n_elem_m','n_elem_h');

%% ---------------------------------------------------------------------
function G = build_smoothness_prior(centroids, category, kappa_lookup, lambda)
% Gamma_ij = kappa_i*kappa_j * exp(-|x_i-x_j|^2/(2 lambda^2)), block form
% (zero covariance between categories), i.e. paper eq (5.1)/(5.2).
n = size(centroids,1);
G = zeros(n,n);
K = max(category);
for a = 1:K
    ia = find(category==a);
    for b = a:K
        kab = kappa_lookup(a,b);
        if kab == 0, continue; end
        ib = find(category==b);
        xa = centroids(ia,:); xb = centroids(ib,:);
        D2 = sum(xa.^2,2) + sum(xb.^2,2).' - 2*(xa*xb.');
        blk = kab^2*exp(-max(D2,0)/(2*lambda^2));
        G(ia,ib) = blk;
        if b ~= a, G(ib,ia) = blk.'; end
    end
end
end
