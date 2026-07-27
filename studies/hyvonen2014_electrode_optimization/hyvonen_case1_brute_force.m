%% Hyvonen et al. 2014 -- Case 1 brute-force validation (M=4)
%
% Mirrors studies/optimal_sensors_bayesian_approach/brute_force_4x4_2d.m:
% exhaustive grid search over the 4 (unordered) electrode angles, then a
% local polish of the best grid point, to check whether the gradient-based
% optimizers (Algorithm 1 and fminunc, in hyvonen_case1.m) find the true
% global optimum -- exactly the paper's own Case-1 validation strategy
% ("evaluate psi on a reasonably large set... brute force global minimizers
% are compared to the outputs of Algorithm 1", sec 5 / Fig. 2).
%
% Since the 4 electrodes are unordered on the circle, use
% nchoosek(1:grid_n,4) combinations (not grid_n^4 permutations), reject
% combinations whose arc-length gaps would go negative (electrodes
% touching/crossing), and evaluate psi (cost only, no gradient) on the
% survivors.
%
% Run "quick_test = true; hyvonen_case1_brute_force" for a fast smoke test
% (small grid_n).

fullpath = mfilename('fullpath');
script_folder = fileparts(fullpath);
cd(script_folder);
addpath(fullfile(script_folder,'lib'));

grandparent_folder = fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));
prepare_workspace(script_folder);

rng(1);

if ~exist('quick_test','var'), quick_test = false; end

%% Must match hyvonen_case1.m exactly (same mesh, prior, noise) or the
%% comparison is meaningless
shape = struct('type','disk','radius',1.0);
maxsz = 1/20;

M = 4;
width = pi/16;
z_contact = 1.0;
sigma_star = 1.0;
alpha_rep = 1e-4;

lambda = 0.5;
Omega_prime_center = [-0.55, 0.55];
Omega_prime_radius = 0.25;
kappa_in  = 0.4;
kappa_out = 0.03;

grid_n = 20;   % 18deg steps -> nchoosek(20,4) = 4845 combinations
if quick_test
    maxsz = 1/10;
    grid_n = 10;
end

opt_mode = 'a-opt';   % paper's Case 1 brute force is shown for A-optimality (Fig. 2)

ctx = build_mesh_ctx(shape, maxsz);
theta_init = 2*pi*(0:M-1)'/M;

centroids = element_centroids(ctx);
d_to_center = sqrt(sum((centroids - Omega_prime_center).^2,2));
category = 2*ones(ctx.n_elem,1);
category(d_to_center <= Omega_prime_radius) = 1;
kappa_lookup = [kappa_in, 0; 0, kappa_out];
Gamma_prior = prior_covariance(ctx, lambda, category, kappa_lookup);

cfg = struct('width',width,'z_contact',z_contact,'sigma_star',sigma_star, ...
    'Gamma_prior',Gamma_prior,'alpha_rep',alpha_rep);
cfg.Gamma_noise = calibrate_noise_covariance(ctx, theta_init, cfg);

cg = @(th) costgrad_electrodes(ctx, th, cfg, opt_mode);
gaps_fn = @(th) electrode_gaps(th, width, shape);

phi_even = cg(theta_init);
fprintf('phi_even = %.6e\n', phi_even);

%% fminunc from even spacing (for reference)
options = optimoptions('fminunc','Algorithm','quasi-newton', ...
    'SpecifyObjectiveGradient',true,'Display','iter', ...
    'MaxIterations',60,'OptimalityTolerance',1e-6,'StepTolerance',1e-9);
[theta_fmin, phi_fmin] = fminunc(cg, theta_init, options);
fprintf('phi_fminunc = %.6e\n', phi_fmin);

%% Exhaustive grid search
angles = 2*pi*(0:grid_n-1)'/grid_n;
combos = nchoosek(1:grid_n, M);   % [n_combos x 4], unordered
n_combos = size(combos,1);
fprintf('Evaluating %i combinations (grid_n=%i)...\n', n_combos, grid_n);

phi_grid = inf(n_combos,1);
min_gap_abs = 0.1*width/shape.radius;   % reject near-colliding combos outright

t0 = tic;
for c = 1:n_combos
    th = angles(combos(c,:));
    g = gaps_fn(th);
    if min(g) < min_gap_abs, continue; end
    try
        phi_grid(c) = cg(th);
    catch
        phi_grid(c) = inf;
    end
end
fprintf('Grid search done in %.1f s\n', toc(t0));

[phi_bf, idx_bf] = min(phi_grid);
theta_bf = angles(combos(idx_bf,:));
fprintf('phi_brute_force (grid) = %.6e at theta = %s deg\n', phi_bf, mat2str(rad2deg(theta_bf),4));

%% Polish the best grid point (removes grid discretization error)
[theta_bf_polished, phi_bf_polished] = fminunc(cg, theta_bf, options);
fprintf('phi_brute_force_polished = %.6e at theta = %s deg\n', ...
    phi_bf_polished, mat2str(rad2deg(mod(theta_bf_polished,2*pi)),4));

%% Comparison
fprintf('\n=== Case 1 brute-force validation (%s) ===\n', opt_mode);
fprintf('phi_even                 = %.6e\n', phi_even);
fprintf('phi_fminunc               = %.6e\n', phi_fmin);
fprintf('phi_brute_force (grid)    = %.6e\n', phi_bf);
fprintf('phi_brute_force_polished  = %.6e\n', phi_bf_polished);
gap_pct = 100*(phi_fmin - phi_bf_polished)/phi_even;
fprintf('gap (fminunc vs polished brute force) = %.4f%% of phi_even\n', gap_pct);

save(fullfile(script_folder,'data','hyvonen_case1_brute_force_results.mat'), ...
    'ctx','theta_init','cfg','category','phi_even','theta_fmin','phi_fmin', ...
    'combos','angles','phi_grid','theta_bf','phi_bf','theta_bf_polished','phi_bf_polished','grid_n');

fprintf('\nSaved to data/hyvonen_case1_brute_force_results.mat\n');
