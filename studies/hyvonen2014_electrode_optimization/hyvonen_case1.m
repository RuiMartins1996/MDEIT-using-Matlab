%% Hyvonen, Seppanen & Staboulis (2014) -- Case 1: circular subdomain with
%% high prior variance, M=4 electrodes (paper sec 5, Fig. 2/3)
%
% Reproduces the paper's motivating validation case: a unit disk with a
% smaller "region of interest" disk Omega' carrying much higher prior
% variance than the background. Four electrodes are optimized for both
% A- and D-optimality; with only 4 electrodes the result can (separately,
% in hyvonen_case1_brute_force.m) be checked against an exhaustive grid
% search, exactly as the paper does.
%
% See PLAN_implementation.md for the full derivation and the two
% deliberate deviations from the paper (elementwise conductivity, fixed
% mesh + continuous arc electrodes instead of remeshing).
%
% Run "quick_test = true; hyvonen_case1" for a fast smoke test.

fullpath = mfilename('fullpath');
script_folder = fileparts(fullpath);
cd(script_folder);
addpath(fullfile(script_folder,'lib'));

grandparent_folder = fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));
prepare_workspace(script_folder);

rng(1);

if ~exist('quick_test','var'), quick_test = false; end

%% Model / prior / optimization parameters (paper sec 4-5)
shape = struct('type','disk','radius',1.0);
maxsz = 1/20;

M = 4;
width = pi/16;         % paper's motivating example electrode width
z_contact = 1.0;        % paper sec 5: unit contact impedances
sigma_star = 1.0;       % paper sec 5: homogeneous prior mean

alpha_rep = 1e-4;       % paper sec 4: alpha in eq (4.3)

lambda = 0.5;           % correlation length, Gamma(0.5,kappa) (paper sec 5)
Omega_prime_center = [-0.55, 0.55];
Omega_prime_radius = 0.25;
kappa_in  = 0.4;        % paper eq (5.2)
kappa_out = 0.03;

max_iters_alg1 = 60;
max_iters_fminunc = 60;

if quick_test
    maxsz = 1/10;
    max_iters_alg1 = 8;
    max_iters_fminunc = 8;
end

opt_modes = {'a-opt','d-opt'};

%% Mesh (fixed for the whole run -- Deviation 2, no remeshing per iterate)
ctx = build_mesh_ctx(shape, maxsz);

%% Electrode initial configuration: evenly spaced
theta_init = 2*pi*(0:M-1)'/M;

%% Prior covariance (paper eq (5.1)-(5.2)): block-diagonal in (Omega', background)
centroids = element_centroids(ctx);
d_to_center = sqrt(sum((centroids - Omega_prime_center).^2,2));
category = 2*ones(ctx.n_elem,1);
category(d_to_center <= Omega_prime_radius) = 1;   % 1 = inside Omega'

kappa_lookup = [kappa_in, 0; 0, kappa_out];
Gamma_prior = prior_covariance(ctx, lambda, category, kappa_lookup);

fprintf('Case 1: %i elements, %i in Omega'' (%.1f%% of prior trace)\n', ...
    ctx.n_elem, nnz(category==1), 100*sum(diag(Gamma_prior).*(category==1))/sum(diag(Gamma_prior)));

%% Noise covariance, calibrated once at the initial configuration
cfg = struct('width',width,'z_contact',z_contact,'sigma_star',sigma_star, ...
    'Gamma_prior',Gamma_prior,'alpha_rep',alpha_rep);
cfg.Gamma_noise = calibrate_noise_covariance(ctx, theta_init, cfg);

%% Optional smoke-test gradient check (do this before any optimization run)
if quick_test
    for om = opt_modes
        fprintf('\nGradient check (%s):\n', om{1});
        cg = @(th) costgrad_electrodes(ctx, th, cfg, om{1});
        check_gradient_fd(cg, theta_init + 0.05*randn(M,1), 3, 1e-3);
    end
end

%% Baseline (even spacing)
[post_var_init, phi_a_init, phi_d_init] = posterior_variance_diag(ctx, theta_init, cfg);
fprintf('\nInitial (even spacing): phi_A = %.6e | phi_D = %.6e\n', phi_a_init, phi_d_init);

%% Optimize with BOTH the paper's Algorithm 1 and fminunc (for comparison)
results = struct();
shape_for_gaps = shape;
gaps_fn = @(th) electrode_gaps(th, width, shape_for_gaps);

for iom = 1:numel(opt_modes)
    om = opt_modes{iom};
    fprintf('\n=== Optimizing %s (Case 1, M=%i) ===\n', om, M);
    cg = @(th) costgrad_electrodes(ctx, th, cfg, om);

    tic
    [theta_alg1, phi_alg1, info_alg1] = steepest_descent_algorithm1(cg, gaps_fn, theta_init, ...
        struct('max_iters',max_iters_alg1,'verbose',true));
    t_alg1 = toc;

    options = optimoptions('fminunc','Algorithm','quasi-newton', ...
        'SpecifyObjectiveGradient',true,'Display','iter', ...
        'MaxIterations',max_iters_fminunc,'OptimalityTolerance',1e-6,'StepTolerance',1e-9);
    tic
    [theta_fmin, phi_fmin] = fminunc(cg, theta_init, options);
    t_fmin = toc;

    [post_var_alg1] = posterior_variance_diag(ctx, theta_alg1, cfg);
    [post_var_fmin]  = posterior_variance_diag(ctx, theta_fmin, cfg);

    results.(strrep(om,'-','_')) = struct( ...
        'theta_alg1',theta_alg1,'phi_alg1',phi_alg1,'t_alg1',t_alg1,'info_alg1',info_alg1, ...
        'theta_fmin',theta_fmin,'phi_fmin',phi_fmin,'t_fmin',t_fmin, ...
        'post_var_alg1',post_var_alg1,'post_var_fmin',post_var_fmin);

    fprintf('%s: phi_init=%.6e | Algorithm1: phi=%.6e (%.1fs) | fminunc: phi=%.6e (%.1fs)\n', ...
        om, iif(strcmp(om,'a-opt'),phi_a_init,phi_d_init), phi_alg1, t_alg1, phi_fmin, t_fmin);
end

%% Show figures




%% Reporting (comparable format to studies/optimal_sensors_bayesian_approach)
roi = category==1;
fprintf('\n=================== Case 1 summary ===================\n');
fprintf('theta_init (deg): %s\n', mat2str(rad2deg(mod(theta_init,2*pi)),4));

for iom = 1:numel(opt_modes)
    om = opt_modes{iom};
    R = results.(strrep(om,'-','_'));
    fprintf('\n--- %s ---\n', om);
    if strcmp(om,'a-opt')
        phi_init = phi_a_init;
        pct_alg1 = 100*(phi_init-R.phi_alg1)/phi_init;
        pct_fmin = 100*(phi_init-R.phi_fmin)/phi_init;
        fprintf('  phi_init = %.6e\n', phi_init);
        fprintf('  Algorithm1: phi = %.6e  (%.2f%% trace reduction)\n', R.phi_alg1, pct_alg1);
        fprintf('  fminunc:    phi = %.6e  (%.2f%% trace reduction)\n', R.phi_fmin, pct_fmin);
    else
        phi_init = phi_d_init;
        nats_alg1 = (phi_init - R.phi_alg1)/2;
        nats_fmin = (phi_init - R.phi_fmin)/2;
        fprintf('  phi_init = %.6e\n', phi_init);
        fprintf('  Algorithm1: phi = %.6e  (%.4f nats / %.4f bits info gain)\n', R.phi_alg1, nats_alg1, nats_alg1/log(2));
        fprintf('  fminunc:    phi = %.6e  (%.4f nats / %.4f bits info gain)\n', R.phi_fmin, nats_fmin, nats_fmin/log(2));
    end
    fprintf('  ROI posterior-variance ratio (Algorithm1/init): %.4f\n', ...
        sum(R.post_var_alg1(roi))/sum(post_var_init(roi)));
    fprintf('  ROI posterior-variance ratio (fminunc/init):    %.4f\n', ...
        sum(R.post_var_fmin(roi))/sum(post_var_init(roi)));
    fprintf('  theta_alg1 (deg): %s\n', mat2str(rad2deg(mod(R.theta_alg1,2*pi)),4));
    fprintf('  theta_fmin (deg): %s\n', mat2str(rad2deg(mod(R.theta_fmin,2*pi)),4));
end

%% Figures: posterior point-variance map with electrodes (paper Fig. 2 style)
plot_case1_figure(ctx, shape, theta_init, post_var_init, 'Initial (even spacing)', ...
    fullfile(script_folder,'figures','case1_init.png'));
for iom = 1:numel(opt_modes)
    om = opt_modes{iom};
    R = results.(strrep(om,'-','_'));
    plot_case1_figure(ctx, shape, R.theta_alg1, R.post_var_alg1, ...
        sprintf('%s optimum (Algorithm 1)', om), ...
        fullfile(script_folder,'figures',sprintf('case1_%s_alg1.png',strrep(om,'-','_'))));
end

save(fullfile(script_folder,'data','hyvonen_case1_results.mat'), ...
    'ctx','theta_init','cfg','category','post_var_init','phi_a_init','phi_d_init','results');

fprintf('\nSaved results to data/hyvonen_case1_results.mat and figures/case1_*.png\n');

%% ------------------------------------------------------------------
function out = iif(cond,a,b)
if cond, out = a; else, out = b; end
end

function plot_case1_figure(ctx, shape, theta_minus, post_var, ttl, outfile)
figure('Visible','off');
patch('Faces',ctx.elems,'Vertices',ctx.nodes,'FaceVertexCData',post_var, ...
    'FaceColor','flat','EdgeColor','none');
axis equal off; colorbar; colormap(hot); title(ttl);
hold on;
M = numel(theta_minus);
width = pi/16; % only used for drawing extent, cosmetic
[theta_plus,~] = solve_right_endpoint(theta_minus, width, shape);
for m = 1:M
    tt = linspace(theta_minus(m), theta_plus(m), 10);
    [g,~,~] = boundary_curve(tt, shape);
    plot(g(:,1)*1.02, g(:,2)*1.02, 'r-', 'LineWidth', 3);
end
saveas(gcf, outfile);
close(gcf);
end
