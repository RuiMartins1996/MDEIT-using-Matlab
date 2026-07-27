%% Hyvonen, Seppanen & Staboulis (2014) -- Case 2: semidisks with different
%% prior variances, M=12 electrodes (paper sec 5, Fig. 4)
%
% Unit disk, upper half (x2>=0) has LOW prior variance (kappa=0.03,
% "almost known"), lower half (x2<0) has HIGH prior variance (kappa=0.4),
% decoupled (kappa=0 across halves) -- paper eq (5.3). With M=12
% electrodes, exhaustive brute force is infeasible (this is the paper's
% point in choosing this case): the goal is simply to check that the
% qualitative output makes sense (electrodes should migrate toward the
% uncertain lower half) for both A- and D-optimality, in a
% higher-dimensional design space than Case 1.
%
% Run "quick_test = true; hyvonen_case2" for a fast smoke test.

fullpath = mfilename('fullpath');
script_folder = fileparts(fullpath);
cd(script_folder);
addpath(fullfile(script_folder,'lib'));

grandparent_folder = fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));
prepare_workspace(script_folder);

rng(1);

if ~exist('quick_test','var'), quick_test = false; end

shape = struct('type','disk','radius',1.0);
maxsz = 1/20;

M = 12;
width = pi/48;         % narrower electrodes to fit 12 comfortably
z_contact = 1.0;
sigma_star = 1.0;
alpha_rep = 1e-4;

lambda = 0.5;
kappa_upper = 0.03;    % paper eq (5.3): x2>=0 -> low variance
kappa_lower = 0.4;     % x2<0  -> high variance

max_iters = 60;
n_starts = 1;
if quick_test
    maxsz = 1/10;
    max_iters = 8;
end

opt_modes = {'a-opt','d-opt'};

ctx = build_mesh_ctx(shape, maxsz);
theta_init = 2*pi*(0:M-1)'/M;

centroids = element_centroids(ctx);
category = ones(ctx.n_elem,1);            % 1 = upper (x2>=0)
category(centroids(:,2) < 0) = 2;          % 2 = lower (x2<0)
kappa_lookup = [kappa_upper, 0; 0, kappa_lower];
Gamma_prior = prior_covariance(ctx, lambda, category, kappa_lookup);

fprintf('Case 2: %i elements, %i in lower half (%.1f%% of prior trace)\n', ...
    ctx.n_elem, nnz(category==2), 100*sum(diag(Gamma_prior).*(category==2))/sum(diag(Gamma_prior)));

cfg = struct('width',width,'z_contact',z_contact,'sigma_star',sigma_star, ...
    'Gamma_prior',Gamma_prior,'alpha_rep',alpha_rep);
cfg.Gamma_noise = calibrate_noise_covariance(ctx, theta_init, cfg);

if quick_test
    for om = opt_modes
        fprintf('\nGradient check (%s):\n', om{1});
        cg = @(th) costgrad_electrodes(ctx, th, cfg, om{1});
        check_gradient_fd(cg, theta_init + 0.02*randn(M,1), 3, 1e-3);
    end
end

[post_var_init, phi_a_init, phi_d_init] = posterior_variance_diag(ctx, theta_init, cfg);
fprintf('\nInitial (even spacing): phi_A = %.6e | phi_D = %.6e\n', phi_a_init, phi_d_init);

options = optimoptions('fminunc','Algorithm','quasi-newton', ...
    'SpecifyObjectiveGradient',true,'Display','iter', ...
    'MaxIterations',max_iters,'OptimalityTolerance',1e-6,'StepTolerance',1e-9);

results = struct();
roi = category==2;

for iom = 1:numel(opt_modes)
    om = opt_modes{iom};
    fprintf('\n=== Optimizing %s (Case 2, M=%i) ===\n', om, M);
    cg = @(th) costgrad_electrodes(ctx, th, cfg, om);

    tic
    [theta_opt, phi_opt] = fminunc(cg, theta_init, options);
    t_opt = toc;

    post_var_opt = posterior_variance_diag(ctx, theta_opt, cfg);

    results.(strrep(om,'-','_')) = struct('theta_opt',theta_opt,'phi_opt',phi_opt, ...
        't_opt',t_opt,'post_var_opt',post_var_opt);

    if strcmp(om,'a-opt')
        pct = 100*(phi_a_init-phi_opt)/phi_a_init;
        fprintf('%s: phi_init=%.6e -> phi_opt=%.6e (%.2f%% trace reduction, %.1fs)\n', ...
            om, phi_a_init, phi_opt, pct, t_opt);
    else
        nats = (phi_d_init-phi_opt)/2;
        fprintf('%s: phi_init=%.6e -> phi_opt=%.6e (%.4f nats / %.4f bits, %.1fs)\n', ...
            om, phi_d_init, phi_opt, nats, nats/log(2), t_opt);
    end
    fprintf('  lower-half posterior-variance ratio (opt/init): %.4f\n', ...
        sum(post_var_opt(roi))/sum(post_var_init(roi)));
    fprintf('  theta_opt (deg): %s\n', mat2str(rad2deg(mod(theta_opt,2*pi)),4));

    plot_case2_figure(ctx, shape, theta_opt, post_var_opt, sprintf('%s optimum', om), ...
        fullfile(script_folder,'figures',sprintf('case2_%s.png',strrep(om,'-','_'))));
end

plot_case2_figure(ctx, shape, theta_init, post_var_init, 'Initial (even spacing)', ...
    fullfile(script_folder,'figures','case2_init.png'));

save(fullfile(script_folder,'data','hyvonen_case2_results.mat'), ...
    'ctx','theta_init','cfg','category','post_var_init','phi_a_init','phi_d_init','results');

fprintf('\nSaved results to data/hyvonen_case2_results.mat and figures/case2_*.png\n');

function plot_case2_figure(ctx, shape, theta_minus, post_var, ttl, outfile)
figure('Visible','off');
patch('Faces',ctx.elems,'Vertices',ctx.nodes,'FaceVertexCData',post_var, ...
    'FaceColor','flat','EdgeColor','none');
axis equal off; colorbar; colormap(hot); title(ttl);
hold on;
M = numel(theta_minus);
width = pi/48;
[theta_plus,~] = solve_right_endpoint(theta_minus, width, shape);
for m = 1:M
    tt = linspace(theta_minus(m), theta_plus(m), 10);
    g = boundary_curve(tt, shape);
    plot(g(:,1)*1.02, g(:,2)*1.02, 'r-', 'LineWidth', 3);
end
saveas(gcf, outfile);
close(gcf);
end
