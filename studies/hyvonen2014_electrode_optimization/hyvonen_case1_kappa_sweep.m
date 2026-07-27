%% Hyvonen et al. 2014 -- Case 1, kappa_out sweep (paper Fig. 3)
%
% Sweeps the background covariance factor kappa_out in eq (5.2) over
% {0.05, 0.10, 0.15, 0.20} with the A-optimality criterion, to reproduce
% the paper's Fig. 3: as the background uncertainty grows relative to
% Omega', the optimal electrode positions spread out and become
% non-symmetric with respect to Omega' (two mirror-image optima appear,
% and the run may land on either one depending on the optimizer path).
%
% Run "quick_test = true; hyvonen_case1_kappa_sweep" for a fast smoke test.

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

M = 4;
width = pi/16;
z_contact = 1.0;
sigma_star = 1.0;
alpha_rep = 1e-4;

lambda = 0.5;
Omega_prime_center = [-0.55, 0.55];
Omega_prime_radius = 0.25;
kappa_in = 0.4;
kappa_out_values = [0.05, 0.10, 0.15, 0.20];

max_iters = 60;
if quick_test
    maxsz = 1/10;
    max_iters = 10;
    kappa_out_values = [0.05, 0.20];
end

ctx = build_mesh_ctx(shape, maxsz);
theta_init = 2*pi*(0:M-1)'/M;

centroids = element_centroids(ctx);
d_to_center = sqrt(sum((centroids - Omega_prime_center).^2,2));
category = 2*ones(ctx.n_elem,1);
category(d_to_center <= Omega_prime_radius) = 1;

gaps_fn = @(th) electrode_gaps(th, width, shape);
options = optimoptions('fminunc','Algorithm','quasi-newton', ...
    'SpecifyObjectiveGradient',true,'Display','off', ...
    'MaxIterations',max_iters,'OptimalityTolerance',1e-6,'StepTolerance',1e-9);

n_k = numel(kappa_out_values);
theta_opts = cell(1,n_k);
phi_opts = zeros(1,n_k);

for ik = 1:n_k
    kappa_out = kappa_out_values(ik);
    kappa_lookup = [kappa_in, 0; 0, kappa_out];
    Gamma_prior = prior_covariance(ctx, lambda, category, kappa_lookup);

    cfg = struct('width',width,'z_contact',z_contact,'sigma_star',sigma_star, ...
        'Gamma_prior',Gamma_prior,'alpha_rep',alpha_rep);
    cfg.Gamma_noise = calibrate_noise_covariance(ctx, theta_init, cfg);

    cg = @(th) costgrad_electrodes(ctx, th, cfg, 'a-opt');
    [theta_opt, phi_opt] = fminunc(cg, theta_init, options);

    theta_opts{ik} = theta_opt;
    phi_opts(ik) = phi_opt;

    fprintf('kappa_out=%.2f: phi_opt=%.6e, theta=%s deg\n', ...
        kappa_out, phi_opt, mat2str(rad2deg(mod(theta_opt,2*pi)),4));

    plot_kappa_figure(ctx, shape, theta_opt, category, kappa_out, ...
        fullfile(script_folder,'figures',sprintf('case1_kappa_sweep_%02d.png',ik)));
end

save(fullfile(script_folder,'data','hyvonen_case1_kappa_sweep_results.mat'), ...
    'kappa_out_values','theta_opts','phi_opts','category','ctx');

fprintf('\nSaved to data/hyvonen_case1_kappa_sweep_results.mat\n');

function plot_kappa_figure(ctx, shape, theta_minus, category, kappa_out, outfile)
figure('Visible','off');
patch('Faces',ctx.elems,'Vertices',ctx.nodes,'FaceVertexCData',double(category), ...
    'FaceColor','flat','EdgeColor','none');
axis equal off; colormap(hot);
title(sprintf('kappa_{out} = %.2f', kappa_out));
hold on;
M = numel(theta_minus);
width = pi/16;
[theta_plus,~] = solve_right_endpoint(theta_minus, width, shape);
for m = 1:M
    tt = linspace(theta_minus(m), theta_plus(m), 10);
    g = boundary_curve(tt, shape);
    plot(g(:,1)*1.02, g(:,2)*1.02, 'r-', 'LineWidth', 3);
end
saveas(gcf, outfile);
close(gcf);
end
