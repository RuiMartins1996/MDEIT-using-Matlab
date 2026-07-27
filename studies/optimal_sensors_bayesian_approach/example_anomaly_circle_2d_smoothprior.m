%% EXAMPLE (2D, SMOOTHNESS PRIOR): Bayesian sensor-position optimization (MDEIT)
%
% ####################################################################
% VARIANT of example_anomaly_circle_2d.m. The ONLY substantive change is
% the PRIOR COVARIANCE: a dense Gaussian smoothness prior (Hyvonen,
% Seppanen & Staboulis 2014, eq (5.1)-(5.2)) instead of the diagonal
% (white) prior used by the original script.
%
% WHY: with a diagonal prior, the A-optimality cost reduction achievable
% by ANY sensor arrangement is bounded above by
%      rho <= (sum of the d_eff largest eigenvalues of Gamma_prior)
%             / trace(Gamma_prior)
% (Ky Fan; see studies/hyvonen2014_electrode_optimization/
% analysis_aopt_reduction_ceiling.m for the derivation). A diagonal prior
% on n_elem elements has n_elem equal-sized eigenvalues inside the ROI, so
% with only d_eff ~ 10 noise-resolvable data modes that ceiling is a few
% percent -- which is why the original script's A-opt reduction is small.
% It is an information ceiling set by the PRIOR, not a deficiency of the
% optimizer or of MDEIT physics. A smoothness prior with correlation
% length comparable to the ROI has effective rank of order 1-10, lifting
% the ceiling to ~90-99%.
%
% Three code changes were needed relative to the original:
%   1. the prior block (below, "SMOOTHNESS PRIOR");
%   2. the noise calibration whitens with chol(Gamma_prior) instead of
%      sqrt(diag(Gamma_prior)) -- otherwise d_target is miscalibrated;
%   3. costgrad_theta_z's gradient weight uses Yd = Y*Gamma_prior instead
%      of the diagonal shortcut Yd = Y.*d.' (the original is only valid
%      for a diagonal prior; using it with a dense prior silently yields a
%      WRONG gradient). Verified by finite differences -- keep
%      do_gradient_check = true when changing anything here.
% ####################################################################
%
% Same study as example_anomaly_circle.m, but the domain is a flat 2D
% circular region (a disk) instead of a 3D cylindrical tank -- i.e. the
% conductivity mesh is a triangular mesh (fmdl.elems is Nx3, fmdl.nodes is
% Nx2), not a tetrahedral extrusion. Sensors are still 3D points (they sit
% at height zs above the plane; the Biot-Savart integral is over the flat
% current sheet z=0, which is exactly how EIDORS's own 2D MDEIT pipeline
% treats it -- see functions/fwd_model/compute_r_matrices.m,
% compute_r_matrices_quadrature_2d).
%
% A conductive circular (disk) anomaly sits off-center. The prior encodes
% "we know roughly WHERE the anomaly may be, but not its value": large
% prior variance in a disk around the suspected anomaly location, small
% variance in the background (cf. Hyvonen, Seppanen & Staboulis, Case 1).
% Because the informative region is off-center, the optimal sensors are
% NOT evenly spaced.
%
% Design criteria (Bz-only MDEIT, linearized around the prior mean):
%   A-optimality: minimize trace(H^-1),      H = J'*inv(Gn)*J + inv(Gpr)
%   D-optimality: minimize -log det(H)  (=  log det posterior covariance)
% Only the z-channel (Bz) magnetometer measurements enter the design and
% the reconstruction; J is the z-only Jacobian dBz/dsigma.
%
% The gradients are analytic: dphi = -2*<W, dJ/dp>_F with
%   W_A = inv(Gn)*J*H^-2,   W_D = inv(Gn)*J*H^-1,
% and dJ/dp assembled from (i) extra adjoint solves for dlambda/dp and
% (ii) the explicit derivative of the Biot-Savart kernel integrals -- the
% derivation is identical to the 3D case; only the element quadrature
% (triangles + 37-point TOMS rule instead of tets + 35-point rule) and the
% Biot-Savart source geometry (flat sheet at z=0 instead of a 3D tet mesh)
% change. G.Gz is identically zero for a 2D mesh (EIDORS convention,
% functions/compute_gradient_matrix.m), so the shared calc_Jz/grad_z code
% needs no branching: the z=0 terms just vanish on their own.
%
% Run "quick_test = true; example_anomaly_circle_2d" for a fast smoke test
% (few iterations + finite-difference gradient check).

%% Prepare workspace
fullpath = mfilename('fullpath');
script_folder = fileparts(fullpath);
cd(script_folder);

grandparent_folder = fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));

model_folder = prepare_workspace(script_folder);

rng(1);

if ~exist('quick_test','var'), quick_test = false; end

%% Model parameters
z0 = 0.0058;  %(Ohm m^2) contact impedance
l0 = 40e-3;   %(m) tank radius (characteristic length)
I0 = 2.4e-3;  %(A) injected current magnitude

sigma0 = l0/z0; %(S/m) characteristic conductivity

model_parameters = create_kai_2d_model_parameters(l0, z0, sigma0, I0);

model_parameters.maxsz = model_parameters.radius/30;
model_parameters.height = 0.6;    % NOT a mesh dimension (mesh is flat, z=0);
                                   % only used below to set the sensors'
                                   % height above the plane (zs = height/2),
                                   % kept for numerical continuity with the
                                   % 3D example.

% Sweep hooks: set any of these in the base workspace before running to
% override the defaults (same idiom as quick_test). Used by
% sweep_smoothprior_configs.m to explore how large the A-optimality
% reduction can be made.
if ~exist('sw_n_electrodes','var'), sw_n_electrodes = 4; end
if ~exist('sw_n_sensors','var'),    sw_n_sensors    = 4; end
if ~exist('sw_d_target','var'),     sw_d_target     = 4; end
if ~exist('sw_anom_offset','var'),  sw_anom_offset  = 0.5; end   % x radius
if ~exist('sw_lambda_fac','var'),   sw_lambda_fac   = 1.0; end   % x ROI radius
if ~exist('sw_bg_std','var'),       sw_bg_std       = 0.03; end  % prior_std_background

model_parameters.numOfElectrodesPerRing = sw_n_electrodes;
model_parameters.numOfRings = 1;

% Conductive circular (disk) anomaly, off-center
model_parameters.material.name   = 'disk';
model_parameters.material.type   = 'cylindrical';   % 2D analog of a sphere: a disk
model_parameters.material.radius = model_parameters.radius/5;
model_parameters.material.position = ...
    [-sw_anom_offset*model_parameters.radius/sqrt(2), ...
     -sw_anom_offset*model_parameters.radius/sqrt(2), 0];

% Sensors: M evenly spaced on a circle of radius rs, at height zs above
% the (flat, z=0) domain
n_sensors = sw_n_sensors;
rs = 1.05*model_parameters.radius;        % circle radius (fixed, close to domain)
zs = model_parameters.height/2;           % sensor height above the plane (fixed)

theta_even = 2*pi*(0:n_sensors-1)'/n_sensors;

model_parameters.numOfSensors = n_sensors;
model_parameters.sensorRadius = rs;
model_parameters.sensorPositions = theta_to_locations(theta_even,rs,zs);

background_conductivity = 1.0;
anomaly_conductivity = 5*background_conductivity;  % conductive anomaly

%% Simulation parameters
opt_modes = {'a-opt','d-opt'};
opt_modes = {'a-opt'};

max_iterations = 30;
n_starts = 1;                  % 1 = start from the even configuration only
do_gradient_check = true;      % ON by default here: this variant changed the
                               % gradient (dense-prior Yd), so it must be
                               % FD-verified on every run.

% See example_anomaly_circle.m for the full discussion of d_target vs. ROI
% size and why this noise calibration matters for design sensitivity.
d_target = sw_d_target;      % desired # of data-dominated modes

prior_std_background = sw_bg_std;
prior_std_roi        = 0.4;
roi_radius_factor    = 1.0;    % ROI disk radius = factor * anomaly radius

% Correlation length of the smoothness prior. This is THE parameter that
% controls the A-optimality ceiling: the prior's effective rank falls as
% prior_lambda grows, and the ceiling rises accordingly. Set to the ROI
% radius here (a physically honest statement of "the anomaly is a blob of
% roughly this size"), which already collapses the effective rank to ~1-3.
prior_lambda = sw_lambda_fac*model_parameters.material.radius;

alpha_rep = 1e-5;                 % pairwise sensor repulsion weight (optional)

if quick_test
    max_iterations = 3;
    n_starts = 1;
    do_gradient_check = true;
end

%% Make forward models (2D: fmdl.nodes is Nx2, fmdl.elems is Nx3 triangles)
[~,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_fwd = fmdls{1};

% model_parameters_recon = model_parameters;
% model_parameters_recon.maxsz = model_parameters_recon.radius/20;
% model_parameters_recon = rmfield(model_parameters_recon,'anomaly');
% model_parameters_recon = rmfield(model_parameters_recon,'material');
% [~,fmdls] = mk_mdeit_model(model_parameters_recon,model_folder);
% fmdl_recon = fmdls{1};

fmdl_recon = fmdl_fwd;

% Homogeneous image = prior mean (the OED linearization point)
imgh_fwd = mk_image_mdeit(fmdl_fwd,background_conductivity);
imgh_recon = mk_image_mdeit(fmdl_recon,background_conductivity);

% Image with the anomaly (ground truth, used for visualization only)
imgi_fwd = add_material_properties(imgh_fwd,[background_conductivity,anomaly_conductivity]);

n_stim  = numel(fmdl_recon.stimulation);
n_elem  = size(fmdl_recon.elems,1);
n_nodes = size(fmdl_recon.nodes,1);

fprintf('Model (2D): %i elements, %i nodes, %i stims, %i sensors (Bz only)\n',...
    n_elem,n_nodes,n_stim,n_sensors);

A = @(x) M(imgh_recon,x);

%% Prior covariance: informative about WHERE the anomaly may be
anomaly_center = model_parameters.material.position;   % [x,y,~] (z unused in 2D)
roi_radius = roi_radius_factor*model_parameters.material.radius;

centroid_dist = sqrt(sum((fmdl_recon.elem_centroids - anomaly_center(1:2)).^2,2));
roi = centroid_dist <= roi_radius;

prior_variance = prior_std_background^2*ones(n_elem,1);
prior_variance(roi) = prior_std_roi^2;

fprintf('Building smooth prior\n');

str = sprintf("(x-%2.2f).^2+(y-%2.2f).^2<%2.2f^2",...
    model_parameters.anomaly.position(1),model_parameters.anomaly.position(2),model_parameters.anomaly.radius);

select_fcn = inline(str,'x','y','z');
idx = elem_select(imgh_recon.fwd_model, select_fcn);
idx(idx>0) = 1;
in_block  = idx * idx.';  % both in
out_block = (~idx) * (~idx).';% both out

% Assemble kappa
kappa = prior_std_roi * sparse(in_block) + prior_std_background * sparse(out_block);


function Gamma = smooth_prior(img, lambda, kappa)

    % Number of elements
    n_elem = size(img.fwd_model.elems,1);
    n_nodes = size(img.fwd_model.nodes,1);
    % centroids = img.
    % --- Compute squared Euclidean distance matrix efficiently ---
    % ||x_i - x_j||^2 = ||x_i||^2 + ||x_j||^2 - 2 x_i·x_j

    % Compute centroids
    nodes = img.fwd_model.nodes; 
    elems = img.fwd_model.elems;

    % -------------------------------------------------------------
    % 1) Compute element centroids
    % -------------------------------------------------------------
    % Works for triangles (2D) and tetrahedra (3D)

    dim = size(nodes,2);
    n_vert = size(elems,2);

    centroids = zeros(n_elem, dim);
    for n = 1:n_vert
        centroids = centroids + nodes(elems(:,n), :);
    end
    centroids = centroids / n_vert;

    % -------------------------------------------------------------
    % 2) Build smoothness matrix
    % -------------------------------------------------------------

    if issparse(kappa)
        % Sparse-efficient version (recommended)

        [i,j,val] = find(kappa);

        xi = centroids(i,:);
        xj = centroids(j,:);

        D2 = sum((xi - xj).^2, 2);

        weights = val.^2 .* exp(-D2/(2*lambda^2));

        Gamma = sparse(i,j,weights,n_elem,n_elem);

    else
        % Dense version (only for small problems)

        x2 = sum(centroids.^2, 2);
        D2 = bsxfun(@plus, x2, x2') - 2*(centroids*centroids');
        D2 = max(D2,0);  % numerical safety

        G = exp(-D2/(2*lambda^2));
        Gamma = kappa .* G;
    end
end

% Gamma_smooth = smooth_prior(imgh_recon,0.99,kappa);
% fprintf('\t Inverting smooth prior (might take a while)\n');
% Gamma_prior = Gamma_smooth;
% % Gamma_prior = Gamma_prior + 1e-6 * speye(size(Gamma_prior)); %need regularization because of ill-conditioning
% inv_Gamma_prior = Gamma_prior \ speye(size(Gamma_prior));

% ---------------------------------------------------------------------
% SMOOTHNESS PRIOR (this is the ONLY substantive change vs.
% example_anomaly_circle_2d.m -- see the header comment).
%   Gamma_ij = kappa_i*kappa_j * exp(-|x_i-x_j|^2/(2*lambda^2))
% with kappa = prior_std_roi inside the ROI, prior_std_background outside,
% and ZERO covariance across the two regions -- i.e. exactly Hyvonen,
% Seppanen & Staboulis (2014) eq (5.1)-(5.2). Built blockwise (dense) so
% the zero cross-block costs nothing.
% ---------------------------------------------------------------------
centroids_pr = fmdl_recon.elem_centroids;
kappa_vec = prior_std_background*ones(n_elem,1);
kappa_vec(roi) = prior_std_roi;

Gamma_prior = zeros(n_elem,n_elem);
for blk = {roi, ~roi}
    ib = find(blk{1});
    xb = centroids_pr(ib,:);
    D2 = sum(xb.^2,2) + sum(xb.^2,2).' - 2*(xb*xb.');
    Gamma_prior(ib,ib) = (kappa_vec(ib)*kappa_vec(ib).') .* exp(-max(D2,0)/(2*prior_lambda^2));
end
Gamma_prior = (Gamma_prior+Gamma_prior.')/2;
inv_Gamma_prior = [];  %#ok<NASGU>  (never needed: Woodbury path only)

prior_variance = diag(Gamma_prior);   % keep downstream figures consistent

mu_pr = sort(eig(Gamma_prior),'descend'); mu_pr = max(mu_pr,0);
fprintf('ROI: %i of %i elements (%.1f%% of prior trace)\n',nnz(roi),n_elem,...
    100*sum(diag(Gamma_prior).*roi)/sum(diag(Gamma_prior)));
fprintf('Smoothness prior: lambda = %.4f (= %.2f x radius), effective rank = %.1f\n',...
    prior_lambda, prior_lambda/model_parameters.radius, sum(mu_pr)^2/sum(mu_pr.^2));
fprintf('A-opt reduction CEILING (top-%i eigenvalues / trace) = %.2f%%\n',...
    d_target, 100*sum(mu_pr(1:min(d_target,n_elem)))/sum(mu_pr));

figure
img_temp = imgh_recon;
img_temp.elem_data = diag(Gamma_prior);
show_fem(img_temp);

drawnow

%% Noise covariance (white), calibrated from the whitened Jacobian spectrum
imgh_fwd = assign_sensor_locations(imgh_fwd,theta_to_locations(theta_even,rs,zs));
Bh = fwd_solve_mdeit(imgh_fwd);

% Giving the same result
% R1 = compute_R_2d(ctx,theta_to_locations(theta_even,rs,zs));
% [Rx,Ry,Rz,fmdl] = compute_r_matrices_local(fmdl_recon,theta_to_locations(theta_even,rs,zs));

imgi_fwd = assign_sensor_locations(imgi_fwd,theta_to_locations(theta_even,rs,zs));
Bi = fwd_solve_mdeit(imgi_fwd);

% Design uses ONLY the z-channel (Bz) measurements and Jacobian.
dB_even = Bi.Bz(:) - Bh.Bz(:);

max_B = max(abs(Bh.Bz(:)));

% Precompute everything that does NOT depend on the sensor positions
ctx = build_ctx(imgh_recon);

J0 = calc_Jz(ctx,theta_to_locations(theta_even,rs,zs));
% Whiten with the FULL prior square root (the diagonal shortcut
% "J0 .* sqrt(prior_variance)'" of the original script is only correct for
% a diagonal prior; with the dense smoothness prior it would badly
% miscalibrate d_target and hence the noise level).
L_pr = chol(Gamma_prior + 1e-12*trace(Gamma_prior)/n_elem*eye(n_elem),'lower');
[~,S,Vs] = svd(J0*L_pr,'econ');                       % prior-whitened spectrum
s = diag(S);

noise_std = s(d_target);
variance_noise = noise_std^2;

n_data = n_stim*n_sensors;                 % z-channel only
Gamma_noise = (variance_noise)*speye(n_data,n_data);
inv_Gamma_noise = (1/variance_noise)*speye(n_data,n_data); %#ok<NASGU>

d_modes = sum(s.^2/variance_noise > 1);
roi_energy = mean(sum(Vs(roi,1:d_modes).^2,1));  % ROI alignment of the modes
fprintf('noise std = %.3e (= max|B|/%.1f)\n',noise_std,max_B/noise_std);
fprintf('whitened spectrum s_i^2/noise_var: max = %.2e, #>1 = %i of %i data / %i params\n',...
    s(1)^2/variance_noise,d_modes,n_data,n_elem);
fprintf('mean ROI energy of the %i data-dominated modes: %.2f (want close to 1)\n',...
    d_modes,roi_energy);

%% Cost/gradient handles (theta is the only free variable)

costgrad = cell(1,numel(opt_modes));
for iom = 1:numel(opt_modes)
    costgrad{iom} = @(theta) costgrad_theta_z(ctx,theta,rs,zs,...
        Gamma_prior,Gamma_noise,opt_modes{iom},alpha_rep);
end

%% Optional: finite-difference gradient check
if do_gradient_check
    for iom = 1:numel(opt_modes)
        fprintf('\nGradient check (%s):\n',opt_modes{iom});
        check_gradient_fd(costgrad{iom},theta_even + 0.1*randn(n_sensors,1),3,1e-4);
    end
end

%% Baseline: evenly spaced sensors
phi_even = zeros(1,numel(opt_modes));
for iom = 1:numel(opt_modes)
    phi_even(iom) = costgrad{iom}(theta_even);
end

%% Optimize with fminunc (quasi-Newton, analytic gradients)
options = optimoptions('fminunc',...
    'Algorithm','quasi-newton','SpecifyObjectiveGradient',true,...
    'Display','iter','MaxIterations',max_iterations,...
    'OptimalityTolerance',1e-7,'StepTolerance',1e-9);

theta_opt = cell(1,numel(opt_modes));
phi_opt = inf(1,numel(opt_modes));

for iom = 1:numel(opt_modes)
    fprintf('\n=== Optimizing %s (Bz-only MDEIT, 2D domain) ===\n',opt_modes{iom});

    for istart = 1:n_starts
        if istart == 1
            th0 = theta_even;
        else
            th0 = 2*pi*rand(n_sensors,1);
        end

        tstart = tic;
        [th,fval] = fminunc(costgrad{iom},th0,options);
        fprintf('start %i: phi = %.6e (%.1f s)\n',istart,fval,toc(tstart));

        if fval < phi_opt(iom)
            phi_opt(iom) = fval;
            theta_opt{iom} = th;
        end
    end
end

%% Report cost reduction
fprintf('\n================ RESULTS ================\n');
for iom = 1:numel(opt_modes)
    fprintf('\n%s:\n',opt_modes{iom});
    fprintf('  phi(even spacing) = %.6e\n',phi_even(iom));
    fprintf('  phi(optimized)    = %.6e\n',phi_opt(iom));
    switch opt_modes{iom}
        case 'a-opt'
            fprintf('  trace reduction   = %.1f%%\n',...
                100*(1 - phi_opt(iom)/phi_even(iom)));
        case 'd-opt'
            dnats = phi_even(iom) - phi_opt(iom);
            fprintf('  extra information = %.2f nats (%.2f bits)\n',...
                dnats,dnats/log(2));
    end
end

%% Reconstruct with MDEIT for both sensor configurations

disp('Doing reconstruction')

img_prior = imgh_recon;
img_even = assign_sensor_locations(imgh_recon,theta_to_locations(theta_even,rs,zs));
img_aopt = assign_sensor_locations(imgh_recon,theta_to_locations(theta_opt{1},rs,zs));

% Jacobian (z-channel only, shared by cost and gradient)
J_even = calc_Jz(ctx,theta_to_locations(theta_even,rs,zs));
J_opt = calc_Jz(ctx,theta_to_locations(theta_opt{1},rs,zs));

% Compute posterior variances for both optimized and even sensor
% configurations

Gamma_post_even = Gamma_prior - (Gamma_prior*J_even.')*((Gamma_noise+J_even*Gamma_prior*J_even.')\(J_even*Gamma_prior));
Gamma_post_opt = Gamma_prior - (Gamma_prior*J_opt.')*((Gamma_noise+J_opt*Gamma_prior*J_opt.')\(J_opt*Gamma_prior));


% Compute data for optimal config
imgh_fwd = assign_sensor_locations(imgh_fwd,theta_to_locations(theta_opt{1},rs,zs));
Bh_opt = fwd_solve_mdeit(imgh_fwd);
imgi_fwd = assign_sensor_locations(imgi_fwd,theta_to_locations(theta_opt{1},rs,zs));
Bi_opt = fwd_solve_mdeit(imgi_fwd);

dB_opt = Bi_opt.Bz(:) - Bh_opt.Bz(:);

% Add measurement noise (every channel is a Bz channel now)
noise_even = sqrt(variance_noise)*randn(size(dB_even));
dB_even_noisy = dB_even + noise_even;

noise_opt = sqrt(variance_noise)*randn(size(dB_opt));
dB_opt_noisy = dB_opt + noise_opt;

x0 = background_conductivity*ones(n_elem,1);
x0(roi) = anomaly_conductivity;

dx0 = x0-background_conductivity*ones(n_elem,1);

posterior_mean_even = x0+Gamma_prior*J_even.'*((J_even*Gamma_prior*J_even.'+Gamma_noise)\(dB_even_noisy-J_even*dx0));
posterior_mean_opt = x0+Gamma_prior*J_opt.'*((J_opt*Gamma_prior*J_opt.'+Gamma_noise)\(dB_opt_noisy-J_opt*dx0));

%%
figure
subplot(1,3,1)
hold on
plot(dB_even_noisy)
plot(dB_even,'r')
title('Even config')
legend('noisy','true')
hold off

subplot(1,3,2)
hold on
plot(dB_opt_noisy)
plot(dB_opt,'r')
title('Opt config')
legend('noisy','true')
hold off

subplot(1,3,3)
hold on
plot(dB_opt)
plot(dB_even)
legend('opt','even')

% Display SNR for even and opt config

snr_even = 10*log10(rms(dB_even)^2/noise_std^2);
snr_opt = 10*log10(rms(dB_opt)^2/noise_std^2);

fprintf('SNR (even config) : %2.2g (dB)\n',snr_even)
fprintf('SNR (opt config) : %2.2g (dB)\n',snr_opt)

img_prior.elem_data = diag(Gamma_prior);
img_even.elem_data = diag(Gamma_post_even);
img_aopt.elem_data = diag(Gamma_post_opt);

figure('Name','Sensor positions');
subplot(1,3,1);
show_fem(imgi_fwd);axis off

hold on
% ROI disk outline (2D analog of the 3D ROI sphere)
% t = linspace(0,2*pi,100);
% fill(roi_radius*cos(t)+anomaly_center(1), roi_radius*sin(t)+anomaly_center(2), ...
%     [0 0.45 0.74],'FaceAlpha',0.3,'EdgeColor','none');
% title('Blue disk is ROI')
axis equal; box on; grid on;

subplot(1,3,2)
show_fem(img_even)
plot_sensors(img_even,false,'k','o');
axis equal; box on; grid on;axis off

subplot(1,3,3)
show_fem(img_aopt)
plot_sensors(img_aopt,false,'g','s');
axis equal; box on; grid on;axis off
%%
close all;
img_prior.calc_colours.npoints = 512;
img_even.calc_colours.npoints = 512;
img_aopt.calc_colours.npoints = 512;

figure('Name','Reconstruction');
subplot(1,4,1)
show_fem(img_prior)
colorbar
title('Prior variance')
subplot(1,4,2)
show_fem(img_even)
colorbar
title('Posterior variance (even)')
subplot(1,4,3)
show_fem(img_aopt)
colorbar
title('Posterior variance (a-opt)')

subplot(1,4,4)
hold on
a = diag(Gamma_prior);
b = diag(Gamma_post_even);
c = diag(Gamma_post_opt);

plot(a(roi),'r-');
plot(b(roi),'g-');
plot(c(roi),'b-');
axis square
title('Prior and Post variances in ROI')
legend('Prior','Even','Opt','Location','southeast');

%% Posterior variance diagnostics (prior vs even vs A-optimized)
post_var_even = compute_posterior_variance_diag(ctx,theta_even,rs,zs,...
    Gamma_prior,Gamma_noise);
post_var_aopt = compute_posterior_variance_diag(ctx,theta_opt{1},rs,zs,...
    Gamma_prior,Gamma_noise);

fprintf('\nPosterior variance sums (ROI | background):\n');
fprintf('  prior      : %.3e | %.3e\n',sum(prior_variance(roi)),sum(prior_variance(~roi)));
fprintf('  even       : %.3e | %.3e\n',sum(post_var_even(roi)),sum(post_var_even(~roi)));
fprintf('  a-optimized: %.3e | %.3e\n',sum(post_var_aopt(roi)),sum(post_var_aopt(~roi)));
fprintf('  ROI posterior variance: optimized/even = %.3f\n',...
    sum(post_var_aopt(roi))/sum(post_var_even(roi)));

%% Save results
results = struct();
results.theta_even = theta_even;
results.theta_opt = {theta_opt};
results.opt_modes = {opt_modes};
results.phi_even = phi_even;
results.phi_opt = phi_opt;
results.post_var_even = post_var_even;
results.post_var_aopt = post_var_aopt;
results.prior_variance = prior_variance;
results.roi = roi;
results.model_parameters = model_parameters;
results.noise_std = noise_std;
% NB: distinct filename from example_anomaly_circle_2d.m's output, so this
% variant never clobbers the baseline study's saved results.
save(fullfile(script_folder,'data','example_anomaly_circle_2d_smoothprior_results.mat'),'results');

% Sweep hook: when driven by sweep_smoothprior_configs.m, drop a one-row
% summary keyed by sw_tag so the driver can aggregate across separate
% MATLAB processes (this script shares variable names with its caller --
% e.g. it assigns 'c' at the reconstruction step -- so running the sweep
% in a single workspace is not safe).
if exist('sw_tag','var')
    sweep_row = struct('tag',sw_tag, ...
        'ceiling',100*sum(mu_pr(1:min(d_target,numel(mu_pr))))/sum(mu_pr), ...
        'phi_even',phi_even(1),'phi_opt',phi_opt(1), ...
        'reduction',100*(1-phi_opt(1)/phi_even(1)), ...
        'n_elec',model_parameters.numOfElectrodesPerRing,'n_sens',n_sensors, ...
        'd_target',d_target,'lambda',prior_lambda,'eff_rank',sum(mu_pr)^2/sum(mu_pr.^2), ...
        'roi_ratio',sum(post_var_aopt(roi))/sum(post_var_even(roi)));
    save(fullfile(script_folder,'data',['sweep_' sw_tag '.mat']),'sweep_row');
end

%% Figures
figures_folder = fullfile(script_folder,'figures');
if ~exist(figures_folder,'dir'), mkdir(figures_folder); end

% --- Sensor positions on the FEM model ---
img_even = assign_sensor_locations(imgh_recon,theta_to_locations(theta_even,rs,zs));
img_aopt = assign_sensor_locations(imgh_recon,theta_to_locations(theta_opt{1},rs,zs));

figure('Name','Sensor positions');
hold on
show_fem(imgi_fwd);
plot_sensors(img_even,false,'b','s');
plot_sensors(img_aopt,false,'r','o');
hold off
axis([-1.1*rs 1.1*rs -1.1*rs 1.1*rs]); axis equal;
box off; grid on;axis off
% title('black o: even | green s: a-opt');
savefig(fullfile(figures_folder,'example_anomaly_circle_2d_sensors.fig'));
saveas(gcf,fullfile(figures_folder,'example_anomaly_circle_2d_sensors.png'));

%% 
try
    figure('Name','Sensor angles');
    polarplot(theta_even,ones(n_sensors,1),'ko','DisplayName','even'); hold on
    try polarplot(mod(theta_opt{1},2*pi),1*ones(n_sensors,1),'gs','DisplayName','a-opt'); end
    polarplot([0 0],[0 1.2],'b-','DisplayName','anomaly azimuth');
    legend show
    title('Sensor azimuths (anomaly at azimuth 0)');
    savefig(fullfile(figures_folder,'example_anomaly_circle_2d_angles.fig'));
    saveas(gcf,fullfile(figures_folder,'example_anomaly_circle_2d_angles.png'));
catch err
    fprintf('WARNING: angle figure failed: %s\n',err.message);
end

% --- Posterior std over the (flat, single-slice) domain ---
try
    xc = fmdl_recon.elem_centroids(:,1);
    yc = fmdl_recon.elem_centroids(:,2);
    cmax = sqrt(max(prior_variance));

    figure('Name','Posterior std','Position',[100 100 1200 400]);
    vals = {sqrt(prior_variance),sqrt(post_var_even),sqrt(post_var_aopt)};
    names = {'prior','posterior (even)','posterior (a-opt)'};
    for ip = 1:3
        subplot(1,3,ip);
        scatter(xc,yc,25,vals{ip},'filled');
        axis equal tight; colorbar; clim([0 cmax]);
        title(names{ip}); xlabel('x'); ylabel('y');
    end
    savefig(fullfile(figures_folder,'example_anomaly_circle_2d_posterior.fig'));
    saveas(gcf,fullfile(figures_folder,'example_anomaly_circle_2d_posterior.png'));
catch err
    fprintf('WARNING: posterior figure failed: %s\n',err.message);
end

fprintf('\nDone.\n');

%% ======================== LOCAL FUNCTIONS ========================

%% Map circle angles to 3D sensor locations (sensors are 3D points; the
%% mesh they sense is a flat z=0 disk)
function sensor_locations = theta_to_locations(theta,rs,zs)
theta = theta(:);
sensor_locations = [rs*cos(theta), rs*sin(theta), zs*ones(numel(theta),1)];
end

%% Cost and gradient in the angle parametrization (Bz-only MDEIT)
% Uses the precomputed context CTX (see build_ctx): sigma is frozen at the
% prior mean, so the EIT forward solve, the factored CEM operator and the
% mesh geometry are all loop invariants and are never recomputed here.
% Identical to the 3D version -- calc_Jz/grad_z are dimension-agnostic
% (G.Gz is zero for a 2D mesh, so the z-terms just vanish).
function [phi,dphi] = costgrad_theta_z(ctx,theta,rs,zs,Gamma_prior,Gamma_noise,opt_mode,alpha_rep)
assert(ismember(opt_mode,{'a-opt','d-opt'}));

sensor_locations = theta_to_locations(theta,rs,zs);

% z-channel Jacobian (shared by cost and gradient)
J = calc_Jz(ctx,sensor_locations);

% Data-space (Woodbury) form: never forms the n x n posterior covariance.
%   H^-1 = Gpr - Gpr J' S^-1 J Gpr,   S = Gn + J Gpr J'   (nd x nd, SPD)
P  = J * Gamma_prior;                  % [nd x n]
S  = Gamma_noise + P * J.';            % [nd x nd]
Ls = chol(S,'lower');                  % factor once (reused for Y and logdet)
Y  = Ls.' \ (Ls \ P);                  % S^-1 P  [nd x n]
d  = diag(Gamma_prior);                % prior variances [n x 1]

switch opt_mode
    case 'a-opt'
        phi_data = sum(d) - sum(sum(P .* Y));          % trace(H^-1)
    case 'd-opt'
        logdet_S = 2*sum(log(diag(Ls)));
        phi_data = sum(log(diag(Gamma_noise))) ...
                 - sum(log(d)) - logdet_S;
end

% Optional pairwise sensor repulsion, scaled by the objective value so the
% penalty is a fixed *fraction* of the data cost regardless of its
% magnitude (alpha_rep is now scale-free).
if alpha_rep > 0
    Grep = repulsion_cost(sensor_locations);
    phi  = phi_data*(1 + alpha_rep*Grep);
else
    phi  = phi_data;
end

if nargout < 2
    return
end

% ---- Gradient ----
% dphi = -2*<W, dJ/dp>_F  with  W_A = Gn^-1*J*H^-2,  W_D = Gn^-1*J*H^-1.
switch opt_mode
    case 'a-opt'
        % Y*Gpr. The original script uses the diagonal shortcut
        % "Yd = Y .* d.'", which is WRONG for the dense smoothness prior
        % used by this variant -- must be a full matrix product.
        Yd = Y * Gamma_prior;        % Y*Gpr   [nd x n]
        W  = Yd - (Yd*J.')*Y;        % Gn^-1*J*H^-2
    case 'd-opt'
        W  = Y;                      % Gn^-1*J*H^-1
end

dphidp = grad_z(ctx,sensor_locations,W);             % [3 x n_sensors], p=x,y,z

% Chain rule to the angle parametrization (evaluated at the CURRENT theta)
dphi_data = -rs*sin(theta(:)).*dphidp(1,:)' + rs*cos(theta(:)).*dphidp(2,:)';

if alpha_rep > 0
    [dGx,dGy,~] = repulsion_grad_cartesian(sensor_locations);
    dGrep = -rs*sin(theta(:)).*dGx + rs*cos(theta(:)).*dGy;
    dphi = dphi_data*(1 + alpha_rep*Grep) + phi_data*alpha_rep*dGrep;
else
    dphi = dphi_data;
end

end

%% Precompute everything independent of the sensor positions (sigma fixed)
function ctx = build_ctx(img)
fmdl = img.fwd_model;

ctx.G         = fmdl.G;               % Gz is identically zero (2D mesh)
ctx.mu_factor = fmdl.mu0/(4*pi);
ctx.sigma     = img.elem_data(:);
ctx.n_elem    = size(fmdl.elems,1);
ctx.n_nodes   = size(fmdl.nodes,1);
ctx.num_stim  = numel(fmdl.stimulation);
ctx.elemV     = fmdl.elem_volume(:);  % element AREAS in 2D (field name kept
                                       % for parity with the 3D script)

% EIT forward solve (depends only on sigma) and its element gradients
u = fwd_solve(img);
ctx.u   = u.volt;                     % [n_nodes x num_stim]
ctx.GxU = ctx.G.Gx*ctx.u;             % [n_elem  x num_stim]
ctx.GyU = ctx.G.Gy*ctx.u;
ctx.GzU = ctx.G.Gz*ctx.u;             % identically zero

% Factor the CEM Schur-complement operator ONCE. It is SPD up to a 1-D
% constant null space (the floating potential): ground the last node.
Amat = M(img,ctx.sigma);
Amat = (Amat + Amat.')/2;             % remove roundoff asymmetry for LDL
ctx.free = 1:size(Amat,1)-1;
ctx.Afac = decomposition(Amat(ctx.free,ctx.free),'ldl');

% Mesh quadrature geometry (depends only on the mesh): triangle quadrature
% points (numQuadPts x 2 x numElements) and element areas.
[ctx.xi,ctx.area,ctx.weights] = precompute_quad_geometry_2d(fmdl);
end

%% Solve Amat*X = RHS for all columns at once (grounded node -> up to a const)
function X = ctx_solve(ctx,RHS)
X = zeros(size(RHS,1),size(RHS,2));
X(ctx.free,:) = ctx.Afac \ RHS(ctx.free,:);
end

%% z-channel Jacobian dBz/dsigma (adjoint method), sigma frozen via CTX
% Assumes identity sensor axes (z measurement axis = e_z), as in this study.
function J = calc_Jz(ctx,sensor_locations)
G = ctx.G;  n_elem = ctx.n_elem;  num_stim = ctx.num_stim;
numSensors = size(sensor_locations,1);

R = compute_R_2d(ctx,sensor_locations);

Sigma = spdiags(ctx.sigma,0,n_elem,n_elem);
Cz     = -R.Ry*Sigma*G.Gx + R.Rx*Sigma*G.Gy;   % [numSensors x n_nodes]
Gamma3 = ctx.mu_factor*Cz;

lambda = ctx_solve(ctx,-Gamma3.');             % [n_nodes x numSensors]

GxL = reshape(G.Gx*lambda,[n_elem 1 numSensors]);
GyL = reshape(G.Gy*lambda,[n_elem 1 numSensors]);
GzL = reshape(G.Gz*lambda,[n_elem 1 numSensors]);   % zero

GxU = reshape(ctx.GxU,[n_elem num_stim 1]);
GyU = reshape(ctx.GyU,[n_elem num_stim 1]);
GzU = reshape(ctx.GzU,[n_elem num_stim 1]);         % zero

elemV = reshape(ctx.elemV,[n_elem 1 1]);
dfdx  = elemV .* ( GxL.*GxU + GyL.*GyU + GzL.*GzU );

Rx_ = reshape(R.Rx.',[n_elem 1 numSensors]);
Ry_ = reshape(R.Ry.',[n_elem 1 numSensors]);
dCzdp = ( -Ry_.*GxU + Rx_.*GyU );              % z measurement axis
dfdp  = ctx.mu_factor*dCzdp;

dfd = dfdx + dfdp;                             % [n_elem x num_stim x numSensors]
dfd = permute(dfd,[3 2 1]);                    % [numSensors x num_stim x n_elem]
J   = reshape(dfd,numSensors*num_stim,n_elem);
end

%% dphi/dp accumulation for the z measurement axis (contracts W with dJ/dp)
function dphidp = grad_z(ctx,sensor_locations,W)
G = ctx.G;  n_elem = ctx.n_elem;  num_stim = ctx.num_stim; %#ok<NASGU>
numSensors = size(sensor_locations,1);
sigma = ctx.sigma;  mu = ctx.mu_factor;
GxU = ctx.GxU;  GyU = ctx.GyU;  GzU = ctx.GzU;    % [n_elem x num_stim]
elemV = ctx.elemV(:).';                           % [1 x n_elem]
block_size = numSensors*num_stim;

% Biot-Savart kernel derivatives (once per j): dim=3 needs dR_y and dR_x
dRy = compute_dR_2d(ctx,sensor_locations,2);
dRx = compute_dR_2d(ctx,sensor_locations,1);

dphidp = zeros(3,numSensors);
for p = 1:3
    % adjoint-derivative source (all sensors), then a single block solve
    rhs = mu*( (dRy{p}.*sigma.')*G.Gx - (dRx{p}.*sigma.')*G.Gy );  % [nS x n_nodes]
    dlambda = ctx_solve(ctx,rhs.');               % [n_nodes x numSensors]
    dlGx = G.Gx*dlambda;                          % [n_elem x numSensors]
    dlGy = G.Gy*dlambda;
    dlGz = G.Gz*dlambda;                          % zero
    for m = 1:numSensors
        ids = m:numSensors:block_size;
        Wm  = W(ids,:);                           % [num_stim x n_elem]
        tmp = dlGx(:,m).*GxU + dlGy(:,m).*GyU + dlGz(:,m).*GzU;   % [n_elem x num_stim]
        dJ1 = tmp.' .* elemV;                     % [num_stim x n_elem]
        dJ2 = -mu*( dRx{p}(m,:).*GyU.' - dRy{p}(m,:).*GxU.' );    % [num_stim x n_elem]
        dphidp(p,m) = -2*sum(sum(Wm.*(dJ1 - dJ2)));
    end
end
end

%% Biot-Savart element integrals R for the current sensor positions (CTX geom)
% 2D triangular mesh, flat at z=0: quadrature points xi are numQuadPts x 2
% x numElements (x,y only); the sensor's out-of-plane height is added
% explicitly as a constant third coordinate (dz = sensor z, since the
% current sheet sits at z=0).
function R = compute_R_2d(ctx,sensor_locations)
xi = ctx.xi;  area = ctx.area;
numElements = size(xi,3);
numSensors  = size(sensor_locations,1);
numQuadPts  = size(xi,1);

Rx = zeros(numSensors,numElements);
Ry = zeros(numSensors,numElements);
Rz = zeros(numSensors,numElements);
w  = reshape(ctx.weights,[numQuadPts,1,1]);

for m = 1:numSensors
    rm = sensor_locations(m,:);
    dxy = rm(1:2) - xi;                            % numQuadPts x 2 x numElements
    dz  = rm(3)*ones(numQuadPts,1,numElements);     % quad points are at z=0
    dm  = cat(2,dxy,dz);                            % numQuadPts x 3 x numElements
    nrm3 = sum(dm.^2,2).^(3/2);
    f = dm ./ nrm3;
    integ = sum(f.*w,1);                          % 1 x 3 x numElements
    Rx(m,:) = squeeze(integ(1,1,:))'.*area;
    Ry(m,:) = squeeze(integ(1,2,:))'.*area;
    Rz(m,:) = squeeze(integ(1,3,:))'.*area;
end
R.Rx = Rx;  R.Ry = Ry;  R.Rz = Rz;
end

%% Derivatives dR_j/dp of the Biot-Savart integrals (p = x,y,z), CTX geometry
% Same closed-form derivative as the 3D case: it only depends on the
% 3-component vector from the (2D, z=0) quadrature point to the (3D)
% sensor, not on the dimensionality of the source mesh.
function dR = compute_dR_2d(ctx,sensor_locations,j)
xi = ctx.xi;  area = ctx.area;
numElements = size(xi,3);
numSensors  = size(sensor_locations,1);
numQuadPts  = size(xi,1);

dRdx = zeros(numSensors,numElements);
dRdy = zeros(numSensors,numElements);
dRdz = zeros(numSensors,numElements);
w = reshape(ctx.weights,[numQuadPts,1,1]);

for m = 1:numSensors
    rm = sensor_locations(m,:);
    dxy = rm(1:2) - xi;
    dz  = rm(3)*ones(numQuadPts,1,numElements);
    dm_vec = cat(2,dxy,dz);                       % numQuadPts x 3 x numElements
    dm_norm2 = sum(dm_vec.^2,2);
    dm_j = dm_vec(:,j,:);
    for p = 1:3
        dm_p = dm_vec(:,p,:);
        funvals = ((j==p)*dm_norm2 - 3*dm_j.*dm_p) ./ (dm_norm2.^(5/2));
        vals = squeeze(sum(funvals.*w,1))'.*area;
        switch p
            case 1, dRdx(m,:) = vals;
            case 2, dRdy(m,:) = vals;
            case 3, dRdz(m,:) = vals;
        end
    end
end
dR = {dRdx,dRdy,dRdz};
end

%% Mesh-only quadrature geometry (quad points xi, element areas) -- triangles
function [xi,area,weights] = precompute_quad_geometry_2d(fmdl)
[coord,weights] = quad_rule_tri37();

numElements = size(fmdl.elems,1);
numQuadPts  = length(weights);

V = reshape(fmdl.nodes(fmdl.elems',:),[3,numElements,2]);
v1 = squeeze(V(1,:,:)); v2 = squeeze(V(2,:,:)); v3 = squeeze(V(3,:,:));

a = v2-v1; b = v3-v1;
area = 0.5*abs(a(:,1).*b(:,2) - a(:,2).*b(:,1))';   % 1 x numElements

Jm = zeros(2,2,numElements);
Jm(:,1,:) = permute(a,[2 3 1]);
Jm(:,2,:) = permute(b,[2 3 1]);

xi = reshape(v1',[2,1,numElements]) + ...
    Jm(:,1,:).*reshape(coord(:,1),[1,numQuadPts,1]) + ...
    Jm(:,2,:).*reshape(coord(:,2),[1,numQuadPts,1]);
xi = permute(xi,[2 1 3]);                        % numQuadPts x 2 x numElements
end

%% Diagonal of the posterior covariance for a given configuration (z-only)
function post_var = compute_posterior_variance_diag(ctx,theta,rs,zs,...
    Gamma_prior,Gamma_noise)

J = calc_Jz(ctx,theta_to_locations(theta,rs,zs));

P = J*Gamma_prior;                          % [nd x n]
S = Gamma_noise + P*J.';                     % [nd x nd]
Y = S \ P;                                   % S^-1 P
post_var = diag(Gamma_prior) - sum(P.*Y,1)'; % diag(H^-1)
end

%% Finite-difference gradient check (central differences)
function check_gradient_fd(costgrad,theta0,n_check,h)
n = numel(theta0);
[~,g] = costgrad(theta0);
idx = randperm(n,min(n_check,n));
for i = idx
    e = zeros(n,1); e(i) = 1;
    fp = costgrad(theta0 + h*e);
    fm = costgrad(theta0 - h*e);
    g_fd = (fp - fm)/(2*h);
    rel_err = abs(g(i)-g_fd)/max(abs(g_fd),eps);
    fprintf('  dphi/dtheta(%i): analytic = %+.6e | FD = %+.6e | rel err = %.2e\n',...
        i,g(i),g_fd,rel_err);
end
end

%% Schur complement of the CEM system matrix (SPD up to a constant null space)
function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac - Ae*(Ad\Ae');
end

%% Biot-Savart element integrals (used by assign_sensor_locations, arbitrary
%% 2D fmdl -- same quadrature as compute_R_2d, but computed fresh from an
%% arbitrary fmdl rather than the precomputed CTX geometry)
function [Rx,Ry,Rz,fmdl] = compute_r_matrices_local(fmdl,sensor_locations)

[coord,weights] = quad_rule_tri37();

numElements = size(fmdl.elems,1);
numSensors  = size(sensor_locations,1);
numQuadPts  = length(weights);

V = reshape(fmdl.nodes(fmdl.elems',:),[3,numElements,2]);
v1 = squeeze(V(1,:,:)); v2 = squeeze(V(2,:,:)); v3 = squeeze(V(3,:,:));

a = v2-v1; b = v3-v1;
area = 0.5*abs(a(:,1).*b(:,2) - a(:,2).*b(:,1))';   % 1 x numElements

Jm = zeros(2,2,numElements);
Jm(:,1,:) = permute(a,[2 3 1]);
Jm(:,2,:) = permute(b,[2 3 1]);

xi = reshape(v1',[2,1,numElements]) + ...
    Jm(:,1,:).*reshape(coord(:,1),[1,numQuadPts,1]) + ...
    Jm(:,2,:).*reshape(coord(:,2),[1,numQuadPts,1]);
xi = permute(xi,[2 1 3]);                % numQuadPts x 2 x numElements

Rx = zeros(numSensors,numElements);
Ry = zeros(numSensors,numElements);
Rz = zeros(numSensors,numElements);

w = reshape(weights,[numQuadPts,1,1]);

for m = 1:numSensors
    rm = sensor_locations(m,:);
    dxy = rm(1:2) - xi;
    dz  = rm(3)*ones(numQuadPts,1,numElements);
    dm  = cat(2,dxy,dz);                 % numQuadPts x 3 x numElements
    nrm3 = sum(dm.^2,2).^(3/2);
    f = dm ./ nrm3;
    integ = sum(f.*w,1);                 % 1 x 3 x numElements
    Rx(m,:) = squeeze(integ(1,1,:))'.*area;
    Ry(m,:) = squeeze(integ(1,2,:))'.*area;
    Rz(m,:) = squeeze(integ(1,3,:))'.*area;
end

fmdl.R.Rx = Rx;
fmdl.R.Ry = Ry;
fmdl.R.Rz = Rz;

end

%% Pairwise sensor repulsion (optional penalty)
function Grep = repulsion_cost(sensor_locations)

n_sensors = size(sensor_locations,1);

X = sensor_locations(:,1);
Y = sensor_locations(:,2);
Z = sensor_locations(:,3);

DX = X - X.'; DY = Y - Y.'; DZ = Z - Z.';

D2 = DX.^2 + DY.^2 + DZ.^2 + 1e-12;

invD2 = 1 ./ D2;
invD2(1:n_sensors+1:end) = 0;

Grep = 0.5*sum(invD2(:));
end

function [dGx,dGy,dGz] = repulsion_grad_cartesian(sensor_locations)

n_sensors = size(sensor_locations,1);

X = sensor_locations(:,1);
Y = sensor_locations(:,2);
Z = sensor_locations(:,3);

DX = X - X.'; DY = Y - Y.'; DZ = Z - Z.';

D2 = DX.^2 + DY.^2 + DZ.^2 + 1e-12;

invD2 = 1 ./ D2;
invD2(1:n_sensors+1:end) = 0;

C = -2*invD2.^2;
C(1:n_sensors+1:end) = 0;

dGx = sum(C.*DX,2);
dGy = sum(C.*DY,2);
dGz = sum(C.*DZ,2);
end

%% Assign sensor locations to the forward model
function img = assign_sensor_locations(img,sensor_locations)
assert(numel(img.fwd_model.sensors) == size(sensor_locations,1));
for id = 1:numel(img.fwd_model.sensors)
    img.fwd_model.sensors(id).position = sensor_locations(id,:);
end

[Rx,Ry,Rz,fmdl] = compute_r_matrices_local(img.fwd_model,sensor_locations);

img.fwd_model = fmdl;
end

%% 37-point TOMS706 quadrature rule on the reference triangle (0,0)-(1,0)-(0,1)
% Same rule used by EIDORS's own 2D Biot-Savart quadrature, see
% functions/fwd_model/compute_r_matrices.m (get_toms_37).
function [coord,weights] = quad_rule_tri37()
coord = [
    0.333333333333333333333333333333  0.333333333333333333333333333333;
    0.950275662924105565450352089520  0.024862168537947217274823955239;
    0.024862168537947217274823955239  0.950275662924105565450352089520;
    0.024862168537947217274823955239  0.024862168537947217274823955239;
    0.171614914923835347556304795551  0.414192542538082326221847602214;
    0.414192542538082326221847602214  0.171614914923835347556304795551;
    0.414192542538082326221847602214  0.414192542538082326221847602214;
    0.539412243677190440263092985511  0.230293878161404779868453507244;
    0.230293878161404779868453507244  0.539412243677190440263092985511;
    0.230293878161404779868453507244  0.230293878161404779868453507244;
    0.772160036676532561750285570113  0.113919981661733719124857214943;
    0.113919981661733719124857214943  0.772160036676532561750285570113;
    0.113919981661733719124857214943  0.113919981661733719124857214943;
    0.009085399949835353883572964740  0.495457300025082323058213517632;
    0.495457300025082323058213517632  0.009085399949835353883572964740;
    0.495457300025082323058213517632  0.495457300025082323058213517632;
    0.062277290305886993497083640527  0.468861354847056503251458179727;
    0.468861354847056503251458179727  0.062277290305886993497083640527;
    0.468861354847056503251458179727  0.468861354847056503251458179727;
    0.022076289653624405142446876931  0.851306504174348550389457672223;
    0.022076289653624405142446876931  0.126617206172027096933163647918;
    0.851306504174348550389457672223  0.022076289653624405142446876931;
    0.851306504174348550389457672223  0.126617206172027096933163647918;
    0.126617206172027096933163647918  0.022076289653624405142446876931;
    0.126617206172027096933163647918  0.851306504174348550389457672223;
    0.018620522802520968955913511549  0.689441970728591295496647976487;
    0.018620522802520968955913511549  0.291937506468887771754472382212;
    0.689441970728591295496647976487  0.018620522802520968955913511549;
    0.689441970728591295496647976487  0.291937506468887771754472382212;
    0.291937506468887771754472382212  0.018620522802520968955913511549;
    0.291937506468887771754472382212  0.689441970728591295496647976487;
    0.096506481292159228736516560903  0.635867859433872768286976979827;
    0.096506481292159228736516560903  0.267625659273967961282458816185;
    0.635867859433872768286976979827  0.096506481292159228736516560903;
    0.635867859433872768286976979827  0.267625659273967961282458816185;
    0.267625659273967961282458816185  0.096506481292159228736516560903;
    0.267625659273967961282458816185  0.635867859433872768286976979827];

weights  = [
    0.051739766065744133555179145422
    0.008007799555564801597804123460
    0.008007799555564801597804123460
    0.008007799555564801597804123460
    0.046868898981821644823226732071
    0.046868898981821644823226732071
    0.046868898981821644823226732071
    0.046590940183976487960361770070
    0.046590940183976487960361770070
    0.046590940183976487960361770070
    0.031016943313796381407646220131
    0.031016943313796381407646220131
    0.031016943313796381407646220131
    0.010791612736631273623178240136
    0.010791612736631273623178240136
    0.010791612736631273623178240136
    0.032195534242431618819414482205
    0.032195534242431618819414482205
    0.032195534242431618819414482205
    0.015445834210701583817692900053
    0.015445834210701583817692900053
    0.015445834210701583817692900053
    0.015445834210701583817692900053
    0.015445834210701583817692900053
    0.015445834210701583817692900053
    0.017822989923178661888748319485
    0.017822989923178661888748319485
    0.017822989923178661888748319485
    0.017822989923178661888748319485
    0.017822989923178661888748319485
    0.017822989923178661888748319485
    0.037038683681384627918546472190
    0.037038683681384627918546472190
    0.037038683681384627918546472190
    0.037038683681384627918546472190
    0.037038683681384627918546472190
    0.037038683681384627918546472190];

end
