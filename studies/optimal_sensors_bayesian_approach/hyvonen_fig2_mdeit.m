%% MDEIT analogue of Hyvonen, Seppanen & Staboulis (2014) Fig. 5.1 ("Fig. 2")
%
% Compute script: 4-electrode / 4-magnetometer flat 2D disk, Bz-only
% MDEIT, brute-force grid search vs. gradient-based (fminunc) optimization
% of the magnetometer azimuths ONLY (electrodes stay evenly spaced -- see
% the plan header for why: EIDORS bakes electrodes into the Netgen mesh,
% so a joint search would require remeshing every iterate).
%
% This is a COPY of brute_force_4x4_2d.m with four changes:
%   1. Hyvonen smoothness prior (5.1)-(5.2) instead of the diagonal prior
%      (see PLAN doc; with a diagonal prior the A-opt reduction is ~1%
%      and all six panels of the figure would look identical).
%   2. Two dense-prior correctness fixes in costgrad_theta_z, required
%      once the prior is no longer diagonal:
%        (a) A-opt gradient weight Yd = Y*Gamma_prior (was Y.*d.')
%        (b) D-opt cost = logdet(Gpr) + logdet(Gn) - logdet(S)
%            (was sum(log(diag(Gamma_noise))) - sum(log(d)) - logdet_S,
%            which has both a wrong sign and a diag-instead-of-logdet bug
%            -- a theta-independent error that never moved the minimizer
%            but corrupted the reported "nats" numbers)
%        (c) D-opt repulsion: alpha_rep = 0 for d-opt (phi_data < 0 there,
%            so the multiplicative penalty used for a-opt would instead
%            REWARD sensor clustering)
%   3. Separate reconstruction/design mesh at radius/15 (not the radius/30
%      forward mesh with the anomaly baked in -- see PLAN doc "mesh
%      sizing" for the inverse-crime / cost-blowup reasons).
%   4. Runs BOTH 'a-opt' and 'd-opt' through the same sweep, saving
%      everything build_hyvonen_fig2_mdeit.m needs for the 2x3 figure.
%
% Run "quick_test = true; hyvonen_fig2_mdeit" for a fast smoke test.
% Requires MATLAB R2024a (Optimization Toolbox / fminunc).

%% Prepare workspace
fullpath = mfilename('fullpath');
script_folder = fileparts(fullpath);
cd(script_folder);

grandparent_folder = fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));
addpath(fullfile(grandparent_folder,'studies','hyvonen2014_electrode_optimization','lib'));

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
                                   % height above the plane (zs = height/2).

model_parameters.numOfElectrodesPerRing = 4;
model_parameters.numOfRings = 1;

% Conductive circular (disk) anomaly, off-center (same placement as
% example_anomaly_circle_2d_smoothprior.m's default sw_anom_offset=0.5)
model_parameters.material.name   = 'disk';
model_parameters.material.type   = 'cylindrical';
model_parameters.material.radius = model_parameters.radius/5;
model_parameters.material.position = ...
    [-0.5*model_parameters.radius/sqrt(2), -0.5*model_parameters.radius/sqrt(2), 0];

% Sensors: M evenly spaced on a circle of radius rs, at height zs above
% the (flat, z=0) domain
n_sensors = 3;
rs = 1.05*model_parameters.radius;        % circle radius (fixed, close to domain)
zs = 0;           % sensor height is 0 in 2D

theta_even = 2*pi*(0:n_sensors-1)'/n_sensors;

model_parameters.numOfSensors = n_sensors;
model_parameters.sensorRadius = rs;
model_parameters.sensorPositions = theta_to_locations(theta_even,rs,zs);

background_conductivity = 1.0;
anomaly_conductivity = 5*background_conductivity;  % conductive anomaly

%% Simulation parameters
% opt_modes = {'a-opt','d-opt'};

rng(1);
opt_modes = {'a-opt'};

n_starts = 4;                  % 1 = start from the even configuration only


d_target = 4;                  % desired # of data-dominated modes

prior_std_background = 0.03;
prior_std_roi        = 0.4;
roi_radius_factor    = 1.0;    % ROI disk radius = factor * anomaly radius
prior_lambda = 1.0*model_parameters.material.radius;

% D-opt: phi_data < 0 (logdet of a small covariance), so the multiplicative
% repulsion penalty phi*(1+alpha*Grep) used for a-opt would REWARD
% clustering there -- disable it for d-opt (see header).
alpha_rep_by_mode = struct('a_opt',1e-5,'d_opt',0);

if ~exist('grid_n','var'), grid_n = 24; end   % angle grid resolution (15 deg steps)
max_iterations = 30;                          % for the fminunc comparison runs

if quick_test
    max_iterations = 3;
    grid_n = 10;                % coarse/fast grid for a smoke test
end

%% Make forward models (2D: fmdl.nodes is Nx2, fmdl.elems is Nx3 triangles)
[~,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_fwd = fmdls{1};

model_parameters_recon = model_parameters;
model_parameters_recon.maxsz = model_parameters_recon.radius/15;
% model_parameters_recon = rmfield(model_parameters_recon,'anomaly');
% model_parameters_recon = rmfield(model_parameters_recon,'material');
[~,fmdls] = mk_mdeit_model(model_parameters_recon,model_folder);
fmdl_recon = fmdls{1};

% Homogeneous image = prior mean (the OED linearization point)
imgh_fwd = mk_image_mdeit(fmdl_fwd,background_conductivity);
imgh_recon = mk_image_mdeit(fmdl_recon,background_conductivity);

% Image with the anomaly (ground truth, used only for noise calibration)
imgi_fwd = add_material_properties(imgh_fwd,[background_conductivity,anomaly_conductivity]);

n_stim  = numel(fmdl_recon.stimulation);
n_elem  = size(fmdl_recon.elems,1);
n_nodes = size(fmdl_recon.nodes,1);

fprintf('Model (2D): %i elements, %i nodes, %i stims, %i sensors (Bz only)\n',...
    n_elem,n_nodes,n_stim,n_sensors);
fprintf('Recon mesh: maxsz = radius/15 (%i elements)\n',n_elem);

%% Prior covariance: Hyvonen smoothness prior (5.1)-(5.2)
%   Gamma_ij = kappa_i*kappa_j * exp(-|x_i-x_j|^2/(2*lambda^2))
% with kappa = prior_std_roi inside the ROI, prior_std_background outside,
% and ZERO cross-block covariance. Built blockwise (dense) so the zero
% cross-block costs nothing.
anomaly_center = model_parameters.material.position;   % [x,y,~] (z unused in 2D)
roi_radius = roi_radius_factor*model_parameters.material.radius;

centroid_dist = sqrt(sum((fmdl_recon.elem_centroids - anomaly_center(1:2)).^2,2));
roi = centroid_dist <= roi_radius;

kappa_vec = prior_std_background*ones(n_elem,1);
kappa_vec(roi) = prior_std_roi;

Gamma_prior = zeros(n_elem,n_elem);
for blk = {roi, ~roi}
    ib = find(blk{1});
    xb = fmdl_recon.elem_centroids(ib,:);
    D2 = sum(xb.^2,2) + sum(xb.^2,2).' - 2*(xb*xb.');
    Gamma_prior(ib,ib) = (kappa_vec(ib)*kappa_vec(ib).') .* exp(-max(D2,0)/(2*prior_lambda^2));
end
Gamma_prior = (Gamma_prior+Gamma_prior.')/2;

prior_variance = diag(Gamma_prior);

mu_pr = sort(eig(Gamma_prior),'descend'); mu_pr = max(mu_pr,0);
fprintf('ROI: %i of %i elements (%.1f%% of prior trace)\n',nnz(roi),n_elem,...
    100*sum(prior_variance(roi))/sum(prior_variance));
fprintf('Smoothness prior: lambda = %.4f (= %.2f x anomaly radius), effective rank = %.1f\n',...
    prior_lambda, prior_lambda/model_parameters.material.radius, sum(mu_pr)^2/sum(mu_pr.^2));
fprintf('A-opt reduction CEILING (top-%i eigenvalues / trace) = %.2f%%\n',...
    d_target, 100*sum(mu_pr(1:min(d_target,n_elem)))/sum(mu_pr));

%% Noise covariance (white), calibrated from the prior-whitened Jacobian spectrum
imgh_fwd = assign_sensor_locations(imgh_fwd,theta_to_locations(theta_even,rs,zs));
Bh = fwd_solve_mdeit(imgh_fwd);

imgi_fwd = assign_sensor_locations(imgi_fwd,theta_to_locations(theta_even,rs,zs));
Bi = fwd_solve_mdeit(imgi_fwd);

dB_even = Bi.Bz(:) - Bh.Bz(:); %#ok<NASGU>
max_B = max(abs(Bh.Bz(:)));

% Precompute everything that does NOT depend on the sensor positions
ctx = build_ctx(imgh_recon);

J0 = calc_Jz(ctx,theta_to_locations(theta_even,rs,zs));
% Whiten with the FULL prior square root: the diagonal shortcut
% "J0 .* sqrt(prior_variance)'" is only correct for a diagonal prior; with
% the dense smoothness prior it would badly miscalibrate d_target.
L_pr = chol(Gamma_prior + 1e-12*trace(Gamma_prior)/n_elem*eye(n_elem),'lower');
[~,S,~] = svd(J0*L_pr,'econ');   % prior-whitened spectrum
s = diag(S);

noise_std = s(d_target);
variance_noise = noise_std^2;

n_data = n_stim*n_sensors;                 % z-channel only
Gamma_noise = (variance_noise)*speye(n_data,n_data);

fprintf('noise std = %.3e (= max|B|/%.1f)\n',noise_std,max_B/noise_std);

%% Cost/gradient handles (theta is the only free variable), per opt_mode
costgrad = cell(1,numel(opt_modes));
for iom = 1:numel(opt_modes)
    alpha_rep = alpha_rep_by_mode.(strrep(opt_modes{iom},'-','_'));
    costgrad{iom} = @(theta) costgrad_theta_z(ctx,theta,rs,zs,...
        Gamma_prior,Gamma_noise,opt_modes{iom},alpha_rep);
end

%% Gradient check, BOTH modes -- gate everything else on this passing
fprintf('\n=== Finite-difference gradient checks ===\n');
for iom = 1:numel(opt_modes)
    fprintf('\nGradient check (%s):\n',opt_modes{iom});
    check_gradient_fd(costgrad{iom},theta_even + 0.1*randn(n_sensors,1),4,1e-4);
end

%% Independent cross-check against lib/oed_criterion.m (Step 2 fix verification)
% NB: oed_criterion computes the pure OED criterion with no repulsion
% penalty, so the comparison here uses alpha_rep = 0 regardless of mode
% (costgrad_theta_z's phi for a-opt otherwise includes the
% phi_data*(1+alpha_rep*Grep) penalty and would not match).
fprintf('\n=== Cross-check vs. hyvonen2014_electrode_optimization/lib/oed_criterion.m ===\n');
J_even = calc_Jz(ctx,theta_to_locations(theta_even,rs,zs));
for iom = 1:numel(opt_modes)
    phi_here = costgrad_theta_z(ctx,theta_even,rs,zs,Gamma_prior,Gamma_noise,opt_modes{iom},0);
    phi_ref  = oed_criterion(J_even,Gamma_prior,Gamma_noise,opt_modes{iom});
    fprintf('  %s: phi(here) = %.10e | phi(oed_criterion) = %.10e | rel diff = %.2e\n',...
        opt_modes{iom},phi_here,phi_ref,abs(phi_here-phi_ref)/abs(phi_ref));
end

%% Posterior variance consistency check (A-opt = trace of H^-1)
post_var_even_check = compute_posterior_variance_diag(ctx,theta_even,rs,zs,Gamma_prior,Gamma_noise);
phi_a_even_direct = costgrad_theta_z(ctx,theta_even,rs,zs,Gamma_prior,Gamma_noise,'a-opt',0);
fprintf('\nPosterior variance consistency: sum(post_var_even) = %.10e | phi_a(even) = %.10e | rel diff = %.2e\n',...
    sum(post_var_even_check),phi_a_even_direct,...
    abs(sum(post_var_even_check)-phi_a_even_direct)/abs(phi_a_even_direct));

%% Baseline: evenly spaced sensors
phi_even = zeros(1,numel(opt_modes));
for iom = 1:numel(opt_modes)
    phi_even(iom) = costgrad{iom}(theta_even);
end

%% Gradient-based (local) optimization, starting from the even configuration
options = optimoptions('fminunc',...
    'Algorithm','quasi-newton','SpecifyObjectiveGradient',true,...
    'Display','iter','MaxIterations',max_iterations,...
    'OptimalityTolerance',1e-7,'StepTolerance',1e-8);

theta_fminunc = cell(n_starts,numel(opt_modes));
phi_fminunc = zeros(n_starts,numel(opt_modes));
for iom = 1:numel(opt_modes)
    fprintf('\n=== Gradient-based optimization (fminunc, from even start, %s) ===\n',opt_modes{iom});

    for istart = 1:n_starts
        fprintf('\n=== Start %i ===\n',istart);

        if istart == 1
            th0 = theta_even;
        else
            th0 = 2*pi*rand(n_sensors,1);
        end

        tstart = tic;
        [theta_fminunc{istart,iom},phi_fminunc(istart,iom)] = fminunc(costgrad{iom},th0,options);
        fprintf('phi(fminunc) = %.6e (%.1f s)\n',phi_fminunc(iom),toc(tstart));
    end

end

%% Keep only best result

theta_fminunc_final = cell(1,numel(opt_modes));
phi_fminunc_final = zeros(1,numel(opt_modes));

for iom = 1:numel(opt_modes)
    %Find minimum
    id = find(phi_fminunc(:,iom) == min(phi_fminunc(:,iom)));
    % Assin
    phi_fminunc_final(iom) = phi_fminunc(id,iom);
    theta_fminunc_final{iom} = theta_fminunc{id,iom};
end

theta_fminunc = theta_fminunc_final;
phi_fminunc =  phi_fminunc_final;

%% Brute-force exhaustive grid search over sensor angles (both modes, same sweep)
% Sensors are unordered -> nchoosek(1:grid_n,n_sensors) combinations
% (not grid_n^n_sensors permutations) cover every distinct configuration
% at this grid resolution.
angle_grid = 2*pi*(0:grid_n-1)/grid_n;
combos = nchoosek(1:grid_n,n_sensors);
n_combos = size(combos,1);

fprintf('\n=== Brute force: %i grid angles (%.1f deg step), %i combinations x %i modes ===\n',...
    grid_n,360/grid_n,n_combos,numel(opt_modes));

phi_grid = inf(n_combos,numel(opt_modes));
tstart = tic;
for k = 1:n_combos
    theta_k = angle_grid(combos(k,:))';

    % Shared Jacobian/Cholesky factor: evaluate both criteria from the
    % same factorization (marginal cost of D-opt is one extra logdet).
    sensor_locations = theta_to_locations(theta_k,rs,zs);
    J = calc_Jz(ctx,sensor_locations);
    P  = J*Gamma_prior;
    Sm = Gamma_noise + P*J.';
    Sm = (Sm+Sm.')/2;
    Ls = chol(Sm + 1e-12*trace(Sm)/n_data*speye(n_data),'lower');
    Y  = Ls.' \ (Ls \ P);
    d  = diag(Gamma_prior);

    phi_grid(k,1) = sum(d) - sum(sum(P.*Y));                              % a-opt
    logdet_S = 2*sum(log(diag(Ls)));
    phi_grid(k,2) = logdet_Gpr_const(L_pr) + sum(log(diag(Gamma_noise))) - logdet_S;  % d-opt

    if mod(k,floor(n_combos/10)) == 0
        fprintf('  %i/%i combinations (%.1f s elapsed)\n',k,n_combos,toc(tstart));
    end
end
fprintf('Brute force done: %i combinations x %i modes in %.1f s\n',n_combos,numel(opt_modes),toc(tstart));

[phi_bf,k_best] = min(phi_grid,[],1);
theta_bf = cell(1,numel(opt_modes));
for iom = 1:numel(opt_modes)
    theta_bf{iom} = angle_grid(combos(k_best(iom),:))';
end

%% Polish the best grid point with a local search (removes grid discretization)
% The grid is only 15 deg coarse; polishing with fminunc removes the
% discretization error and is what makes the "brute force" column a fair
% estimate of the true global minimum rather than of the best grid vertex.
theta_bf_polished = cell(1,numel(opt_modes));
phi_bf_polished = zeros(1,numel(opt_modes));
for iom = 1:numel(opt_modes)
    fprintf('\n=== Local polish of the brute-force optimum (%s) ===\n',opt_modes{iom});
    [theta_bf_polished{iom},phi_bf_polished(iom)] = fminunc(costgrad{iom},theta_bf{iom},options);
    assert(phi_bf_polished(iom) <= phi_bf(iom) + 1e-8*abs(phi_bf(iom)),...
        'fminunc polish diverged for %s: phi_bf_polished > phi_bf',opt_modes{iom});
end

%% Posterior variance diagonals for every panel of the figure
post_var_even = compute_posterior_variance_diag(ctx,theta_even,rs,zs,Gamma_prior,Gamma_noise);

post_var_bf_polished = cell(1,numel(opt_modes));
post_var_fminunc = cell(1,numel(opt_modes));
for iom = 1:numel(opt_modes)
    post_var_bf_polished{iom} = compute_posterior_variance_diag(ctx,theta_bf_polished{iom},rs,zs,Gamma_prior,Gamma_noise);
    post_var_fminunc{iom} = compute_posterior_variance_diag(ctx,theta_fminunc{iom},rs,zs,Gamma_prior,Gamma_noise);

end

%% Pure (repulsion-free) tr(Gamma)/logdet(Gamma) for reporting and titling.
% phi_even/phi_fminunc/phi_bf/phi_bf_polished are the OPTIMIZED objective:
% for a-opt that includes the alpha_rep repulsion penalty
% (phi = phi_data*(1+alpha_rep*Grep)), which is what fminunc actually
% minimizes but is NOT tr(Gamma). Recompute the alpha_rep=0 value at each
% theta so the figure/report show the paper's actual criterion, not the
% penalized surrogate. (d-opt already runs with alpha_rep=0, so
% phi_pure(2,:) == phi(2,:) exactly.)
phi_pure_even     = zeros(1,numel(opt_modes));
phi_pure_fminunc  = zeros(1,numel(opt_modes));
phi_pure_bf       = zeros(1,numel(opt_modes));
phi_pure_bf_polished = zeros(1,numel(opt_modes));
for iom = 1:numel(opt_modes)
    om = opt_modes{iom};
    phi_pure_even(iom)        = costgrad_theta_z(ctx,theta_even,rs,zs,Gamma_prior,Gamma_noise,om,0);
    phi_pure_fminunc(iom)     = costgrad_theta_z(ctx,theta_fminunc{iom},rs,zs,Gamma_prior,Gamma_noise,om,0);
    phi_pure_bf(iom)          = costgrad_theta_z(ctx,theta_bf{iom},rs,zs,Gamma_prior,Gamma_noise,om,0);
    phi_pure_bf_polished(iom) = costgrad_theta_z(ctx,theta_bf_polished{iom},rs,zs,Gamma_prior,Gamma_noise,om,0);
end

%% Report (same style as hyvonen_case1.m)
fprintf('\n================ BRUTE FORCE vs GRADIENT-BASED RESULTS (MDEIT, 2D) ================\n');
fprintf('theta_even (deg) = %s\n',mat2str(mod(theta_even,2*pi)*180/pi,4));

for iom = 1:numel(opt_modes)
    om = opt_modes{iom};
    fprintf('\n--- %s ---\n',om);
    fprintf('  phi(even)                = %.6e  [pure: %.6e]\n',phi_even(iom),phi_pure_even(iom));
    fprintf('  phi(fminunc)              = %.6e  [pure: %.6e]\n',phi_fminunc(iom),phi_pure_fminunc(iom));
    fprintf('  phi(brute-force grid)     = %.6e  (%.1f deg grid)  [pure: %.6e]\n',phi_bf(iom),360/grid_n,phi_pure_bf(iom));
    fprintf('  phi(brute-force polished) = %.6e  [pure: %.6e]\n',phi_bf_polished(iom),phi_pure_bf_polished(iom));
    fprintf('  gap: fminunc vs brute-force-polished = %+.3e (%.4f%% of phi_even)\n',...
        phi_fminunc(iom)-phi_bf_polished(iom),100*(phi_fminunc(iom)-phi_bf_polished(iom))/phi_even(iom));
    switch om
        case 'a-opt'
            fprintf('  trace reduction (fminunc,  pure tr(Gamma))          = %.1f%%\n',...
                100*(phi_pure_even(iom)-phi_pure_fminunc(iom))/phi_pure_even(iom));
            fprintf('  trace reduction (brute-force pol., pure tr(Gamma))  = %.1f%%\n',...
                100*(phi_pure_even(iom)-phi_pure_bf_polished(iom))/phi_pure_even(iom));
        case 'd-opt'
            dnats_fmin = (phi_pure_even(iom)-phi_pure_fminunc(iom))/2;
            dnats_bf   = (phi_pure_even(iom)-phi_pure_bf_polished(iom))/2;
            fprintf('  info gain (fminunc)          = %.2f nats (%.2f bits)\n',dnats_fmin,dnats_fmin/log(2));
            fprintf('  info gain (brute-force pol.)  = %.2f nats (%.2f bits)\n',dnats_bf,dnats_bf/log(2));
    end
    fprintf('  theta_fminunc        (deg) = %s\n',mat2str(mod(theta_fminunc{iom},2*pi)*180/pi,4));
    fprintf('  theta_brute_polished (deg) = %s\n',mat2str(mod(theta_bf_polished{iom},2*pi)*180/pi,4));
    fprintf('  ROI posterior variance sums:\n');
    fprintf('    even        : %.4e\n',sum(post_var_even(roi)));
    fprintf('    fminunc     : %.4e  (ratio vs even = %.4f)\n',...
        sum(post_var_fminunc{iom}(roi)),sum(post_var_fminunc{iom}(roi))/sum(post_var_even(roi)));
    fprintf('    brute force : %.4e  (ratio vs even = %.4f)\n',...
        sum(post_var_bf_polished{iom}(roi)),sum(post_var_bf_polished{iom}(roi))/sum(post_var_even(roi)));
end

%% Save results
results = struct();
results.opt_modes = {opt_modes};
results.grid_n = grid_n;
results.n_starts = n_starts;
results.n_combos = n_combos;
results.theta_even = theta_even;
results.theta_fminunc = {theta_fminunc};
results.theta_bf = {theta_bf};
results.theta_bf_polished = {theta_bf_polished};
results.phi_even = phi_even;
results.phi_fminunc = phi_fminunc;
results.phi_bf = phi_bf;
results.phi_bf_polished = phi_bf_polished;
% "Pure" = alpha_rep=0, i.e. the actual tr(Gamma)/logdet(Gamma) criterion
% with no repulsion penalty (see comment above where these are computed).
% For d-opt these equal phi_* exactly (alpha_rep_by_mode.d_opt = 0).
results.phi_pure_even = phi_pure_even;
results.phi_pure_fminunc = phi_pure_fminunc;
results.phi_pure_bf = phi_pure_bf;
results.phi_pure_bf_polished = phi_pure_bf_polished;
results.phi_grid = phi_grid;
results.combos = combos;
results.angle_grid = angle_grid;
results.roi = roi;
results.prior_variance = prior_variance;
results.post_var_even = post_var_even;
results.post_var_bf_polished = {post_var_bf_polished};
results.post_var_fminunc = {post_var_fminunc};
results.Gamma_noise_scalar = variance_noise;
results.noise_std = noise_std;
results.model_parameters = model_parameters;
results.rs = rs;
results.zs = zs;

fmdl_recon_fig = struct();
fmdl_recon_fig.nodes = fmdl_recon.nodes;
fmdl_recon_fig.elems = fmdl_recon.elems;
fmdl_recon_fig.elem_centroids = fmdl_recon.elem_centroids;
fmdl_recon_fig.electrode_nodes = {fmdl_recon.electrode(:).nodes};
results.fmdl_recon = fmdl_recon_fig;

data_folder = fullfile(script_folder,'data');
if ~exist(data_folder,'dir'), mkdir(data_folder); end
save(fullfile(data_folder,'hyvonen_fig2_mdeit_results.mat'),'results');

fprintf('\nSaved results to data/hyvonen_fig2_mdeit_results.mat\n');
fprintf('Done.\n');

%% ======================== LOCAL FUNCTIONS ========================
% (identical to brute_force_4x4_2d.m / example_anomaly_circle_2d.m -- kept
% self-contained here so this script does not depend on locals that are
% not exported by those files -- except for costgrad_theta_z, which has
% the two dense-prior bug fixes described in the file header.)

%% Constant logdet(Gamma_prior), from the L_pr already computed for noise calib
function v = logdet_Gpr_const(L_pr)
v = 2*sum(log(diag(L_pr)));
end

%% Map circle angles to 3D sensor locations
function sensor_locations = theta_to_locations(theta,rs,zs)
theta = theta(:);
sensor_locations = [rs*cos(theta), rs*sin(theta), zs*ones(numel(theta),1)];
end

%% Cost and gradient in the angle parametrization (Bz-only MDEIT)
function [phi,dphi] = costgrad_theta_z(ctx,theta,rs,zs,Gamma_prior,Gamma_noise,opt_mode,alpha_rep)
assert(ismember(opt_mode,{'a-opt','d-opt'}));

sensor_locations = theta_to_locations(theta,rs,zs);

% z-channel Jacobian (shared by cost and gradient)
J = calc_Jz(ctx,sensor_locations);

% Data-space (Woodbury) form: never forms the n x n posterior covariance.
%   H^-1 = Gpr - Gpr J' S^-1 J Gpr,   S = Gn + J Gpr J'   (nd x nd, SPD)
P  = J * Gamma_prior;                  % [nd x n]
S  = Gamma_noise + P * J.';            % [nd x nd]
S  = (S+S.')/2;                        % roundoff asymmetry hurts chol/Y precision
n_data_local = size(S,1);
Ls = chol(S + 1e-12*trace(S)/n_data_local*speye(n_data_local),'lower');  % factor once (reused for Y and logdet)
Y  = Ls.' \ (Ls \ P);                  % S^-1 P  [nd x n]
d  = diag(Gamma_prior);                % prior variances [n x 1]

switch opt_mode
    case 'a-opt'
        phi_data = sum(d) - sum(sum(P .* Y));          % trace(H^-1)
    case 'd-opt'
        % logdet(Gamma_post) = logdet(Gpr) + logdet(Gn) - logdet(S).
        % (Fixes a bug present in brute_force_4x4_2d.m, which computed
        % sum(log(diag(Gamma_noise))) - sum(log(d)) - logdet_S: both a
        % wrong sign and diag-instead-of-logdet for a dense prior. The
        % sign error is theta-independent so it never moved the minimizer
        % or the gradient, but it corrupted the reported "nats" numbers.)
        logdet_S   = 2*sum(log(diag(Ls)));
        n_elem     = size(Gamma_prior,1);
        L_pr_local = chol(Gamma_prior + 1e-12*trace(Gamma_prior)/n_elem*eye(n_elem),'lower');
        logdet_Gpr = 2*sum(log(diag(L_pr_local)));
        phi_data   = logdet_Gpr + sum(log(diag(Gamma_noise))) - logdet_S;
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
        Yd = Y * Gamma_prior;         % Y*Gpr   [nd x n] -- NOT the diagonal
                                      % shortcut Y.*d.', which is only valid
                                      % for a diagonal prior.
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
ctx.elemV     = fmdl.elem_volume(:);  % element AREAS in 2D

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
% explicitly as a constant third coordinate.
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
