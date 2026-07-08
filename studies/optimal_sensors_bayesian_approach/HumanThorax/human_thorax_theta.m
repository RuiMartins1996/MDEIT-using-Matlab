%% TASK A: Bayesian optimization of sensor AZIMUTHS on the human thorax (MDEIT)
%
% Sensors are constrained to a circle of FIXED radius R0 at half the height of
% a 3D extruded adult-male thorax model. Only the azimuth theta_m of each
% sensor is free (radius R0 and height zs = height/2 are fixed). A diagonal
% "Hyvonen Case 1" prior encodes WHERE the anomaly may be: elevated std in a
% ball (the ROI) around a suspected off-centre location, small std in the
% background. Because the thorax cross-section is strongly non-circular AND the
% ROI is off-centre, the fixed-radius circle has an azimuth-dependent standoff
% from the body and a broken rotational symmetry, so the optimal azimuths are
% expected to be genuinely structured (unlike the near-rotation-invariant
% A-opt result on the cylinder).
%
% Design criteria (3-axis MDEIT, linearized around the prior mean = the
% homogeneous-lung thorax image):
%   A-optimality: minimize trace(H^-1),   H = J'*inv(Gn)*J + inv(Gpr)
%   D-optimality: minimize -log det(H)
%
% Optimization is driven by the SHARED, validated optimizer
% functions/optimize_sensor_configuration.m (bug-fixed dlambda indexing,
% correct Cholesky solve order, chain-rule Jacobian re-evaluated at the current
% iterate -- see HANDOFF_human_thorax.m Sec 2.1), fed the theta coordinate
% maps (map_q_to_x_theta / map_x_to_q_theta / jac_coord_transf_theta).
%
% This script ALSO contains self-contained cost/gradient helpers (copied from
% the FD-validated cylinder example example_anomaly_circle.m) used ONLY for:
%   (i)  the mandatory central finite-difference gradient check of the theta
%        chain rule (jac_coord_transf_theta), the arbiter per Sec 5;
%   (ii) computing phi(even)/phi(opt)/phi(boundary) and the posterior-variance
%        and roi_energy diagnostics for the report.
% Their phi(even) is cross-checked against the shared optimizer's own printed
% f0 (must match to all digits).
%
% Traps of human_thorax_example.m removed (Sec 2.3): the ROI prior is built
% unconditionally (NOT gated on the existence of a results file), the noise is
% set by the d_target calibration (NOT an ad-hoc max_B*1e-1), and all outputs
% use task-specific filenames (thorax_theta_*).
%
% Workspace-overridable knobs (set in the base workspace before running):
%   quick_test, max_iterations, n_starts, d_target,
%   prior_std_background_factor, prior_std_roi_factor, roi_radius_factor.
%
% Run a fast smoke test with:
%   quick_test = true; human_thorax_theta

%% Prepare workspace
clc; close all
fullpath = mfilename('fullpath');
script_folder = fileparts(fullpath);
cd(script_folder);

grandgrandparent_folder = fileparts(fileparts(fileparts(script_folder)));
addpath(genpath(fullfile(grandgrandparent_folder,'functions')));

model_folder = prepare_workspace(script_folder);
addpath(genpath(fullfile(script_folder,'model')));

rng(1);

if ~exist('quick_test','var'), quick_test = false; end

%% Characteristic units (as in human_thorax_example.m)
z_0 = 0.0058;   %(Ohm m^2) contact impedance (58 Ohm cm^2)
l0  = 40e-3;    %(m) tank radius (characteristic length)
I0  = 2.4e-3;   %(A) injected current magnitude

sigma0 = l0/z_0; %(S/m) characteristic conductivity

%% Geometry / model parameters (MUST match the cached model_*.mat, Sec 2.3.4)
height  = 0.5;
mu0     = 1.0;

num_sensors = 16;
R0 = 1.5;               % FIXED sensor-circle radius (Task A)
zs = height/2;          % FIXED sensor-circle height

inj  = [0 3];
meas = [0 3];
current_amplitude = 2.4e-3/I0;

num_of_electrodes_per_ring = 12;
electrode_radius = 0.15;
maxsz_mesh = 1;         % coarse mesh
maxsz_electrode = 0.15;

% Suspected anomaly / ROI centre (off-centre; azimuth = atan2(0.1,0.1) = 45 deg)
anomaly.position = [0.1, 0.1, zs];
anomaly.radius   = 0.1;

% Even (baseline == optimizer start) azimuths
theta_even = 2*pi*(0:num_sensors-1)'/num_sensors;
sensor_positions_0 = map_q_to_x_theta_locations(theta_even,R0,zs);

%% Simulation knobs (workspace-overridable)
opt_modes = {'a-opt','d-opt'};

if ~exist('max_iterations','var'), max_iterations = 40; end
if ~exist('n_starts','var'),       n_starts = 1;  end
do_gradient_check = false;

% Prior / noise regime (Sec 3). Cylinder final: bg=0.005, roi=1.0, d_target=100.
% roi_radius_factor raised 1.0 -> 2.0 for the thorax: on this coarse mesh a
% factor-1 ball captures only 14 elements (roi_energy pathologically low ~0.11);
% a factor-2 ball captures 129 elements -- comparable to the cylinder's ~100 ROI
% elements -- restoring roi_energy ~0.33 and a ~26% even-config ROI variance
% reduction (tuning sweep, see HANDOFF UPDATE). d_target=100 kept as on cylinder.
if ~exist('d_target','var'),                   d_target = 100;   end
if ~exist('prior_std_background_factor','var'),prior_std_background_factor = 0.005; end
if ~exist('prior_std_roi_factor','var'),       prior_std_roi_factor = 1.0; end
if ~exist('roi_radius_factor','var'),          roi_radius_factor = 2.0; end

delta_boundary = 0.05;  % standoff for the boundary ("as close as allowed") config

if quick_test
    max_iterations = 3;
    n_starts = 1;
    do_gradient_check = true;
end

%% Load (or build) the cached FEM models
% Homogeneous-lung model = OED linearization point (prior mean).
model_hom_file   = fullfile(script_folder,'model','model_homogeneous.mat');
model_recon_file = fullfile(script_folder,'model','model_recon.mat');

if exist(model_hom_file,'file')
    var = load(model_hom_file);  fmdl_hom = var.fmdl;
else
    fmdl_hom = build_thorax_model(height,{'thorax','rlung','llung'},maxsz_mesh,...
        num_of_electrodes_per_ring,electrode_radius,maxsz_electrode,zs,mu0,sensor_positions_0);
    fmdl = fmdl_hom; save(model_hom_file,'fmdl');
end

if exist(model_recon_file,'file')
    var = load(model_recon_file);  fmdl_recon = var.fmdl;
else
    fmdl_recon = build_thorax_model(height,{'thorax'},maxsz_mesh,...
        num_of_electrodes_per_ring,electrode_radius,maxsz_electrode,zs,mu0,sensor_positions_0);
    fmdl = fmdl_recon; save(model_recon_file,'fmdl');
end

%% Stimulation
stimulation = mk_stim_patterns(num_of_electrodes_per_ring,1,inj,meas,{'meas_current'},current_amplitude);
fmdl_hom.stimulation   = stimulation;
fmdl_recon.stimulation = stimulation;

%% Images
% Prior-mean image (background = 1.0, lungs = 0.3), 3 material regions
imgh = mk_image_mdeit(fmdl_hom,1.0);
imgh = add_material_properties(imgh,[1.0,0.3,0.3]);
imgh = assign_sensor_locations(imgh,sensor_positions_0);

background_sigma = 1.0;  % background elem_data value (prior std scales with this)

n_stim  = numel(fmdl_hom.stimulation);
n_elem  = size(fmdl_hom.elems,1);
n_nodes = size(fmdl_hom.nodes,1);
fprintf('Model: %i elements, %i nodes, %i stims, %i sensors (3-axis)\n',...
    n_elem,n_nodes,n_stim,num_sensors);

A = @(x) M(imgh,x);

% Jacobian at the even (start) configuration -- independent of prior/noise, so
% computed once and reused by the noise calibration, the roi_energy diagnostic
% and the (optional) tuning sweep.
img0 = assign_sensor_locations(imgh,sensor_positions_0);
J0 = calc_jacobian_3axis_local(img0,A);

centroids = fmdl_hom.elem_centroids;
centroid_dist = sqrt(sum((centroids - anomaly.position).^2,2));

%% Optional tuning sweep (Sec 3): pick (roi_radius_factor, bg, d_target)
% Set tune_sweep=true in the base workspace to print the diagnostics table and
% exit BEFORE optimizing (roi_energy + even-config ROI variance reduction).
if exist('tune_sweep','var') && tune_sweep
    if ~exist('rr_list','var'), rr_list = [1 2 3 4]; end
    if ~exist('bg_list','var'), bg_list = [0.005 0.02]; end
    if ~exist('dt_list','var'), dt_list = [20 50 100]; end
    fprintf('\n===== TUNING SWEEP (roi_std_factor=%.2f) =====\n',prior_std_roi_factor);
    fprintf('roi_rf  bg_fac  nROI  d_tgt  noise_std   d_modes  roi_energy  evenROI/prior\n');
    for rr = rr_list
        roi_s = centroid_dist <= rr*anomaly.radius;
        nroi = nnz(roi_s);
        for bg = bg_list
            pv = (bg*background_sigma)^2*ones(n_elem,1);
            pv(roi_s) = (prior_std_roi_factor*background_sigma)^2;
            [~,Ss,Vss] = svd(J0 .* sqrt(pv)','econ');  ss = diag(Ss);
            for dt = dt_list
                dtc = min(dt,numel(ss));
                vn = ss(dtc)^2;
                dm = sum(ss.^2/vn > 1);
                re = mean(sum(Vss(roi_s,1:max(dm,1)).^2,1));
                Hs = full(J0.'*(J0/vn) + spdiags(1./pv,0,n_elem,n_elem));
                Ls = chol(Hs,'lower'); Xs = Ls\speye(n_elem);
                pvar = sum(Xs.^2,1)';
                red = sum(pvar(roi_s))/sum(pv(roi_s));
                fprintf('%5.1f  %6.3f  %4i  %5i  %9.3e  %6i    %6.3f      %6.3f\n',...
                    rr,bg,nroi,dtc,sqrt(vn),dm,re,red);
            end
        end
    end
    fprintf('===== END TUNING SWEEP =====\n');
    return
end

%% Prior covariance: diagonal, informative about WHERE the anomaly may be
% (Hyvonen Case 1). Built UNCONDITIONALLY (trap 2.3.1 removed).
roi_radius = roi_radius_factor*anomaly.radius;
roi = centroid_dist <= roi_radius;

prior_std_background = prior_std_background_factor*background_sigma;
prior_std_roi        = prior_std_roi_factor*background_sigma;

prior_variance = prior_std_background^2*ones(n_elem,1);
prior_variance(roi) = prior_std_roi^2;

Gamma_prior     = spdiags(prior_variance,0,n_elem,n_elem);
inv_Gamma_prior = spdiags(1./prior_variance,0,n_elem,n_elem);

fprintf('ROI: %i of %i elements (%.1f%% of prior trace)\n',nnz(roi),n_elem,...
    100*sum(prior_variance(roi))/sum(prior_variance));

%% Noise covariance (white), calibrated from the whitened Jacobian spectrum
% noise_std = s(d_target), the d_target-th singular value of J0*sqrt(Gpr) at the
% even configuration -> ~d_target modes are data-dominated (trap 2.3.2 removed).
% J0 was computed once above; reuse it here.
[~,S,Vs] = svd(J0 .* sqrt(prior_variance)','econ');
s = diag(S);

d_target = min(d_target,numel(s));
noise_std = s(d_target);
variance_noise = noise_std^2;

n_data = 3*n_stim*num_sensors;
inv_Gamma_noise = (1/variance_noise)*speye(n_data,n_data);

d_modes = sum(s.^2/variance_noise > 1);
roi_energy = mean(sum(Vs(roi,1:max(d_modes,1)).^2,1));
fprintf('noise std = %.3e\n',noise_std);
fprintf('whitened spectrum: s1^2/var = %.2e, #data-dominated = %i of %i data / %i params\n',...
    s(1)^2/variance_noise,d_modes,n_data,n_elem);
fprintf('mean ROI energy of the %i data-dominated modes: %.2f (want close to 1)\n',...
    d_modes,roi_energy);

%% Coordinate-transform handles for the SHARED optimizer (theta parametrization)
q_to_x_theta          = @(q) map_q_to_x_theta(q,R0,zs);
x_to_q_theta          = @(x) map_x_to_q_theta(x);
jac_coord_transf_theta= @(s) compute_jac_coord_transf_theta(s);

% Round-trip assert (Sec 4.2): x0 -> q -> x0
x0 = sensor_locations_to_vector_cartesian(sensor_positions_0);
assert(all(abs(x0 - q_to_x_theta(x_to_q_theta(x0))) < 1e-9),...
    'theta coordinate round-trip failed');

%% Self-contained cost/grad handles (for FD check, phi values, diagnostics)
costgrad = cell(1,numel(opt_modes));
for iom = 1:numel(opt_modes)
    costgrad{iom} = @(theta) costgrad_theta(imgh,theta,R0,zs,...
        inv_Gamma_prior,inv_Gamma_noise,A,opt_modes{iom},jac_coord_transf_theta);
end

%% Mandatory FD gradient check (Sec 5): even start AND a perturbed config
if do_gradient_check
    idx_pool = unique([1 round(num_sensors/2) num_sensors ...
                       randperm(num_sensors,min(3,num_sensors))]);
    for iom = 1:numel(opt_modes)
        fprintf('\n=== Gradient check (%s), even start ===\n',opt_modes{iom});
        check_gradient_fd(costgrad{iom},theta_even,idx_pool,1e-4);
        fprintf('=== Gradient check (%s), perturbed config ===\n',opt_modes{iom});
        check_gradient_fd(costgrad{iom},theta_even+0.3*randn(num_sensors,1),idx_pool,1e-4);
    end
end

%% Baseline phi(even)
phi_even = zeros(1,numel(opt_modes));
for iom = 1:numel(opt_modes)
    phi_even(iom) = costgrad{iom}(theta_even);
    fprintf('phi(even, %s) = %.10e\n',opt_modes{iom},phi_even(iom));
end

%% Boundary ("as close as allowed") baseline config
sensor_positions_boundary = find_sensor_positions_boundary(fmdl_hom,[0,0,zs],delta_boundary);
phi_boundary = zeros(1,numel(opt_modes));
for iom = 1:numel(opt_modes)
    phi_boundary(iom) = cost_at_positions(imgh,sensor_positions_boundary,...
        inv_Gamma_prior,inv_Gamma_noise,A,opt_modes{iom});
end

%% Parallel pool + EIDORS on workers (Sec 2.4)
pool = gcp('nocreate');
if isempty(pool), pool = parpool('local'); end
pctRunOnAll(sprintf('setupEidors(''%s'');',script_folder));

%% Optimize with the SHARED optimizer (theta maps), both criteria + multistart
options = optimoptions('fminunc',...
    'Algorithm','quasi-newton','SpecifyObjectiveGradient',true,...
    'Display','iter','MaxIterations',max_iterations,...
    'OptimalityTolerance',1e-8,'StepTolerance',1e-12,'UseParallel',true);

theta_opt = cell(1,numel(opt_modes));
phi_opt   = inf(1,numel(opt_modes));
theta_starts = cell(1,numel(opt_modes));   % all multistart results

for iom = 1:numel(opt_modes)
    fprintf('\n================ Optimizing %s (theta, R0=%.2f) ================\n',...
        opt_modes{iom},R0);
    theta_starts{iom} = nan(num_sensors,n_starts);

    for istart = 1:n_starts
        if istart == 1
            th0 = theta_even;
        else
            th0 = 2*pi*rand(num_sensors,1);
        end
        pos0 = map_q_to_x_theta_locations(th0,R0,zs);
        imgstart = assign_sensor_locations(imgh,pos0);

        tstart = tic;
        img_out = optimize_sensor_configuration(imgstart,...
            inv_Gamma_prior,inv_Gamma_noise,...
            jac_coord_transf_theta,q_to_x_theta,x_to_q_theta,...
            opt_modes{iom},'mdeit3',3,options);

        pos_out = fetch_sensor_locations(img_out);
        th = atan2(pos_out(:,2),pos_out(:,1));
        fval = cost_at_positions(imgh,pos_out,inv_Gamma_prior,inv_Gamma_noise,A,opt_modes{iom});
        fprintf('start %i: phi = %.10e (%.1f s)\n',istart,fval,toc(tstart));

        theta_starts{iom}(:,istart) = th;
        if fval < phi_opt(iom)
            phi_opt(iom) = fval;
            theta_opt{iom} = th;
        end

        % Checkpoint after every (mode,start): a crash/shutdown cannot waste
        % the completed optimizations (deterministic, rng(1)).
        ckpt = struct('opt_modes',{opt_modes},'theta_even',theta_even,...
            'theta_opt',{theta_opt},'phi_opt',phi_opt,'phi_even',phi_even,...
            'theta_starts',{theta_starts},'done_mode',iom,'done_start',istart,...
            'R0',R0,'zs',zs);
        save(fullfile(script_folder,'data','thorax_theta_checkpoint.mat'),'ckpt');
    end
end

%% Report
fprintf('\n================ RESULTS (Task A: theta free, R0=%.2f) ================\n',R0);
for iom = 1:numel(opt_modes)
    fprintf('\n%s:\n',opt_modes{iom});
    fprintf('  phi(even)     = %.6e\n',phi_even(iom));
    fprintf('  phi(opt)      = %.6e\n',phi_opt(iom));
    fprintf('  phi(boundary) = %.6e\n',phi_boundary(iom));
    switch opt_modes{iom}
        case 'a-opt'
            fprintf('  trace reduction vs even     = %.2f%%\n',100*(1-phi_opt(iom)/phi_even(iom)));
            fprintf('  trace reduction vs boundary = %.2f%%\n',100*(1-phi_opt(iom)/phi_boundary(iom)));
        case 'd-opt'
            fprintf('  info gain vs even     = %.3f nats (%.3f bits)\n',...
                phi_even(iom)-phi_opt(iom),(phi_even(iom)-phi_opt(iom))/log(2));
            fprintf('  info gain vs boundary = %.3f nats (%.3f bits)\n',...
                phi_boundary(iom)-phi_opt(iom),(phi_boundary(iom)-phi_opt(iom))/log(2));
    end
end

%% Posterior-variance diagnostics (prior vs even vs a-opt)
post_var_prior = prior_variance;
post_var_even  = compute_posterior_variance_diag(imgh,map_q_to_x_theta_locations(theta_even,R0,zs),...
    inv_Gamma_prior,inv_Gamma_noise,A);
post_var_aopt  = compute_posterior_variance_diag(imgh,map_q_to_x_theta_locations(theta_opt{1},R0,zs),...
    inv_Gamma_prior,inv_Gamma_noise,A);

fprintf('\nPosterior variance sums (ROI | background):\n');
fprintf('  prior      : %.3e | %.3e\n',sum(post_var_prior(roi)),sum(post_var_prior(~roi)));
fprintf('  even       : %.3e | %.3e\n',sum(post_var_even(roi)),sum(post_var_even(~roi)));
fprintf('  a-optimized: %.3e | %.3e\n',sum(post_var_aopt(roi)),sum(post_var_aopt(~roi)));
fprintf('  ROI posterior variance: opt/even = %.4f, even/prior = %.4f\n',...
    sum(post_var_aopt(roi))/sum(post_var_even(roi)),...
    sum(post_var_even(roi))/sum(post_var_prior(roi)));

%% Save results
results = struct();
results.theta_even   = theta_even;
results.theta_opt    = theta_opt;
results.theta_starts = theta_starts;
results.opt_modes    = {opt_modes};
results.phi_even     = phi_even;
results.phi_opt      = phi_opt;
results.phi_boundary = phi_boundary;
results.sensor_positions_boundary = sensor_positions_boundary;
results.post_var_prior = post_var_prior;
results.post_var_even  = post_var_even;
results.post_var_aopt  = post_var_aopt;
results.prior_variance = prior_variance;
results.roi = roi;
results.R0 = R0; results.zs = zs; results.height = height;
results.anomaly = anomaly;
results.d_target = d_target; results.noise_std = noise_std;
results.d_modes = d_modes; results.roi_energy = roi_energy;
results.prior_std_background_factor = prior_std_background_factor;
results.prior_std_roi_factor = prior_std_roi_factor;
results.roi_radius_factor = roi_radius_factor;
save(fullfile(script_folder,'data','thorax_theta_results.mat'),'results');

%% Figures
make_figures(script_folder,imgh,R0,zs,height,theta_even,theta_opt,...
    sensor_positions_boundary,anomaly,fmdl_hom,roi,...
    post_var_prior,post_var_even,post_var_aopt);

fprintf('\nDone (Task A).\n');

%% ======================== LOCAL FUNCTIONS ========================

%% ---- Coordinate maps (theta free, radius R0 fixed, height zs fixed) ----
function locations = map_q_to_x_theta_locations(theta,R0,zs)
theta = theta(:);
locations = [R0*cos(theta), R0*sin(theta), zs*ones(numel(theta),1)];
end

function x = map_q_to_x_theta(q,R0,zs)
% q = theta (M x 1) -> x = [x1..xM; y1..yM; z1..zM] (coordinate-major, Sec 2.2)
q = q(:);
x = [R0*cos(q); R0*sin(q); zs*ones(numel(q),1)];
end

function q = map_x_to_q_theta(x)
assert(mod(numel(x),3)==0);
M = numel(x)/3;
xm = x(1:M); ym = x(M+1:2*M);
q = atan2(ym,xm);
q = q(:);
end

function Jt = compute_jac_coord_transf_theta(sensor_locations)
% Returns qmax x 3 x M with entry (1,dim,m) = d p_dim(m)/d theta_m.
% r_m, theta_m are recomputed from the CURRENT sensor_locations (Sec 2.1).
M = size(sensor_locations,1);
rm = sqrt(sensor_locations(:,1).^2 + sensor_locations(:,2).^2);
thm = atan2(sensor_locations(:,2),sensor_locations(:,1));
Jt = zeros(1,3,M);
Jt(1,1,:) = -rm.*sin(thm);
Jt(1,2,:) =  rm.*cos(thm);
Jt(1,3,:) =  0;
end

function x = sensor_locations_to_vector_cartesian(sensor_locations)
M = size(sensor_locations,1);
x = [sensor_locations(:,1); sensor_locations(:,2); sensor_locations(:,3)];
end

%% ---- Self-contained cost + gradient in theta (FD-validated machinery) ----
function [phi,dphi] = costgrad_theta(img,theta,R0,zs,...
    inv_Gamma_prior,inv_Gamma_noise,A,opt_mode,jac_coord_transf)

assert(ismember(opt_mode,{'a-opt','d-opt'}));
n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);
assert(numel(theta)==n_sensors);

sensor_locations = map_q_to_x_theta_locations(theta,R0,zs);
img = assign_sensor_locations(img,sensor_locations);

J = calc_jacobian_3axis_local(img,A);

H = full(J.'*inv_Gamma_noise*J + inv_Gamma_prior);
L = chol(H,'lower');

switch opt_mode
    case 'a-opt'
        X = L\eye(n_elem);
        phi = sum(X(:).^2);          % trace(H^-1)
    case 'd-opt'
        phi = -2*sum(log(diag(L)));  % -log det(H)
end

if nargout < 2, return; end

switch opt_mode
    case 'a-opt'
        invH = L'\X; Z = L\invH; Hpow = L'\Z;   % H^-2
    case 'd-opt'
        X = L\eye(n_elem); Hpow = L'\X;         % H^-1
end

W = inv_Gamma_noise*(J*Hpow);
dJ = compute_dJxyz_xyz(img);

dphidp = zeros(3,n_sensors);
block_size = n_sensors*n_stim;
for p = 1:3
    for m = 1:n_sensors
        ids_local = m : n_sensors : block_size;
        acc = 0;
        for dim = 1:3
            ids = ids_local + (dim-1)*block_size;
            dJm = reshape(dJ{dim,p}(m,:,:),[n_stim,n_elem]);
            Wm  = W(ids,:);
            acc = acc + sum(Wm(:).*dJm(:));
        end
        dphidp(p,m) = -2*acc;
    end
end

% Chain rule via the SAME handle passed to the shared optimizer (validates it)
Jt = jac_coord_transf(sensor_locations);   % 1 x 3 x M
dphi = zeros(n_sensors,1);
for m = 1:n_sensors
    dphi(m) = Jt(1,1,m)*dphidp(1,m) + Jt(1,2,m)*dphidp(2,m) + Jt(1,3,m)*dphidp(3,m);
end
end

%% ---- Cost at arbitrary sensor positions (for even/opt/boundary phi) ----
function phi = cost_at_positions(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,opt_mode)
n_elem = size(img.fwd_model.elems,1);
img = assign_sensor_locations(img,sensor_locations);
J = calc_jacobian_3axis_local(img,A);
H = full(J.'*inv_Gamma_noise*J + inv_Gamma_prior);
L = chol(H,'lower');
switch opt_mode
    case 'a-opt'
        X = L\eye(n_elem); phi = sum(X(:).^2);
    case 'd-opt'
        phi = -2*sum(log(diag(L)));
end
end

%% ---- Posterior covariance diagonal ----
function post_var = compute_posterior_variance_diag(img,sensor_locations,...
    inv_Gamma_prior,inv_Gamma_noise,A)
n_elem = size(img.fwd_model.elems,1);
img = assign_sensor_locations(img,sensor_locations);
J = calc_jacobian_3axis_local(img,A);
H = full(J.'*inv_Gamma_noise*J + inv_Gamma_prior);
L = chol(H,'lower');
X = L\eye(n_elem);
post_var = sum(X.^2,1)';
end

%% ---- Central finite-difference gradient check ----
function check_gradient_fd(costgrad,theta0,idx,h)
n = numel(theta0);
[~,g] = costgrad(theta0);
idx = idx(idx>=1 & idx<=n);
for i = idx(:)'
    e = zeros(n,1); e(i) = 1;
    fp = costgrad(theta0 + h*e);
    fm = costgrad(theta0 - h*e);
    g_fd = (fp - fm)/(2*h);
    rel_err = abs(g(i)-g_fd)/max(abs(g_fd),eps);
    fprintf('  dphi/dtheta(%2i): analytic = %+.6e | FD = %+.6e | rel err = %.2e\n',...
        i,g(i),g_fd,rel_err);
end
end

%% ---- Stacked 3-axis Jacobian (recomputes R from sensor positions) ----
function J = calc_jacobian_3axis_local(img,A)
Jx = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,1);
Jy = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,2);
Jz = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,3);
J = [Jx;Jy;Jz];
end

function out = M(img,sigma)
numNodes = size(img.fwd_model.nodes,1);
img.elem_data = sigma;
s_mat = system_mat_1st_order(img);
Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);
out = Ac-Ae*inv(Ad)*Ae';
end

function img = compute_gamma_matrices_local(img)
mu_factor = img.fwd_model.mu0/(4*pi);
num_sensors = numel(img.fwd_model.sensors);
sensor_locations = zeros(num_sensors,3);
for i = 1:num_sensors
    sensor_locations(i,:) = img.fwd_model.sensors(i).position;
end
[Rx,Ry,Rz,fmdl] = compute_r_matrices_local(img.fwd_model,sensor_locations);
img.fwd_model = fmdl;   % CRUCIAL: keep R in sync with sensor positions
R.Rx = Rx; R.Ry = Ry; R.Rz = Rz;
G = img.fwd_model.G;
Sigma = spdiags(img.elem_data(:),0,length(img.elem_data),length(img.elem_data));
g = zeros(num_sensors,3,3);
for m = 1:num_sensors
    g(m,:,:) = [img.fwd_model.sensors(m).axes.axis1;
                img.fwd_model.sensors(m).axes.axis2;
                img.fwd_model.sensors(m).axes.axis3];
end
Cx = ( -R.Rz * Sigma * G.Gy +  R.Ry * Sigma * G.Gz );
Cy = ( -R.Rx * Sigma * G.Gz +  R.Rz * Sigma * G.Gx );
Cz = ( -R.Ry * Sigma * G.Gx +  R.Rx * Sigma * G.Gy );
img.Gamma1 = mu_factor*(g(:,1,1).*Cx + g(:,1,2).*Cy + g(:,1,3).*Cz);
img.Gamma2 = mu_factor*(g(:,2,1).*Cx + g(:,2,2).*Cy + g(:,2,3).*Cz);
img.Gamma3 = mu_factor*(g(:,3,1).*Cx + g(:,3,2).*Cy + g(:,3,3).*Cz);
end

function [Rx,Ry,Rz,fmdl] = compute_r_matrices_local(fmdl,sensor_locations)
[coord,weights] = quad_rule_35();
numElements = size(fmdl.elems,1);
numSensors  = size(sensor_locations,1);
numQuadPts  = length(weights);
V = reshape(fmdl.nodes(fmdl.elems',:),[4,numElements,3]);
v1 = squeeze(V(1,:,:)); v2 = squeeze(V(2,:,:));
v3 = squeeze(V(3,:,:)); v4 = squeeze(V(4,:,:));
a = v2-v1; b = v3-v1; c = v4-v1;
detJ = abs(sum(a .* cross(b,c,2),2))';
Jm = zeros(3,3,numElements);
Jm(:,1,:) = permute(a,[2 3 1]);
Jm(:,2,:) = permute(b,[2 3 1]);
Jm(:,3,:) = permute(c,[2 3 1]);
xi = reshape(v1',[3,1,numElements]) + ...
    Jm(:,1,:).*reshape(coord(:,1),[1,numQuadPts,1]) + ...
    Jm(:,2,:).*reshape(coord(:,2),[1,numQuadPts,1]) + ...
    Jm(:,3,:).*reshape(coord(:,3),[1,numQuadPts,1]);
xi = permute(xi,[2 1 3]);
Rx = zeros(numSensors,numElements);
Ry = zeros(numSensors,numElements);
Rz = zeros(numSensors,numElements);
w = reshape(weights,[numQuadPts,1,1]);
for m = 1:numSensors
    rm = sensor_locations(m,:);
    dm = rm - xi;
    nrm3 = sum(dm.^2,2).^(3/2);
    f = dm ./ nrm3;
    integ = sum(f.*w,1);
    Rx(m,:) = squeeze(integ(1,1,:))'.*(detJ/6);
    Ry(m,:) = squeeze(integ(1,2,:))'.*(detJ/6);
    Rz(m,:) = squeeze(integ(1,3,:))'.*(detJ/6);
end
fmdl.R.Rx = Rx; fmdl.R.Ry = Ry; fmdl.R.Rz = Rz;
end

function dR = dRmkj_xyz(fmdl,j)
[coord,weights] = quad_rule_35();
numElements = size(fmdl.elems,1);
numSensors  = numel(fmdl.sensors);
numQuadPts  = length(weights);
V = reshape(fmdl.nodes(fmdl.elems',:),[4,numElements,3]);
v1 = squeeze(V(1,:,:)); v2 = squeeze(V(2,:,:));
v3 = squeeze(V(3,:,:)); v4 = squeeze(V(4,:,:));
a = v2-v1; b = v3-v1; c = v4-v1;
detJ = abs(sum(a .* cross(b,c,2),2))';
Jm = zeros(3,3,numElements);
Jm(:,1,:) = permute(a,[2 3 1]);
Jm(:,2,:) = permute(b,[2 3 1]);
Jm(:,3,:) = permute(c,[2 3 1]);
xi = reshape(v1',[3,1,numElements]) + ...
    Jm(:,1,:).*reshape(coord(:,1),[1,numQuadPts,1]) + ...
    Jm(:,2,:).*reshape(coord(:,2),[1,numQuadPts,1]) + ...
    Jm(:,3,:).*reshape(coord(:,3),[1,numQuadPts,1]);
xi = permute(xi,[2 1 3]);
dRdx = zeros(numSensors,numElements);
dRdy = zeros(numSensors,numElements);
dRdz = zeros(numSensors,numElements);
w = reshape(weights,[numQuadPts,1,1]);
for m = 1:numSensors
    rm = fmdl.sensors(m).position;
    dm_vec = rm - xi;
    dm_norm2 = sum(dm_vec.^2,2);
    dm_j = dm_vec(:,j,:);
    for p = 1:3
        dm_p = dm_vec(:,p,:);
        funvals = ((j==p)*dm_norm2 - 3*dm_j.*dm_p) ./ (dm_norm2.^(5/2));
        vals = squeeze(sum(funvals.*w,1))'.*(detJ/6);
        switch p
            case 1, dRdx(m,:) = vals;
            case 2, dRdy(m,:) = vals;
            case 3, dRdz(m,:) = vals;
        end
    end
end
dR = {dRdx,dRdy,dRdz};
end

function [dlambda,dR1t,dR2t] = compute_dlambdaxyz_xyz(img)
num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);
dlambda = cell(3,3); dR1t = cell(3,3); dR2t = cell(3,3);
mu_factor = img.fwd_model.mu0/(4*pi);
G = img.fwd_model.G;
sigma = img.elem_data;
A_matrix = M(img,sigma);
d = sqrt(diag(A_matrix));
Mfun = @(x) x ./ d; Nfun = @(x) x ./ d;
for dim = 1:3
    switch dim
        case 1, dR1 = dRmkj_xyz(img.fwd_model,3); dR2 = dRmkj_xyz(img.fwd_model,2); G1 = G.Gy; G2 = G.Gz;
        case 2, dR1 = dRmkj_xyz(img.fwd_model,1); dR2 = dRmkj_xyz(img.fwd_model,3); G1 = G.Gz; G2 = G.Gx;
        case 3, dR1 = dRmkj_xyz(img.fwd_model,2); dR2 = dRmkj_xyz(img.fwd_model,1); G1 = G.Gx; G2 = G.Gy;
    end
    rhs1 = mu_factor*( (dR1{1} .* sigma.')*G1 - (dR2{1} .* sigma.')*G2 );
    rhs2 = mu_factor*( (dR1{2} .* sigma.')*G1 - (dR2{2} .* sigma.')*G2 );
    rhs3 = mu_factor*( (dR1{3} .* sigma.')*G1 - (dR2{3} .* sigma.')*G2 );
    dlambdaX = zeros(n_nodes,num_sensors);
    dlambdaY = zeros(n_nodes,num_sensors);
    dlambdaZ = zeros(n_nodes,num_sensors);
    parfor m = 1:num_sensors
        [dlambdaX(:,m),~,~] = pcg(A_matrix,rhs1(m,:)',1e-9,numel(sigma),Mfun,Nfun);
        [dlambdaY(:,m),~,~] = pcg(A_matrix,rhs2(m,:)',1e-9,numel(sigma),Mfun,Nfun);
        [dlambdaZ(:,m),~,~] = pcg(A_matrix,rhs3(m,:)',1e-9,numel(sigma),Mfun,Nfun);
    end
    dlambda{dim,1} = dlambdaX; dlambda{dim,2} = dlambdaY; dlambda{dim,3} = dlambdaZ;
    dR1t{dim,1} = dR1{1}; dR1t{dim,2} = dR1{2}; dR1t{dim,3} = dR1{3};
    dR2t{dim,1} = dR2{1}; dR2t{dim,2} = dR2{2}; dR2t{dim,3} = dR2{3};
end
end

function dJ = compute_dJxyz_xyz(img)
num_sensors = numel(img.fwd_model.sensors);
n_elem      = size(img.fwd_model.elems,1);
num_stim    = numel(img.fwd_model.stimulation);
mu_factor = img.fwd_model.mu0/(4*pi);
G = img.fwd_model.G;
elemV = img.fwd_model.elem_volume(:);
u = fwd_solve(img); u = u.volt;
GxU = G.Gx*u; GyU = G.Gy*u; GzU = G.Gz*u;
[dlambda,dR1dp,dR2dp] = compute_dlambdaxyz_xyz(img);
G1U = cell(1,3); G2U = cell(1,3);
for dim = 1:3
    switch dim
        case 1, G1U{dim} = GyU; G2U{dim} = GzU;
        case 2, G1U{dim} = GzU; G2U{dim} = GxU;
        case 3, G1U{dim} = GxU; G2U{dim} = GyU;
    end
end
dJ = cell(3,3);
dJ1 = zeros(num_sensors,num_stim,n_elem);
dJ2 = zeros(num_sensors,num_stim,n_elem);
for dim = 1:3
    for p = 1:3
        for m = 1:num_sensors
            dlGx = G.Gx*dlambda{dim,p}(:,m);
            dlGy = G.Gy*dlambda{dim,p}(:,m);
            dlGz = G.Gz*dlambda{dim,p}(:,m);
            tmp = dlGx.*GxU + dlGy.*GyU + dlGz.*GzU;
            dJ1(m,:,:) = tmp.' .* elemV(:).';
            dJ2(m,:,:) = -mu_factor * ( ...
                dR2dp{dim,p}(m,:) .* G2U{dim}.' - dR1dp{dim,p}(m,:) .* G1U{dim}.' );
        end
        dJ{dim,p} = dJ1 - dJ2;
    end
end
end

function J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,select_sensor_axis)
mu0 = img.fwd_model.mu0;
n_elem = size(img.fwd_model.elems,1);
num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);
img = compute_gamma_matrices_local(img);
R = img.fwd_model.R;
G = img.fwd_model.G;
switch select_sensor_axis
    case 1, Gamma = img.Gamma1;
    case 2, Gamma = img.Gamma2;
    case 3, Gamma = img.Gamma3;
    otherwise, error('Invalid sensor axis');
end
u = fwd_solve(img); u = u.volt;
A_matrix = A(img.elem_data);
GammaT = Gamma.';
lambda = A_matrix \ (-GammaT);
Gx_times_lambda = G.Gx*lambda;
Gy_times_lambda = G.Gy*lambda;
Gz_times_lambda = G.Gz*lambda;
Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;
mu_factor = mu0/(4*pi);
elemV = reshape(img.fwd_model.elem_volume(:),[1 1 n_elem]);
GxL = reshape(Gx_times_lambda.',[num_sensors 1 n_elem]);
GyL = reshape(Gy_times_lambda.',[num_sensors 1 n_elem]);
GzL = reshape(Gz_times_lambda.',[num_sensors 1 n_elem]);
GxU = reshape(Gx_times_u.',[1 num_stim n_elem]);
GyU = reshape(Gy_times_u.',[1 num_stim n_elem]);
GzU = reshape(Gz_times_u.',[1 num_stim n_elem]);
dfdx = elemV .* ( GxL.*GxU + GyL.*GyU + GzL.*GzU );
Rx_ = reshape(R.Rx,[num_sensors 1 n_elem]);
Ry_ = reshape(R.Ry,[num_sensors 1 n_elem]);
Rz_ = reshape(R.Rz,[num_sensors 1 n_elem]);
dCxdp = ( -Rz_.*GyU + Ry_.*GzU );
dCydp = ( -Rx_.*GzU + Rz_.*GxU );
dCzdp = ( -Ry_.*GxU + Rx_.*GyU );
g = zeros(num_sensors,3,3);
for m = 1:num_sensors
    g(m,:,:) = [img.fwd_model.sensors(m).axes.axis1;
                img.fwd_model.sensors(m).axes.axis2;
                img.fwd_model.sensors(m).axes.axis3];
end
gx = reshape(g(:,select_sensor_axis,1),[num_sensors 1 1]);
gy = reshape(g(:,select_sensor_axis,2),[num_sensors 1 1]);
gz = reshape(g(:,select_sensor_axis,3),[num_sensors 1 1]);
dfdp = mu_factor*( gx.*dCxdp + gy.*dCydp + gz.*dCzdp );
dfd = dfdx + dfdp;
J = reshape(dfd,num_sensors*num_stim,n_elem);
end

%% ---- Sensor bookkeeping ----
function img = assign_sensor_locations(img,sensor_locations)
assert(numel(img.fwd_model.sensors) == size(sensor_locations,1));
for id = 1:numel(img.fwd_model.sensors)
    img.fwd_model.sensors(id).position = sensor_locations(id,:);
end
end

function sensor_locations = fetch_sensor_locations(img)
n = numel(img.fwd_model.sensors);
sensor_locations = zeros(n,3);
for m = 1:n
    sensor_locations(m,:) = img.fwd_model.sensors(m).position;
end
end

function fmdl = assign_magnetometers(fmdl,sensor_positions)
num_sensors = size(sensor_positions,1);
sensor_axes = repmat(struct('axis1',[],'axis2',[],'axis3',[]),1,num_sensors);
for m = 1:num_sensors
    sensor_axes(m).axis1 = [1,0,0];
    sensor_axes(m).axis2 = [0,1,0];
    sensor_axes(m).axis3 = [0,0,1];
end
sensors = repmat(struct('position',[],'axes',[]),1,num_sensors);
for m = 1:num_sensors
    sensors(m).position = sensor_positions(m,:);
    sensors(m).axes = sensor_axes(m);
end
fmdl.sensors = sensors;
end

%% ---- Boundary ("as close as allowed") config ----
function sensor_position_boundary = find_sensor_positions_boundary(fmdl,ray_center,delta)
num_sensors = numel(fmdl.sensors);
dtheta = 2*pi/num_sensors;
thetas = 0:dtheta:2*pi-dtheta;
sensor_position_boundary = zeros(num_sensors,3);
for m = 1:numel(thetas)
    [hit_point,d] = find_hit_point_fmdl(fmdl,ray_center,thetas(m));
    sensor_position_boundary(m,:) = hit_point + delta*d;
end
end

function [hit_point,d] = find_hit_point_fmdl(fmdl,ray_center,theta)
assert(all(size(ray_center)==[1,3]))
boundary = fmdl.boundary;
d = [cos(theta), sin(theta), 0]; d = d / norm(d);
t_min = inf; hit_point = [];
for i = 1:size(boundary,1)
    tri = fmdl.nodes(boundary(i,:), :);
    v1 = tri(1,:); v2 = tri(2,:); v3 = tri(3,:);
    e1 = v2 - v1; e2 = v3 - v1;
    h = cross(d, e2); a = dot(e1, h);
    if abs(a) < 1e-12, continue; end
    f = 1.0 / a; s = ray_center - v1; u = f * dot(s, h);
    if u < 0 || u > 1, continue; end
    q = cross(s, e1); v = f * dot(d, q);
    if v < 0 || (u + v) > 1, continue; end
    t = f * dot(e2, q);
    if t > 0 && t < t_min
        hit_point = ray_center + t * d; break
    end
end
if isempty(hit_point), hit_point = ray_center; end
end

%% ---- Model builder (only used if a cache is missing) ----
function fmdl = build_thorax_model(height,contour_names,maxsz_mesh,...
    num_elec,electrode_radius,maxsz_electrode,zplane,mu0,sensor_positions_0)
contours = cell(1,numel(contour_names));
for i = 1:numel(contour_names)
    switch contour_names{i}
        case 'thorax', contours{i} = shape_library('get','adult_male','boundary');
        case 'rlung',  contours{i} = shape_library('get','adult_male','right_lung');
        case 'llung',  contours{i} = shape_library('get','adult_male','left_lung');
    end
end
shape = {height, contours, 1, maxsz_mesh};
elec_pos = [num_elec, 1, zplane]';
elec_shape = [electrode_radius, 0, maxsz_electrode]';
fmdl = ng_mk_extruded_model(shape, elec_pos, elec_shape);
fmdl.mu0 = mu0;
fmdl = assign_magnetometers(fmdl,sensor_positions_0);
fmdl = compute_geometry_matrices(fmdl,mu0);
end

%% ---- Figures ----
function make_figures(script_folder,imgh,R0,zs,height,theta_even,theta_opt,...
    sensor_positions_boundary,anomaly,fmdl_hom,roi,...
    post_var_prior,post_var_even,post_var_aopt) %#ok<INUSD>

pos_even = map_q_to_x_theta_locations(theta_even,R0,zs);
pos_a    = map_q_to_x_theta_locations(theta_opt{1},R0,zs);
pos_d    = map_q_to_x_theta_locations(theta_opt{2},R0,zs);
roi_azimuth = atan2(anomaly.position(2),anomaly.position(1));

% --- Sensor maps overlaid on the thorax slice (top view) ---
try
    figure('Name','Thorax theta: sensor configs','Position',[80 80 1300 350]);
    cfgs = {pos_even,pos_a,pos_d,sensor_positions_boundary};
    ttls = {'even','a-opt','d-opt','boundary'};
    for k = 1:4
        subplot(1,4,k); hold on
        himg = show_fem(imgh); set(findobj(himg,'Type','Patch'),'EdgeAlpha',0.08);
        plot3(cfgs{k}(:,1),cfgs{k}(:,2),cfgs{k}(:,3)+0.01,'r.','MarkerSize',16);
        plot3(anomaly.position(1),anomaly.position(2),anomaly.position(3)+0.01,'bp','MarkerFaceColor','b','MarkerSize',10);
        view(0,90); axis equal; axis([-1.2*R0 1.2*R0 -1.2*R0 1.2*R0]); box on
        title(ttls{k},'Interpreter','latex');
    end
    savefig(fullfile(script_folder,'thorax_theta_sensor_configs.fig'));
    saveas(gcf,fullfile(script_folder,'thorax_theta_sensor_configs.png'));
catch err
    fprintf('WARNING: sensor-config figure failed: %s\n',err.message);
end

% --- Polar plot of optimal azimuths vs ROI azimuth ---
try
    figure('Name','Thorax theta: azimuths');
    polarplot(mod(theta_even,2*pi),1.00*ones(numel(theta_even),1),'ko','DisplayName','even'); hold on
    polarplot(mod(theta_opt{1},2*pi),1.05*ones(numel(theta_opt{1}),1),'gs','DisplayName','a-opt');
    polarplot(mod(theta_opt{2},2*pi),1.10*ones(numel(theta_opt{2}),1),'rd','DisplayName','d-opt');
    polarplot([roi_azimuth roi_azimuth],[0 1.2],'b-','LineWidth',1.5,'DisplayName','ROI azimuth');
    legend show; title('Sensor azimuths (ROI at 45 deg)');
    savefig(fullfile(script_folder,'thorax_theta_azimuths.fig'));
    saveas(gcf,fullfile(script_folder,'thorax_theta_azimuths.png'));
catch err
    fprintf('WARNING: azimuth figure failed: %s\n',err.message);
end

% --- Mid-height posterior-std slice: prior vs even vs a-opt ---
try
    C = fmdl_hom.elem_centroids;
    slice = abs(C(:,3)-zs) < 0.6*height;   % thin extruded model: take all
    xc = C(slice,1); yc = C(slice,2);
    cmax = sqrt(max(post_var_prior(slice)));
    figure('Name','Thorax theta: posterior std','Position',[100 100 1200 380]);
    vals = {sqrt(post_var_prior(slice)),sqrt(post_var_even(slice)),sqrt(post_var_aopt(slice))};
    names = {'prior','posterior (even)','posterior (a-opt)'};
    for ip = 1:3
        subplot(1,3,ip);
        scatter(xc,yc,18,vals{ip},'filled');
        hold on; plot(anomaly.position(1),anomaly.position(2),'rp','MarkerFaceColor','r');
        axis equal tight; colorbar; clim([0 cmax]);
        title(names{ip}); xlabel('x'); ylabel('y');
    end
    savefig(fullfile(script_folder,'thorax_theta_posterior.fig'));
    saveas(gcf,fullfile(script_folder,'thorax_theta_posterior.png'));
catch err
    fprintf('WARNING: posterior figure failed: %s\n',err.message);
end
end

%% ---- 35-point quadrature on the reference tetrahedron ----
function [coord,weights] = quad_rule_35()
coord = [    0 0 0
    0 0 1
    0 1 0
    1 0 0
    0 0 0.75
    0 0 0.25
    0 0.75 0
    0 0.75 0.25
    0 0.25 0
    0 0.25 0.75
    0.75 0 0
    0.75 0 0.25
    0.75 0.25 0
    0.25 0 0
    0.25 0 0.75
    0.25 0.75 0
    0.5 0.5 0
    0.5 0 0.5
    0.5 0 0
    0 0.5 0.5
    0 0.5 0
    0 0 0.5
    0.25 0.25 0
    0.25 0.25 0.5
    0.25 0 0.25
    0.25 0 0.5
    0.25 0.5 0.25
    0.25 0.5 0
    0 0.25 0.25
    0 0.25 0.5
    0 0.5 0.25
    0.5 0.25 0.25
    0.5 0.25 0
    0.5 0 0.25
    0.25 0.25 0.25];
weights = [ -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.3047619047619048];
end
