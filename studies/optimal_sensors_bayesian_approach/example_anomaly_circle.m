%% EXAMPLE: Bayesian optimization of sensor positions on a circle (MDEIT)
%
% Sensors are constrained to a circle of radius rs at half the height of a
% cylindrical tank. A conductive spherical anomaly sits off-center at half
% height. The prior encodes "we know roughly WHERE the anomaly may be, but
% not its value": large prior variance in a ball around the suspected
% anomaly location, small variance in the background (cf. Hyvonen,
% Seppanen & Staboulis, Case 1). Because the informative region is
% off-center, the optimal sensors are NOT evenly spaced: they migrate
% toward the anomaly side, and the A-/D-optimality costs drop considerably
% compared to the evenly spaced configuration.
%
% Design criteria (3-axis MDEIT, linearized around the prior mean):
%   A-optimality: minimize trace(H^-1),      H = J'*inv(Gn)*J + inv(Gpr)
%   D-optimality: minimize -log det(H)  (=  log det posterior covariance)
%
% The gradients are analytic: dphi = -2*<W, dJ/dp>_F with
%   W_A = inv(Gn)*J*H^-2,   W_D = inv(Gn)*J*H^-1,
% and dJ/dp assembled from (i) extra adjoint solves for dlambda/dp and
% (ii) the explicit derivative of the Biot-Savart kernel integrals.
%
% Differences w.r.t. main.m (intentional):
%   * compute_r_matrices_local is vectorized over elements (much faster);
%   * the triangular solves for H^-1/H^-2 use the correct order
%     (H^-1 = L'\(L\I) for H = L*L');
%   * compute_dJxyz_xyz indexes dlambda{dim,p} (main.m uses dlambda{p});
%   * cost and gradient share one Jacobian evaluation per iterate;
%   * the cylindrical chain rule is evaluated at the CURRENT iterate.
%
% Run "quick_test = true; example_anomaly_circle" for a fast smoke test
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

model_parameters = create_kai_3d_model_parameters(l0, z0, sigma0, I0);

model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/8;
model_parameters.numOfElectrodesPerRing = 16;
model_parameters.numOfRings = 1;      % single electrode ring at half height

% Conductive spherical anomaly, off-center, at half height
model_parameters.material.name   = 'sphere';
model_parameters.material.type   = 'spherical';
model_parameters.material.radius = model_parameters.radius/3;
model_parameters.material.position = ...
    [model_parameters.radius/2, 0, model_parameters.height/2];

% Sensors: M evenly spaced on a circle of radius rs at half height
n_sensors = 8;
rs = 1.5*model_parameters.radius;         % circle radius (fixed)
zs = model_parameters.height/2;           % circle height (fixed)

theta_even = 2*pi*(0:n_sensors-1)'/n_sensors;

model_parameters.numOfSensors = n_sensors;
model_parameters.sensorRadius = rs;
model_parameters.sensorPositions = theta_to_locations(theta_even,rs,zs);

background_conductivity = 3.28e-1/sigma0;
anomaly_conductivity = 5*background_conductivity;  % conductive anomaly

%% Simulation parameters
opt_modes = {'a-opt','d-opt'};

max_iterations = 40;
n_starts = 1;                  % 1 = start from the even configuration only
                                % (multistart with n_starts=3 was checked: all
                                % starts converge to the same a-opt minimum;
                                % d-opt starts agree to within 0.03%)
do_gradient_check = false;

% The noise level is calibrated from the whitened Jacobian spectrum so
% that d_target eigenmodes are data-dominated. This guarantees the design
% problem sits in the regime where the data (and hence the sensor
% positions) actually matter: with too much noise the prior dominates and
% the cost is insensitive to the positions; with too little noise every
% configuration recovers the ROI almost perfectly and the criteria
% saturate.
% Note the interplay between d_target and the ROI size: the data can
% constrain at most ~d_target directions of the prior, so if the ROI
% contains many more elements than d_target, NO sensor arrangement can
% reduce the A-cost substantially and the design becomes insensitive.
d_target = 100;                % desired # of data-dominated modes

prior_std_background = 0.005*background_conductivity;
prior_std_roi        = 1.00*background_conductivity;
roi_radius_factor    = 1.0;    % ROI ball radius = factor * anomaly radius

alpha_rep = 0;                 % pairwise sensor repulsion weight (optional)

if quick_test
    max_iterations = 3;
    n_starts = 1;
    do_gradient_check = true;
end

%% Make forward model
[~,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl = fmdls{1};

% Homogeneous image = prior mean (the OED linearization point)
imgh = mk_image_mdeit(fmdl,background_conductivity);

% Image with the anomaly (ground truth, used for visualization only)
imgi = add_material_properties(imgh,[background_conductivity,anomaly_conductivity]);

n_stim  = numel(fmdl.stimulation);
n_elem  = size(fmdl.elems,1);
n_nodes = size(fmdl.nodes,1);

fprintf('Model: %i elements, %i nodes, %i stims, %i sensors (3-axis)\n',...
    n_elem,n_nodes,n_stim,n_sensors);

A = @(x) M(imgh,x);

%% Prior covariance: informative about WHERE the anomaly may be
anomaly_center = model_parameters.material.position;
roi_radius = roi_radius_factor*model_parameters.material.radius;

centroid_dist = sqrt(sum((fmdl.elem_centroids - anomaly_center).^2,2));
roi = centroid_dist <= roi_radius;

prior_variance = prior_std_background^2*ones(n_elem,1);
prior_variance(roi) = prior_std_roi^2;

Gamma_prior     = spdiags(prior_variance,0,n_elem,n_elem);
inv_Gamma_prior = spdiags(1./prior_variance,0,n_elem,n_elem);

fprintf('ROI: %i of %i elements (%.1f%% of prior trace)\n',nnz(roi),n_elem,...
    100*sum(prior_variance(roi))/sum(prior_variance));

%% Noise covariance (white), calibrated from the whitened Jacobian spectrum
% The singular values s_i of J0*Gpr^(1/2) measure how strongly the data
% constrain the prior eigenmodes. Setting the noise std to s(d_target)
% makes exactly ~d_target modes data-dominated (s_i^2/noise_var > 1).
img0 = assign_sensor_locations(imgh,theta_to_locations(theta_even,rs,zs));
B0 = fwd_solve_mdeit(img0);
max_B = max([abs(B0.Bx(:));abs(B0.By(:));abs(B0.Bz(:))]);

J0 = calc_jacobian_3axis_local(img0,A);
[~,S,Vs] = svd(J0 .* sqrt(prior_variance)','econ');   % unwhitened spectrum
s = diag(S);

noise_std = s(d_target);
variance_noise = noise_std^2;

n_data = 3*n_stim*n_sensors;
inv_Gamma_noise = (1/variance_noise)*speye(n_data,n_data);

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
    costgrad{iom} = @(theta) costgrad_theta_3axis(imgh,theta,rs,zs,...
        inv_Gamma_prior,inv_Gamma_noise,A,opt_modes{iom},alpha_rep);
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
    'OptimalityTolerance',1e-8,'StepTolerance',1e-12);

theta_opt = cell(1,numel(opt_modes));
phi_opt = inf(1,numel(opt_modes));

for iom = 1:numel(opt_modes)
    fprintf('\n=== Optimizing %s (3-axis MDEIT) ===\n',opt_modes{iom});

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

%% Posterior variance diagnostics (prior vs even vs A-optimized)
post_var_even = compute_posterior_variance_diag(imgh,theta_even,rs,zs,...
    inv_Gamma_prior,inv_Gamma_noise,A);
post_var_aopt = compute_posterior_variance_diag(imgh,theta_opt{1},rs,zs,...
    inv_Gamma_prior,inv_Gamma_noise,A);

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
save(fullfile(script_folder,'data','example_anomaly_circle_results.mat'),'results');

%% Figures
figures_folder = fullfile(script_folder,'figures');
if ~exist(figures_folder,'dir'), mkdir(figures_folder); end

% --- Sensor positions on the FEM model ---
try
    img_even = assign_sensor_locations(imgi,theta_to_locations(theta_even,rs,zs));
    img_aopt = assign_sensor_locations(imgi,theta_to_locations(theta_opt{1},rs,zs));
    img_dopt = assign_sensor_locations(imgi,theta_to_locations(theta_opt{2},rs,zs));

    figure('Name','Sensor positions');
    hold on
    show_fem(imgi);
    plot_sensors(img_even,false,'k','o');
    plot_sensors(img_aopt,false,'g','s');
    h = plot_sensors(img_dopt,false,'r','d');
    tt = linspace(0,2*pi,200);
    plot3(h,rs*cos(tt),rs*sin(tt),zs*ones(size(tt)),'b--');
    axis([-1.1*rs 1.1*rs -1.1*rs 1.1*rs 0 model_parameters.height])
    box on; grid on; view(3);
    title('black o: even | green s: a-opt | red d: d-opt');
    savefig(fullfile(figures_folder,'example_anomaly_circle_sensors.fig'));
    saveas(gcf,fullfile(figures_folder,'example_anomaly_circle_sensors.png'));
catch err
    fprintf('WARNING: sensor-position figure failed: %s\n',err.message);
end

% --- Sensor angles (polar view; anomaly sits at azimuth 0) ---
try
    figure('Name','Sensor angles');
    polarplot(theta_even,ones(n_sensors,1),'ko','DisplayName','even'); hold on
    polarplot(mod(theta_opt{1},2*pi),1.05*ones(n_sensors,1),'gs','DisplayName','a-opt');
    polarplot(mod(theta_opt{2},2*pi),1.10*ones(n_sensors,1),'rd','DisplayName','d-opt');
    polarplot([0 0],[0 1.2],'b-','DisplayName','anomaly azimuth');
    legend show
    title('Sensor azimuths (anomaly at 0)');
    savefig(fullfile(figures_folder,'example_anomaly_circle_angles.fig'));
    saveas(gcf,fullfile(figures_folder,'example_anomaly_circle_angles.png'));
catch err
    fprintf('WARNING: angle figure failed: %s\n',err.message);
end

% --- Posterior std at the mid-height slice ---
try
    slice = abs(fmdl.elem_centroids(:,3)-zs) < model_parameters.maxsz;
    xc = fmdl.elem_centroids(slice,1);
    yc = fmdl.elem_centroids(slice,2);
    cmax = sqrt(max(prior_variance(slice)));

    figure('Name','Posterior std, mid slice','Position',[100 100 1200 400]);
    vals = {sqrt(prior_variance(slice)),sqrt(post_var_even(slice)),sqrt(post_var_aopt(slice))};
    names = {'prior','posterior (even)','posterior (a-opt)'};
    for ip = 1:3
        subplot(1,3,ip);
        scatter(xc,yc,25,vals{ip},'filled');
        axis equal tight; colorbar; clim([0 cmax]);
        title(names{ip}); xlabel('x'); ylabel('y');
    end
    savefig(fullfile(figures_folder,'example_anomaly_circle_posterior.fig'));
    saveas(gcf,fullfile(figures_folder,'example_anomaly_circle_posterior.png'));
catch err
    fprintf('WARNING: posterior figure failed: %s\n',err.message);
end

fprintf('\nDone.\n');

%% ======================== LOCAL FUNCTIONS ========================

%% Map circle angles to 3D sensor locations
function sensor_locations = theta_to_locations(theta,rs,zs)
theta = theta(:);
sensor_locations = [rs*cos(theta), rs*sin(theta), zs*ones(numel(theta),1)];
end

%% Cost and gradient in the angle parametrization (3-axis MDEIT)
function [phi,dphi] = costgrad_theta_3axis(img,theta,rs,zs,...
    inv_Gamma_prior,inv_Gamma_noise,A,opt_mode,alpha_rep)

assert(ismember(opt_mode,{'a-opt','d-opt'}));

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);
assert(numel(theta) == n_sensors);

sensor_locations = theta_to_locations(theta,rs,zs);
img = assign_sensor_locations(img,sensor_locations);

% Jacobian (shared by cost and gradient)
J = calc_jacobian_3axis_local(img,A);

H = full(J.'*inv_Gamma_noise*J + inv_Gamma_prior);
L = chol(H,'lower');   % H = L*L'

switch opt_mode
    case 'a-opt'
        X = L\eye(n_elem);            % X = inv(L)
        phi = sum(X(:).^2);           % trace(H^-1) = ||inv(L)||_F^2
    case 'd-opt'
        phi = -2*sum(log(diag(L)));   % -log det(H)
end

if alpha_rep > 0
    phi = phi + alpha_rep*repulsion_cost(sensor_locations);
end

if nargout < 2
    return
end

% ---- Gradient ----
% H^-1 = inv(L)'*inv(L);  computed with the correct solve order
switch opt_mode
    case 'a-opt'
        invH  = L'\X;                 % H^-1
        Z     = L\invH;
        Hpow  = L'\Z;                 % H^-2
    case 'd-opt'
        X    = L\eye(n_elem);
        Hpow = L'\X;                  % H^-1
end

W = inv_Gamma_noise*(J*Hpow);         % [3*n_stim*n_sensors x n_elem]

% Derivative of the stacked Jacobian w.r.t. sensor coordinates
dJ = compute_dJxyz_xyz(img);          % dJ{dim,p}: [n_sensors x n_stim x n_elem]

dphidp = zeros(3,n_sensors);          % rows: p = x,y,z
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

% Chain rule to the angle parametrization (evaluated at the CURRENT theta)
dphi = -rs*sin(theta(:)).*dphidp(1,:)' + rs*cos(theta(:)).*dphidp(2,:)';

if alpha_rep > 0
    [dGx,dGy,~] = repulsion_grad_cartesian(sensor_locations);
    dphi = dphi + alpha_rep*(-rs*sin(theta(:)).*dGx + rs*cos(theta(:)).*dGy);
end

end

%% Diagonal of the posterior covariance for a given configuration
function post_var = compute_posterior_variance_diag(img,theta,rs,zs,...
    inv_Gamma_prior,inv_Gamma_noise,A)

n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,theta_to_locations(theta,rs,zs));
J = calc_jacobian_3axis_local(img,A);

H = full(J.'*inv_Gamma_noise*J + inv_Gamma_prior);
L = chol(H,'lower');
X = L\eye(n_elem);                    % inv(L)
post_var = sum(X.^2,1)';              % diag(H^-1) = diag(inv(L)'*inv(L))
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

%% Stacked 3-axis Jacobian
function J = calc_jacobian_3axis_local(img,A)
Jx = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,1);
Jy = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,2);
Jz = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,3);
J = [Jx;Jy;Jz];
end

%% Schur complement of the CEM system matrix (SPD operator)
function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
end

%% Gamma (measurement) matrices for the current sensor positions
function img = compute_gamma_matrices_local(img)

mu_factor = img.fwd_model.mu0/(4*pi);

num_sensors = numel(img.fwd_model.sensors);

sensor_locations = zeros(num_sensors,3);
for i = 1:num_sensors
    sensor_locations(i,:) = img.fwd_model.sensors(i).position;
end

[Rx,Ry,Rz,fmdl] = compute_r_matrices_local(img.fwd_model,sensor_locations);
img.fwd_model = fmdl; % CRUCIAL: keep R in sync with the sensor positions

R.Rx = Rx; R.Ry = Ry; R.Rz = Rz;
G = img.fwd_model.G;

Sigma = spdiags(img.elem_data(:),0,length(img.elem_data),length(img.elem_data));

% g_{dl}^m: components of measurement axis d of sensor m in the R^3 basis
g = zeros(num_sensors,3,3);
for m = 1:num_sensors
    g(m,:,:) = [...
        img.fwd_model.sensors(m).axes.axis1;
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

%% Biot-Savart element integrals (vectorized over elements and quad points)
function [Rx,Ry,Rz,fmdl] = compute_r_matrices_local(fmdl,sensor_locations)

[coord,weights] = quad_rule_35();

numElements = size(fmdl.elems,1);
numSensors  = size(sensor_locations,1);
numQuadPts  = length(weights);

V = reshape(fmdl.nodes(fmdl.elems',:),[4,numElements,3]);
v1 = squeeze(V(1,:,:)); v2 = squeeze(V(2,:,:));
v3 = squeeze(V(3,:,:)); v4 = squeeze(V(4,:,:));

a = v2-v1; b = v3-v1; c = v4-v1;
detJ = abs(sum(a .* cross(b,c,2),2))';   % 1 x numElements

% Quadrature points in all elements: numQuadPts x 3 x numElements
Jm = zeros(3,3,numElements);
Jm(:,1,:) = permute(a,[2 3 1]);
Jm(:,2,:) = permute(b,[2 3 1]);
Jm(:,3,:) = permute(c,[2 3 1]);

xi = reshape(v1',[3,1,numElements]) + ...
    Jm(:,1,:).*reshape(coord(:,1),[1,numQuadPts,1]) + ...
    Jm(:,2,:).*reshape(coord(:,2),[1,numQuadPts,1]) + ...
    Jm(:,3,:).*reshape(coord(:,3),[1,numQuadPts,1]);
xi = permute(xi,[2 1 3]);                % numQuadPts x 3 x numElements

Rx = zeros(numSensors,numElements);
Ry = zeros(numSensors,numElements);
Rz = zeros(numSensors,numElements);

w = reshape(weights,[numQuadPts,1,1]);

for m = 1:numSensors
    rm = sensor_locations(m,:);
    dm = rm - xi;                        % numQuadPts x 3 x numElements
    nrm3 = sum(dm.^2,2).^(3/2);
    f = dm ./ nrm3;
    integ = sum(f.*w,1);                 % 1 x 3 x numElements
    Rx(m,:) = squeeze(integ(1,1,:))'.*(detJ/6);
    Ry(m,:) = squeeze(integ(1,2,:))'.*(detJ/6);
    Rz(m,:) = squeeze(integ(1,3,:))'.*(detJ/6);
end

fmdl.R.Rx = Rx;
fmdl.R.Ry = Ry;
fmdl.R.Rz = Rz;

end

%% Derivatives dR_j/dp of the Biot-Savart integrals, all p = x,y,z
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
    dm_vec = rm - xi;                       % numQuadPts x 3 x numElements
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

%% Derivatives of the adjoint solutions w.r.t. all sensor coordinates
% dlambda{dim,p} = d(lambda for measurement axis dim)/d(coordinate p)
function [dlambda,dR1t,dR2t] = compute_dlambdaxyz_xyz(img)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);

dlambda = cell(3,3);
dR1t = cell(3,3);
dR2t = cell(3,3);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

sigma = img.elem_data;
A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));
Mfun = @(x) x ./ d;
Nfun = @(x) x ./ d;

for dim = 1:3
    switch dim
        case 1
            dR1 = dRmkj_xyz(img.fwd_model,3);
            dR2 = dRmkj_xyz(img.fwd_model,2);
            G1 = G.Gy; G2 = G.Gz;
        case 2
            dR1 = dRmkj_xyz(img.fwd_model,1);
            dR2 = dRmkj_xyz(img.fwd_model,3);
            G1 = G.Gz; G2 = G.Gx;
        case 3
            dR1 = dRmkj_xyz(img.fwd_model,2);
            dR2 = dRmkj_xyz(img.fwd_model,1);
            G1 = G.Gx; G2 = G.Gy;
    end

    rhs1 = mu_factor*( (dR1{1} .* sigma.')*G1 - (dR2{1} .* sigma.')*G2 ); % p = x
    rhs2 = mu_factor*( (dR1{2} .* sigma.')*G1 - (dR2{2} .* sigma.')*G2 ); % p = y
    rhs3 = mu_factor*( (dR1{3} .* sigma.')*G1 - (dR2{3} .* sigma.')*G2 ); % p = z

    dlambdaX = zeros(n_nodes,num_sensors);
    dlambdaY = zeros(n_nodes,num_sensors);
    dlambdaZ = zeros(n_nodes,num_sensors);

    parfor m = 1:num_sensors
        [dlambdaX(:,m),~,~] = pcg(A_matrix,rhs1(m,:)',1e-9,numel(sigma),Mfun,Nfun);
        [dlambdaY(:,m),~,~] = pcg(A_matrix,rhs2(m,:)',1e-9,numel(sigma),Mfun,Nfun);
        [dlambdaZ(:,m),~,~] = pcg(A_matrix,rhs3(m,:)',1e-9,numel(sigma),Mfun,Nfun);
    end

    dlambda{dim,1} = dlambdaX;
    dlambda{dim,2} = dlambdaY;
    dlambda{dim,3} = dlambdaZ;

    dR1t{dim,1} = dR1{1}; dR1t{dim,2} = dR1{2}; dR1t{dim,3} = dR1{3};
    dR2t{dim,1} = dR2{1}; dR2t{dim,2} = dR2{2}; dR2t{dim,3} = dR2{3};
end

end

%% Derivative of the stacked 3-axis Jacobian w.r.t. sensor coordinates
function dJ = compute_dJxyz_xyz(img)

num_sensors = numel(img.fwd_model.sensors);
n_elem      = size(img.fwd_model.elems,1);
num_stim    = numel(img.fwd_model.stimulation);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

elemV = img.fwd_model.elem_volume(:);

u = fwd_solve(img);
u = u.volt;

GxU = G.Gx*u;
GyU = G.Gy*u;
GzU = G.Gz*u;

[dlambda,dR1dp,dR2dp] = compute_dlambdaxyz_xyz(img);

G1U = cell(1,3);
G2U = cell(1,3);
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
            % --- dJ1: contribution from dlambda/dp ---
            % NOTE: dlambda must be indexed {dim,p} (main.m uses {p})
            dlGx = G.Gx*dlambda{dim,p}(:,m);
            dlGy = G.Gy*dlambda{dim,p}(:,m);
            dlGz = G.Gz*dlambda{dim,p}(:,m);

            tmp = dlGx.*GxU + dlGy.*GyU + dlGz.*GzU;
            dJ1(m,:,:) = tmp.' .* elemV(:).';

            % --- dJ2: contribution from dR/dp ---
            dJ2(m,:,:) = -mu_factor * ( ...
                dR2dp{dim,p}(m,:) .* G2U{dim}.' - dR1dp{dim,p}(m,:) .* G1U{dim}.' ...
                );
        end
        dJ{dim,p} = dJ1 - dJ2;
    end
end

end

%% Jacobian dB/dsigma for one measurement axis (adjoint method)
function J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,select_sensor_axis)

mu0 = img.fwd_model.mu0;

n_nodes = size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

img = compute_gamma_matrices_local(img);

R = img.fwd_model.R;
G = img.fwd_model.G;

switch select_sensor_axis
    case 1
        Gamma = img.Gamma1;
        R1 = R.Rz.'; R2 = R.Ry.';
    case 2
        Gamma = img.Gamma2;
        R1 = R.Rx.'; R2 = R.Rz.';
    case 3
        Gamma = img.Gamma3;
        R1 = R.Ry.'; R2 = R.Rx.';
    otherwise
        error('Invalid sensor axis');
end

u = fwd_solve(img);
u = u.volt;

lambda = zeros(n_nodes,num_sensors);

A_matrix = A(img.elem_data);

d = sqrt(diag(A_matrix));
Mfun = @(x) x ./ d;
Nfun = @(x) x ./ d;

GammaT = Gamma.';

num_elements = numel(img.elem_data);
parfor m = 1:num_sensors
    [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-8,num_elements,Mfun,Nfun);
end

Gx_times_lambda = G.Gx*lambda;
Gy_times_lambda = G.Gy*lambda;
Gz_times_lambda = G.Gz*lambda;

Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);
elemV = reshape(elemV,[n_elem 1 1]);

GxL = reshape(Gx_times_lambda,[n_elem 1 num_sensors]);
GyL = reshape(Gy_times_lambda,[n_elem 1 num_sensors]);
GzL = reshape(Gz_times_lambda,[n_elem 1 num_sensors]);

GxU = reshape(Gx_times_u,[n_elem num_stim 1]);
GyU = reshape(Gy_times_u,[n_elem num_stim 1]);
GzU = reshape(Gz_times_u,[n_elem num_stim 1]);

dfdx = elemV .* ( GxL.*GxU + GyL.*GyU + GzL.*GzU );

Rx_ = reshape(R.Rx.',[n_elem 1 num_sensors]);
Ry_ = reshape(R.Ry.',[n_elem 1 num_sensors]);
Rz_ = reshape(R.Rz.',[n_elem 1 num_sensors]);

dCxdp = ( -Rz_.*GyU + Ry_.*GzU );
dCydp = ( -Rx_.*GzU + Rz_.*GxU );
dCzdp = ( -Ry_.*GxU + Rx_.*GyU );

g = zeros(num_sensors,3,3);
for m = 1:num_sensors
    g(m,:,:) = [...
        img.fwd_model.sensors(m).axes.axis1;
        img.fwd_model.sensors(m).axes.axis2;
        img.fwd_model.sensors(m).axes.axis3];
end

gx = reshape(g(:,select_sensor_axis,1),[1 1 num_sensors]);
gy = reshape(g(:,select_sensor_axis,2),[1 1 num_sensors]);
gz = reshape(g(:,select_sensor_axis,3),[1 1 num_sensors]);

dfdp = mu_factor*( gx.*dCxdp + gy.*dCydp + gz.*dCzdp );

dfd = dfdx + dfdp;             % [n_elem x num_stim x num_sensors]

dfd = permute(dfd,[3 2 1]);    % [num_sensors x num_stim x n_elem]
J = reshape(dfd,num_sensors*num_stim,n_elem);

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
end

%% 35-point quadrature rule on the reference tetrahedron
function [coord,weights] = quad_rule_35()
coord = [    0.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  1.0000000000000000
    0.0000000000000000  1.0000000000000000  0.0000000000000000
    1.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.7500000000000000
    0.0000000000000000  0.0000000000000000  0.2500000000000000
    0.0000000000000000  0.7500000000000000  0.0000000000000000
    0.0000000000000000  0.7500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.7500000000000000
    0.7500000000000000  0.0000000000000000  0.0000000000000000
    0.7500000000000000  0.0000000000000000  0.2500000000000000
    0.7500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.7500000000000000
    0.2500000000000000  0.7500000000000000  0.0000000000000000
    0.5000000000000000  0.5000000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.5000000000000000
    0.5000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.5000000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.2500000000000000  0.5000000000000000
    0.2500000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.5000000000000000  0.2500000000000000
    0.2500000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.2500000000000000  0.2500000000000000];

weights = [  -0.0119047619047619
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
