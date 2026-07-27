clc;clear all;close all

% This script optimizes sensor configuration in a ring at half height. The
% radius and the angle of the sensors is free!

%% Prepare workspace
% Get the full path of the current script
fullpath = mfilename('fullpath');
% Extract just the folder
script_folder = fileparts(fullpath);
cd(script_folder);

% Have to add the functions path manually so prepare_workspace runs
grandgrandparent_folder =fileparts(fileparts(fileparts(script_folder)));
addpath(genpath(fullfile(grandgrandparent_folder,'functions')));

model_folder = prepare_workspace(script_folder);

% Add the local model folder to the path
addpath(genpath(fullfile(script_folder,'model')));

rng(1);

colors = [228,26,28;
55,126,184;
77,175,74;
152,78,163;
255,127,0;
255,255,51;
166,86,40;
247,129,191]/255;

markers = {'+','o','d','s','p','x','^' , '>' , '<' };


% Make sure your parallel pool is started
pool = gcp('nocreate');  % get current pool, do not create
if isempty(pool)
    pool = parpool('local');  % only create if none exists
end
% Run EIDORS startup on all workers
pctRunOnAll(sprintf('setupEidors(''%s'');', script_folder));

%% Model parameters 
z_0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z_0*I0/(l0^2); %(V)
sigma0 = l0/z_0; %(S/m)
J0 = I0/(l0^2);

R0 = 1.5;

model_parameters = create_default_3d_model_parameters(l0, z_0, sigma0, I0);

model_parameters.maxsz = min(model_parameters.height,model_parameters.radius)/3;
model_parameters.height = 2;
model_parameters.numOfElectrodesPerRing = 8;
model_parameters.numOfRings = 3;
model_parameters.numOfSensors = model_parameters.numOfElectrodesPerRing*model_parameters.numOfRings;
model_parameters.sensorRadius = R0;
model_parameters.material.type = 'spherical';

model_parameters.material.position(1) = 0; model_parameters.material.position(2) = 0;
model_parameters.material.radius = model_parameters.radius/5;
% model_parameters.material.position(1) = 0.95*(model_parameters.radius-model_parameters.material.radius);
model_parameters.material.position(3) = 0.5*model_parameters.height;

background_conductivity = 3.28e-1/sigma0;
anomaly_conductivity = background_conductivity*1.1;

inj = [0 3]; %skip 0 pattern (pg 172)
meas = [0 3]; 
current_amplitude = 2.4e-3/I0;

mu0 = 1.0;

rmin = 1.05*sqrt(model_parameters.height^2+model_parameters.radius^2); % rhom is measured from (0,0,0), so the minimaml r cant just be the height
rmax = 2*rmin;
%% Make forward models

% Make forward homogeneous/inhomogeneous forward models with material
[~,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl = fmdls{1};

% Make model for reconstruction
model_parameters_recon = model_parameters;
model_parameters_recon.maxsz = min(model_parameters.height,model_parameters.radius)/10;

model_parameters_recon = rmfield(model_parameters_recon,'material');

[~,fmdls] = mk_mdeit_model(model_parameters_recon,model_folder);
fmdl_recon = fmdls{1};

%% Assign stimulation pattern
stimulation = mk_stim_patterns(model_parameters.numOfElectrodesPerRing,model_parameters.numOfRings,inj,meas,{'meas_current'},current_amplitude);

fmdl.stimulation = stimulation;
fmdl.stimulation = stimulation;
fmdl_recon.stimulation = stimulation;

%% Make images

imgh = mk_image_mdeit(fmdl,1.0);
imgh = add_material_properties(imgh, [background_conductivity,background_conductivity]);

imgi = mk_image_mdeit(fmdl,1.0);
imgi = add_material_properties(imgi,  [background_conductivity,anomaly_conductivity]);

img_recon = mk_image_mdeit(fmdl_recon,1.0);

imgh.calc_colours.ref_level = background_conductivity;
imgh.calc_colours.greylev = -0.01;
imgh.calc_colours.transparency_thresh = 0.01;

imgi.calc_colours.ref_level = background_conductivity;
imgi.calc_colours.greylev = -0.01;
imgi.calc_colours.transparency_thresh = 0.01;

% For visualization
figure
subplot(1,3,1)
show_fem_transparent_edges(imgh);
plot_sensors(imgh);
view(0,90)
subplot(1,3,2)
show_fem_transparent_edges(imgi);
plot_sensors(imgi);
view(0,90)
subplot(1,3,3)
show_fem_transparent_edges(img_recon);
plot_sensors(img_recon);
view(0,90)

%% Define prior and noise covariance matrices

num_sensors = model_parameters.numOfSensors;

Bi = fwd_solve_mdeit(imgi);
Bh = fwd_solve_mdeit(imgh);
dB = [Bi.Bx(:)-Bh.Bx(:);Bi.By(:)-Bh.By(:);Bi.Bz(:)-Bh.Bz(:)];
max_B = max(dB(:));

% Set the noise variance with respect to the data magnitude
noise_std_deviation = max_B*1e-1; %when this is times 1e-5 I get a singular matrix later, problems!
variance_noise = noise_std_deviation^2;

n_stim = numel(imgh.fwd_model.stimulation);

Gamma_noise_3_axis = variance_noise.*speye(3*n_stim*num_sensors,3*n_stim*num_sensors);
inv_Gamma_noise_3_axis = inv(Gamma_noise_3_axis);

inv_Gamma_noise = inv_Gamma_noise_3_axis;

%% Determine the prior variance
% such that lambda = sigma_epsilon^2/sigma_p^2 is the optimal tykhonov regularization 
% paramter according to L-Curve method for the initial sensor configuration

imgtemp = img_recon;

sensor_positions_0 = fetch_sensor_locations(imgh);

theta0 = atan2(sensor_positions_0(:,2),sensor_positions_0(:,1));
z0 = sensor_positions_0(:,3);

imgtemp = assign_sensor_locations(imgtemp,sensor_positions_0);

A = @(x) M(img_recon,x);

[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtemp,A);
J = [Jx;Jy;Jz];

n_elem = size(J,2);

% Add noise to measurements
dy_opt = dB + sqrt(variance_noise)*randn(size(dB,1),1);

% lambda_vector = logspace(-15,2,50);
% [lambda_opt,dx] = l_curve_method_new(J,dy_opt,lambda_vector);
% 
% variance_prior = variance_noise/lambda_opt;

%% Build smooth prior encoding the anomaly information

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

    % Sparse-efficient version (recommended)

    [i,j,val] = find(kappa);

    xi = centroids(i,:);
    xj = centroids(j,:);

    D2 = sum((xi - xj).^2, 2);

    weights = val .* exp(-D2/(2*lambda^2));

    Gamma = sparse(i,j,weights,n_elem,n_elem);

    %Symmetrize (since Gamma is supposed to be symmetric positive definite)
    Gamma = 0.5*(Gamma+Gamma.');

    % Regularize to guarante positive defintness
    Gamma = Gamma + 1e-6*speye(n_elem);
    end

fprintf('Building smooth prior\n');

anomaly_elem_ids = imgh.fwd_model.mat_idx{2}(:);
anomaly_idx = false(size(imgh.fwd_model.elems,1),1);
anomaly_idx(anomaly_elem_ids) = true;

in_block = anomaly_idx*anomaly_idx.';

bg_elem_ids = imgh.fwd_model.mat_idx{1}(:);
bg_idx = false(size(imgh.fwd_model.elems,1),1);
bg_idx(bg_elem_ids) = true;

out_block = (bg_idx)*(bg_idx).';

% Assemble kappa
kappa = 0.8 * sparse(in_block) + 0.03 * sparse(out_block);

Gamma_smooth = smooth_prior(imgh,0.5,kappa);

% Check if Gamma is block diagonal under permutation
% perm = [bg_elem_ids; lung_elem_ids];  % proposed block ordering
% Gamma_perm = Gamma_smooth(perm, perm);
% spy(Gamma_perm);  % should show two clean diagonal blocks if truly block diagonal

% tic
% fprintf('\t Inverting smooth prior (might take a while)\n');
% Gamma_prior = variance_prior*Gamma_smooth;
% Gamma_prior = Gamma_prior + 1e-6 * speye(size(Gamma_prior)); %need regularization because of ill-conditioning
% inv_Gamma_prior = Gamma_prior \ speye(size(Gamma_prior));
% inv_Gamma_prior_3_axis_2 = inv_Gamma_prior;
% disp(toc);

% This is faster inversion since it takes advantage of the block diagonal
% structure under an appropriate permutation
tic
fprintf('\t Inverting smooth prior (block diagonal structure)\n');

Gamma_prior = Gamma_smooth;

% Gamma_prior = variance_prior*Gamma_smooth;
% Gamma_prior = Gamma_prior + 1e-6 * speye(size(Gamma_prior));

% --- Permute to block diagonal form ---
perm = [bg_elem_ids; anomaly_elem_ids];
[~, iperm] = sort(perm);  % inverse permutation to restore original ordering later

Gamma_perm = Gamma_prior(perm, perm);

% --- Extract blocks ---
n_bg   = numel(bg_elem_ids);
n_lung = numel(anomaly_elem_ids);

idx_bg   = 1:n_bg;
idx_lung = (n_bg+1):(n_bg+n_lung);

G_bg   = Gamma_perm(idx_bg,   idx_bg);
G_lung = Gamma_perm(idx_lung, idx_lung);

% --- Invert each block independently ---
inv_G_bg   = G_bg   \ speye(n_bg);
inv_G_lung = G_lung \ speye(n_lung);

% --- Reassemble block diagonal inverse ---
inv_Gamma_perm = blkdiag(inv_G_bg, inv_G_lung);

% --- Restore original element ordering ---
inv_Gamma_prior = inv_Gamma_perm(iperm, iperm);
inv_Gamma_prior_3_axis = inv_Gamma_prior;
disp(toc);

%%
figure;
img_prior = imgh;
img_prior.elem_data = diag(Gamma_prior);

img_prior.calc_colours.ref_level = min(diag(Gamma_prior));
imgi.calc_colours.greylev = -0.01;
imgi.calc_colours.transparency_thresh = 0.01;

show_fem_transparent_edges(img_prior)
title('Smooth prior auto-correlation')
drawnow;


%% Define jacobian of coordinate transformations

function jacobian_coordinate_transformation = compute_jacobian_coordinate_transformation_spherical(sensor_locations,rmin,rmax)

n_sensors = size(sensor_locations,1);

jacobian_coordinate_transformation = zeros(3,3,n_sensors);

% Get spherical coordinates
rhom = sqrt(sensor_locations(:,1).^2+sensor_locations(:,2).^2+sensor_locations(:,3).^2);
rm = sqrt(sensor_locations(:,1).^2+sensor_locations(:,2).^2);
thetam = acos(sensor_locations(:,3)./rhom);
phim = sign(sensor_locations(:,2)).*acos(sensor_locations(:,1)./rm);

xim = log((rhom-rmin)./(rmax-rhom));

sigmoid = @(x) 1./(1+exp(-x));
dsigmoid = @(x) (sigmoid(x)).*(1-sigmoid(x));

drhomdxim = (rmax-rmin)*dsigmoid(xim);

jacobian_coordinate_transformation_xi = zeros(3,n_sensors);
jacobian_coordinate_transformation_xi(1,:) = sin(thetam).*cos(phim).*drhomdxim;
jacobian_coordinate_transformation_xi(2,:) = sin(thetam).*sin(phim).*drhomdxim;
jacobian_coordinate_transformation_xi(3,:) = cos(thetam).*drhomdxim;

jacobian_coordinate_transformation_theta = zeros(3,n_sensors);
jacobian_coordinate_transformation_theta(1,:) = rhom.*cos(thetam).*cos(phim);
jacobian_coordinate_transformation_theta(2,:) = rhom.*cos(thetam).*sin(phim);
jacobian_coordinate_transformation_theta(3,:) = -rm.*sin(thetam);

jacobian_coordinate_transformation_phi = zeros(3,n_sensors);
jacobian_coordinate_transformation_phi(1,:) = -rhom.*sin(thetam).*sin(phim);
jacobian_coordinate_transformation_phi(2,:) = rhom.*sin(thetam).*cos(phim);
jacobian_coordinate_transformation_phi(3,:) = zeros(1,n_sensors);

jacobian_coordinate_transformation(1,1,:) = jacobian_coordinate_transformation_xi(1,:);
jacobian_coordinate_transformation(1,2,:) = jacobian_coordinate_transformation_xi(2,:);
jacobian_coordinate_transformation(1,3,:) = jacobian_coordinate_transformation_xi(3,:);

jacobian_coordinate_transformation(2,1,:) = jacobian_coordinate_transformation_theta(1,:);
jacobian_coordinate_transformation(2,2,:) = jacobian_coordinate_transformation_theta(2,:);
jacobian_coordinate_transformation(2,3,:) = jacobian_coordinate_transformation_theta(3,:);

jacobian_coordinate_transformation(3,1,:) = jacobian_coordinate_transformation_phi(1,:);
jacobian_coordinate_transformation(3,2,:) = jacobian_coordinate_transformation_phi(2,:);
jacobian_coordinate_transformation(3,3,:) = jacobian_coordinate_transformation_phi(3,:);

end

jac_coord_transf_spherical = @(s)  compute_jacobian_coordinate_transformation_spherical(s,rmin,rmax);

%% Define mappings between cartesian and parametric coordinates

% Map sensor locations to vector in cartesian coordinates
function x = sensor_locations_to_vector_cartesian(sensor_locations)

n_sensors = size(sensor_locations,1);

x = zeros(3*n_sensors,1);

x(1:n_sensors) =  sensor_locations(:,1);
x(n_sensors+1:2*n_sensors) = sensor_locations(:,2);
x(2*n_sensors+1:3*n_sensors) = sensor_locations(:,3);

end

function sensor_locations = vector_to_sensor_locations_cartesian(x)

assert(mod(numel(x),3)==0);
n_sensors = numel(x)/3;

sensor_locations = zeros(n_sensors,3);

sensor_locations(:,1) = x(1:n_sensors);
sensor_locations(:,2) = x(n_sensors+1:2*n_sensors);
sensor_locations(:,3) = x(2*n_sensors+1:3*n_sensors);

end

% Map x-vector to generalized coordinates q (xi,theta) and back
function x = map_q_to_x_spherical(q,rmin,rmax)

assert(mod(numel(q),3)==0);
n_sensors = numel(q)/3;

sigmoid = @(x) 1./(1+exp(-x));

xim = q(1:n_sensors);
thetam = q(n_sensors+1:2*n_sensors);
phim = q(2*n_sensors+1:3*n_sensors);

rhom = rmin + (rmax-rmin).*sigmoid(xim);

x = [...
    rhom(:).*sin(thetam(:)).*cos(phim(:));...
    rhom(:).*sin(thetam(:)).*sin(phim(:));...
    rhom(:).*cos(thetam(:))];

end

function q = map_x_to_q_spherical(x,rmin,rmax)

assert(mod(numel(x),3)==0);
n_sensors = numel(x)/3;

xm = x(1:n_sensors);
ym = x(n_sensors+1:2*n_sensors);
zm = x(2*n_sensors+1:3*n_sensors);

% Get spherical coordinates
rhom = sqrt(xm.^2+ym.^2+zm.^2);
rm = sqrt(xm.^2+ym.^2);
thetam = acos(zm./rhom);
phim = sign(ym).*acos(xm./rm);

xim = log((rhom-rmin)./(rmax-rhom)); 

q = [xim(:);thetam(:);phim(:)];
end

x_to_q_spherical = @(x) map_x_to_q_spherical(x,rmin,rmax);
q_to_x_spherical = @(q) map_q_to_x_spherical(q,rmin,rmax);

%% Define vector_to_sensor_locations and vice-versa

x0 = sensor_locations_to_vector_cartesian(sensor_positions_0);

assert(norm(x0-[sensor_positions_0(:,1);sensor_positions_0(:,2);sensor_positions_0(:,3)])<1e-5)

q0 = x_to_q_spherical(x0);

% Full map from q coordinates to sensor locations and back for
% unconstrained optimization
vector_to_sensor_locations_spherical = @(q) vector_to_sensor_locations_cartesian(q_to_x_spherical(q));
sensor_locations_to_vector_spherical = @(sensor_locations) x_to_q_spherical(sensor_locations_to_vector_cartesian(sensor_locations));

%% Run optimization

max_iterations = 30;

options = optimoptions('fminunc',...
    'Display','iter','MaxIterations',max_iterations,...
    'StepTolerance',1e-7,'OptimalityTolerance',1e-7,...
    'Algorithm','quasi-newton','HessianApproximation','lbfgs',...
    'SpecifyObjectiveGradient',true,'UseParallel',true);


% Launch optimization from different initial conditions
% Use a grid based multi-start

num_of_sensors = model_parameters.numOfSensors;

if ~exist('data/sensor_positions_opt.mat','file')
    
    imgsh{1} = imgh;
    imgsh{1} = assign_sensor_locations(imgsh{1},sensor_positions_0);

    img_out = optimize_sensor_configuration(imgsh{1},inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,...
        jac_coord_transf_spherical,q_to_x_spherical,x_to_q_spherical,'a-opt','mdeit3',3,options);
    imgout{1} = img_out;

else
    imgsh{1} = imgh;

    imgtemp = imgh;
    var = load(fullfile(script_folder,'data/sensor_positions_opt.mat'));

    sensor_positions_opt = var.sensor_position_bfgs;
    imgtemp = assign_sensor_locations(imgtemp,sensor_positions_opt);
    imgout{1} = imgtemp;
end

%% Find sensor positions close to boundary

delta = 0.05;
sensor_positions_boundary = find_sensor_positions_boundary(fmdl,delta);

%%
figure

subplot(1,3,1)
show_fem_transparent_edges(imgsh{1});
plot_sensors(imgsh{1});
view(0,90)
axis([-rmax rmax -rmax rmax])
box on
title('Initial Condition for first model','Interpreter','latex')

subplot(1,3,2)
show_fem_transparent_edges(imgout{1});
plot_sensors(imgout{1});
view(0,90)
axis([-rmax rmax -rmax rmax])
box on
title(sprintf('Opt config %i',1),'Interpreter','latex')

subplot(1,3,3)
imgtemp = imgh;
imgtemp = assign_sensor_locations(imgtemp,sensor_positions_boundary);
show_fem_transparent_edges(imgtemp);
plot_sensors(imgtemp);
title('Boundary config','Interpreter','latex')
view(0,90)
axis([-rmax rmax -rmax rmax])
box on

drawnow
%% Plots

figure

sensor_position_init = zeros(num_of_sensors,3);
sensor_position_bfgs =  zeros(num_of_sensors,3);

for n = 1:num_of_sensors
    sensor_position_init(n,:) = imgsh{1}.fwd_model.sensors(n).position;
    sensor_position_bfgs(n,:) = imgout{1}.fwd_model.sensors(n).position;
end

hold on;
plot3(sensor_position_init(:,1),sensor_position_init(:,2),sensor_position_init(:,3),'r.','MarkerSize',10);
plot3(sensor_position_bfgs(:,1),sensor_position_bfgs(:,2),sensor_position_bfgs(:,3),'b.','MarkerSize',10);
plot3(sensor_positions_boundary(:,1),sensor_positions_boundary(:,2),sensor_positions_boundary(:,3),'g.','MarkerSize',10);
show_fem_transparent_edges(imgsh{1});

view(0,90)
axis([-rmax rmax -rmax rmax])
grid on; grid minor;
xlabel('x','Interpreter','latex');ylabel('y','Interpreter','latex')
axis([-rmax rmax -rmax rmax])
box on
legend({'Initial Condition','Optimal configuration BFGS','Boundary configuration'},'Interpreter','latex','Location','northwest')

drawnow

save(fullfile(script_folder,"\data\sensor_positions_init"),"sensor_position_init");
save(fullfile(script_folder,"\data\sensor_positions_opt"),"sensor_position_bfgs");

%% Compute forward solutions for initial config and optimal config

% Assign sensor locations and compute geometry matrices
imghinit = imgh;
imghinit = assign_sensor_locations(imghinit,sensor_position_init);
imghinit.fwd_model = compute_geometry_matrices(imghinit.fwd_model,mu0);

imgiinit = imgi;
imgiinit = assign_sensor_locations(imgiinit,sensor_position_init);
imgiinit.fwd_model = compute_geometry_matrices(imgiinit.fwd_model,mu0);

imghopt = imgh;
imghopt = assign_sensor_locations(imghopt,sensor_position_bfgs);
imghopt.fwd_model = compute_geometry_matrices(imghopt.fwd_model,mu0);

imgiopt = imgi;
imgiopt = assign_sensor_locations(imgiopt,sensor_position_bfgs);
imgiopt.fwd_model = compute_geometry_matrices(imgiopt.fwd_model,mu0);


imghbdy = imgh;
imghbdy = assign_sensor_locations(imghbdy,sensor_positions_boundary);
imghbdy.fwd_model = compute_geometry_matrices(imghbdy.fwd_model,mu0);

imgibdy = imgi;
imgibdy = assign_sensor_locations(imgibdy,sensor_positions_boundary);
imgibdy.fwd_model = compute_geometry_matrices(imgibdy.fwd_model,mu0);


Bi_init = fwd_solve_mdeit(imgiinit);
Bh_init = fwd_solve_mdeit(imghinit);
dB_init = [Bi_init.Bx(:)-Bh_init.Bx(:);Bi_init.By(:)-Bh_init.By(:);Bi_init.Bz(:)-Bh_init.Bz(:)];

Bi_opt = fwd_solve_mdeit(imgiopt);
Bh_opt = fwd_solve_mdeit(imghopt);
dB_opt = [Bi_opt.Bx(:)-Bh_opt.Bx(:);Bi_opt.By(:)-Bh_opt.By(:);Bi_opt.Bz(:)-Bh_opt.Bz(:)];

Bi_bdy = fwd_solve_mdeit(imgibdy);
Bh_bdy = fwd_solve_mdeit(imghbdy);
dB_bdy = [Bi_bdy.Bx(:)-Bh_bdy.Bx(:);Bi_bdy.By(:)-Bh_bdy.By(:);Bi_bdy.Bz(:)-Bh_bdy.Bz(:)];


% Add measurement noise
dy_init = add_measurement_noise_difference_absolute(...
    [Bi_init.Bx(:);Bi_init.By(:);Bi_init.Bz(:)],[Bh_init.Bx(:);Bh_init.By(:);Bh_init.Bz(:)],...
    noise_std_deviation);
dy_opt = add_measurement_noise_difference_absolute([Bi_opt.Bx(:);Bi_opt.By(:);Bi_opt.Bz(:)],...
    [Bh_opt.Bx(:);Bh_opt.By(:);Bh_opt.Bz(:)],...
    noise_std_deviation);

dy_bdy = add_measurement_noise_difference_absolute([Bi_bdy.Bx(:);Bi_bdy.By(:);Bi_bdy.Bz(:)],...
    [Bh_bdy.Bx(:);Bh_bdy.By(:);Bh_bdy.Bz(:)],...
    noise_std_deviation);

%% Plots
figure
subplot(1,2,1)
hold on
plot(dB_init,'r-')
plot(dB_opt,'b--')
plot(dB_bdy,'g-.')
xlabel('Measurement id','Interpreter','latex');
ylabel('dB','Interpreter','latex')
legend('Init','Opt','Bdy')

subplot(1,2,2)
hold on
plot(dy_init-dB_init,'r-')
plot(dy_opt-dB_opt,'b--')
plot(dy_bdy-dB_bdy,'g-.')
xlabel('Measurement id','Interpreter','latex');
ylabel('Noise','Interpreter','latex')
grid on;
grid minor;
box on;
legend('Init','Opt','Bdy')

%% Compute Jacobians for recon img

imgtemph = img_recon;
A = @(x) M(imgtemph,x);

imgtemph = assign_sensor_locations(imgtemph,sensor_position_init);
fprintf('Jacobian: Initial config\n')
[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtemph,A);
J_3axis_no_opt = [Jx;Jy;Jz];

imgtemph = assign_sensor_locations(imgtemph,sensor_position_bfgs);
fprintf('Jacobian: Opt config\n')
[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtemph,A);
J_3axis_opt = [Jx;Jy;Jz];

imgtemph = assign_sensor_locations(imgtemph,sensor_positions_boundary);
fprintf('Jacobian: Bdy config\n')
[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtemph,A);
J_3axis_bdy = [Jx;Jy;Jz];

%% Reconstruct posterior means pre and post sensor optimization



% fprintf('Reconstructing posterior means pre and post sensor optimization\n');
% 
% % We need only the diagonals of the posterior covariance matrix, so avoid
% % computing inv(...), it's very expensive. Its faster than the alternative
% % that chat gpt suggested tho...
% 
% % Is this form faster computationally????????????
% % Check!!!!!!!!!!!!!!!!!!!!!!!! (NOT SURE!!!!!!!! MAYBE CHANGE BACK TO Gamma_post_no_opt = inv(J_3axis_no_opt.'*inv_Gamma_noise_3_axis*J_3axis_no_opt + inv_Gamma_prior_3_axis);
% % Γpost = Γ_pr −Γ_pr*J^T*(J*Γ_pr*J^T +Γ_noise)^{−1}J*Γ_pr.
% 


imgtemph = imgh;
A = @(x) M(imgtemph,x);

imgtemph = assign_sensor_locations(imgtemph,sensor_position_init);
fprintf('Jacobian: Initial config\n')
[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtemph,A);
J_3axis_no_opt_post = [Jx;Jy;Jz];

imgtemph = assign_sensor_locations(imgtemph,sensor_position_bfgs);
fprintf('Jacobian: Opt config\n')
[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtemph,A);
J_3axis_opt_post = [Jx;Jy;Jz];

imgtemph = assign_sensor_locations(imgtemph,sensor_positions_boundary);
fprintf('Jacobian: Bdy config\n')
[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtemph,A);
J_3axis_bdy_post = [Jx;Jy;Jz];

fprintf('\t Recon post Gamma/mean no opt\n');
% Gamma_post_no_opt = inv(J_3axis_no_opt.'*inv_Gamma_noise_3_axis*J_3axis_no_opt + inv_Gamma_prior_3_axis);
Gamma_post_no_opt = Gamma_prior-(Gamma_prior*J_3axis_no_opt_post.')/(J_3axis_no_opt_post*Gamma_prior*J_3axis_no_opt_post.'+Gamma_noise_3_axis)*J_3axis_no_opt_post*Gamma_prior;

% Prior contribution is zero because mean of noise is 0. Avoid it so
% computation is faster
posterior_mean_no_opt = Gamma_post_no_opt*(...
    J_3axis_no_opt_post.'*inv_Gamma_noise*dy_init);

fprintf('\t Recon post Gamma/mean opt\n');

% Gamma_post_opt = inv(J_3axis_opt.'*inv_Gamma_noise_3_axis*J_3axis_opt + inv_Gamma_prior_3_axis);
Gamma_post_opt = Gamma_prior-(Gamma_prior*J_3axis_opt_post.')/(J_3axis_opt_post*Gamma_prior*J_3axis_opt_post.'+Gamma_noise_3_axis)*J_3axis_opt_post*Gamma_prior;

posterior_mean_opt = Gamma_post_opt*(...
    J_3axis_opt_post.'*inv_Gamma_noise*dy_opt);

fprintf('\t Recon post Gamma/mean boundary config\n');

Gamma_post_bdy = Gamma_prior-(Gamma_prior*J_3axis_bdy_post.')/(J_3axis_bdy_post*Gamma_prior*J_3axis_bdy_post.'+Gamma_noise_3_axis)*J_3axis_bdy_post*Gamma_prior;

posterior_mean_bdy = Gamma_post_bdy*(...
    J_3axis_bdy_post.'*inv_Gamma_noise*dy_bdy);

%%

imgtemp = imgh;
figure('Position',[100,100,800,500])
subplot(1,3,1);

imgtemp.elem_data = posterior_mean_no_opt;
imgtemp.calc_colours.ref_level = 0;
imgtemp.calc_colours.transparency_thresh = 0.1;

show_fem_transparent_edges(imgtemp)
plot_sensors(imghinit);
view(0,90)
axis([-R0 R0 -R0 R0])
title('Post Mean - Initial configuration','Interpreter','latex')
subplot(1,3,2);
imgtemp.elem_data = posterior_mean_opt;
show_fem_transparent_edges(imgtemp)
plot_sensors(imghopt);
view(0,90)
axis([-R0 R0 -R0 R0])
title('Post Mean - Optimal configuration','Interpreter','latex')
subplot(1,3,3);
imgtemp.elem_data = posterior_mean_bdy;
show_fem_transparent_edges(imgtemp)
plot_sensors(imghbdy);
view(0,90)
axis([-R0 R0 -R0 R0])
title('Post Mean - Boundary configuration','Interpreter','latex')

figure

mean_bdy_background = mean(posterior_mean_bdy(imgtemp.fwd_model.mat_idx{1}));
mean_opt_background = mean(posterior_mean_opt(imgtemp.fwd_model.mat_idx{1}));
mean_no_opt_background = mean(posterior_mean_no_opt(imgtemp.fwd_model.mat_idx{1}));


mean_bdy_anomaly = mean(posterior_mean_bdy(imgtemp.fwd_model.mat_idx{2}));
mean_opt_anomaly = mean(posterior_mean_opt(imgtemp.fwd_model.mat_idx{2}));
mean_no_opt_anomaly = mean(posterior_mean_no_opt(imgtemp.fwd_model.mat_idx{2}));


fprintf('Conductivity change:\n')
fprintf('\t Expected| Anomaly: %g, Background: %g\n',anomaly_conductivity-background_conductivity,0);
fprintf('\t Bdy| Anomaly %g, Background: %g\n',mean_bdy_anomaly,mean_bdy_background);
fprintf('\t Opt| Anomaly %g, Background: %g\n',mean_opt_anomaly,mean_opt_background);
fprintf('\t No-opt| Anomaly %g, Background: %g\n',mean_no_opt_anomaly,mean_no_opt_background);

error_bdy_anomaly = abs(mean_bdy_anomaly - anomaly_conductivity+background_conductivity)/abs(anomaly_conductivity-background_conductivity)*100;
error_opt_anomaly = abs(mean_opt_anomaly - anomaly_conductivity+background_conductivity)/abs(anomaly_conductivity-background_conductivity)*100;
error_no_opt_anomaly = abs(mean_no_opt_anomaly- anomaly_conductivity+background_conductivity)/abs(anomaly_conductivity-background_conductivity)*100;

error_bdy_background = abs(mean_bdy_background);
error_opt_background = abs(mean_opt_background);
error_no_opt_background = abs(mean_no_opt_background);


fprintf('Errors (%%):\n')
fprintf('\t Bdy| Anomaly %g %%, Background: %g\n',error_bdy_anomaly,error_bdy_background);
fprintf('\t Opt| Anomaly %g %%, Background: %g\n',error_opt_anomaly,error_opt_background);
fprintf('\t No-opt| Anomaly %g %%, Background: %g\n',error_no_opt_anomaly,error_no_opt_background);

hold on
plot(posterior_mean_bdy)
plot(posterior_mean_opt)
plot(posterior_mean_no_opt)
hold off
legend('bdy','opt','no-opt');
xlabel('elem id','Interpreter','latex')
ylabel('$\sigma$','Interpreter','latex')

%% Jacobian SVDs, need for later functions
fprintf('\t Jacobian SVD: Initial config\n')
[U_opt,S_opt,V_opt] = svd(J_3axis_opt,'econ');
fprintf('\t Jacobian SVD: Opt config\n')
[U_no_opt,S_no_opt,V_no_opt] = svd(J_3axis_no_opt,'econ');
fprintf('\t Jacobian SVD: Bdy config\n')
[U_bdy,S_bdy,V_bdy] = svd(J_3axis_bdy,'econ');

%% Reconstruct with GCV
function [lambda_opt,sigma] = generalized_cross_validation(J,data,lambda_vector,U,S,V)

m = size(J,1);
n = length(lambda_vector);

if nargin <4
    [U,S,V] = svds(J,m);
end

assert(exist('U','var'));
assert(exist('S','var'))
assert(exist('V','var'))

sigma = diag(S);             % singular values
Uy = U' * data;              % coordinates of data in U-basis

V_lambda = zeros(n,1);

for i = 1:n
    lambda = lambda_vector(i);

    gamma = sigma.^2 ./ (sigma.^2 + m*lambda);    % shrinkage factors

    one_minus_gamma = 1 - gamma;

    numerator = (1/m) * sum( (one_minus_gamma .* Uy).^2 );

    denominator = ( (1/m) * sum(one_minus_gamma) )^2;

    V_lambda(i) = numerator / denominator;
end

optimal_id = find(V_lambda == min(V_lambda));

lambda_opt = m*lambda_vector(optimal_id);

temp = diag(S)+lambda_opt./diag(S);
sigma = V*diag(1./temp)*U'*data;

end

lambda_vector = logspace(-20,0,20);
fprintf('\t For optimal config\n');
[lambda_opt,sigma_gcv_opt] = generalized_cross_validation(J_3axis_opt,dy_opt,lambda_vector,U_opt,S_opt,V_opt);
fprintf('\t For initial config\n');
[lambda_no_opt,sigma_gcv_no_opt] = generalized_cross_validation(J_3axis_no_opt,dy_init,lambda_vector,U_no_opt,S_no_opt,V_no_opt);
fprintf('\t For boundary config\n');
[lambda_bdy,sigma_gcv_bdy] = generalized_cross_validation(J_3axis_bdy,dy_bdy,lambda_vector,U_bdy,S_bdy,V_bdy);

%% Setup noise correction 
fprintf('Doing noise correction\n')


fprintf('\t Computing geometry matrices once\n')
imgh_nbc = compute_gamma_matrices_local_optimized(imgh);
imgi_nbc = compute_gamma_matrices_local_optimized(imgi);

function dB = noisy_data_generator_mdeit_3(imgh,imgi,noise_amplitude,num_noise_repetitions)
% Forward solve
[datah,~] = fwd_solve_mdeit(imgh);
Bh = [datah.Bx(:);datah.By(:);datah.Bz(:)];

% Forward solve
[datai,~] = fwd_solve_mdeit(imgi);
Bi = [datai.Bx(:);datai.By(:);datai.Bz(:)];

dB = zeros(size(Bi,1),num_noise_repetitions);

for t = 1:num_noise_repetitions
    dB(:,t) = add_measurement_noise_difference_absolute(Bi,Bh,noise_amplitude);
end

end

%% Noise correction 

num_noise_repetitions = 100;

func_mdeit = @(imgh,imgi,num_noise_repetitions) ...
    noisy_data_generator_mdeit_3(imgh_nbc,imgi_nbc,noise_std_deviation,num_noise_repetitions);

fprintf('\t Running NBC for opt config\n')
std_sigma_mdeit_opt = noise_correction(imgh,imgi,J_3axis_opt,lambda_opt,func_mdeit,num_noise_repetitions,U_opt,S_opt,V_opt);
fprintf('\t Running NBC for initial config\n')
std_sigma_mdeit_no_opt = noise_correction(imgh,imgi,J_3axis_no_opt,lambda_no_opt,func_mdeit,num_noise_repetitions,U_no_opt,S_no_opt,V_no_opt);
fprintf('\t Running NBC for boundary config\n')
std_sigma_mdeit_bdy = noise_correction(imgh,imgi,J_3axis_bdy,lambda_bdy,func_mdeit,num_noise_repetitions,U_bdy,S_bdy,V_bdy);

%% Plots

imgtemp = img_recon;
figure('Position',[100,100,800,500])
subplot(3,3,1);
imgtemp.elem_data = sigma_gcv_no_opt;
imgtemp.calc_colours.ref_level = mean(sigma_gcv_no_opt);
imgtemp.calc_colours.transparency_thresh = 0.2;

show_fem_transparent_edges(imgtemp)
plot_sensors(imghinit);
view(0,90)
title('GCV - Initial configuration','Interpreter','latex')
subplot(3,3,2);
imgtemp.elem_data = sigma_gcv_opt;
show_fem_transparent_edges(imgtemp)
plot_sensors(imghopt);
view(0,90)
title('GCV - Optimal configuration','Interpreter','latex')
subplot(3,3,3);
imgtemp.elem_data = sigma_gcv_bdy;
show_fem_transparent_edges(imgtemp)
plot_sensors(imghbdy);
view(0,90)
title('GCV - Boundary configuration','Interpreter','latex')

subplot(3,3,4);
imgtemp.elem_data = sigma_gcv_no_opt./std_sigma_mdeit_no_opt;
show_fem_transparent_edges(imgtemp)
% plot_sensors(imghinit);
view(0,90)
grid on;box on;
title('NBC - Initial configuration','Interpreter','latex')
subplot(3,3,5);
imgtemp.elem_data = sigma_gcv_opt./std_sigma_mdeit_opt;
show_fem_transparent_edges(imgtemp)
% plot_sensors(imghopt);
view(0,90)
grid on;box on;
title('NBC - Optimal configuration','Interpreter','latex')
subplot(3,3,6);
imgtemp.elem_data = sigma_gcv_bdy./std_sigma_mdeit_bdy;
show_fem_transparent_edges(imgtemp)
% plot_sensors(imghbdy);
view(0,90)
grid on;box on;
title('NBC - Boundary configuration','Interpreter','latex')


subplot(3,3,7);
show_fem_transparent_edges(imgi)
view(0,90)
title('Ground Truth','Interpreter','latex')

%% 

figure('Position',[100,100,1400,300])
subplot(1,4,1);
imgtemp.elem_data = sigma_gcv_no_opt;
show_fem_transparent_edges(imgtemp)
plot_sensors(imghinit);
view(0,90)
axis([-R0 R0 -R0 R0])
grid on;box on;
title('GCV - Initial configuration','Interpreter','latex')
subplot(1,4,2);
imgtemp.elem_data = sigma_gcv_opt;
show_fem_transparent_edges(imgtemp)
plot_sensors(imghopt);
view(0,90)
axis([-R0 R0 -R0 R0])
grid on;box on;
title('GCV - Optimal configuration','Interpreter','latex')
subplot(1,4,3);
imgtemp.elem_data = sigma_gcv_bdy;
show_fem_transparent_edges(imgtemp)
plot_sensors(imghbdy);
view(0,90)
axis([-R0 R0 -R0 R0])
grid on;box on;
title('GCV - Boundary configuration','Interpreter','latex')
subplot(1,4,4);
show_fem_transparent_edges(imgi)
view(0,90)
axis([-R0 R0 -R0 R0])
title('Ground Truth','Interpreter','latex')


figure('Position',[100,100,1400,300])
subplot(1,4,1);
imgtemp.elem_data = sigma_gcv_no_opt./std_sigma_mdeit_no_opt;
show_fem_transparent_edges(imgtemp)
plot_sensors(imghinit);
view(0,90)
axis([-R0 R0 -R0 R0])
grid on;box on;
title('NBC - Initial configuration','Interpreter','latex')
subplot(1,4,2);
imgtemp.elem_data = sigma_gcv_opt./std_sigma_mdeit_opt;
show_fem_transparent_edges(imgtemp)
plot_sensors(imghopt);
view(0,90)
axis([-R0 R0 -R0 R0])
grid on;box on;
title('NBC - Optimal configuration','Interpreter','latex')
subplot(1,4,3);
imgtemp.elem_data = sigma_gcv_bdy./std_sigma_mdeit_bdy;
show_fem_transparent_edges(imgtemp)
plot_sensors(imghbdy);
view(0,90)
axis([-R0 R0 -R0 R0])
grid on;box on;
title('NBC - Boundary configuration','Interpreter','latex')
subplot(1,4,4);
show_fem_transparent_edges(imgi)
view(0,90)
axis([-R0 R0 -R0 R0])
title('Ground Truth','Interpreter','latex')

%% Plot functions

function show_fem_transparent_edges(img)

hh = show_fem(img);                % draw the model (hh may be a handle or array)
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
% Apply transparency
set(patches, 'EdgeAlpha', 0.1);

% Disable legend participation
hh.Annotation.LegendInformation.IconDisplayStyle = 'off';
hh.HandleVisibility = 'off';

end

%% Helper functions

%function: Assign sensor locations
function img = assign_sensor_locations(img,sensor_locations)
assert(numel(img.fwd_model.sensors) == size(sensor_locations,1));
for id = 1: numel(img.fwd_model.sensors)
    img.fwd_model.sensors(id).position = sensor_locations(id,:);
end
end

function sensor_locations = fetch_sensor_locations(mdl)

assert(isfield(mdl,'type'));

switch mdl.type
    case 'image'
        fmdl = mdl.fwd_model;
    case 'fwd_model'
        fmdl = mdl;
end

n_sensors = numel(fmdl.sensors);
sensor_locations = zeros(n_sensors,3);

for m = 1: numel(fmdl.sensors)
    sensor_locations(m,:) = fmdl.sensors(m).position;
end

end

%% FUNCTIONS
function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
end

function [Rx,Ry,Rz,fmdl] = compute_r_matrices_local_optimized(fmdl,sensor_locations)

numElements = size(fmdl.elems,1);
numSensors = size(sensor_locations,1);

% ------------------------------------------------------------
% Persistent quadrature rule (avoid reallocation every call)
% ------------------------------------------------------------
persistent ref_pts weights nq;

if isempty(ref_pts)

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
    ref_pts = coord.';      % 3 × nq
    nq = size(ref_pts,2);
end

nodes = fmdl.nodes;
elems = fmdl.elems;

% ------------------------------------------------------------
% Precompute element geometry (done once!)
% ------------------------------------------------------------

v1  = zeros(3,numElements);
Jall = zeros(3,3,numElements);
detJ = zeros(numElements,1);

for k = 1:numElements

    verts = nodes(elems(k,:),:);

    v1(:,k) = verts(1,:).';

    J = [ ...
        (verts(2,:) - verts(1,:))' , ...
        (verts(3,:) - verts(1,:))' , ...
        (verts(4,:) - verts(1,:))' ];

    Jall(:,:,k) = J;
    detJ(k)     = abs(det(J));
end

% ------------------------------------------------------------
% Allocate outputs
% ------------------------------------------------------------

Rx = zeros(numSensors,numElements);
Ry = zeros(numSensors,numElements);
Rz = zeros(numSensors,numElements);

% ------------------------------------------------------------
% Main loop over sensors (parallelizable)
% ------------------------------------------------------------
parfor m = 1:numSensors   % <-- change to "for" if no Parallel Toolbox

    rm = sensor_locations(m,:).';

    Rx_local = zeros(1,numElements);
    Ry_local = zeros(1,numElements);
    Rz_local = zeros(1,numElements);

    for k = 1:numElements

        % Map quadrature points to physical tetrahedron
        Xi = v1(:,k) + Jall(:,:,k) * ref_pts;    % 3 × nq

        diff = rm - Xi;                         % 3 × nq
        r2   = sum(diff.^2,1);                  % 1 × nq

        kernel = diff ./ (r2.^1.5);             % 3 × nq

        R = (detJ(k)/6) * (kernel * weights);  % 3 × 1

        Rx_local(k) = R(1);
        Ry_local(k) = R(2);
        Rz_local(k) = R(3);
    end

    Rx(m,:) = Rx_local;
    Ry(m,:) = Ry_local;
    Rz(m,:) = Rz_local;
end

% ------------------------------------------------------------
% Store
% ------------------------------------------------------------

fmdl.R.Rx = Rx;
fmdl.R.Ry = Ry;
fmdl.R.Rz = Rz;

end

function img = compute_gamma_matrices_local_optimized(img)

mu_factor = img.fwd_model.mu0/(4*pi);

num_sensors = numel(img.fwd_model.sensors);

sensor_locations = zeros(numel(img.fwd_model.sensors),3);

for i = 1: numel(img.fwd_model.sensors)
    sensor_locations(i,:) = img.fwd_model.sensors(i).position;
end

[Rx,Ry,Rz,fmdl] = compute_r_matrices_local_optimized(img.fwd_model,sensor_locations);

img.fwd_model = fmdl; % THIS IS CRUCIAL!!!!! so the R matrices are updated whenever compute_gamma_matrices_local is called, otherwise they will be the same as the initial img.fwd_model, which might be wrong!

% Convenience handles
R.Rx = Rx;
R.Ry = Ry;
R.Rz = Rz;
G = img.fwd_model.G;

% Sigma = sparse(1:length(img.elem_data), 1:length(img.elem_data), img.elem_data);
Sigma = spdiags(img.elem_data(:), 0, length(img.elem_data), length(img.elem_data));

% NEW: The matrix g_{dl}^m, gives the components of the measurement axis of
% sensor m on the canonical R^3 basis

g = zeros(num_sensors,3,3);
for m = 1:numel(img.fwd_model.sensors)
    g(m,:,:) = [...
        img.fwd_model.sensors(m).axes.axis1;
        img.fwd_model.sensors(m).axes.axis2;
        img.fwd_model.sensors(m).axes.axis3];
end

Cx = ( -R.Rz * Sigma * G.Gy +  R.Ry * Sigma * G.Gz );
Cy = ( -R.Rx * Sigma * G.Gz +  R.Rz * Sigma * G.Gx );
Cz = ( -R.Ry * Sigma * G.Gx +  R.Rx * Sigma * G.Gy );

Gamma1 = mu_factor*(g(:,1,1).*Cx + g(:,1,2).*Cy + g(:,1,3).*Cz);
Gamma2 = mu_factor*(g(:,2,1).*Cx + g(:,2,2).*Cy + g(:,2,3).*Cz);
Gamma3 = mu_factor*(g(:,3,1).*Cx + g(:,3,2).*Cy + g(:,3,3).*Cz);

img.Gamma1 = Gamma1;
img.Gamma2 = Gamma2;
img.Gamma3 = Gamma3;

end

function [Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A)

mu0 = img.fwd_model.mu0;

n_nodes =  size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices_local_optimized(img);

R = img.fwd_model.R;
G = img.fwd_model.G;

Gamma1 = img.Gamma1;
Gamma2 = img.Gamma2;
Gamma3 = img.Gamma3;

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

A_matrix = A(img.elem_data);

Gamma1T = Gamma1.';
Gamma2T = Gamma2.';
Gamma3T = Gamma3.';

% Solve the adjoint problem for each sensor to get lambda vectors

lambdaX = A_matrix \ (-Gamma1T);
lambdaY = A_matrix \ (-Gamma2T);
lambdaZ = A_matrix \ (-Gamma3T);

Gx_times_lambda_X = G.Gx*lambdaX;
Gy_times_lambda_X = G.Gy*lambdaX;
Gz_times_lambda_X = G.Gz*lambdaX;

Gx_times_lambda_Y = G.Gx*lambdaY;
Gy_times_lambda_Y = G.Gy*lambdaY;
Gz_times_lambda_Y = G.Gz*lambdaY;

Gx_times_lambda_Z = G.Gx*lambdaZ;
Gy_times_lambda_Z = G.Gy*lambdaZ;
Gz_times_lambda_Z = G.Gz*lambdaZ;

Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% We want to broadcast arrays into num_sensors*num_stim*num_elems, in that order, so we avoid a permute in dfd

% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [1 1 n_elem]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u.', [1 num_stim n_elem]); % [: × numStim × 1]
GyU = reshape(Gy_times_u.', [1 num_stim n_elem]);
GzU = reshape(Gz_times_u.', [1 num_stim n_elem]);

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx, [num_sensors 1 n_elem]);
Ry_ = reshape(R.Ry, [num_sensors 1 n_elem]);
Rz_ = reshape(R.Rz, [num_sensors 1 n_elem]);

% These are the derivatives with respect to sigma of the C components of the Gamma matrix,
dCxdp = ( -Rz_.*GyU + Ry_.*GzU );
dCydp = ( -Rx_.*GzU + Rz_.*GxU );
dCzdp = ( -Ry_.*GxU + Rx_.*GyU );

% The g matrix does not depend on sigma.
g = zeros(num_sensors,3,3);
for m = 1:numel(img.fwd_model.sensors)
    g(m,:,:) = [...
        img.fwd_model.sensors(m).axes.axis1;
        img.fwd_model.sensors(m).axes.axis2;
        img.fwd_model.sensors(m).axes.axis3];
end

for select_sensor_axis = 1:3

    switch select_sensor_axis
        case 1
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_X.', [num_sensors 1 n_elem]); 
            GyL = reshape(Gy_times_lambda_X.', [num_sensors 1 n_elem]);
            GzL = reshape(Gz_times_lambda_X.', [num_sensors 1 n_elem]);
        case 2
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_Y.', [num_sensors 1 n_elem]); 
            GyL = reshape(Gy_times_lambda_Y.', [num_sensors 1 n_elem]);
            GzL = reshape(Gz_times_lambda_Y.', [num_sensors 1 n_elem]);
        case 3
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_Z.', [num_sensors 1 n_elem]); 
            GyL = reshape(Gy_times_lambda_Z.', [num_sensors 1 n_elem]);
            GzL = reshape(Gz_times_lambda_Z.', [num_sensors 1 n_elem]);
    end

    % Compute all dfdx for all sensors+stim
    dfdx = elemV .* ( ...
        GxL.*GxU + ...
        GyL.*GyU + ...
        GzL.*GzU );



    % g: [num_sensors × 3 × 3]
    gx = reshape(g(:,select_sensor_axis,1), [num_sensors 1 1 ]);
    gy = reshape(g(:,select_sensor_axis,2), [num_sensors 1 1]);
    gz = reshape(g(:,select_sensor_axis,3), [num_sensors 1 1]);

    dfdp = mu_factor*(...
        gx.*dCxdp +...
        gy.*dCydp +...
        gz.*dCzdp);

    dfd = dfdx + dfdp;   % size: [numStim × numSensors × numElems]

    % Now reshape to match J(ids,:)

    % collapse first 2 dims → [numSensors*numStim × numElems]
    J = reshape(dfd, num_sensors*num_stim, n_elem);

    switch select_sensor_axis
        case 1
            Jx = J;
        case 2
            Jy = J;
        case 3
            Jz = J;
    end

end

return
end

function fmdl = assign_magnetometers(fmdl,sensor_positions)

num_sensors = size(sensor_positions,1);

sensor_axes = ...
    repmat(struct('axis1', [], 'axis2', [],'axis3',[]), 1, num_sensors);

for m = 1:num_sensors
    sensor_axes(m).axis1 = [1,0,0];
    sensor_axes(m).axis2 = [0,1,0];
    sensor_axes(m).axis3 = [0,0,1];
end

% Assign sensorLocations and sensorAxes to fmdl
sensors = repmat(struct('position', [], 'axes', []), 1, num_sensors);

for m = 1:num_sensors
    sensors(m).position = sensor_positions(m,:);
    sensors(m).axes = sensor_axes(m);
end

fmdl.sensors = sensors;

end


function data_noisy = add_measurement_noise_difference_absolute(datai,datah,noise_amplitude)

if isempty(noise_amplitude)
    data_noisy = datai-datah;
    return;
else
    d = datai-datah;

    data_noisy = d+noise_amplitude*randn(size(d));
end
end




function sensor_position_boundary = find_sensor_positions_boundary(fmdl,delta)
    
    num_sensors = numel(fmdl.sensors);
    sensor_positions = fetch_sensor_locations(fmdl);
    thetas = atan2(sensor_positions(:,2),sensor_positions(:,1));
    z = sensor_positions(:,3);

    sensor_position_boundary = zeros(num_sensors,3);
    
    for m = 1:numel(thetas)
        ray_center = [0,0,z(m)];
        [hit_point,d] = find_hit_point_fmdl(fmdl,ray_center,thetas(m));
        
        sensor_position_boundary(m,:) = hit_point + delta*d;
    end

    
end

function [hit_point,d] = find_hit_point_fmdl(fmdl,ray_center,theta)

    assert(all(size(ray_center)==[1,3]))
    boundary = fmdl.boundary;

    % figure
    % show_fem_transparent_edges(fmdl);
    
    % Direction (horizontal plane)
    d = [cos(theta), sin(theta), 0];

    % Normalize direction (good practice)
    d = d / norm(d);
    
    % % Plot
    % % Choose a length for visualization
    % t_max = 1.2;  % adjust depending on mesh size

    % % Parametric ray
    % t = linspace(0, t_max, 100);
    % ray = ray_center + t.' * d;

    % hold on; grid on;
    % plot3(ray(:,1), ray(:,2), ray(:,3), 'r-', 'LineWidth', 2);
    % % Plot origin
    % plot3(ray_center(1), ray_center(2), ray_center(3), 'ro', 'MarkerFaceColor', 'r');

    t_min = inf;
    hit_point = [];

    for i = 1:size(boundary,1)
        tri = fmdl.nodes(boundary(i,:), :);
        v1 = tri(1,:); v2 = tri(2,:); v3 = tri(3,:);

        e1 = v2 - v1;
        e2 = v3 - v1;

        h = cross(d, e2);
        a = dot(e1, h);

        if abs(a) < 1e-12
            continue; % parallel
        end

        f = 1.0 / a;
        s = ray_center - v1;
        u = f * dot(s, h);

        if u < 0 || u > 1
            continue;
        end

        q = cross(s, e1);
        v = f * dot(d, q);

        if v < 0 || (u + v) > 1
            continue;
        end

        t = f * dot(e2, q);

        if t > 0 && t < t_min
            hit_point = ray_center + t * d;
            % plot3(hit_point(1), hit_point(2), hit_point(3), 'b.', 'MarkerSize',15);
            break
        end
    end


end




%% Test


function cost = compute_cost_function_a_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A)

n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);
% [Jx,Jy,Jz] = calc_jacobian_local(img,A);
J = [Jx;Jy;Jz];

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% This is generally better than the alternatives proposed by chatgpt
L = chol(H,'lower');
Hinv = L'\(L\eye(n_elem));
cost = trace(Hinv);

end