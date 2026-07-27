clc;clear all;close all
%% Prepare workspace
% Get the full path of the current script
fullpath = mfilename('fullpath');
% Extract just the folder
script_folder = fileparts(fullpath);
cd(script_folder);

% Have to add the functions path manually so prepare_workspace runs
grandparent_folder =fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));

model_folder = prepare_workspace(script_folder);

rng(1);

%% Model parameters 
z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);

model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/8;
model_parameters.numOfElectrodesPerRing = 8;
model_parameters.numOfRings = 4;
model_parameters.numOfSensors = model_parameters.numOfElectrodesPerRing*model_parameters.numOfRings;
model_parameters.sensorRadius = model_parameters.radius*1.5;
model_parameters.material.type = 'spherical';
model_parameters.material.position(1) = 0.95*(model_parameters.radius-model_parameters.material.radius);
model_parameters.material.position(3) = 0.5*model_parameters.height;

background_conductivity = 3.28e-1/sigma0;
anomaly_conductivity = background_conductivity*1.1;

%% Simulation parameters
dim = 3;            %if doing 1-axis, which dimmension?

%Optimizer parameters
max_iterations = 100;
use_parallel = true;

rmax = 2*model_parameters.radius;
rmin = model_parameters.radius*1.01;

zmax = model_parameters.height;
zmin = 0;

R0 = model_parameters.radius*1.5;

alpha = 0; %controls the force of the repulsion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

n_trial = 1; % number of initial conditions

snr = 30;
%% Make forward model with material

% Forward model with material 
model_parameters.maxsz = model_parameters.radius/10;
[~,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl = fmdls{1};

% Make homogeneous image
imgh = mk_image_mdeit(fmdl,background_conductivity);
% Add plastic cylinder
imgi = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

%% Make forward model without material for reconstruction
model_parameters_temp = model_parameters;
model_parameters_temp.maxsz = model_parameters.radius/2;

model_parameters_temp.material = struct();

% Forward model without material 
[~,fmdls] = mk_mdeit_model(model_parameters_temp,model_folder);
fmdl_recon = fmdls{1};

imgh_recon = mk_image_mdeit(fmdl_recon,background_conductivity);



%%
n_sensors = numel(imgh_recon.fwd_model.sensors);
n_stim = numel(imgh_recon.fwd_model.stimulation);
n_elem = size(imgh_recon.fwd_model.elems,1);
n_nodes = size(fmdl.nodes,1);

A = @(x) M(imgh_recon,x);





%% Initialize
r0 = R0*ones(1,n_sensors);

all_sensor_locations_0 = zeros(n_sensors,3,n_trial);

% Set multiple initial conditions
for i = 1:n_trial
    theta_new = 2*pi*rand(n_sensors,1);
    z_new = model_parameters.height/5 + 4/5*model_parameters.height*rand(1,n_sensors);

    sensor_locations_0 = [(r0(:).*cos(theta_new(:))),(r0(:).*sin(theta_new(:))),z_new(:)];

    all_sensor_locations_0(:,:,i) = sensor_locations_0;
end

theta0 = 0:2*pi/n_sensors:2*pi-2*pi/n_sensors;
zed0 = model_parameters.material.position(3)*ones(n_sensors,1);
zed1 = model_parameters.height/10*ones(n_sensors,1);

sensor_locations_circle = [r0(:).*cos(theta0(:)),r0(:).*sin(theta0(:)),zed0];
sensor_locations_circle_height = [r0(:).*cos(theta0(:)),r0(:).*sin(theta0(:)),zed1];
%% Define prior and noise covariance matrices

B = fwd_solve_mdeit(imgi);
max_B = max([abs(B.Bx(:));abs(B.By(:));abs(B.Bz(:))]);

% Set the noise variance with respect to the data magnitude
noise_std_deviation = max_B/3;
variance_noise = noise_std_deviation^2;

% Jdim = calc_jacobian_1axis_direct_fully_vectorized_local(imgh_recon,A,dim);
% 
% % This should be big enough so #rank eigenvalues are significant, that is,
% % we want d = rank(Jdim) and d_3axis = rank(J_3axis)
% coeff = 100;
% variance_prior = coeff*variance_noise/eigs(Jdim'*Jdim,1);
% 
% % Check how many eigenvectors of J'J are in the data-dominated
% % regime, as opposed to the prior-dominated regime
% d = sum(eigs(Jdim'*Jdim,n_elem).*variance_prior/variance_noise>1);
% 
% fprintf('# lambda_i*alpha:beta>1 = %i (%i)\n',d,n_elem);
% 
% % White noise with zero mean and \mu variance
% Gamma_noise = variance_noise.*speye(n_stim*n_sensors,n_stim*n_sensors);
% inv_Gamma_noise = inv(Gamma_noise);
% 
% Gamma_prior = variance_prior.*speye(n_elem,n_elem);
% inv_Gamma_prior = inv(Gamma_prior);
% 
% [Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgh_recon,A);
% 
% J_3axis_opt = [Jx;Jy;Jz];
% 
% variance_prior_3_axis = coeff*variance_noise/eigs(J_3axis_opt'*J_3axis_opt,1);
% 
% d_3axis = sum(eigs(J_3axis_opt'*J_3axis_opt,n_elem).*variance_prior_3_axis/variance_noise>1);
% 
% fprintf('# lambda_i*alpha:beta>1 = %i (%i)\n',d_3axis,n_elem);
% 
% Gamma_prior_3_axis = variance_prior_3_axis.*speye(n_elem,n_elem);
% inv_Gamma_prior_3_axis = inv(Gamma_prior_3_axis);

Gamma_noise_3_axis = variance_noise.*speye(3*n_stim*n_sensors,3*n_stim*n_sensors);
inv_Gamma_noise_3_axis = inv(Gamma_noise_3_axis);

%% Smoothness prior covariance matrix

str = sprintf("(x-%2.2f).^2+(y-%2.2f).^2+(z-%2.2f).^2<%2.2f^2",...
    model_parameters.anomaly.position(1),model_parameters.anomaly.position(2),model_parameters.anomaly.position(3),model_parameters.anomaly.radius);

select_fcn = inline(str,'x','y','z');
idx = elem_select(imgh_recon.fwd_model, select_fcn);
idx(idx>0) = 1;
in_block  = idx * idx.';  % both in
out_block = (~idx) * (~idx).';% both out

% Assemble kappa
kappa = 0.4 * sparse(in_block) + 0.03 * sparse(out_block);

kappa = variance_noise * sparse(in_block) + variance_noise/10 * sparse(out_block);


Gamma_smooth = smooth_prior(imgh_recon,0.5,kappa);

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

        weights = val .* exp(-D2/(2*lambda^2));

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

% Gamma_prior = variance_prior_3_axis*Gamma_smooth;
Gamma_prior = Gamma_smooth;
Gamma_prior = Gamma_prior + 1e-6 * speye(size(Gamma_prior)); %need regularization because of ill-conditioning
inv_Gamma_prior = Gamma_prior \ speye(size(Gamma_prior));
inv_Gamma_prior_3_axis = inv_Gamma_prior;

figure
img_temp = imgh_recon;
img_temp.elem_data = diag(Gamma_prior);
show_fem(img_temp);

drawnow

% % DEBUG 
% variance_prior_3_axis = variance_noise*1000;
% Gamma_prior_3_axis = variance_prior_3_axis.*speye(n_elem,n_elem);
% inv_Gamma_prior_3_axis = inv(Gamma_prior_3_axis);


%% Define jacobian of coordinate transformations

function jacobian_coordinate_transformation = compute_jacobian_coordinate_transformation_cylindrical(sensor_locations)

n_sensors = size(sensor_locations,1);

jacobian_coordinate_transformation = zeros(3,3,n_sensors);

rm = sqrt(sensor_locations(:,1).^2+sensor_locations(:,2).^2);
thetam = atan2(sensor_locations(:,2),sensor_locations(:,1));
zm = sensor_locations(:,3);

jacobian_coordinate_transformation_r = zeros(3,n_sensors);
jacobian_coordinate_transformation_r(1,:) = cos(thetam);
jacobian_coordinate_transformation_r(2,:) = sin(thetam);
jacobian_coordinate_transformation_r(3,:) = zeros(1,n_sensors);

jacobian_coordinate_transformation_theta = zeros(3,n_sensors);
jacobian_coordinate_transformation_theta(1,:) = -rm'.*sin(thetam)';
jacobian_coordinate_transformation_theta(2,:) = rm'.*cos(thetam)';
jacobian_coordinate_transformation_theta(3,:) = zeros(1,n_sensors);

jacobian_coordinate_transformation_z = zeros(3,n_sensors);
jacobian_coordinate_transformation_z(1,:) = zeros(1,n_sensors);
jacobian_coordinate_transformation_z(2,:) = zeros(1,n_sensors);
jacobian_coordinate_transformation_z(3,:) = ones(1,n_sensors);

jacobian_coordinate_transformation(1,:,:) = jacobian_coordinate_transformation_r;
jacobian_coordinate_transformation(2,:,:) = jacobian_coordinate_transformation_theta;
jacobian_coordinate_transformation(3,:,:) = jacobian_coordinate_transformation_z;
end

% DEBUG!
%%%%%%%%
function jacobian_coordinate_transformation = compute_jacobian_coordinate_transformation_test(sensor_locations)

n_sensors = size(sensor_locations,1);

jacobian_coordinate_transformation = zeros(3,3,n_sensors);

jacobian_coordinate_transformation_x = zeros(3,n_sensors);
jacobian_coordinate_transformation_x(1,:) = ones(1,n_sensors);
jacobian_coordinate_transformation_x(2,:) = zeros(1,n_sensors);
jacobian_coordinate_transformation_x(3,:) = zeros(1,n_sensors);

jacobian_coordinate_transformation_y = zeros(3,n_sensors);
jacobian_coordinate_transformation_y(1,:) = zeros(1,n_sensors);
jacobian_coordinate_transformation_y(2,:) = ones(1,n_sensors);
jacobian_coordinate_transformation_y(3,:) = zeros(1,n_sensors);

jacobian_coordinate_transformation_z = zeros(3,n_sensors);
jacobian_coordinate_transformation_z(1,:) = zeros(1,n_sensors);
jacobian_coordinate_transformation_z(2,:) = zeros(1,n_sensors);
jacobian_coordinate_transformation_z(3,:) = ones(1,n_sensors);

jacobian_coordinate_transformation(1,1,:) = jacobian_coordinate_transformation_x(1,:);
jacobian_coordinate_transformation(1,2,:) = jacobian_coordinate_transformation_x(2,:);
jacobian_coordinate_transformation(1,3,:) = jacobian_coordinate_transformation_x(3,:);

jacobian_coordinate_transformation(2,1,:) = jacobian_coordinate_transformation_y(1,:);
jacobian_coordinate_transformation(2,2,:) = jacobian_coordinate_transformation_y(2,:);
jacobian_coordinate_transformation(2,3,:) = jacobian_coordinate_transformation_y(3,:);

jacobian_coordinate_transformation(3,1,:) = jacobian_coordinate_transformation_z(1,:);
jacobian_coordinate_transformation(3,2,:) = jacobian_coordinate_transformation_z(2,:);
jacobian_coordinate_transformation(3,3,:) = jacobian_coordinate_transformation_z(3,:);

end
%%%%%%%%

function jacobian_coordinate_transformation = compute_jacobian_coordinate_transformation_z(sensor_locations)
n_sensors = size(sensor_locations,1);

jacobian_coordinate_transformation = zeros(1,3,n_sensors);

jacobian_coordinate_transformation_z = zeros(3,n_sensors);
jacobian_coordinate_transformation_z(1,:) = zeros(1,n_sensors);
jacobian_coordinate_transformation_z(2,:) = zeros(1,n_sensors);
jacobian_coordinate_transformation_z(3,:) = ones(1,n_sensors);

jacobian_coordinate_transformation(1,1,:) = jacobian_coordinate_transformation_z(1,:);
jacobian_coordinate_transformation(1,2,:) = jacobian_coordinate_transformation_z(2,:);
jacobian_coordinate_transformation(1,3,:) = jacobian_coordinate_transformation_z(3,:);
end

function jacobian_coordinate_transformation = compute_jacobian_coordinate_transformation_r(sensor_locations)

n_sensors = size(sensor_locations,1);

jacobian_coordinate_transformation = zeros(1,3,n_sensors);

thetam = atan2(sensor_locations(:,2),sensor_locations(:,1));

jacobian_coordinate_transformation_r = zeros(3,n_sensors);
jacobian_coordinate_transformation_r(1,:) = cos(thetam);
jacobian_coordinate_transformation_r(2,:) = sin(thetam);
jacobian_coordinate_transformation_r(3,:) = zeros(1,n_sensors);

jacobian_coordinate_transformation(1,1,:) = jacobian_coordinate_transformation_r(1,:);
jacobian_coordinate_transformation(1,2,:) = jacobian_coordinate_transformation_r(2,:);
jacobian_coordinate_transformation(1,3,:) = jacobian_coordinate_transformation_r(3,:);
end

function jacobian_coordinate_transformation = compute_jacobian_coordinate_transformation_xi(sensor_locations,rmin,rmax)

n_sensors = size(sensor_locations,1);

jacobian_coordinate_transformation = zeros(1,3,n_sensors);

rm = sqrt(sensor_locations(:,1).^2+sensor_locations(:,2).^2);
thetam = atan2(sensor_locations(:,2),sensor_locations(:,1));
xim = log((rm-rmin)./(rmax-rm));

sigmoid = @(x) 1./(1+exp(-x));
dsigmoid = @(x) (sigmoid(x)).*(1-sigmoid(x));

drmdxim = (rmax-rmin)*dsigmoid(xim);

jacobian_coordinate_transformation_xi = zeros(3,n_sensors);
jacobian_coordinate_transformation_xi(1,:) = cos(thetam).*drmdxim;
jacobian_coordinate_transformation_xi(2,:) = sin(thetam).*drmdxim;
jacobian_coordinate_transformation_xi(3,:) = zeros(1,n_sensors);

jacobian_coordinate_transformation(1,1,:) = jacobian_coordinate_transformation_xi(1,:);
jacobian_coordinate_transformation(1,2,:) = jacobian_coordinate_transformation_xi(2,:);
jacobian_coordinate_transformation(1,3,:) = jacobian_coordinate_transformation_xi(3,:);
end

function jacobian_coordinate_transformation = compute_jacobian_coordinate_transformation_xi_eta(sensor_locations,rmin,rmax,zmin,zmax)

n_sensors = size(sensor_locations,1);

jacobian_coordinate_transformation = zeros(2,3,n_sensors);

rm = sqrt(sensor_locations(:,1).^2+sensor_locations(:,2).^2);
thetam = atan2(sensor_locations(:,2),sensor_locations(:,1));
zm = sensor_locations(:,3);

xim = log((rm-rmin)./(rmax-rm));
etam = log((zm-zmin)./(zmax-zm));

sigmoid = @(x) 1./(1+exp(-x));
dsigmoid = @(x) (sigmoid(x)).*(1-sigmoid(x));

drmdxim = (rmax-rmin)*dsigmoid(xim);
dzmdetam = (zmax-zmin)*dsigmoid(etam);

jacobian_coordinate_transformation_xi = zeros(3,n_sensors);
jacobian_coordinate_transformation_xi(1,:) = cos(thetam).*drmdxim;
jacobian_coordinate_transformation_xi(2,:) = sin(thetam).*drmdxim;
jacobian_coordinate_transformation_xi(3,:) = zeros(1,n_sensors);

jacobian_coordinate_transformation_eta = zeros(3,n_sensors);
jacobian_coordinate_transformation_eta(1,:) = zeros(1,n_sensors);
jacobian_coordinate_transformation_eta(2,:) = zeros(1,n_sensors);
jacobian_coordinate_transformation_eta(3,:) = dzmdetam;

jacobian_coordinate_transformation(1,1,:) = jacobian_coordinate_transformation_xi(1,:);
jacobian_coordinate_transformation(1,2,:) = jacobian_coordinate_transformation_xi(2,:);
jacobian_coordinate_transformation(1,3,:) = jacobian_coordinate_transformation_xi(3,:);

jacobian_coordinate_transformation(2,1,:) = jacobian_coordinate_transformation_eta(1,:);
jacobian_coordinate_transformation(2,2,:) = jacobian_coordinate_transformation_eta(2,:);
jacobian_coordinate_transformation(2,3,:) = jacobian_coordinate_transformation_eta(3,:);
end

% Transform gradient of cost function w.r.t. cartesian coordinates to
% w.r.t. q coordinates
function dphidq = dphidp_to_dphidq(sensor_locations,dphidp,jacobian_coordinate_transformation)

n_sensors = size(sensor_locations,1);

assert(size(jacobian_coordinate_transformation,2) == 3);
assert(size(jacobian_coordinate_transformation,3) == n_sensors);


qmax = size(jacobian_coordinate_transformation,1);

dphidq = zeros(qmax,n_sensors);

for q = 1:qmax
    temp = zeros(1,n_sensors);
    for dim = 1:3
        temp = temp + squeeze(jacobian_coordinate_transformation(q,dim,:)).'.*dphidp(dim,:);
    end
    dphidq(q,:) = temp;
end

end

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


% Map cartesian coordinates to cylindrical coordinates
function q = x_to_q_cylindrical(x)

assert(mod(numel(x),3)==0);
n_sensors = numel(x)/3;

rm = sqrt(x(1:n_sensors).^2+x(n_sensors+1:2*n_sensors).^2);
thetam = atan2(x(n_sensors+1:2*n_sensors),x(1:n_sensors));
zm = x(2*n_sensors+1:3*n_sensors);

q = [rm(:);thetam(:);zm(:)];
end

function x = q_to_x_cylindrical(q)

assert(mod(numel(q),3)==0);
n_sensors = numel(q)/3;

rm = q(1:n_sensors);
thetam = q(n_sensors+1:2*n_sensors);
zm = q(2*n_sensors+1:3*n_sensors);

x = [rm(:).*cos(thetam(:));rm(:).*sin(thetam(:));zm(:)];
end



% Map cartesian coordinates to parametrized region (region contained
% between two cylinders of radius rmin and rmax and height [0,3])
function q = map_x_to_q_cyl_region(x,rmin,rmax,zmin,zmax)
assert(mod(numel(x),3)==0);
n_sensors = numel(x)/3;

rm = sqrt(x(1:n_sensors).^2+x(n_sensors+1:2*n_sensors).^2);
thetam = atan2(x(n_sensors+1:2*n_sensors),x(1:n_sensors));
zm = x(2*n_sensors+1:3*n_sensors);


xim = log((rm-rmin)./(rmax-rm));
etam = log((zm-zmin)./(zmax-zm));

q = [xim(:);thetam(:);etam(:)];
end

function x = map_q_to_x_cyl_region(q,rmin,rmax,zmin,zmax)

sigmoid = @(y) 1./(1+exp(-y));

assert(mod(numel(q),3)==0);
n_sensors = numel(q)/3;

q = q(:);
xim = q(1:n_sensors);
thetam = q(n_sensors+1:2*n_sensors);
etam = q(2*n_sensors+1:3*n_sensors);

rm = rmin + (rmax-rmin).*sigmoid(xim);
zm = zmin + (zmax-zmin).*sigmoid(etam);

x = [rm(:).*cos(thetam(:));rm(:).*sin(thetam(:));zm(:)];
end

x_to_q_cyl_region = @(x) map_x_to_q_cyl_region(x,rmin,rmax,zmin,zmax);

q_to_x_cyl_region = @(q) map_q_to_x_cyl_region(q,rmin,rmax,zmin,zmax);



% Map cartesian coordinates to parametrized region (cylinder shell at radius r0)
function q = map_x_to_q_cyl_shell(x,r0,zmin,zmax)
assert(mod(numel(x),2)==0);
n_sensors = numel(x)/2;

thetam = atan2(x(n_sensors+1:2*n_sensors),x(1:n_sensors));
zm = x(2*n_sensors+1:3*n_sensors);

etam = log((zm-zmin)./(zmax-zm));

q = [thetam(:);etam(:)];
end

function x = map_q_to_x_cyl_shell(q,r0,zmin,zmax)
sigmoid = @(y) 1./(1+exp(-y));

assert(mod(numel(q),2)==0);
n_sensors = numel(q)/2;

q = q(:);
thetam = q(1:n_sensors);
etam = q(n_sensors+1:2*n_sensors);

zm = zmin + (zmax-zmin).*sigmoid(etam);

x = [r0.*cos(thetam(:));r0.*sin(thetam(:));zm(:)];
end

x_to_q_cyl_shell = @(x) map_x_to_q_cyl_shell(x,r0,zmin,zmax);
q_to_x_cyl_shell = @(q) map_q_to_x_cyl_shell(q,r0,zmin,zmax);



function x = map_q_to_x_z(q,r0,theta0)
n_sensors = numel(q);

zm = q(:);

x = [r0.'.*cos(theta0(:));r0.'.*sin(theta0(:));zm(:)];
end

function q = map_x_to_q_z(x)

assert(mod(numel(x),3)==0);
n_sensors = numel(x)/3;

zm = x(2*n_sensors+1:3*n_sensors);

q = [zm(:)];
end

x_to_q_z = @(x) map_x_to_q_z(x);
q_to_x_z = @(q) map_q_to_x_z(q,r0,theta0);



function x = map_q_to_x_xi_r(q,theta0,zed0,rmin,rmax)
n_sensors = numel(q);

xim = q(:);

sigmoid = @(x) 1./(1+exp(-x));

rm = rmin + (rmax-rmin).*sigmoid(xim);

x = [rm(:).*cos(theta0(:));rm(:).*sin(theta0(:));zed0(:)];
end

function q = map_x_to_q_xi_r(x,rmin,rmax)

assert(mod(numel(x),3)==0);
n_sensors = numel(x)/3;

xm = x(0*n_sensors+1:1*n_sensors);
ym = x(1*n_sensors+1:2*n_sensors);
rm = sqrt(xm.^2+ym.^2);

xim = log((rm-rmin)./(rmax-rm));

q = [xim(:)];
end

x_to_q_xi_r = @(x) map_x_to_q_xi_r(x,rmin,rmax);
q_to_x_xi_r = @(q) map_q_to_x_xi_r(q,theta0,zed0,rmin,rmax);



function x = map_q_to_x_xi_r_z(q,theta0,rmin,rmax,zmin,zmax)
assert(mod(numel(q),2)==0)
n_sensors = numel(q)/2;

xim = q(1:n_sensors);
etam = q(n_sensors+1:2*n_sensors);

sigmoid = @(x) 1./(1+exp(-x));

rm = rmin + (rmax-rmin).*sigmoid(xim);

zm = zmin + (zmax-zmin).*sigmoid(etam);

x = [rm(:).*cos(theta0(:));rm(:).*sin(theta0(:));zm(:)];
end

function q = map_x_to_q_xi_r_z(x,rmin,rmax,zmin,zmax)

assert(mod(numel(x),3)==0);
n_sensors = numel(x)/3;

xm = x(0*n_sensors+1:1*n_sensors);
ym = x(1*n_sensors+1:2*n_sensors);
zm = x(2*n_sensors+1:3*n_sensors);

rm = sqrt(xm.^2+ym.^2);

xim = log((rm-rmin)./(rmax-rm));
etam = log((zm-zmin)./(zmax-zm));

q = [xim(:);etam(:)];
end

x_to_q_xi_eta = @(x) map_x_to_q_xi_r_z(x,rmin,rmax,zmin,zmax);
q_to_x_xi_eta = @(q) map_q_to_x_xi_r_z(q,theta0,rmin,rmax,zmin,zmax);



% DEBUG!
x_to_q_test = @(x) x;
q_to_x_test = @(q) q;


%% Define vector_to_sensor_locations and vice-versa

x0 = sensor_locations_to_vector_cartesian(sensor_locations_0);


vector_to_sensor_locations_test = @(q) vector_to_sensor_locations_cartesian(q_to_x_test(q));
sensor_locations_to_vector_test = @(sensor_locations) x_to_q_test(sensor_locations_to_vector_cartesian(sensor_locations));



% Full map from q coordinates to sensor locations and back for cylindrical
% coordinates
vector_to_sensor_locations_con = @(q) vector_to_sensor_locations_cartesian(q_to_x_cylindrical(q));
sensor_locations_to_vector_con = @(sensor_locations) x_to_q_cylindrical(sensor_locations_to_vector_cartesian(sensor_locations));

% Full map from q coordinates to sensor locations and back for z-coordinates
vector_to_sensor_locations_z = @(q) vector_to_sensor_locations_cartesian(q_to_x_z(q));
sensor_locations_to_vector_z = @(sensor_locations) x_to_q_z(sensor_locations_to_vector_cartesian(sensor_locations));

% Full map from q coordinates to sensor locations and back for xi-coordinates
vector_to_sensor_locations_xi = @(q) vector_to_sensor_locations_cartesian(q_to_x_xi_r(q));
sensor_locations_to_vector_xi = @(sensor_locations) x_to_q_xi_r(sensor_locations_to_vector_cartesian(sensor_locations));

% Full map from q coordinates to sensor locations and back for (xi,eta)-coordinates
vector_to_sensor_locations_xi_eta = @(q) vector_to_sensor_locations_cartesian(q_to_x_xi_eta(q));
sensor_locations_to_vector_xi_eta = @(sensor_locations) x_to_q_xi_eta(sensor_locations_to_vector_cartesian(sensor_locations));


% Full map from q coordinates to sensor locations and back for
% unconstrained optimization
vector_to_sensor_locations_unc = @(q) vector_to_sensor_locations_cartesian(q_to_x_cyl_region(q));
sensor_locations_to_vector_unc = @(sensor_locations) x_to_q_cyl_region(sensor_locations_to_vector_cartesian(sensor_locations));

% % Sanity check
% q = x_to_q_cylindrical(x0);
% x = q_to_x_cylindrical(q);
% assert(norm(x-x0)<1e-5,'Unexpected');
%
% % Sanity check
% q0 = sensor_locations_to_vector_con(sensor_locations_0);
% sensor_locations_new = vector_to_sensor_locations_con(q0);
% assert(norm(sensor_locations_new-sensor_locations_0)<1e-5,'Unexpected');
%
% % Sanity check
% sensor_locations = vector_to_sensor_locations_con(q0);
% q_new = sensor_locations_to_vector_con(sensor_locations);
% assert(norm(q0-q_new)<1e-5,'Unexpected');


%% Functions to compute cost and cost gradient
function out = f(img,q,inv_Gamma_prior,inv_Gamma_noise,A,vector_to_sensor_locations,opt_mode,mode,dim)

allowed_opt_modes = {'a-opt','d-opt'};
assert(ismember(opt_mode,allowed_opt_modes));

allowed_modes = {'mdeit3','mdeit1'};
assert(ismember(mode,allowed_modes));

if strcmp(mode,'mdeit3') && nargin<9
    dim = 'default';
end

sensor_locations = vector_to_sensor_locations(q);

switch mode
    case 'mdeit1'
        switch opt_mode
            case 'a-opt'
                phi_oed = compute_cost_function_a_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
            case 'd-opt'
                phi_oed = compute_cost_function_d_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
        end
    case 'mdeit3'
        switch opt_mode
            case 'a-opt'
                phi_oed = compute_cost_function_a_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
            case 'd-opt'
                phi_oed = compute_cost_function_d_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
        end
end

out = phi_oed;
end

function dphi = g(img,q,inv_Gamma_prior,inv_Gamma_noise,A,...
    vector_to_sensor_locations,jac_coor_transf,opt_mode,mode,dim)

n_sensors = numel(img.fwd_model.sensors);

allowed_opt_modes = {'a-opt','d-opt'};
assert(ismember(opt_mode,allowed_opt_modes));

allowed_modes = {'mdeit3','mdeit1'};
assert(ismember(mode,allowed_modes));

if strcmp(mode,'mdeit3') && nargin<10
    dim = 'default';
end

sensor_locations = vector_to_sensor_locations(q);

switch mode
    case 'mdeit1'
        switch opt_mode
            case 'a-opt'
                dphidp = compute_cost_function_gradient_a_opt_optimized(...
                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
            case 'd-opt'
                dphidp = compute_cost_function_gradient_d_opt_optimized(...
                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
        end
    case 'mdeit3'
        switch opt_mode
            case 'a-opt'
                dphidp = compute_cost_function_gradient_a_opt_optimized_3_axis(...
                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
            case 'd-opt'
                dphidp = compute_cost_function_gradient_d_opt_optimized_3_axis(...
                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
        end
end

% Convert cartesian derivatives to other coordinate derivatives (
% generalized)
jacobian_coordinate_transformation = jac_coor_transf(sensor_locations);
dphidq = dphidp_to_dphidq(sensor_locations,dphidp,jacobian_coordinate_transformation);

dphi = reshape(dphidq.', 1, []);

end

jac_coor_transf_cylindrical = @(s) compute_jacobian_coordinate_transformation_cylindrical(s);
jac_coor_transf_z = @(s) compute_jacobian_coordinate_transformation_z(s);
jac_coor_transf_r = @(s) compute_jacobian_coordinate_transformation_r(s);
jac_coor_transf_xi = @(s) compute_jacobian_coordinate_transformation_xi(s,rmin,rmax);
jac_coor_transf_xi_eta = @(s) compute_jacobian_coordinate_transformation_xi_eta(s,rmin,rmax,zmin,zmax);


% Compute gradient w.r.t. z: use jac_coor_transf_z
f_a_opt_mdeit_dim_z = @(q) f(imgh_recon,q,inv_Gamma_prior,inv_Gamma_noise,A,vector_to_sensor_locations_z,'a-opt','mdeit1',dim);
g_a_opt_mdeit_dim_z  = @(q) g(imgh_recon,q,inv_Gamma_prior,inv_Gamma_noise,A,...
    vector_to_sensor_locations_z,jac_coor_transf_z,'a-opt','mdeit1',dim);

% Compute gradient w.r.t. r: use jac_coor_transf_r
f_a_opt_mdeit_dim_r = @(q) f(imgh_recon,q,inv_Gamma_prior,inv_Gamma_noise,A,vector_to_sensor_locations_xi,'a-opt','mdeit1',dim);
g_a_opt_mdeit_dim_r  = @(q) g(imgh_recon,q,inv_Gamma_prior,inv_Gamma_noise,A,...
    vector_to_sensor_locations_xi,jac_coor_transf_r,'a-opt','mdeit1',dim);

% 3D-mdeit for z
f_a_opt_mdeit3_z = @(q) f(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,vector_to_sensor_locations_z,'a-opt','mdeit3',dim);
g_a_opt_mdeit3_z  = @(q) g(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,...
    vector_to_sensor_locations_z,jac_coor_transf_z,'a-opt','mdeit3',dim);

% 3D-mdeit for xi
f_d_opt_mdeit3_xi = @(q) f(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,vector_to_sensor_locations_xi,'d-opt','mdeit3',dim);
g_d_opt_mdeit3_xi  = @(q) g(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,...
    vector_to_sensor_locations_xi,jac_coor_transf_xi,'d-opt','mdeit3',dim);

% 3D-mdeit for (xi,eta) - d opt
f_d_opt_mdeit3_xi_eta = @(q) f(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,vector_to_sensor_locations_xi_eta,'d-opt','mdeit3',dim);
g_d_opt_mdeit3_xi_eta  = @(q) g(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,...
    vector_to_sensor_locations_xi_eta,jac_coor_transf_xi_eta,'d-opt','mdeit3',dim);
 
% 3D-mdeit for (xi,eta) - a opt
f_a_opt_mdeit3_xi_eta = @(q) f(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,vector_to_sensor_locations_xi_eta,'a-opt','mdeit3',dim);
g_a_opt_mdeit3_xi_eta  = @(q) g(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,...
    vector_to_sensor_locations_xi_eta,jac_coor_transf_xi_eta,'a-opt','mdeit3',dim);


% 3D-mdeit for cylindrical coordinates
f_d_opt_mdeit3_cyl = @(q) f(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,vector_to_sensor_locations_con,'d-opt','mdeit3',dim);
g_d_opt_mdeit3_cyl  = @(q) g(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,...
    vector_to_sensor_locations_con,jac_coor_transf_cylindrical,'d-opt','mdeit3',dim);

%% PlotFcns
% function stop = outfun_d_opt(x,optimValues,state)
%
% switch state
%     case 'iter'
%         % Make updates to plot or guis as needed
%         this_axis = gca;
%         sensor_locations_k = vector_to_sensor_locations(x);
%         img_k = assign_sensor_locations(imgi,sensor_locations_k);
%         plot_sensors(img_k,false,'r','s',this_axis);
%         axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
%         box on;grid on;
%         view(3)
%         drawnow
%     case 'interrupt'
%         % Probably no action here. Check conditions to see
%         % whether optimization should quit.
%     case 'init'
%         hold on
%         show_fem(imgi);
%         % camlight;lighting gouraud
%     case 'done'
%         % Cleanup of plots, guis, or final plot
%         this_axis = gca;
%         sensor_locations_k = vector_to_sensor_locations(x);
%         img_k = assign_sensor_locations(imgi,sensor_locations_k);
%         plot_sensors(img_k,false,'b','s',this_axis);
%         axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
%         box on;grid on;
%         view(3)
%         drawnow
%     otherwise
% end
%
% stop = false; %continue
% end
%
% function stop = outfun_a_opt(x,optimValues,state)
% switch state
%     case 'iter'
%
%         % Make updates to plot or guis as needed
%         this_axis = gca;
%         sensor_locations_k = vector_to_sensor_locations(x);
%         img_k = assign_sensor_locations(imgi,sensor_locations_k);
%         plot_sensors(img_k,false,'r','o',this_axis);
%         axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
%         box on;grid on;
%         view(3)
%         drawnow
%     case 'interrupt'
%         % Probably no action here. Check conditions to see
%         % whether optimization should quit.
%     case 'init'
%         hold on
%         show_fem(imgi);
%         % camlight; lighting gouraud
%
%     case 'done'
%         % Cleanup of plots, guis, or final plot
%         this_axis = gca;
%         sensor_locations_k = vector_to_sensor_locations(x);
%         img_k = assign_sensor_locations(imgi,sensor_locations_k);
%         plot_sensors(img_k,false,'b','x',this_axis);
%         axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
%         box on;grid on;
%         view(3)
%         drawnow
%     otherwise
% end
%
% stop = false; %continue
% end

%% OutFcns
function [xsol,fval,history] = runfmincon(funcgrad,q0,lb,ub,options)
% Set up shared variables with outfun
history.x = [];
history.fval = [];

options.OutputFcn = @outfun;

    function stop = outfun(x,optimValues,state)
        stop = false;

        switch state
            case 'init'
                history.fval = [history.fval, optimValues.fval];
                history.x = [history.x, x];
            case 'iter'
                history.fval = [history.fval, optimValues.fval];
                history.x = [history.x, x];
            case 'done'
            otherwise
        end
    end

[xsol,fval] = fmincon(funcgrad,q0,[],[],[],[],lb,ub,[],options);

end


function [xsol,fval,history] = runfminunc(funcgrad,q0,options)
% Set up shared variables with outfun
history.x = [];
history.fval = [];

options.OutputFcn = @outfun;

    function stop = outfun(x,optimValues,state)
        stop = false;

        switch state
            case 'init'
                history.fval = [history.fval, optimValues.fval];
                history.x = [history.x, x];
            case 'iter'
                history.fval = [history.fval, optimValues.fval];
                history.x = [history.x, x];
            case 'done'
            otherwise
        end
    end

[xsol,fval] = fminunc(funcgrad,q0,options);

end
%% Define function+gradient function

function [func,grad] = funcwithgrad(q,f_impl,g_impl)
% Calculate objective f
func = f_impl(q);

if nargout > 1 % gradient required
    grad =  g_impl(q);
end
end

%% Define homotopy

imgtempi = assign_sensor_locations(imgh_recon,all_sensor_locations_0(:,:,1));
[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtempi,A);

% THIS DOES NOT WORK WELL, the cost function change is smaller than without
% this
% f0 = f_d_opt_mdeit3_cyl(sensor_locations_to_vector_con(all_sensor_locations_0(:,:,1)));
% surrogate0 = -norm([Jx;Jy;Jz],'fro')^2;

function phi0 = compute_surrogate_objective(img,A,q,vector_to_sensor_locations)
sensor_locations = vector_to_sensor_locations(q);

img = assign_sensor_locations(img,sensor_locations);

[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);
J_3axis = [Jx;Jy;Jz];

phi0 = -norm(J_3axis,'fro')^2;
end

function [phi0,dphi] = compute_obj_gradient_surrogate(img,A,q,jac_coor_transf,vector_to_sensor_locations)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

sensor_locations = vector_to_sensor_locations(q);

img = assign_sensor_locations(img,sensor_locations);

[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);

J = [Jx;Jy;Jz];

phi0 = -norm(J,'fro')^2;

% dJxds2 = compute_dJ_xyz(img,1);
dJ = compute_dJxyz_xyz_optimized(img);

dJxds = cell(1,3);
dJyds = cell(1,3);
dJzds = cell(1,3);

for p = 1:3
    dJxds{p} = dJ{1,p};
    dJyds{p} = dJ{2,p};
    dJzds{p} = dJ{3,p};
end

W = J;

gphi0 = zeros(3,n_sensors);

block_size = n_sensors * n_stim;

for p = 1:3
    for m = 1:n_sensors
        % indices inside one block (x, y, or z)
        ids_local = m : n_sensors : block_size;

        % --- X component ---
        ids_x = ids_local;
        dJx_m = reshape(dJxds{p}(m,:,:), [n_stim, n_elem]);
        Wx_m  = W(ids_x,:);
        term_x = sum(Wx_m(:) .* dJx_m(:));

        % --- Y component ---
        ids_y = ids_local + block_size;
        dJy_m = reshape(dJyds{p}(m,:,:), [n_stim, n_elem]);
        Wy_m  = W(ids_y,:);
        term_y = sum(Wy_m(:) .* dJy_m(:));

        % --- Z component ---
        ids_z = ids_local + 2*block_size;
        dJz_m = reshape(dJzds{p}(m,:,:), [n_stim, n_elem]);
        Wz_m  = W(ids_z,:);
        term_z = sum(Wz_m(:) .* dJz_m(:));

        gphi0(p,m) = -2 * (term_x + term_y + term_z);
    end
end

% gphi0 = reshape(gphi0.', 1, []);

% Convert cartesian derivatives to other coordinate derivatives (
% generalized)
jacobian_coordinate_transformation = jac_coor_transf(sensor_locations);
dphidq = dphidp_to_dphidq(sensor_locations,gphi0,jacobian_coordinate_transformation);

dphi = reshape(dphidq.', 1, []);
end

surrogate_obj = @(q) compute_surrogate_objective(imgh_recon,A,q,vector_to_sensor_locations_con);
surrogate_obj_grad = @(q) compute_obj_gradient_surrogate(imgh_recon,A,q,jac_coor_transf_cylindrical,vector_to_sensor_locations_con);

surrogate_obj_xi = @(q) compute_surrogate_objective(imgh_recon,A,q,vector_to_sensor_locations_xi);
surrogate_obj_grad_xi = @(q) compute_obj_gradient_surrogate(imgh_recon,A,q,jac_coor_transf_xi,vector_to_sensor_locations_xi);

surrogate_obj_xi_eta = @(q) compute_surrogate_objective(imgh_recon,A,q,vector_to_sensor_locations_xi_eta);
surrogate_obj_grad_xi_eta = @(q) compute_obj_gradient_surrogate(imgh_recon,A,q,jac_coor_transf_xi,vector_to_sensor_locations_xi_eta);

function phi = homotopy_func(q,alpha,obj,surrogate_obj)

if abs(alpha-0.0)<1e-9
    phi0 = surrogate_obj(q);   
    phi = phi0;
    return
elseif abs(alpha-1.0)<1e-9
    phi = obj(q);
    return
end

phi0 = surrogate_obj(q);
phi = alpha*obj(q) + (1-alpha)*phi0;
end

function dphi = grad_homotopy_func(q,alpha,grad,surrogate_obj_grad)

if abs(alpha-0.0)<1e-9
    [phi0,gphi0] = surrogate_obj_grad(q);
    dphi = gphi0;
    return
elseif abs(alpha-1.0)<1e-9
    dphi = grad(q);
    % disp(dphi);
    return
end

[phi0,gphi0] = surrogate_obj_grad(q);
dphi = alpha*grad(q) + (1-alpha)*gphi0;

return

end





%% TEST IF GRADIENT COMPUTATION IS THE SAME AS FD result (se central differences)
% jac_coor_transf_test = @(s) compute_jacobian_coordinate_transformation_xi_eta(s,rmin,rmax,zmin,zmax);
% 
% f_a_opt_mdeit3_test = @(q) f(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,vector_to_sensor_locations_xi_eta,'a-opt','mdeit3',dim);
% g_a_opt_mdeit3_test  = @(q) g(imgh_recon,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,...
%     vector_to_sensor_locations_xi_eta,jac_coor_transf_test,'a-opt','mdeit3',dim);
% 
% 
% % CHECK if gradient is correct
% q0 = sensor_locations_to_vector_xi_eta(sensor_locations_circle_height);
% 
% g_0 = g_a_opt_mdeit3_test (q0);
% f_impl = f_a_opt_mdeit3_test;
% 
% deltas = [1e-3 1e-5 1e-9];
% 
% grad_fd = zeros(numel(deltas),numel(q0));
% 
% for i = 1:numel(deltas)
%     delta = 1e-7;
% 
%     for k = 1:length(q0)
%         fprintf('Iteration %i/%i\n',k,numel(q0));
%         qplus = q0;
%         qminus = q0;
%         qplus(k) = qplus(k) + delta;
%         qminus(k) = qminus(k) - delta;
% 
%         grad_fd(i,k) = (f_impl(qplus)-f_impl(qminus))/(2*delta);
%     end
% end
% 
% figure;
% hold on
% markers = {'.','o','s'};
% for i = 1:numel(deltas)
%     plot(grad_fd(i,:),markers{i})
% end
% plot(g_0)
% hold off
% 
% msg = {};
% for i = 1:numel(deltas)
%     msg{i} = sprintf('FD: d = %2.2g',deltas(i));
% end
% msg{i+1} = 'Analytic';
% legend(msg)

% disp('Finite-Differences')
% disp(grad_fd(1:5));
% disp('Analytical')
% disp(g_0(1:5))



%% Optimize 3-axis A-optimality with fmincon in (xi,eta)-coordinate of sensors - fminunc
fprintf(' (xi,eta) 3-axis A-optimality OED - fminunc\n')

options = optimoptions('fminunc',...
    'Display','iter','MaxIterations',max_iterations,...
    'StepTolerance',1e-5,'OptimalityTolerance',1e-5,...
    'Algorithm','quasi-newton','HessianApproximation','lbfgs',...
    'SpecifyObjectiveGradient',true,'UseParallel',false);

fun_grad_fminunc_a_opt_mdeit_3 = @(q) funcwithgrad(q,f_a_opt_mdeit3_xi_eta,g_a_opt_mdeit3_xi_eta);

q0_unc = sensor_locations_to_vector_xi_eta(sensor_locations_circle_height);

%Sanity check
s = vector_to_sensor_locations_xi_eta(q0_unc);
assert(all(abs(s(:)-sensor_locations_circle_height(:))<1e-10));
sensor_locations_n = sensor_locations_circle_height;

% A PROBLEM HAPPENS WHEN WE PICK THIS INITIAL CONDITION. THE OPTIMIZATION
% DOES OT EVEN ITERATE!!!!!
% q0_unc = sensor_locations_to_vector_xi_eta(sensor_locations_circle);
% sensor_locations_n = sensor_locations_circle;


alphas = [1.0];

all_x_a_opt_con = zeros(numel(alphas),length(q0_unc));

history = cell(1,length(alphas));

figure;
show_fem(imgh_recon);
imgtempi = assign_sensor_locations(imgh_recon,sensor_locations_circle_height);
plot_sensors(imgtempi,[],'r','s');
drawnow
axis([-2 2 -2 2 -0.1 3.1]);

for n = 1:numel(alphas)
    % Initial conditions
    q0_unc = sensor_locations_to_vector_xi_eta(sensor_locations_n);

    cla;
    hold on
    show_fem(imgh_recon);
    imgtempi = assign_sensor_locations(imgh_recon,sensor_locations_circle_height);
    plot_sensors(imgtempi,[],'r','s');
    imgtempi = assign_sensor_locations(imgh_recon,sensor_locations_n);
    plot_sensors(imgtempi,[],'b','.');
    axis([-2 2 -2 2 -0.1 3.1]);
    drawnow
    hold off

    func_homotopy = @(q) homotopy_func(q,alphas(n),fun_grad_fminunc_a_opt_mdeit_3,surrogate_obj_xi_eta);
    grad_homotopy = @(q) grad_homotopy_func(q,alphas(n),g_a_opt_mdeit3_xi_eta,surrogate_obj_grad_xi_eta);

    func_grad = @(q) funcwithgrad(q,func_homotopy,grad_homotopy);

    [q_a_opt_unc,fval,history{n}] = runfminunc(func_grad,q0_unc,options);

    all_x_a_opt_con(n,:) = q_a_opt_unc;
    img_a_opt_unc{n} = assign_sensor_locations(imgh_recon,vector_to_sensor_locations_xi_eta(q_a_opt_unc));

    sensor_locations_n = vector_to_sensor_locations_xi_eta(q_a_opt_unc);
end

%Change in cost function:

f0 = fun_grad_fminunc_a_opt_mdeit_3(sensor_locations_to_vector_xi_eta(sensor_locations_circle_height));
fend = fun_grad_fminunc_a_opt_mdeit_3(q_a_opt_unc);

fprintf('Change in objective : %2.2f %%\n',(fend-f0)/f0*100)


%% Optimize 3-axis D-optimality with fmincon in (xi,eta)-coordinate of sensors - fminunc
% fprintf(' (xi,eta) 3-axis D-optimality OED - fminunc\n')
% 
% options = optimoptions('fminunc',...
%     'Display','iter','MaxIterations',max_iteratons,...
%     'StepTolerance',1e-5,'OptimalityTolerance',1e-5,...
%     'Algorithm','quasi-newton','HessianApproximation','lbfgs',...
%     'SpecifyObjectiveGradient',true,'UseParallel',false);
% 
% fun_grad_fminunc_d_opt_mdeit_3 = @(q) funcwithgrad(q,f_d_opt_mdeit3_xi_eta,g_d_opt_mdeit3_xi_eta);
% 
% q0_unc = sensor_locations_to_vector_xi_eta(sensor_locations_circle_height);
% 
% %Sanity check
% s = vector_to_sensor_locations_xi_eta(q0_unc);
% assert(all(abs(s(:)-sensor_locations_circle_height(:))<1e-10));
% 
% sensor_locations_n = sensor_locations_circle_height;
% alphas = [1.0];
% 
% all_x_a_opt_con = zeros(numel(alphas),length(q0_unc));
% 
% history = cell(1,length(alphas));
% 
% figure;
% show_fem(imgh_recon);
% imgtempi = assign_sensor_locations(imgh_recon,sensor_locations_circle_height);
% plot_sensors(imgtempi,[],'r','s');
% drawnow
% axis([-2 2 -2 2 -0.1 3.1]);
% 
% for n = 1:numel(alphas)
%     % Initial conditions
%     q0_unc = sensor_locations_to_vector_xi_eta(sensor_locations_n);
% 
%     cla;
%     hold on
%     show_fem(imgh_recon);
%     imgtempi = assign_sensor_locations(imgh_recon,sensor_locations_circle_height);
%     plot_sensors(imgtempi,[],'r','s');
%     imgtempi = assign_sensor_locations(imgh_recon,sensor_locations_n);
%     plot_sensors(imgtempi,[],'b','.');
%     axis([-2 2 -2 2 -0.1 3.1]);
%     drawnow
%     hold off
% 
%     func_homotopy = @(q) homotopy_func(q,alphas(n),fun_grad_fminunc_d_opt_mdeit_3,surrogate_obj_xi_eta);
%     grad_homotopy = @(q) grad_homotopy_func(q,alphas(n),g_d_opt_mdeit3_xi_eta,surrogate_obj_grad_xi_eta);
% 
%     func_grad = @(q) funcwithgrad(q,func_homotopy,grad_homotopy);
% 
%     [q_a_opt_unc,fval,history{n}] = runfminunc(func_grad,q0_unc,options);
% 
%     all_x_a_opt_con(n,:) = q_a_opt_unc;
%     img_a_opt_unc{n} = assign_sensor_locations(imgh_recon,vector_to_sensor_locations_xi_eta(q_a_opt_unc));
% 
%     sensor_locations_n = vector_to_sensor_locations_xi_eta(q_a_opt_unc);
% end
% 
% %Change in cost function:
% 
% f0 = fun_grad_fminunc_d_opt_mdeit_3(sensor_locations_to_vector_xi_eta(sensor_locations_circle_height));
% fend = fun_grad_fminunc_d_opt_mdeit_3(q_a_opt_unc);
% 
% fprintf('Change in objective : %2.2f %%\n',(fend-f0)/f0*100)
% 
% % Plot results figure
% 
% [Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img_a_opt_unc{1},A);
% J_3axis_opt = [Jx;Jy;Jz];

%% Simulate forward data

fprintf('Solving fmdl for initial sensors config\n');
% Reconstruct with old sensor postions
imgtempi = imgi;
imgtempi = assign_sensor_locations(imgtempi,sensor_locations_circle_height);
imgtemph = imgh;
imgtemph = assign_sensor_locations(imgtemph,sensor_locations_circle_height);

yi = fwd_solve_mdeit_local(imgtempi);
yi = [yi.Bx(:);yi.By(:);yi.Bz(:)];
yh = fwd_solve_mdeit_local(imgtemph);
yh = [yh.Bx(:);yh.By(:);yh.Bz(:)];
dy_old = add_measurement_noise_difference(yi,yh,snr);

fprintf('Solving fmdl for optimal sensor config\n');

% Reconstruct with new sensor postions
imgtempi = imgi;
imgtempi = assign_sensor_locations(imgtempi,sensor_locations_n);
imgtemph = imgh;
imgtemph = assign_sensor_locations(imgtemph,sensor_locations_n);

yi = fwd_solve_mdeit_local(imgtempi);
yi = [yi.Bx(:);yi.By(:);yi.Bz(:)];
yh = fwd_solve_mdeit_local(imgtemph);
yh = [yh.Bx(:);yh.By(:);yh.Bz(:)];
dy = add_measurement_noise_difference(yi,yh,snr);

%% Reconstruct posterior means pre and post sensor optimization

% profile on 

fprintf('Reconstructing posterior means pre and post sensor optimization\n');

% We need only the diagonals of the posterior covariance matrix, so avoid
% computing inv(...), it's very expensive. Its faster than the alternative
% that chat gpt suggested tho...

% Is this form faster computationally????????????
% Check!!!!!!!!!!!!!!!!!!!!!!!! (NOT SURE!!!!!!!! MAYBE CHANGE BACK TO Gamma_post_no_opt = inv(J_3axis_no_opt.'*inv_Gamma_noise_3_axis*J_3axis_no_opt + inv_Gamma_prior_3_axis);
% Γpost = Γ_pr −Γ_pr*J^T*(J*Γ_pr*J^T +Γ_noise)^{−1}J*Γ_pr.

imgtemph = imgh_recon;
imgtemph = assign_sensor_locations(imgtemph,sensor_locations_circle_height);

[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtemph,A);
J_3axis_no_opt = [Jx;Jy;Jz];

Gamma_post_no_opt = Gamma_prior-Gamma_prior*J_3axis_no_opt.'*inv(J_3axis_no_opt*Gamma_prior*J_3axis_no_opt.'+Gamma_noise_3_axis)*J_3axis_no_opt*Gamma_prior;
% Gamma_post_no_opt = inv(J_3axis_no_opt.'*inv_Gamma_noise_3_axis*J_3axis_no_opt + inv_Gamma_prior_3_axis);

imgtemph = imgh_recon;
imgtemph = assign_sensor_locations(imgtemph,sensor_locations_n);

[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtemph,A);
J_3axis_opt = [Jx;Jy;Jz];

% Gamma_post_opt = inv(J_3axis_opt.'*inv_Gamma_noise_3_axis*J_3axis_opt + inv_Gamma_prior_3_axis);
Gamma_post_opt = Gamma_prior-Gamma_prior*J_3axis_opt.'*inv(J_3axis_opt*Gamma_prior*J_3axis_opt.'+Gamma_noise_3_axis)*J_3axis_opt*Gamma_prior;

% Prior contribution is zero because mean of noise is 0. Avoid it so
% computation is faster
posterior_mean_no_opt = Gamma_post_no_opt*(...
    J_3axis_no_opt.'*inv_Gamma_noise_3_axis*(dy_old-zeros(3*n_stim*n_sensors,1)));

posterior_mean_opt = Gamma_post_opt*(...
    J_3axis_opt.'*inv_Gamma_noise_3_axis*(dy-zeros(3*n_stim*n_sensors,1)));

img_temp = assign_sensor_locations(imgh_recon,sensor_locations_circle_height);
img_temp.calc_colours.ref_level = 0;
img_temp.calc_colours.greylev = -0.05;
img_temp.calc_colours.cb_shrink_move = [0.3,0.6,0.02];


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

% Reconstruct with generalized cross validation and TSVD

fprintf('Doing reconstruction with generalized cross-validation\n');
lambda_vector = logspace(-20,0,20);
[lambda_opt,sigma_gcv] = generalized_cross_validation(J_3axis_opt,dy,lambda_vector);

% profile viewer
%%
figure
subplot(3,2,1)
img_temp.elem_data = diag(Gamma_prior);
show_fem(img_temp);
eidors_colourbar(img_temp);
plot_sensors(img_temp)
title('Prior point-variance','Interpreter','latex')

subplot(3,2,2)
img_temp.elem_data = diag(Gamma_post_opt);
show_fem(img_temp);
eidors_colourbar(img_temp);
plot_sensors(img_a_opt_unc{1})
title('Posteriror point-variance','Interpreter','latex')

subplot(3,2,3)
img_temp.elem_data = posterior_mean_no_opt;
show_fem(img_temp);
eidors_colourbar(img_temp);
plot_sensors(img_temp)
title('Posterior mean - no opt','Interpreter','latex')

subplot(3,2,4)
img_temp.elem_data = posterior_mean_opt;
show_fem(img_temp);
eidors_colourbar(img_temp);
plot_sensors(img_a_opt_unc{1})
title('Posterior mean - opt','Interpreter','latex')

subplot(3,2,5)
img_temp.elem_data = sigma_gcv;
show_fem(img_temp);
eidors_colourbar(img_temp);
plot_sensors(img_a_opt_unc{1})
title('GCV Tykhonov Prior - opt','Interpreter','latex')

subplot(3,2,6)
show_fem(imgi);
eidors_colourbar(imgi);
title('ground-truth','Interpreter','latex')




%% Optimize 3-axis D-optimality with fmincon in xi-coordinate of sensors - fminunc
fprintf(' (xi) 3-axis D-optimality OED - fminunc\n')

options = optimoptions('fminunc',...
    'Display','iter','MaxIterations',max_iterations,...
    'StepTolerance',1e-5,'OptimalityTolerance',1e-5,...
    'Algorithm','quasi-newton','HessianApproximation','lbfgs',...
    'SpecifyObjectiveGradient',true,'UseParallel',false);

fun_grad_fminunc_d_opt_mdeit_3 = @(q) funcwithgrad(q,f_d_opt_mdeit3_xi,g_d_opt_mdeit3_xi);

q0_unc = sensor_locations_to_vector_xi(sensor_locations_circle);

%Sanity check
s = vector_to_sensor_locations_xi(q0_unc);
assert(all(abs(s(:)-sensor_locations_circle(:))<1e-10));

sensor_locations_n = sensor_locations_circle;
alphas = [1.0];

all_x_a_opt_con = zeros(numel(alphas),length(q0_unc));


history = cell(1,length(alphas));

figure;
show_fem(imgh_recon);
imgtempi = assign_sensor_locations(imgh_recon,sensor_locations_circle);
plot_sensors(imgtempi,[],'r','s');
drawnow
axis([-2 2 -2 2 -0.1 3.1]);

for n = 1:numel(alphas)
    % Initial conditions
    q0_unc = sensor_locations_to_vector_xi(sensor_locations_n);
    
    cla;
    hold on
    show_fem(imgh_recon);
    imgtempi = assign_sensor_locations(imgh_recon,sensor_locations_circle);
    plot_sensors(imgtempi,[],'r','s');
    imgtempi = assign_sensor_locations(imgh_recon,sensor_locations_n);
    plot_sensors(imgtempi,[],'b','.');
    axis([-2 2 -2 2 -0.1 3.1]);
    drawnow
    hold off
    
    func_homotopy = @(q) homotopy_func(q,alphas(n),fun_grad_fminunc_d_opt_mdeit_3,surrogate_obj_xi);
    grad_homotopy = @(q) grad_homotopy_func(q,alphas(n),g_d_opt_mdeit3_xi,surrogate_obj_grad_xi);

    func_grad = @(q) funcwithgrad(q,func_homotopy,grad_homotopy);

    [q_a_opt_unc,fval,history{n}] = runfminunc(func_grad,q0_unc,options);

    all_x_a_opt_con(n,:) = q_a_opt_unc;
    img_a_opt_unc{n} = assign_sensor_locations(imgh_recon,vector_to_sensor_locations_xi(q_a_opt_unc));

    sensor_locations_n = vector_to_sensor_locations_xi(q_a_opt_unc);
end

%Change in cost function:

f0 = fun_grad_fminunc_d_opt_mdeit_3(sensor_locations_to_vector_xi(sensor_locations_circle));
fend = fun_grad_fminunc_d_opt_mdeit_3(q_a_opt_unc);

fprintf('Change in objective : %2.2f %%\n',(fend-f0)/f0*100)



%% Optimize 3-axis D-optimality with fmincon in z-coordinate of sensors - fmincon
fprintf(' (Cylinder Region) 3-axis A-optimality OED - fmincon\n')

options = optimoptions('fmincon',...
    'Display','iter','MaxIterations',max_iterations,...
    'StepTolerance',1e-5,'OptimalityTolerance',1e-5,...
    'Algorithm','interior-point','HessianApproximation','lbfgs',...
    'SpecifyObjectiveGradient',true,'UseParallel',use_parallel);

fun_grad_fmincon_d_opt_mdeit_3 = @(q) funcwithgrad(q,f_d_opt_mdeit3_cyl,g_d_opt_mdeit3_cyl);

lb = [ rmin*ones(1,n_sensors),   -inf*ones(1,n_sensors),   0*ones(1,n_sensors)];
ub = [ rmax*ones(1,n_sensors),   +inf*ones(1,n_sensors),   model_parameters.height*ones(1,n_sensors)];

q0_con = sensor_locations_to_vector_con(all_sensor_locations_0(:,:,1));
all_x_a_opt_con = zeros(n_trial,length(q0_con));
history = cell(1,n_trial);

sensor_locations_n = all_sensor_locations_0(:,:,1);

alphas = [0 0.1 0.2 0.3 0.4 0.5 0.75 1];

figure;
show_fem(imgh_recon);
imgtempi = assign_sensor_locations(imgh_recon,all_sensor_locations_0(:,:,1));
plot_sensors(imgtempi,[],'r','s');
drawnow
axis([-2 2 -2 2 -0.1 3.1]);

for n = 1:numel(alphas)

    % Initial conditions
    q0_con = sensor_locations_to_vector_con(sensor_locations_n);
    
    cla;
    hold on
    show_fem(imgh_recon);
    imgtempi = assign_sensor_locations(imgh_recon,all_sensor_locations_0(:,:,1));
    plot_sensors(imgtempi,[],'r','s');
    imgtempi = assign_sensor_locations(imgh_recon,sensor_locations_n);
    plot_sensors(imgtempi,[],'b','.');
    axis([-2 2 -2 2 -0.1 3.1]);
    drawnow
    hold off
    
    func_homotopy = @(q) homotopy_func(q,alphas(n),fun_grad_fmincon_d_opt_mdeit_3,surrogate_obj);
    grad_homotopy = @(q) grad_homotopy_func(q,alphas(n),g_d_opt_mdeit3_cyl,surrogate_obj_grad);

    func_grad = @(q) funcwithgrad(q,func_homotopy,grad_homotopy);

    [q_a_opt_con,fval,history{n}] = runfmincon(func_grad,q0_con,lb,ub,options);

    all_x_a_opt_con(n,:) = q_a_opt_con;
    img_a_opt_con{n} = assign_sensor_locations(imgh_recon,vector_to_sensor_locations_con(q_a_opt_con));

    sensor_locations_n = vector_to_sensor_locations_con(q_a_opt_con);
end

%Change in cost function:

f0 = fun_grad_fmincon_d_opt_mdeit_3(sensor_locations_to_vector_con(all_sensor_locations_0(:,:,1)));
fend = fun_grad_fmincon_d_opt_mdeit_3(q_a_opt_con);

fprintf('Change in objective : %2.2f %%\n',(fend-f0)/f0*100)

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

function lambda = compute_lambda_xyz(img)

n_nodes = size(img.fwd_model.nodes,1);
num_sensors = numel(img.fwd_model.sensors);
n_elem = size(img.fwd_model.elems,1);

sigma = img.elem_data;

% Compute Gamma matrices
img = compute_gamma_matrices_local(img);
GammaX = img.Gamma1;
GammaY = img.Gamma2;
GammaZ = img.Gamma3;

% Solve the adjoint problem for each sensor to get lambda vectors
A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

GammaXT = GammaX.';
GammaYT = GammaY.';
GammaZT = GammaZ.';

lambdaX = zeros(n_nodes,num_sensors);
lambdaY = zeros(n_nodes,num_sensors);
lambdaZ = zeros(n_nodes,num_sensors);

% Solve the adjoint problem for each sensor to get lambda vectors
for m = 1:num_sensors

    [lambdaX(:,m),~,~] = pcg(A_matrix,-GammaXT(:,m),1e-9,n_elem,Mfun,Nfun);
    [lambdaY(:,m),~,~] = pcg(A_matrix,-GammaYT(:,m),1e-9,n_elem,Mfun,Nfun);
    [lambdaZ(:,m),~,~] = pcg(A_matrix,-GammaZT(:,m),1e-9,n_elem,Mfun,Nfun);
end

lambda{1} = lambdaX;
lambda{2} = lambdaY;
lambda{3} = lambdaZ;


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

function img = compute_gamma_matrices_local(img)

mu_factor = img.fwd_model.mu0/(4*pi);

num_sensors = numel(img.fwd_model.sensors);

sensor_locations = zeros(numel(img.fwd_model.sensors),3);

for i = 1: numel(img.fwd_model.sensors)
    sensor_locations(i,:) = img.fwd_model.sensors(i).position;
end

[Rx,Ry,Rz,fmdl] = compute_r_matrices_local(img.fwd_model,sensor_locations);

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

function [data,u] = fwd_solve_mdeit_local(img,u)

if img.fwd_solve.get_all_meas ~=1
    error('get_all_meas should be set to 1 for u.volt field to exist')
end

if nargin<2
    %% Forward solve the EIT model
    u = fwd_solve(img);
end

%% Output
img = compute_gamma_matrices_local_optimized(img);

Bx = img.Gamma1*u.volt;
By = img.Gamma2*u.volt;
Bz = img.Gamma3*u.volt;

data = struct('Bx',Bx,'By',By,'Bz',Bz);
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

function [Rx,Ry,Rz,fmdl] = compute_r_matrices_local(fmdl,sensor_locations)
numElements = size(fmdl.elems,1);
numSensors = size(sensor_locations,1);

Rx = zeros(numSensors,numElements);
Ry = zeros(numSensors,numElements);
Rz = zeros(numSensors,numElements);

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

for m = 1:numSensors
    rm = sensor_locations(m,:);
    fun = @(x,y,z) (rm - [x,y,z])./ (sum((rm - [x, y, z]).^2)^1.5);
    for k = 1:numElements

        %Find the vertices of the tetrahedron
        v = fmdl.nodes(fmdl.elems(k,:),:);

        J = [(v(2,:)-v(1,:))',(v(3,:)-v(1,:))',(v(4,:)-v(1,:))'];

        detJ = abs(det(J));

        R = 0;
        for i = 1:length(weights)

            r = coord(i,1);
            s = coord(i,2);
            t = coord(i,3);

            xi = v(1,:)' + J * [r; s; t];

            R = R + weights(i)*fun(xi(1),xi(2),xi(3));
        end

        R =  (detJ / 6) * R;
        Rx(m,k) = R(1);
        Ry(m,k) = R(2);
        Rz(m,k) = R(3);
    end
end

fmdl.R.Rx = Rx;
fmdl.R.Ry = Ry;
fmdl.R.Rz = Rz;

end

function dRdp = dRmkj_dp(fmdl,j,p)
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

numElements = size(fmdl.elems,1);
numSensors  = numel(fmdl.sensors);
numQuadPts  = length(weights);

% --- Precompute element vertices ---
% V: 4 x 3 x numElements
% Stack vertices: 4 x numElements x 3
V = reshape(fmdl.nodes(fmdl.elems',:), [4, numElements, 3]);

v1 = squeeze(V(1,:,:));  % numElements x 3
v2 = squeeze(V(2,:,:));
v3 = squeeze(V(3,:,:));
v4 = squeeze(V(4,:,:));

% Initialize J
J = zeros(3,3,numElements);  % 3 x 3 x numElements

% Fill J using transpose of row vectors
J(:,1,:) = permute(v2 - v1, [2 3 1]);  % 3 x 1 x numElements
J(:,2,:) = permute(v3 - v1, [2 3 1]);
J(:,3,:) = permute(v4 - v1, [2 3 1]);

% Determinants of all J
detJ = zeros(1,numElements);
for k = 1:numElements
    detJ(k) = abs(det(J(:,:,k)));
end

% --- Compute dRdp for all sensors ---
dRdp = zeros(numSensors,numElements);

for m = 1:numSensors
    rm = fmdl.sensors(m).position;  % 1 x 3
    dm = @(X) rm - X;                % vectorized over X: N x 3

    % Evaluate quadrature points in all elements
    % xi: numQuadPts x 3 x numElements
    xi = reshape(v1', [3,1,numElements]) + ...
        J(:,1,:).*reshape(coord(:,1),[1,numQuadPts,1]) + ...
        J(:,2,:).*reshape(coord(:,2),[1,numQuadPts,1]) + ...
        J(:,3,:).*reshape(coord(:,3),[1,numQuadPts,1]); % 3 x numQuadPts x numElements

    xi = permute(xi,[2 1 3]); % numQuadPts x 3 x numElements

    % dm_vec = rm - xi, same size
    dm_vec = rm - xi; % numQuadPts x 3 x numElements

    % Compute fun values at all quadrature points
    % ((j==p)*norm(dm)^2 - 3*dm_j*dm_p) / norm(dm)^5
    dm_norm2 = sum(dm_vec.^2,2);           % numQuadPts x 1 x numElements
    dm_j    = dm_vec(:,j,:);               % numQuadPts x 1 x numElements
    dm_p    = dm_vec(:,p,:);

    funvals = ((j==p)*dm_norm2 - 3*dm_j.*dm_p) ./ (dm_norm2.^(5/2)); % numQuadPts x 1 x numElements

    % Integrate over quadrature points using weights
    dRdp(m,:) = squeeze(sum(funvals .* reshape(weights, [numQuadPts,1,1]),1))' .* (detJ/6);
end


end

function dR = dRmkj_xyz(fmdl,j)
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

numElements = size(fmdl.elems,1);
numSensors  = numel(fmdl.sensors);
numQuadPts  = length(weights);

% --- Precompute element vertices ---
% V: 4 x 3 x numElements
% Stack vertices: 4 x numElements x 3
V = reshape(fmdl.nodes(fmdl.elems',:), [4, numElements, 3]);

v1 = squeeze(V(1,:,:));  % numElements x 3
v2 = squeeze(V(2,:,:));
v3 = squeeze(V(3,:,:));
v4 = squeeze(V(4,:,:));

% Initialize J
J = zeros(3,3,numElements);  % 3 x 3 x numElements

% Fill J using transpose of row vectors
J(:,1,:) = permute(v2 - v1, [2 3 1]);  % 3 x 1 x numElements
J(:,2,:) = permute(v3 - v1, [2 3 1]);
J(:,3,:) = permute(v4 - v1, [2 3 1]);

% Determinants of all J
detJ = zeros(1,numElements);
for k = 1:numElements
    detJ(k) = abs(det(J(:,:,k)));
end

% --- Compute dRdp for all sensors, and for all p ---
dRdx = zeros(numSensors,numElements);
dRdy = zeros(numSensors,numElements);
dRdz = zeros(numSensors,numElements);



for m = 1:numSensors
    rm = fmdl.sensors(m).position;  % 1 x 3

    % Evaluate quadrature points in all elements
    % xi: numQuadPts x 3 x numElements
    xi = reshape(v1', [3,1,numElements]) + ...
        J(:,1,:).*reshape(coord(:,1),[1,numQuadPts,1]) + ...
        J(:,2,:).*reshape(coord(:,2),[1,numQuadPts,1]) + ...
        J(:,3,:).*reshape(coord(:,3),[1,numQuadPts,1]); % 3 x numQuadPts x numElements

    xi = permute(xi,[2 1 3]); % numQuadPts x 3 x numElements

    % dm_vec = rm - xi, same size
    dm_vec = rm - xi; % numQuadPts x 3 x numElements

    dm_norm2 = sum(dm_vec.^2,2);           % numQuadPts x 1 x numElements
    dm_j    = dm_vec(:,j,:);               % numQuadPts x 1 x numElements

    for p = 1:3
        % Compute fun values at all quadrature points
        % ((j==p)*norm(dm)^2 - 3*dm_j*dm_p) / norm(dm)^5

        dm_p    = dm_vec(:,p,:);

        funvals = ((j==p)*dm_norm2 - 3*dm_j.*dm_p) ./ (dm_norm2.^(5/2)); % numQuadPts x 1 x numElements

        % Integrate over quadrature points using weights
        switch p
            case 1
                dRdx(m,:) = squeeze(sum(funvals .* reshape(weights, [numQuadPts,1,1]),1))' .* (detJ/6);
            case 2
                dRdy(m,:) = squeeze(sum(funvals .* reshape(weights, [numQuadPts,1,1]),1))' .* (detJ/6);
            case 3
                dRdz(m,:) = squeeze(sum(funvals .* reshape(weights, [numQuadPts,1,1]),1))' .* (detJ/6);
        end
    end
end

dR = {dRdx,dRdy,dRdz};

end

function dR = dRmkj_xyz_optimized(fmdl,j)

persistent ref_pts weights nq

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

numElements = size(fmdl.elems,1);
numSensors  = numel(fmdl.sensors);

nodes = fmdl.nodes;
elems = fmdl.elems;

% ---------------------------------------------------------
% Precompute element geometry (ONCE)
% ---------------------------------------------------------
v1  = zeros(3,numElements);
J   = zeros(3,3,numElements);

for k = 1:numElements
    verts = nodes(elems(k,:),:);

    v1(:,k) = verts(1,:).';

    J(:,:,k) = [ ...
        (verts(2,:) - verts(1,:))' , ...
        (verts(3,:) - verts(1,:))' , ...
        (verts(4,:) - verts(1,:))' ];
end

% Vectorized determinant (much faster than loop)
detJ = abs( ...
      J(1,1,:).*(J(2,2,:).*J(3,3,:) - J(2,3,:).*J(3,2,:)) ...
    - J(1,2,:).*(J(2,1,:).*J(3,3,:) - J(2,3,:).*J(3,1,:)) ...
    + J(1,3,:).*(J(2,1,:).*J(3,2,:) - J(2,2,:).*J(3,1,:)) );

detJ = reshape(detJ,numElements,1);

% ---------------------------------------------------------
% Precompute physical quadrature points for ALL elements
% Xi: 3 × nq × numElements
% ---------------------------------------------------------

Xi = zeros(3,nq,numElements);

for k = 1:numElements
    Xi(:,:,k) = v1(:,k) + J(:,:,k)*ref_pts;
end


% ---------------------------------------------------------
% Allocate output
% ---------------------------------------------------------

dRdx = zeros(numSensors,numElements);
dRdy = zeros(numSensors,numElements);
dRdz = zeros(numSensors,numElements);

w = reshape(weights,1,nq,1);   % 1 × nq × 1

sensors = fmdl.sensors;

for m = 1:numSensors

    rm = sensors(m).position';  % 1 x 3

    diff = rm - Xi;              % 3 × nq × numElements

    r2 = sum(diff.^2,1);         % 1 × nq × numElements
    r5 = r2.^(5/2);

    dm_j = diff(j,:,:);

    for p = 1:3

        dm_p = diff(p,:,:);

        num = (j==p)*r2 - 3*dm_j.*dm_p;

        integrand = num ./ r5;

        val = squeeze(sum(integrand .* w,2));  % 1 × numElements

        switch p
            case 1
                dRdx(m,:) = val .* (detJ/6);
            case 2
                dRdy(m,:) = val .* (detJ/6);
            case 3
                dRdz(m,:) = val .* (detJ/6);
        end
    end
end

dR = {dRdx,dRdy,dRdz};

end

function [dlambda,dR1,dR2] = compute_dlambda_ds(img,dim,p)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);

dlambda = zeros(n_nodes,num_sensors);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

sigma = img.elem_data;
A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

switch dim
    case 1
        dR1 = dRmkj_dp(img.fwd_model, 3, p);
        G1 = G.Gy;

        dR2 = dRmkj_dp(img.fwd_model, 2, p);
        G2 = G.Gz;
    case 2
        dR1 = dRmkj_dp(img.fwd_model, 1, p);
        dR2 = dRmkj_dp(img.fwd_model, 3, p);

        G1 = G.Gz;
        G2 = G.Gx;
    case 3
        dR1 = dRmkj_dp(img.fwd_model, 2, p);
        dR2 = dRmkj_dp(img.fwd_model, 1, p);

        G1 = G.Gx;
        G2 = G.Gy;
    otherwise
        error('Invalid dimension');
end


% Multiply elementwise by sigma (diagonal) without creating sparse diagonal matrix
% rhs = mu_factor*( dR1 * Sigma * G1 -  dR2* Sigma * G2);
rhs = mu_factor*( (dR1 .* sigma.')*G1 - (dR2 .* sigma.')*G2 );


parfor m = 1:num_sensors
    [dlambda(:,m),~,~] = pcg(A_matrix,rhs(m,:)',1e-9,numel(sigma),Mfun,Nfun);
end

end

function [dlambda,dR1,dR2] = compute_dlambda_xyz(img,dim)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);

dlambdaX = zeros(n_nodes,num_sensors);
dlambdaY = zeros(n_nodes,num_sensors);
dlambdaZ = zeros(n_nodes,num_sensors);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

switch dim
    case 1
        dR1 = dRmkj_xyz(img.fwd_model, 3);
        dR2 = dRmkj_xyz(img.fwd_model, 2);

        G1 = G.Gy;
        G2 = G.Gz;
    case 2
        dR1  = dRmkj_xyz(img.fwd_model, 1);
        dR2 = dRmkj_xyz(img.fwd_model, 3);

        G1 = G.Gz;
        G2 = G.Gx;
    case 3
        dR1 = dRmkj_xyz(img.fwd_model, 2);
        dR2 = dRmkj_xyz(img.fwd_model, 1);

        G1 = G.Gx;
        G2 = G.Gy;
    otherwise
        error('Invalid dimension');
end

sigma = img.elem_data;
A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

% Multiply elementwise by sigma (diagonal) without creating sparse diagonal matrix
% rhs = mu_factor*( dR1 * Sigma * G1 -  dR2* Sigma * G2);
rhs1 = mu_factor*( (dR1{1} .* sigma.')*G1 - (dR2{1} .* sigma.')*G2 );
rhs2 = mu_factor*( (dR1{2} .* sigma.')*G1 - (dR2{2} .* sigma.')*G2 );
rhs3 = mu_factor*( (dR1{3} .* sigma.')*G1 - (dR2{3} .* sigma.')*G2 );

parfor m = 1:num_sensors

    [dlambdaX(:,m),~,~] = pcg(A_matrix,rhs1(m,:)',1e-9,numel(sigma),Mfun,Nfun);
    [dlambdaY(:,m),~,~] = pcg(A_matrix,rhs2(m,:)',1e-9,numel(sigma),Mfun,Nfun);
    [dlambdaZ(:,m),~,~] = pcg(A_matrix,rhs3(m,:)',1e-9,numel(sigma),Mfun,Nfun);
end

dlambda = {dlambdaX,dlambdaY,dlambdaZ};

end




function [dlambda,dR1t,dR2t] = compute_dlambdaxyz_xyz_optimized(img)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);

dlambda = cell(3,3);

dlambda_dim_dx = zeros(n_nodes,num_sensors);
dlambda_dim_dy = zeros(n_nodes,num_sensors);
dlambda_dim_dz = zeros(n_nodes,num_sensors);

dR1t = cell(3,3);
dR2t = cell(3,3);


mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

sigma = img.elem_data;
A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

dR_all{1} = dRmkj_xyz_optimized(img.fwd_model,1);
dR_all{2} = dRmkj_xyz_optimized(img.fwd_model,2);
dR_all{3} = dRmkj_xyz_optimized(img.fwd_model,3);

sigmaT = sigma.';

n_elem = numel(sigma);

for dim = 1:3

    switch dim
        case 1
            dR1 = dR_all{3};
            dR2 = dR_all{2};

            G1 = G.Gy;
            G2 = G.Gz;
        case 2
            dR1 = dR_all{1};
            dR2 = dR_all{3};

            G1 = G.Gz;
            G2 = G.Gx;
        case 3
            dR1 = dR_all{2};
            dR2 = dR_all{1};

            G1 = G.Gx;
            G2 = G.Gy;
    end

    % Multiply elementwise by sigma (diagonal) without creating sparse diagonal matrix
    % rhs = mu_factor*( dR1 * Sigma * G1 -  dR2* Sigma * G2);

    rhsx = mu_factor*( (dR1{1} .* sigmaT)*G1 - (dR2{1} .* sigmaT)*G2 ); %for p=1
    rhsy = mu_factor*( (dR1{2} .* sigmaT)*G1 - (dR2{2} .* sigmaT)*G2 ); %for p=2
    rhsz = mu_factor*( (dR1{3} .* sigmaT)*G1 - (dR2{3} .* sigmaT)*G2 ); %for p=3

    % Solve the adjoint problem for each sensor to get lambda vectors
    parfor m = 1:num_sensors

        [dlambda_dim_dx(:,m),~,~] = pcg(A_matrix,rhsx(m,:)',1e-9,n_elem,Mfun,Nfun);
        [dlambda_dim_dy(:,m),~,~] = pcg(A_matrix,rhsy(m,:)',1e-9,n_elem,Mfun,Nfun);
        [dlambda_dim_dz(:,m),~,~] = pcg(A_matrix,rhsz(m,:)',1e-9,n_elem,Mfun,Nfun);
    end

    dlambda{dim,1} = dlambda_dim_dx;
    dlambda{dim,2} = dlambda_dim_dy;
    dlambda{dim,3} = dlambda_dim_dz;

    dR1t{dim,1} = dR1{1};
    dR1t{dim,2} = dR1{2};
    dR1t{dim,3} = dR1{3};

    dR2t{dim,1} = dR2{1};
    dR2t{dim,2} = dR2{2};
    dR2t{dim,3} = dR2{3};
end

end

function [dlambda,dR1t,dR2t] = compute_dlambdaxyz_xyz(img)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);

dlambda = cell(3,3);

dlambdaX = zeros(n_nodes,num_sensors);
dlambdaY = zeros(n_nodes,num_sensors);
dlambdaZ = zeros(n_nodes,num_sensors);

dR1t = cell(3,3);
dR2t = cell(3,3);


mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

sigma = img.elem_data;
A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

for dim = 1:3

    switch dim
        case 1
            dR1 = dRmkj_xyz(img.fwd_model, 3);
            dR2 = dRmkj_xyz(img.fwd_model, 2);

            G1 = G.Gy;
            G2 = G.Gz;
        case 2
            dR1  = dRmkj_xyz(img.fwd_model, 1);
            dR2 = dRmkj_xyz(img.fwd_model, 3);

            G1 = G.Gz;
            G2 = G.Gx;
        case 3
            dR1 = dRmkj_xyz(img.fwd_model, 2);
            dR2 = dRmkj_xyz(img.fwd_model, 1);

            G1 = G.Gx;
            G2 = G.Gy;
    end

    % Multiply elementwise by sigma (diagonal) without creating sparse diagonal matrix
    % rhs = mu_factor*( dR1 * Sigma * G1 -  dR2* Sigma * G2);

    rhs1 = mu_factor*( (dR1{1} .* sigma.')*G1 - (dR2{1} .* sigma.')*G2 ); %for p=1
    rhs2 = mu_factor*( (dR1{2} .* sigma.')*G1 - (dR2{2} .* sigma.')*G2 ); %for p=2
    rhs3 = mu_factor*( (dR1{3} .* sigma.')*G1 - (dR2{3} .* sigma.')*G2 ); %for p=3



    % Solve the adjoint problem for each sensor to get lambda vectors
    parfor m = 1:num_sensors

        [dlambdaX(:,m),~,~] = pcg(A_matrix,rhs1(m,:)',1e-9,numel(sigma),Mfun,Nfun);
        [dlambdaY(:,m),~,~] = pcg(A_matrix,rhs2(m,:)',1e-9,numel(sigma),Mfun,Nfun);
        [dlambdaZ(:,m),~,~] = pcg(A_matrix,rhs3(m,:)',1e-9,numel(sigma),Mfun,Nfun);
    end

    dlambda{dim,1} = dlambdaX;
    dlambda{dim,2} = dlambdaY;
    dlambda{dim,3} = dlambdaZ;

    dR1t{dim,1} = dR1{1};
    dR1t{dim,2} = dR1{2};
    dR1t{dim,3} = dR1{3};

    dR2t{dim,1} = dR2{1};
    dR2t{dim,2} = dR2{2};
    dR2t{dim,3} = dR2{3};
end

end



function dB = compute_dB_ds(img,dim,p)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);
num_stim = numel(img.fwd_model.stimulation);

mu_factor = img.fwd_model.mu0/(4*pi);

R = img.fwd_model.R;
G = img.fwd_model.G;

% sigma = img.elem_data;

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

switch dim
    case 1
        dRzdp = dRmkj_dp(img.fwd_model, 3, p);

        dRydp = dRmkj_dp(img.fwd_model, 2, p);
    otherwise
        error('Here');
end

Sigma = spdiags(img.elem_data(:), 0, length(img.elem_data), length(img.elem_data));

dB = zeros(num_sensors,num_stim);

for j = 1:num_stim
    % Gx_times_u = G.Gx*u(:,j);

    Gy_times_u = G.Gy*u(:,j);
    Gz_times_u = G.Gz*u(:,j);

    dB(:,j) = -mu_factor*( -dRydp*Sigma*Gz_times_u + dRzdp*Sigma*Gy_times_u);
end

end



function [dJ,dJ1,dJ2] = compute_dJ_ds(img,dim,p)
% Computes derivative of the Jacobian w.r.t. sensor positions
% Vectorized over elements and stimulations, loop only over sensors

num_sensors = numel(img.fwd_model.sensors);
n_elem      = size(img.fwd_model.elems,1);
num_stim    = numel(img.fwd_model.stimulation);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

% Element volumes
elemV = img.fwd_model.elem_volume(:);  % n_elem x 1

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% G*u for all stim, size: n_elem x num_stim
GxU = G.Gx*u;
GyU = G.Gy*u;
GzU = G.Gz*u;

% Compute derivative of lambda w.r.t sensor positions (output dR1 and dR2
% to avoid recomputing)
[dlambda,dR1dp,dR2dp] = compute_dlambda_ds(img,dim,p);

% Compute derivatives of R w.r.t sensor positions
switch dim
    case 1
        % dR1dp = dRmkj_dp(img.fwd_model, 3, p);
        % dR2dp = dRmkj_dp(img.fwd_model, 2, p);

        G1U = GyU;
        G2U = GzU;
    case 2
        % dR1dp = dRmkj_dp(img.fwd_model, 1, p);
        % dR2dp = dRmkj_dp(img.fwd_model, 3, p);

        G1U = GzU;
        G2U = GxU;
    case 3
        % dR1dp = dRmkj_dp(img.fwd_model, 2, p);
        % dR2dp = dRmkj_dp(img.fwd_model, 1, p);

        G1U = GxU;
        G2U = GyU;
    otherwise
        error('Dimension not supported.');
end


% Preallocate
dJ1 = zeros(num_sensors,num_stim,n_elem);
dJ2 = zeros(num_sensors,num_stim,n_elem);



%% Loop only over sensors
for m = 1:num_sensors
    %% --- dJ1: contribution from dlambda ---
    % dlambda * G matrices, size: n_elem x 1
    dlGx = G.Gx*dlambda(:,m);
    dlGy = G.Gy*dlambda(:,m);
    dlGz = G.Gz*dlambda(:,m);

    % Elementwise multiplication and sum over components
    % tmp: n_elem x num_stim
    tmp = dlGx .* GxU + dlGy .* GyU + dlGz .* GzU;

    % Multiply by element volumes, permute to [num_stim x n_elem]
    dJ1(m,:,:) = tmp.' .* elemV(:).';

    %% --- dJ2: contribution from dR/dp ---
    % dRydp, dRzdp: 1 x n_elem
    % Multiply by G*u per stimulation, permute to [num_stim x n_elem]
    dJ2(m,:,:) = -mu_factor * ( ...
        dR2dp(m,:) .* G2U.' - dR1dp(m,:) .* G1U.' ...
        );
end

%% Total derivative
dJ = dJ1 - dJ2;

end

function [dJ,dJ1,dJ2] = compute_dJ_xyz(img,dim)
% Computes derivative of the Jacobian w.r.t. sensor positions
% Vectorized over elements and stimulations, loop only over sensors

num_sensors = numel(img.fwd_model.sensors);
n_elem      = size(img.fwd_model.elems,1);
num_stim    = numel(img.fwd_model.stimulation);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

% Element volumes
elemV = img.fwd_model.elem_volume(:);  % n_elem x 1

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% G*u for all stim, size: n_elem x num_stim
GxU = G.Gx*u;
GyU = G.Gy*u;
GzU = G.Gz*u;

% Compute derivative of lambda w.r.t sensor positions (output dR1 and dR2
% to avoid recomputing)
[dlambda,dR1dp,dR2dp] = compute_dlambda_xyz(img,dim);

% Compute derivatives of R w.r.t sensor positions
switch dim
    case 1
        G1U = GyU;
        G2U = GzU;
    case 2
        G1U = GzU;
        G2U = GxU;
    case 3
        G1U = GxU;
        G2U = GyU;
    otherwise
        error('Dimension not supported.');
end

% Preallocate
dJ = cell(1,3);
dJ1 = zeros(num_sensors,num_stim,n_elem);
dJ2 = zeros(num_sensors,num_stim,n_elem);

%% Loop only over sensors
for p = 1:3
    for m = 1:num_sensors
        %% --- dJ1: contribution from dlambda ---
        % dlambda * G matrices, size: n_elem x 1
        dlGx = G.Gx*dlambda{p}(:,m);
        dlGy = G.Gy*dlambda{p}(:,m);
        dlGz = G.Gz*dlambda{p}(:,m);

        % Elementwise multiplication and sum over components
        % tmp: n_elem x num_stim
        tmp = dlGx .* GxU + dlGy .* GyU + dlGz .* GzU;

        % Multiply by element volumes, permute to [num_stim x n_elem]
        dJ1(m,:,:) = tmp.' .* elemV(:).';

        %% --- dJ2: contribution from dR/dp ---
        % dRydp, dRzdp: 1 x n_elem
        % Multiply by G*u per stimulation, permute to [num_stim x n_elem]
        dJ2(m,:,:) = -mu_factor * ( ...
            dR2dp{p}(m,:) .* G2U.' - dR1dp{p}(m,:) .* G1U.' ...
            );
    end

    %% Total derivative
    dJ{p} = dJ1 - dJ2;
end

end


% DEBUGGING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function test_dR_finite_difference(fmdl,j)
% TEST_DR_FINITE_DIFFERENCE
% Verifies the derivative computed by dRmkj_xyz_optimized using central
% finite differences w.r.t sensor positions.
%
% INPUTS:
%   fmdl - FEM model structure with nodes, elems, sensors
%   j    - the component (1=x, 2=y, 3=z) of the R derivative

numSensors = numel(fmdl.sensors);
eps_val    = 1e-6;      % small perturbation
rel_err    = zeros(numSensors,3);

% Compute analytic derivative
dR_analytic = dRmkj_xyz_optimized(fmdl,j);

% Loop over all sensors and coordinate directions
for m = 1:numSensors
    for p = 1:3
        % Perturb sensor m along coordinate p
        fmdl_plus  = fmdl;
        fmdl_minus = fmdl;
        
        fmdl_plus.sensors(m).position(p)  = fmdl_plus.sensors(m).position(p) + eps_val;
        fmdl_minus.sensors(m).position(p) = fmdl_minus.sensors(m).position(p) - eps_val;
        
        sensor_locations_plus = fetch_sensor_locations(fmdl_plus);
        sensor_locations_minus = fetch_sensor_locations(fmdl_minus);

        % Compute R matrices for perturbed sensors
        
        [Rx_plus,Ry_plus,Rz_plus]  = compute_r_matrices_local_optimized(fmdl,sensor_locations_plus);
        [Rx_minus,Ry_minus,Rz_minus] = compute_r_matrices_local_optimized(fmdl,sensor_locations_minus);
        
        % Finite difference approximation
        switch j
            case 1
                dR_fd = (Rx_plus(m,:)-Rx_minus(m,:))/(2*eps_val);
            case 2
                dR_fd = (Ry_plus(m,:)-Ry_minus(m,:))/(2*eps_val);

            case 3
                dR_fd = (Rz_plus(m,:)-Rz_minus(m,:))/(2*eps_val);
        end

        % Extract the corresponding analytic derivative
        %

        % Compare each component
        err = norm(dR_fd-dR_analytic{p}(m,:)) / norm(dR_fd);
        
        
        rel_err(m,p) = err;
        fprintf('Sensor %d, coord %d: relative error = %.3e\n', m, p, rel_err(m,p));
    end
end

figure;
imagesc(rel_err); colorbar;
xlabel('Coordinate (x=1,y=2,z=3)');
ylabel('Sensor index');
title('Relative error of dRmkj derivative vs finite difference');
end


% SEEMS LIKE THE DERIVATIVES FOR THE FIRST SENSOR COORDINATE ARE CORRECT;
% BUT NOT FOR THE SECOND AND THIRD SENSOR COORDINATES!!! THIS IS WHAT I
% HAVE FOUND UNTIL NOW.
function test_dlambda_finite_difference(img)
% Test the derivative of lambda w.r.t sensor positions using finite differences
% img: EIT image structure with fwd_model

sensor_locations = reshape([img.fwd_model.sensors.position],3,[])';
n_sensors = size(sensor_locations,1);
eps = 1e-6;  % small perturbation

% Compute analytic derivatives
[dlambda,~,~] = compute_dlambdaxyz_xyz_optimized(img);

fprintf('Testing dlambda derivatives with finite differences...\n');

for m = 1:n_sensors
    for p = 1:3
        % perturb sensor m along coordinate p
        sensor_plus  = sensor_locations;
        sensor_minus = sensor_locations;
        sensor_plus(m,p)  = sensor_plus(m,p) + eps;
        sensor_minus(m,p) = sensor_minus(m,p) - eps;

        img_plus  = assign_sensor_locations(img,sensor_plus);
        img_minus = assign_sensor_locations(img,sensor_minus);

        % recompute lambda
        lambda_plus  = compute_lambda_xyz(img_plus);
        lambda_minus = compute_lambda_xyz(img_minus);

        lambda_plus_x = lambda_plus{1};
        lambda_plus_y = lambda_plus{2};
        lambda_plus_z = lambda_plus{3};

        lambda_minus_x = lambda_minus{1};
        lambda_minus_y = lambda_minus{2};
        lambda_minus_z = lambda_minus{3};
        
        % finite difference
        dlambda_x_dp_fd = (lambda_plus_x(:,m) - lambda_minus_x(:,m))/(2*eps);
        dlambda_y_dp_fd = (lambda_plus_y(:,m) - lambda_minus_y(:,m))/(2*eps);
        dlambda_z_dp_fd = (lambda_plus_z(:,m) - lambda_minus_z(:,m))/(2*eps);


        % analytic derivative
        dlambda_x_dp_analytic = dlambda{1,p}(:,m);
        dlambda_y_dp_analytic = dlambda{2,p}(:,m);
        dlambda_z_dp_analytic = dlambda{3,p}(:,m);

        rel_err_x = norm(dlambda_x_dp_fd - dlambda_x_dp_analytic)/norm(dlambda_x_dp_fd);
        rel_err_y = norm(dlambda_y_dp_fd - dlambda_y_dp_analytic)/norm(dlambda_y_dp_fd);
        rel_err_z = norm(dlambda_z_dp_fd - dlambda_z_dp_analytic)/norm(dlambda_z_dp_fd);

        fprintf('Sensor %d, coordinate %d: relative error = %.3e\n', m, p, max([rel_err_x,rel_err_y,rel_err_z]));
    end
end


end

function test_single_sensor_derivative(img,A)

sensor_locations = reshape([img.fwd_model.sensors.position],3,[])';
n_sensors = size(sensor_locations,1);

[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);
J0 = [Jx;Jy;Jz];

dJ = compute_dJxyz_xyz_optimized(img);

for m = 1:n_sensors
    for p = 1:3
        % m = 1;     % choose sensor
        % p = 1;     % coordinate

        eps = 1e-6;

        sensor_plus  = sensor_locations;
        sensor_minus = sensor_locations;

        sensor_plus(m,p)  = sensor_plus(m,p)  + eps;
        sensor_minus(m,p) = sensor_minus(m,p) - eps;

        img_plus  = assign_sensor_locations(img,sensor_plus);
        img_minus = assign_sensor_locations(img,sensor_minus);

        [Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img_plus,A);
        J_plus = [Jx;Jy;Jz];

        [Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img_minus,A);
        J_minus = [Jx;Jy;Jz];

        dJ_fd = (J_plus - J_minus)/(2*eps);

        n_stim = numel(img.fwd_model.stimulation);
        block_size = n_sensors * n_stim;

        ids = m:n_sensors:block_size;

        dJ_pred = zeros(size(J0));

        dJ_pred(ids,:)              = squeeze(dJ{1,p}(m,:,:));
        dJ_pred(ids+block_size,:)   = squeeze(dJ{2,p}(m,:,:));
        dJ_pred(ids+2*block_size,:) = squeeze(dJ{3,p}(m,:,:));

        rel_err = norm(dJ_fd-dJ_pred,'fro')/norm(dJ_fd,'fro');

        fprintf('Relative error = %.3e\n',rel_err)

    end
end


end


function test_dJ_finite_difference(img,A)

sensor_locations = reshape([img.fwd_model.sensors.position],3,[])';

n_sensors = size(sensor_locations,1);
n_stim = numel(img.fwd_model.electrode);
n_elem = size(img.fwd_model.elems,1);

% Compute analytic derivatives
dJ = compute_dJxyz_xyz_optimized(img);

eps = 1e-6;

for m = 1:n_sensors
    for p = 1:3
        % perturb sensor m along coordinate p
        sensor_plus  = sensor_locations;
        sensor_minus = sensor_locations;
        sensor_plus(m,p)  = sensor_plus(m,p) + eps;
        sensor_minus(m,p) = sensor_minus(m,p) - eps;

        img_plus  = assign_sensor_locations(img,sensor_plus);
        img_minus = assign_sensor_locations(img,sensor_minus);

        [Jx_plus,Jy_plus,Jz_plus] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img_plus,A);

        [Jx_minus,Jy_minus,Jz_minus] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img_minus,A);

        % Central difference change
        dJx_dp_fd = (Jx_plus - Jx_minus) / (2*eps);
        dJy_dp_fd = (Jy_plus - Jy_minus) / (2*eps);
        dJz_dp_fd = (Jz_plus - Jz_minus) / (2*eps);
        
        % dJ{1,p}(1,:,1) derivatives for Jx, w.r.t. p coordinates: 1 -
        % sensor 1, : - any stim pattern, 1 - element 1

        % dJx_dp_fd_final(:,1) derivatives for Jx w.r.t. p coordinaes: :
        % any stim pattern, 1 - element 1 , for sensor m

        % What we have to compare is dJx_dp_fd_final(:,k) with
        % dJ{1,p}(m,:,k) at each m

        % Fetch the entries which are not null
        dJx_dp_fd_final = dJx_dp_fd(m:n_stim:n_stim*n_sensors,:);
        dJy_dp_fd_final = dJy_dp_fd(m:n_stim:n_stim*n_sensors,:);
        dJz_dp_fd_final = dJz_dp_fd(m:n_stim:n_stim*n_sensors,:);
        
        errors_x = zeros(1,n_elem);
        errors_y = zeros(1,n_elem);
        errors_z = zeros(1,n_elem);

        for k = 1:n_elem
            errors_x(k) = norm(dJ{1,p}(m,:,k) - dJx_dp_fd_final(:,k)')/norm(dJx_dp_fd_final(:,k));
            errors_y(k) = norm(dJ{2,p}(m,:,k) - dJy_dp_fd_final(:,k)')/norm(dJy_dp_fd_final(:,k));
            errors_z(k) = norm(dJ{3,p}(m,:,k) - dJz_dp_fd_final(:,k)')/norm(dJz_dp_fd_final(:,k));
        end

        fprintf('S(%d), coord(%d): rerrx = %.3e\n', m, p,max(errors_x));
        fprintf('S(%d), coord(%d): rerry = %.3e\n', m, p,max(errors_y));
        fprintf('S(%d), coord(%d): rerrz = %.3e\n', m, p,max(errors_z));
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dJ,dJ1,dJ2] = compute_dJxyz_xyz_optimized(img)
% Computes derivative of the Jacobian w.r.t. sensor positions
% Vectorized over elements and stimulations, loop only over sensors

n_sensors = numel(img.fwd_model.sensors);
n_elem      = size(img.fwd_model.elems,1);
num_stim    = numel(img.fwd_model.stimulation);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

% Element volumes
elemV = img.fwd_model.elem_volume(:);  
elemV = reshape(elemV,[1 1 n_elem]);

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% G*u for all stim, size: n_elem x num_stim
Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

% Expand u-terms to [1 x n_stim x n_elem]
GxU = reshape(Gx_times_u.', [1 num_stim n_elem]); 
GyU = reshape(Gy_times_u.', [1 num_stim n_elem]);
GzU = reshape(Gz_times_u.', [1 num_stim n_elem]);


% Compute derivative of lambda w.r.t sensor positions (output dR1 and dR2
% to avoid recomputing)
[dlambda,dR1dp,dR2dp] = compute_dlambdaxyz_xyz_optimized(img);

% THIS WAS WRONG BEFORE ;(. Quite hard to catch the bug!
% Precompute G*dlambda for all p. 
% 
% for dim = 1:3
%     DLx{dim} = reshape((G.Gx * dlambda{dim}).',[n_sensors 1 n_elem]);   
%     DLy{dim} = reshape((G.Gy * dlambda{dim}).',[n_sensors 1 n_elem]);
%     DLz{dim} = reshape((G.Gz * dlambda{dim}).',[n_sensors 1 n_elem]);
% end

% Careful here: dlambda is a 3 x 3 cell of dim x coord. dlambdaXdx =
% dlambda{1,1}
for dim = 1:3
    for p = 1:3
        DLx{dim,p} = reshape((G.Gx * dlambda{dim,p}).',[n_sensors 1 n_elem]);
        DLy{dim,p} = reshape((G.Gy * dlambda{dim,p}).',[n_sensors 1 n_elem]);
        DLz{dim,p} = reshape((G.Gz * dlambda{dim,p}).',[n_sensors 1 n_elem]);
    end
end

% Precompute G1U and G2U
G1U = {GyU, GzU, GxU};
G2U = {GzU, GxU, GyU};

% Preallocate
dJ = cell(3,3);

%% Loop only over sensors
for dim = 1:3
    for p = 1:3
        % --- dJ1 term (fully vectorized over sensors) ---
        

        tmp = DLx{dim,p} .* GxU + ...
            DLy{dim,p} .* GyU + ...
            DLz{dim,p} .* GzU;   % num_sensors xnum_stim x n_elem 

        % THIS WAS WRONG BEFORE ;(. Quite hard to catch the bug!
        % tmp = DLx{p} .* GxU + ...
        %     DLy{p} .* GyU + ...
        %     DLz{p} .* GzU;   % num_sensors xnum_stim x n_elem 

        % Expand element volumes
        dJ1 = tmp .* elemV;

        % --- dJ2 term ---

        term2 = ...
            reshape(dR2dp{dim,p},[n_sensors 1 n_elem]) .* G2U{dim} - ...
            reshape(dR1dp{dim,p},[n_sensors 1 n_elem]) .* G1U{dim};

        dJ2 = -mu_factor * term2;

        % Total
        dJ{dim,p} = dJ1 - dJ2;
    end
end
end

function dJ = compute_dJxyz_xyz(img)
% Computes derivative of the Jacobian w.r.t. sensor positions
% Vectorized over elements and stimulations, loop only over sensors

num_sensors = numel(img.fwd_model.sensors);
n_elem      = size(img.fwd_model.elems,1);
num_stim    = numel(img.fwd_model.stimulation);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

% Element volumes
elemV = img.fwd_model.elem_volume(:);  % n_elem x 1

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% G*u for all stim, size: n_elem x num_stim
GxU = G.Gx*u;
GyU = G.Gy*u;
GzU = G.Gz*u;

% Compute derivative of lambda w.r.t sensor positions (output dR1 and dR2
% to avoid recomputing)
[dlambda,dR1dp,dR2dp] = compute_dlambdaxyz_xyz(img);

% Compute derivatives of R w.r.t sensor positions
G1U = cell(1,3);
G2U = cell(1,3);

for dim = 1:3
    switch dim
        case 1
            G1U{dim} = GyU;
            G2U{dim} = GzU;
        case 2
            G1U{dim} = GzU;
            G2U{dim} = GxU;
        case 3
            G1U{dim} = GxU;
            G2U{dim} = GyU;
        otherwise
            error('Dimension not supported.');
    end
end

% Preallocate
dJ = cell(3,3);

dJ1 = zeros(num_sensors,num_stim,n_elem);
dJ2 = zeros(num_sensors,num_stim,n_elem);


%% Loop only over sensors
for dim = 1:3
    for p = 1:3
        for m = 1:num_sensors
            %% --- dJ1: contribution from dlambda ---
            % dlambda * G matrices, size: n_elem x 1
            dlGx = G.Gx*dlambda{p}(:,m);
            dlGy = G.Gy*dlambda{p}(:,m);
            dlGz = G.Gz*dlambda{p}(:,m);

            % Elementwise multiplication and sum over components
            % tmp: n_elem x num_stim
            tmp = dlGx .* GxU + dlGy .* GyU + dlGz .* GzU;

            % Multiply by element volumes, permute to [num_stim x n_elem]
            dJ1(m,:,:) = tmp.' .* elemV(:).';

            %% --- dJ2: contribution from dR/dp ---
            switch dim
                case 1
                    dJ2(m,:,:) = -mu_factor * ( ...
                        dR2dp{dim,p}(m,:) .* G2U{dim}.' - dR1dp{dim,p}(m,:) .* G1U{dim}.' ...
                        );
                case 2
                    dJ2(m,:,:) = -mu_factor * ( ...
                        dR2dp{dim,p}(m,:) .* G2U{dim}.' - dR1dp{dim,p}(m,:) .* G1U{dim}.' ...
                        );
                case 3
                    dJ2(m,:,:) = -mu_factor * ( ...
                        dR2dp{dim,p}(m,:) .* G2U{dim}.' - dR1dp{dim,p}(m,:) .* G1U{dim}.' ...
                        );
            end


        end

        %% Total derivative
        dJ{dim,p} = dJ1 - dJ2;
    end
end
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


function [Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local(img,A)

mu0 = img.fwd_model.mu0;

n_nodes =  size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices_local(img);

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

% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [n_elem 1 1]);
% Later this will broadcast to [numElems × numStim × numSensors]

% Expand u-terms to 3D
GxU = reshape(Gx_times_u, [n_elem num_stim 1]); % [: × numStim × 1]
GyU = reshape(Gy_times_u, [n_elem num_stim 1]);
GzU = reshape(Gz_times_u, [n_elem num_stim 1]);

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx.', [n_elem 1 num_sensors]);
Ry_ = reshape(R.Ry.', [n_elem 1 num_sensors]);
Rz_ = reshape(R.Rz.', [n_elem 1 num_sensors]);

% These are the derivatives with respect to sigma of the C components of the Gamma matrix,
dCxdp = ( -Rz_.*GyU + Ry_.*GzU );
dCydp = ( -Rx_.*GzU + Rz_.*GxU );
dCzdp = ( -Ry_.*GxU + Rx_.*GyU );

for select_sensor_axis = 1:3

    switch select_sensor_axis
        case 1
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_X, [n_elem 1 num_sensors]); % [: × 1 × numSensors]
            GyL = reshape(Gy_times_lambda_X, [n_elem 1 num_sensors]);
            GzL = reshape(Gz_times_lambda_X, [n_elem 1 num_sensors]);
        case 2
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_Y, [n_elem 1 num_sensors]); % [: × 1 × numSensors]
            GyL = reshape(Gy_times_lambda_Y, [n_elem 1 num_sensors]);
            GzL = reshape(Gz_times_lambda_Y, [n_elem 1 num_sensors]);
        case 3
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_Z, [n_elem 1 num_sensors]); % [: × 1 × numSensors]
            GyL = reshape(Gy_times_lambda_Z, [n_elem 1 num_sensors]);
            GzL = reshape(Gz_times_lambda_Z, [n_elem 1 num_sensors]);
    end

    % Compute all dfdx for all sensors+stim
    dfdx = elemV .* ( ...
        GxL.*GxU + ...
        GyL.*GyU + ...
        GzL.*GzU );

    % The g matrix does not depend on sigma.
    g = zeros(num_sensors,3,3);
    for m = 1:numel(img.fwd_model.sensors)
        g(m,:,:) = [...
            img.fwd_model.sensors(m).axes.axis1;
            img.fwd_model.sensors(m).axes.axis2;
            img.fwd_model.sensors(m).axes.axis3];
    end

    % g: [num_sensors × 3 × 3]
    gx = reshape(g(:,select_sensor_axis,1), [1 1 num_sensors]);
    gy = reshape(g(:,select_sensor_axis,2), [1 1 num_sensors]);
    gz = reshape(g(:,select_sensor_axis,3), [1 1 num_sensors]);

    dfdp = mu_factor*(...
        gx.*dCxdp +...
        gy.*dCydp +...
        gz.*dCzdp);

    dfd = dfdx + dfdp;   % size: [numElems × numStim × numSensors]

    % Now reshape to match J(ids,:)
    % permute to [numSensors × numStim × numElems]
    dfd = permute(dfd, [3 2 1]);

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


function J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,select_sensor_axis)

mu0 = img.fwd_model.mu0;

n_nodes =  size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices_local(img);

R = img.fwd_model.R;
G = img.fwd_model.G;

switch select_sensor_axis
    case 1
        Gamma = img.Gamma1;
        R1 = R.Rz.';
        R2 = R.Ry.';
    case 2
        Gamma = img.Gamma2;
        R1 = R.Rx.';
        R2 = R.Rz.';
    case 3
        Gamma = img.Gamma3;
        R1 = R.Ry.';
        R2 = R.Rx.';
    otherwise
        error('here')
end

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% Solve the adjoint problem for each sensor to get lambda vectors
lambda = zeros(n_nodes,num_sensors);

A_matrix = A(img.elem_data);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

GammaT = Gamma.';

% Incomplete Cholesky factorization preconditioner seems to be a bit faster
% than Jacobi preconditioner. However, it breaks down when the
% conductivities become negative
% R = ichol(A_matrix);
% Rt = R';

num_elements = numel(img.elem_data);
for m = 1:num_sensors
    [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-6,num_elements,Mfun,Nfun);
end

Gx_times_lambda = G.Gx*lambda;
Gy_times_lambda = G.Gy*lambda;
Gz_times_lambda = G.Gz*lambda;

Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [n_elem 1 1]);
% Later this will broadcast to [numElems × numStim × numSensors]

% Expand lambda and R terms to 3D
GxL = reshape(Gx_times_lambda, [n_elem 1 num_sensors]); % [: × 1 × numSensors]
GyL = reshape(Gy_times_lambda, [n_elem 1 num_sensors]);
GzL = reshape(Gz_times_lambda, [n_elem 1 num_sensors]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u, [n_elem num_stim 1]); % [: × numStim × 1]
GyU = reshape(Gy_times_u, [n_elem num_stim 1]);
GzU = reshape(Gz_times_u, [n_elem num_stim 1]);

% Compute all dfdx for all sensors+stim
dfdx = elemV .* ( ...
    GxL.*GxU + ...
    GyL.*GyU + ...
    GzL.*GzU );

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx.', [n_elem 1 num_sensors]);
Ry_ = reshape(R.Ry.', [n_elem 1 num_sensors]);
Rz_ = reshape(R.Rz.', [n_elem 1 num_sensors]);

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

% g: [num_sensors × 3 × 3]
gx = reshape(g(:,select_sensor_axis,1), [1 1 num_sensors]);
gy = reshape(g(:,select_sensor_axis,2), [1 1 num_sensors]);
gz = reshape(g(:,select_sensor_axis,3), [1 1 num_sensors]);

dfdp = mu_factor*(...
    gx.*dCxdp +...
    gy.*dCydp +...
    gz.*dCzdp);

dfd = dfdx + dfdp;   % size: [numElems × numStim × numSensors]

% Now reshape to match J(ids,:)
% permute to [numSensors × numStim × numElems]
dfd = permute(dfd, [3 2 1]);

% collapse first 2 dims → [numSensors*numStim × numElems]
J = reshape(dfd, num_sensors*num_stim, n_elem);

return
end

% For 1 axis-MDEIT
function cost = compute_cost_function_d_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

L = chol(H,'lower');
logdetH = 2*sum(log(diag(L)));

cost = -logdetH;
end

function cost = compute_cost_function_a_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% L = chol(H,'lower');
% cost = sum(1./diag(L).^2);

cost = trace(inv(H));
end


% For 3 axis-MDEIT
function cost = compute_cost_function_d_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A)
% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);

J = [Jx;Jy;Jz];

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

L = chol(H,'lower');
logdetH = 2*sum(log(diag(L)));

cost = -logdetH;
end


function cost = compute_cost_function_a_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A)

n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);

J = [Jx;Jy;Jz];

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;


% cost = trace(inv(H));

% CHANGED HERE!
L = chol(H,'lower');
Hinv = L'\(L\eye(n_elem));
cost = trace(Hinv);

end





function G = repulsion_cost(sensor_locations)

n_sensors = size(sensor_locations,1);

% --- pairwise differences ---
% Expand coordinates
X = sensor_locations(:,1);
Y = sensor_locations(:,2);
Z = sensor_locations(:,3);

DX = X - X.';    % [n_sensors x n_sensors]
DY = Y - Y.';
DZ = Z - Z.';

% Squared distances
D2 = DX.^2 + DY.^2 + DZ.^2+ 1e-12; %smooth correction

% --- Cost ---
invD2 = 1 ./ D2;

% remove diagonal
invD2(1:n_sensors+1:end) = 0;

% sum only i<j. We sum repeated pairs (i,j) and (j,i), so we multiply by
% 0.5
G = 0.5 * sum(invD2(:));

end

function [dGx,dGy,dGz] = repulsion_grad_cartesian(sensor_locations)

n_sensors = size(sensor_locations,1);

% --- pairwise differences ---
% Expand coordinates
X = sensor_locations(:,1);
Y = sensor_locations(:,2);
Z = sensor_locations(:,3);

DX = X - X.';    % [n_sensors x n_sensors]
DY = Y - Y.';
DZ = Z - Z.';

% Squared distances
D2 = DX.^2 + DY.^2 + DZ.^2+ 1e-12; %smooth correction

% --- Cost ---
invD2 = 1 ./ D2;

% remove diagonal
invD2(1:n_sensors+1:end) = 0;

% We need 1/d^4
invD4 = invD2.^2;

% Compute force matrix coefficients
C = -2 * invD4;

% zero diagonal
C(1:n_sensors+1:end) = 0;

% Each gradient component is:
% dG_i = sum_j C_ij * (p_i - p_j)

dGx = sum(C .* DX, 2);
dGy = sum(C .* DY, 2);
dGz = sum(C .* DZ, 2);

end



% For 1-axis MDEIT
function dphidp = compute_cost_function_gradient_d_opt_optimized(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);
dJds = compute_dJ_xyz(img,dim);

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves

% Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
% Y = L'\eye(n_elem);
% Hinv = L\Y;

Hinv = L'\(L\eye(n_elem));

W = Gamma_noise_inv*J*Hinv;

dphidp = zeros(3,n_sensors);

for p = 1:3
    for m = 1:n_sensors
        ids = m : n_sensors : n_sensors*n_stim;

        dJm = reshape(dJds{p}(m,:,:), [n_stim, n_elem]);
        Wm  = W(ids,:);

        dphidp(p,m) = -2 * sum(Wm(:) .* dJm(:));
    end
end

end

function dphidp = compute_cost_function_gradient_a_opt_optimized(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);
dJds = compute_dJ_xyz(img,dim);

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves


% % Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
% Y = L'\eye(n_elem);
% X = L\Y;
% W = L'\X;
% inv_H2 = L\W;

Hinv = L'\(L\eye(n_elem));
inv_H2 = Hinv*Hinv;

W = Gamma_noise_inv*J*inv_H2;

dphidp = zeros(3,n_sensors);

for p = 1:3
    for m = 1:n_sensors
        ids = m : n_sensors : n_sensors*n_stim;

        dJm = reshape(dJds{p}(m,:,:), [n_stim, n_elem]);
        Wm  = W(ids,:);

        dphidp(p,m) = -2 * sum(Wm(:) .* dJm(:));
    end
end

end

% For 3-axis MDEIT
function dphidp = compute_cost_function_gradient_a_opt_optimized_3_axis(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);

% DEBUG!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_single_sensor_derivative(img,A)

% test_dJ_finite_difference(img,A) % THE ERROR WAS HERE. CHECK THE RELEVANT
% FUNCTION (compute_dJ_xyz_xyz) FOR THE CHANGES!

% test_dlambda_finite_difference(img) % Its correct as well, so the error
% must be above

% These are correct!!!!!! So the error must be somewhere above
% test_dR_finite_difference(img.fwd_model,1);
% test_dR_finite_difference(img.fwd_model,2);
% test_dR_finite_difference(img.fwd_model,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dJxds2 = compute_dJ_xyz(img,1);
dJ = compute_dJxyz_xyz_optimized(img);

dJxds = cell(1,3);
dJyds = cell(1,3);
dJzds = cell(1,3);

for p = 1:3
    dJxds{p} = dJ{1,p};
    dJyds{p} = dJ{2,p};
    dJzds{p} = dJ{3,p};
end

J = [Jx;Jy;Jz];

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves

% THIS WAS WRONG!!!!!!!!!!!
% Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
% Y = L'\eye(n_elem);
% X = L\Y;
% W = L'\X;
% inv_H2 = L\W;

Hinv = L'\(L\eye(n_elem));
inv_H2 = Hinv*Hinv;

W = Gamma_noise_inv*J*inv_H2;

dphidp = zeros(3,n_sensors);

block_size = n_sensors * n_stim;

for p = 1:3
    for m = 1:n_sensors
        % indices inside one block (x, y, or z)
        ids_local = m : n_sensors : block_size;

        % --- X component ---
        ids_x = ids_local;
        dJx_m = reshape(dJxds{p}(m,:,:), [n_stim, n_elem]);
        Wx_m  = W(ids_x,:);
        term_x = sum(Wx_m(:) .* dJx_m(:));

        % --- Y component ---
        ids_y = ids_local + block_size;
        dJy_m = reshape(dJyds{p}(m,:,:), [n_stim, n_elem]);
        Wy_m  = W(ids_y,:);
        term_y = sum(Wy_m(:) .* dJy_m(:));

        % --- Z component ---
        ids_z = ids_local + 2*block_size;
        dJz_m = reshape(dJzds{p}(m,:,:), [n_stim, n_elem]);
        Wz_m  = W(ids_z,:);
        term_z = sum(Wz_m(:) .* dJz_m(:));

        dphidp(p,m) = -2 * (term_x + term_y + term_z);
    end
end

% % Precompute index offsets
% ids_base = reshape(1:block_size, n_sensors, n_stim);
% 
% for p = 1:3
% 
%     dJx = squeeze(dJ{1,p});   % [n_sensors x n_stim x n_elem]
%     dJy = squeeze(dJ{2,p});
%     dJz = squeeze(dJ{3,p});
% 
%     for m = 1:n_sensors
% 
%         ids_local = ids_base(m,:);
% 
%         Wx = W(ids_local, :);
%         Wy = W(ids_local + block_size, :);
%         Wz = W(ids_local + 2*block_size, :);
% 
%         % Use Frobenius inner products (faster than elementwise sum)
%         term_x = sum(sum(Wx .* squeeze(dJx(m,:,:))));
%         term_y = sum(sum(Wy .* squeeze(dJy(m,:,:))));
%         term_z = sum(sum(Wz .* squeeze(dJz(m,:,:))));
% 
%         dphidp(p,m) = -2 * (term_x + term_y + term_z);
%     end
% end


end

function dphidp = compute_cost_function_gradient_d_opt_optimized_3_axis(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);

% dJxds2 = compute_dJ_xyz(img,1);
dJ = compute_dJxyz_xyz_optimized(img);

dJxds = cell(1,3);
dJyds = cell(1,3);
dJzds = cell(1,3);

for p = 1:3
    dJxds{p} = dJ{1,p};
    dJyds{p} = dJ{2,p};
    dJzds{p} = dJ{3,p};
end

J = [Jx;Jy;Jz];

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves

% Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
% Y = L'\eye(n_elem);
% Hinv = L\Y;


Hinv = L'\(L\eye(n_elem));
% inv_H2 = Hinv*Hinv;


W = Gamma_noise_inv*J*Hinv;

dphidp = zeros(3,n_sensors);

block_size = n_sensors * n_stim;

for p = 1:3
    for m = 1:n_sensors
        % indices inside one block (x, y, or z)
        ids_local = m : n_sensors : block_size;

        % --- X component ---
        ids_x = ids_local;
        dJx_m = reshape(dJxds{p}(m,:,:), [n_stim, n_elem]);
        Wx_m  = W(ids_x,:);
        term_x = sum(Wx_m(:) .* dJx_m(:));

        % --- Y component ---
        ids_y = ids_local + block_size;
        dJy_m = reshape(dJyds{p}(m,:,:), [n_stim, n_elem]);
        Wy_m  = W(ids_y,:);
        term_y = sum(Wy_m(:) .* dJy_m(:));

        % --- Z component ---
        ids_z = ids_local + 2*block_size;
        dJz_m = reshape(dJzds{p}(m,:,:), [n_stim, n_elem]);
        Wz_m  = W(ids_z,:);
        term_z = sum(Wz_m(:) .* dJz_m(:));

        dphidp(p,m) = -2 * (term_x + term_y + term_z);
    end
end
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

if isfield(mdl,'fwd_model')
    fmdl = mdl.fwd_model;
else
    fmdl = mdl;
end
    

n_sensors = numel(fmdl.sensors);
sensor_locations = zeros(n_sensors,3);

for m = 1: numel(fmdl.sensors)
    sensor_locations(m,:) = fmdl.sensors(m).position;
end

end

function data_noisy = add_measurement_noise_difference(datai,datah,SNRdb)

if isempty(SNRdb)
    data_noisy = datai-datah;
    return;
else
    d = datai-datah;
    d_amplitude = max(abs(d-mean(d)));

    noise_level= d_amplitude/10^(SNRdb/20);

    data_noisy = d+noise_level*randn(size(d));
end
end