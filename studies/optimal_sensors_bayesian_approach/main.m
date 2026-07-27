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

model_parameters = create_kai_3d_model_parameters(l0, z0, sigma0, I0);

model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/8;
model_parameters.numOfElectrodesPerRing = 8;
model_parameters.numOfRings = 4;
model_parameters.numOfSensors = 4;
model_parameters.sensorRadius = model_parameters.radius*1.5;
model_parameters.material.name = 'sphere';
model_parameters.material.radius = model_parameters.radius/3;
model_parameters.material.type = 'spherical';
model_parameters.material.position(3) = 3/4*model_parameters.height;

background_conductivity = 3.28e-1/sigma0;
anomaly_conductivity = background_conductivity/10;

%% Simulation parameters

do_3_axis = false;
dim = 3;            %if doing 1-axis, which dimmension?

%Optimizer parameters
max_iteratons = 50;
hessian_approximation = 'lbfgs';
use_parallel = true;
algorithm = 'quasi-newton';

rmax = 2*model_parameters.radius;
rmin = model_parameters.radius*1.1;

zmax = model_parameters.height;
zmin = 0;

r0 = 1.5*model_parameters.radius; %radius of cylinder shell to perform optimizaiton on

alpha = 0; %controls the force of the repulsion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%% Make forward model

[~,fmdls] = mk_mdeit_model(model_parameters,model_folder);

fmdl = fmdls{1};

sensorLocations = zeros(numel(fmdl.sensors),3);
for i = 1: numel(fmdl.sensors)
    sensorLocations(i,:) = fmdl.sensors(i).position;
end

% Make homogeneous image
imgh = mk_image_mdeit(fmdl,background_conductivity);

% Add plastic cylinder
imgi = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

% show_fem(imgi);
% plot_sensors(imgi);

n_sensors = numel(imgi.fwd_model.sensors);
n_stim = numel(imgi.fwd_model.stimulation);
n_elem = size(imgi.fwd_model.elems,1);
n_nodes = size(fmdl.nodes,1);

A = @(x) M(imgi,x);

%% Initialize

% Random starting positions in radius r = 1.5;
R0 = model_parameters.radius*1.5;
r_new = R0*ones(1,n_sensors);

z_new = model_parameters.height*rand(1,n_sensors);
theta_new = 2*pi*rand(1,n_sensors);

sensor_locations_0 = [(r_new.*cos(theta_new))',(r_new.*sin(theta_new))',z_new'];

img_init = assign_sensor_locations(imgi,sensor_locations_0);

R0 = max(sqrt(sensor_locations_0(:,1).^2+sensor_locations_0(:,2).^2));
z0 = max(sensor_locations_0(:,3));

%% Define prior and noise covariance matrices

B = fwd_solve_mdeit(imgi);
max_B = max([abs(B.Bx(:));abs(B.By(:));abs(B.Bz(:))]);

% Set the noise variance with respect to the data magnitude
noise_std_deviation = max_B/10;
variance_noise = noise_std_deviation^2;



Jdim = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,dim);

coeff = 50;
variance_prior = coeff*variance_noise/eigs(Jdim'*Jdim,1);

% Check how many eigenvectors of J'J are in the data-dominated
% regime, as opposed to the prior-dominated regime
d = sum(eigs(Jdim'*Jdim,n_elem).*variance_prior/variance_noise>1);

fprintf('# lambda_i*alpha:beta>1 = %i (%i)\n',d,n_elem);

% White noise with zero mean and \mu variance
Gamma_noise = variance_noise.*speye(n_stim*n_sensors,n_stim*n_sensors);
inv_Gamma_noise = inv(Gamma_noise);

Gamma_prior = variance_prior.*speye(n_elem,n_elem);
inv_Gamma_prior = inv(Gamma_prior);


Jx = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,1);
Jy = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,2);
Jz = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,3);

J_3axis = [Jx;Jy;Jz];

coeff = 50;
variance_prior_3_axis = coeff*variance_noise/eigs(J_3axis'*J_3axis,1);

d_3axis = sum(eigs(J_3axis'*J_3axis,n_elem).*variance_prior_3_axis/variance_noise>1);

fprintf('# lambda_i*alpha:beta>1 = %i (%i)\n',d_3axis,n_elem);

Gamma_prior_3_axis = variance_prior_3_axis.*speye(n_elem,n_elem);
inv_Gamma_prior_3_axis = inv(Gamma_prior_3_axis);

Gamma_noise_3_axis = variance_noise.*speye(3*n_stim*n_sensors,3*n_stim*n_sensors);
inv_Gamma_noise_3_axis = inv(Gamma_noise_3_axis);


%% Define functions
jacobian_coordinate_transformation_cylindrical = compute_jacobian_coordinate_transformation_cylindrical(sensor_locations_0);

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



x0 = sensor_locations_to_vector_cartesian(sensor_locations_0);

% Full map from q coordinates to sensor locations and back for constrained
% optimization
vector_to_sensor_locations_con = @(q) vector_to_sensor_locations_cartesian(q_to_x_cylindrical(q));
sensor_locations_to_vector_con = @(sensor_locations) x_to_q_cylindrical(sensor_locations_to_vector_cartesian(sensor_locations));


% Full map from q coordinates to sensor locations and back for
% unconstrained optimization
vector_to_sensor_locations_unc = @(q) vector_to_sensor_locations_cartesian(q_to_x_cyl_region(q));
sensor_locations_to_vector_unc = @(sensor_locations) x_to_q_cyl_region(sensor_locations_to_vector_cartesian(sensor_locations));

% Sanity check
q = x_to_q_cylindrical(x0);
x = q_to_x_cylindrical(q);
assert(norm(x-x0)<1e-5,'Unexpected');

% Sanity check
q0 = sensor_locations_to_vector_con(sensor_locations_0);
sensor_locations_new = vector_to_sensor_locations_con(q0);
assert(norm(sensor_locations_new-sensor_locations_0)<1e-5,'Unexpected');

% Sanity check
sensor_locations = vector_to_sensor_locations_con(q0);
q_new = sensor_locations_to_vector_con(sensor_locations);
assert(norm(q0-q_new)<1e-5,'Unexpected');

% Transform gradient of cost function w.r.t. cartesian coordinates to
% w.r.t. q coordinates
function dphidq = dphidp_to_dphidq(sensor_locations,dphidp,jacobian_coordinate_transformation)

n_sensors = size(sensor_locations,1);

assert(all(size(jacobian_coordinate_transformation) == [3,3,n_sensors]));

dphidq = zeros(3,n_sensors);

for q = 1:3
    temp = zeros(1,n_sensors);
    for dim = 1:3
        temp = temp + squeeze(jacobian_coordinate_transformation(q,dim,:)).'.*dphidp(dim,:);
    end
    dphidq(q,:) = temp;
end

end


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
    vector_to_sensor_locations,jacobian_coordinate_transformation,opt_mode,mode,dim)

n_sensors = numel(img.fwd_model.sensors);

allowed_opt_modes = {'a-opt','d-opt'};
assert(ismember(opt_mode,allowed_opt_modes));

allowed_modes = {'mdeit3','mdeit1'};
assert(ismember(mode,allowed_modes));

if strcmp(mode,'mdeit3') && nargin<10
    dim = 'default';
end

sensor_locations = vector_to_sensor_locations(q);

% The cylindrical chain-rule Jacobian depends on the current (r,theta) of
% each sensor, not just the initial configuration: it must be re-evaluated
% at the current iterate (the passed-in jacobian_coordinate_transformation
% was previously stale, captured once at sensor_locations_0).
jacobian_coordinate_transformation = ...
    compute_jacobian_coordinate_transformation_cylindrical(sensor_locations);

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
dphidq = dphidp_to_dphidq(sensor_locations,dphidp,jacobian_coordinate_transformation);

dphi = [dphidq(1,:),dphidq(2,:),dphidq(3,:)];


end

f_a_opt_mdeit_dim_con = @(q) f(imgi,q,inv_Gamma_prior,inv_Gamma_noise,A,vector_to_sensor_locations_con,'a-opt','mdeit1',dim);
g_a_opt_mdeit_dim_con  = @(q) g(imgi,q,inv_Gamma_prior,inv_Gamma_noise,A,...
    vector_to_sensor_locations_con,jacobian_coordinate_transformation_cylindrical,'a-opt','mdeit1',dim);

f_d_opt_mdeit_dim_con = @(q) f(imgi,q,inv_Gamma_prior,inv_Gamma_noise,A,vector_to_sensor_locations_con,'d-opt','mdeit1',dim);
g_d_opt_mdeit_dim_con  = @(q) g(imgi,q,inv_Gamma_prior,inv_Gamma_noise,A,...
    vector_to_sensor_locations_con,jacobian_coordinate_transformation_cylindrical,'d-opt','mdeit1',dim);

f_a_opt_mdeit3_con = @(q) f(imgi,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,vector_to_sensor_locations_con,'a-opt','mdeit3',dim);
g_a_opt_mdeit3_con  = @(q) g(imgi,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,...
    vector_to_sensor_locations_con,jacobian_coordinate_transformation_cylindrical,'a-opt','mdeit3',dim);

f_d_opt_mdeit3_con = @(q) f(imgi,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,vector_to_sensor_locations_con,'d-opt','mdeit3',dim);
g_d_opt_mdeit3_con  = @(q) g(imgi,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,...
    vector_to_sensor_locations_con,jacobian_coordinate_transformation_cylindrical,'d-opt','mdeit3',dim);

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

%% Define function+gradient function

function [func,grad] = funcwithgrad(q,f_impl,g_impl)
% Calculate objective f
func = f_impl(q);

if nargout > 1 % gradient required
    grad =  g_impl(q);
end
end

%% Optimize 3-axis A-optimality with fmincon in cylinder shell
fprintf('(Shell) 3-axis A-optimality OED - fmincon\n')

options = optimoptions('fmincon',...
    'OptimalityTolerance',1e-9,'Display','iter','MaxIterations',max_iteratons, ...
    'Algorithm','interior-point','HessianApproximation',hessian_approximation,...
    'SpecifyObjectiveGradient',true,'UseParallel',use_parallel);

% Initial conditionss
q0_con = sensor_locations_to_vector_con(sensor_locations_0);

fun_a_opt_mdeit_1 = @(q) funcwithgrad(q,f_a_opt_mdeit3_con,g_a_opt_mdeit3_con);

lb = [ r0*ones(1,n_sensors),   -inf(1,n_sensors),   0*ones(1,n_sensors)];
ub = [ r0*ones(1,n_sensors),    inf(1,n_sensors),   model_parameters.height*ones(1,n_sensors) ];

x_a_opt_con_shell = fmincon(fun_a_opt_mdeit_1,q0_con,[],[],[],[],lb,ub,[],options);
img_a_opt_con_shell = assign_sensor_locations(imgi,vector_to_sensor_locations_con(x_a_opt_con_shell));


%% Optimize 3-axis D-optimality with fmincon in cylinder shell
fprintf('(Shell) 3-axis A-optimality OED - fmincon\n')

options = optimoptions('fmincon',...
    'OptimalityTolerance',1e-9,'Display','iter','MaxIterations',max_iteratons, ...
    'Algorithm','interior-point','HessianApproximation',hessian_approximation,...
    'SpecifyObjectiveGradient',true,'UseParallel',use_parallel);

% Initial conditionss
q0_con = sensor_locations_to_vector_con(sensor_locations_0);

fun_d_opt_mdeit_1 = @(q) funcwithgrad(q,f_d_opt_mdeit3_con,g_d_opt_mdeit3_con);

lb = [ r0*ones(1,n_sensors),   -inf(1,n_sensors),   0*ones(1,n_sensors)];
ub = [ r0*ones(1,n_sensors),    inf(1,n_sensors),   model_parameters.height*ones(1,n_sensors) ];

x_d_opt_con_shell = fmincon(fun_d_opt_mdeit_1,q0_con,[],[],[],[],lb,ub,[],options);
img_d_opt_con_shell = assign_sensor_locations(imgi,vector_to_sensor_locations_con(x_d_opt_con_shell));



%% Optimize 1-axis A-optimality with fmincon in cylinder region
fprintf('(Region) 1-axis A-optimality OED - fmincon\n')

options = optimoptions('fmincon',...
    'OptimalityTolerance',1e-9,'Display','iter','MaxIterations',max_iteratons, ...
    'Algorithm','interior-point','HessianApproximation',hessian_approximation,...
    'SpecifyObjectiveGradient',true,'UseParallel',use_parallel);

% Initial conditions
q0_con = sensor_locations_to_vector_con(sensor_locations_0);

fun_a_opt_mdeit_1 = @(q) funcwithgrad(q,f_a_opt_mdeit_dim_con,g_a_opt_mdeit_dim_con);

lb = [ rmin*ones(1,n_sensors),   -inf(1,n_sensors),   0*ones(1,n_sensors)];
ub = [ rmax*ones(1,n_sensors),    inf(1,n_sensors),   model_parameters.height*ones(1,n_sensors) ];

x_a_opt_con = fmincon(fun_a_opt_mdeit_1,q0_con,[],[],[],[],lb,ub,[],options);
img_a_opt_con_region = assign_sensor_locations(imgi,vector_to_sensor_locations_con(x_a_opt_con));

%% Optimize 1-axis D-optimality with fmincon in cylinder region
fprintf('(Region) 1-axis D-optimality OED - fmincon\n')

options = optimoptions('fmincon',...
    'OptimalityTolerance',1e-9,'Display','iter','MaxIterations',max_iteratons, ...
    'Algorithm','interior-point','HessianApproximation',hessian_approximation,...
    'SpecifyObjectiveGradient',true,'UseParallel',use_parallel);

% Initial conditions
q0_con = sensor_locations_to_vector_con(sensor_locations_0);


fun_d_opt_mdeit_1 = @(q) funcwithgrad(q,f_d_opt_mdeit_dim_con,g_d_opt_mdeit_dim_con);

lb = [ rmin*ones(1,n_sensors),   -inf(1,n_sensors),   0*ones(1,n_sensors)];
ub = [ rmax*ones(1,n_sensors),    inf(1,n_sensors),   model_parameters.height*ones(1,n_sensors) ];

x_d_opt_con = fmincon(fun_d_opt_mdeit_1,q0_con,[],[],[],[],lb,ub,[],options);
img_d_opt_con_region = assign_sensor_locations(imgi,vector_to_sensor_locations_con(x_d_opt_con));

%% Optimize 3-axis A-optimality with fmincon in cylinder region
fprintf('(Region) 3-axis A-optimality OED - fmincon\n')

options = optimoptions('fmincon',...
    'OptimalityTolerance',1e-9,'Display','iter','MaxIterations',max_iteratons, ...
    'Algorithm','interior-point','HessianApproximation',hessian_approximation,...
    'SpecifyObjectiveGradient',true,'UseParallel',use_parallel);

% Initial conditions
q0_con = sensor_locations_to_vector_con(sensor_locations_0);

fun_a_opt_mdeit_3 = @(q) funcwithgrad(q,f_a_opt_mdeit3_con,g_a_opt_mdeit3_con);

lb = [ rmin*ones(1,n_sensors),   -inf(1,n_sensors),   0*ones(1,n_sensors)];
ub = [ rmax*ones(1,n_sensors),    inf(1,n_sensors),   model_parameters.height*ones(1,n_sensors) ];

x_a_opt_con = fmincon(fun_a_opt_mdeit_3,q0_con,[],[],[],[],lb,ub,[],options);
img_a_opt_con_mdeit3_region = assign_sensor_locations(imgi,vector_to_sensor_locations_con(x_a_opt_con));

%% Optimize 3-axis D-optimality with fmincon in cylinder region
fprintf('(Region) 3-axis D-optimality OED - fmincon\n')

options = optimoptions('fmincon',...
    'OptimalityTolerance',1e-9,'Display','iter','MaxIterations',max_iteratons, ...
    'Algorithm','interior-point','HessianApproximation',hessian_approximation,...
    'SpecifyObjectiveGradient',true,'UseParallel',use_parallel);

% Initial conditions
q0_con = sensor_locations_to_vector_con(sensor_locations_0);

fun_d_opt_mdeit_3 = @(q) funcwithgrad(q,f_d_opt_mdeit3_con,g_d_opt_mdeit3_con);

lb = [ rmin*ones(1,n_sensors),   -inf(1,n_sensors),   0*ones(1,n_sensors)];
ub = [ rmax*ones(1,n_sensors),    inf(1,n_sensors),   model_parameters.height*ones(1,n_sensors) ];

x_d_opt_con = fmincon(fun_d_opt_mdeit_3,q0_con,[],[],[],[],lb,ub,[],options);
img_d_opt_con_mdeit3_region = assign_sensor_locations(imgi,vector_to_sensor_locations_con(x_d_opt_con));

%% Plot the optimal sensor locations
theta = linspace(0,2*pi,100);
[x, y, z] = cylinder(rmin, 50);
[x2, y2, z2] = cylinder(rmax, 50);
z = z * 3;
z2 = z2*3;

figure
hold on
show_fem(imgi);
plot_sensors(img_init);
% plot_sensors(img_a_opt,false,'g','s');
% plot_sensors(img_d_opt,false,'r','o');
plot_sensors(img_d_opt_con_region,false,'r','s');
plot_sensors(img_a_opt_con_region,false,'g','d');
plot_sensors(img_a_opt_con_mdeit3_region,false,'b','o');
h = plot_sensors(img_d_opt_con_mdeit3_region,false,'y','+');

hold on
hSurf = surf(h,x, y, z);
set(hSurf, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', 'red')
hSurf2 = surf(h,x2, y2, z2);
set(hSurf2, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', 'blue')
plot3(h,R0*cos(theta),R0*sin(theta),model_parameters.height/2*ones(size(theta)),'LineStyle','--','Color','b')
axis([-1.1*rmax 1.1*rmax -1.1*rmax 1.1*rmax 0 model_parameters.height])
hold off
box on;grid on;
camlight; lighting gouraud
drawnow








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

function lambda = compute_lambda_x(img)

n_nodes = size(img.fwd_model.nodes,1);
num_sensors = numel(img.fwd_model.sensors);

sigma = img.elem_data;

% Compute Gamma matrices
img = compute_gamma_matrices_local(img);
Gamma = img.Gamma1;

% Solve the adjoint problem for each sensor to get lambda vectors
lambda = zeros(n_nodes,num_sensors);

A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

GammaT = Gamma.';

parfor m = 1:num_sensors
    [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-6,numel(sigma),Mfun,Nfun);
end

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

function [Rx,Ry,Rz,fmdl] = compute_r_matrices_local(fmdl,sensor_locations)
% Vectorized over elements and quadrature points (loop only over sensors).
% Replaces the previous interpreted triple loop (sensors x elements x 35
% quad points), which is ~100x slower for meshes of a few thousand elements.
numElements = size(fmdl.elems,1);
numSensors = size(sensor_locations,1);

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

numQuadPts = length(weights);

V = reshape(fmdl.nodes(fmdl.elems',:),[4,numElements,3]);
v1 = squeeze(V(1,:,:)); v2 = squeeze(V(2,:,:));
v3 = squeeze(V(3,:,:)); v4 = squeeze(V(4,:,:));

a = v2-v1; b = v3-v1; c = v4-v1;
detJ = abs(sum(a .* cross(b,c,2),2))';   % 1 x numElements

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
            % NOTE: dlambda must be indexed {dim,p} (dlambda is a 3x3 cell;
            % linear indexing {p} previously gave the wrong adjoint
            % derivative for most (dim,p) combinations)
            dlGx = G.Gx*dlambda{dim,p}(:,m);
            dlGy = G.Gy*dlambda{dim,p}(:,m);
            dlGz = G.Gz*dlambda{dim,p}(:,m);

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
Jx = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,1);
Jy = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,2);
Jz = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,3);

J = [Jx;Jy;Jz];

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

L = chol(H,'lower');
logdetH = 2*sum(log(diag(L)));

cost = -logdetH;
end

function cost = compute_cost_function_a_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A)

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
Jx = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,1);
Jy = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,2);
Jz = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,3);

J = [Jx;Jy;Jz];

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% L = chol(H,'lower');
% cost = sum(1./diag(L).^2);

cost = trace(inv(H));
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



function dphidp = compute_cost_function_gradient_d_opt(...
    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim,p)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);

% Compute the Jacobian derivative w.r.t sensor positions
dJds = compute_dJ_ds(img,dim,p);

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% Intermediate variable W (different in a-opt to d-opt)
W =  inv_Gamma_noise*J*(inv(H)); %[n_stim*n_sensors,n_elem]

% Compute the derivative of the cost w.r.t. the p-coordinate of
% the sensor
dphidp = zeros(1,n_sensors);
for m = 1:n_sensors
    temp = 0.0;
    for l = 1:n_stim
        id = m+(l-1)*n_sensors;
        temp = temp -2*sum(W(id,:) .* squeeze(dJds(m,l,:)).'); % Frobenius inner product
    end
    dphidp(m) = temp;
end

end

function dphidp = compute_cost_function_gradient_a_opt(...
    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim,p)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);

% Compute the Jacobian derivative w.r.t sensor positions
dJds = compute_dJ_ds(img,dim,p);

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% Intermediate variable W (different in a-opt to d-opt)
W =  inv_Gamma_noise*J*(inv(H)^2); %[n_stim*n_sensors,n_elem]

% Compute the derivative of the cost w.r.t. the p-coordinate of
% the sensor
dphidp = zeros(1,n_sensors);
for m = 1:n_sensors
    temp = 0.0;
    for l = 1:n_stim
        id = m+(l-1)*n_sensors;
        temp = temp -2*sum(W(id,:) .* squeeze(dJds(m,l,:)).'); % Frobenius inner product
    end
    dphidp(m) = temp;
end

end


% For 1-axis MDEIT
function dphidp = compute_cost_function_gradient_d_opt_optimized(...
    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);
dJds = compute_dJ_xyz(img,dim);

H = J.'*inv_Gamma_noise*J + inv_Gamma_prior;

L = chol(H,'lower');

% Compute H^{-1} through Cholesky factorization and linear system solves
% H = L*L' -> H^{-1} = L'^{-1}*L^{-1}: solve L first, then L'
% (previous order L'\eye then L\... computed inv(L*L') = inv(L'*L), wrong)
X = L\eye(n_elem);
invH = L'\X;

W = inv_Gamma_noise*J*invH;

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
    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);
dJds = compute_dJ_xyz(img,dim);

H = J.'*inv_Gamma_noise*J + inv_Gamma_prior;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves
% H = L*L' -> H^{-1} = L'^{-1}*L^{-1}: solve L first, then L' (twice)
% (previous order L'\eye first computed inv(L*L') = inv(L'*L), wrong)
X = L\eye(n_elem);
invH = L'\X;
Z = L\invH;
inv_H2 = L'\Z;

W = inv_Gamma_noise*J*inv_H2;

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
    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

Jx = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,1);
Jy = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,2);
Jz = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,3);

% dJxds2 = compute_dJ_xyz(img,1);
dJ = compute_dJxyz_xyz(img);

dJxds = cell(1,3);
dJyds = cell(1,3);
dJzds = cell(1,3);

for p = 1:3
    dJxds{p} = dJ{1,p};
    dJyds{p} = dJ{2,p};
    dJzds{p} = dJ{3,p};
end

J = [Jx;Jy;Jz];

H = J.'*inv_Gamma_noise*J + inv_Gamma_prior;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves
% H = L*L' -> H^{-1} = L'^{-1}*L^{-1}: solve L first, then L' (twice)
% (previous order L'\eye first computed inv(L*L') = inv(L'*L), wrong)
X = L\eye(n_elem);
invH = L'\X;
Z = L\invH;
inv_H2 = L'\Z;

W = inv_Gamma_noise*J*inv_H2;

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


function dphidp = compute_cost_function_gradient_d_opt_optimized_3_axis(...
    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

Jx = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,1);
Jy = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,2);
Jz = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,3);

% dJxds2 = compute_dJ_xyz(img,1);
dJ = compute_dJxyz_xyz(img);

dJxds = cell(1,3);
dJyds = cell(1,3);
dJzds = cell(1,3);

for p = 1:3
    dJxds{p} = dJ{1,p};
    dJyds{p} = dJ{2,p};
    dJzds{p} = dJ{3,p};
end

J = [Jx;Jy;Jz];

H = J.'*inv_Gamma_noise*J + inv_Gamma_prior;

L = chol(H,'lower');

% Compute H^{-1} through Cholesky factorization and linear system solves
% H = L*L' -> H^{-1} = L'^{-1}*L^{-1}: solve L first, then L'
% (previous order L'\eye then L\... computed inv(L*L') = inv(L'*L), wrong)
X = L\eye(n_elem);
invH = L'\X;

W = inv_Gamma_noise*J*invH;

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

function sensor_locations = fetch_sensor_locations(img)
n_sensors = numel(img.fwd_model.sensors);
sensor_locations = zeros(n_sensors,3);

for m = 1: numel(img.fwd_model.sensors)
    sensor_locations(m,:) = img.fwd_model.sensors(m).position;
end

end
