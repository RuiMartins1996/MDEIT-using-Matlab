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
grandparent_folder =fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));

model_folder = prepare_workspace(script_folder);

rng(1);

% Make sure your parallel pool is started
pool = gcp('nocreate');  % get current pool, do not create
if isempty(pool)
    pool = parpool('local');  % only create if none exists
end
% Run EIDORS startup on all workers
pctRunOnAll(sprintf('setupEidors(''%s'');', script_folder));


colors = [228,26,28;
55,126,184;
77,175,74;
152,78,163;
255,127,0;
255,255,51;
166,86,40;
247,129,191]/255;

markers = {'+','o','d','s','p','x','^' , '>' , '<' };
%% Model parameters 
z_0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z_0*I0/(l0^2); %(V)
sigma0 = l0/z_0; %(S/m)
J0 = I0/(l0^2);

model_parameters = create_default_3d_model_parameters(l0, z_0, sigma0, I0);

model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/8;
model_parameters.numOfElectrodesPerRing = 3;
model_parameters.numOfRings = 1;
model_parameters.numOfSensors = 2;
model_parameters.material.type = 'spherical';
model_parameters.material.position(1) = 0.95*(model_parameters.radius-model_parameters.material.radius);
model_parameters.material.position(3) = 0.5*model_parameters.height;

background_conductivity = 3.28e-1/sigma0;
anomaly_conductivity = background_conductivity*1.1;

current_amplitude = 2.4e-3/I0;

inj = [0 1]; %skip 0 pattern (pg 172)
meas = [0 1]; 
%% Simulation parameters
dim = 3;            %if doing 1-axis, which dimmension?

rmax = 2*model_parameters.radius;
rmin = model_parameters.radius*1.01;

R0 = model_parameters.radius*1.5;

snr = 5;
num_noise_repetitions = 30;
%% Make forward model with material

% Forward model with material 
model_parameters.maxsz = model_parameters.radius/3;
[~,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl = fmdls{1};

stimulation = ...
    mk_stim_patterns(...
    model_parameters.numOfElectrodesPerRing,...
    model_parameters.numOfRings,...
    inj,meas,{'meas_current'},current_amplitude);

fmdl.stimulation = stimulation;


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

n_sensors = numel(imgh_recon.fwd_model.sensors);
n_stim = numel(imgh_recon.fwd_model.stimulation);
n_elem = size(imgh_recon.fwd_model.elems,1);
n_nodes = size(fmdl.nodes,1);

A = @(x) M(imgh_recon,x);


%% Initialize
num_of_sensors_per_ring = n_sensors/model_parameters.numOfRings;

r_0 = R0*ones(1,model_parameters.numOfRings).*ones(num_of_sensors_per_ring,1);
r_0 = r_0(:);

dtheta = 2*pi/(num_of_sensors_per_ring);
theta_ring = 0:dtheta:2*pi-dtheta;
theta_0 = ones(1,model_parameters.numOfRings).*theta_ring.';
theta_0 = theta_0(:);

fmdl_height = max(fmdl_recon.nodes(:,3))-min(fmdl_recon.nodes(:,3));
dh = fmdl_height/(model_parameters.numOfRings+1);
sensor_heights = dh:dh:fmdl_height-dh;

z_0 = sensor_heights.*ones(num_of_sensors_per_ring,1);
z_0 = z_0(:);

% Same number of rings as electrodes
sensor_locations_0 = [(r_0(:).*cos(theta_0(:))),(r_0(:).*sin(theta_0(:))),z_0(:)];

imgh_recon = assign_sensor_locations(imgh_recon,sensor_locations_0);

show_fem(imgh_recon)
plot_sensors(imgh_recon);

%% Define prior and noise covariance matrices

Bi = fwd_solve_mdeit(imgi);
Bh = fwd_solve_mdeit(imgh);
dB = [Bi.Bx(:)-Bh.Bx(:);Bi.By(:)-Bh.By(:);Bi.Bz(:)-Bh.Bz(:)];
max_B = max(dB(:));

% Set the noise variance with respect to the data magnitude
noise_std_deviation = max_B*1e-1;
variance_noise = noise_std_deviation^2;

Gamma_noise_3_axis = variance_noise.*speye(3*n_stim*n_sensors,3*n_stim*n_sensors);
inv_Gamma_noise_3_axis = inv(Gamma_noise_3_axis);

inv_Gamma_noise = inv_Gamma_noise_3_axis;

%% Determine the prior variance
% such that lambda = sigma_epsilon^2/sigma_p^2 is the optimal tykhonov regularization 
% paramter according to L-Curve method for the initial sensor configuration

imgtemp = imgh_recon;
imgtemp = assign_sensor_locations(imgtemp,sensor_locations_0);

[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtemp,A);
J = [Jx;Jy;Jz];

n_elem = size(J,2);

% Add noise to measurements
dy = dB + sqrt(variance_noise)*randn(size(dB,1),1);

lambda_vector = logspace(-20,1,50);
[lambda_opt,dx] = l_curve_method_new(J,dy,lambda_vector);

variance_prior = variance_noise/lambda_opt;

%% Diagonal prior covariance matrix

% imgtemp = imgh_recon;
% imgtemp = assign_sensor_locations(imgtemp,sensor_locations_0);
% 
% [Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgtemp,A);
% J = [Jx;Jy;Jz];
% variance_prior = variance_noise/(trace(J'*J)/size(J,2));

% THIS CONSTANT IS A REGULARIZATION PARAMETER. HOW SHOULD WE PICK IT???????
% variance_prior = 1; 
%this is data dominated regime, since Gamma_noise = sigma_noise*I, Gamma_prior = sigma_prior*I, 
% then Gamma_post = (1/sigma_noise J'*J + 1/sigma_prior*I)^-1 =
% sigma_noise*(J'J sigma_noise/sigma_prior I)^-1, then
% sigma_noise/sigma_prior is small compared to eigenvalues of J'J

% MAYBE WE CAN PICK THIS REGULARIZATION TERM WITH L-CURVE SINCE THIS IS
% TIKHONOV!!!!!!!!!

% Gamma_prior = variance_prior.*speye(n_elem,n_elem);
% inv_Gamma_prior = inv(Gamma_prior);
% inv_Gamma_prior_3_axis = inv_Gamma_prior;

fprintf('Building smooth prior\n');

str = sprintf("(x-%2.2f).^2+(y-%2.2f).^2+(z-%2.2f).^2<%2.2f^2",...
    model_parameters.anomaly.position(1),model_parameters.anomaly.position(2),model_parameters.anomaly.position(3),model_parameters.anomaly.radius);

select_fcn = inline(str,'x','y','z');
idx = elem_select(imgh_recon.fwd_model, select_fcn);
idx(idx>0) = 1;
in_block  = idx * idx.';  % both in
out_block = (~idx) * (~idx).';% both out

% Assemble kappa
kappa = (0.4 * sparse(in_block) + 0.03 * sparse(out_block));

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
fprintf('\t Inverting smooth prior (might take a while)\n');
Gamma_prior = variance_prior*Gamma_smooth;
Gamma_prior = Gamma_prior + 1e-6 * speye(size(Gamma_prior)); %need regularization because of ill-conditioning
inv_Gamma_prior = Gamma_prior \ speye(size(Gamma_prior));
inv_Gamma_prior_3_axis = inv_Gamma_prior;

img_prior = imgh_recon;
img_prior.elem_data = diag(Gamma_prior);
show_fem_transparent_edges(img_prior);
drawnow;

%% Define jacobian of coordinate transformations

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


jac_coord_transf_xi = @(s)  compute_jacobian_coordinate_transformation_xi(s,rmin,rmax);

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
function x = map_q_to_x_xi(q,rmin,rmax,theta_0,z_0)
n_sensors = numel(q);

xim = q(1:n_sensors);
thetam = theta_0;
zm = z_0;

sigmoid = @(x) 1./(1+exp(-x));

rm = rmin + (rmax-rmin).*sigmoid(xim);

x = [rm(:).*cos(thetam(:));rm(:).*sin(thetam(:));zm(:)];
end

function q = map_x_to_q_xi(x,rmin,rmax)

assert(mod(numel(x),3)==0);
n_sensors = numel(x)/3;

xm = x(1:n_sensors);
ym = x(n_sensors+1:2*n_sensors);

rm = sqrt(xm.^2+ym.^2);
xim = log((rm-rmin)./(rmax-rm));

q = [xim(:)];
end

x_to_q_xi = @(x) map_x_to_q_xi(x,rmin,rmax);
q_to_x_xi = @(q) map_q_to_x_xi(q,rmin,rmax,theta_0,z_0);

%% Define vector_to_sensor_locations and vice-versa

x0 = sensor_locations_to_vector_cartesian(sensor_locations_0);

% Full map from q coordinates to sensor locations and back for
% unconstrained optimization
vector_to_sensor_locations_xi = @(q) vector_to_sensor_locations_cartesian(q_to_x_xi(q));
sensor_locations_to_vector_xi = @(sensor_locations) x_to_q_xi(sensor_locations_to_vector_cartesian(sensor_locations));

assert(all(all(sensor_locations_0 == vector_to_sensor_locations_xi(sensor_locations_to_vector_xi(sensor_locations_0)))));
%% Run optimization

max_iterations = 100;

options = optimoptions('fminunc',...
    'Display','iter','MaxIterations',max_iterations,...
    'StepTolerance',1e-9,'OptimalityTolerance',1e-9,...
    'Algorithm','quasi-newton','HessianApproximation','lbfgs',...
    'SpecifyObjectiveGradient',true,'UseParallel',true);


% Launch optimization from different initial conditions
% Use a grid based multi-start

num_of_sensors = size(z_0,1);
num_initial_conditions = 1;

imgsh = cell(num_initial_conditions,1);
imgout = cell(num_initial_conditions,1);


%% DEBUG
r_rand = rmin+(rmax-rmin)*rand(2,1);
sensor_locations_0 = [(r_rand(:).*cos(theta_0(:))),(r_rand(:).*sin(theta_0(:))),z_0(:)];

imgtemp = imgh_recon;
imgtemp  = assign_sensor_locations(imgtemp ,sensor_locations_0);

Ln = chol(Gamma_noise_3_axis, 'lower');
Lp = chol(Gamma_prior, 'lower');

c0 = compute_cost_function_a_opt_3_axis(imgtemp,sensor_locations_0,inv_Gamma_prior,inv_Gamma_noise_3_axis,A);

c1 =  compute_cost_function_a_opt_3_axis_chol(imgtemp ,sensor_locations_0,inv_Gamma_prior,inv_Gamma_noise_3_axis,A);

c2 = compute_cost_function_a_opt_3_axis_no_inverses_gpt(imgtemp, sensor_locations_0, Lp, Ln, A);

c3 =  compute_cost_function_a_opt_3_axis_no_inverses_no_cholesky(...
    imgtemp, sensor_locations_0, Gamma_prior, Gamma_noise_3_axis, A);

%%
n = 1;
for k = 1:num_initial_conditions
        
        r_rand = rmin+(rmax-rmin)*rand(2,1);

        sensor_locations_0 = [(r_rand(:).*cos(theta_0(:))),(r_rand(:).*sin(theta_0(:))),z_0(:)];

        imgsh{n} = imgh_recon;
        imgsh{n} = assign_sensor_locations(imgsh{n},sensor_locations_0);

        img_out = optimize_sensor_configuration(imgsh{n},inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,...
            jac_coord_transf_xi,q_to_x_xi,x_to_q_xi,'a-opt','mdeit3',3,options);

        imgout{n} = img_out;
        n=n+1;
end

figure
for n = 1:num_initial_conditions
    subplot(1,num_initial_conditions,n)
    show_fem_transparent_edges(imgout{n});
    plot_sensors(imgout{n});
end
drawnow

%% Optimize sensor positions with brute force
function [sensor_locations_brute_force,min_cost,R1,R2,C] = ...
    optimize_sensor_positions_brute_force(dq,imgh_recon,inv_Gamma_prior,inv_Gamma_noise,A,rmin,rmax,theta_0,z_0,N)

num_of_sensors = numel(z_0);

r = linspace(rmin,rmax,N);

index_combinations = nchoosek(1:N, num_of_sensors);
% Add pairs of (r1,r2) where r1 = r2;
index_combinations = [index_combinations; [(1:N).',(1:N).']];

num_of_combinations = size(index_combinations,1);

cost  = zeros(num_of_combinations,1);

afterEach(dq, @(~) updateProgress(num_of_combinations));

% Function that each worker uses to output progress
    function updateProgress(num_of_combinations, reset)
        persistent completed t_global
        if nargin > 1 && reset
            completed = 0;
            t_global = tic;
            return
        end

        if isempty(completed)
            completed = 0;
            t_global = tic;
        end

        completed = completed + 1;
        elapsed = toc(t_global);

        avg_time_per_worker = elapsed / completed;

        eta = avg_time_per_worker * (num_of_combinations - completed);
        fprintf('ETA (s): %d\n', ceil(eta));
    end

updateProgress(0, true);

%% Loop through all the sensor positions and compute cost
parfor k = 1:num_of_combinations
    
    sensor_locations_brute_force = zeros(num_of_sensors,3);
    for m = 1:num_of_sensors
        rm = r(index_combinations(k,m));
        sensor_locations_brute_force(m,:) = [rm.*cos(theta_0(m)),rm.*sin(theta_0(m)),z_0(m)];
    end

    cost(k) = ...
        compute_cost_function_a_opt_3_axis(imgh_recon,sensor_locations_brute_force,inv_Gamma_prior,inv_Gamma_noise,A);

    send(dq, k);
end

fprintf('Done\n');

%% Find minimal cost 


%% DEBUG ( mesh the cost function w.r.t. theta1, theta2)

[R1,R2] = meshgrid(r,r);
C = zeros(size(R1));

for i = 1:size(index_combinations,1)
    ids = index_combinations(i,:);
    
    % Symmetry of changing sensor 1 with sensor 2
    C(ids(1),ids(2)) = cost(i);
    C(ids(2),ids(1)) = cost(i);
end


%%

[min_cost, linear_idx] = min(cost(:));  % minimum value and its linear index

r_min = zeros(num_of_sensors,1);

for m = 1:num_of_sensors
    r_min(m) = r(index_combinations(linear_idx,m));
end

fprintf('Minimal cost: %2.2f\n', min_cost);
fprintf('Radius: \n \t t1=%.4f\n', r_min);

sensor_locations_brute_force = zeros(num_of_sensors,3);
for m = 1:num_of_sensors
    rm = r_min(m);
    sensor_locations_brute_force(m,:) = [rm.*cos(theta_0(m)),rm.*sin(theta_0(m)),z_0(m)];
end

end

%% Do brute force for several values of N

% --- Load existing results if available ---
if isfile('data/sensor_optimization_results.mat')
    load('data/sensor_optimization_results.mat','resultsTable');
else
    resultsTable = table([], [], [], 'VariableNames', {'N','SensorLocations','MinCost'});
end

N_vec = [20];

sensor_locations_brute_force = cell(1,numel(N_vec));
min_cost = zeros(1,numel(N_vec));

for n = 1:numel(N_vec)
    N = N_vec(n);

    % Check if N is already in the saved table
    if any(resultsTable.N == N)
        fprintf('Skipping N=%d, already computed.\n', N);
        idx = find(resultsTable.N == N, 1);
        sensor_locations_brute_force{n} = resultsTable.SensorLocations{idx};
        min_cost(n) = resultsTable.MinCost(idx);
        continue
    end

    dq = parallel.pool.DataQueue;
    
    [sensor_locations_brute_force{n},min_cost(n),R1,R2,C] = ...
        optimize_sensor_positions_brute_force(dq,imgh_recon,inv_Gamma_prior,inv_Gamma_noise,A,rmin,rmax,theta_0,z_0,N);
    
    %% DEBUG

    figure
    hold on
    surf(R1,R2,C);
    shading interp
    colormap hot

    xlabel('$\theta_1$','Interpreter','latex')
    ylabel('$\theta_2$','Interpreter','latex')
    zlabel('$C(\theta_1,\theta_2)$','Interpreter','latex')
    title('Cost function')

    % Get all local minima
    Cmin = imregionalmin(C);

    % Get indices of minima
    [row, col] = find(Cmin);

    % Map to actual coordinates
    r1_min = R1(sub2ind(size(R1), row, col));
    r2_min = R2(sub2ind(size(R2), row, col));
    C_values   = C(sub2ind(size(C), row, col));

    % Local minima
    for i = 1:numel(r1_min)
        plot3(r1_min(i),r2_min(i),C_values(i),'r.','MarkerSize',5)
    end
    % Global minimum
    ids = find(C_values == min(C_values));
    plot3(r1_min(ids),r2_min(ids),C_values(ids),'cs','MarkerSize',10)
    
    min_cost_sopt = inf;
    min_cost_id = inf;

    for id = 1:num_initial_conditions
        s0 = fetch_sensor_locations(imgsh{id});
        % cost_s0 = ...
        %     compute_cost_function_d_opt_3_axis(imgsh{id},fetch_sensor_locations(imgsh{id}),inv_Gamma_prior,inv_Gamma_noise,A);
        % 
        cost_s0 = ...
            compute_cost_function_a_opt_3_axis(imgsh{id},fetch_sensor_locations(imgsh{id}),inv_Gamma_prior,inv_Gamma_noise,A);
        
        r0 = sqrt(s0(:,1).^2+s0(:,2).^2);
        plot3(r0(1),r0(2),cost_s0,'k.','MarkerSize',5)

        sopt = fetch_sensor_locations(imgout{id});
        % cost_sopt = ...
        %     compute_cost_function_d_opt_3_axis(imgout{id},sopt,inv_Gamma_prior,inv_Gamma_noise,A);

        cost_sopt = ...
            compute_cost_function_a_opt_3_axis(imgout{id},sopt,inv_Gamma_prior,inv_Gamma_noise,A);
        
        if cost_sopt<min_cost_sopt
            min_cost_id = id;
            min_cost_sopt = cost_sopt;
        end

        r_opt = sqrt(sopt(:,2).^2+sopt(:,1).^2);
        
        plot3(r_opt(1),r_opt(2),cost_sopt,'y.','MarkerSize',15)
    end
    
    r_bf = sqrt(...
        sensor_locations_brute_force{n}(:,2).^2+...
        sensor_locations_brute_force{n}(:,1).^2);

    plot3(r_bf(1),r_bf(2),min_cost(n),'rh','MarkerSize',10)


    axis([rmin,rmax,rmin,rmax])
    
    legend('Cost function','Initial configuration','Optimal configuration')
    hold off

    %% --- Append new result to the table ---
    newRow = {N, sensor_locations_brute_force{n}, min_cost(n)};
    resultsTable = [resultsTable; newRow];
    save('data/sensor_optimization_results.mat', 'resultsTable');
end

%% Show brute-force and optimization solutions side by side
figure

function h = plot_sensors_local(ax,img,color,marker)

if nargin<3
    color = 'b';
end
if nargin <4
    marker = '.';
end

sensor_positions = zeros(numel(img.fwd_model.sensors),3);

for m = 1:numel(img.fwd_model.sensors)
    sensor_positions(m,:) = img.fwd_model.sensors(m).position;
end

h = plot3(ax,sensor_positions(:,1),sensor_positions(:,2),sensor_positions(:,3),'LineStyle','none','Color',color,'Marker',marker);

end

ax = axes;
hold(ax, 'on')
all_handles = [];

ct = 1;
sopt = fetch_sensor_locations(imgout{min_cost_id});
% min_cost_opt = ...
%     compute_cost_function_d_opt_3_axis(imgout{n},sopt,inv_Gamma_prior,inv_Gamma_noise,A);

min_cost_opt = ...
    compute_cost_function_a_opt_3_axis(imgout{min_cost_id},sopt,inv_Gamma_prior,inv_Gamma_noise,A);


h_opt = plot_sensors_local(ax,imgout{min_cost_id},colors(n,:),'h');
all_handles = [all_handles h_opt];
all_msgs{ct} = sprintf('BFGS - cost = %2.2f',min_cost_opt);
ct = ct+1;

if isfile('data/sensor_optimization_results.mat')
    load('data/sensor_optimization_results.mat','resultsTable');
else
    error('No results found');
end

[sortedN,ids] = sort(resultsTable.N);

for n = 1:numel(resultsTable.N)
    imgtemp = assign_sensor_locations(imgh_recon,resultsTable.SensorLocations{ids(n)});
    
    if mod(n,numel(markers)) == 0
        marker = markers{1};
    else
        marker = markers{mod(n,numel(markers))};
    end

    if mod(n,size(colors,1)) == 0
        color = colors(1,:);
    else
        color = colors(mod(n,size(colors,1)),:);
    end
    
    h_brute = plot_sensors_local(ax,imgtemp,color,marker);
    all_handles = [all_handles h_brute];
    all_msgs{ct} = sprintf('BF| N = %i |cost = %2.3f',resultsTable.N(ids(n)),resultsTable.MinCost(ids(n)));
    ct = ct+1;
end

grid on;
grid minor;
box on;

show_fem_transparent_edges(img_prior)
axis([-rmax rmax -rmax rmax 0 3])

legend(all_handles,all_msgs,...
    'Location','southeast','interpreter','latex');


diff_sensors_locations = zeros(1,numel(sensor_locations_brute_force));

for n = 1:size(resultsTable,1)
    diff_sensors_locations(n) = norm(sopt-resultsTable.SensorLocations{ids(n)},'fro');
    fprintf('||s_1-s_2||_f: %2.2f\n',diff_sensors_locations(n));
end

figure

plot(sortedN,diff_sensors_locations,'k.-')
xlabel('N','Interpreter','latex')
ylabel('$||s_{BFGS}-s_{BF}||_{F}$','Interpreter','latex')
set(gca,'YScale','log')
grid on;
grid minor;
box on;

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

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);

J = [Jx;Jy;Jz];

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% L = chol(H,'lower');
% cost = sum(1./diag(L).^2);

cost = trace(inv(H));
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


function data_noisy = add_measurement_noise_difference_absolute(datai,datah,noise_amplitude)

if isempty(noise_amplitude)
    data_noisy = datai-datah;
    return;
else
    d = datai-datah;

    data_noisy = d+noise_amplitude*randn(size(d));
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



%% Functions
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

function cost = compute_cost_function_a_opt_3_axis_chol(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A)

n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
[Jx,Jy,Jz]= calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);
J = [Jx;Jy;Jz];

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% CHANGED HERE!
% cost = trace(inv(H));

L = chol(H,'lower');
Hinv = L'\(L\eye(n_elem));
cost = trace(Hinv);

end


function cost = compute_cost_function_a_opt_3_axis_no_inverses_gpt( ...
    img, sensor_locations, Lp, Ln, A)

n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img, sensor_locations);

% Compute Jacobian
[Jx, Jy, Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img, A);
J = [Jx; Jy; Jz];

% Whiten J (no inverse formed)
Jw   = Ln \ J;        % solves Ln * Jw = J
Jhat = Jw * Lp;       % consistent with Lp Lp' = Gamma_prior

% Build and factor Hhat
Hhat = Jhat' * Jhat + speye(n_elem);
Lh   = chol(Hhat, 'lower');

% Compute trace(H^{-1})
Z = Lh \ Lp';         % solves Lh * Z = Lp'
cost = sum(Z(:).^2); % Frobenius norm squared

end


function cost =  compute_cost_function_a_opt_3_axis_no_inverses_no_cholesky(...
    img, sensor_locations, Gamma_prior, Gamma_noise, A)

n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img, sensor_locations);

% Compute Jacobian
[Jx, Jy, Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img, A);
J = [Jx; Jy; Jz];

% Precompute
GJ = J * Gamma_prior;
M  = GJ * J' + Gamma_noise;

% trace term 1
t1 = trace(Gamma_prior);

% Solve M X = J*Gamma_prior 
X = M \ (J*Gamma_prior);

t2 = trace(Gamma_prior*J.'*X);

cost = t1-t2;

end

