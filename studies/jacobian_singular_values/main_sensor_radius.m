clc;clear all;close all;

%% Prepare workspace
% Get the full path of the current script
fullpath = mfilename('fullpath');
% Extract just the folder
script_folder = fileparts(fullpath);
cd(script_folder);

% Have to add the functions path manually so prepare_workspace runs
parent_folder = fileparts(script_folder);
grandparent_folder =fileparts(parent_folder);
addpath(genpath(fullfile(grandparent_folder,'functions')));

model_folder = prepare_workspace(script_folder);

data_folder = '.\data_sensor_radius';

%% Define test parameters

min_sensor_radius = 1.01;
max_sensor_radius = 1.2;
n_steps = 50;

sweep_mode = 'linear';

%Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);

maxsz_reconstruction = 4.5e-3/l0;
background_conductivity =  3.28e-1/sigma0; 

%Stimulation pattern was not the default. For now, edit it manually
current_amplitude = 2.4e-3/I0;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used

%% Create template model_parameters to prepare for sweep
model_parameters = create_kai_3d_model_parameters(l0, z0, sigma0, I0);

% model_parameters.maxsz = 4.5e-3/l0;
% model_parameters.numOfRings = 3;
% model_parameters.radius = 40e-3/l0;

model_parameters.maxsz = 2*model_parameters.maxsz;
model_parameters.radius = 80e-3/l0;
model_parameters.height = model_parameters.radius;
model_parameters.numOfSensors = model_parameters.numOfRings*model_parameters.numOfElectrodesPerRing;


% Get default parameters name so we can test other default cases
model_name = model_parameters_to_file_name(model_parameters,model_folder);

%model_name has the fullfile format into the forlder models. Extract only
%the name of the file and forget about the directory

model_name = regexp(model_name,'([^\\\/]+)(?=\.[^.]+$)','match');
model_name = model_name{1};

% Set file names according to model name
if ~exist(fullfile(script_folder,data_folder,model_name),'dir')
    mkdir(fullfile(script_folder,data_folder,model_name));
end


file_name_eit = fullfile(data_folder,model_name,'/singular_values_eit.mat');
file_name_mdeit_x = fullfile(script_folder,data_folder,model_name,'/singular_values_mdeit_x.mat');
file_name_mdeit_y = fullfile(script_folder,data_folder,model_name,'/singular_values_mdeit_y.mat');
file_name_mdeit_z = fullfile(script_folder,data_folder,model_name,'/singular_values_mdeit_z.mat');
file_name_mdeit_3 = fullfile(script_folder,data_folder,model_name,'/singular_values_mdeit_3.mat');

file_name_model_parameters = fullfile(data_folder,model_name,'/model_parameters.mat');

save(file_name_model_parameters,"model_parameters");


%% Sweep model_parameters over magnetometer ring radius
model_parameters_array = ...
    sweep_model_parameters({'sensorRadius'},min_sensor_radius,max_sensor_radius,n_steps,model_parameters,sweep_mode);


%% Create models for forward simulation and assign stimulation pattern (no material here, homogeneous conductivity)

[model_parameters,fmdls] = mk_mdeit_model(model_parameters_array,model_folder);

for i = 1:length(fmdls)

    num_of_rings = model_parameters(i).numOfRings;
    
    stimulation = ...
        mk_stim_patterns(numel(fmdls{i}.electrode)/num_of_rings,num_of_rings,inj,meas,{'meas_current'},current_amplitude);

    fmdls{i}.stimulation = stimulation;
end

%% Compute singular values of Jacobian for all the models

num_of_elements_array = zeros(n_steps,1);
random_seed_array = zeros(n_steps,1);

num_of_measurements_eit_array = zeros(n_steps,1);
num_of_measurements_mdeit_array = zeros(n_steps,1);

rank_eit_array = zeros(n_steps,1);
rank_mdeit_x_array = zeros(n_steps,1);
rank_mdeit_y_array = zeros(n_steps,1);
rank_mdeit_z_array = zeros(n_steps,1);
rank_mdeit_3_array = zeros(n_steps,1);

% NUmber of measurements is constant when varying sensor radius
num_of_sensors = model_parameters_array(1).numOfSensors;
num_of_injections = numel(stimulation);
num_of_measurements_eit_per_pattern = size(stimulation(1).meas_pattern,1);

num_of_measurements_eit = num_of_measurements_eit_per_pattern*num_of_injections;
num_of_measurements_mdeit = num_of_sensors*num_of_injections;

% Arrays for storing singular values
s_eit_array = zeros(num_of_measurements_eit,n_steps);
s_mdeit_x_array = zeros(num_of_measurements_mdeit,n_steps);
s_mdeit_y_array = zeros(num_of_measurements_mdeit,n_steps);
s_mdeit_z_array = zeros(num_of_measurements_mdeit,n_steps);
s_mdeit_3_array = zeros(3*num_of_measurements_mdeit,n_steps);

       
function r = compute_rank(s,J)
% tolerance for computing rank according to MATLAB's rank function
tol = max(size(J)) * eps(norm(J));
r = sum(s > tol);
end

data_already_exists = false;
if exist(fullfile(data_folder,model_name,'model_parameters.mat'),'file')
    model_file = load(fullfile(data_folder,model_name,'model_parameters.mat'));
    model_parameters_loaded = model_file.model_parameters;
    data_already_exists = true;
end

for i = 1:length(fmdls)
    
    fprintf('Running model %i of %i\n',i,length(fmdls));
    
    fmdl = fmdls{i};

    num_of_elements_array(i) = size(fmdl.elems,1);
    random_seed_array(i) = model_parameters_array(i).randomConductivitySeed;
    
    % Sanity check ( see if the number of eit measurements is the same as
    % 1-axis mdeit measurements)
    assert(num_of_measurements_eit ==  num_of_measurements_mdeit, 'Expected these values to be the same');

    num_of_measurements_eit_array(i) = num_of_measurements_eit;
    num_of_measurements_mdeit_array(i) = num_of_measurements_mdeit;
    sensor_radius = sqrt(fmdl.sensors(1).position(1)^2+fmdl.sensors(1).position(2)^2);

    % Check if this entry is already on the data matrices
    % id1 = find(num_of_elements_array(i)== data_file.num_of_elements);
    % id2 = find(num_of_measurements_eit_array(i)== data_file.num_of_measurements);
    % id3 = find(random_seed_array(i) == data_file.random_seed);
    % id4 = find(sensor_radius-data_file.sensor_radius<1e-6);
    
    % TODO: Don't repeat calculations of data that already exists. Must
    % think of a smart way to detect that
    is_match = false;

    if ~is_match
    
        tic;
    % Make homogeneous image
    imgh = mk_image_mdeit(fmdl,background_conductivity);
    
    % Compute jacobians
    J_eit = calc_jacobian(imgh);

    lambdatimesdAdp = @(lambda) computeLambdaTimesDaDp(imgh,lambda);
    A = @(sigma) M(imgh,sigma);

    J_mdeit_x = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',1);
    J_mdeit_y = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',2);
    J_mdeit_z = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',3);

    J_mdeit_3 = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit3');
    
    % Using a custom function compute_rank avoids doing svds twice! 

    % Compute singular values
    s_eit = svds(J_eit,size(J_eit,1));
    s_mdeit_x = svds(J_mdeit_x,size(J_mdeit_x,1));
    s_mdeit_y = svds(J_mdeit_y,size(J_mdeit_y,1));
    s_mdeit_z = svds(J_mdeit_z,size(J_mdeit_z,1));
    s_mdeit_3 = svds(J_mdeit_3,size(J_mdeit_3,1));

    % Compute ranks
    rank_eit_array(i) = compute_rank(s_eit,J_eit);
    rank_mdeit_x_array(i) = compute_rank(s_mdeit_x,J_mdeit_x);
    rank_mdeit_y_array(i) = compute_rank(s_mdeit_y,J_mdeit_y);
    rank_mdeit_z_array(i) = compute_rank(s_mdeit_z,J_mdeit_z);
    rank_mdeit_3_array(i) = compute_rank(s_mdeit_3,J_mdeit_3);

    s_eit = s_eit(1:rank_eit_array(i));
    s_mdeit_x = s_mdeit_x(1:rank_mdeit_x_array(i));
    s_mdeit_y = s_mdeit_y(1:rank_mdeit_y_array(i));
    s_mdeit_z = s_mdeit_z(1:rank_mdeit_z_array(i));
    s_mdeit_3 = s_mdeit_3(1:rank_mdeit_3_array(i));
    
    % Store them in these arrays
    s_eit_array(1:rank_eit_array(i),i) = s_eit;
    s_mdeit_x_array(1:rank_mdeit_x_array(i),i) = s_mdeit_x;
    s_mdeit_y_array(1:rank_mdeit_y_array(i),i) = s_mdeit_y;
    s_mdeit_z_array(1:rank_mdeit_z_array(i),i) = s_mdeit_z;
    s_mdeit_3_array(1:rank_mdeit_3_array(i),i) = s_mdeit_3;
    
    %% Save data
    
    data_eit = struct(...
        'num_of_elements',num_of_elements_array(i),...
        'num_of_measurements',num_of_measurements_eit_array(i), ...
        'rank',rank_eit_array(i), ...
        'random_seed',model_parameters_array(i).randomConductivitySeed,...
        'sensor_radius',model_parameters_array(i).sensorRadius,...
        'singular_values',s_eit_array(:,i));

    data_mdeit_x = struct(...
        'num_of_elements',num_of_elements_array(i),...
        'num_of_measurements',num_of_measurements_mdeit_array(i), ...
        'rank',rank_mdeit_x_array(i), ...
        'random_seed',model_parameters_array(i).randomConductivitySeed,...
        'sensor_radius',model_parameters_array(i).sensorRadius,...
        'singular_values',s_mdeit_x_array(:,i));

    data_mdeit_y = struct(...
        'num_of_elements',num_of_elements_array(i),...
        'num_of_measurements',num_of_measurements_mdeit_array(i), ...
        'rank',rank_mdeit_y_array(i), ...
        'random_seed',model_parameters_array(i).randomConductivitySeed,...
        'sensor_radius',model_parameters_array(i).sensorRadius,...
        'singular_values',s_mdeit_y_array(:,i));

    data_mdeit_z = struct(...
        'num_of_elements',num_of_elements_array(i),...
        'num_of_measurements',num_of_measurements_mdeit_array(i), ...
        'rank',rank_mdeit_z_array(i), ...
        'random_seed',model_parameters_array(i).randomConductivitySeed,...
        'sensor_radius',model_parameters_array(i).sensorRadius,...
        'singular_values',s_mdeit_z_array(:,i));

    data_mdeit_3 = struct(...
        'num_of_elements',num_of_elements_array(i),...
        'num_of_measurements',3*num_of_measurements_mdeit_array(i), ...
        'rank',rank_mdeit_3_array(i), ...
        'random_seed',model_parameters_array(i).randomConductivitySeed,...
        'sensor_radius',model_parameters_array(i).sensorRadius,...
        'singular_values',s_mdeit_3_array(:,i));


    data_updated_eit = update_data(data_eit,file_name_eit);
    data_updated_mdeit_x = update_data(data_mdeit_x,file_name_mdeit_x);
    data_updated_mdeit_y = update_data(data_mdeit_y,file_name_mdeit_y);
    data_updated_mdeit_z = update_data(data_mdeit_z,file_name_mdeit_z);
    data_updated_mdeit_3 = update_data(data_mdeit_3,file_name_mdeit_3);
    
 
    save(file_name_eit,'data_updated_eit');
    save(file_name_mdeit_x,'data_updated_mdeit_x');
    save(file_name_mdeit_y,'data_updated_mdeit_y');
    save(file_name_mdeit_z,'data_updated_mdeit_z');
    save(file_name_mdeit_3,'data_updated_mdeit_3');

    t_end = toc;

    elapsed_time_per_iteration(i) = t_end;
    expected_time = (length(fmdls)-i)*mean(elapsed_time_per_iteration);
    fprintf('Expected time to completion: %2.2f (min)\n',expected_time/60);

    % all parameters are the same
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