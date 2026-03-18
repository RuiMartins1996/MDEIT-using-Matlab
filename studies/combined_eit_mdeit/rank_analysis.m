clc; clear all; close all;

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

data_folder = fullfile(script_folder,'data');


file_name_eit = strcat(data_folder,'/singular_values_eit.mat');
file_name_mdeit_x = strcat(data_folder,'/singular_values_mdeit_x.mat');
file_name_mdeit_y = strcat(data_folder,'/singular_values_mdeit_y.mat');
file_name_mdeit_z = strcat(data_folder,'/singular_values_mdeit_z.mat');
file_name_mdeit_3 = strcat(data_folder,'/singular_values_mdeit_3.mat');
file_name_aug = strcat(data_folder,'/singular_values_aug.mat');

%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);
H0 = I0/l0; 
B0 = 1.25663706127e-6*I0/l0; 

beta = V0/H0;

% Set the minimum sensor radius
sensor_radius_0 = 1.01;
rmin = sensor_radius_0;
rmax = 3;

%% Problem parameters

% Default 3D model
model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);
% with some changes
model_parameters.maxsz = model_parameters.radius/5;
model_height = model_parameters.height/2;
model_parameters.numOfElectrodesPerRing = 8;
model_parameters.numOfRings = 4;
model_parameters.numOfSensors = 8*4;
model_parameters.anomaly.type = 'spherical';
model_parameters.material.type = 'spherical';
model_parameters.material.name = 'sphere_anomaly';
model_parameters.sensorRadius = 3.5*model_parameters.radius;


% model_parameters.material.position = [0.5 0 model_parameters.height/2];
model_parameters.material.position = [0 0 model_parameters.height/2];
% model_parameters.material.radius = model_parameters.material.radius/2;

background_conductivity = 3.28e-1/sigma0;
anomaly_conductivity = 1e-12/sigma0;

maxsz_data =  model_parameters.radius/10;
maxsz_reconstruction = model_parameters.radius/7;

current_amplitude = 2.4e-3/I0;

% Set noise type and snr
noise_type = 'white';

snr = 1;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used

%% Generate coarse forward model for reconstruction (different mesh than the data)
model_parameters = rmfield(model_parameters,'material');
model_parameters.maxsz = maxsz_reconstruction;

num_of_rings_vector = [1,2,3,4];
num_of_electrodes_per_ring_vector = [4:16];

[R, E] = ndgrid(num_of_rings_vector, num_of_electrodes_per_ring_vector);
parameter_matrix = [R(:)'; E(:)'];

num_of_sensors = parameter_matrix(1,:).*parameter_matrix(2,:);
parameter_matrix(3,:) = num_of_sensors;

% Remove some triples from parameter matrix because Netgen has problems 

idx = all(parameter_matrix == [3;15;45], 1);
parameter_matrix(:, idx) = [];
% idx = all(parameter_matrix == [2;4;8], 1);
% parameter_matrix(:, idx) = [];

model_parameters_array = sweep_model_parameters({'numOfRings','numOfElectrodesPerRing','numOfSensors'},parameter_matrix,[],[],model_parameters,[]);

%% Create models
[model_parameters_array,fmdl_array] = mk_mdeit_model(model_parameters_array,model_folder);

%% Assign stimulation patterns

max_num_of_measurements_eit = -inf;

for i = 1:numel(fmdl_array)
    stimulation = ...
        mk_stim_patterns(...
        model_parameters_array(i).numOfElectrodesPerRing,...
        model_parameters_array(i).numOfRings,...
        inj,meas,{'meas_current'},current_amplitude);

    fmdl_array{i}.stimulation = stimulation;


    % Compute maximum number of measurements
    num_of_injections = numel(stimulation);
    num_of_measurements_eit_per_pattern = size(stimulation(1).meas_pattern,1);
    num_of_measurements_eit = num_of_measurements_eit_per_pattern*num_of_injections;

    max_num_of_measurements_eit = max(max_num_of_measurements_eit,num_of_measurements_eit);

end


%% Make homogeneous images
imgs = cell(1,numel(fmdl_array));

for i = 1:numel(fmdl_array)
    imgs{i} = mk_image_mdeit(fmdl_array{i},sigma0,sprintf('Model %i',i));
end


%% Compute singular values of Jacobian for all the models

num_of_elements_array = zeros(numel(model_parameters_array),1);
random_seed_array = zeros(numel(model_parameters_array),1);

num_of_measurements_eit_array = zeros(numel(model_parameters_array),1);
num_of_measurements_mdeit_array = zeros(numel(model_parameters_array),1);

rank_eit_array = zeros(numel(model_parameters_array),1);
rank_mdeit_x_array = zeros(numel(model_parameters_array),1);
rank_mdeit_y_array = zeros(numel(model_parameters_array),1);
rank_mdeit_z_array = zeros(numel(model_parameters_array),1);
rank_mdeit_3_array = zeros(numel(model_parameters_array),1);
rank_aug = zeros(numel(model_parameters_array),1);

% Arrays for storing singular values

s_eit_array = zeros(max_num_of_measurements_eit,numel(model_parameters_array));
s_mdeit_x_array = zeros(max_num_of_measurements_eit,numel(model_parameters_array));
s_mdeit_y_array = zeros(max_num_of_measurements_eit,numel(model_parameters_array));
s_mdeit_z_array = zeros(max_num_of_measurements_eit,numel(model_parameters_array));
s_mdeit_3_array = zeros(3*max_num_of_measurements_eit,numel(model_parameters_array));
s_aug_array = zeros(4*max_num_of_measurements_eit,numel(model_parameters_array));

%% Compute Jacobians


for i = 1:numel(imgs)
    fprintf('Running model %i of %i\n',i,length(fmdl_array));

    fmdl = fmdl_array{i};

    num_of_elements_array(i) = size(fmdl.elems,1);
    random_seed_array(i) = model_parameters_array(i).randomConductivitySeed;

    num_of_sensors = model_parameters_array(i).numOfSensors;
    num_of_injections = numel(fmdl.stimulation);
    num_of_measurements_eit_per_pattern = size(fmdl.stimulation(1).meas_pattern,1);

    num_of_measurements_eit = num_of_measurements_eit_per_pattern*num_of_injections;
    num_of_measurements_mdeit = num_of_sensors*num_of_injections;

    % Sanity check ( see if the number of eit measurements is the same as
    % 1-axis mdeit measurements)
    assert(num_of_measurements_eit ==  num_of_measurements_mdeit, 'Expected these values to be the same');
    
    num_of_measurements_eit_array(i) = num_of_measurements_eit;
    num_of_measurements_mdeit_array(i) = num_of_measurements_mdeit;


    % Make homogeneous image
    imgh = mk_image_mdeit(fmdl,background_conductivity);

    A = @(sigma) M(imgh,sigma);

    % Compute jacobians
    J_eit = calc_jacobian(imgh);
    J_mdeit_x = calc_jacobian_mdeit(imgh,imgh.elem_data,[],A,'mdeit1',1);
    J_mdeit_y = calc_jacobian_mdeit(imgh,imgh.elem_data,[],A,'mdeit1',2);
    J_mdeit_z = calc_jacobian_mdeit(imgh,imgh.elem_data,[],A,'mdeit1',3);
    J_mdeit_3 = calc_jacobian_mdeit(imgh,imgh.elem_data,[],A,'mdeit3');  
    J_aug = [J_eit; sqrt(beta) * J_mdeit_3]; 

    % Compute ranks
    rank_eit_array(i) = rank(J_eit);
    rank_mdeit_x_array(i) = rank(J_mdeit_x);
    rank_mdeit_y_array(i) = rank(J_mdeit_y);
    rank_mdeit_z_array(i) = rank(J_mdeit_z);
    rank_mdeit_3_array(i) = rank(J_mdeit_3);
    rank_aug(i) = rank(J_aug); 

    % Compute singular values
    s_eit = svds(J_eit,rank_eit_array(i));
    s_mdeit_x = svds(J_mdeit_x,rank_mdeit_x_array(i));
    s_mdeit_y = svds(J_mdeit_y,rank_mdeit_y_array(i));
    s_mdeit_z = svds(J_mdeit_z,rank_mdeit_z_array(i));
    s_mdeit_3 = svds(J_mdeit_3,rank_mdeit_3_array(i));
    s_aug = svds(J_aug,rank_aug(i));

    % Store them in these arrays
    s_eit_array(1:rank_eit_array(i),i) = s_eit;
    s_mdeit_x_array(1:rank_mdeit_x_array(i),i) = s_mdeit_x;
    s_mdeit_y_array(1:rank_mdeit_y_array(i),i) = s_mdeit_y;
    s_mdeit_z_array(1:rank_mdeit_z_array(i),i) = s_mdeit_z;
    s_mdeit_3_array(1:rank_mdeit_3_array(i),i) = s_mdeit_3;
    s_aug_array(1:rank_aug(i),i) = s_aug;

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

    data_aug = struct(...
        'num_of_elements',num_of_elements_array(i),...
        'num_of_measurements',4*num_of_measurements_mdeit_array(i), ...
        'rank',rank_aug(i), ...
        'random_seed',model_parameters_array(i).randomConductivitySeed,...
        'sensor_radius',model_parameters_array(i).sensorRadius,...
        'singular_values',s_aug_array(:,i));


    data_updated_eit = update_data(data_eit,file_name_eit);
    data_updated_mdeit_x = update_data(data_mdeit_x,file_name_mdeit_x);
    data_updated_mdeit_y = update_data(data_mdeit_y,file_name_mdeit_y);
    data_updated_mdeit_z = update_data(data_mdeit_z,file_name_mdeit_z);
    data_updated_mdeit_3 = update_data(data_mdeit_3,file_name_mdeit_3);
    data_updated_aug = update_data(data_aug,file_name_aug);

    save(file_name_eit,'data_updated_eit');
    save(file_name_mdeit_x,'data_updated_mdeit_x');
    save(file_name_mdeit_y,'data_updated_mdeit_y');
    save(file_name_mdeit_z,'data_updated_mdeit_z');
    save(file_name_mdeit_3,'data_updated_mdeit_3');
    save(file_name_aug,'data_updated_aug');

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
%% FUNCTIONS
function show_fem_transparent_edges(img,sensors_plot)

if nargin<2
    sensors_plot = false;
end

hh = show_fem(img);                % draw the model (hh may be a handle or array)
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);

if sensors_plot
    hold on
    plot_sensors(img);

    axis([-3 3 -3 3])
    view(0,90)
end
end