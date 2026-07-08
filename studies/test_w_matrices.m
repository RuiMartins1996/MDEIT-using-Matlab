clc; clear all; close all;

%% Prepare workspace
% Get the full path of the current script
fullpath = mfilename('fullpath');
% Extract just the folder
script_folder = fileparts(fullpath);
cd(script_folder);

% Have to add the functions path manually so prepare_workspace runs
parent_folder =fileparts(script_folder);
addpath(genpath(fullfile(parent_folder,'functions')));

model_folder = prepare_workspace(script_folder);

data_folder = fullfile(script_folder,'data');


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
model_parameters = create_kai_3d_model_parameters(l0, z0, sigma0, I0);
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

max_sizes = [model_parameters.radius/2];

maxsz_data =  model_parameters.radius/10;

current_amplitude = 2.4e-3/I0;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used

%% Generate two forward models with different number of elements

fmdl_array = cell(numel(max_sizes),1);

for i = 1:numel(max_sizes)
    model_parameters_temp = rmfield(model_parameters,'material');
    model_parameters_temp.maxsz = max_sizes(i);
    
    [~,fmdl] = mk_mdeit_model(model_parameters_temp,model_folder);
    
    fmdl_array{i} = fmdl{1};
end
%% Assign stimulation patterns

stimulation = ...
    mk_stim_patterns(...
    model_parameters.numOfElectrodesPerRing,...
    model_parameters.numOfRings,...
    inj,meas,{'meas_current'},current_amplitude);

for i = 1:numel(max_sizes)
fmdl_array{i}.stimulation = stimulation;
end

%% Test
% [Ix,Iy,Iz] = compute_w_matrices(fmdl_array{1},fmdl_array{1}.mu0,3);

img = mk_image_mdeit(fmdl_array{1},background_conductivity,'Model 1');

% [data] = fwd_solve_mdeit(img);
% [m,u] = fwd_solve_mdeit_alternative(img);

A = @(sigma) M(img,sigma);
J1 = calc_jacobian_mdeit(img,img.elem_data,[],A,'mdeit1',3);

J2 = calc_jacobian_mdeit_alternative(img,'mdeit1',3);

function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
end
