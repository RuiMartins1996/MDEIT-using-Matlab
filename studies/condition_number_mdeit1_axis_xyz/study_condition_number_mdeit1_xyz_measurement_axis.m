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


%% Setup EIDORS
clc;

rng(1)

data_folder = strcat(script_folder ,'\data');

%% Build/load multiple forward models of different sensor radius

model_parameters.maxsz = 0.5;
model_parameters.numOfRings = 4;
model_parameters.numOfElectrodesPerRing = 4;
model_parameters.numOfSensors = 16;
model_parameters.isCylindrical = true;

% Create file name
file_name = create_file_name(data_folder,model_parameters);

r_min = 1.01;
r_max = 3;
n_steps = 10;

% Sweep model_parameters over multiple cylindrical radius
model_parameters_array = ...
    sweep_model_parameters({'sensorRadius'},r_min,r_max,n_steps,model_parameters);

% Build models
[model_parameters_array,fmdl_array] = mk_mdeit_model(model_parameters_array,model_folder,[]);

s.r_min = r_min;
s.r_max = r_max;
s.n_steps = n_steps;
s.model_parameters = model_parameters;
%% Set anonymous functions needed by calc_jacobian_mdeit
imgh = mk_image_mdeit(fmdl_array{1},1.0);

lambdatimesdAdp = @(lambda) computeLambdaTimesDaDp(imgh,lambda);
A = @(sigma) M(imgh,sigma);

%% Compute the jacobian for each of them

condition_number_array = nan(length(fmdl_array),3);

for d = 1:3
    select_sensor_axis = d;

    for i = 1:length(fmdl_array)
        fprintf('Model %i of %i\n',i,length(fmdl_array));

        % Make homogeneous image
        imgh = mk_image_mdeit(fmdl_array{i},1.0);

        J = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',select_sensor_axis);

        singular_values = svds(J,rank(J));

        condition_number_array(i,d) = singular_values(1)/singular_values(end);
        
        s.condition_number_array = condition_number_array;
        save(file_name,"s");
    end
end

%% FUNCTIONS
function file_name = create_file_name(data_folder,model_parameters)
    name = sprintf('data_E_%i_R_%i_M_%i',...
        model_parameters.numOfElectrodesPerRing,...
        model_parameters.numOfRings,...
        model_parameters.numOfSensors);
    file_name = strcat(data_folder,"\",name);
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