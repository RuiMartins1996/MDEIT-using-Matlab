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

data_folder = strcat(script_folder ,'\data');

clc;
rng(1)

%%
lambda_vec = logspace(-15,-5,5);


%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);

maxsz_reconstruction = 5e-3/l0;

%% Assign the parameters for several models (should create utility functions for this)

SNRdb = 100;

background_conductivity = 3.28e-1/sigma0; 

prior = 'tykhonov';
% prior = 'noser';

%% Build/load forward model

model_parameters.maxsz = 3e-3/l0;

% Tank dimensions in characteristic units
model_parameters.isCylindrical = true;
model_parameters.height=70e-3/l0;  
model_parameters.radius=40e-3/l0; 

%Electrode configuration and properties
model_parameters.numOfRings = 2;
model_parameters.numOfElectrodesPerRing = 8; 

% what's the size of the electrodes? He mentions EEG for the tank experiment, but what for the simulated experiment?
model_parameters.electrodeRadius = 8e-3/l0; %(mm) unknown from the thesis, but good estimate for Ag/AgCl EEG electrode 
assert(model_parameters.electrodeRadius*model_parameters.numOfElectrodesPerRing < 2*pi*model_parameters.radius)
model_parameters.electrodeContactImpedance = 0.0058/z0; 

% Magnetic sensor configuration and properties
model_parameters.numOfSensors = 16; %(mm) pg 172
model_parameters.sensorRadius = 70e-3/l0;  %(mm) pg 172
model_parameters.mu0 = 1; %With this definition, we're computing the magnetic field intensity H in units of I0/l0;

anomaly_conductivity = 1e-12/l0; % (0 according to Kai's thesis, but check notes)
anomaly_position = [0,0,35]*1e-3/l0; %  [20,0,35]*1e-3/l0; % pg 172
anomaly_radius = 25/2*1e-3/l0; % 25mm of diameter, pg 172

% except for the anomaly conductivity
model_parameters.anomaly = struct(...
    'type','spherical',...
    'conductivity',anomaly_conductivity,...
    'radius',anomaly_radius,...
    'position',anomaly_position);

% A material defines the geometry of that material inside the domain, not
% its properties, like conductivity.
model_parameters.material = struct(...
    'type','spherical',...
    'name','my_sphere',...
    'radius',anomaly_radius,...
    'position',anomaly_position);

% Stimulation pattern was not the default. For now, edit it manually
current_amplitude = 2.4e-3/I0;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used
options = {};

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

fmdl = fmdls{1};
%% Assign the stimulation pattern to this forward model

stimulation = mk_stim_patterns(numel(fmdl.electrode),1,inj,meas,options,current_amplitude);
fmdl.stimulation = stimulation;
%% Make images from forward models and forward solve

% Make homogeneous image
imgh = mk_image_mdeit(fmdl,background_conductivity);

% Add plastic cylinder
imgi = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

%% Generate coarse forward model for reconstruction (different mesh than the data)
model_parameters.material = struct();
model_parameters.maxsz = maxsz_reconstruction;

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_reconstruction = fmdls{1};

% Don't forget to assignt the same stimulation pattern
fmdl_reconstruction.stimulation = stimulation;

n_elem = size(fmdl_reconstruction.elems,1);
%% Run L-curve method to determine good lambda
s_mdeit1 = l_curve_method(fmdl_reconstruction,imgh,imgi,'mdeit1',lambda_vec);
s_mdeit3 = l_curve_method(fmdl_reconstruction,imgh,imgi,'mdeit3',lambda_vec);

file_name = create_file_name(data_folder,model_parameters,SNRdb);

model_struct = struct();
model_struct.imgi = imgi;
model_struct.model_parameters = model_parameters;

l_curve_struct = struct();
l_curve_struct.lambda_vec = lambda_vec;
l_curve_struct.s_mdeit1 = s_mdeit1;
l_curve_struct.s_mdeit3 = s_mdeit3;

save(file_name,'model_struct','l_curve_struct');

%% FUNCTION
function file_name = create_file_name(data_folder,model_parameters,SNRdb)
    if isempty(SNRdb)
        name = sprintf('data_E_%i_R_%i_M_%i_noiseless',...
        model_parameters.numOfElectrodesPerRing,...
        model_parameters.numOfRings,...
        model_parameters.numOfSensors);
    else
    name = sprintf('data_E_%i_R_%i_M_%i_SNR_%s',...
        model_parameters.numOfElectrodesPerRing,...
        model_parameters.numOfRings,...
        model_parameters.numOfSensors,...
        SNRdb);
    end
    file_name = strcat(data_folder,'\',name);
end