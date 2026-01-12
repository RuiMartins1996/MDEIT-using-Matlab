
clc; clear all; close all;

%% Prepare workspace

% Get the full path of the current script
fullpath = mfilename('fullpath');

% Extract just the folder
script_folder = fileparts(fullpath);

cd(script_folder);

% Set or create data folder
data_folder = strcat(script_folder ,'\data');
if ~exist(data_folder, 'dir')
    mkdir(data_folder);
end
addpath(data_folder);

cd("..\..\");

addpath(genpath("functions"));
addpath(genpath("libraries"));

run("globalParameters.m")

% Set or create model folder
model_folder = './models';
if ~exist(model_folder, 'dir')
    mkdir(model_folder);
end
addpath(genpath("models"));

%% Setup EIDORS
eidors_folder = setupEidors(cd);
clc;
rng(1)

%% Experiment parameters

num_of_repetitions = 3;
tolerances = 10.^(-4:-1:-7);

%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);
%% Assign the parameters for several models (should create utility functions for this)

% For small SNR, noise level might be bigger in magnitude than difference
% data. Does that measurement even make sense?
SNRdb = 100;

% minsz_mesh_convergence = 0.05;
% maxsz_mesh_convergence = 0.2;
% num_meshes_mesh_convergence = 5;

% maxsz_reconstruction = 0.05;

minsz_mesh_convergence = 0.05;
maxsz_mesh_convergence = 0.05;
num_meshes_mesh_convergence = 1;

maxsz_reconstruction = 0.05;

background_conductivity = 3.28e-1/sigma0;  %page 163 mentions a saline solution (NaCl + water) at 0.2% mass concentration, but can't find data for that conductivity, check notes

%% Define forward model (2D real tank experiment)

% maxsz is determined by mesh convergence study
model_parameters.maxsz = 10e-3/l0;

% This data came from Kai's Thesis
model_parameters.is2D = true;

% Tank dimensions in characteristic units
model_parameters.isCylindrical = true; % page 163
model_parameters.height=70e-3/l0; %(m) page 163 
model_parameters.radius=40e-3/l0; %(m) page 163 

%Electrode configuration and properties
model_parameters.numOfRings = 1;
model_parameters.numOfElectrodesPerRing = 16; 

% what's the size of the electrodes? He mentions EEG for the tank experiment, but what for the simulated experiment?
model_parameters.electrodeRadius = 8e-3/l0; %(mm) unknown from the thesis, but good estimate for Ag/AgCl EEG electrode 
assert(model_parameters.electrodeRadius*model_parameters.numOfElectrodesPerRing < 2*pi*model_parameters.radius)
model_parameters.electrodeContactImpedance = 0.0058/z0; 

% Magnetic sensor configuration and properties
model_parameters.numOfSensors = 25; %(mm) pg 172
model_parameters.sensorRadius = 70e-3/l0;  %(mm) pg 172
model_parameters.mu0 = 1; %With this definition, we're computing the magnetic field intensity H in units of I0/l0;


anomaly_conductivity = 1e-12/l0; % (0 according to Kai's thesis, but check notes)
anomaly_position = [20,0,35]*1e-3/l0; %  [0,0,35]*1e-3/l0; % pg 172
anomaly_radius = 25/2*1e-3/l0; % 25mm of diameter, pg 172

% This is deprecated, the .material property overrides this behaviour,
% except for the anomaly conductivity
model_parameters.anomaly = struct(...
    'type','spherical',...
    'conductivity',anomaly_conductivity,...
    'radius',anomaly_radius,...
    'position',anomaly_position);

% A material defines the geometry of that material inside the domain, not
% its properties, like conductivity.
model_parameters.material = struct(...
    'type','cylindrical',...
    'name','plastic_cylinder',...
    'radius',anomaly_radius,...
    'position',anomaly_position);

% electrode height 35 mm is correct given that the ring is placed in the
% middle of the height of the cylinder even though there is no parameter of
% model_parameters controlling it

% Stimulation pattern was not the default. For now, edit it manually
current_amplitude = 2.4e-3/I0;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used
options = {};
%% Mesh convergence study to obtain the "correct" forward model

[fmdl_converged,u_norms,b_norms,n_elems,n_nodes] = ...
    mesh_convergence(model_parameters,model_folder,minsz_mesh_convergence,maxsz_mesh_convergence,num_meshes_mesh_convergence,'multiple_stimulation');

fmdl = fmdl_converged;

dim = size(fmdl.nodes,2);

%% Assign the stimulation pattern to this forward model

stimulation = mk_stim_patterns(numel(fmdl.electrode),1,inj,meas,options,current_amplitude);
fmdl.stimulation = stimulation;
%% Make images from forward models and forward solve

% Make homogeneous image
imgh = mk_image_mdeit(fmdl,background_conductivity);

% Add plastic cylinder
imgi = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

%% Generate data from these forward models, but consider a different one for reconstruction

% Forward solve
[datah,uh] = fwd_solve_mdeit(imgh);
Bzh = datah.Bz(:);

% Forward solve
[datai,ui] = fwd_solve_mdeit(imgi);
Bzi = datai.Bz(:);

Bzh_real = Bzh;
Bzi_real = Bzi;

ui_real = ui;
uh_real = uh;

%% Add measurement noise
% Bzh = add_measurement_noise(Bzh,SNRdb);
% Bzi = add_measurement_noise(Bzi,SNRdb);

% uh.meas = add_measurement_noise(uh.meas,SNRdb);
% ui.meas= add_measurement_noise(ui.meas,SNRdb);

dB = add_measurement_noise_difference(Bzi,Bzh,SNRdb);
du = add_measurement_noise_difference(ui.meas,uh.meas,SNRdb);

Bzi = Bzh+dB;
ui.meas = uh.meas+du;

%% Generate coarse forward model for reconstruction (different mesh than the data)
model_parameters.material = struct();
model_parameters.maxsz = maxsz_reconstruction;

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_reconstruction = fmdls{1};

% Don't forget to assignt the same stimulation pattern
fmdl_reconstruction.stimulation = stimulation;

n_elem = size(fmdl_reconstruction.elems,1);

%% Make inverse model
imdl_mdeit_1= eidors_obj('inv_model','my_inverse_model');

% 1-axis imdl
imdl_mdeit_1.recon_mode = 'mdeit1';
imdl_mdeit_1.select_sensor_axis = 3; % z-axis only page 172


imdl_mdeit_1.fwd_model = fmdl_reconstruction; %Use a different forward model for reconstruction
imdl_mdeit_1.jacobian_bkgnd = struct('value',background_conductivity);
% imdl_mdeit_1.solver = 'gn';
imdl_mdeit_1.solver = 'lm';
imdl_mdeit_1.recon_type = 'absolute';
imdl_mdeit_1.hyperparameter.value = 1e-1;
imdl_mdeit_1.max_iterations = 3;
imdl_mdeit_1.x0 = background_conductivity*ones(size(fmdl_reconstruction.elems,1),1);


imdl_mdeit_1.verbose = true; % print debug info

imdl_mdeit_2= imdl_mdeit_1;

%% NEED to expose this functions

img = mk_image_mdeit(fmdl_reconstruction,background_conductivity);

recon_mode = 'mdeit1';
select_sensor_axis = imdl_mdeit_1.select_sensor_axis;

lambdatimesdAdp = @(lambda) nan; %dummy, its no longer needed 
A = @(sigma) M(img,sigma);

function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
end

res = @(x) calc_residual_mdeit(img, x,Bzi,recon_mode,select_sensor_axis);
jac = @(x) calc_jacobian_mdeit(img, x,lambdatimesdAdp,A,recon_mode,select_sensor_axis);
%% Compare the J0'J0 regularization

x0 = background_conductivity*ones(size(fmdl_reconstruction.elems,1),1);
j0 = jac(x0);
B = sum(j0.^2, 1)';

imdl_mdeit_1.RtR_prior = @(x,J) sum(J.^2, 1)'.*x;
imdl_mdeit_2.RtR_prior = @(x,J) B.*x;

%% Perform reconstruction

% Reconstruct (difference) for 1-MDEIT
img_output_mdeit_1 = inv_solve_mdeit(imdl_mdeit_1,Bzi);
img_output_mdeit_2 = inv_solve_mdeit(imdl_mdeit_2,Bzi);

%% Plots
subplot(1,2,1)
show_fem(img_output_mdeit_1)
eidors_colourbar(img_output_mdeit_1)
subplot(1,2,2)
show_fem(img_output_mdeit_2)
eidors_colourbar(img_output_mdeit_2)


%% Functions
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


