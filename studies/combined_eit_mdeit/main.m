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

%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);

% Set the minimum sensor radius
sensor_radius_0 = 1.01;
rmin = sensor_radius_0;
rmax = 3;

% Number of noise samples to gather for noise based correction
num_noise_repetitions = 15;

%% Problem parameters
background_conductivity = 3.28e-1/sigma0;
anomaly_conductivity = 1e-12/sigma0;

maxsz_reconstruction = 5e-3/l0;

current_amplitude = 2.4e-3/I0;

% Set noise type and snr
noise_type = 'white';

snr = 20;

% Default 3D model
model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);
% with some changes
model_height = model_parameters.height/2;
model_parameters.numOfElectrodesPerRing = 8;
model_parameters.numOfRings = 4;
model_parameters.numOfSensors = 8*4;
model_parameters.anomaly.type = 'spherical';
model_parameters.material.type = 'spherical';
model_parameters.material.name = 'sphere_anomaly';
%% Stimulation pattern is not the default. For now, edit it manually

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used
stimulation = mk_stim_patterns(model_parameters.numOfElectrodesPerRing,model_parameters.numOfRings,inj,meas,{},current_amplitude);

%% Create model for data generation
[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

fmdl = fmdls{1};
fmdl.stimulation = stimulation;

% Make homogeneous image
imgh = mk_image_mdeit(fmdl,background_conductivity);

imgi = add_material_properties(imgh,[background_conductivity,anomaly_conductivity]);

figure
show_fem_transparent_edges(imgi);
drawnow

%% Generate coarse forward model for reconstruction (different mesh than the data)
model_parameters.material = struct();
model_parameters.maxsz = maxsz_reconstruction;

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_reconstruction = fmdls{1};

% Don't forget to assignt the same stimulation pattern
fmdl_reconstruction.stimulation = stimulation;

n_elem = size(fmdl_reconstruction.elems,1);

%% Make inverse models
imdl_mdeit_3= eidors_obj('inv_model','my_inverse_model');

% 1-axis imdl
imdl_mdeit_3.recon_mode = 'mdeit3';

imdl_mdeit_3.fwd_model = fmdl_reconstruction; %Use a different forward model for reconstruction
imdl_mdeit_3.jacobian_bkgnd = struct('value',background_conductivity);
imdl_mdeit_3.solver = 'gn';
imdl_mdeit_3.RtR_prior = @(x,Jx) x; %tykhonov
imdl_mdeit_3.recon_type = 'difference';
imdl_mdeit_3.verbose = true; % print debug info


imdl_eit = mk_common_model('a2c2',8); % Will replace most fields
imdl_eit.jacobian_bkgnd = struct('value',background_conductivity); %If I
% use this, the solution blows up
imdl_eit.fwd_model = fmdl_reconstruction; %Use a different forward model for reconstruction
imdl_eit.recon_mode = 'eit';
imdl_eit.reconst_type = 'difference';
imdl_eit.RtR_prior = @prior_tikhonov;
imdl_eit.inv_solve_core.do_pcg = true;

%% Generate data
options.noise_structure.type = noise_type;
options.noise_structure.snr = snr;
options.noise_structure.sampling_rate = 1000;
options.noise_structure.epoch_time = 0.01;

[data_b,data_u,snr,noise_b,noise_u] = generate_data(options,imgh,imgi);

% Since this is difference noisy data, it outputs the noisy difference
% measurements, we have to do the following:

% Forward solve
[datah_signal,uh_signal] = fwd_solve_mdeit(imgh);

Bxh_signal = datah_signal.Bx(:);
Byh_signal = datah_signal.By(:);
Bzh_signal = datah_signal.Bz(:);

[datai_signal,ui_signal] = fwd_solve_mdeit(imgi);

Bxi_signal = datai_signal.Bx(:);
Byi_signal = datai_signal.By(:);
Bzi_signal = datai_signal.Bz(:);

% Correction
Bxi_noisy = Bxh_signal + data_b.dBx;
Byi_noisy = Byh_signal + data_b.dBy;
Bzi_noisy = Bzh_signal + data_b.dBz;

ui_noisy = ui_signal;
ui_noisy.meas = uh_signal.meas + data_u.meas_du;

% % 
% figure
% hold on
% plot(Bzi_signal,'r');
% plot(Bzi_noisy,'b.');
% hold off
% 
% figure
% hold on
% plot(ui_noisy.meas,'r');
% plot(ui_signal.meas,'b.');
% hold off

%% Find optimal regularization parameter with GCV
lambda_vector = logspace(-17,3,30);

[lambda_optimal_mdeit1,optimal_id_mdeit_1,V_mu_mdeit1,~] = generalized_cross_validation(imdl_mdeit_3,[data_b.dBx; data_b.dBy; data_b.dBz],lambda_vector);
imdl_mdeit_3.hyperparameter = struct('value',lambda_optimal_mdeit1);


[lambda_optimal_eit,optimal_id_eit,V_mu_eit,~] = generalized_cross_validation(imdl_eit,data_u.meas_du,lambda_vector);
% Since EIDORS uses hp^2*RtR, but the generalized_cross_validation function
% uses hp*RtR, we have to correct for that.
imdl_eit.hyperparameter = struct('value',sqrt(lambda_optimal_eit));

%% Perform reconstruction

% Reconstruct (difference) for 1-MDEIT
Bh = [Bxh_signal;Byh_signal;Bzh_signal];
Bi_noisy = [Bxi_noisy;Byi_noisy;Bzi_noisy];

img_output_mdeit_3 = inv_solve_mdeit(imdl_mdeit_3,Bh,Bi_noisy);

% Reconstruct (difference) for EIT
img_output_eit = inv_solve(imdl_eit,uh_signal,ui_noisy);

%%
figure
subplot(1,2,1)
show_fem_transparent_edges(img_output_mdeit_3);
subplot(1,2,2)
show_fem_transparent_edges(img_output_eit);


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