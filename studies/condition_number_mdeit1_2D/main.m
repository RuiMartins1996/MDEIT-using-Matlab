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

rng(1)
%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);
%% Assign the parameters for several models (should create utility functions for this)


global SNRdb
SNRdb = 30;

maxsz_simulation = 10e-3/l0;

maxsz_reconstruction = 10e-3/l0;

background_conductivity = 3.28e-1/sigma0;  %page 163 mentions a saline solution (NaCl + water) at 0.2% mass concentration, but can't find data for that conductivity, check notes

min_sensor_radius = 45e-3/l0;
max_sensor_radius = 90e-3/l0;

n_steps = 10;

select_sensor_axis = 3;


global num_noise_repetitions
num_noise_repetitions = 100;

lambda_vector = logspace(-15,3,30);

%% Define template forward model (2D real tank experiment)

% maxsz is determined by mesh convergence study
model_parameters.maxsz = maxsz_simulation;

% Tank dimensions in characteristic units
model_parameters.isCylindrical = true; % page 163
model_parameters.height= 120e-3/l0; %(m) page 163 
model_parameters.radius= 40e-3/l0; %(m) page 163 

%Electrode configuration and properties
model_parameters.numOfRings = 4;
model_parameters.numOfElectrodesPerRing = 4; 

% what's the size of the electrodes? He mentions EEG for the tank experiment, but what for the simulated experiment?
model_parameters.electrodeRadius = 10e-3/l0; %(mm) unknown from the thesis, but good estimate for Ag/AgCl EEG electrode 
assert(model_parameters.electrodeRadius*model_parameters.numOfElectrodesPerRing < 2*pi*model_parameters.radius)
model_parameters.electrodeContactImpedance = 0.0058/z0; 

% Magnetic sensor configuration and properties
model_parameters.numOfSensors = 16; %(mm) pg 172
model_parameters.sensorRadius = 70e-3/l0;  %(mm) pg 172
model_parameters.mu0 = 1; %With this definition, we're computing the magnetic field intensity H in units of I0/l0;

anomaly_conductivity = 1e-12/sigma0; % (0 according to Kai's thesis, but check notes)
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
    'type','spherical',...
    'name','my_anomaly',...
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


%% Sweep model_parameters over multiple cylindrical radius
model_parameters_array = ...
    sweep_model_parameters({'sensorRadius'},min_sensor_radius,max_sensor_radius,n_steps,model_parameters,'log');
%% Make inverse models (fwd_models are assigned lated)
imdl_mdeit_1_min= eidors_obj('inv_model','my_inverse_model');

% 1-axis imdl
imdl_mdeit_1_min.recon_mode = 'mdeit1';
imdl_mdeit_1_min.select_sensor_axis = 3; % z-axis only page 172
imdl_mdeit_1_min.tol = 1e-5;

imdl_mdeit_1_min.fwd_model = []; %Use a different forward model for reconstruction
imdl_mdeit_1_min.jacobian_bkgnd = struct('value',background_conductivity);
imdl_mdeit_1_min.solver = 'gn';
imdl_mdeit_1_min.RtR_prior = @(x,Jx) x; %tykhonov
imdl_mdeit_1_min.recon_type = 'difference';

imdl_mdeit_1_min.verbose = true; % print debug info
% imdl_mdeit_1_min.solver = 'tsvd';

%Same for max sensor radius
imdl_mdeit_1_max = imdl_mdeit_1_min;

% MDEIT3
imdl_mdeit_3_min = imdl_mdeit_1_min;
imdl_mdeit_3_min = rmfield(imdl_mdeit_3_min,'select_sensor_axis');
imdl_mdeit_3_min.recon_mode = 'mdeit3';

imdl_mdeit_3_max = imdl_mdeit_1_max;
imdl_mdeit_3_max = rmfield(imdl_mdeit_3_max,'select_sensor_axis');
imdl_mdeit_3_max.recon_mode = 'mdeit3';

% eit imdl
% eidors_cache('clear_all'); %For safety, because EIDORS is using the cache forward solve from a different forward model.
imdl_eit = mk_common_model('a2c2',8); % Will replace most fields
imdl_eit.jacobian_bkgnd = struct('value',background_conductivity); %If I
% use this, the solution blows up
imdl_eit.fwd_model = []; %Use a different forward model for reconstruction
imdl_eit.recon_mode = 'eit';
imdl_eit.reconst_type = 'difference';
imdl_eit.RtR_prior = @prior_tikhonov;

%% Make a stimulation pattern, same for all models
stimulation = mk_stim_patterns(model_parameters.numOfElectrodesPerRing*model_parameters.numOfElectrodesPerRing,1,inj,meas,options,current_amplitude);

%% Compute the condition number of Jacobian for MDEIT and EIT
condition_number_array_eit = nan(numel(model_parameters_array),1);
condition_number_array_mdeit = nan(numel(model_parameters_array),1);
condition_number_array_mdeit_3= nan(numel(model_parameters_array),1);
r_vector = nan(numel(model_parameters_array),1);

%EIT does not depend on sensor radius
model_parameters = model_parameters_array(1);

% Edit these fields so the model has no knowledge of the anomaly
model_parameters.material = struct();
model_parameters.maxsz = maxsz_reconstruction;

condition_number_array_eit(:) = compute_jacobian_condition_number_r(...
    model_parameters,model_folder,stimulation,background_conductivity,'eit',[]);

for n = 1:numel(model_parameters_array)
    %% Create model
    model_parameters = model_parameters_array(n);
    
    % Edit these fields so the model has no knowledge of the anomaly
    model_parameters.material = struct();
    model_parameters.maxsz = maxsz_reconstruction;
    
    r_vector(n) = model_parameters.sensorRadius;

    condition_number_array_mdeit(n) = compute_jacobian_condition_number_r(...
        model_parameters,model_folder,stimulation,background_conductivity,'mdeit1',select_sensor_axis);

    condition_number_array_mdeit_3(n) = compute_jacobian_condition_number_r(...
        model_parameters,model_folder,stimulation,background_conductivity,'mdeit3',[]);
end

%% Find minimum of condition number and reconstruct for that
id_min_mdeit1 = find(condition_number_array_mdeit == min(condition_number_array_mdeit));
id_max_mdeit1 = find(r_vector == min(r_vector));

id_min_mdeit3 = find(condition_number_array_mdeit_3 == min(condition_number_array_mdeit_3));
id_max_mdeit3 = find(r_vector == min(r_vector));

%% Plots condition number

figure
hold on
plot(r_vector,condition_number_array_eit(:),'b.-')
plot(r_vector,condition_number_array_mdeit(:),'r.-')
plot(r_vector,condition_number_array_mdeit_3(:),'g.-')

plot(r_vector(id_min_mdeit1),condition_number_array_mdeit(id_min_mdeit1),'r*')
plot(r_vector(id_max_mdeit1),condition_number_array_mdeit(id_max_mdeit1),'b*')

grid on;grid minor;
box on;

set(gca,'YScale','log')
min_y = min([condition_number_array_eit(:);condition_number_array_mdeit(:);condition_number_array_mdeit_3(:)]);
max_y = max([condition_number_array_eit(:);condition_number_array_mdeit(:);condition_number_array_mdeit_3(:)]);

ylim([min_y*0.5,1.1*max_y])

xlabel('$r_{sensors}$ (units of $l_0$)','Interpreter','latex');
ylabel('$\kappa$','Interpreter','latex')
legend('EIT','MDEIT-$1$ axis','MDEIT-$3$ axis','interpreter','latex')
%% Generate simulated data for min/max condition number situation

model_parameters = model_parameters_array(1);
[du,ui,uh,imgi] = generate_simulated_data(model_parameters,model_folder,stimulation,background_conductivity,anomaly_conductivity,'eit');

model_parameters = model_parameters_array(1);
model_parameters.sensorRadius = r_vector(id_min_mdeit1);

[dBz_min,Bzi_min,Bzh_min] = generate_simulated_data(model_parameters,model_folder,stimulation,background_conductivity,anomaly_conductivity,'mdeit1');
[dB_min,Bi_min,Bh_min] = generate_simulated_data(model_parameters,model_folder,stimulation,background_conductivity,anomaly_conductivity,'mdeit3');

model_parameters = model_parameters_array(1);
model_parameters.sensorRadius = r_vector(id_max_mdeit1);

[dBz_max,Bzi_max,Bzh_max] = generate_simulated_data(model_parameters,model_folder,stimulation,background_conductivity,anomaly_conductivity,'mdeit1');
[dB_max,Bi_max,Bh_max] = generate_simulated_data(model_parameters,model_folder,stimulation,background_conductivity,anomaly_conductivity, ...
    'mdeit3');

%% Plot generated data

figure
subplot(2,3,1)
plot(du)
subplot(2,3,2)
plot(dBz_min)
subplot(2,3,3)
plot(dBz_max)
subplot(2,3,4)
plot(dB_min)
subplot(2,3,5)
plot(dB_max)

%% Generate coarse forward model for reconstruction (different mesh than the data)

model_parameters = model_parameters_array(1);
model_parameters.material = struct();
model_parameters.maxsz = maxsz_reconstruction;
model_parameters.sensorRadius = r_vector(id_min_mdeit1);

[~,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_reconstruction_min = fmdls{1};
fmdl_reconstruction_min.stimulation = stimulation;

model_parameters = model_parameters_array(1);
model_parameters.material = struct();
model_parameters.maxsz = maxsz_reconstruction;
model_parameters.sensorRadius = r_vector(id_max_mdeit1);

[~,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_reconstruction_max = fmdls{1};
fmdl_reconstruction_max.stimulation = stimulation;

%% Assign fmdl_reconstruction to inverse models
imdl_eit.fwd_model = fmdl_reconstruction_min; %Use a different forward model for reconstruction

imdl_mdeit_1_min.fwd_model = fmdl_reconstruction_min; %Use a different forward model for reconstruction
imdl_mdeit_1_max.fwd_model = fmdl_reconstruction_max; %Use a different forward model for reconstruction

imdl_mdeit_3_min.fwd_model = fmdl_reconstruction_min;
imdl_mdeit_3_max.fwd_model = fmdl_reconstruction_max;


%%
model_parameters = model_parameters_array(1);
model_parameters.material = struct();
model_parameters.maxsz = maxsz_reconstruction;
model_parameters.sensorRadius = r_vector(id_min_mdeit1);

% DEBUG: Are these two different solves with same hyperparameter giving the
% same results??
% Answer: Yes, the different results is because I'm using different noisy
% data!

[img_output_mdeit_1_min,imdl_mdeit_1_min] = inverse_solve(imdl_mdeit_1_min,lambda_vector,Bzi_min,Bzh_min);
[img_output_mdeit_1_max,imdl_mdeit_1_max] = inverse_solve(imdl_mdeit_1_max,lambda_vector,Bzi_min,Bzh_min);

[img_output_mdeit_3_min,imdl_mdeit_3_min] = inverse_solve(imdl_mdeit_3_min,lambda_vector,Bi_min,Bh_min);
[img_output_mdeit_3_max,imdl_mdeit_3_max] = inverse_solve(imdl_mdeit_3_max,lambda_vector,Bi_min,Bh_min);

a = img_output_mdeit_1_min.elem_data-img_output_mdeit_1_max.elem_data;
b = img_output_mdeit_3_min.elem_data-img_output_mdeit_3_max.elem_data;
figure
hold on
subplot(1,2,1)
plot(img_output_mdeit_1_min.elem_data);
plot(a)
subplot(1,2,2)
plot(img_output_mdeit_3_min.elem_data);
plot(b)
hold off

[img_output_eit,imdl_eit] = inverse_solve(imdl_eit,lambda_vector,ui,uh);
[img_output_mdeit_1_min,imdl_mdeit_1_min] = inverse_solve(imdl_mdeit_1_min,lambda_vector,Bzi_min,Bzh_min);
[img_output_mdeit_1_max,imdl_mdeit_1_max] = inverse_solve(imdl_mdeit_1_max,lambda_vector,Bzi_max,Bzh_max);

[img_output_mdeit_3_min,imdl_mdeit_3_min] = inverse_solve(imdl_mdeit_3_min,lambda_vector,Bi_min,Bh_min);
[img_output_mdeit_3_max,imdl_mdeit_3_max] = inverse_solve(imdl_mdeit_3_max,lambda_vector,Bi_max,Bh_max);

%%

figure
subplot(2,3,1)
show_fem(imgi)
title('Ground Truth','Interpreter','latex')
subplot(2,3,2)
show_fem(img_output_eit)
title('EIT','Interpreter','latex')

subplot(2,3,3)
show_fem(img_output_mdeit_1_min)
title('MDEIT min','Interpreter','latex')

subplot(2,3,4)
show_fem(img_output_mdeit_1_max)
title('MDEIT max','Interpreter','latex')

subplot(2,3,5)
show_fem(img_output_mdeit_3_min)
title('MDEIT3 min','Interpreter','latex')

subplot(2,3,6)
show_fem(img_output_mdeit_3_max)
title('MDEIT3 max','Interpreter','latex')

%% Apply KAI's noise correction
sigma_std_deviation_eit = noise_correction(img_output_eit,imdl_eit,background_conductivity,ui.meas-uh.meas);

sigma_std_deviation_mdeit1_min = noise_correction(img_output_mdeit_1_min,imdl_mdeit_1_min,background_conductivity,dBz_min);
sigma_std_deviation_mdeit1_max = noise_correction(img_output_mdeit_1_max,imdl_mdeit_1_max,background_conductivity,dBz_max);

sigma_std_deviation_mdeit3_min = noise_correction(img_output_mdeit_3_min,imdl_mdeit_3_min,background_conductivity,dB_min);
sigma_std_deviation_mdeit3_max = noise_correction(img_output_mdeit_3_max,imdl_mdeit_3_max,background_conductivity,dB_max);

% Threshold

max_elem_data = max(abs(img_output_eit.elem_data));
ids = abs(img_output_eit.elem_data)<0.5*max_elem_data;
img_output_eit.elem_data(not(ids)) = 0;

max_elem_data = max(abs(img_output_mdeit_1_min.elem_data));
ids = abs(img_output_mdeit_1_min.elem_data)<0.5*max_elem_data;
img_output_mdeit_1_min.elem_data(not(ids)) = 0;

max_elem_data = max(abs(img_output_mdeit_1_max.elem_data));
ids = abs(img_output_mdeit_1_max.elem_data)<0.5*max_elem_data;
img_output_mdeit_1_max.elem_data(not(ids)) = 0;

max_elem_data = max(abs(img_output_mdeit_3_min.elem_data));
ids = abs(img_output_mdeit_3_min.elem_data)<0.5*max_elem_data;
img_output_mdeit_3_min.elem_data(not(ids)) = 0;

max_elem_data = max(abs(img_output_mdeit_3_max.elem_data));
ids = abs(img_output_mdeit_3_max.elem_data)<0.5*max_elem_data;
img_output_mdeit_3_max.elem_data(not(ids)) = 0;

% Noise correct
img_output_eit_corrected = img_output_eit;
img_output_eit_corrected.elem_data = img_output_eit.elem_data./sigma_std_deviation_eit;

img_output_mdeit_1_min_corrected = img_output_mdeit_1_min;
img_output_mdeit_1_min_corrected.elem_data = img_output_mdeit_1_min.elem_data./sigma_std_deviation_mdeit1_min;

img_output_mdeit_1_max_corrected = img_output_mdeit_1_max;
img_output_mdeit_1_max_corrected.elem_data = img_output_mdeit_1_max.elem_data./sigma_std_deviation_mdeit1_max;

img_output_mdeit_3_min_corrected = img_output_mdeit_3_min;
img_output_mdeit_3_min_corrected.elem_data = img_output_mdeit_3_min.elem_data./sigma_std_deviation_mdeit3_min;

img_output_mdeit_3_max_corrected = img_output_mdeit_3_max;
img_output_mdeit_3_max_corrected.elem_data = img_output_mdeit_3_max.elem_data./sigma_std_deviation_mdeit3_max;

%% Show reconstructions

posn= [inf,inf,model_parameters_array(1).material.position(3),1,1];
npoints = 128;

figure('Position',[100,100,1000,500])
subplot(2,6,1)
imgi.calc_colours.npoints= npoints;
show_slices(imgi, posn);
title('Ground Truth','Interpreter','latex')

subplot(2,6,2)
img_output_eit.calc_colours.npoints= npoints;
show_slices(img_output_eit, posn);
title('EIT','Interpreter','latex')

subplot(2,6,3)
img_output_mdeit_1_min.calc_colours.npoints= npoints;
show_slices(img_output_mdeit_1_min, posn);
title('MDEIT: Min condition number','Interpreter','latex')
subplot(2,6,4)
img_output_mdeit_1_max.calc_colours.npoints= npoints;
show_slices(img_output_mdeit_1_max, posn);
title('MDEIT: Max condition number','Interpreter','latex')

subplot(2,6,5)
img_output_mdeit_3_min.calc_colours.npoints= npoints;
show_slices(img_output_mdeit_3_min, posn);
title('MDEIT3: Min condition number','Interpreter','latex')
subplot(2,6,6)
img_output_mdeit_3_max.calc_colours.npoints= npoints;
show_slices(img_output_mdeit_3_max, posn);
title('MDEIT3: Max condition number','Interpreter','latex')


subplot(2,6,7)
show_slices(imgi, posn);
title('Ground Truth','Interpreter','latex')

subplot(2,6,8)
img_output_eit_corrected.calc_colours.npoints= npoints;
show_slices(img_output_eit_corrected,posn)
title('EIT corrected','Interpreter','latex')

subplot(2,6,9)
img_output_mdeit_1_min_corrected.calc_colours.npoints= npoints;
show_slices(img_output_mdeit_1_min_corrected,posn);
title('MDEIT: Min condition number, corrected','Interpreter','latex')

subplot(2,6,10)
img_output_mdeit_1_max_corrected.calc_colours.npoints= npoints;
show_slices(img_output_mdeit_1_max_corrected,posn);
title('MDEIT: Max condition number, corrected','Interpreter','latex')

subplot(2,6,11)
img_output_mdeit_3_min_corrected.calc_colours.npoints= npoints;
show_slices(img_output_mdeit_3_min_corrected,posn);
title('MDEIT3: Min condition number,corrected','Interpreter','latex')

subplot(2,6,12)
img_output_mdeit_3_max_corrected.calc_colours.npoints= npoints;
show_slices(img_output_mdeit_3_max_corrected,posn);
title('MDEIT3: Max condition number,corrected','Interpreter','latex')


figure

subplot(2,3,1)
show_fem(imgi);
title('Ground Truth','Interpreter','latex')

subplot(2,3,2)
show_fem(img_output_eit_corrected)
title('EIT corrected','Interpreter','latex')

subplot(2,3,3)
img_output_mdeit_1_min_corrected.calc_colours.npoints= npoints;
show_fem(img_output_mdeit_1_min_corrected);
title('MDEIT: Min condition number, corrected','Interpreter','latex')

subplot(2,3,4)
img_output_mdeit_1_max_corrected.calc_colours.npoints= npoints;
show_fem(img_output_mdeit_1_max_corrected);
title('MDEIT: Max condition number, corrected','Interpreter','latex')

subplot(2,3,5)
img_output_mdeit_3_min_corrected.calc_colours.npoints= npoints;
show_fem(img_output_mdeit_3_min_corrected);
title('MDEIT3: Min condition number,corrected','Interpreter','latex')

subplot(2,3,6)
img_output_mdeit_3_max_corrected.calc_colours.npoints= npoints;
show_fem(img_output_mdeit_3_max_corrected);
title('MDEIT3: Max condition number,corrected','Interpreter','latex')

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
%% FUNCTION
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


%% FUNCTION
function condition_number = compute_jacobian_condition_number_r(model_parameters,model_folder,stimulation,background_conductivity,recon_mode,select_sensor_axis)

valid_modes = {'mdeit1', 'mdeit3','eit'};

if ~ismember(recon_mode, valid_modes)
    error('invalid recon_mode "%s". Must be''mdeit1'', or ''mdeit3''.', recon_mode);
end

if nargin<6
    select_sensor_axis = [];
end

%% Create/load model

[~,fmdls] = ...
    mk_mdeit_model(model_parameters,model_folder);

fmdl = fmdls{1};
fmdl.stimulation = stimulation;

imgh = mk_image_mdeit(fmdl,background_conductivity);

if strcmp(recon_mode,'eit')
    J = calc_jacobian(imgh);
end

% Compute MDEIT Jacobians
lambdatimesdAdp = @(lambda) computeLambdaTimesDaDp(imgh,lambda);
A = @(sigma) M(imgh,sigma);

if strcmp(recon_mode,'mdeit1')
    J = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',select_sensor_axis);

elseif strcmp(recon_mode,'mdeit3')
    J = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit3');

end

fprintf('Doing svd\n')
singular_values = svds(J,rank(J));
condition_number = singular_values(1)/singular_values(end);

end

%% FUNCTIONS

function [dM,Mi,Mh,imgi] = generate_simulated_data(model_parameters,model_folder,stimulation,background_conductivity,anomaly_conductivity,recon_mode)

global SNRdb

valid_modes = {'mdeit1', 'mdeit3','eit'};

if ~ismember(recon_mode, valid_modes)
    error('invalid recon_mode "%s". Must be''mdeit1'', or ''mdeit3''.', recon_mode);
end

options = [];

[~,fmdls] = ...
    mk_mdeit_model(model_parameters,model_folder,options);

fmdl = fmdls{1};
fmdl.stimulation = stimulation;

% For min sensor radius
imgh = mk_image_mdeit(fmdl,background_conductivity);
imgi = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

% Forward solve
[datah,uh] = fwd_solve_mdeit(imgh);
[datai,ui] = fwd_solve_mdeit(imgi);
Bzh = datah.Bz(:);
Bzi = datai.Bz(:);
Bh = [datah.Bx(:);datah.By(:);datah.Bz(:)];
Bi = [datai.Bx(:);datai.By(:);datai.Bz(:)];

if strcmp(recon_mode,'eit')
    du = add_measurement_noise_difference(ui.meas,uh.meas,SNRdb);
    ui.meas = uh.meas+du;
    dM = du;
    Mi = ui;
    Mh = uh;
    return
elseif strcmp(recon_mode,'mdeit1')
    dBz = add_measurement_noise_difference(Bzi,Bzh,SNRdb);
    Bzi = Bzh+dBz;
    
    dM = dBz;
    Mi = Bzi;
    Mh = Bzh;
    return;
elseif strcmp(recon_mode,'mdeit3')
    dB = add_measurement_noise_difference(Bi,Bh,SNRdb);
    Bi = Bh+dB;

    dM = dB;
    Mi = Bi;
    Mh = Bh;
end


end


%% FUNCTION
function [img_output,imdl] = inverse_solve(imdl,lambda_vector,datai,datah)
    
    if isempty(imdl.fwd_model) || not(isfield(imdl,'fwd_model'))
        error('here')
    end
    
    if isstruct(datai) && isstruct(datah)
        dM = datai.meas-datah.meas;
    else
        dM = datai-datah;
    end

    % GCV
    [lambda_optimal,optimal_id,V_mu,~] = generalized_cross_validation(imdl,dM,lambda_vector);
    if strcmp(imdl.recon_mode,'eit')
        imdl.hyperparameter  = struct('value',sqrt(lambda_optimal));
    else
        imdl.hyperparameter = struct('value',lambda_optimal);
    end

    % Inverse solve
    fprintf('Solving\n');
    if strcmp(imdl.recon_mode,'eit')
        %% DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % EIT SOLVE IS TAKING TO LONG FOR LARGE JACOBIAN
        % img_output = inv_solve(imdl,datah,datai);
        img_output = mk_image(imdl);
    else
        img_output = inv_solve_mdeit(imdl,datah,datai);
    end

end

%% FUNCTION: Noise correction
function sigma_std_deviation = noise_correction(img,imdl,background_conductivity,diff_data)

global num_noise_repetitions SNRdb

if strcmp(imdl.recon_mode,'mdeit1') || strcmp(imdl.recon_mode,'mdeit3')
    hyperparameter = imdl.hyperparameter.value;
    if isfield(img,'jacobian')
        J = img.jacobian;
    else
        error('here');

    end
elseif strcmp(imdl.recon_mode,'eit')
    imgh = mk_image(imdl.fwd_model,background_conductivity);
    J = calc_jacobian(imgh);
    hyperparameter = imdl.hyperparameter.value.^2;
end

[U,S,V] = svd(J,'econ');

n_elem = numel(img.elem_data);

M = zeros(n_elem,num_noise_repetitions);

noise_level= max(abs(diff_data-mean(diff_data)))/10^(SNRdb/20);

for t = 1:num_noise_repetitions
    
    fprintf('Solving for noise realization %i\n',t);

    data_noisy = noise_level*randn(size(diff_data,1),1);

    sv = diag(S)+hyperparameter./diag(S);
    sigma_noisy = V*diag(1./sv)*U'*data_noisy;
    
    M(:,t) = sigma_noisy;
end

sigma_std_deviation = std(M,[],2);

end
