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
%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);
H0 = I0/l0; %With this definition, we're computing the magnetic field intensity H in units of I0/l0;

%% Assign the parameters for several models (should create utility functions for this)

SNRdb = 40;

num_noise_repetitions = 100;

minsz_mesh_convergence = 0.05;
maxsz_mesh_convergence = 0.2;
num_meshes_mesh_convergence = 5;

maxsz_reconstruction = 0.05;

background_conductivity = 3.28e-1/sigma0;  %page 163 mentions a saline solution (NaCl + water) at 0.2% mass concentration, but can't find data for that conductivity, check notes

% lambda_vector = logspace(-15,1,5);
lambda_vector = 10.^linspace(-12,2,10);


%% Define forward model (2D real tank experiment)

% maxsz is determined by mesh convergence study
model_parameters.maxsz = 5e-3/l0;

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


%DEBUG!!!!!!!!!
% model_parameters.sensorRadius = 45e-3/l0;  %(mm) pg 172


model_parameters.mu0 = 1; % setting mu0 to 1 is the same as computing H in units of I0/l0


anomaly_conductivity = 1e-12/l0; % (0 according to Kai's thesis, but check notes)
anomaly_position = [20,0,35]*1e-3/l0; % pg 172
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

%% Visualize pattern 1 and conductivity anomaly

% Conductivity

% Electric potential
img_visualize_ui = rmfield(imgi, 'elem_data');
img_visualize_ui.node_data = ui.volt(:,1);

img_visualize_uh = rmfield(imgh, 'elem_data');
img_visualize_uh.node_data = uh.volt(:,1);

img_visualize_J = rmfield(imgi, 'elem_data');

% Current density norm
Jx = -imgi.elem_data.*imgi.fwd_model.G.Gx*ui.volt(:,1);
Jy = -imgi.elem_data.*imgi.fwd_model.G.Gy*ui.volt(:,1);
Jz = -imgi.elem_data.*imgi.fwd_model.G.Gz*ui.volt(:,1);

normJ = sqrt(Jx.^2+Jy.^2+Jz.^2);

img_visualize_J.elem_data = normJ;

% Set up colorbars

% Compute colors and scaling limits
imgi.calc_colours.ref_level = background_conductivity;  %  centre of the colour scale
imgi.calc_colours.clim = abs(background_conductivity-anomaly_conductivity); %  max diff from ref_level

img_visualize_ui.calc_colours.ref_level = mean(ui.volt(:,1));
img_visualize_ui.calc_colours.clim = max(ui.volt(:,1)-mean(ui.volt(:,1)));

img_visualize_uh.calc_colours.ref_level = mean(uh.volt(:,1));
img_visualize_uh.calc_colours.clim = max(uh.volt(:,1)-mean(uh.volt(:,1)));

img_visualize_J.calc_colours.ref_level = mean(normJ);
img_visualize_J.calc_colours.clim = max(normJ);

% calc_colours
figure('Position',[200,100,1200,600]);
hold on
subplot(3,2,1)

hh = show_fem(imgi);                % draw the model (hh may be a handle or array)
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);

plot_sensors(fmdl_converged)
if dim ==2
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius])
    axis square
else
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius...
        0 model_parameters.height])
end
eidors_colourbar(imgi);
title('$\sigma_{ground \: truth}$ - units $z_0^{-1}l_0$','Interpreter','latex')
box on;
subplot(3,2,2)
hold on

hh = show_fem(img_visualize_ui);              
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);

plot_sensors(fmdl_converged)
if dim ==2
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius])
    axis square
else
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius...
        0 model_parameters.height])
end
hold off
eidors_colourbar(img_visualize_ui);
box on;
title('$u_i$ - injection pattern $1$ - units $z_0 l_0^{-2} I_0$','Interpreter','latex')
hold off

subplot(3,2,3)
hold on

hh = show_fem(img_visualize_uh);            
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);

plot_sensors(fmdl_converged)
if dim ==2
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius])
    axis square
else
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius...
        0 model_parameters.height])
end
cb = eidors_colourbar(img_visualize_uh);

hold off
box on;
title('$u_h$ - injection pattern $1$ - units $z_0 l_0^{-2} I_0$','Interpreter','latex')
hold off


subplot(3,2,4)
hold on

hh = show_fem(img_visualize_J);            
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);

plot_sensors(fmdl_converged)
if dim ==2
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius])
    axis square
else
    axis(...
        [-model_parameters.sensorRadius model_parameters.sensorRadius ...
        -model_parameters.sensorRadius model_parameters.sensorRadius...
        0 model_parameters.height])
end
cb = eidors_colourbar(img_visualize_J);

hold off
box on;
title('$||\vec{J}||$ - injection pattern $1$ - units $I_0l_0^{-2}$','Interpreter','latex')
hold off

subplot(3,2,5)
hold on
plot(1:numel(Bzh),Bzh);

xlabel('measurement index','Interpreter','latex');
ylabel('$H^z$','Interpreter','latex');

hold off
box on;
title('$H^z_{homogeneous}$ - units $I_0l_0^{-1}$','Interpreter','latex')
hold off
axis square

pause(1e-5);

%% Add measurement noise
Bzh = add_measurement_noise(Bzh,SNRdb);
Bzi = add_measurement_noise(Bzi,SNRdb);

uh.meas = add_measurement_noise(uh.meas,SNRdb);
ui.meas = add_measurement_noise(ui.meas,SNRdb);

%% Plot real data vs noisy data

dB = Bzi-Bzh;
du = ui.meas-uh.meas;

dB_real = Bzi_real-Bzh_real;
du_real = ui_real.meas-uh_real.meas;

noise_mdeit = dB-dB_real;
noise_eit = du-du_real;

figure
subplot(3,2,1)
hold on
plot(Bzi_real,'b');
plot(Bzh_real,'g-.');
plot(Bzi,'r.');

ylabel('$H^z_i$','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
legend('Real inh','Real h','Noisy inh')
title('$H^z_i$','Interpreter','latex');

subplot(3,2,2)
hold on
plot(ui_real.meas,'b');
plot(uh_real.meas,'g-.');

plot(ui.meas,'r.');
ylabel('$u_i$','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
legend('Real inh','Real h','Noisy inh')
title('$u_i$','Interpreter','latex');


subplot(3,2,3)
hold on
plot(dB_real,'b');
plot(dB,'r.');

id_mdeit = find(abs(noise_mdeit) == max(abs(noise_mdeit)));
line = [id_mdeit,dB_real(id_mdeit);id_mdeit,dB(id_mdeit)];
plot(line(:,1),line(:,2),'k-','LineWidth',1,'LineStyle','--');
plot(line(1,1),line(1,2),'k*');

ylabel('$\Delta H^z$','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
legend('Real','Noisy','Max Noise')
title('$H^z_i-H^z_h$','Interpreter','latex');

subplot(3,2,4)
hold on
plot(du_real,'b');
plot(du,'r.');

id_eit = find(abs(noise_eit) == max(abs(noise_eit)));
line = [id_eit,du_real(id_eit);id_eit,du(id_eit)];
plot(line(:,1),line(:,2),'k-','LineWidth',1,'LineStyle','--');
plot(line(1,1),line(1,2),'k*');


ylabel('$\Delta u$','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
legend('Real','Noisy','Max Noise')

title('$u_i-u_h$','Interpreter','latex');

subplot(3,2,5)
hold on
plot(noise_mdeit,'b');

ylabel('noise','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
title('noise','Interpreter','latex');

subplot(3,2,6)
hold on
plot(noise_eit,'b');

ylabel('noise','Interpreter','latex');
xlabel('Measurement Index','Interpreter','latex')
title('noise','Interpreter','latex');

relative_change_eit = abs(du_real)./(abs(uh_real.meas));
relative_change_mdeit = abs(dB_real)./(abs(Bzh_real));

% Print some diagonostics
fprintf('------------------------------------------\n')
max_noise_eit = abs(noise_eit(id_eit));
max_noise_mdeit = abs(noise_mdeit(id_mdeit));

fprintf('Max noise EIT [%%]: %.3g %% of max signal\n',...
    100*max_noise_eit/max(du_real));
fprintf('Max noise MDEIT [%%]: %.3g %% of max signal\n',...
    100*max_noise_mdeit/max(dB_real));


pause(1e-5)
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
imdl_mdeit_1.solver = 'gn';
imdl_mdeit_1.RtR_prior = @(x,Jx) x; %tykhonov
imdl_mdeit_1.recon_type = 'difference';

imdl_mdeit_1.verbose = true; % print debug info
imdl_mdeit_1.tol = 1e-12;

% 3-axis imdl
imdl_mdeit_3 = rmfield(imdl_mdeit_1, 'select_sensor_axis');
imdl_mdeit_3.recon_mode = 'mdeit3';

% eit imdl
% eidors_cache('clear_all'); %For safety, because EIDORS is using the cache forward solve from a different forward model.
imdl_eit = mk_common_model('a2c2',8); % Will replace most fields
imdl_eit.jacobian_bkgnd = struct('value',background_conductivity); %If I
% use this, the solution blows up
imdl_eit.fwd_model = fmdl_reconstruction; %Use a different forward model for reconstruction
imdl_eit.recon_mode = 'eit';
imdl_eit.reconst_type = 'difference';
imdl_eit.RtR_prior = @prior_tikhonov;
imdl_eit.tol = 1e-12;

%% Find optimal regularization parameter with GCV 
[lambda_optimal_mdeit1,optimal_id_mdeit_1,V_mu_mdeit1,~] = generalized_cross_validation(imdl_mdeit_1,dB,lambda_vector);
[lambda_optimal_eit,optimal_id_eit,V_mu_eit,~] = generalized_cross_validation(imdl_eit,du,lambda_vector);

%% Show the generalized cross validation results
figure;
subplot(1,2,1);
hold on
plot(lambda_vector,V_mu_eit);
plot(lambda_vector(optimal_id_eit),V_mu_eit(optimal_id_eit),'r.');
box on;grid on;grid minor;
set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel('$\lambda$','Interpreter','latex')
ylabel('$V(\lambda)$','Interpreter','latex')
title('GCV EIT','Interpreter','latex')

subplot(1,2,2);
hold on
plot(lambda_vector,V_mu_mdeit1);
plot(lambda_vector(optimal_id_mdeit_1),V_mu_mdeit1(optimal_id_mdeit_1),'r.');
set(gca,'YScale','log');
set(gca,'XScale','log');box on;grid on;grid minor;
xlabel('$\lambda$','Interpreter','latex')
ylabel('$V(\lambda)$','Interpreter','latex')
title('GCV MDEIT','Interpreter','latex')

%% Perform reconstruction for multiple regularization parameters

%Jacobian is the same in all iterations, so in the first iteration compute
%the jacobian, but on the others re-use the pre-computed one.

fprintf('Iteration %i\n',1);

% Set hyperparameters
imdl_mdeit_1.hyperparameter = struct('value',lambda_vector(1));
% Since EIDORS uses hp^2*RtR, but the generalized_cross_validation function
% uses hp*RtR, we have to correct for that.
imdl_eit.hyperparameter = struct('value',sqrt(lambda_vector(1)));

% Reconstruct
img_output_mdeit_1 = inv_solve_mdeit(imdl_mdeit_1,Bzh,Bzi);
img_output_eit = inv_solve(imdl_eit,uh,ui);

if isfield(img_output_mdeit_1,'jacobian')
    J_mdeit = img_output_mdeit_1.jacobian;
else
    error('Expected to have jacobian here');
end

imgs_reconstructed_mdeit{1} = img_output_mdeit_1; 
imgs_reconstructed_eit{1} = img_output_eit;

for i = 2:numel(lambda_vector)
    fprintf('Iteration %i\n',i);
    % Set hyperparameters
    imdl_mdeit_1.hyperparameter = struct('value',lambda_vector(i));
    imdl_eit.hyperparameter = struct('value',sqrt(lambda_vector(i)));

    % Reconstruct with pre-computed jacobian
    img_output_mdeit_1 = inv_solve_mdeit(imdl_mdeit_1,Bzh,Bzi,J_mdeit);
    img_output_eit = inv_solve(imdl_eit,uh,ui);
    
    imgs_reconstructed_mdeit{i} = img_output_mdeit_1; 
    imgs_reconstructed_eit{i} = img_output_eit;
end

%Reconstruct for optimal GCV hyperparameter too
fprintf('Reconstructing for GCV hyperparameter\n');
% Set hyperparameters
imdl_mdeit_1.hyperparameter = struct('value',lambda_optimal_mdeit1);
imdl_eit.hyperparameter = struct('value',sqrt(lambda_optimal_eit));

% imdl_eit.hyperparameter = struct('value',sqrt(1e-5));
% imdl_mdeit_1.hyperparameter = struct('value',1e-5);

img_output_eit = inv_solve(imdl_eit,uh,ui);
img_output_mdeit_1 = inv_solve_mdeit(imdl_mdeit_1,Bzh,Bzi,J_mdeit);


imgs_reconstructed_mdeit{numel(lambda_vector)+1} = img_output_mdeit_1;
imgs_reconstructed_eit{numel(lambda_vector)+1} = img_output_eit;
%% Show reconstructions

sq = floor(sqrt(numel(imgs_reconstructed_mdeit)))^2; %closest smaller perfect square

if sq == 1
    m = 1;
    n = numel(imgs_reconstructed_mdeit);
else
    m = sqrt(sq)+1;
    n = sqrt(sq);
end

figure('Position',[200,200,1000,500],'Name','EIT Reconstruction');
for i = 1:numel(lambda_vector)
    subplot(m,n,i)
    show_fem_transparent_edges(imgs_reconstructed_eit{i});
    title_str = strcat('$\lambda = $',num2str(lambda_vector(i)));
    title(title_str,'Interpreter','latex')
    box on
end

subplot(m,n,numel(lambda_vector)+1)
show_fem_transparent_edges(imgs_reconstructed_eit{numel(lambda_vector)+1});
title_str = strcat('(GCV) $\lambda = $',num2str(lambda_optimal_eit));
title(title_str,'Interpreter','latex')
box on


figure('Position',[200,200,1000,500],'Name','MDEIT Reconstruction');
for i = 1:numel(lambda_vector)
    subplot(m,n,i)
    show_fem_transparent_edges(imgs_reconstructed_mdeit{i});
    title_str = strcat('$\lambda = $',num2str(lambda_vector(i)));
    title(title_str,'Interpreter','latex')
    box on
end

subplot(m,n,numel(lambda_vector)+1)
show_fem_transparent_edges(imgs_reconstructed_mdeit{numel(lambda_vector)+1});
title_str = strcat('(GCV) $\lambda = $',num2str(lambda_optimal_mdeit1));
title(title_str,'Interpreter','latex')
box on

%% FUNCTIONS

function data_noisy = add_measurement_noise(data,SNRdb)
assert(isnumeric(data) && isvector(data),'data must be a numerical vector');
assert(size(data,1)>=size(data,2),'data must be a column vector');

if isempty(SNRdb)
    data_noisy = data;
    return;
else
%     data_amplitude = max(abs(data-mean(data)));
%     noise_level= data_amplitude/10^(SNRdb/20);
%     data_noisy = data+noise_level*randn(size(data));

    data_amplitudes = abs(data-mean(data));
    noise_levels = data_amplitudes/10^(SNRdb/20);
    
    data_noisy = data;
    for i = 1:length(noise_levels)
        data_noisy(i) = data_noisy(i)+noise_levels(i)*randn();
    end
end
end

function show_fem_transparent_edges(img)

hh = show_fem(img);                % draw the model (hh may be a handle or array)
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);

end