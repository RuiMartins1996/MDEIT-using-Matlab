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


SNR_vector = [20];
noise_type = 'white';

% For optimiziation
sensor_radius_0 = 1.01;
rmin = sensor_radius_0;
rmax = 3;

num_noise_repetitions = 30;

%% Problem parameters
background_conductivity = 3.28e-1/sigma0;
anomaly_conductivity = 1e-12/sigma0;

maxsz_reconstruction = 5e-3/l0;

% model_parameters = create_kai_2d_model_parameters(l0, z0, sigma0, I0);
% model_height = 0;

model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);
model_height = model_parameters.height/2;
%% Stimulation pattern was not the default. For now, edit it manually
current_amplitude = 2.4e-3/I0;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used
stimulation = mk_stim_patterns(model_parameters.numOfElectrodesPerRing,model_parameters.numOfRings,inj,meas,{},current_amplitude);

%% Create modelsensor_radius_0;

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

model_parameters.sensorRadius = sensor_radius_0;
theta = linspace(0,2*pi,model_parameters.numOfSensors);
model_parameters.sensorPositions = ...
    [sensor_radius_0*cos(theta)',sensor_radius_0*sin(theta)',model_height*ones(model_parameters.numOfSensors,1)];

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

fmdl = fmdls{1};
fmdl.stimulation = stimulation;

% Make homogeneous image
imgh = mk_image_mdeit(fmdl,background_conductivity);

imgi = add_material_properties(imgh,[background_conductivity,anomaly_conductivity]);

figure
show_fem(imgi);
pause(1e-10)
%% Optimization
num_of_sensors = model_parameters.numOfSensors;
x0 = model_parameters.sensorPositions;
hold on

% Avoid giving 3rd coordinate
sensor_locations_vector_0 = [x0(:,1);x0(:,2)];

%% Set up worker enviroments

% Start pool
if isempty(gcp('nocreate'))
    parpool('local', 6);

end

% Initialize EIDORS once per worker
parfevalOnAll(@prepare_workspace, 0, script_folder);

%% Try particle swarm

opts = optimoptions('particleswarm', ...
    'Display','iter', ...
    'UseParallel', true, ...
    'SwarmSize', 10, ...        % tuneable
    'MaxIterations', 10);      % tuneable


condition_number_at_x0 = compute_jacobian_mdeit1_condition_number_v2(imgh,sensor_locations_vector_0);

% First compute optimal radius at which condition number is minimal
lb = rmin;
ub = rmax;
obj_r = @(r) compute_cost_function_r(imgh,r,rmax);

if exist('r_pswarm_file.mat','file')
    var = load('r_pswarm_file.mat');

    r_pswarm = var.r_pswarm;
    f_r_pswarm = var.f_r_pswarm;
    condition_number_at_r_pswarm = var.condition_number_at_r_pswarm;

    theta = linspace(0,2*pi,num_of_sensors);
    sensor_positions = ...
        [r_pswarm*cos(theta)',r_pswarm*sin(theta)',model_height*ones(num_of_sensors,1)];
    sensor_locations_r_pswarm = [sensor_positions(:,1);sensor_positions(:,2)];

    sensor_locations_opt_r_pswarm = reshape([sensor_locations_r_pswarm(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);
else
    [r_pswarm, f_r_pswarm] = particleswarm(obj_r,1, lb, ub, opts);

    theta = linspace(0,2*pi,num_of_sensors);
    sensor_positions = ...
        [r_pswarm*cos(theta)',r_pswarm*sin(theta)',model_height*ones(num_of_sensors,1)];
    sensor_locations_r_pswarm = [sensor_positions(:,1);sensor_positions(:,2)];
    condition_number_at_r_pswarm = compute_jacobian_mdeit1_condition_number_v2(imgh,sensor_locations_r_pswarm);

    sensor_locations_opt_r_pswarm = reshape([sensor_locations_r_pswarm(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);

    save(fullfile(data_folder,'r_pswarm_file.mat'),'r_pswarm','f_r_pswarm','condition_number_at_r_pswarm');
end

%% Use a different method than particle swarm

obj_r = @(r) compute_cost_function_r(imgh,r,rmax);

if exist('r_fminbnd_file.mat','file')
    var = load('r_fminbnd_file.mat');

    r_fminbnd = var.r_fminbnd;
    f_r_pswarm = var.f_r_fminbnd;
    condition_number_at_r_fminbnd = var.condition_number_at_r_fminbnd;

    theta = linspace(0,2*pi,num_of_sensors);
    sensor_positions = ...
        [r_fminbnd*cos(theta)',r_fminbnd*sin(theta)',model_height*ones(num_of_sensors,1)];
    sensor_locations_r_fminbnd = [sensor_positions(:,1);sensor_positions(:,2)];

    sensor_locations_opt_r_fminbnd = reshape([sensor_locations_r_fminbnd(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);
else
    options = optimset('Display','iter');
    [r_fminbnd,f_r_fminbnd] = fminbnd(obj_r,rmin,rmax,options);

    theta = linspace(0,2*pi,num_of_sensors);
    sensor_positions = [r_fminbnd*cos(theta)',r_fminbnd*sin(theta)',model_height*ones(num_of_sensors,1)];
    sensor_locations_r_fminbnd = [sensor_positions(:,1);sensor_positions(:,2)];
    condition_number_at_r_fminbnd = compute_jacobian_mdeit1_condition_number_v2(imgh,sensor_locations_r_fminbnd);

    save(fullfile(data_folder,'r_fminbnd_file.mat'),'r_fminbnd','f_r_fminbnd','condition_number_at_r_fminbnd');

end

%% Compute the condition number and cost function for several r for plotting
r_vec = logspace(log10(rmin),log10(rmax)-eps(log10(rmax)),20);

if exist('r_cond_file.mat','file')
    var = load('r_cond_file.mat');

    condition_number_at_r = var.condition_number_at_r;
else

    for i = 1:length(r_vec)

        fprintf('Index %i of %i\n',i,length(r_vec));
        theta = linspace(0,2*pi,num_of_sensors);
        sensor_positions = ...
            [r_vec(i)*cos(theta)',r_vec(i)*sin(theta)',model_height*ones(num_of_sensors,1)];
        sensor_locations_vector_r = [sensor_positions(:,1);sensor_positions(:,2)];

        condition_number_at_r(i) = compute_jacobian_mdeit1_condition_number_v2(imgh,sensor_locations_vector_r);
    end

    save(fullfile(data_folder,'r_cond_file.mat'),'condition_number_at_r');
end

%%
opts = optimoptions('particleswarm', ...
    'Display','iter', ...
    'UseParallel', true, ...
    'SwarmSize', 20, ...        % tuneable
    'MaxIterations', 20);      % tuneable

lb = -rmax*ones(length(sensor_locations_vector_0),1);
ub =  rmax*ones(length(sensor_locations_vector_0),1);
obj = @(sensor_locations_vector) compute_cost_function_v2(imgh,sensor_locations_vector,rmax);

objective_at_x0 = obj(sensor_locations_vector_0);

if exist('x_pswarm_file.mat','file')
    var = load('x_pswarm_file.mat');

    x_pswarm = var.x_pswarm;
    f_x_pswarm = var.f_x_pswarm;
    condition_number_at_x_pswarm = var.condition_number_at_x_pswarm;

    sensor_locations_opt_x_pswarm = reshape([x_pswarm(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);

else

    [x_pswarm, f_x_pswarm] = particleswarm(obj,length(sensor_locations_vector_0), lb, ub, opts);

    sensor_locations_opt_x_pswarm = reshape([x_pswarm(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);

    condition_number_at_x_pswarm = compute_jacobian_mdeit1_condition_number_v2(imgh,x_pswarm);

    save(fullfile(data_folder,'x_pswarm_file.mat'),'x_pswarm', 'f_x_pswarm','condition_number_at_x_pswarm');
end


%%
cost_at_r = compute_cost(condition_number_at_r,r_vec,1,rmax);
vert_line_cost = linspace(min(cost_at_r),max(cost_at_r));
vert_line_condition_number = linspace(min(condition_number_at_r),max(condition_number_at_r));

id_min_cond_number = find(condition_number_at_r==min(condition_number_at_r));
id_cost = find(cost_at_r==min(cost_at_r));


figure

% Subplot 1
subplot(2,2,1)
legendStr = [];
hold on
plot(r_vec,condition_number_at_r)
plot(r_vec(id_min_cond_number),condition_number_at_r(id_min_cond_number),'rx','MarkerSize',10)

legendStr = [legendStr,{'plot','find() minimum'}];

if exist('r_pswarm')
    x = r_pswarm*ones(size(vert_line_cost));
    plot(r_pswarm,condition_number_at_r_pswarm,'r.','MarkerSize',10)
    plot(x,vert_line_condition_number,'--')
    legendStr = [legendStr,'p\_swarm() solution'];
end

axis([min(r_vec) max(r_vec) 0.9*min(condition_number_at_r) 1.1*max(condition_number_at_r)])
legend(legendStr)
set(gca,'YScale','log')
set(gca,'XScale','log')
grid on;grid minor;box on;
title('Condition number w.r.t sensor radius','Interpreter','latex')
hold off

subplot(2,2,3)
legendStr = [];
hold on
plot(r_vec,cost_at_r)
plot(r_vec(id_cost ),cost_at_r(id_cost ),'rx','MarkerSize',10)

legendStr = [legendStr,{'plot','find() minimum'}];

if exist('r_pswarm')
    x = r_pswarm*ones(size(vert_line_cost));
    cost_r_pswarm = compute_cost(condition_number_at_r_pswarm,r_pswarm,1,rmax);

    plot(r_pswarm,cost_r_pswarm,'b.','MarkerSize',10);
    plot(x,vert_line_condition_number,'--')
    legendStr = [legendStr,'p\_swarm() solution'];
end

axis([min(r_vec) max(r_vec) 0.9*min(cost_at_r) 1.1*max(cost_at_r)])
legend(legendStr)
set(gca,'YScale','log')
set(gca,'XScale','log')
grid on;grid minor;box on;

title('Cost w.r.t sensor radius','Interpreter','latex')
hold off

subplot(2,2,2)
legendStr = [];
hold on
plot(r_vec,condition_number_at_r);
plot(r_vec(id_min_cond_number),condition_number_at_r(id_min_cond_number),'rx','MarkerSize',10)
legendStr = {'plot','find() minimum'};

if exist('r_pswarm')
    x = r_fminbnd*ones(size(vert_line_cost));
    plot(r_fminbnd,condition_number_at_r_fminbnd,'b.','MarkerSize',10);
    plot(x,vert_line_condition_number,'--')
    legendStr = [legendStr,'fminbnd() solution'];
end

axis([min(r_vec) max(r_vec) 0.9*min(condition_number_at_r) 1.1*max(condition_number_at_r)])
legend(legendStr)
set(gca,'YScale','log')
set(gca,'XScale','log')
grid on;grid minor;box on;

title('Condition number w.r.t sensor radius','Interpreter','latex')
hold off

subplot(2,2,4)
legendStr = [];
hold on
plot(r_vec,cost_at_r)
plot(r_vec(id_cost ),cost_at_r(id_cost ),'rx','MarkerSize',10)

legendStr = [legendStr,{'plot','find() minimum'}];

if exist('r_fminbnd')
    x = r_fminbnd*ones(size(vert_line_cost));
    cost_r_fminbnd = compute_cost(condition_number_at_r_fminbnd,r_fminbnd,1,rmax);

    plot(r_fminbnd,cost_r_fminbnd,'b.','MarkerSize',10);
    plot(x,vert_line_condition_number,'--')
    legendStr = [legendStr,'fminbnd() solution'];
end

axis([min(r_vec) max(r_vec) 0.9*min(cost_at_r) 1.1*max(cost_at_r)])
legend(legendStr)
set(gca,'YScale','log')
set(gca,'XScale','log')
grid on;grid minor;box on;

title('Cost w.r.t sensor radius','Interpreter','latex')
hold off

% Render now
drawnow;

%% Plot sensor positions

figure

subplot(1,3,1)
hold on
show_fem(imgh)
plot3(sensor_locations_opt_r_pswarm(:,1),sensor_locations_opt_r_pswarm(:,2),sensor_locations_opt_r_pswarm(:,3),'b.','MarkerSize',10)
hold off
title('p\_swarm radius','Interpreter','latex')
axis('square')
box on
view(2)

subplot(1,3,2)
hold on
show_fem(imgh)
plot3(sensor_locations_opt_r_fminbnd(:,1),sensor_locations_opt_r_fminbnd(:,2),sensor_locations_opt_r_fminbnd(:,3),'b.','MarkerSize',10)
hold off
title('fminbnd radius','Interpreter','latex')
axis('square')
box on
view(2)

subplot(1,3,3)
hold on
show_fem(imgh)
plot3(sensor_locations_opt_x_pswarm(:,1),sensor_locations_opt_x_pswarm(:,2),sensor_locations_opt_x_pswarm(:,3),'b.','MarkerSize',10)
hold off
title('p\_swarm free','Interpreter','latex')
axis('square')
box on
view(2)

%% Create several models for forward simulation and for reconstruction

theta = linspace(0,2*pi,model_parameters.numOfSensors);

%Create model for eit

% Need sensors field to avoid errors in generate_data, but won't use the
% simulated magnetic field
model_parameters.sensorRadius = sensor_radius_0;
model_parameters.sensorPositions = ...
    [sensor_radius_0*cos(theta)',sensor_radius_0*sin(theta)',model_height*ones(model_parameters.numOfSensors,1)];

% Can't use {'no_geometry_matrices','no_magnetometers'} as options,
% generate_data does not expect fmdl to have no 'sensors' field
[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_eit = fmdls{1};
fmdl_eit.stimulation = stimulation;

imgh_eit = mk_image_mdeit(fmdl_eit,background_conductivity);
imgi_eit = add_material_properties(imgh_eit, [background_conductivity,anomaly_conductivity]);

%Create model for sensor_radius_0
if exist('sensor_radius_0','var')
    model_parameters.sensorRadius = sensor_radius_0;
    model_parameters.sensorPositions = ...
        [sensor_radius_0*cos(theta)',sensor_radius_0*sin(theta)',model_height*ones(model_parameters.numOfSensors,1)];
    [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

    fmdl_r_0 = fmdls{1};
    fmdl_r_0.stimulation = stimulation;

    imgh_r_0 = mk_image_mdeit(fmdl_r_0,background_conductivity);
    imgi_r_0 = add_material_properties(imgh_r_0, [background_conductivity,anomaly_conductivity]);
end

%Create model for r_pswarm
if exist('sensor_locations_opt_r_pswarm','var')
    model_parameters.sensorRadius = [];
    model_parameters.sensorPositions = sensor_locations_opt_r_pswarm;
    [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

    fmdl_r_pswarm = fmdls{1};
    fmdl_r_pswarm.stimulation = stimulation;

    imgh_r_pswarm = mk_image_mdeit(fmdl_r_pswarm,background_conductivity);
    imgi_r_pswarm = add_material_properties(imgh_r_pswarm, [background_conductivity,anomaly_conductivity]);
end

%Create model for r_fminbnd
if exist('sensor_locations_opt_r_fminbnd','var')

    model_parameters.sensorRadius = [];
    model_parameters.sensorPositions = sensor_locations_opt_r_fminbnd;
    [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

    fmdl_r_fminbnd = fmdls{1};
    fmdl_r_fminbnd.stimulation = stimulation;

    imgh_r_fminbnd = mk_image_mdeit(fmdl_r_fminbnd,background_conductivity);
    imgi_r_fminbnd = add_material_properties(imgh_r_fminbnd, [background_conductivity,anomaly_conductivity]);
end

%Create model for sensor_locations_opt_pswarm
if exist('sensor_locations_opt_x_pswarm','var')

    model_parameters.sensorRadius = [];
    model_parameters.sensorPositions = sensor_locations_opt_x_pswarm;
    [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

    fmdl_x_pswarm = fmdls{1};
    fmdl_x_pswarm.stimulation = stimulation;

    imgh_x_pswarm = mk_image_mdeit(fmdl_x_pswarm,background_conductivity);
    imgi_x_pswarm = add_material_properties(imgh_x_pswarm, [background_conductivity,anomaly_conductivity]);
end

%% Generate coarse forward model for reconstruction (different mesh than the data)
model_parameters.material = struct();
model_parameters.maxsz = maxsz_reconstruction;

% Create model for eit
model_parameters.sensorRadius = [];
%MUST HAVE field sensorPositions to avoid errors later, but won't use it
model_parameters.sensorPositions = reshape([sensor_locations_vector_0(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);;
[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder,{'no_geometry_matrices','no_magnetometers'});

fmdl_reconstruction_eit = fmdls{1};
fmdl_reconstruction_eit.stimulation = stimulation;

% Create model with default sensor positions
if exist('fmdl_r_0','var')
    model_parameters.sensorPositions = reshape([sensor_locations_vector_0(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);
    [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
    fmdl_reconstruction_r_0 = fmdls{1};
    fmdl_reconstruction_r_0.stimulation = stimulation;
end

%Create model for r_pswarm
if exist('fmdl_r_pswarm','var')
    model_parameters.sensorPositions = sensor_locations_opt_r_pswarm;
    [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
    fmdl_reconstruction_r_pswarm = fmdls{1};
    fmdl_reconstruction_r_pswarm.stimulation = stimulation;
end

%Create model for r_fminbnd
if exist('fmdl_r_fminbnd','var')
    model_parameters.sensorPositions = sensor_locations_opt_r_fminbnd;
    [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
    fmdl_reconstruction_r_fminbnd = fmdls{1};
    fmdl_reconstruction_r_fminbnd.stimulation = stimulation;
end

%  Create model for x_pswarm
if exist('fmdl_x_pswarm','var')
    model_parameters.sensorPositions = sensor_locations_opt_x_pswarm;
    [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
    fmdl_reconstruction_x_pswarm = fmdls{1};
    fmdl_reconstruction_x_pswarm.stimulation = stimulation;
end

%% Make inverse model
imdl_mdeit= eidors_obj('inv_model','my_inverse_model');

% 1-axis imdl
imdl_mdeit.recon_mode = 'mdeit1';
imdl_mdeit.select_sensor_axis = 3; % z-axis only page 172
imdl_mdeit.jacobian_bkgnd = struct('value',background_conductivity);
imdl_mdeit.solver = 'gn';
imdl_mdeit.RtR_prior = @(x,Jx) x; %tykhonov
imdl_mdeit.recon_type = 'difference';

imdl_mdeit.verbose = true; % print debug info

% eit imdl
imdl_eit = mk_common_model('a2c2',8); % Will replace most fields
imdl_eit.jacobian_bkgnd = struct('value',background_conductivity); %If I
% use this, the solution blows up
imdl_eit.recon_mode = 'eit';
imdl_eit.reconst_type = 'difference';
imdl_eit.RtR_prior = @prior_tikhonov;
imdl_eit.inv_solve_core.do_pcg = true;

%% Noise generator functions
function data_noisy = noisy_data_generator_mdeit(imgh,imgi,options)
[out_data_b,out_data_u,snr,noise_b,noise_u] =  generate_data(options,imgh,imgi);
data_noisy = noise_b.noise_b_z;
end

function data_noisy = noisy_data_generator_eit(imgh,imgi,options)
[out_data_b,out_data_u,snr,noise_b,noise_u] =  generate_data(options,imgh,imgi);
data_noisy = noise_u;
end

%% Reconstruction
for i = 1:length(SNR_vector)
    %Generate data
    options.noise_structure.type = noise_type;
    options.noise_structure.snr = SNR_vector(i);
    options.noise_structure.epoch_time = 50e-3; %50 milliseconds
    options.noise_structure.sampling_rate = 1000;

    % To find optimal regularization parameter with GCV
    lambda_vector = logspace(-17,3,30);
    
    % For noise correction
    func_eit = @(imgh,imgi) noisy_data_generator_eit(imgh,imgi,options);
    func_mdeit = @(imgh,imgi) noisy_data_generator_mdeit(imgh,imgi,options);

    % EIT
    [~,data_u_eit,snr] = generate_data(options,imgh_eit,imgi_eit);
    uh_eit = fwd_solve(imgh_eit);
    % Bzh_r_0 = datah_r_0.Bz(:);
    % Bzi_r_0 = Bzh_r_0 + data_b_r_0.dBz;
    ui_eit = fwd_solve(imgi_eit);
    ui_eit.meas = uh_eit.meas + data_u_eit.meas_du;
    
    imdl_eit.fwd_model = fmdl_reconstruction_eit; %Use a different forward model for reconstruction
    [lambda_optimal_eit,optimal_id_eit,V_mu_eit,~]  = generalized_cross_validation(imdl_eit,data_u_eit.meas_du,lambda_vector);
    % Since EIDORS uses hp^2*RtR, but the generalized_cross_validation function
    % uses hp*RtR, we have to correct for that.
    imdl_eit.hyperparameter = struct('value',sqrt(lambda_optimal_eit));
    lambda_eit = sqrt(lambda_optimal_eit);

    % Reconstruct (difference) for EIT
    img_output_eit = inv_solve(imdl_eit,uh_eit,ui_eit);

    Jh_eit = calc_jacobian(mk_image(fmdl_reconstruction_eit,background_conductivity));
    sigma_std_eit = noise_correction(imgh_eit,imgi_eit,Jh_eit,lambda_optimal_eit,func_eit,num_noise_repetitions);

    % MDEIT
    n_plots = 0;
    img_outputs = [];
    sigma_stds = [];
    figure_names = [];
    condition_numbers = [];

    % Data for default model
    if exist("fmdl_reconstruction_r_0",'var')
        [data_b_r_0,data_u_r_0,snr] = generate_data(options,imgh_r_0,imgi_r_0);
        [datah_r_0,uh_r_0] = fwd_solve_mdeit(imgh_r_0);
        Bzh_r_0 = datah_r_0.Bz(:);
        Bzi_r_0 = Bzh_r_0 + data_b_r_0.dBz;
        [~,ui_r_0] = fwd_solve_mdeit(imgi_r_0);
        ui_r_0.meas = uh_r_0.meas +data_u_r_0.meas_du;

        imdl_mdeit_r_0 = imdl_mdeit;
        imdl_mdeit_r_0.fwd_model = fmdl_reconstruction_r_0;
        [lambda_optimal_mdeit_r_0,optimal_id_mdeit_r_0,V_mu_mdeit_r_0,~] = generalized_cross_validation(imdl_mdeit_r_0,data_b_r_0.dBz,lambda_vector);
        
        imdl_mdeit_r_0.hyperparameter = struct('value',lambda_optimal_mdeit_r_0);

        img_output_mdeit_r_0 = inv_solve_mdeit(imdl_mdeit_r_0,Bzh_r_0,Bzi_r_0);

        % for k = 1:length(lambda_vector)
        %     imdl_mdeit_r_0.hyperparameter = struct('value',lambda_vector(k));
        %     img_output_mdeit_r_0 = inv_solve_mdeit(imdl_mdeit_r_0,Bzh_r_0,Bzi_r_0);
        %     figure
        % 
        %     subplot(1,2,1)
        %     img_temp = img_output_mdeit_r_0;
        %     img_temp.calc_colours.npoints =128;
        %     img_temp.show_slices.do_colourbar = true;
        %     show_slices(img_temp,[inf,inf,1.5])
        %     title(strcat('lambda = ',num2str(lambda_vector(k))));
        %     subplot(1,2,2)
        %     show_fem_transparent_edges(img_temp)
        %     drawnow
        % end

        Jh_r_0 = img_output_mdeit_r_0.jacobian;
        sigma_std_mdeit_r_0 = noise_correction(imgh_r_0,imgi_r_0,Jh_r_0,lambda_optimal_mdeit_r_0,func_mdeit,num_noise_repetitions);

        n_plots = n_plots +1;
        sigma_stds = [sigma_stds {sigma_std_mdeit_r_0}];
        img_outputs = [img_outputs img_output_mdeit_r_0];
        figure_names = [figure_names,"r\_0"];
        condition_numbers = [condition_numbers condition_number_at_x0];
    end

    % Data for r_pswarm model
    if exist("fmdl_reconstruction_r_pswarm",'var')
        [data_b_r_pswarm,data_u_r_pswarm,snr] = generate_data(options,imgh_r_pswarm,imgi_r_pswarm);
        [datah_r_pswarm,uh_r_pswarm] = fwd_solve_mdeit(imgh_r_pswarm);
        Bzh_r_pswarm = datah_r_pswarm.Bz(:);
        Bzi_r_pswarm = Bzh_r_pswarm + data_b_r_pswarm.dBz;
        [~,ui_r_pswarm] = fwd_solve_mdeit(imgi_r_pswarm);
        ui_r_pswarm.meas = uh_r_pswarm.meas +data_u_r_pswarm.meas_du;

        imdl_mdeit_r_pswarm = imdl_mdeit;
        imdl_mdeit_r_pswarm.fwd_model = fmdl_reconstruction_r_pswarm;
        [lambda_optimal_mdeit_r_pswarm,optimal_id_mdeit_r_pswarm,V_mu_mdeit_r_pswarm,~] = generalized_cross_validation(imdl_mdeit_r_pswarm,data_b_r_pswarm.dBz,lambda_vector);
        
        imdl_mdeit_r_pswarm.hyperparameter = struct('value',lambda_optimal_mdeit_r_pswarm);

        img_output_mdeit_r_pswarm = inv_solve_mdeit(imdl_mdeit_r_pswarm,Bzh_r_pswarm,Bzi_r_pswarm);

        Jh_r_pswarm = img_output_mdeit_r_pswarm.jacobian;
        sigma_std_mdeit_r_pswarm = noise_correction(imgh_r_pswarm,imgi_r_pswarm,Jh_r_pswarm,lambda_optimal_mdeit_r_pswarm,func_mdeit,num_noise_repetitions);

        n_plots = n_plots +1;
        sigma_stds = [sigma_stds {sigma_std_mdeit_r_pswarm}];

        img_outputs = [img_outputs img_output_mdeit_r_pswarm];
        figure_names = [figure_names,"r\_pswarm"];
        condition_numbers = [condition_numbers condition_number_at_r_pswarm];
    end

    % Data for r_fminbnd model
    if exist("fmdl_reconstruction_r_fminbnd",'var')
        [data_b_r_fminbnd,data_u_r_fminbnd,snr] = generate_data(options,imgh_r_fminbnd,imgi_r_fminbnd);
        [datah_r_fminbnd,uh_r_fminbnd] = fwd_solve_mdeit(imgh_r_fminbnd);
        Bzh_r_fminbnd = datah_r_fminbnd.Bz(:);
        Bzi_r_fminbnd = Bzh_r_fminbnd + data_b_r_fminbnd.dBz;
        [~,ui_r_fminbnd] = fwd_solve_mdeit(imgi_r_fminbnd);
        ui_r_fminbnd.meas = uh_r_fminbnd.meas +data_u_r_fminbnd.meas_du;

        imdl_mdeit_r_fminbnd = imdl_mdeit;
        imdl_mdeit_r_fminbnd.fwd_model = fmdl_reconstruction_r_fminbnd;
        [lambda_optimal_mdeit_r_fminbnd,optimal_id_mdeit_r_fminbnd,V_mu_mdeit_r_fminbnd,~] = generalized_cross_validation(imdl_mdeit_r_fminbnd,data_b_r_fminbnd.dBz,lambda_vector);
    
        imdl_mdeit_r_fminbnd.hyperparameter = struct('value',lambda_optimal_mdeit_r_fminbnd);

        img_output_mdeit_r_fminbnd = inv_solve_mdeit(imdl_mdeit_r_fminbnd,Bzh_r_fminbnd,Bzi_r_fminbnd);
        
        Jh_r_fminbnd = img_output_mdeit_r_fminbnd.jacobian;
        sigma_std_mdeit_r_fminbnd = noise_correction(imgh_r_fminbnd,imgi_r_fminbnd,Jh_r_fminbnd,lambda_optimal_mdeit_r_fminbnd,func_mdeit,num_noise_repetitions);

        n_plots = n_plots +1;
        sigma_stds = [sigma_stds {sigma_std_mdeit_r_fminbnd}];

        img_outputs = [img_outputs img_output_mdeit_r_fminbnd];
        figure_names = [figure_names,"r\_fminbnd"];
        condition_numbers = [condition_numbers condition_number_at_r_fminbnd];
    end

    % Data for x_pswarm model
    if exist("fmdl_reconstruction_x_pswarm",'var')
        [data_b_x_pswarm,data_u_x_pswarm,snr] = generate_data(options,imgh_x_pswarm,imgi_x_pswarm);
        [datah_x_pswarm,uh_x_pswarm] = fwd_solve_mdeit(imgh_x_pswarm);
        Bzh_x_pswarm = datah_x_pswarm.Bz(:);
        Bzi_x_pswarm = Bzh_x_pswarm + data_b_x_pswarm.dBz;
        [~,ui_x_pswarm] = fwd_solve_mdeit(imgi_x_pswarm);
        ui_x_pswarm.meas = uh_x_pswarm.meas +data_u_x_pswarm.meas_du;

        imdl_mdeit_x_pswarm = imdl_mdeit;
        imdl_mdeit_x_pswarm .fwd_model = fmdl_reconstruction_x_pswarm;
        [lambda_optimal_mdeit_x_pswarm,optimal_id_mdeit_x_pswarm,V_mu_mdeit_x_pswarm,~] = generalized_cross_validation(imdl_mdeit_x_pswarm,data_b_x_pswarm.dBz,lambda_vector);
        
        imdl_mdeit_x_pswarm.hyperparameter = struct('value',lambda_optimal_mdeit_x_pswarm);

        img_output_mdeit_x_pswarm = inv_solve_mdeit(imdl_mdeit_x_pswarm,Bzh_x_pswarm,Bzi_x_pswarm);
        
        Jh_x_pswarm = img_output_mdeit_x_pswarm.jacobian;
        sigma_std_mdeit_x_pswarm = noise_correction(imgh_x_pswarm,imgi_x_pswarm,Jh_x_pswarm,lambda_optimal_mdeit_x_pswarm,func_mdeit,num_noise_repetitions);

        n_plots = n_plots +1;
        sigma_stds = [sigma_stds {sigma_std_mdeit_x_pswarm}];

        img_outputs = [img_outputs img_output_mdeit_x_pswarm];
        figure_names = [figure_names,"x\_pswarm"];
        
        condition_numbers = [condition_numbers condition_number_at_x_pswarm];
    end
    
    %% Plots 
    figure
    % The slice definition (horizontal at half height)
    n_points = 128;
    half_height = model_parameters.height/2;
    level.centre = [0,0,half_height];
    level.normal_angle = [0,0,1];
    
    % To draw a red circle where the anomaly is supposed to be
    theta = linspace(0, 2*pi, 200);
    pixel_center = [1,1]*0.5*n_points;
    pixel_anomaly_radius = n_points*model_parameters.anomaly.radius/(2*model_parameters.radius);
    pixel_anomaly_x0 = pixel_center(1)+n_points*model_parameters.anomaly.position(1)/(2*model_parameters.radius);
    pixel_anomaly_y0 = pixel_center(2)+n_points*model_parameters.anomaly.position(2)/(2*model_parameters.radius);

    x = pixel_anomaly_x0  + pixel_anomaly_radius*cos(theta);
    y = pixel_anomaly_y0 + pixel_anomaly_radius*sin(theta);

    for n = 1:n_plots
        subplot(4,n_plots+1,n)
        img_temp = img_outputs(n);
        img_temp.calc_colours.npoints = n_points;
        img_temp.show_slices.do_colourbar = true;

        hold on
        show_slices(img_temp,[inf,inf,1.5])
        plot(x, y, 'r--', 'LineWidth', 0.5)
        hold off;
        
        condition_number = condition_numbers(n);
        this_title = strcat(figure_names(n),', $\log(\kappa) =',num2str(log10(condition_number)),'$');
        title(this_title,'Interpreter','latex');
        box on;grid on;grid minor;
        
        subplot(4,n_plots+1,n+n_plots+1)
        show_fem_transparent_edges(img_temp)
        zlim([half_height-model_parameters.height/5,half_height+model_parameters.height/5])
        view(3)
        box on;grid on;grid minor;
    
        subplot(4,n_plots+1,n+2*(n_plots+1))
        img_temp.elem_data = img_temp.elem_data./sigma_stds{n};
        show_fem_transparent_edges(img_temp)
        zlim([half_height-model_parameters.height/5,half_height+model_parameters.height/5])
        view(3)
        box on;grid on;grid minor;

        subplot(4,n_plots+1,n+3*(n_plots+1))
        img_temp = img_outputs(n);
        img_temp.elem_data = img_temp.elem_data./sigma_stds{n};
        img_temp.calc_colours.npoints = n_points;
        img_temp.show_slices.do_colourbar = true;

        hold on
        show_slices(img_temp,[inf,inf,1.5])
        plot(x, y, 'r--', 'LineWidth', 0.5)
        hold off;

        box on;grid on;grid minor;
    end
    
    n = n_plots+1;
    subplot(4,n_plots+1,n)
    img_temp = img_output_eit;
    img_temp.calc_colours.npoints = n_points;
    img_temp.show_slices.do_colourbar = true;
    
    hold on
    show_slices(img_temp,[inf,inf,1.5])
    plot(x, y, 'r--', 'LineWidth', 0.5)
    hold off;

    title('eit','Interpreter','latex');
    box on;grid on;grid minor;

    subplot(4,n_plots+1,n+n_plots+1)
    show_fem_transparent_edges(img_temp)
    zlim([half_height-model_parameters.height/5,half_height+model_parameters.height/5])
    view(3)
    box on;grid on;grid minor;

    subplot(4,n_plots+1,n+2*(n_plots+1))
    img_temp.elem_data = img_temp.elem_data./sigma_std_eit;
    show_fem_transparent_edges(img_temp)
    zlim([half_height-model_parameters.height/5,half_height+model_parameters.height/5])
    view(3)
    box on;grid on;grid minor;

    subplot(4,n_plots+1,n+3*(n_plots+1))
    img_temp = img_output_eit;
    img_temp.elem_data = img_temp.elem_data./sigma_std_eit;
    img_temp.calc_colours.npoints = n_points;
    img_temp.show_slices.do_colourbar = true;

    hold on
    show_slices(img_temp,[inf,inf,1.5])
    plot(x, y, 'r--', 'LineWidth', 0.5)
    hold off;

    box on;grid on;grid minor;

    % %% Plots
    % figure
    % subplot(1,3,1)
    % hold on
    % plot(lambda_vector,V_mu_mdeit_r_0);
    % plot(lambda_vector(optimal_id_mdeit_r_0),V_mu_mdeit_r_0(optimal_id_mdeit_r_0),'r.');
    % set(gca,'XScale','log')
    % set(gca,'YScale','log')
    % xlabel('$\lambda$','Interpreter','latex')
    % hold off
    %
    % subplot(1,3,2)
    % hold on
    % plot(lambda_vector,V_mu_mdeit_r_pswarm);
    % plot(lambda_vector(optimal_id_mdeit_r_pswarm),V_mu_mdeit_r_pswarm(optimal_id_mdeit_r_pswarm),'r.');
    % set(gca,'XScale','log')
    % set(gca,'YScale','log')
    % xlabel('$\lambda$','Interpreter','latex')
    % hold off
    %
    % subplot(1,3,3)
    % hold on
    % plot(lambda_vector,V_mu_mdeit_x_pswarm);
    % plot(lambda_vector(optimal_id_mdeit_x_pswarm),V_mu_mdeit_x_pswarm(optimal_id_mdeit_x_pswarm),'r.');
    % set(gca,'XScale','log')
    % set(gca,'YScale','log')
    % xlabel('$\lambda$','Interpreter','latex')
    % hold off
end

%% Function
function cost = compute_cost(condition_number,r,rmin,rmax)
max_log_condition_number = max(log10(condition_number(:)));
cost = log10(condition_number(:)) + max_log_condition_number*sqrt((r(:)-rmin).^2/(rmax-rmin)^2);
end

%% Function
function condition_number = compute_jacobian_mdeit1_condition_number_v2(imgh,sensor_locations_vector)
if any(sensor_locations_vector < -5) || any(sensor_locations_vector > 5)
    error('Optimizer passed values outside bounds!');
end

if size(sensor_locations_vector,1)==1
    sensor_locations_vector = sensor_locations_vector(:);
end
num_of_sensors = numel(imgh.fwd_model.sensors);


%% Parse sensor_locations_vector
dim = size(imgh.fwd_model.nodes,2);
if dim == 2
    model_height = 0;
elseif dim ==3
    model_height = max(imgh.fwd_model.nodes(:,3))-min(imgh.fwd_model.nodes(:,3));
    model_height = model_height/2;
else
    error('Here');
end
sensor_locations = reshape([sensor_locations_vector(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);

select_sensor_axis = 3;

%% Edit fwd_model so sensors field is changed
for m = 1:num_of_sensors
    imgh.fwd_model.sensors(m).position = sensor_locations(m,:);
end
%% Recompute geometry matrices and reassign to fwd_model
[Rx,Ry,Rz] = compute_r_matrices(imgh.fwd_model);
R = struct('Rx',Rx,'Ry',Ry,'Rz',Rz);

imgh.fwd_model.R = R;
%% Compute MDEIT Jacobian
lambdatimesdAdp = @(lambda) computeLambdaTimesDaDp(imgh,lambda);
A = @(sigma) M(imgh,sigma);

J = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',select_sensor_axis);

%% Do singular value decomposition
singular_values = svds(J,rank(J));
condition_number = singular_values(1)/singular_values(end);
end

%% Function
function cost = compute_cost_function_v2(imgh,sensor_locations_vector,rmax)

dim = size(imgh.fwd_model.nodes,2);
if dim == 2
    model_height = 0;
elseif dim ==3
    model_height = max(imgh.fwd_model.nodes(:,3))-min(imgh.fwd_model.nodes(:,3));
    model_height = model_height/2;
else
    error('Here');
end

condition_number = compute_jacobian_mdeit1_condition_number_v2(imgh,sensor_locations_vector);
num_of_sensors = numel(imgh.fwd_model.sensors);
sensor_locations = reshape([sensor_locations_vector(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);

%% Build cost function
rmin = max(sqrt(imgh.fwd_model.nodes(:,1).^2+imgh.fwd_model.nodes(:,2).^2));

radius = sqrt(sensor_locations(:,1).^2+sensor_locations(:,2).^2);
avg_radius = sum(radius)/num_of_sensors;

cost = compute_cost(condition_number,avg_radius,rmin,rmax);
end

%% FUNCITON
function cost = compute_cost_function_r(imgh,r,rmax)
num_of_sensors = numel(imgh.fwd_model.sensors);

dim = size(imgh.fwd_model.nodes,2);
if dim == 2
    model_height = 0;
elseif dim ==3
    model_height = max(imgh.fwd_model.nodes(:,3))-min(imgh.fwd_model.nodes(:,3));
    model_height = model_height/2;
else
    error('Here');
end

theta = linspace(0,2*pi,num_of_sensors);
sensor_positions = ...
    [r*cos(theta)',r*sin(theta)',model_height*ones(num_of_sensors,1)];
sensor_locations_vector = [sensor_positions(:,1);sensor_positions(:,2)];

condition_number = compute_jacobian_mdeit1_condition_number_v2(imgh,sensor_locations_vector);
sensor_locations = reshape([sensor_locations_vector(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);

%% Build cost function
rmin = max(sqrt(imgh.fwd_model.nodes(:,1).^2+imgh.fwd_model.nodes(:,2).^2));

radius = sqrt(sensor_locations(:,1).^2+sensor_locations(:,2).^2);
avg_radius = sum(radius)/num_of_sensors;

cost = compute_cost(condition_number,avg_radius,rmin,rmax);
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