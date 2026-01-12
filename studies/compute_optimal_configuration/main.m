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


%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);


SNR_vector = [50];
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
    'MaxIterations', 2);      % tuneable


condition_number_at_x0 = compute_jacobian_mdeit1_condition_number_v2(imgh,sensor_locations_vector_0);

% First compute optimal radius at which condition number is minimal
lb = rmin;
ub = rmax;
obj_r = @(r) compute_cost_function_r(imgh,r,rmax);

if exist('ropt_pswarm_file.mat','file')
    var = load('ropt_pswarm_file.mat');

    ropt_pswarm = var.ropt_pswarm;
    fropt_pswarm = var.fropt_pswarm;
    condition_number_at_ropt_pswarm = var.condition_number_at_ropt_pswarm;
    
    theta = linspace(0,2*pi,num_of_sensors);
    sensor_positions = ...
        [ropt_pswarm*cos(theta)',ropt_pswarm*sin(theta)',model_height*ones(num_of_sensors,1)];
    sensor_locations_r_pswarm = [sensor_positions(:,1);sensor_positions(:,2)];
    
    sensor_locations_opt_r_pswarm = reshape([sensor_locations_r_pswarm(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);
else
    [ropt_pswarm, fropt_pswarm] = particleswarm(obj_r,1, lb, ub, opts);
    
    theta = linspace(0,2*pi,num_of_sensors);
    sensor_positions = ...
        [ropt_pswarm*cos(theta)',ropt_pswarm*sin(theta)',model_height*ones(num_of_sensors,1)];
    sensor_locations_r_pswarm = [sensor_positions(:,1);sensor_positions(:,2)];
    condition_number_at_ropt_pswarm = compute_jacobian_mdeit1_condition_number_v2(imgh,sensor_locations_r_pswarm);
    
    sensor_locations_opt_r_pswarm = reshape([sensor_locations_r_pswarm(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);

    save('ropt_pswarm_file.mat','ropt_pswarm','fropt_pswarm','condition_number_at_ropt_pswarm');
end



%% Plot the condition number and cost function for several r
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

    save('r_cond_file.mat','condition_number_at_r');
end



%%
cost_r = compute_cost(condition_number_at_r,r_vec,1,rmax);
id = find(cost_r==min(cost_r));
y1 = linspace(min(cost_r),max(cost_r));
y2 = linspace(min(condition_number_at_r),max(condition_number_at_r));


figure
subplot(1,2,1)
hold on
plot(r_vec,condition_number_at_r)

if exist('ropt_pswarm')
    x = ropt_pswarm*ones(size(y1));
    plot(ropt_pswarm,condition_number_at_ropt_pswarm,'r.','MarkerSize',10)
    plot(x,y2,'--')
end
set(gca,'YScale','log')
set(gca,'XScale','log')

title('Condition number')
hold off
subplot(1,2,2)
hold on
plot(r_vec,cost_r);
plot(r_vec(id),cost_r(id),'rx','MarkerSize',10)

if exist('ropt_pswarm')
    x = ropt_pswarm*ones(size(y1));
    plot(ropt_pswarm,fropt_pswarm,'b.','MarkerSize',10);
    plot(x,y1,'--')
end
set(gca,'YScale','log')
set(gca,'XScale','log')

hold off
title('Cost Function')

%%
opts = optimoptions('particleswarm', ...
    'Display','iter', ...
    'UseParallel', true, ...
    'SwarmSize', 20, ...        % tuneable
    'MaxIterations', 2);      % tuneable

lb = -rmax*ones(length(sensor_locations_vector_0),1);
ub =  rmax*ones(length(sensor_locations_vector_0),1);
obj = @(sensor_locations_vector) compute_cost_function_v2(imgh,sensor_locations_vector,rmax);

objective_at_x0 = obj(sensor_locations_vector_0);

if exist('xopt_pswarm_file.mat','file')
    var = load('xopt_pswarm_file.mat');

    xopt_pswarm = var.xopt_pswarm;
    fopt_pswarm = var.fopt_pswarm;
    condition_number_at_xopt = var.condition_number_at_xopt;

    sensor_locations_opt_pswarm = reshape([xopt_pswarm(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);

else

    [xopt_pswarm, fopt_pswarm] = particleswarm(obj,length(sensor_locations_vector_0), lb, ub, opts);
    
    sensor_locations_opt_pswarm = reshape([xopt_pswarm(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);
    
    condition_number_at_xopt = compute_jacobian_mdeit1_condition_number_v2(imgh,xopt_pswarm);
    
    save('xopt_pswarm_file.mat','xopt_pswarm', 'fopt_pswarm','condition_number_at_xopt');
end

%% Create three models
%Create model for sensor_locations_0
[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
model_parameters.sensorRadius = sensor_radius_0;
theta = linspace(0,2*pi,model_parameters.numOfSensors);
model_parameters.sensorPositions = ...
    [sensor_radius_0*cos(theta)',sensor_radius_0*sin(theta)',model_height*ones(model_parameters.numOfSensors,1)];

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

fmdl_0 = fmdls{1};
fmdl_0.stimulation = stimulation;

imgh_0 = mk_image_mdeit(fmdl_0,background_conductivity);
imgi_0 = add_material_properties(imgh_0, [background_conductivity,anomaly_conductivity]);


%Create model for sensor_sensor_locations_opt_r_pswarm
[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

model_parameters.sensorRadius = [];
%DEBUGGG!!!! Create model of optimal radius not optimal
%config!!!!!!!!!!!!!!!
% model_parameters.sensorPositions = sensor_locations_opt_pswarm;

model_parameters.sensorPositions = sensor_locations_opt_r_pswarm;

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

fmdl_1 = fmdls{1};
fmdl_1.stimulation = stimulation;

imgh_1 = mk_image_mdeit(fmdl_1,background_conductivity);
imgi_1 = add_material_properties(imgh_1, [background_conductivity,anomaly_conductivity]);


%Create model for sensor_locations_opt_pswarm
[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

model_parameters.sensorRadius = [];
model_parameters.sensorPositions = sensor_locations_opt_pswarm;

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

fmdl_2 = fmdls{1};
fmdl_2.stimulation = stimulation;

imgh_2 = mk_image_mdeit(fmdl_2,background_conductivity);
imgi_2 = add_material_properties(imgh_2, [background_conductivity,anomaly_conductivity]);

%% Generate coarse forward model for reconstruction (different mesh than the data)
model_parameters.material = struct();
model_parameters.maxsz = maxsz_reconstruction;

% Create model with default sensor positions
model_parameters.sensorPositions = reshape([sensor_locations_vector_0(:);model_height*ones(num_of_sensors,1)],[num_of_sensors,3]);
[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_reconstruction_r_0 = fmdls{1};
fmdl_reconstruction_r_0.stimulation = stimulation;

%Create model of optimal radius config
model_parameters.sensorPositions = sensor_locations_opt_r_pswarm;
[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_reconstruction_r_1 = fmdls{1};
fmdl_reconstruction_r_1.stimulation = stimulation;

%  Create model with optimal config
model_parameters.sensorPositions = sensor_locations_opt_pswarm;
[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_reconstruction_r_2 = fmdls{1};
fmdl_reconstruction_r_2.stimulation = stimulation;

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

% Different inverse models ( with different fwd models) for each case
imdl_mdeit_0 = imdl_mdeit;
imdl_mdeit_1 = imdl_mdeit;
imdl_mdeit_2 = imdl_mdeit;

imdl_eit = mk_common_model('a2c2',8); % Will replace most fields
imdl_eit.jacobian_bkgnd = struct('value',background_conductivity); %If I
% use this, the solution blows up
imdl_eit.recon_mode = 'eit';
imdl_eit.reconst_type = 'difference';
imdl_eit.RtR_prior = @prior_tikhonov;
imdl_eit.inv_solve_core.do_pcg = true;

function data_noisy = noisy_data_generator_mdeit(imgh,imgi,options)
[out_data_b,out_data_u,snr,noise_b,noise_u] =  generate_data(options,imgh,imgi);
data_noisy = noise_b.noise_b_z;
end

function data_noisy = noisy_data_generator_eit(imgh,imgi,options)
[out_data_b,out_data_u,snr,noise_b,noise_u] =  generate_data(options,imgh,imgi);
data_noisy = noise_u;
end

for i = 1:length(SNR_vector)
    %% Generate data
    options.noise_structure.type = noise_type;
    options.noise_structure.snr = SNR_vector(i);
    options.noise_structure.epoch_time = 50e-3; %50 milliseconds
    options.noise_structure.sampling_rate = 1000;

    % Data for standard model
    [data_b_0,data_u_0,snr] = generate_data(options,imgh_0,imgi_0);

    % Forward solve
    [datah_0,uh_0] = fwd_solve_mdeit(imgh_0);
    Bzh_0 = datah_0.Bz(:);

    % Correction
    Bzi_0 = Bzh_0 + data_b_0.dBz;

    [~,ui_0] = fwd_solve_mdeit(imgi_0);
    ui_0.meas = uh_0.meas +data_u_0.meas_du;
    
    % Data for optimal radius model
    [data_b_1,data_u_1,snr] = generate_data(options,imgh_1,imgi_1);
    
    % Forward solve
    [datah_1,uh_1] = fwd_solve_mdeit(imgh_1);
    Bzh_1 = datah_1.Bz(:);

    % Correction
    Bzi_1 = Bzh_1 + data_b_1.dBz;

    [~,ui_1] = fwd_solve_mdeit(imgi_1);
    ui_1.meas = uh_1.meas +data_u_1.meas_du;

    % Data for optimal config model
    [data_b_2,data_u_2,snr] = generate_data(options,imgh_2,imgi_2);
    
    % Forward solve
    [datah_2,uh_2] = fwd_solve_mdeit(imgh_2);
    Bzh_2 = datah_2.Bz(:);

    % Correction
    Bzi_2 = Bzh_2 + data_b_2.dBz;

    [~,ui_2] = fwd_solve_mdeit(imgi_2);
    ui_2.meas = uh_2.meas +data_u_2.meas_du;

    %% Find optimal regularization parameter with GCV
    lambda_vector = logspace(-17,3,30);

    imdl_mdeit_0.fwd_model = fmdl_reconstruction_r_0; %Use a different forward model for reconstruction

    [lambda_optimal_mdeit_0,optimal_id_mdeit_0,V_mu_mdeit_0,~] = generalized_cross_validation(imdl_mdeit_0,data_b_0.dBz,lambda_vector);
    imdl_mdeit_0.hyperparameter = struct('value',lambda_optimal_mdeit_0);
    lambda_mdeit_0 = lambda_optimal_mdeit_0;
    
    imdl_mdeit_1.fwd_model = fmdl_reconstruction_r_1; %Use a different forward model for reconstruction

    [lambda_optimal_mdeit_1,optimal_id_mdeit_1,V_mu_mdeit_1,~] = generalized_cross_validation(imdl_mdeit_1,data_b_1.dBz,lambda_vector);
    imdl_mdeit_1.hyperparameter = struct('value',lambda_optimal_mdeit_1);
    lambda_mdeit_1 = lambda_optimal_mdeit_1;

    imdl_mdeit_2.fwd_model = fmdl_reconstruction_r_2; %Use a different forward model for reconstruction

    [lambda_optimal_mdeit_2,optimal_id_mdeit_2,V_mu_mdeit_2,~] = generalized_cross_validation(imdl_mdeit_2,data_b_2.dBz,lambda_vector);
    imdl_mdeit_2.hyperparameter = struct('value',lambda_optimal_mdeit_2);
    lambda_mdeit_2 = lambda_optimal_mdeit_2;
    
    % EIT
    imdl_eit.fwd_model = fmdl_reconstruction_r_0; %Use a different forward model for reconstruction
    [lambda_optimal_eit_0,optimal_id_eit_0,V_mu_eit_0,~]  = generalized_cross_validation(imdl_eit,data_u_0.meas_du,lambda_vector);
    % Since EIDORS uses hp^2*RtR, but the generalized_cross_validation function
    % uses hp*RtR, we have to correct for that.
    imdl_eit.hyperparameter = struct('value',sqrt(lambda_optimal_eit_0));
    
    lambda_eit_0 = sqrt(lambda_optimal_eit_0);
    %% Plots
    subplot(1,3,1)

    hold on
    plot(lambda_vector,V_mu_mdeit_0);
    plot(lambda_vector(optimal_id_mdeit_0),V_mu_mdeit_0(optimal_id_mdeit_0),'r.');
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlabel('$\lambda$','Interpreter','latex')
    hold off
    subplot(1,3,2)
    
    hold on
    plot(lambda_vector,V_mu_mdeit_1);
    plot(lambda_vector(optimal_id_mdeit_1),V_mu_mdeit_1(optimal_id_mdeit_1),'r.');
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlabel('$\lambda$','Interpreter','latex')
    hold off
    
    subplot(1,3,3)
    hold on
    plot(lambda_vector,V_mu_mdeit_2);
    plot(lambda_vector(optimal_id_mdeit_2),V_mu_mdeit_2(optimal_id_mdeit_2),'r.');
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    xlabel('$\lambda$','Interpreter','latex')
    hold off
    
    %% Perform reconstruction

    % Reconstruct (difference) for 1-MDEIT
    img_output_mdeit_0 = inv_solve_mdeit(imdl_mdeit_0,Bzh_0,Bzi_0);

    img_output_mdeit_1 = inv_solve_mdeit(imdl_mdeit_1,Bzh_1,Bzi_1);

    img_output_mdeit_2 = inv_solve_mdeit(imdl_mdeit_2,Bzh_2,Bzi_2);

    % Reconstruct (difference) for EIT

    img_output_eit_0 = inv_solve(imdl_eit,uh_0,ui_0);
    
    %% Noise correction
    func_eit = @(imgh,imgi) noisy_data_generator_eit(imgh,imgi,options);
    func_mdeit = @(imgh,imgi) noisy_data_generator_mdeit(imgh,imgi,options);

    Jh_0 = img_output_mdeit_0.jacobian;
    sigma_std_mdeit_0 = noise_correction(imgh_0,imgi_0,Jh_0,lambda_mdeit_0,func_mdeit,num_noise_repetitions);

    Jh_1 = img_output_mdeit_1.jacobian;
    sigma_std_mdeit_1 = noise_correction(imgh_1,imgi_1,Jh_1,lambda_mdeit_1,func_mdeit,num_noise_repetitions);

    Jh_2 = img_output_mdeit_2.jacobian;
    sigma_std_mdeit_2 = noise_correction(imgh_2,imgi_2,Jh_2,lambda_mdeit_2,func_mdeit,num_noise_repetitions);
        
    % I LEFT WORK HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    % DEBUG!!!!!!!!!!!!!!!!!!!!!!

    img_eit_h_0 = mk_image(fmdl_reconstruction_r_0,background_conductivity);
    Jh = calc_jacobian(img_eit_h_0);
    sigma_std_eit_0 = noise_correction(imgh_0,imgi_0,Jh,lambda_eit_0,func_eit,num_noise_repetitions);
    % sigma_std_eit_0 = ones(size(img_eit_h_0.elem_data));

    img_eit_h_1 = mk_image(fmdl_reconstruction_r_1,background_conductivity);
    % Jh = calc_jacobian(img_eit_h_1);
    % sigma_std_eit_1 = noise_correction(imgh_1,imgi_1,Jh,lambda_eit_1,func_eit,num_noise_repetitions);
    sigma_std_eit_1 = ones(size(img_eit_h_1.elem_data));


    pause(1e-20)
    %% Plots (2D)
    half_height = model_parameters.height/2;

    n_points = 128; %Images will be 128x128 pixels, and will correspond to a window of  [-radius,radius]^2 in real space

    % Draw a red circle where the anomaly is supposed to be
    theta = linspace(0, 2*pi, 200);
    pixel_center = [1,1]*0.5*n_points;
    pixel_anomaly_radius = n_points*model_parameters.anomaly.radius/(2*model_parameters.radius);
    pixel_anomaly_x0 = pixel_center(1)+n_points*model_parameters.anomaly.position(1)/(2*model_parameters.radius);
    pixel_anomaly_y0 = pixel_center(2)+n_points*model_parameters.anomaly.position(2)/(2*model_parameters.radius);
    
    x = pixel_anomaly_x0  + pixel_anomaly_radius*cos(theta);
    y = pixel_anomaly_y0 + pixel_anomaly_radius*sin(theta);


    % The slice definition (horizontal at half height)
    level.centre = [0,0,half_height];
    level.normal_angle = [0,0,1];

    figure_name = strcat('2D SNR = ',num2str(SNR_vector(i),'%i'));
    figure('Name', figure_name);
    
    subplot(2,2,1)
    img_temp = img_output_eit_0;
    img_temp.elem_data = img_temp.elem_data./sigma_std_eit_0;
    img_temp.calc_colours.npoints = n_points;
    img_temp.show_slices.do_colourbar = true;
    
    show_slices(img_temp,[inf,inf,1.5])
    hold on
    plot(x, y, 'r', 'LineWidth', 1.5)
    hold off;

    title('EIT','Interpreter','latex')

    subplot(2,2,2)
    img_temp = img_output_mdeit_0;
    img_temp.elem_data = img_temp.elem_data./sigma_std_eit_0;
    img_temp.calc_colours.npoints = n_points;
    img_temp.show_slices.do_colourbar = true;
    img_temp.calc_colours.clim = 100;
    show_slices(img_temp,level)
    hold on
    plot(x, y, 'r', 'LineWidth', 1.5)
    hold off;

    title_str = strcat('MDEIT min sensor radius $\kappa = ',num2str(condition_number_at_x0, '%.1g'),'$');
    title(title_str,'interpreter','latex')
    
    subplot(2,2,3)
    img_temp = img_output_mdeit_1;
    img_temp.elem_data = img_temp.elem_data./sigma_std_mdeit_1;
    img_temp.calc_colours.npoints = n_points;
    img_temp.show_slices.do_colourbar = true;

    show_slices(img_temp,level)
    hold on
    plot(x, y, 'r', 'LineWidth', 1.5)
    hold off;

    title_str = strcat('MDEIT opt radius $\kappa = ',num2str(condition_number_at_ropt_pswarm, '%.1g'),'$');
    title(title_str,'interpreter','latex')

    subplot(2,2,4)
    img_temp = img_output_mdeit_2;
    img_temp.elem_data = img_temp.elem_data./sigma_std_mdeit_2;
    img_temp.calc_colours.npoints = n_points;
    img_temp.show_slices.do_colourbar = true;

    show_slices(img_temp,level)
    hold on
    plot(x, y, 'r', 'LineWidth', 1.5)
    hold off;

    title_str = strcat('MDEIT opt config $\kappa = ',num2str(condition_number_at_xopt, '%.1g'),'$');
    title(title_str,'interpreter','latex')


    %% Plots (3D)
    figure_name = strcat('With NBC, SNR = ',num2str(SNR_vector(i),'%i'));
    figure('Name', figure_name);
    
    subplot(2,2,1)
    img_temp = img_output_eit_0;
    img_temp.elem_data = img_temp.elem_data./sigma_std_eit_0;
    show_fem_transparent_edges(img_temp)
    title('EIT min sensor radius')

    subplot(2,2,2)
    img_temp = img_output_mdeit_0;
    img_temp.elem_data = img_temp.elem_data./sigma_std_mdeit_0;
    show_fem_transparent_edges(img_temp)

    title_str = strcat('MDEIT min sensor radius $\kappa = ',num2str(condition_number_at_x0, '%.1g'),'$');
    title(title_str,'interpreter','latex')


    subplot(2,2,3)
    img_temp = img_output_mdeit_1;
    img_temp.elem_data = img_temp.elem_data./sigma_std_mdeit_1;
    show_fem_transparent_edges(img_temp)
    title_str = strcat('MDEIT opt radius $\kappa = ',num2str(condition_number_at_ropt_pswarm, '%.1g'),'$');
    title(title_str,'interpreter','latex')

    subplot(2,2,4)
    img_temp = img_output_mdeit_2;
    img_temp.elem_data = img_temp.elem_data./sigma_std_mdeit_2;
    show_fem_transparent_edges(img_temp)

    title_str = strcat('MDEIT opt config $\kappa = ',num2str(condition_number_at_xopt, '%.1g'),'$');
    title(title_str,'interpreter','latex')
    
    %% Plots (3D with sensors)
    figure_name = strcat('No NBC, SNR = ',num2str(SNR_vector(i),'%i'));
    figure('Name', figure_name);
    subplot(2,2,1)
    show_fem_transparent_edges(img_output_eit_0,true)
    title('EIT min sensor radius')
    subplot(2,2,2)
    show_fem_transparent_edges(img_output_mdeit_0,true)
    title('MDEIT min sensor radius')

    subplot(2,2,3)
    show_fem_transparent_edges(img_output_mdeit_1,true)
    title('MDEIT opt radius')
    subplot(2,2,4)
    show_fem_transparent_edges(img_output_mdeit_2,true)
    title('MDEIT opt config')

    pause(1e-10)

end

%%
% figure
% hold on
% show_fem(imgh)
% plot3(x0(1:num_of_sensors),x0(num_of_sensors+1:2*num_of_sensors),model_height*ones(num_of_sensors,1),'b.','MarkerSize',5)
% 
% legend_axis = [];
% legend_str = {};
% ct = 0;
% if exist('sensor_locations_opt_pswarm')
%     h = plot3(sensor_locations_opt_pswarm(1:num_of_sensors),sensor_locations_opt_pswarm(num_of_sensors+1:2*num_of_sensors),model_height*ones(num_of_sensors,1),'r.','MarkerSize',10);
%     legend_axis = [legend_axis h];
%     legend_1 = strcat('pswarm, $\kappa = ',num2str(fopt_pswarm),'$');
%     legend_str{ct} =  legend_1;
%     ct = ct+1;
% end
% if exist('sensor_locations_opt_psearch')
%     q = plot3(sensor_locations_opt_psearch(1:num_of_sensors),sensor_locations_opt_psearch(num_of_sensors+1:2*num_of_sensors),model_height*ones(num_of_sensors,1),'y.','MarkerSize',10);
%     legend_axis = [legend_axis q];
%     legend_2 = strcat('psearch, $\kappa = ',num2str(fopt_psearch),'$');
% 
%     legend_str{ct} =  legend_2;
%     ct = ct+1;
% end
% if exist('sensor_locations_opt_ga')
%     p = plot3(sensor_locations_opt_ga(1:num_of_sensors),sensor_locations_opt_ga(num_of_sensors+1:2*num_of_sensors),model_height*ones(num_of_sensors,1),'g.','MarkerSize',10);
%     legend_axis = [legend_axis p];
%     legend_3 = strcat('ga, $\kappa = ',num2str(fopt_ga),'$');
% 
%     legend_str{ct} =  legend_3;
%     ct = ct+1;
% end
% 
% hold off
% 
% legend(legend_axis,legend_str,'Interpreter','latex');

%% Function
function cost = compute_cost(condition_number,r,rmin,rmax)
    cost = log10(condition_number(:)) + (r(:)-rmin).^2/(rmax-rmin)^2;
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