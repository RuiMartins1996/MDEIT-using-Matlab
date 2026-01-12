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

background_conductivity = 3.28e-1/sigma0;

%%  Noise
SNRdb = 50;
num_noise_repetitions = 30;

%% Build/load multiple forward models of different sensor radius
model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

model_parameters.numOfRings = 4;
model_parameters.numOfElectrodesPerRing = 4;
model_parameters.numOfSensors = 4*4;
model_parameters.sensorRadius = model_parameters.radius*1.1;

anomaly_conductivity = 5*background_conductivity;
anomaly_position = [0.3 0 model_parameters.height/2];
anomaly_radius = 0.3*model_parameters.radius;

model_parameters.material = struct( ...
    'type', 'spherical', ...
    'name', 'sphere_anomaly', ...
    'radius', anomaly_radius, ...
    'position', anomaly_position);
model_parameters.anomaly = struct(...
    'type','spherical',...
    'conductivity',anomaly_conductivity,...
    'radius',anomaly_radius);

% Build models
[model_parameters,fmdl_array] = mk_mdeit_model(model_parameters,model_folder,[]);

%% Set anonymous functions needed by calc_jacobian_mdeit
imgh = mk_image_mdeit(fmdl_array{1},1.0);

lambdatimesdAdp = @(lambda) computeLambdaTimesDaDp(imgh,lambda);
A = @(sigma) M(imgh,sigma);

% Make homogeneous image
imgh = mk_image_mdeit(fmdl_array{1},1.0);
data1 = fwd_solve_mdeit(imgh);

% Place anomaly
imgi = add_material_properties(imgh,[background_conductivity,anomaly_conductivity]);
data2 = fwd_solve_mdeit(imgi);

figure;
show_fem(imgi);

pause(1e-10);

%% Make inverse model to test leave one out cross validation

imdl= eidors_obj('inv_model','my_inverse_model');

imdl.fwd_model = fmdl_array{1};
imdl.jacobian_bkgnd = struct('value',1.0);
imdl.solver = 'gn';
imdl.hyperparameter = struct('value',1e-6);
imdl.RtR_prior = @(x,Jx) speye(numel(x)); %tykhonov

imdl.recon_mode = 'mdeit1';
imdl.recon_type = 'difference';
imdl.select_sensor_axis = 1;

%% Run generalized cross validation ( works for linearized GN ) 

mu_vector = logspace(-20,-3,50);

function data_noisy = add_measurement_noise(data,SNRdb)
assert(isnumeric(data) && isvector(data),'data must be a numerical vector');
assert(size(data,1)>=size(data,2),'data must be a column vector');

if isempty(SNRdb)
    data_noisy = data;
    return;
else
    noise_level= std(data)/10^(SNRdb/20);
    data_noisy = data+noise_level*randn(size(data));
end
end

if strcmp(imdl.recon_mode,'mdeit3')
    datai = [data2.Bx(:);data2.By(:);data2.Bz(:)];
    datah = [data1.Bx(:);data1.By(:);data1.Bz(:)];

    datai_noisy = add_measurement_noise(datai,SNRdb);
    datah_noisy = add_measurement_noise(datah,SNRdb);

    data = datai_noisy-datah_noisy;
elseif strcmp(imdl.recon_mode,'mdeit1')
    datai = data2.Bx(:);
    datah = data1.Bx(:);

    datai_noisy = add_measurement_noise(datai,SNRdb);
    datah_noisy = add_measurement_noise(datah,SNRdb);

    data = datai_noisy-datah_noisy;
else
    error('here')
end

opts.reconstruct = true;
[mu_min,optimal_id,V_mu,dx] = generalized_cross_validation(imdl,data,mu_vector,'reconstruct',true);

figure
img_output = mk_image_mdeit(imdl.fwd_model,1);
img_output.elem_data = dx;
show_fem(img_output);
%% Plots
figure
hold on
plot(mu_vector,V_mu,'b.');
plot(mu_vector(optimal_id),V_mu(optimal_id),'r.','MarkerSize',15)
hold off
set(gca,'YScale','log');
set(gca,'XScale','log');

xlabel('$\lambda$','Interpreter','latex')
ylabel('$V(\lambda)$','Interpreter','latex')
grid on; grid minor;
legend('$V(\lambda)$','interpreter','latex')

%% Run leave-one-out cross validation
% leave_one_out_CV(imdl,data.Bx(:),1e-5)

%% Save results

model_struct = struct();
model_struct.imdl = imdl;
model_struct.imgi = imgi;
model_struct.imgh = imgh;
model_struct.model_parameters = model_parameters;


gcv_struct = struct();
gcv_struct.V_mu = V_mu;
gcv_struct.mu_min = mu_min;
gcv_struct.mu_vector = mu_vector;
gcv_struct.conductivity_optimal_reg_par = dx;
gcv_struct.optimal_id = optimal_id;

loocv_struct = struct();

file_name = create_file_name(data_folder,model_parameters,SNRdb,imdl.recon_mode);
save(file_name,'model_struct','gcv_struct','loocv_struct');


%% FUNCTIONS: Leave one out
function leave_one_out_CV(imdl,data,lambda)
%% Check if inputs are valid
assert(isfield(imdl,'recon_mode'),'No imdl.recon_mode field');
recon_mode = imdl.recon_mode;

allowed_recon_modes = {'mdeit1','mdeit3'};
assert(ismember(recon_mode,allowed_recon_modes),...
    'Invalid reconstruction mode: %s. Allowed modes are: %s', ...
    recon_mode, strjoin(allowed_recon_modes, ', '))

num_sensors = numel(imdl.fwd_model.sensors);
num_stim = numel(imdl.fwd_model.stimulation);

% if strcmp(recon_mode,'eit')
%     assert(isfield(imdl.fwd_model.stimulation.meas_pattern),'No fwd_model.stimulation.meas_pattern defined');
% 
%     num_measurements = 0;
%     for i = 1:num_stim
%         num_measurements = num_measurements + size(imdl.fwd_model.stimulation(i).meas_pattern,1);
%     end
% 
%     assert(all(size(data)==[num_measurements,1]),'Size of data is not consistent with imdl.recon_mode')
% end

if strcmp(recon_mode,'mdeit1')
    assert(all(size(data)==[num_stim*num_sensors,1]),'Size of data is not consistent with imdl.recon_mode')
    
    assert(isfield(imdl,'select_sensor_axis'),'Need to define imdl.select_sensor_axis for this reconstruction mode');
    
    select_sensor_axis = imdl.select_sensor_axis;
    
elseif strcmp(recon_mode,'mdeit3')
    assert(all(size(data)==[3*num_stim*num_sensors,1]),'Size of data is not consistent with imdl.recon_mode')
end

assert(isfield(imdl,'RtR_prior'),'imdl.RtR_prior must be defined');
RtR = imdl.RtR_prior;

n_elem = size(imdl.fwd_model.elems,1);
tol = 1e-6;
max_iterations = 5;

%% Setup inverse solve 
% have to do this since its incovenient to use inv_solve_mdeit while leaving one data point out


% a dummy image is enough, since its elem_data is change inside these
% functions or they just need other img info
img = mk_image_mdeit(imdl.fwd_model,ones(n_elem,1));

lambdatimesdAdp = @(lambda) computeLambdaTimesDaDp(img,lambda);
A = @(sigma) M(img,sigma);

% Define the residual and jacobian anonymous functions

res = @(x) calc_residual_mdeit(img, x, data,recon_mode,select_sensor_axis);
jac = @(x) calc_jacobian_mdeit(img, x,lambdatimesdAdp,A,recon_mode,select_sensor_axis);

%% Do leave-one-out cross validation

x0 = ones(n_elem,1);

cross_validation_error = zeros(length(lambda),1);

for n = 1:length(lambda)
    
    e = zeros(length(data),1);

    for i = 1:length(data)
        % residual and jacobian functions that leave one out in the data
        % (equivalent to leave one out in the data)
        res_i = @(x) res_leave_one_out(x,res,i);
        jac_i = @(x) jac_leave_one_out(x,jac,i);

        % solve the inverse problem
        x_i = gn_solve(res_i, jac_i, x0, tol, max_iterations, RtR, lambda(n));
        
        % predict the left out measurement
        img = mk_image_mdeit(imdl.fwd_model,x_i);
        data_i = fwd_solve_mdeit(img);
        
        if strcmp(recon_mode,'mdeit1')
            switch select_sensor_axis
                case 1
                    y_i = data_i.Bx(:);
                case 2
                    y_i = data_i.By(:);
                case 3
                    y_i = data_i.Bz(:);
            end
            e(i) = (y_i(i)-data(i))^2;
        else
            y_i = [data_i.Bx(:);data_i.By(:);data_i.Bz(:)];
            e(i) = (y_i(i)-data(i))^2;
        end
    end

    cross_validation_error(n) = avg(e);
end

end

%% FUNCTIONS
function r = res_leave_one_out(x,res,i)
r = res(x);
r(i) = [];
end

%% FUNCTIONS
function j = jac_leave_one_out(x,jac,i)
j = jac(x);
j(i,:) = [];
end

%% FUNCTIONS
function file_name = create_file_name(data_folder,model_parameters,SNRdb,recon_mode)
    if isempty(SNRdb)
        name = sprintf('data_E_%i_R_%i_M_%i_mode_%s_noiseless',...
        model_parameters.numOfElectrodesPerRing,...
        model_parameters.numOfRings,...
        model_parameters.numOfSensors, ...
        recon_mode);
    else
    name = sprintf('data_E_%i_R_%i_M_%i_mode_%s_SNR_%d',...
        model_parameters.numOfElectrodesPerRing,...
        model_parameters.numOfRings,...
        model_parameters.numOfSensors,...
        recon_mode,...
        SNRdb);
    end
    file_name = strcat(data_folder,'\',name);
end