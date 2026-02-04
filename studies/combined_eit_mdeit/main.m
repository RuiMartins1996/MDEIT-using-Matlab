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
H0 = I0/l0; 
B0 = 1.25663706127e-6*I0/l0; 

% Set the minimum sensor radius
sensor_radius_0 = 1.01;
rmin = sensor_radius_0;
rmax = 3;

% Number of noise samples to gather for noise based correction
num_noise_repetitions = 15;

%% Problem parameters

% Default 3D model
model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);
% with some changes
model_parameters.maxsz = model_parameters.radius/15;
model_height = model_parameters.height/2;
model_parameters.numOfElectrodesPerRing = 8;
model_parameters.numOfRings = 4;
model_parameters.numOfSensors = 8*4;
model_parameters.anomaly.type = 'spherical';
model_parameters.material.type = 'spherical';
model_parameters.material.name = 'sphere_anomaly';

% model_parameters.material.position = [0.5 0 model_parameters.height/2];
model_parameters.material.position = [0 0 model_parameters.height/2];
% model_parameters.material.radius = model_parameters.material.radius/2;

background_conductivity = 3.28e-1/sigma0;
anomaly_conductivity = 1e-12/sigma0;

maxsz_reconstruction = model_parameters.maxsz;

current_amplitude = 2.4e-3/I0;

% Set noise type and snr
noise_type = 'white';

snr = 20;

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

%% Find optimal regularization parameter with GCV, and perform tsvd reconstruction
lambda_vector = logspace(-17,3,30);

options = struct('reconstruct',true);

[lambda_optimal_mdeit3,optimal_id_mdeit_1,V_mu_mdeit1,dx_mdeit] = generalized_cross_validation(imdl_mdeit_3,[data_b.dBx; data_b.dBy; data_b.dBz],lambda_vector,'reconstruct',true);
imdl_mdeit_3.hyperparameter = struct('value',lambda_optimal_mdeit3);


[lambda_optimal_eit,optimal_id_eit,V_mu_eit,dx_eit] = generalized_cross_validation(imdl_eit,data_u.meas_du,lambda_vector,'reconstruct',true);
% Since EIDORS uses hp^2*RtR, but the generalized_cross_validation function
% uses hp*RtR, we have to correct for that.
imdl_eit.hyperparameter = struct('value',sqrt(lambda_optimal_eit));

%% Show results
img_tsvd_eit =  mk_image(imdl_eit.fwd_model,imdl_eit.jacobian_bkgnd.value);
img_tsvd_eit.elem_data = dx_eit;


img_tsvd_mdeit =  mk_image(imdl_eit.fwd_model,imdl_eit.jacobian_bkgnd.value);
img_tsvd_mdeit.elem_data = dx_mdeit;

figure
subplot(1,2,1)
show_fem_transparent_edges(img_tsvd_mdeit);
subplot(1,2,2)
show_fem_transparent_edges(img_tsvd_eit);

drawnow

%% Perform reconstruction (DO THIS ANYWAY TO COMPUTE THE JACOBIANS. CHANGE LATER!!!)

% Reconstruct (difference) for 1-MDEIT
Bh = [Bxh_signal;Byh_signal;Bzh_signal];
Bi_noisy = [Bxi_noisy;Byi_noisy;Bzi_noisy];

imgh_r = mk_image_mdeit(fmdl_reconstruction,background_conductivity);
A = @(sigma) M(imgh_r,sigma);
J_mdeit = calc_jacobian_mdeit(imgh_r, imgh_r.elem_data,[],A,'mdeit3');

J_eit = calc_jacobian(imgh_r);



% img_output_mdeit_3 = inv_solve_mdeit(imdl_mdeit_3,Bh,Bi_noisy);
% J_mdeit = img_output_mdeit_3.jacobian;

% Reconstruct (difference) for EIT
% img_output_eit = inv_solve(imdl_eit,uh_signal,ui_noisy);

%% Try to reconstruct with both information
fprintf('Attempting EIT+MDEIT reconstruction: \n');

% Functional to minimize will be :
% 0.5||r_eit||^2+beta/2||r_mdeit||^2+lambda/2||dsigma||^2

% lambda is the regularization parameter and beta is a weighting factor, to
% make both residuals of the same order. A good choice might be beta =
% V0/H0;

% Normal equations will be: (J_eit'*J_eit' + beta*J_mdeit'*J_mdeit +
% lambda*eye)*dsigma = -J_eit*r_eit - beta*J_mdeit*r_mdeit

% Residual for both EIT and MDEIT:
residual_eit =  uh_signal.meas(:)-ui_noisy.meas(:);
residual_mdeit =  Bh(:)-Bi_noisy(:);

beta = V0/H0;

% Matrix aglomerating both jacobians  
% A_aug = [J_eit; sqrt(beta) * J_mdeit; sqrt(lambda) * eye(n_elem)]; 
% r_aug = [residual_eit; sqrt(beta) * residual_mdeit; zeros(n_elem,1)];  % n = number of unknowns

A_aug = [J_eit; sqrt(beta) * J_mdeit]; 
r_aug = [residual_eit; sqrt(beta) * residual_mdeit];  % n = number of unknowns

% Generalized cross-validation
fprintf('Doing SVD: \n')
tic
[U,S,V] = svd(A_aug,'econ');

fprintf('Took %i (s)\n',toc);

%% Generalized cross validation
% sigma = diag(S);             % singular values
% Uy = U' * r_aug;              % coordinates of data in U-basis
% m = size(A_aug,1);
% n = length(lambda_vector);
% 
% V_lambda = zeros(n,1);
% for i = 1:n
%     lambda = lambda_vector(i);
% 
%     gamma = sigma.^2 ./ (sigma.^2 + m*lambda);    % shrinkage factors
% 
%     one_minus_gamma = 1 - gamma;
% 
%     numerator = (1/m) * sum( (one_minus_gamma .* Uy).^2 );
% 
%     denominator = ( (1/m) * sum(one_minus_gamma) )^2;
% 
%     V_lambda(i) = numerator / denominator;
% end

% %version without m*lambda, I think it is wrong tho
% V_lambda = zeros(n,1);
% for i = 1:n
%     lambda = lambda_vector(i);
%     one_minus_gamma = 1 - sigma.^2 ./ (sigma.^2 + lambda);
%     V_lambda(i) = sum((one_minus_gamma .* Uy).^2) / (sum(one_minus_gamma)^2);
% end
% 
% optimal_id = find(V_lambda == min(V_lambda));
% 
% lambda_optimal_combined = m*lambda_vector(optimal_id);

%% 
% figure
% hold on
% plot(lambda_vector,V_lambda);
% plot(lambda_vector(optimal_id),V_lambda(optimal_id),'r.','MarkerSize',10)
% hold off
% grid on;grid minor;box on;
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% xlabel('$\lambda$','Interpreter','latex')
% ylabel('$V(\lambda)$','Interpreter','latex')

%% L-curve method (for some reason, works better than GCV here)

for i = 1:length(lambda_vector)
    
    lambda = lambda_vector(i);

    sigma = diag(S);             % singular values
    Uy = U' * r_aug;             % coordinates in U-basis
    z = (sigma ./ (sigma.^2 + lambda)) .* Uy;  % correct if sigma and Uy are same length
    d_sigma = V * z;

    residual_norms(i) = norm(A_aug*d_sigma-r_aug,2); 
    x_norms(i) = norm(d_sigma,2);
end

% Given data
r = residual_norms(:);
xnorm = x_norms(:);

% Sort by residual norm (important for monotone x)
[rs, idx] = sort(r);
xs = xnorm(idx);

% Log-log scale for L-curve
lr = log(rs);
lx = log(xs);

% --- Use monotone-preserving piecewise cubic interpolation ---
lr_dense = linspace(min(lr), max(lr), 200);
lx_dense = interp1(lr, lx, lr_dense, 'pchip');

% --- Numerical derivatives using finite differences ---
d1 = gradient(lx_dense, lr_dense);          % first derivative
d2 = gradient(d1, lr_dense);               % second derivative

% --- Curvature formula: κ = |y''| / (1 + y'^2)^(3/2)
kappa = abs(d2) ./ (1 + d1.^2).^(3/2);

% --- Interpolate curvature back to original lr points if needed
kappa_interp = interp1(lr_dense, kappa, lr, 'pchip');

% Find maximum curvature (corner of L-curve)
[~, imax] = max(kappa_interp);
opt_lr = lr(imax);
opt_lx = lx(imax);

% Assuming lambda_vector corresponds to lr_dense
lambda_opt = lambda_vector(imax);

% Convert back to linear scale
opt_r = exp(opt_lr);
opt_x = exp(opt_lx);

fprintf('Optimal residual norm = %.4e\n', opt_r);
fprintf('Optimal solution norm = %.4e\n', opt_x);
fprintf('Maximum curvature = %.4e\n', kappa_interp(imax));
fprintf('Optimal hyperparameter = %.4e\n', lambda_opt);

%% figure
figure
hold on
plot(exp(lr_dense),exp(lx_dense),'b.');
plot(residual_norms,x_norms);
plot(opt_r,opt_x,'r.','MarkerSize',10)
text(opt_r+0.1,opt_x,strcat('$\lambda = ',num2str(lambda_opt),'$'),'Interpreter','latex')
hold off
grid on;grid minor;box on;
set(gca,'YScale','log')
set(gca,'XScale','log')


%% Solve with TSVD
lambda =  lambda_opt;

sigma = diag(S);             % singular values
Uy = U' * r_aug;             % coordinates in U-basis
z = (sigma ./ (sigma.^2 + lambda)) .* Uy;  % correct if sigma and Uy are same length
dx = V * z;

img_combined = imgh_r;
img_combined.elem_data = dx;

%% Noise correction

function data_noisy = noisy_data_generator_mdeit(imgh,imgi,options)
[~,~,~,noise_b,~] =  generate_data(options,imgh,imgi);
noise_b_vec = [noise_b.noise_b_x;noise_b.noise_b_y;noise_b.noise_b_z];
data_noisy = noise_b_vec;
end

function data_noisy = noisy_data_generator_eit(imgh,imgi,options)
[~,~,~,~,noise_u] =  generate_data(options,imgh,imgi);
data_noisy = noise_u;
end

function data_noisy = noisy_data_generator_eit_mdeit(imgh,imgi,options)
assert(isfield(options,'V0H0'),'Need field');
V0H0 = options.V0H0;

[~,~,~,noise_b,noise_u] =  generate_data(options,imgh,imgi);
noise_b_vec = [noise_b.noise_b_x;noise_b.noise_b_y;noise_b.noise_b_z];
data_noisy = [noise_u;V0H0*noise_b_vec];
end


options.V0H0 = V0/H0;
options.noise_structure.type = 'white';
options.noise_structure.snr = snr;
options.noise_structure.epoch_time = 50e-3; %50 milliseconds
options.noise_structure.sampling_rate = 1000;
options.noise_structure.B0 = 1.2566370612720e-6*J0*l0;
options.noise_structure.measurement_protocol = 'default_rui';

func_eit = @(imgh,imgi) noisy_data_generator_eit(imgh,imgi,options);
func_mdeit = @(imgh,imgi) noisy_data_generator_mdeit(imgh,imgi,options);
funct_eit_mdeit = @(imgh,imgi) noisy_data_generator_eit_mdeit(imgh,imgi,options);

lambda_eit = imdl_eit.hyperparameter.value^2;
lambda_mdeit =imdl_mdeit_3.hyperparameter.value;

sigma_std_eit = noise_correction(imgh,imgi,J_eit,lambda_eit,func_eit,15);
sigma_std_mdeit = noise_correction(imgh,imgi,J_mdeit,lambda_mdeit,func_mdeit,15);
sigma_std_eit_mdeit = noise_correction(imgh,imgi,A_aug,lambda_opt,funct_eit_mdeit,15);


%%
figure
subplot(2,3,1)
show_fem_transparent_edges(img_tsvd_mdeit);
title('MDEIT')
subplot(2,3,2)
show_fem_transparent_edges(img_tsvd_eit);
title('EIT')
subplot(2,3,3)
show_fem_transparent_edges(img_combined);
title('MDEIT+EIT')

img_eit_n = img_tsvd_eit;
img_eit_n.elem_data = img_eit_n.elem_data./sigma_std_eit;
img_mdeit_n = img_tsvd_mdeit;
img_mdeit_n.elem_data = img_mdeit_n.elem_data./sigma_std_mdeit;
img_combined_n = img_combined;
img_combined_n.elem_data = img_combined_n.elem_data./sigma_std_eit_mdeit;

subplot(2,3,4)
show_fem_transparent_edges(img_mdeit_n);
title('MDEIT')
subplot(2,3,5)
show_fem_transparent_edges(img_eit_n);
title('EIT')
subplot(2,3,6)
show_fem_transparent_edges(img_combined_n);
title('MDEIT+EIT')
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