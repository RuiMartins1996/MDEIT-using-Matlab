clc;clear all;close all

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

rng(1);

%% First of all, lets find the number of measurements that maximizes the number of discared singular values

% Do this by hand


%% Model parameters 
z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);

model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

% Geometry
model_parameters.height = 3;

% Electrodes
model_parameters.numOfElectrodesPerRing = 6;
model_parameters.numOfRings = 4;
% Sensors
model_parameters.numOfSensors = model_parameters.numOfElectrodesPerRing*model_parameters.numOfRings;

model_parameters.isCylindrical = false;
model_parameters.sensorRadius = 3.5*model_parameters.radius;

% CHANGED HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% model_parameters.isCylindrical = true;
% model_parameters.sensorRadius = 1.5*model_parameters.radius;



% Material/anomaly
model_parameters.material.type = 'spherical';
% model_parameters.material.position = [0*model_parameters.radius 0 0.5*model_parameters.height];

model_parameters.material.position = [0.5*model_parameters.radius 0 0.5*model_parameters.height];

% model_parameters.material.position(1) = 0.95*(model_parameters.radius-model_parameters.material.radius);
% model_parameters.material.position(3) = 0.5*model_parameters.height;


background_conductivity = 3.28e-1/sigma0;
anomaly_conductivity = background_conductivity*1.1;
current_amplitude = 2.4e-3/I0;

minsz_mesh_convergence = 0.1;
maxsz_mesh_convergence = 0.1;
num_meshes_mesh_convergence = 1;

maxsz_reconstruction = 0.1;


inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used
options = {};

num_noise_repetitions = 30;
SNRdb = 10;
%% Simulation parameters
dim = 3;            %if doing 1-axis, which dimmension?

rmax = 2*model_parameters.radius;
rmin = model_parameters.radius*1.01;

zmax = model_parameters.height;
zmin = 0;

R0 = model_parameters.radius*1.5;

%% Mesh convergence study to obtain the "correct" forward model
model_parameters.maxsz = minsz_mesh_convergence;
[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);

fmdl = fmdls{1};
%% Assign the stimulation pattern to this forward model

stimulation = mk_stim_patterns(numel(fmdl.electrode),1,inj,meas,options,current_amplitude);
fmdl.stimulation = stimulation;
for i = 1:numel(fmdl.electrode)
    fmdl.electrode(i).z_contact = z0/z0;
end

%% Make images from forward models and forward solve

% Make homogeneous image
imgh = mk_image_mdeit(fmdl,background_conductivity);

% Add plastic cylinder
imgi = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

figure
show_fem(imgi);
drawnow
%% Generate data from these forward models

% Forward solve
[datah,uh] = fwd_solve_mdeit(imgh);
Bxh = datah.Bx(:);
Byh = datah.By(:);
Bzh = datah.Bz(:);
Bh = [datah.Bx(:);datah.By(:);datah.Bz(:)];

% Forward solve
[datai,ui] = fwd_solve_mdeit(imgi);
Bxi = datai.Bx(:);
Byi = datai.By(:);
Bzi = datai.Bz(:);
Bi = [datai.Bx(:);datai.By(:);datai.Bz(:)];

dB = add_measurement_noise_difference(Bi,Bh,SNRdb);
dBx = add_measurement_noise_difference(Bxi,Bxh,SNRdb);
dBy = add_measurement_noise_difference(Byi,Byh,SNRdb);
dBz = add_measurement_noise_difference(Bzi,Bzh,SNRdb);

du = add_measurement_noise_difference(ui.meas,uh.meas,SNRdb);


%% Generate coarse forward model for reconstruction (different mesh than the data)
model_parameters.material = struct();
model_parameters.maxsz = maxsz_reconstruction;

[model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
fmdl_reconstruction = fmdls{1};

% Don't forget to assignt the same stimulation pattern
fmdl_reconstruction.stimulation = stimulation;

n_elem = size(fmdl_reconstruction.elems,1);

%% 

imgh_reconstruction = mk_image_mdeit(fmdl_reconstruction,background_conductivity);

A = @(x) M(imgh_reconstruction,x);

[Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(imgh_reconstruction,A);

J = [Jx;Jy;Jz];

J_eit = calc_jacobian(imgh_reconstruction);

%% Singular value decomposition of Jacobians
fprintf('Computing SVDs\n')

% fprintf('SVD mdeit-x\n')
% [Ux,Sx,Vx] = svds(Jx,size(Jx,1));
% s_x = diag(Sx);
% 
% fprintf('SVD mdeit-y\n')
% [Uy,Sy,Vy] = svds(Jy,size(Jy,1));
% s_y = diag(Sy);

fprintf('SVD mdeit-z\n')
[Uz,Sz,Vz] = svds(Jz,size(Jz,1));
s_z = diag(Sz);

fprintf('SVD mdeit-3\n')
[U,S,V] = svds(J,size(J,1));
s_xyz = diag(S);

fprintf('SVD eit\n')
[U_eit,S_eit,V_eit] = svds(J_eit,size(J_eit,1));
s_eit = diag(S_eit);


%%
% rank_x = sum(s_x>1e-10);
% rank_y = sum(s_y>1e-10);
rank_z = sum(s_z>1e-10);
rank_xyz = sum(s_xyz>1e-10);
rank_eit =  sum(s_eit>1e-10);

if size(J,2) == rank_xyz
    warning('Number of elements is too small for the current number of measurements! The rank is the number of elements!')
end

% condition_number_x = s_x(1)/s_x(rank_x);
% condition_number_y = s_y(1)/s_y(rank_y);
condition_number_z = s_z(1)/s_z(rank_z);
condition_number_xyz = s_xyz(1)/s_xyz(rank_xyz);
condition_number_eit = s_eit(1)/s_eit(rank_eit);

%% Truncate the singular values in multiple ways:

ns = length(s_xyz);

% Find the singular value of s_xyz closest to the the least singular value
% of s_z
id1 = find(abs(s_xyz-s_z(rank_z)) == min(abs(s_xyz-s_z(rank_z))));

% Find the index of the sv of mdeit3 such that 50% of them are discarded
id2 = truncate_at_percentage(s_xyz,50);

% Find the index of the sv of mdeit3 such that 75% of them are discarded
id3 = truncate_at_percentage(s_xyz,75);



function id = truncate_at_percentage(s_mdeit3,target_percentange)

get_rank = @(s) sum(s>1e-10);

rank_mdeit3 = get_rank(s_mdeit3);

left_limit = 1;
right_limit = rank_mdeit3;

% Use bissection method 
id = (left_limit-1)+floor((right_limit-left_limit)/2);

while 1
    s_truncated = s_mdeit3(1:id);

    rank_truncated = get_rank(s_truncated);
    
    num_of_discarded_sv = rank_mdeit3-rank_truncated;

    rank_percentage = 100*(num_of_discarded_sv)/(rank_mdeit3);

    if rank_percentage >= target_percentange
        left_limit = id;
    else
        right_limit = id;
    end

    id_new = (left_limit-1) + floor((right_limit-left_limit)/2);
    
    %fprintf('id: %i,id_new: %i\n',id,id_new);

    % Stopping criterion
    if abs(id_new-id) == 1 || abs(id_new-id) == 0
        break
    end

    id = id_new;

end

end

%% Plot the truncation lines
figure
hold on
plot(s_z,'b*','MarkerSize',1);
plot(s_xyz,'.');
% Truncation line
plot(1:ns,s_xyz(id1)*ones(ns,1),'b--');
plot(1:ns,s_xyz(id2)*ones(ns,1),'r--');
plot(1:ns,s_xyz(id3)*ones(ns,1),'--','Color',[52,87,49]/256);
hold off
set(gca,'YScale','log');
legend('mdeit-z svs','mdeit-3 svs','Least mdeit-z sv','Discard 50% mdeit3 sv','Discard 75% mdeit3 sv')
ylabel('$\sigma$','Interpreter','latex');

axis([1 ns 1e-9 1e-2])
grid on;grid minor;box on;
drawnow

%% Truncate the singular values at ids

% Truncate at id corresponding to least singular value of z-axis
Un1 = U(:,1:id1);
Vn1 = V(:,1:id1);
Sn1 = S(1:id1,1:id1);

J_xyz_truncated_1 = Un1*Sn1*Vn1.';

s_xyz_truncated_1 = diag(Sn1);
rank_xyz_truncated_1 = sum(s_xyz_truncated_1>1e-10);

condition_number_xyz_truncated = s_xyz_truncated_1(1)/s_xyz_truncated_1(rank_xyz_truncated_1);

fprintf('Discarded SV MDEIT3 (id = %i):\n\t %i (%2.2f %%) \n',id1,rank_xyz-rank_xyz_truncated_1,100*(rank_xyz-rank_xyz_truncated_1)/(rank_xyz));

% Truncate at id corresponding to 50% of discarded SV for mdeit3
Un2 = U(:,1:id2);
Vn2 = V(:,1:id2);
Sn2 = S(1:id2,1:id2);

J_xyz_truncated_2 = Un2*Sn2*Vn2.';

s_xyz_truncated_2 = diag(Sn2);
rank_xyz_truncated_2 = sum(s_xyz_truncated_2>1e-10);
condition_number_xyz_truncated = s_xyz_truncated_2(1)/s_xyz_truncated_2(rank_xyz_truncated_2);

fprintf('Discarded SV MDEIT3 (id = %i):\n\t %i (%2.2f %%) \n',id2,rank_xyz-rank_xyz_truncated_2,100*(rank_xyz-rank_xyz_truncated_2)/(rank_xyz));
% Find corresponding id for z-axis mdeit
id2z = find(abs(s_xyz_truncated_2(end)-s_z) == min(abs(s_xyz_truncated_2(end)-s_z)));

Un2z = Uz(:,1:id2z);
Vn2z = Vz(:,1:id2z);
Sn2z = Sz(1:id2z,1:id2z);

J_z_truncated_2 =  Un2z*Sn2z*Vn2z.';

s_z_truncated_2 = diag(Sn2z);

rank_z_truncated_2 = sum(s_z_truncated_2>1e-10);
condition_number_z_truncated_2 = s_z_truncated_2(1)/s_z_truncated_2(rank_z_truncated_2);

fprintf('Discarded SV MDEIT1 (id = %i):\n\t %i (%2.2f %%) \n',id2z,rank_z-rank_z_truncated_2,100*(rank_z-rank_z_truncated_2)/(rank_z));


% Truncate at id corresponding to 75% of discarded SV for mdeit3
Un3 = U(:,1:id3);
Vn3 = V(:,1:id3);
Sn3 = S(1:id3,1:id3);

J_xyz_truncated_3 = Un3*Sn3*Vn3.';

s_xyz_truncated_3 = diag(Sn3);
rank_xyz_truncated_3 = sum(s_xyz_truncated_3>1e-10);
condition_number_xyz_truncated_3 = s_xyz_truncated_3(1)/s_xyz_truncated_3(rank_xyz_truncated_3);

fprintf('Discarded SV MDEIT3 (id = %i):\n\t %i (%2.2f %%) \n',id3,rank_xyz-rank_xyz_truncated_3,100*(rank_xyz-rank_xyz_truncated_3)/(rank_xyz));
% Find corresponding id for z-axis mdeit

id3z = find(abs(s_xyz_truncated_3(end)-s_z) == min(abs(s_xyz_truncated_3(end)-s_z)));

Un3z = Uz(:,1:id3z);
Vn3z = Vz(:,1:id3z);
Sn3z = Sz(1:id3z,1:id3z);

J_z_truncated_3 =  Un3z*Sn3z*Vn3z.';

s_z_truncated_3 = diag(Sn3z);

rank_z_truncated_3 = sum(s_z_truncated_3>1e-10);
condition_number_z_truncated_3 = s_z_truncated_3(1)/s_z_truncated_3(rank_z_truncated_3);

fprintf('Discarded SV MDEIT1 (id = %i):\n\t %i (%2.2f %%) \n',id3z,rank_z-rank_z_truncated_3,100*(rank_z-rank_z_truncated_3)/(rank_z));

%% Doing L-curve method

lambda_vector = logspace(-10,0,10);

fprintf('(SKIP) L-curve method\n')

% fprintf('Mdeit-x\n')
% [lambda_opt_mdeit_x,sigma_mdeit_x_lc] =...
%     l_curve_method(Jx,dBx,lambda_vector,Ux,Sy,Vx);
% fprintf('Mdeit-y\n')
% [lambda_opt_mdeit_y,sigma_mdeit_y_lc] =...
%     l_curve_method(Jy,dBy,lambda_vector,Uy,Sy,Vy);
% fprintf('Mdeit-z\n')
% [lambda_opt_mdeit_z,sigma_mdeit_z_lc] =...
%     l_curve_method(Jz,dBz,lambda_vector,Uz,Sz,Vz);
% fprintf('Mdeit-3 truncated\n')
% [lambda_opt_mdeit_3_truncated,sigma_mdeit_3_truncated_lc] =...
%     l_curve_method(J_truncated,dB,lambda_vector,Un,Sn,Vn);
% fprintf('Mdeit-3\n')
% [lambda_opt_mdeit_3,sigma_mdeit_3_lc] =...
%     l_curve_method(J,dB,lambda_vector,U,S,V);

function [lambda_opt,sigma] = l_curve_method(J,data,lambda_vector,U,S,V)

res = @(x,y) norm(J*x-y,2); 

% Allocate for residual and solution norms
residual_norms = zeros(1,length(lambda_vector));
x_norms = zeros(1,length(lambda_vector));

for i = 1:length(lambda_vector)
    
    % Solve using TSVD
    temp = diag(S)+lambda_vector(i)./diag(S);
    sigma = V*diag(1./temp)*U'*data;
    
    % Compute residual norm and solution norm
    residual_norms(i) = res(sigma,data);
    x_norms(i) = norm(sigma,2);
end

%% Compute the curvature by fitting a smoothing cubic spline to graph

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

temp = diag(S)+lambda_opt./diag(S);
sigma = V*diag(1./temp)*U'*data;

% fprintf('Optimal residual norm = %.4e\n', opt_r);
% fprintf('Optimal solution norm = %.4e\n', opt_x);
% fprintf('Maximum curvature = %.4e\n', kappa_interp(imax));
% fprintf('Optimal hyperparameter = %.4e\n', lambda_opt);

% %%
% hold on
% plot(xs,rs);
% plot(opt_x,opt_r,'r.')
% text(opt_x,opt_r,sprintf('lambda = %d',lambda_opt))
% hold off

end
%% Cross validation to find optimal hyperparameter and reconstruct with TSVD

fprintf('Doing cross validation\n')

function [lambda_opt,sigma] = generalized_cross_validation(J,data,lambda_vector,U,S,V)

m = size(J,1);
n = length(lambda_vector);

if nargin <4
    [U,S,V] = svds(J,m);
else
    assert(exist('U','var'));
    assert(exist('S','var'))
    assert(exist('V','var'))
end

sigma = diag(S);             % singular values
Uy = U' * data;              % coordinates of data in U-basis

V_lambda = zeros(n,1);

for i = 1:n
    lambda = lambda_vector(i);

    gamma = sigma.^2 ./ (sigma.^2 + m*lambda);    % shrinkage factors

    one_minus_gamma = 1 - gamma;

    numerator = (1/m) * sum( (one_minus_gamma .* Uy).^2 );

    denominator = ( (1/m) * sum(one_minus_gamma) )^2;

    V_lambda(i) = numerator / denominator;
end

optimal_id = find(V_lambda == min(V_lambda));

lambda_opt = m*lambda_vector(optimal_id);

temp = diag(S)+lambda_opt./diag(S);
sigma = V*diag(1./temp)*U'*data;

end


% fprintf('Mdeit-x\n')
% [lambda_opt_mdeit_x,sigma_mdeit_x_gcv] =...
%     generalized_cross_validation(Jx,dBx,lambda_vector,Ux,Sx,Vx);
% fprintf('Mdeit-y\n')
% [lambda_opt_mdeit_y,sigma_mdeit_y_gcv] =...
%     generalized_cross_validation(Jy,dBy,lambda_vector,Uy,Sy,Vy);

fprintf('Mdeit-z\n')
[lambda_opt_mdeit_z,sigma_mdeit_z_gcv] =...
    generalized_cross_validation(Jz,dBz,lambda_vector,Uz,Sz,Vz);

fprintf('Mdeit-1 truncated (id2z)\n')
[lambda_opt_mdeit_z_truncated_2,sigma_mdeit_z_truncated_gcv_2] =...
    generalized_cross_validation(J_z_truncated_2,dBz,lambda_vector,Un2z,Sn2z,Vn2z);

fprintf('Mdeit-1 truncated (id3z)\n')
[lambda_opt_mdeit_z_truncated_3,sigma_mdeit_z_truncated_gcv_3] =...
    generalized_cross_validation(J_z_truncated_3,dBz,lambda_vector,Un3z,Sn3z,Vn3z);

fprintf('Mdeit-3 truncated (id1)\n')
[lambda_opt_mdeit_3_truncated_1,sigma_mdeit_3_truncated_gcv_1] =...
    generalized_cross_validation(J_xyz_truncated_1,dB,lambda_vector,Un1,Sn1,Vn1);

fprintf('Mdeit-3 truncated (id2)\n')
[lambda_opt_mdeit_3_truncated_2,sigma_mdeit_3_truncated_gcv_2] =...
    generalized_cross_validation(J_xyz_truncated_2,dB,lambda_vector,Un2,Sn2,Vn2);

fprintf('Mdeit-3 truncated (id3)\n')
[lambda_opt_mdeit_3_truncated_3,sigma_mdeit_3_truncated_gcv_3] =...
    generalized_cross_validation(J_xyz_truncated_3,dB,lambda_vector,Un3,Sn3,Vn3);

fprintf('Mdeit-3\n')
[lambda_opt_mdeit_3,sigma_mdeit_3_gcv] =...
    generalized_cross_validation(J,dB,lambda_vector,U,S,V);


%% Reconstruct with only TSVD (without Tikhonov and wihout NB

fprintf('Doing cross validation\n')

function sigma = reconstruct_tsvd(U,S,V,data)
sigma = V*diag(1./(diag(S)))*U'*data;
end




fprintf('Mdeit-z\n')
sigma_mdeit_z_tsvd = reconstruct_tsvd(Uz,Sz,Vz,dBz);

fprintf('Mdeit-1 truncated (id2z)\n')
sigma_mdeit_z_truncated_tsvd_2 =  reconstruct_tsvd(Un2z,Sn2z,Vn2z,dBz);

fprintf('Mdeit-1 truncated (id3z)\n')
sigma_mdeit_z_truncated_tsvd_3 = reconstruct_tsvd(Un3z,Sn3z,Vn3z,dBz);

fprintf('Mdeit-3 truncated (id1)\n')
sigma_mdeit_3_truncated_tsvd_1 = reconstruct_tsvd(Un1,Sn1,Vn1,dB);

fprintf('Mdeit-3 truncated (id2)\n')
sigma_mdeit_3_truncated_tsvd_2 = reconstruct_tsvd(Un2,Sn2,Vn2,dB);

fprintf('Mdeit-3 truncated (id3)\n')
sigma_mdeit_3_truncated_tsvd_3 =  reconstruct_tsvd(Un3,Sn3,Vn3,dB);

fprintf('Mdeit-3\n')
sigma_mdeit_3_tsvd= reconstruct_tsvd(U,S,V,dB);

%% Noise correction

function dBx = noisy_data_generator_mdeit_x(imgh,imgi,SNRdb,num_noise_repetitions)
% Forward solve
[datah,~] = fwd_solve_mdeit(imgh);
Bxh = datah.Bx(:);

% Forward solve
[datai,~] = fwd_solve_mdeit(imgi);
Bxi = datai.Bx(:);

dBx = zeros(size(Bxi,1),num_noise_repetitions);

for t = 1:num_noise_repetitions
    dBx(:,t) = add_measurement_noise_difference(Bxi,Bxh,SNRdb);
end

end

function dBy = noisy_data_generator_mdeit_y(imgh,imgi,SNRdb,num_noise_repetitions)
% Forward solve
[datah,~] = fwd_solve_mdeit(imgh);
Byh = datah.By(:);

% Forward solve
[datai,~] = fwd_solve_mdeit(imgi);
Byi = datai.By(:);

dBy = zeros(size(Byi,1),num_noise_repetitions);

for t = 1:num_noise_repetitions
    dBy(:,t) = add_measurement_noise_difference(Byi,Byh,SNRdb);
end

end

function dBz = noisy_data_generator_mdeit_z(imgh,imgi,SNRdb,num_noise_repetitions)
% Forward solve
[datah,~] = fwd_solve_mdeit(imgh);
Bzh = datah.Bz(:);

% Forward solve
[datai,~] = fwd_solve_mdeit(imgi);
Bzi = datai.Bz(:);

dBz = zeros(size(Bzi,1),num_noise_repetitions);

for t = 1:num_noise_repetitions
    dBz(:,t) = add_measurement_noise_difference(Bzi,Bzh,SNRdb);
end


end

function dB = noisy_data_generator_mdeit_3(imgh,imgi,SNRdb,num_noise_repetitions)
% Forward solve
[datah,uh] = fwd_solve_mdeit(imgh);
Bxh = datah.Bx(:);
Byh = datah.By(:);
Bzh = datah.Bz(:);
Bh = [datah.Bx(:);datah.By(:);datah.Bz(:)];

% Forward solve
[datai,ui] = fwd_solve_mdeit(imgi);
Bxi = datai.Bx(:);
Byi = datai.By(:);
Bzi = datai.Bz(:);
Bi = [datai.Bx(:);datai.By(:);datai.Bz(:)];

dB = zeros(size(Bi,1),num_noise_repetitions);

for t = 1:num_noise_repetitions
    dB(:,t) = add_measurement_noise_difference(Bi,Bh,SNRdb);
end

end

func_mdeit_x = @(imgh,imgi,num_noise_repetitions) noisy_data_generator_mdeit_x(imgh,imgi,SNRdb,num_noise_repetitions);
func_mdeit_y = @(imgh,imgi,num_noise_repetitions) noisy_data_generator_mdeit_y(imgh,imgi,SNRdb,num_noise_repetitions);
func_mdeit_z = @(imgh,imgi,num_noise_repetitions) noisy_data_generator_mdeit_z(imgh,imgi,SNRdb,num_noise_repetitions);
func_mdeit_3 = @(imgh,imgi,num_noise_repetitions) noisy_data_generator_mdeit_3(imgh,imgi,SNRdb,num_noise_repetitions);

fprintf('Doing noise correction\n');

% fprintf('Mdeit-x\n')
% std_sigma_mdeit_x = noise_correction(imgh,imgi,Jx,lambda_opt_mdeit_x,func_mdeit_x,num_noise_repetitions,Ux,Sx,Vx);
% fprintf('Mdeit-y\n')
% std_sigma_mdeit_y = noise_correction(imgh,imgi,Jy,lambda_opt_mdeit_y,func_mdeit_y,num_noise_repetitions,Uy,Sy,Vy);

fprintf('Mdeit-z\n')
std_sigma_mdeit_z = noise_correction(imgh,imgi,Jz,lambda_opt_mdeit_z,func_mdeit_z,num_noise_repetitions,Uz,Sz,Vz);
fprintf('Mdeit-3 truncated (id1)\n')
std_sigma_mdeit_3_truncated_1 = noise_correction(imgh,imgi,J_xyz_truncated_1,lambda_opt_mdeit_3_truncated_1,func_mdeit_3,num_noise_repetitions,Un1,Sn1,Vn1);
fprintf('Mdeit-3 truncated (id2)\n')
std_sigma_mdeit_3_truncated_2 = noise_correction(imgh,imgi,J_xyz_truncated_2,lambda_opt_mdeit_3_truncated_2,func_mdeit_3,num_noise_repetitions,Un2,Sn2,Vn2);
fprintf('Mdeit-1 truncated (id2z)\n')
std_sigma_mdeit_z_truncated_2 = noise_correction(imgh,imgi,J_z_truncated_2,lambda_opt_mdeit_z_truncated_2,func_mdeit_z,num_noise_repetitions,Un2z,Sn2z,Vn2z);
fprintf('Mdeit-3 truncated (id3)\n')
std_sigma_mdeit_3_truncated_3 = noise_correction(imgh,imgi,J_xyz_truncated_3,lambda_opt_mdeit_3_truncated_3,func_mdeit_3,num_noise_repetitions,Un3,Sn3,Vn3);
fprintf('Mdeit-1 truncated (id3z)\n')
std_sigma_mdeit_z_truncated_3 = noise_correction(imgh,imgi,J_z_truncated_3,lambda_opt_mdeit_z_truncated_3,func_mdeit_z,num_noise_repetitions,Un3z,Sn3z,Vn3z);


fprintf('Mdeit-3\n')
std_sigma_mdeit_3 = noise_correction(imgh,imgi,J,lambda_opt_mdeit_3,func_mdeit_3,num_noise_repetitions,U,S,V);

fprintf('Done\n');

%% Solve also with Gauss-Newton

% fprintf('Solving with gn\n')
% function sigma = gauss_newton_solve(J,data,lambda)
% %Form normal equation    
% lhs = J.'*J + lambda*eye(size(J,2));
% rhs = -J.'*data;
% 
% d = sqrt(diag(lhs));        % vector of diagonal entries
% 
% Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
% Nfun = @(x) x ./ d;              % right preconditioner
% 
% tol = 1e-5;
% maxit = 1000;
% sigma = pcg(lhs,rhs,tol,maxit,Mfun,Nfun);
% 
% end

% sigma_gn_mdeit_z = gauss_newton_solve(Jz,-dBz,lambda_opt_mdeit_1);
% sigma_gn_mdeit_3_truncated = gauss_newton_solve(J_truncated,-dB,lambda_opt_mdeit_3_truncated);
% sigma_gn_mdeit_3 = gauss_newton_solve(J,-dB,lambda_opt_mdeit_3);
%% Figure

figure('Position',[100,100,1300,600])

num_of_subplots_per_row = 8;

z_level = model_parameters.height/2;
npoints= 128;

% Function to plot a slice of the 3D image
function h = plot_image_slice(z_level,npoints,img)
img.calc_colours.transparency_thresh = -1;
img.calc_colours.npoints = npoints;
h = show_3d_slices(img,z_level,[],[]);
view(2);
axis square;
box on;
end

id_column = 1;
% First row
subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_gcv;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-z','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_truncated_gcv_2;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-z truncated (id2)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_truncated_gcv_3;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-z truncated (id3)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_1;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 truncated (id1)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_2;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 truncated (id2)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_3;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 truncated (id3)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_gcv;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 not truncated','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = imgi;
show_fem_transparent_edges(imgtemp);
view(2)
title('Ground Truth','Interpreter','latex')
id_column = id_column+1;


% Second row
subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_gcv./std_sigma_mdeit_z;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-z','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_truncated_gcv_2./std_sigma_mdeit_z_truncated_2;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-z truncated (id2)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_truncated_gcv_3./std_sigma_mdeit_z_truncated_3;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-z truncated (id3)','Interpreter','latex')
id_column = id_column+1;


subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_1./std_sigma_mdeit_3_truncated_1;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 truncated (id1)','Interpreter','latex')
id_column = id_column+1;


subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_2./std_sigma_mdeit_3_truncated_2;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 truncated (id2)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_3./std_sigma_mdeit_3_truncated_3;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 truncated (id2)','Interpreter','latex')
id_column = id_column+1;


subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_gcv./std_sigma_mdeit_3;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 not truncated','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = imgi;
show_fem_transparent_edges(imgtemp);
view(2);
title('Ground Truth','Interpreter','latex')
id_column = id_column+1;


% Third row
subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_gcv;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-z','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_truncated_gcv_2;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-z truncated (id2)','Interpreter','latex')
id_column = id_column+1;


subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_truncated_gcv_3;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-z truncated (id3)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_1;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-3 truncated (id1)','Interpreter','latex')
id_column = id_column+1;


subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_2;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-3 truncated (id2)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_3;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-3 truncated (id3)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_gcv;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-3 not truncated','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = imgi;
plot_image_slice(z_level,npoints,imgtemp);
title('Ground Truth','Interpreter','latex')
id_column = id_column+1;

% Fourth Row
subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data =sigma_mdeit_z_gcv./std_sigma_mdeit_z;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-z','Interpreter','latex')
id_column = id_column+1;


subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data =sigma_mdeit_z_truncated_gcv_2./std_sigma_mdeit_z_truncated_2;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-z truncated (id2)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data =sigma_mdeit_z_truncated_gcv_3./std_sigma_mdeit_z_truncated_3;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-z truncated (id3)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_1./std_sigma_mdeit_3_truncated_1;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-3 truncated (id1)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_2./std_sigma_mdeit_3_truncated_2;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-3 truncated (id2)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_gcv_3./std_sigma_mdeit_3_truncated_3;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-3 truncated (id3)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_gcv./std_sigma_mdeit_3;
plot_image_slice(z_level,npoints,imgtemp);
title('Mdeit-3 not truncated','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = imgi;
plot_image_slice(z_level,npoints,imgtemp);
title('Ground Truth','Interpreter','latex')
id_column = id_column+1;



%Fifth row
subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_tsvd;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-z','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_truncated_tsvd_2;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-z truncated (id2)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_z_truncated_tsvd_3;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-z truncated (id3)','Interpreter','latex')
id_column = id_column+1;


subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_tsvd_1;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 truncated (id1)','Interpreter','latex')
id_column = id_column+1;


subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_tsvd_2;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 truncated (id2)','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_truncated_tsvd_3;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 truncated (id3)','Interpreter','latex')
id_column = id_column+1;


subplot(5,num_of_subplots_per_row,id_column)
imgtemp = mk_image_mdeit(imgh_reconstruction);
imgtemp.elem_data = sigma_mdeit_3_tsvd;
show_fem_transparent_edges(imgtemp);
view(2)
title('Mdeit-3 not truncated','Interpreter','latex')
id_column = id_column+1;

subplot(5,num_of_subplots_per_row,id_column)
imgtemp = imgi;
show_fem_transparent_edges(imgtemp);
view(2);
title('Ground Truth','Interpreter','latex')
id_column = id_column+1;





%% Functions
function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
end

function [Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A)

mu0 = img.fwd_model.mu0;

n_nodes =  size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices_local_optimized(img);

R = img.fwd_model.R;
G = img.fwd_model.G;

Gamma1 = img.Gamma1;
Gamma2 = img.Gamma2;
Gamma3 = img.Gamma3;

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

A_matrix = A(img.elem_data);

Gamma1T = Gamma1.';
Gamma2T = Gamma2.';
Gamma3T = Gamma3.';

% Solve the adjoint problem for each sensor to get lambda vectors

lambdaX = A_matrix \ (-Gamma1T);
lambdaY = A_matrix \ (-Gamma2T);
lambdaZ = A_matrix \ (-Gamma3T);

Gx_times_lambda_X = G.Gx*lambdaX;
Gy_times_lambda_X = G.Gy*lambdaX;
Gz_times_lambda_X = G.Gz*lambdaX;

Gx_times_lambda_Y = G.Gx*lambdaY;
Gy_times_lambda_Y = G.Gy*lambdaY;
Gz_times_lambda_Y = G.Gz*lambdaY;

Gx_times_lambda_Z = G.Gx*lambdaZ;
Gy_times_lambda_Z = G.Gy*lambdaZ;
Gz_times_lambda_Z = G.Gz*lambdaZ;

Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% We want to broadcast arrays into num_sensors*num_stim*num_elems, in that order, so we avoid a permute in dfd

% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [1 1 n_elem]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u.', [1 num_stim n_elem]); % [: × numStim × 1]
GyU = reshape(Gy_times_u.', [1 num_stim n_elem]);
GzU = reshape(Gz_times_u.', [1 num_stim n_elem]);

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx, [num_sensors 1 n_elem]);
Ry_ = reshape(R.Ry, [num_sensors 1 n_elem]);
Rz_ = reshape(R.Rz, [num_sensors 1 n_elem]);

% These are the derivatives with respect to sigma of the C components of the Gamma matrix,
dCxdp = ( -Rz_.*GyU + Ry_.*GzU );
dCydp = ( -Rx_.*GzU + Rz_.*GxU );
dCzdp = ( -Ry_.*GxU + Rx_.*GyU );

% The g matrix does not depend on sigma.
g = zeros(num_sensors,3,3);
for m = 1:numel(img.fwd_model.sensors)
    g(m,:,:) = [...
        img.fwd_model.sensors(m).axes.axis1;
        img.fwd_model.sensors(m).axes.axis2;
        img.fwd_model.sensors(m).axes.axis3];
end

for select_sensor_axis = 1:3

    switch select_sensor_axis
        case 1
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_X.', [num_sensors 1 n_elem]); % [: × 1 × numSensors]
            GyL = reshape(Gy_times_lambda_X.', [num_sensors 1 n_elem]);
            GzL = reshape(Gz_times_lambda_X.', [num_sensors 1 n_elem]);
        case 2
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_Y.', [num_sensors 1 n_elem]); % [: × 1 × numSensors]
            GyL = reshape(Gy_times_lambda_Y.', [num_sensors 1 n_elem]);
            GzL = reshape(Gz_times_lambda_Y.', [num_sensors 1 n_elem]);
        case 3
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_Z.', [num_sensors 1 n_elem]); % [: × 1 × numSensors]
            GyL = reshape(Gy_times_lambda_Z.', [num_sensors 1 n_elem]);
            GzL = reshape(Gz_times_lambda_Z.', [num_sensors 1 n_elem]);
    end

    % Compute all dfdx for all sensors+stim
    dfdx = elemV .* ( ...
        GxL.*GxU + ...
        GyL.*GyU + ...
        GzL.*GzU );



    % g: [num_sensors × 3 × 3]
    gx = reshape(g(:,select_sensor_axis,1), [num_sensors 1 1 ]);
    gy = reshape(g(:,select_sensor_axis,2), [num_sensors 1 1]);
    gz = reshape(g(:,select_sensor_axis,3), [num_sensors 1 1]);

    dfdp = mu_factor*(...
        gx.*dCxdp +...
        gy.*dCydp +...
        gz.*dCzdp);

    dfd = dfdx + dfdp;   % size: [numStim × numSensors × numElems]

    % Now reshape to match J(ids,:)

    % collapse first 2 dims → [numSensors*numStim × numElems]
    J = reshape(dfd, num_sensors*num_stim, n_elem);

    switch select_sensor_axis
        case 1
            Jx = J;
        case 2
            Jy = J;
        case 3
            Jz = J;
    end

end

return
end

function [Rx,Ry,Rz,fmdl] = compute_r_matrices_local_optimized(fmdl,sensor_locations)

numElements = size(fmdl.elems,1);
numSensors = size(sensor_locations,1);

% ------------------------------------------------------------
% Persistent quadrature rule (avoid reallocation every call)
% ------------------------------------------------------------
persistent ref_pts weights nq;

if isempty(ref_pts)

    coord = [    0.0000000000000000  0.0000000000000000  0.0000000000000000
        0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  1.0000000000000000  0.0000000000000000
        1.0000000000000000  0.0000000000000000  0.0000000000000000
        0.0000000000000000  0.0000000000000000  0.7500000000000000
        0.0000000000000000  0.0000000000000000  0.2500000000000000
        0.0000000000000000  0.7500000000000000  0.0000000000000000
        0.0000000000000000  0.7500000000000000  0.2500000000000000
        0.0000000000000000  0.2500000000000000  0.0000000000000000
        0.0000000000000000  0.2500000000000000  0.7500000000000000
        0.7500000000000000  0.0000000000000000  0.0000000000000000
        0.7500000000000000  0.0000000000000000  0.2500000000000000
        0.7500000000000000  0.2500000000000000  0.0000000000000000
        0.2500000000000000  0.0000000000000000  0.0000000000000000
        0.2500000000000000  0.0000000000000000  0.7500000000000000
        0.2500000000000000  0.7500000000000000  0.0000000000000000
        0.5000000000000000  0.5000000000000000  0.0000000000000000
        0.5000000000000000  0.0000000000000000  0.5000000000000000
        0.5000000000000000  0.0000000000000000  0.0000000000000000
        0.0000000000000000  0.5000000000000000  0.5000000000000000
        0.0000000000000000  0.5000000000000000  0.0000000000000000
        0.0000000000000000  0.0000000000000000  0.5000000000000000
        0.2500000000000000  0.2500000000000000  0.0000000000000000
        0.2500000000000000  0.2500000000000000  0.5000000000000000
        0.2500000000000000  0.0000000000000000  0.2500000000000000
        0.2500000000000000  0.0000000000000000  0.5000000000000000
        0.2500000000000000  0.5000000000000000  0.2500000000000000
        0.2500000000000000  0.5000000000000000  0.0000000000000000
        0.0000000000000000  0.2500000000000000  0.2500000000000000
        0.0000000000000000  0.2500000000000000  0.5000000000000000
        0.0000000000000000  0.5000000000000000  0.2500000000000000
        0.5000000000000000  0.2500000000000000  0.2500000000000000
        0.5000000000000000  0.2500000000000000  0.0000000000000000
        0.5000000000000000  0.0000000000000000  0.2500000000000000
        0.2500000000000000  0.2500000000000000  0.2500000000000000];

    weights = [  -0.0119047619047619
        -0.0119047619047619
        -0.0119047619047619
        -0.0119047619047619
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        -0.0285714285714286
        -0.0285714285714286
        -0.0285714285714286
        -0.0285714285714286
        -0.0285714285714286
        -0.0285714285714286
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.3047619047619048];
    ref_pts = coord.';      % 3 × nq
    nq = size(ref_pts,2);
end

nodes = fmdl.nodes;
elems = fmdl.elems;

% ------------------------------------------------------------
% Precompute element geometry (done once!)
% ------------------------------------------------------------

v1  = zeros(3,numElements);
Jall = zeros(3,3,numElements);
detJ = zeros(numElements,1);

for k = 1:numElements

    verts = nodes(elems(k,:),:);

    v1(:,k) = verts(1,:).';

    J = [ ...
        (verts(2,:) - verts(1,:))' , ...
        (verts(3,:) - verts(1,:))' , ...
        (verts(4,:) - verts(1,:))' ];

    Jall(:,:,k) = J;
    detJ(k)     = abs(det(J));
end

% ------------------------------------------------------------
% Allocate outputs
% ------------------------------------------------------------

Rx = zeros(numSensors,numElements);
Ry = zeros(numSensors,numElements);
Rz = zeros(numSensors,numElements);

% ------------------------------------------------------------
% Main loop over sensors (parallelizable)
% ------------------------------------------------------------
for m = 1:numSensors   % <-- change to "for" if no Parallel Toolbox

    rm = sensor_locations(m,:).';

    Rx_local = zeros(1,numElements);
    Ry_local = zeros(1,numElements);
    Rz_local = zeros(1,numElements);

    for k = 1:numElements

        % Map quadrature points to physical tetrahedron
        Xi = v1(:,k) + Jall(:,:,k) * ref_pts;    % 3 × nq

        diff = rm - Xi;                         % 3 × nq
        r2   = sum(diff.^2,1);                  % 1 × nq

        kernel = diff ./ (r2.^1.5);             % 3 × nq

        R = (detJ(k)/6) * (kernel * weights);  % 3 × 1

        Rx_local(k) = R(1);
        Ry_local(k) = R(2);
        Rz_local(k) = R(3);
    end

    Rx(m,:) = Rx_local;
    Ry(m,:) = Ry_local;
    Rz(m,:) = Rz_local;
end

% ------------------------------------------------------------
% Store
% ------------------------------------------------------------

fmdl.R.Rx = Rx;
fmdl.R.Ry = Ry;
fmdl.R.Rz = Rz;

end

function img = compute_gamma_matrices_local_optimized(img)


mu_factor = img.fwd_model.mu0/(4*pi);

num_sensors = numel(img.fwd_model.sensors);

sensor_locations = zeros(numel(img.fwd_model.sensors),3);

for i = 1: numel(img.fwd_model.sensors)
    sensor_locations(i,:) = img.fwd_model.sensors(i).position;
end

[Rx,Ry,Rz,fmdl] = compute_r_matrices_local_optimized(img.fwd_model,sensor_locations);

img.fwd_model = fmdl; % THIS IS CRUCIAL!!!!! so the R matrices are updated whenever compute_gamma_matrices_local is called, otherwise they will be the same as the initial img.fwd_model, which might be wrong!

% Convenience handles
R.Rx = Rx;
R.Ry = Ry;
R.Rz = Rz;
G = img.fwd_model.G;

% Sigma = sparse(1:length(img.elem_data), 1:length(img.elem_data), img.elem_data);
Sigma = spdiags(img.elem_data(:), 0, length(img.elem_data), length(img.elem_data));

% NEW: The matrix g_{dl}^m, gives the components of the measurement axis of
% sensor m on the canonical R^3 basis

g = zeros(num_sensors,3,3);
for m = 1:numel(img.fwd_model.sensors)
    g(m,:,:) = [...
        img.fwd_model.sensors(m).axes.axis1;
        img.fwd_model.sensors(m).axes.axis2;
        img.fwd_model.sensors(m).axes.axis3];
end

Cx = ( -R.Rz * Sigma * G.Gy +  R.Ry * Sigma * G.Gz );
Cy = ( -R.Rx * Sigma * G.Gz +  R.Rz * Sigma * G.Gx );
Cz = ( -R.Ry * Sigma * G.Gx +  R.Rx * Sigma * G.Gy );

Gamma1 = mu_factor*(g(:,1,1).*Cx + g(:,1,2).*Cy + g(:,1,3).*Cz);
Gamma2 = mu_factor*(g(:,2,1).*Cx + g(:,2,2).*Cy + g(:,2,3).*Cz);
Gamma3 = mu_factor*(g(:,3,1).*Cx + g(:,3,2).*Cy + g(:,3,3).*Cz);

img.Gamma1 = Gamma1;
img.Gamma2 = Gamma2;
img.Gamma3 = Gamma3;

end



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

