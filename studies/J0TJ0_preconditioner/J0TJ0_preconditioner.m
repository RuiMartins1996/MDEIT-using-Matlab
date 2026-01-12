
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

%% Experiment parameters

num_of_repetitions = 1;
tolerances = 10.^[-3,-7];
lambdas = [1e-3,1e-5];

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
SNRdb = [];

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


anomaly_conductivity = 1e-2/background_conductivity; % (0 according to Kai's thesis, but check notes)
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

%% For difference reconstruction

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

lambda_vector = logspace(-15,3,30);

[lambda_optimal_mdeit1,optimal_id_mdeit_1,V_mu_mdeit1,~] = ...
    generalized_cross_validation(imdl_mdeit_1,Bzi-Bzh,lambda_vector);

img = mk_image_mdeit(fmdl_reconstruction,background_conductivity);
sigma0 = img.elem_data;

x0 = img.elem_data;
max_iterations = 50;
recon_mode = 'mdeit1';
select_sensor_axis = 3;
RtR = @(x,Jx) x; % tikhonov prior

lambdatimesdAdp = @(lambda) nan; %dummy, its no longer needed 
A = @(sigma) M(img,sigma);

res = @(x) -dB;
jac = @(x) calc_jacobian_mdeit(img, x,lambdatimesdAdp,A,recon_mode,select_sensor_axis);

j0 = jac(x0);
[U,S,V] = svd(j0);
sigma = diag(S);
diagonal = [sigma;zeros(size(S,2)-numel(sigma),1)];

tol = 1e-10;

% The Moore-Penrose inverse of lhs of normal equations is 
Vt = V';
A_dagger =@(lambda) V * diag(1./(diagonal.^2+lambda)) * Vt; %moore penrose pseudo-inverse
apply_A_dagger = @(b,lambda) V * ((1./(diagonal.^2 + lambda)) .* (Vt * b));

verbose = true;

% [sigma_gn,num_iter,~,it,res1] = gn_solve(...
%     res, jac, x0, tol, 1, RtR, lambda_optimal_mdeit1,[],verbose);
% 
% [sigma_gn_pre,num_of_iterations_5,~,it5,res5] = gn_solve(...
%     res, jac, x0, tol, 1, RtR, lambda_optimal_mdeit1,[],verbose);

% figure
% img_temp = img;
% img_temp.elem_data = sigma_gn;
% subplot(1,2,1)
% show_fem(img_temp)
% img_temp = img;
% img_temp.elem_data = sigma_gn_pre;
% subplot(1,2,2)
% show_fem(img_temp)
% 
% pause(1e-10)

% COMENTS:
% For the difference reconstruction, it is quicker to use the
% preconditioner, because the preconditioner is the exact inverse of the
% lhs of the normal equations. What this proves is that it might be better
% to use TSVD for difference reconstruction, always.
%% For absolute reconstruction 

%Expose the GN method so we can explore how the preconditioner affects the solve:

img = mk_image_mdeit(fmdl_reconstruction,background_conductivity);
sigma0 = img.elem_data;

x0 = img.elem_data;
max_iterations = 50;
recon_mode = 'mdeit1';
select_sensor_axis = 3;
RtR = @(x,Jx) x;

% FOR Absolute reconstruction, the pre-conditioner choice cannot be given
% by GCV, for low noise (SNRdb = 100) !!!
% lambda = imdl_mdeit_1.hyperparameter.value;

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

j0 = jac(x0);
[U,S,V] = svd(j0);
sigma = diag(S);
diagonal = [sigma;zeros(size(S,2)-numel(sigma),1)];

% The Moore-Penrose inverse of lhs of normal equations is 
Vt = V';
A_dagger =@(lambda) V * diag(1./(diagonal.^2+lambda)) * Vt; %moore penrose pseudo-inverse
apply_A_dagger = @(b,lambda) V * ((1./(diagonal.^2 + lambda)) .* (Vt * b));

%Sanity check
% fprintf('Condition number of preconditioner times matrix:\n %.1d\n',abs(cond(A_dagger(lambda)\A)-1));

% tol = max(size(S)) * eps(max(diag(S)));
% A_moore_penrose = @(lambda) pinv(j0'*j0+lambda*speye(size(j0,2)),tol);

l2_norm_ground_truth = computeL2norm_img(imgi);

verbose = true;

apply_A_dagger = @(b,lambda) V * ((1./(diagonal.^2 + lambda)) .* (Vt * b));
apply_A_dagger_gn = @(b) V * ((1./(diagonal.^2 + lambda_optimal_mdeit1)) .* (Vt * b));

% [sigma_gn,num_of_iterations_1,it1,res1,exit_flags_gn] = gn_solve(res, jac, x0, tol, max_iterations, RtR,lambda_optimal_mdeit1,[],verbose);
% [sigma_gn_pre,num_of_iterations_2,it2,res2,exit_flags_gn_pre] = gn_solve(res, jac, x0, tol, max_iterations, RtR,lambda_optimal_mdeit1,apply_A_dagger_gn,verbose);

%%
% figure('Position',[100,100,1000,500])
% subplot(1,2,1)
% img_out = img;
% img_out.elem_data = sigma_gn;
% show_fem_corrected(img_out);
% eidors_colourbar(img_out)
% title('GN - No Pre')
% 
% subplot(1,2,2)
% img_out = img;
% img_out.elem_data = sigma_gn_pre;
% show_fem_corrected(img_out);
% eidors_colourbar(img_out)
% title('GN - Pre')

%%

for i = 1:numel(lambdas)
    lambda = lambdas(i);

    Apre = j0' * j0 + lambda * speye(size(j0,2));

    opts.type = 'ict';        % incomplete Cholesky with threshold
    opts.droptol = 1e-3;      % adjust for accuracy/speed
    R = chol(Apre, 'upper');    % R'*R ≈ Apre
    
    P = @(v) R \ (R' \ v);  % M * v ≈ A^{-1} * v

    for j = 1:numel(tolerances)
        tol = tolerances(j);

        times1 = zeros(num_of_repetitions,1);
        times2 = zeros(num_of_repetitions,1);
        times3 = zeros(num_of_repetitions,1);

        for n = 1:num_of_repetitions
            tic
            [sigma_lm_no_pre,num_of_iterations_1,it1,res1,exit_flags_lm_no_pre(j)] = lm_solve(res, jac, x0, tol, max_iterations, RtR,lambda,[],verbose);
            times1(n) = toc;


            tic
            [sigma_lm,num_of_iterations_2,it2,res2,exit_flags_lm(j)] = lm_solve(res, jac, x0, tol, max_iterations, RtR,lambda,apply_A_dagger,verbose);
            times2(n) = toc;
            
            tic
            [sigma_gn,num_of_iterations_3,it3,res3,exit_flags_gn(j)] = gn_solve(res, jac, x0, tol, max_iterations, RtR,lambda_optimal_mdeit1,[],verbose);
            times3(n) = toc;
            
            tic
            [sigma_gn_pre,num_of_iterations_4,it4,res4,exit_flags_gn_pre(j)]  = gn_solve(res, jac, x0, tol, max_iterations, RtR,lambda_optimal_mdeit1,apply_A_dagger_gn,verbose);
            times4(n) = toc;
        end

        pcg_iterations_lm_no_pre{i,j} = it1;
        pcg_iterations_lm{i,j} = it2;
        pcg_iterations_gn{i,j} = it3;
        pcg_iterations_gn_pre{i,j} = it4;
            
        residual_lm_no_pre(i,j) = norm(res1);
        residual_lm(i,j) = norm(res2);
        residual_gn(i,j) = norm(res3);
        residual_gn_pre(i,j) = norm(res4);

        num_of_iterations_lm_no_pre(i,j) = num_of_iterations_1;
        num_of_iterations_lm(i,j) = num_of_iterations_2;
        num_of_iterations_gn(i,j) = num_of_iterations_3;
        num_of_iterations_gn_pre(i,j) = num_of_iterations_4;

        time_lm_no_pre(i,j) = mean(times1);
        std_lm_no_pre(i,j) = std(times1);
        time_lm(i,j) = mean(times2);
        std_lm(i,j) = std(times2);
        time_gn(i,j) = mean(times3);
        std_gn(i,j) = std(times3);
        time_gn_pre(i,j) = mean(times4);
        std_gn_pre(i,j) = std(times4);
    end
end

%% Time w.r.t tolerances
figure('Position',[100,100,1000,500]);
for i = 1:numel(lambdas)
    subplot(1,numel(lambdas),i)

    hold on
    errorbar(tolerances,time_lm_no_pre(i,:),std_lm_no_pre(i,:),'b.-');
    errorbar(tolerances,time_lm(i,:),std_lm(i,:),'r.-');
    errorbar(tolerances,time_gn(i,:),std_gn(i,:),'g.-');
    errorbar(tolerances,time_gn_pre(i,:),std_gn_pre(i,:),'k.-');
    hold off
    
    set(gca,'YScale','log')
    legend('LM - No Pre','LM - Pre','GN - No Pre','GN - Pre')
    set(gca,'XScale','log');
    grid on;grid minor;
    box on;
    xlabel('tol','Interpreter','latex');
    ylabel('time(s)','Interpreter','latex')
    title_str = strcat('$\lambda = ',num2str(lambdas(i)),'$');
    title(title_str,'Interpreter','latex')
end

%%

figure('Position',[100,100,1000,500]);
for i = 1:numel(lambdas)
        subplot(1,numel(lambdas),i)

    hold on
    plot(residual_lm_no_pre(i,:),time_lm_no_pre(i,:),'b.')
    plot(residual_lm(i,:),time_lm(i,:),'r * ')

    hold off

    legend({'LM - No Pre','LM - Pre'})
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    grid on;grid minor;
    box on;
    ylabel('time(s)','Interpreter','latex')
    xlabel('error','Interpreter','latex')
    title_str = strcat('$\lambda = ',num2str(lambdas(i)),'$');
    title(title_str,'Interpreter','latex')
end

%% Plots of number of PCG iterations per LM iteration
figure('Position',[100,100,1000,500]);

names = {'LM - No Pre','LM - Pre','GN - No Pre','GN - Pre'};
markers = {'p','s','o','p','h'};

all_plots = {pcg_iterations_gn,pcg_iterations_gn_pre};

id = 0;
for i = 1:numel(lambdas)
    id = id+1;
    subplot(numel(lambdas),2,id)
    ct = 1;
    legendStr = {i,j};

    % this_plot = pcg_iterations_lm_no_pre;
    this_plot = all_plots{i};
    hold on

    for j = numel(tolerances):-1:1
        marker = markers{exit_flags_lm_no_pre(j)};
        color = colors(j,:);
        plot(this_plot{i,j},'Color',color,'LineWidth',1,'Marker',marker)
        legendStr{ct} = sprintf('tol = %.1d',tolerances(j));
        ct = ct+1;
    end

    hold off
    legend(legendStr)

    title_str = strcat('LM - No Pre $\lambda = ',num2str(lambdas(i)),'$');
    title(title_str,'Interpreter','latex')    

    box on;
    grid on

    ylabel('number of PCG iterations','Interpreter','latex')
    xlabel('number of LM iterations','Interpreter','latex')
    
    id = id+1;
    subplot(numel(lambdas),2,id)
    ct = 1;
    legendStr = {};

    this_plot = pcg_iterations_lm;
    hold on

    for j = numel(tolerances):-1:1
        marker = markers{exit_flags_lm(j)};
        color = colors(j,:);
        plot(this_plot{i,j},'Color',color,'LineWidth',1,'Marker',marker)
        legendStr{ct} = sprintf('tol = %.1d',tolerances(j));
        ct = ct+1;
    end

    hold off
    legend(legendStr)

    title_str = strcat('LM - Pre $\lambda = ',num2str(lambdas(i)),'$');
    title(title_str,'Interpreter','latex')  

    box on;
    grid on

    ylabel('number of PCG iterations','Interpreter','latex')
    xlabel('number of LM iterations','Interpreter','latex')
end
%% Statistics
i = 2;
for j = 1:numel(tolerances)

    a1 = pcg_iterations_lm_no_pre{i,j};
    a2 = pcg_iterations_lm{i,j};


    fprintf('_____________________________\n');
    fprintf('tol = %.2g\n',tolerances(j));
    fprintf('Total number of PCG iterations LM: %i\n',sum(a1));
    fprintf('Total number of PCG iterations LM - Pre: %i\n',sum(a2));
end

for i = 1:numel(tolerances)
    
    a1 = residual_lm_no_pre(i);
    a2 = residual_lm(i);

    fprintf('_____________________________\n');
    fprintf('tol = %.2g\n',tolerances(i));
    fprintf('Residual LM: %i\n',a1);
    fprintf('Residual LM - Pre: %i\n',a2);
end

for i = 1:numel(tolerances)
    
    a1 = time_lm_no_pre(i);
    a2 = time_lm(i);

    b1 = std_lm_no_pre(i);
    b2 = std_lm(i);

    fprintf('_____________________________\n');
    fprintf('tol = %.2g\n',tolerances(i));
    fprintf('Time LM: %.2f +/- %.2f\n',a1,b1);
    fprintf('Time LM - Pre: %.2f +/- %.2f\n',a2,b2);
end
%%


figure('Position',[100,100,1000,500])
subplot(1,2,1)
img_out = img;
img_out.elem_data = sigma_lm_no_pre;
show_fem_corrected(img_out);
eidors_colourbar(img_out)
title('LM - No Pre')

subplot(1,2,2)
img_out = img;
img_out.elem_data = sigma_lm;
show_fem_corrected(img_out);
eidors_colourbar(img_out)
title('LM - Pre')

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

%% FUNCTIONs
function hh = show_fem_corrected(img)

hh = show_fem(img);                % draw the model (hh may be a handle or array)
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);
axis off
end

function L2norm = computeL2norm_img(img)
% computeL2norm_img - computes the L2 norm of a piecewise constant function
% defined on an EIDORS fwd_model
%
% Input:
%   img - EIDORS img structure (elem_data must be defined)
%
% Output:
%   L2norm - scalar, the L2 norm of the piecewise constant function

if ~isfield(img, 'elem_data') || isempty(img.elem_data)
    error('img.elem_data is empty or missing');
end

f = img.elem_data;           % piecewise constant values
mdl = img.fwd_model;         % underlying mesh

num_elem = size(mdl.elems,1);
dim = size(mdl.nodes,2);

L2sq = 0;

for k = 1:num_elem
    verts = mdl.nodes(mdl.elems(k,:), :);
    
    % Compute element volume/area
    switch dim
        case 1
            vol = abs(verts(end) - verts(1));
        case 2
            vol = polyarea(verts(:,1), verts(:,2));
        case 3
            vol = abs(dot(verts(2,:) - verts(1,:), ...
                          cross(verts(3,:) - verts(1,:), verts(4,:) - verts(1,:))))/6;
        otherwise
            error('Unsupported dimension');
    end
    
    L2sq = L2sq + f(k)^2 * vol;
end

L2norm = sqrt(L2sq);
end
