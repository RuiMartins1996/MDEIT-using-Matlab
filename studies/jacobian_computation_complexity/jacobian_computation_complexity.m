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

file_name = strcat(script_folder,'/data/data.mat');


%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);
%% Assign the parameters for several models (should create utility functions for this)

num_of_repetitions_eit = 15;
num_of_repetitions_mdeit = 5;

maxsz_reconstruction = 0.03;

background_conductivity = 3.28e-1/sigma0;  %page 163 mentions a saline solution (NaCl + water) at 0.2% mass concentration, but can't find data for that conductivity, check notes

%% Define forward model (2D real tank experiment)

model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

% A material defines the geometry of that material inside the domain, not
% its properties, like conductivity.
model_parameters.material = struct();

% Stimulation pattern was not the default. For now, edit it manually
current_amplitude = 2.4e-3/I0;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used
options = {};

% min_maxsz = model_parameters.radius;
% max_maxsz = model_parameters.radius/10;
% 
% min_maxsz = model_parameters.radius/10;
% max_maxsz = model_parameters.radius/15;

% min_maxsz = model_parameters.radius/15;
% max_maxsz = model_parameters.radius/20;

% min_maxsz = model_parameters.radius/20;
% max_maxsz = model_parameters.radius/25;

min_maxsz = model_parameters.radius/25;
max_maxsz = model_parameters.radius/30;

n_steps = 5;

%% Sweep model_parameters over multiple cylindrical radius
model_parameters_array = ...
    sweep_model_parameters({'maxsz'},min_maxsz,max_maxsz,n_steps,model_parameters,'log');

%% First of all, build or load every forward model

% If this does not terminate, then we can do something different

fmdls_all = cell(1,numel(model_parameters_array));

for n = 1:numel(model_parameters_array)
    model_parameters = model_parameters_array(n);
    [model_parameters,fmdls] = ...
        mk_mdeit_model(model_parameters,model_folder,options);
    fmdls_all{n} = fmdls{1};
end

%% 

times_eit = zeros(numel(model_parameters_array),1);
times_mdeit = zeros(numel(model_parameters_array),1);
std_eit = zeros(numel(model_parameters_array),1);
std_mdeit = zeros(numel(model_parameters_array),1);

time_forward_solve_eit = zeros(numel(model_parameters_array),1);
std_forward_solve_eit = zeros(numel(model_parameters_array),1);

time_forward_solve_mdeit = zeros(numel(model_parameters_array),1);
std_forward_solve_mdeit = zeros(numel(model_parameters_array),1);

n_elem_vector = zeros(numel(model_parameters_array),1);

if isfile(file_name)
    loaded = load(file_name, 'T');
    T = loaded.T;  % retrieve the table
else
    % Initialize empty table if it doesn't exist yet
    T = table( ...
        [],[],[],[],[],[],[],[],[],[],[],[],...
        'VariableNames', { ...
        'num_sensors','num_electrodes_per_ring','num_rings', ...
        'times_mdeit','times_eit', ...
        'std_mdeit','std_eit', ...
        'time_forward_solve_eit','std_forward_solve_eit', ...
        'time_forward_solve_mdeit','std_forward_solve_mdeit', ...
        'n_elems'} ...
        );
end

for n = 1:numel(model_parameters_array)

    model_parameters = model_parameters_array(n);

    [model_parameters,fmdls] = ...
        mk_mdeit_model(model_parameters,model_folder,options);
    
    fmdl = fmdls{1};
    stimulation = mk_stim_patterns(numel(fmdl.electrode),1,inj,meas,options,current_amplitude);
    fmdl.stimulation = stimulation;

    % logical index of matching rows
    existing_idx = (T.num_sensors == model_parameters.numOfSensors) & ...
        (T.num_electrodes_per_ring == model_parameters.numOfElectrodesPerRing) & ...
        (T.num_rings == model_parameters.numOfRings) & ...
        (T.n_elems == size(fmdl.elems,1));

    if any(existing_idx)
        fprintf('Entry for [%d, %d, %d,%d] already exists.\n', ...
            model_parameters.numOfSensors, model_parameters.numOfElectrodesPerRing, model_parameters.numOfRings,size(fmdl.elems,1));
        skip_iteration = true;   % flag to skip computation
    else
        skip_iteration = false;
    end


    if not(skip_iteration)
        n_elem_vector(n) = size(fmdl.elems,1);

        imgh = mk_image_mdeit(fmdl,background_conductivity);

        % FORWARD SOLVE: Compute mean execution time and standard deviation
        [time_forward_solve_mdeit(n),std_forward_solve_mdeit(n)] = ...
            compute_execution_time(@fwd_solve_mdeit, num_of_repetitions_mdeit, imgh);

        [time_forward_solve_eit(n),std_forward_solve_eit(n)] = ...
            compute_execution_time(@fwd_solve, num_of_repetitions_eit, imgh);

        % JACOBIAN SOLVE: Compute mean execution time and standard deviation
        [times_eit(n),std_eit(n)] = compute_execution_time(@calc_jacobian, num_of_repetitions_eit, imgh);

        lambdatimesdAdp = @(lambda) computeLambdaTimesDaDp(imgh,lambda);
        A = @(sigma) M(imgh,sigma);

        [times_mdeit(n),std_mdeit(n)] = ...
            compute_execution_time(@ calc_jacobian_mdeit, num_of_repetitions_mdeit, imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',3);

        new_row = table( ...
            model_parameters.numOfSensors,model_parameters.numOfElectrodesPerRing,model_parameters.numOfRings,...
            times_mdeit(n), times_eit(n), ...
            std_mdeit(n), std_eit(n), ...
            time_forward_solve_eit(n), std_forward_solve_eit(n), ...
            time_forward_solve_mdeit(n), std_forward_solve_mdeit(n), ...
            n_elem_vector(n),...
            'VariableNames', { ...
            'num_sensors','num_electrodes_per_ring','num_rings', ...
            'times_mdeit','times_eit', ...
            'std_mdeit','std_eit', ...
            'time_forward_solve_eit','std_forward_solve_eit', ...
            'time_forward_solve_mdeit','std_forward_solve_mdeit', ...
            'n_elems'} ...
            );

        T = [T; new_row];
        
        % Save at each iteration so no progress is lost
        save(file_name, "T");
    end
end

save(file_name, "T");

fprintf('Done!\n')

%% FUNCTION: compute_execution_time
function [mean_time,std_dev] = compute_execution_time(func, repetitions, varargin)
    
    times = zeros(1,repetitions);
    
    tic
    for t = 1:repetitions
        tic
        func(varargin{:});
        times(t) = toc;
    end

    mean_time = sum(times) / repetitions;
    std_dev = std(times);
end




%% PLOTS
% 

% figure('Position',[100,100,1000,500]);
% cla;
% 
% subplot(1,3,1)
% hold on
% errorbar(n_nodes_vector,time_forward_solve_eit,std_forward_solve_eit,'o-','MarkerSize',5,'Color',colors(2,:))
% errorbar(n_nodes_vector,time_forward_solve_mdeit,std_forward_solve_mdeit,'d-','MarkerSize',5,'Color',colors(1,:))
% hold off
% 
% p_f_mdeit = polyfit(...
%     log10(n_nodes_vector),...
%     log10(time_forward_solve_mdeit),...
%     1);
% 
% p_f_eit = polyfit(...
%     log10(n_nodes_vector),...
%     log10(time_forward_solve_eit),...
%     1);
% 
% hold on
% x = linspace(min(n_nodes_vector),max(n_nodes_vector));
% plot(x,10^p_f_eit(2)*x.^p_f_eit(1),'LineStyle','--','Color',colors(2,:))
% plot(x,10^p_f_mdeit(2)*x.^p_f_mdeit(1),'LineStyle','--','Color',colors(1,:))
% hold off
% 
% 
% msg1 = strcat('EIT $(m = ',num2str(p_f_eit(1)),'$)'); 
% msg2 = strcat('MDEIT $(m = ',num2str(p_f_mdeit(1)),'$)'); 
% 
% legend({msg1,msg2},'Interpreter','latex','Location','northwest')
% 
% box on;
% grid on;grid minor;
% 
% set(gca,'YScale','log');
% set(gca,'XScale','log');
% 
% xlabel('N','Interpreter','latex')
% ylabel('t(s)','Interpreter','latex')
% 
% subplot(1,3,2)
% hold on
% errorbar(nnz_A_vector,time_forward_solve_eit,std_forward_solve_eit,'o-','MarkerSize',5,'Color',colors(2,:))
% errorbar(nnz_A_vector,time_forward_solve_mdeit,std_forward_solve_mdeit,'d-','MarkerSize',5,'Color',colors(1,:))
% hold off
% 
% p_A_mdeit = polyfit(...
%     log10(nnz_A_vector),...
%     log10(time_forward_solve_mdeit),...
%     1);
% 
% p_A_eit = polyfit(...
%     log10(nnz_A_vector),...
%     log10(time_forward_solve_eit),...
%     1);
% 
% hold on
% x = linspace(min(nnz_A_vector),max(nnz_A_vector));
% plot(x,10^p_A_eit(2)*x.^p_A_eit(1),'LineStyle','--','Color',colors(2,:))
% plot(x,10^p_A_mdeit(2)*x.^p_A_mdeit(1),'LineStyle','--','Color',colors(1,:))
% hold off
% 
% 
% msg1 = strcat('EIT $(m = ',num2str(p_A_eit(1)),'$)'); 
% msg2 = strcat('MDEIT $(m = ',num2str(p_A_mdeit(1)),'$)'); 
% 
% legend({msg1,msg2},'Interpreter','latex','Location','northwest')
% 
% box on;
% grid on;grid minor;
% 
% set(gca,'YScale','log');
% set(gca,'XScale','log');
% 
% xlabel('nnz(A)','Interpreter','latex')
% ylabel('t(s)','Interpreter','latex')
% 
% subplot(1,3,3)
% 
% hold on
% plot(n_nodes_vector,nnz_A_vector,'o','Color',colors(3,:));
% 
% box on;
% grid on;grid minor;
% 
% ylabel('nnz(A)','Interpreter','latex')
% xlabel('t(s)','Interpreter','latex')
% 
% p = polyfit(n_nodes_vector,nnz_A_vector,1);
% x = linspace(min(n_nodes_vector),max(n_nodes_vector));
% plot(x,p(1).*x+p(2),'Color',colors(3,:),'LineStyle','--')
% 
% msg = strcat('$(m = ',num2str(p(1)),'$)'); 
% legend({msg},'Interpreter','latex','Location','northwest')
% title('6 is the average number of neighbour nodes + 1 for diagonal')

%% PLOTS
figure;
cla;

% fit_data_points = length(n_elem_vector):-1:length(n_elem_vector)-7;

fit_data_points = 1:7;

hold on
errorbar(n_elem_vector,times_eit,std_eit,'o-','MarkerSize',5,'Color',colors(2,:))
errorbar(n_elem_vector,times_mdeit,std_mdeit,'d-','MarkerSize',5,'Color',colors(1,:))
hold off

p_mdeit = polyfit(...
    log10(n_elem_vector(fit_data_points)),...
    log10(times_mdeit(fit_data_points)),...
    1);

p_eit = polyfit(...
    log10(n_elem_vector(fit_data_points)),...
    log10(times_eit(fit_data_points)),...
    1);

% We're checking a fit of the type t = 10^(b)*K^m

hold on
x = linspace(min(n_elem_vector(fit_data_points)),max(n_elem_vector(fit_data_points)));
plot(x,10^p_eit(2)*x.^p_eit(1),'LineStyle','--','Color',colors(2,:))
plot(x,10^p_mdeit(2)*x.^p_mdeit(1),'LineStyle','--','Color',colors(1,:))
hold off

msg1 = strcat('EIT $(m = ',num2str(p_eit(1)),'$)'); 
msg2 = strcat('MDEIT $(m = ',num2str(p_mdeit(1)),'$)'); 

legend({msg1,msg2},'Interpreter','latex','Location','northwest')

box on;
grid on;grid minor;

set(gca,'YScale','log');
set(gca,'XScale','log');

min_x = min(min(n_elem_vector),10^floor(log10(min(n_elem_vector))));
max_x = max(max(n_elem_vector),10^ceil(log10(max(n_elem_vector))));
min_y = 0.5*min([times_mdeit(:);times_eit(:)]);
max_y = 1.1*max([times_mdeit(:);times_eit(:)]);

xlabel('$K$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')

axis([min_x,max_x,min_y,max_y])

title('Jacobian computation time log-scale')
%%
function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
end


