clc; clear all; close all;

% In this script, we see how the computation time of the Jacobian depends
% of the number of magnetometers and the number of electrodes

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


file_name = strcat(script_folder,'/data/data_m_l.mat');



%% Model parameters 
z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);

model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

model_parameters.maxsz =  model_parameters.radius/5;

% A material defines the geometry of that material inside the domain, not
% its properties, like conductivity.
model_parameters.material = struct();

% Stimulation pattern was not the default. For now, edit it manually
current_amplitude = 2.4e-3/I0;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used
options = {};

%% Assign the parameters for several models (should create utility functions for this)
num_of_repetitions_eit = 15;

num_of_repetitions_mdeit = 5;

background_conductivity = 3.28e-1/sigma0;  %page 163 mentions a saline solution (NaCl + water) at 0.2% mass concentration, but can't find data for that conductivity, check notes

%% Sweep model_parameters over multiple cylindrical radius

% num_of_electrodes_per_ring_array = 4:16; %because of skip 2 pattern, we have to start at 4
% num_of_rings_array = 1:6;

num_of_electrodes_per_ring_array = [4:16]; 
num_of_rings_array = 1:8;

% num_of_electrodes_per_ring_array = [4]; 
% num_of_rings_array = [1,2,3,4];

[E, R] = ndgrid(num_of_electrodes_per_ring_array, num_of_rings_array);
combinations = [E(:)'; R(:)'];
num_of_sensors = combinations(1,:).*combinations(2,:);
combinations = [combinations;num_of_sensors];

% for some reason, the following combinations make netgen stuck, so avoid them 
% [16,8,128]


remove ={...
    [16,8,128]};

% remove = {};

for n = 1:numel(remove)
    id = find(all(ismember(combinations,remove{n}'),1));
    combinations(:,id) = [];
end

model_parameters_array = ...
    sweep_model_parameters({'numOfElectrodesPerRing','numOfRings','numOfSensors'}, combinations, [], [],model_parameters);


%% First of all, build or load every forward model

% If this does not terminate, then we can do something different

fmdls_all = cell(1,numel(model_parameters_array));

for n = 1:numel(model_parameters_array)
    model_parameters = model_parameters_array(n);
    fprintf('Working on model : [%i,%i,%i]\n',model_parameters.numOfElectrodesPerRing,model_parameters.numOfRings,model_parameters.numOfSensors);
    [model_parameters,fmdls] = ...
        mk_mdeit_model(model_parameters,model_folder,options);
    fmdls_all{n} = fmdls{1};
end

%% A small plot

figure('Position',[100,100,1000,300])
for n = 1:5
    subplot(1,5,n)
    mdl = fmdls_all{8*n};
    show_fem_transparent_edges(mdl);
    plot_sensors(mdl)
    box on;
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

fprintf('Computing Jacobian times\n');

for n = 1:numel(model_parameters_array)
    fprintf('Iteration %i of %i\n',n,numel(model_parameters_array));
    model_parameters = model_parameters_array(n);

    fmdl = fmdls_all{n};
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
    [time_forward_solve_mdeit(n), std_forward_solve_mdeit(n)] = ...
        compute_execution_time(@fwd_solve_mdeit, num_of_repetitions_mdeit, imgh);

    [time_forward_solve_eit(n), std_forward_solve_eit(n)] = ...
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

electrode_count = combinations(1,:).*combinations(2,:);
sensor_count = combinations(3,:);

colors = [228,26,28;... %Colors for figure representation
    55,126,184;...
    77,175,74;...
    152,78,163;...
    255,127,0;...
    255,255,51;...
    166,86,40;...
    247,129,191;
    202,178,214;
    106,61,154]/255;

%% PLOTS

% SHOULD PACK THIS INTO A FUNCTION!!!!!!!!

figure;
cla;

% The fit for EIT time seems for to be log(t) = m*log(L^2) + b - > t =
% exp(mlog(L^2)+b) = exp(mlog(L^2)exp(b) = exp(log((L^2)^m)exp(b) =
% (L^2)^m*exp(b)
p_eit_alternative = polyfit(log(electrode_count.^2),log(times_eit),1);

x_alt = linspace(min(electrode_count.*electrode_count),max(electrode_count.*electrode_count));
y_plot_alt = x_alt.^(p_eit_alternative(1))*exp(p_eit_alternative(2));

% EIT Jacobian computation time plot
subplot(2,2,1)
hold on
errorbar(electrode_count.^2,times_eit,std_eit,'o','MarkerSize',5,'Color',colors(3,:))

p_eit = polyfit(electrode_count.*electrode_count,times_eit,1);
x = linspace(min(electrode_count.*electrode_count),max(electrode_count.*electrode_count));

y_plot = polyval(p_eit,x);
y_fit = polyval(p_eit,electrode_count.*electrode_count);
y = times_eit;

SS_res = sum((y(:) - y_fit(:)).^2);
SS_tot = sum((y(:) - mean(y(:))).^2);

R2 = 1 - (SS_res / SS_tot);

plot(x,y_plot,'.','Color',colors(4,:),'LineWidth',2);
plot(x_alt,y_plot_alt,'.','Color',colors(6,:),'LineWidth',2);
hold off
msg = strcat('$ R^2 = ',num2str(R2),'$');
legend('EIT',msg,'interpreter','latex','location','southeast');

xlabel('$L^2$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;

title('EIT Jacobian computation','Interpreter','latex')

set(gca,'YScale','log');
set(gca,'XScale','log');

% MDEIT Jacobian computation time plot
subplot(2,2,2)
hold on
errorbar(electrode_count.*sensor_count,times_mdeit,std_mdeit,'d','MarkerSize',5,'Color',colors(1,:))

p_mdeit = polyfit(electrode_count.*sensor_count,times_mdeit,1);
x = linspace(min(electrode_count.*sensor_count),max(electrode_count.*sensor_count));

y_plot = polyval(p_mdeit,x);
y_fit = polyval(p_mdeit,electrode_count.*sensor_count);
y = times_mdeit;

SS_res = sum((y(:) - y_fit(:)).^2);
SS_tot = sum((y(:) - mean(y(:))).^2);

R2 = 1 - (SS_res / SS_tot);

plot(x,y_plot,'.','Color',colors(2,:),'LineWidth',2);
hold off
msg = strcat('$ R^2 = ',num2str(R2),'$');
legend('$1$-axis MDEIT',msg,'interpreter','latex','location','southeast');

set(gca,'YScale','log');
set(gca,'XScale','log');

xlabel('$M \times L$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;

title('$1$-axis MDEIT Jacobian computation','Interpreter','latex')



% EIT forward computation time plot
subplot(2,2,3)
hold on
errorbar(electrode_count,time_forward_solve_eit,std_forward_solve_eit,'o','MarkerSize',5,'Color',colors(3,:))

p_eit = polyfit(electrode_count,time_forward_solve_eit,1);
x = linspace(min(electrode_count),max(electrode_count));

y_plot = polyval(p_eit,x);
y_fit = polyval(p_eit,electrode_count);
y = time_forward_solve_eit;

SS_res = sum((y(:) - y_fit(:)).^2);
SS_tot = sum((y(:) - mean(y(:))).^2);

R2 = 1 - (SS_res / SS_tot);

plot(x,y_plot,'.','Color',colors(4,:),'LineWidth',2);
hold off
msg = strcat('$ R^2 = ',num2str(R2),'$');
legend('EIT',msg,'interpreter','latex','location','southeast');

xlabel('$L$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;

title('EIT forward solve','Interpreter','latex')

set(gca,'YScale','log');
set(gca,'XScale','log');

% MDEIT forward computation time plot
subplot(2,2,4)
hold on
% errorbar(sensor_count,time_forward_solve_mdeit,std_forward_solve_mdeit,'o','MarkerSize',5,'Color',colors(1,:))

p_mdeit = polyfit(sensor_count,time_forward_solve_mdeit,1);
x = linspace(min(sensor_count),max(sensor_count));

y_fit = polyval(p_mdeit,sensor_count);
y = time_forward_solve_mdeit;

% Remove outliers
residual = y(:)-y_fit(:);
outliers = abs(residual) > 3*std(residual);

sensor_count_clean = sensor_count(~outliers);
y_clean = y(~outliers);
std_clean = std_forward_solve_mdeit(~outliers);
p_mdeit_clean = polyfit(sensor_count_clean,y_clean,1);

y_plot = polyval(p_mdeit_clean,x);
y_fit = polyval(p_mdeit_clean,sensor_count_clean);

SS_res = sum((y_clean(:) - y_fit(:)).^2);
SS_tot = sum((y_clean(:) - mean(y_clean(:))).^2);

R2 = 1 - (SS_res / SS_tot);


errorbar(sensor_count_clean,y_clean,std_clean,'o','MarkerSize',5,'Color',colors(1,:))
plot(x,y_plot,'.','Color',colors(2,:),'LineWidth',2);

% plot(sensor_count(out_id),time_forward_solve_mdeit(out_id),'k*')
hold off
msg = strcat('$ R^2 = ',num2str(R2),'$');
legend('$1$-axis MDEIT',msg,'interpreter','latex','location','southeast');

xlabel('$M$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;

title('MDEIT forward solve','Interpreter','latex')
set(gca,'YScale','log');
set(gca,'XScale','log');




%% FUNCTION
function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
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