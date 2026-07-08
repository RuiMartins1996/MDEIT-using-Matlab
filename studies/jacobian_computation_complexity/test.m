clc; clear all; close all;

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

file_name = strcat(script_folder,'/data/data_unfactorized_factorized.mat');


%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);
%% Assign the parameters for several models (should create utility functions for this)

num_of_repetitions_eit = 5;
num_of_repetitions_mdeit = 3;

maxsz_reconstruction = 0.03;

background_conductivity = 3.28e-1/sigma0;  %page 163 mentions a saline solution (NaCl + water) at 0.2% mass concentration, but can't find data for that conductivity, check notes

%% Define forward model (2D real tank experiment)

model_parameters = create_kai_3d_model_parameters(l0, z0, sigma0, I0);

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

% min_maxsz = model_parameters.radius/10;
% max_maxsz = model_parameters.radius/15;

% min_maxsz = model_parameters.radius/15;
% max_maxsz = model_parameters.radius/20;

% min_maxsz = model_parameters.radius/20;
% max_maxsz = model_parameters.radius/25;

% min_maxsz = model_parameters.radius/25;
% max_maxsz = model_parameters.radius/30;

% min_maxsz = model_parameters.radius/35;
% max_maxsz = model_parameters.radius/40;

% Iteration 5 of this run is not worth it, has only 2229153 eleemnts, but
% last data point had 2153817 . Use below

min_maxsz = model_parameters.radius/50;
max_maxsz = model_parameters.radius/55;
n_steps = 1;

% n_steps = 5;

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
times_mdeit_factorized = zeros(numel(model_parameters_array),1);
times_mdeit_factorized_fwd = zeros(numel(model_parameters_array),1);

std_eit = zeros(numel(model_parameters_array),1);
std_mdeit = zeros(numel(model_parameters_array),1);
std_mdeit_factorized = zeros(numel(model_parameters_array),1);
std_mdeit_factorized_fwd = zeros(numel(model_parameters_array),1);


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
        [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],...
        'VariableNames', { ...
        'num_sensors','num_electrodes_per_ring','num_rings', ...
        'times_mdeit','times_mdeit_fact','times_mdeit_fact_fwd','times_eit', ...
        'std_mdeit','std_mdeit_fact','std_mdeit_fact_fwd','std_eit', ...
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
        fprintf('Working on iteration %i\n',n);
        n_elem_vector(n) = size(fmdl.elems,1);

        fprintf('# elements: %i\n',n_elem_vector(n));


        imgh = mk_image_mdeit(fmdl,background_conductivity);
         
        % FORWARD SOLVE: Compute mean execution time and standard deviation
        [time_forward_solve_mdeit(n),std_forward_solve_mdeit(n)] = ...
            compute_execution_time(@fwd_solve_mdeit, num_of_repetitions_mdeit, imgh);

        [time_forward_solve_eit(n),std_forward_solve_eit(n)] = ...
            compute_execution_time(@fwd_solve, num_of_repetitions_eit, imgh);

        % JACOBIAN SOLVE: Compute mean execution time and standard deviation

        A = @(sigma) M(imgh,sigma);

        lambdatimesdAdp = [];
        s_mat = system_mat_1st_order(imgh);
        E = s_mat.E;
        tic
        F = factorise_symmetric(E);
        fprintf('Factorization time: %2.2f\n',toc);
        
        % Factorization precomputed
        [times_mdeit_factorized(n),std_mdeit_factorized(n)] = ...
            compute_execution_time(@ calc_jacobian_mdeit, num_of_repetitions_mdeit, imgh,imgh.elem_data,lambdatimesdAdp,F,'mdeit1',3);
        
        % EIT
        [times_eit(n),std_eit(n)] = compute_execution_time(@calc_jacobian, num_of_repetitions_eit, imgh);

        % No shennanigans
        [times_mdeit(n),std_mdeit(n)] = ...
            compute_execution_time(@ calc_jacobian_mdeit, num_of_repetitions_mdeit, imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',3);
        
        u = fwd_solve(imgh);
        u = u.volt;
        u_struct.Gx_times_u = imgh.fwd_model.G.Gx*u;
        u_struct.Gy_times_u = imgh.fwd_model.G.Gy*u;
        u_struct.Gz_times_u = imgh.fwd_model.G.Gz*u; 
        
        % Factorization and solution gradients precomputed
        [times_mdeit_factorized_fwd(n),std_mdeit_factorized_fwd(n)] = ...
            compute_execution_time(@ calc_jacobian_mdeit, num_of_repetitions_mdeit, imgh,imgh.elem_data,lambdatimesdAdp,F,'mdeit1',3,u_struct);

        new_row = table( ...
            model_parameters.numOfSensors,model_parameters.numOfElectrodesPerRing,model_parameters.numOfRings,...
            times_mdeit(n), times_mdeit_factorized(n), times_mdeit_factorized_fwd(n), times_eit(n), ...
            std_mdeit(n),std_mdeit_factorized(n),std_mdeit_factorized_fwd(n), std_eit(n), ...
            time_forward_solve_eit(n), std_forward_solve_eit(n), ...
            time_forward_solve_mdeit(n), std_forward_solve_mdeit(n), ...
            n_elem_vector(n),...
            'VariableNames', { ...
            'num_sensors','num_electrodes_per_ring','num_rings', ...
            'times_mdeit','times_mdeit_fact','times_mdeit_fact_fwd','times_eit', ...
            'std_mdeit','std_mdeit_fact','std_mdeit_fact_fwd','std_eit', ...
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

    fprintf("Time: %2.2f (s)\n",times(t));

end



%% PLOTS
figure;
cla;

% fit_data_points = length(n_elem_vector):-1:length(n_elem_vector)-7;

fit_data_points = 1:length(times_mdeit);

hold on
errorbar(n_elem_vector,times_mdeit,std_mdeit,'d-','MarkerSize',5,'Color',colors(1,:))
errorbar(n_elem_vector,times_eit,std_eit,'o-','MarkerSize',5,'Color',colors(2,:))
errorbar(n_elem_vector,times_mdeit_factorized,std_mdeit_factorized,'d-','MarkerSize',5,'Color',colors(3,:))
hold off

msg1 = strcat('MDEIT'); 
msg2 = strcat('EIT'); 
msg3 = strcat('MDEIT - factorized'); 

legend({msg1,msg2,msg3},'Interpreter','latex','Location','northwest')

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

title('Jacobian computation time log-scale','Interpreter','latex')





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

%% FUNCTIONS
function F = factorise_symmetric(A)
    F.kind = 'ldl';
    try
        [F.L,F.D,F.P] = ldl(A,'vector'); 
        F.n = size(A,1);
    catch
        [F.L,F.U,F.pv,F.qv] = lu(A,'vector'); 
        F.kind='lu'; 
        F.n   = size(A,1);
    end
end

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
