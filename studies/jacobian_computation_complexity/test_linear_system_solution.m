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

file_name = strcat(script_folder,'/data/data_linear_system_solvers.mat');

check_entries = {'pcg','ldl','backslash','pcg_1_e_10','pcg_1_e_12'};

%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);
%% Assign the parameters for several models (should create utility functions for this)

num_of_repetitions = 3;

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

min_maxsz = model_parameters.radius/10;
max_maxsz = model_parameters.radius/15;

% min_maxsz = model_parameters.radius/15;
% max_maxsz = model_parameters.radius/20;

% min_maxsz = model_parameters.radius/20;
% max_maxsz = model_parameters.radius/25;

% min_maxsz = model_parameters.radius/25;
% max_maxsz = model_parameters.radius/30;

% min_maxsz = model_parameters.radius/35;
% max_maxsz = model_parameters.radius/40;

min_maxsz = model_parameters.radius/15;
max_maxsz = model_parameters.radius/30;

n_steps = 15;

%% Sweep model_parameters over multiple cylindrical radius
model_parameters_array = ...
    sweep_model_parameters({'maxsz'},min_maxsz,max_maxsz,n_steps,model_parameters,'log');

model_parameters_array(9) = []; % THIS FAILS TO MESH
%% First of all, build or load every forward model

% If this does not terminate, then we can do something different

fmdls_all = cell(1,numel(model_parameters_array));

for n = 1:numel(model_parameters_array)
    model_parameters = model_parameters_array(n);

    [model_parameters,fmdls] = ...
        mk_mdeit_model(model_parameters,model_folder,options);
    fmdls_all{n} = fmdls{1};
end

%% Allocate variables

time_pcg = zeros(numel(model_parameters_array),1);
time_ldl = zeros(numel(model_parameters_array),1);
time_backslash = zeros(numel(model_parameters_array),1);
time_pcg_1_e_10 = zeros(numel(model_parameters_array),1);
time_pcg_1_e_12 = zeros(numel(model_parameters_array),1);


time_eidors_left_divide = zeros(numel(model_parameters_array),1);

std_pcg = zeros(numel(model_parameters_array),1);
std_ldl = zeros(numel(model_parameters_array),1);
std_backslash = zeros(numel(model_parameters_array),1);
std_pcg_1_e_10  = zeros(numel(model_parameters_array),1);
std_pcg_1_e_12  = zeros(numel(model_parameters_array),1);


std_eidors_left_divide = zeros(numel(model_parameters_array),1);

n_elem_vector = zeros(numel(model_parameters_array),1);

%% Initialize table
n_models = numel(model_parameters_array);

num_sensors = zeros(n_models,1);
num_electrodes_per_ring = zeros(n_models,1);
num_rings = zeros(n_models,1);
n_elem = zeros(n_models,1);

for n = 1:n_models
    mp = model_parameters_array(n);

    num_sensors(n) = mp.numOfSensors;
    num_electrodes_per_ring(n) = mp.numOfElectrodesPerRing;
    num_rings(n) = mp.numOfRings;

    % assumes fmdls_all already built
    n_elem(n) = size(fmdls_all{n}.elems, 1);

end

T = table( ...
    num_sensors, ...
    num_electrodes_per_ring, ...
    num_rings, ...
    n_elem, ...
    'VariableNames', { ...
    'num_sensors', ...
    'num_electrodes_per_ring', ...
    'num_rings', ...
    'n_elems', ...
    } ...
);

%% Get the data table or copy from existing data
base_names = { ...
    'num_sensors','num_electrodes_per_ring','num_rings', ...
    'n_elems'};

residual_names = strcat('rel_residual_', check_entries);
times_names = strcat('times_', check_entries);
std_names   = strcat('std_',   check_entries);

var_names = [base_names, times_names, std_names,residual_names];

% THIS IS NOT COPYING OLD DATA CORRECTLY, MOST LIKELY!!!!!!!!!!!!!!!!!
% ---------------------------------------------------------
% CASE 1: file exists → preserve old + extend schema
% ---------------------------------------------------------
if isfile(file_name)

    loaded = load(file_name, 'T');
    T_old = loaded.T;
    n = height(T_old);

    T = table();

    % initialize full schema with NaNs
    for i = 1:numel(var_names)
        T.(var_names{i}) = nan(n,1);
    end

    % copy old data where available
    common_vars = intersect(var_names, T_old.Properties.VariableNames, 'stable');
    T(:, common_vars) = T_old(:, common_vars);

% ---------------------------------------------------------
% CASE 2: no file → initialize empty run-ready table
% ---------------------------------------------------------
else

    n = numel(model_parameters_array); % or your planned size

    T = table();

    for i = 1:numel(var_names)
        T.(var_names{i}) = nan(n,1);
    end

    % Fill know entries
    for row = 1:n
        model_parameters = model_parameters_array(row);

        [model_parameters,fmdls] = ...
            mk_mdeit_model(model_parameters,model_folder,options);
        fmdl = fmdls{1};

        % Fill known entries
        T.num_sensors(row) = model_parameters.numOfSensors;
        T.num_electrodes_per_ring(row) = model_parameters.numOfElectrodesPerRing;
        T.num_rings(row) = model_parameters.numOfRings;
        T.n_elems(row) = size(fmdl.elems,1);
    end
end

%% Cull table if there are repeated entries
[~, ia] = unique(T(:, base_names), 'rows', 'stable');
T = T(ia, :);

%% Find missing entries in the table

if ~isempty(T) && height(T) > 0
    prefix = 'times_';

    col_names = strcat(prefix, check_entries);
    col_names = intersect(col_names, T.Properties.VariableNames, 'stable');

    data = T(:, col_names);
    mask = ismissing(data);

    [row_idx, col_idx] = find(mask);

    missing_cols = col_names(col_idx);
    missing_n_elems = T.n_elems(row_idx);
else 
    error('Should not be here')
end



%% Find indices of fmdl_all for which missing_n_elems exists 

for n = 1:numel(fmdls_all)
    n_elem_fmdl(n) = size(fmdls_all{n}.elems,1);
end

row_idx_new = [];
for n = 1:numel(row_idx)
     this_missing_n_elems = T.n_elems(row_idx(n));
     cols(n) = col_names(col_idx(n));
    
     row_idx_new = [row_idx_new ;find(n_elem_fmdl == this_missing_n_elems)];
end

%%% THIS IS WRONG. WHAT I WANT IS THE INDEX IN T WHERE THE MISSING n_elem
%%% IS LOCATED, NOT THE INDEX IN fmdls_all ..., but will do for now because
%%% the n_elems in table are sorted

%% Loop through missing entries 

for n = 1:numel(row_idx_new)
    
    row = row_idx_new(n);
    col = cols{n};
    
    % Load model
    model_parameters = model_parameters_array(row);
    fmdl = fmdls_all{row};
    stimulation = mk_stim_patterns(numel(fmdl.electrode),1,inj,meas,options,current_amplitude);
    fmdl.stimulation = stimulation;

    fprintf('Working on iteration %i | %s\n ',row,col);
    fprintf('# elements: %i\n',size(fmdl.elems,1));

    % % Fill known entries
    % T.num_sensors(row) = model_parameters.numOfSensors;
    % T.num_electrodes_per_ring(row) = model_parameters.numOfElectrodesPerRing;
    % T.num_rings(row) = model_parameters.numOfRings;
    % T.n_elems(row) = n_elem_vector(row);

    % Compute missing entries
    imgh = mk_image_mdeit(fmdl,background_conductivity);
    A = @(sigma) M(imgh,sigma);
    select_sensor_axis = 3;



    % Compute Gamma matrices
    imgh = compute_gamma_matrices(imgh);


    switch col
        case 'times_pcg'
            tol = 1e-5;
            
            % Compute error:
            lambda = solve_pcg(imgh,A,select_sensor_axis,tol);
            rel_residual = compute_residual(imgh,lambda);
            
            % Compute time
            [time_pcg(n),std_pcg(n)] = ...
                compute_execution_time(@solve_pcg, num_of_repetitions, imgh,A,select_sensor_axis,tol);
            
            T.rel_residual_pcg(row) = rel_residual;
            T.times_pcg(row) = time_pcg(n);
            T.std_pcg(row)   = std_pcg(n);
        case 'times_ldl'

            % Compute error:
            lambda = solve_ldl(imgh,select_sensor_axis);
            rel_residual = compute_residual(imgh,lambda);

            [time_ldl(n),std_ldl(n)] = ...
                compute_execution_time(@solve_ldl, num_of_repetitions, imgh,select_sensor_axis);

            T.rel_residual_ldl(row) = rel_residual;
            T.times_ldl(row) = time_ldl(n);
            T.std_ldl(row)   = std_ldl(n);
        case 'times_backslash'

            % Compute error:
            lambda = solve_backslash(imgh,A,select_sensor_axis);
            rel_residual = compute_residual(imgh,lambda);

            [time_backslash(n),std_backslash(n)] = ...
                compute_execution_time(@solve_backslash, num_of_repetitions, imgh,A,select_sensor_axis);

            T.rel_residual_backslash(row) = rel_residual;
            T.times_backslash(row) = time_backslash(n);
            T.std_backslash(row)   = std_backslash(n);

        case 'times_pcg_1_e_10'
            tol = 1e-10;
            
            lambda = solve_pcg(imgh,A,select_sensor_axis,tol);
            rel_residual = compute_residual(imgh,lambda);

            [time_pcg_1_e_10(n),std_pcg_1_e_10(n)] = ...
                compute_execution_time(@solve_pcg, num_of_repetitions, imgh,A,select_sensor_axis,tol);

            T.rel_residual_pcg_1_e_10(row) = rel_residual;
            T.times_pcg_1_e_10(row) = time_pcg_1_e_10(n);
            T.std_pcg_1_e_10(row)   = std_pcg_1_e_10(n);

        case 'times_pcg_1_e_12'
            tol = 1e-12;

            lambda = solve_pcg(imgh,A,select_sensor_axis,tol);
            rel_residual = compute_residual(imgh,lambda);

            [time_pcg_1_e_12(n),std_pcg_1_e_12(n)] = ...
                compute_execution_time(@solve_pcg, num_of_repetitions, imgh,A,select_sensor_axis,tol);

            T.rel_residual_pcg_1_e_12(row) = rel_residual;
            T.times_pcg_1_e_12(row) = time_pcg_1_e_12(n);
            T.std_pcg_1_e_12(row)   = std_pcg_1_e_12(n);

        otherwise
            error('Unknown column');
    end

    % Save at each iteration so no progress is lost
    save(file_name, "T");
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

function rel_residual = compute_residual(img,lambda)

Gamma = img.Gamma3;
GammaT = Gamma.';

 A = @(sigma) M(img,sigma);

residual = A(img.elem_data)*lambda+GammaT;

norm_residual = zeros(size(residual,2),1);
rel_residual = zeros(size(residual,2),1);

for j = 1:size(residual,2)
    norm_residual(j) = norm(residual(:,j),2);
    rel_residual(j) = norm_residual(j)./norm(-GammaT(:,j));
end

rel_residual = max(rel_residual);

end



function [norm_residual,rel_residual,abs_error,rel_error] = compute_error(func,varargin)

% Check if varargin{1} is an EIDORS img
img = varargin{1};
assert(strcmp(img.type,'image'));

% Solve with function
lambda = func(varargin{:});

% Solve with PCG with very accurate tolerance:
A = @(sigma) M(img,sigma);
lambda_exact = solve_pcg(img,A,3,1e-12);

% Compute Gamma matrices
img = compute_gamma_matrices(img);
Gamma = img.Gamma3;

n_nodes = size(img.fwd_model.nodes,1);
num_sensors = numel(img.fwd_model.sensors);
num_electrodes = numel(img.fwd_model.electrode);

GammaT = Gamma.';

residual = A(img.elem_data)*lambda+GammaT;

norm_residual = zeros(size(residual,2),1);
rel_residual = zeros(size(residual,2),1);

for n = 1:size(residual,2)
    norm_residual(n) = norm(residual(:,n),2);
    rel_residual(n) = norm_residual(n)./norm(-GammaT(:,n));
end

lambda_ldl = solve_ldl(img,3);
lambda_backslash = solve_backslash(img,A,3);

residual = A(img.elem_data)*lambda+GammaT;
residual_ldl = A(img.elem_data)*lambda_ldl+GammaT;
residual_backslash = A(img.elem_data)*lambda_backslash+GammaT;



abs_error = norm(lambda(:)-lambda_exact(:),'inf');
rel_error = norm(lambda(:)-lambda_exact(:),'inf')/norm(lambda_exact(:),'inf');

end
%% FUNCTION: solve_pcg

function lambda = solve_pcg(img,A,select_sensor_axis,tol)

n_nodes = size(img.fwd_model.nodes,1);
num_sensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

switch select_sensor_axis
    case 1
        Gamma = img.Gamma1;
    case 2
        Gamma = img.Gamma2;
    case 3
        Gamma = img.Gamma3;
    otherwise
        error('here')
end

% Solve the adjoint problem for each sensor to get lambda vectors
GammaT = Gamma.';

A_matrix = A(img.elem_data);

lambda = zeros(n_nodes,num_sensors);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

num_elements = numel(img.elem_data);
parfor m = 1:num_sensors
    [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),tol,num_elements,Mfun,Nfun);
end

end

function lambda = solve_ldl(img,select_sensor_axis)

n_nodes = size(img.fwd_model.nodes,1);
num_sensors = numel(img.fwd_model.sensors);
num_electrodes = numel(img.fwd_model.electrode);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

% Factorize
s_mat = system_mat_1st_order(img);
E = s_mat.E;
F = factorise_symmetric(E);

switch select_sensor_axis
    case 1
        Gamma = img.Gamma1;
    case 2
        Gamma = img.Gamma2;
    case 3
        Gamma = img.Gamma3;
    otherwise
        error('here')
end

% Solve the adjoint problem for each sensor to get lambda vectors
GammaT = Gamma.';

rhs = sparse(n_nodes+num_electrodes,num_sensors);
rhs(1:end-num_electrodes,:) = -GammaT;
lambda = solve_fact_multiple_rhs(F,rhs);
lambda = lambda(1:end-num_electrodes,:);

end

function lambda = solve_backslash(img,A,select_sensor_axis)

% Compute Gamma matrices
img = compute_gamma_matrices(img);

switch select_sensor_axis
    case 1
        Gamma = img.Gamma1;
    case 2
        Gamma = img.Gamma2;
    case 3
        Gamma = img.Gamma3;
    otherwise
        error('here')
end

% Solve the adjoint problem for each sensor to get lambda vectors
GammaT = Gamma.';

A_matrix = A(img.elem_data);

lambda = A_matrix \ (-GammaT);

end

function lambda = solve_eidors_left_divide(img,select_sensor_axis,tol)
n_nodes = size(img.fwd_model.nodes,1);
num_sensors = numel(img.fwd_model.sensors);
num_electrodes = numel(img.fwd_model.electrode);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

% Factorize
s_mat = system_mat_1st_order(img);
E = s_mat.E;

switch select_sensor_axis
    case 1
        Gamma = img.Gamma1;
    case 2
        Gamma = img.Gamma2;
    case 3
        Gamma = img.Gamma3;
    otherwise
        error('here')
end

% Solve the adjoint problem for each sensor to get lambda vectors
GammaT = Gamma.';

rhs = sparse(n_nodes+num_electrodes,num_sensors);
rhs(1:end-num_electrodes,:) = -GammaT;

lambda = left_divide( E, rhs,tol);

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

%% FUNCTIONS
function F = factorise_symmetric(A)
    F.kind = 'ldl';
    try
        [F.L,F.D,F.P] = ldl(A,'vector'); 
        F.n = size(A,1);
    catch
        error('Couldnt do it')
        % [F.L,F.U,F.pv,F.qv] = lu(A,'vector'); 
        % F.kind='lu'; 
        % F.n   = size(A,1);
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



function X = solve_fact_multiple_rhs(F, rhs)

    switch F.kind

        case 'ldl'
            % Permute RHS (each column independently)
            rp = rhs(F.P, :);

            % LDL solves (all column-wise)
            y  = F.L \ rp;
            z  = F.D \ y;
            w  = F.L' \ z;

            % Allocate full solution matrix
            X = zeros(F.n, size(rhs,2));

            % Unpermute rows
            X(F.P, :) = w;

        case 'lu'
            % Row permutation of RHS
            y = rhs(F.pv, :);

            % Triangular solves
            z = F.L \ y;
            w = F.U \ z;

            % Allocate solution
            X = zeros(F.n, size(rhs,2));

            % Column permutation recovery
            X(F.qv, :) = w;

        otherwise
            error('Unknown factorisation kind.');
    end
end
