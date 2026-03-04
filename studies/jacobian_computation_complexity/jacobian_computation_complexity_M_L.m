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

%% Define the characteristic scales in SI units

%% Model parameters 
z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);

model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

model_parameters.maxsz =  1e-1/l0;

% A material defines the geometry of that material inside the domain, not
% its properties, like conductivity.
model_parameters.material = struct();

% Stimulation pattern was not the default. For now, edit it manually
current_amplitude = 2.4e-3/I0;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used
options = {};

%% Assign the parameters for several models (should create utility functions for this)
num_of_repetitions_eit = 45;

num_of_repetitions_mdeit = 5;

background_conductivity = 3.28e-1/sigma0;  %page 163 mentions a saline solution (NaCl + water) at 0.2% mass concentration, but can't find data for that conductivity, check notes

min_maxsz = 0.5e-3/l0;
max_maxsz = 0.5e-3/l0;
n_steps = 1;


%% Sweep model_parameters over multiple cylindrical radius

% num_of_electrodes_per_ring_array = 4:16; %because of skip 2 pattern, we have to start at 4
% num_of_rings_array = 1:6;

num_of_electrodes_per_ring_array = [4:16]; 
num_of_rings_array = 1:5;

[E, R] = ndgrid(num_of_electrodes_per_ring_array, num_of_rings_array);
combinations = [E(:)'; R(:)'];
num_of_sensors = combinations(1,:).*combinations(2,:);
combinations = [combinations;num_of_sensors];

% for some reason, the following combinations make netgen stuck, so avoid them 
% [15,3,45]
% [13,5,65]

remove ={...
    [15,3,45],...
    [13,5,65]};

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

nnz_A_vector = zeros(numel(model_parameters_array),1);
n_nodes_vector = zeros(numel(model_parameters_array),1);

mesh_created_successfully = true;
maxTime = 2;

% % Start parallel pool if not already running
% pool = gcp('nocreate'); 
% if isempty(pool)
%     pool = parpool; 
% end
% 
% % Run EIDORS startup on all workers
% spmd
%     % Extract just the folder
%     script_folder = 'C:\Users\RuiMartins\Desktop\MDEIT-using-Matlab-main\studies\jacobian_computation_complexity';
% 
%     % Dont output model folder here, re-use the one defined outside par
%     % pool
%     prepare_workspace(script_folder);
% end
% 
% clc;

mesh_created_successfully = true;

fprintf('Computing Jacobian times\n');
for n = 1:numel(model_parameters_array)
    fprintf('Iteration %i of %i\n',n,numel(model_parameters_array));
    model_parameters = model_parameters_array(n);

    if mesh_created_successfully
        fmdl = fmdls_all{n};
        stimulation = mk_stim_patterns(numel(fmdl.electrode),1,inj,meas,options,current_amplitude);
        fmdl.stimulation = stimulation;

        n_elem_vector(n) = size(fmdl.elems,1);

        imgh = mk_image_mdeit(fmdl,background_conductivity);

        temp1 = zeros(1,num_of_repetitions_mdeit);
        temp2 = zeros(1,num_of_repetitions_eit);

        for t = 1:num_of_repetitions_mdeit
            tic
            r1 = fwd_solve_mdeit(imgh);
            temp1(t) = toc;
        end

        for t = 1:num_of_repetitions_eit
            tic
            r2 = fwd_solve(imgh);
            temp2(t) = toc;
        end

        time_forward_solve_mdeit(n) = sum(temp1)/num_of_repetitions_mdeit;
        std_forward_solve_mdeit(n) = std(temp1);
        time_forward_solve_eit(n) = sum(temp2)/num_of_repetitions_eit;
        std_forward_solve_eit(n) = std(temp2);

        t_eit = 0;
        temp = zeros(1,num_of_repetitions_eit);

        for t = 1:num_of_repetitions_eit
            tic;

            J_EIT = calc_jacobian(imgh);
            t_eit = t_eit + toc;
            temp(t) = toc;
        end

        times_eit(n) = t_eit/num_of_repetitions_eit;
        std_eit(n) = std(temp);

        lambdatimesdAdp = @(lambda) computeLambdaTimesDaDp(imgh,lambda);
        A = @(sigma) M(imgh,sigma);

        n_nodes_vector(n) = size(A(imgh.elem_data),1);
        nnz_A_vector(n) = nnz(A(imgh.elem_data));

        t_mdeit = 0;
        temp = zeros(1,num_of_repetitions_mdeit);
        for t = 1:num_of_repetitions_mdeit
            tic;
            J_MDEIT = ...
                calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',3);
            t_mdeit = t_mdeit + toc;
            temp(t) = toc;
        end

        times_mdeit(n) = t_mdeit/num_of_repetitions_mdeit;
        std_mdeit(n) = std(temp);
    end

    % file_name = strcat(script_folder,'/data/data');
    %
    % save(file_name,...
    %     "times_mdeit","times_eit","std_mdeit","std_eit",...
    %     "time_forward_solve_eit","std_forward_solve_eit",...
    %     "time_forward_solve_mdeit","std_forward_solve_mdeit",...
    %     "n_elem_vector","nnz_A_vector","n_nodes_vector");
end

fprintf('Done!\n')

fprintf('Time MDEIT: %.2d +- %.2d\n',times_mdeit(end),std_mdeit(end))


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
figure;
cla;

subplot(1,2,1)
hold on
errorbar(electrode_count.^2,times_eit,std_eit,'o','MarkerSize',5,'Color',colors(3,:))

p_eit = polyfit(electrode_count.*electrode_count,times_eit,1);
x = linspace(min(electrode_count.*sensor_count),max(electrode_count.*sensor_count));

y_plot = polyval(p_eit,x);
y_fit = polyval(p_eit,electrode_count.*electrode_count);
y = times_eit;

SS_res = sum((y(:) - y_fit(:)).^2);
SS_tot = sum((y(:) - mean(y(:))).^2);

R2 = 1 - (SS_res / SS_tot);

plot(x,y_plot,'.','Color',colors(4,:),'LineWidth',2);
hold off
msg = strcat('$ R^2 = ',num2str(R2),'$');
legend('EIT',msg,'interpreter','latex');

xlabel('$L^2$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;

 
% set(gca,'YScale','log');
% set(gca,'XScale','log');

subplot(1,2,2)
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
legend('$1$-axis MDEIT',msg,'interpreter','latex');

% set(gca,'YScale','log');
% set(gca,'XScale','log');

xlabel('$M \times L$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;




%% FUNCTION 



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


%% FUNCTION 
function avg_neigh = average_node_neighbours(fmdl)
    % fmdl.nodes : [n_nodes x dim]
    % fmdl.elems : [n_elems x n_vertex_per_elem]

    elems = fmdl.elems;
    n_nodes = size(fmdl.nodes, 1);

    % Build adjacency list
    neighbours = cell(n_nodes,1);

    for el = 1:size(elems,1)
        verts = elems(el,:);
        % all unique unordered node pairs in this element
        pairs = nchoosek(verts,2);

        % add each pair to adjacency
        for k = 1:size(pairs,1)
            i = pairs(k,1);
            j = pairs(k,2);
            neighbours{i}(end+1) = j;
            neighbours{j}(end+1) = i;
        end
    end

    % remove duplicates and count
    deg = zeros(n_nodes,1);
    for i = 1:n_nodes
        deg(i) = numel(unique(neighbours{i}));
    end

    % compute average number of neighbours
    avg_neigh = mean(deg);

    fprintf('Average number of node neighbours: %.3f\n', avg_neigh);
end