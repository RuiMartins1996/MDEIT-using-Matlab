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
%% Assign the parameters for several models (should create utility functions for this)
num_of_repetitions = 3;

maxsz_reconstruction = 0.03;

background_conductivity = 3.28e-1/sigma0;  %page 163 mentions a saline solution (NaCl + water) at 0.2% mass concentration, but can't find data for that conductivity, check notes

min_maxsz = 0.5e-3/l0;
max_maxsz = 10e-3/l0;
n_steps = 20;
%% Define forward model (2D real tank experiment)

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


anomaly_conductivity = 1e-12/l0; % (0 according to Kai's thesis, but check notes)
anomaly_position = [0,0,35]*1e-3/l0; %  [20,0,35]*1e-3/l0; % pg 172
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

% Sweep model_parameters over multiple cylindrical radius
model_parameters_array = ...
    sweep_model_parameters({'maxsz'},min_maxsz,max_maxsz,n_steps,model_parameters,'log');

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

for n = 1:numel(model_parameters_array)

    model_parameters = model_parameters_array(n);

    [model_parameters,fmdls] = ...
        mk_mdeit_model(model_parameters,model_folder,options);
    
    fmdl = fmdls{1};
    stimulation = mk_stim_patterns(numel(fmdl.electrode),1,inj,meas,options,current_amplitude);
    fmdl.stimulation = stimulation;

    n_elem_vector(n) = size(fmdl.elems,1);

    imgh = mk_image_mdeit(fmdl,background_conductivity);
    
    temp1 = zeros(1,num_of_repetitions);
    temp2 = zeros(1,num_of_repetitions);

    for t = 1:num_of_repetitions
        tic
        r1 = fwd_solve_mdeit(imgh);
        temp1(t) = toc;
        tic
        r2 = fwd_solve(imgh);
        temp2(t) = toc;
    end
    
    time_forward_solve_mdeit(n) = sum(temp1)/num_of_repetitions;
    std_forward_solve_mdeit(n) = std(temp1);
    time_forward_solve_eit(n) = sum(temp2)/num_of_repetitions;
    std_forward_solve_eit(n) = std(temp2);

    t_eit = 0;
    temp = zeros(1,num_of_repetitions);

    for t = 1:num_of_repetitions
        tic;

        J_EIT = calc_jacobian(imgh);
        t_eit = t_eit + toc;
        temp(t) = toc;
    end

    times_eit(n) = t_eit/num_of_repetitions;
    std_eit(n) = std(temp);

    lambdatimesdAdp = @(lambda) computeLambdaTimesDaDp(imgh,lambda);
    A = @(sigma) M(imgh,sigma);
    
    n_nodes_vector(n) = size(A(imgh.elem_data),1);
    nnz_A_vector(n) = nnz(A(imgh.elem_data));
    
    t_mdeit = 0;
    temp = zeros(1,num_of_repetitions);
    for t = 1:num_of_repetitions
        tic;
        J_MDEIT = ...
            calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',3);
        t_mdeit = t_mdeit + toc;
        temp(t) = toc;
    end

    times_mdeit(n) = t_mdeit/num_of_repetitions;
    std_mdeit(n) = std(temp);

    file_name = strcat(script_folder,'/data/data');

    save(file_name,...
        "times_mdeit","times_eit","std_mdeit","std_eit",...
        "time_forward_solve_eit","std_forward_solve_eit",...
        "time_forward_solve_mdeit","std_forward_solve_mdeit",...
        "n_elem_vector","nnz_A_vector","n_nodes_vector");
end

fprintf('Done!\n')

fprintf('Time MDEIT: %.2d +- %.2d\n',times_mdeit(end),std_mdeit(end))


%% PLOTS

figure('Position',[100,100,1000,500]);
cla;

subplot(1,3,1)
hold on
errorbar(n_nodes_vector,time_forward_solve_eit,std_forward_solve_eit,'o-','MarkerSize',5,'Color',colors(2,:))
errorbar(n_nodes_vector,time_forward_solve_mdeit,std_forward_solve_mdeit,'d-','MarkerSize',5,'Color',colors(1,:))
hold off

p_f_mdeit = polyfit(...
    log10(n_nodes_vector),...
    log10(time_forward_solve_mdeit),...
    1);

p_f_eit = polyfit(...
    log10(n_nodes_vector),...
    log10(time_forward_solve_eit),...
    1);

hold on
x = linspace(min(n_nodes_vector),max(n_nodes_vector));
plot(x,10^p_f_eit(2)*x.^p_f_eit(1),'LineStyle','--','Color',colors(2,:))
plot(x,10^p_f_mdeit(2)*x.^p_f_mdeit(1),'LineStyle','--','Color',colors(1,:))
hold off


msg1 = strcat('EIT $(m = ',num2str(p_f_eit(1)),'$)'); 
msg2 = strcat('MDEIT $(m = ',num2str(p_f_mdeit(1)),'$)'); 

legend({msg1,msg2},'Interpreter','latex','Location','northwest')

box on;
grid on;grid minor;

set(gca,'YScale','log');
set(gca,'XScale','log');

xlabel('N','Interpreter','latex')
ylabel('t(s)','Interpreter','latex')

subplot(1,3,2)
hold on
errorbar(nnz_A_vector,time_forward_solve_eit,std_forward_solve_eit,'o-','MarkerSize',5,'Color',colors(2,:))
errorbar(nnz_A_vector,time_forward_solve_mdeit,std_forward_solve_mdeit,'d-','MarkerSize',5,'Color',colors(1,:))
hold off

p_A_mdeit = polyfit(...
    log10(nnz_A_vector),...
    log10(time_forward_solve_mdeit),...
    1);

p_A_eit = polyfit(...
    log10(nnz_A_vector),...
    log10(time_forward_solve_eit),...
    1);

hold on
x = linspace(min(nnz_A_vector),max(nnz_A_vector));
plot(x,10^p_A_eit(2)*x.^p_A_eit(1),'LineStyle','--','Color',colors(2,:))
plot(x,10^p_A_mdeit(2)*x.^p_A_mdeit(1),'LineStyle','--','Color',colors(1,:))
hold off


msg1 = strcat('EIT $(m = ',num2str(p_A_eit(1)),'$)'); 
msg2 = strcat('MDEIT $(m = ',num2str(p_A_mdeit(1)),'$)'); 

legend({msg1,msg2},'Interpreter','latex','Location','northwest')

box on;
grid on;grid minor;

set(gca,'YScale','log');
set(gca,'XScale','log');

xlabel('nnz(A)','Interpreter','latex')
ylabel('t(s)','Interpreter','latex')

subplot(1,3,3)

hold on
plot(n_nodes_vector,nnz_A_vector,'o','Color',colors(3,:));

box on;
grid on;grid minor;

ylabel('nnz(A)','Interpreter','latex')
xlabel('t(s)','Interpreter','latex')

p = polyfit(n_nodes_vector,nnz_A_vector,1);
x = linspace(min(n_nodes_vector),max(n_nodes_vector));
plot(x,p(1).*x+p(2),'Color',colors(3,:),'LineStyle','--')

msg = strcat('$(m = ',num2str(p(1)),'$)'); 
legend({msg},'Interpreter','latex','Location','northwest')
title('6 is the average number of neighbour nodes + 1 for diagonal')

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