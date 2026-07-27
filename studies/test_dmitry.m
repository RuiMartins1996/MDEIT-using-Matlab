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


file_name_eit = strcat(data_folder,'/singular_values_eit.mat');
file_name_mdeit_x = strcat(data_folder,'/singular_values_mdeit_x.mat');
file_name_mdeit_y = strcat(data_folder,'/singular_values_mdeit_y.mat');
file_name_mdeit_z = strcat(data_folder,'/singular_values_mdeit_z.mat');
file_name_mdeit_3 = strcat(data_folder,'/singular_values_mdeit_3.mat');
file_name_aug = strcat(data_folder,'/singular_values_aug.mat');

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

beta = V0/H0;

% Set the minimum sensor radius
sensor_radius_0 = 1.01;
rmin = sensor_radius_0;
rmax = 3;

%% Problem parameters

% Default 3D model
model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);
% with some changes
model_parameters.maxsz = model_parameters.radius/5;
model_height = model_parameters.height/2;
model_parameters.numOfElectrodesPerRing = 8;
model_parameters.numOfRings = 4;
model_parameters.numOfSensors = 8*4;
model_parameters.anomaly.type = 'spherical';
model_parameters.material.type = 'spherical';
model_parameters.material.name = 'sphere_anomaly';
model_parameters.sensorRadius = 3.5*model_parameters.radius;


% model_parameters.material.position = [0.5 0 model_parameters.height/2];
model_parameters.material.position = [0 0 model_parameters.height/2];
% model_parameters.material.radius = model_parameters.material.radius/2;

background_conductivity = 3.28e-1/sigma0;
anomaly_conductivity = 1e-12/sigma0;

max_sizes = [model_parameters.radius/4;model_parameters.radius/5; model_parameters.radius/10;model_parameters.radius/15];

maxsz_data =  model_parameters.radius/10;

current_amplitude = 2.4e-3/I0;

inj = [0 3]; %skip 2 pattern (pg 172)
meas = [0 3]; %for EIT, skip2 measurement protocol was used

%% Generate two forward models with different number of elements

fmdl_array = cell(numel(max_sizes),1);

for i = 1:numel(max_sizes)
    model_parameters_temp = rmfield(model_parameters,'material');
    model_parameters_temp.maxsz = max_sizes(i);
    
    [~,fmdl] = mk_mdeit_model(model_parameters_temp,model_folder);
    
    fmdl_array{i} = fmdl{1};
end
%% Assign stimulation patterns


stimulation = ...
    mk_stim_patterns(...
    model_parameters.numOfElectrodesPerRing,...
    model_parameters.numOfRings,...
    inj,meas,{'meas_current'},current_amplitude);

for i = 1:numel(max_sizes)
fmdl_array{i}.stimulation = stimulation;
end

%% Make homogeneous images

imgs = cell(numel(max_sizes),1);
for i = 1:numel(imgs)
    imgs{i} =  mk_image_mdeit(fmdl_array{i},background_conductivity,sprintf('Model %i',i));
end

%% Compute Jacobians

ranks = zeros(numel(max_sizes),1);
s = cell(numel(max_sizes),1);

for i = 1:numel(max_sizes)
    fprintf('Working on model %i\n',i);

    imgh = imgs{i};
    A = @(sigma) M(imgh,sigma);

    J_mdeit = calc_jacobian_mdeit(imgh,imgh.elem_data,[],A,'mdeit3');
    
    s_mdeit = svd(J_mdeit);
    
    rank_mdeit = sum(s_mdeit>1e-12);
    
    ranks(i) = rank_mdeit;
    s{i} = s_mdeit;
end



%% Plots

miny = 1e-11;
maxy = 1e-1;

figure
hold on
for i = 1:numel(max_sizes)
    plot(s{i});
end
plot(ones(1,100)*ranks(1),linspace(miny,maxy))
text(ranks(1),1e-2,sprintf('rank = %i',ranks(1)))

kappas = zeros(numel(max_sizes),1);

for i = 1:numel(max_sizes)
    kappas(i) = s{i}(1)/s{i}(ranks(i));
end

msgs = cell(numel(max_sizes),1);

for i = 1:numel(max_sizes)
msgs{i} = strcat(...
    sprintf('Model %i',i),'|',...
    sprintf('K = %i',size(fmdl_array{i}.elems,1)),'|',...
    '$\kappa = ',sprintf('%g%',kappas(i)),'$');
end

legend(msgs,'interpreter','latex','Location','southwest')

set(gca,'YScale','log')
grid on;grid minor;box on

ylim([miny,maxy])


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