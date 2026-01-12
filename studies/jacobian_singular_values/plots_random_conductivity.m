clc;clear all;close all;

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

markers = {'x','d','o','s','*','+'};
marker_size = 3;


%% Prepare workspace
% Get the full path of the current script
fullpath = mfilename('fullpath');
% Extract just the folder
script_folder = fileparts(fullpath);
cd(script_folder);

% Have to add the functions path manually so prepare_workspace runs
parent_folder = fileparts(script_folder);
grandparent_folder =fileparts(parent_folder);
addpath(genpath(fullfile(grandparent_folder,'functions')));

model_folder = prepare_workspace(script_folder);

data_folder = './data_random_conductivity';

file_name_eit = strcat(data_folder,'/singular_values_eit.mat');
file_name_mdeit_x = strcat(data_folder,'/singular_values_mdeit_x.mat');
file_name_mdeit_y = strcat(data_folder,'/singular_values_mdeit_y.mat');
file_name_mdeit_z = strcat(data_folder,'/singular_values_mdeit_z.mat');
file_name_mdeit_3 = strcat(data_folder,'/singular_values_mdeit_3.mat');

files = {file_name_eit,file_name_mdeit_x,file_name_mdeit_y,file_name_mdeit_z,file_name_mdeit_3};
%% Check if files are missing

assert(all(cellfun(@(f) exist(f,'file'), files)), 'Files are missing');

%% Load data 

function data = load_data(file_name)

% Load the existing old data if it exists
var = load(file_name);

if numel(fieldnames(var)) ~= 1
    error('Expected number of structure fields to be 1');
end

field_names = fieldnames(var);

data = var.(field_names{1});
end

data_eit = load_data(file_name_eit);
data_mdeit_x = load_data(file_name_mdeit_x);
data_mdeit_y = load_data(file_name_mdeit_y);
data_mdeit_z = load_data(file_name_mdeit_z);
data_mdeit_3 = load_data(file_name_mdeit_3);

%% Compute the condition number of the jacobians

function condition_number_array = compute_condition_number(data)

num_of_data_points = size(data.singular_values,2);

max_singular_values = zeros(num_of_data_points,1);
min_singular_values = zeros(num_of_data_points,1);

for i = 1:num_of_data_points
    max_singular_values(i) = max(data.singular_values(1:data.rank(i),i));
    min_singular_values(i) = min(data.singular_values(1:data.rank(i),i));
end

condition_number_array = max_singular_values./min_singular_values;

end


data_vec = {data_eit,data_mdeit_x,data_mdeit_y,data_mdeit_z,data_mdeit_3};
name_vec = {'eit','mdeit-x','mdeit-y','mdeit-z','mdeit-3'};

%% Plot singular values
hold on
for i = 1:length(data_vec)
    for j = 1:length(data_vec{i}.random_seed)
        rank_ij = data_vec{i}.rank(j);
        if j==1
        plot(1:rank_ij,data_vec{i}.singular_values(1:rank_ij,j),...
            'Marker',markers{i},'MarkerSize',marker_size,'Color',colors(i,:));
        else
            plot(1:rank_ij,data_vec{i}.singular_values(1:rank_ij,j),...
            'Marker',markers{i},'MarkerSize',marker_size,'Color',colors(i,:),'HandleVisibility','off');
        end
    end
end
hold off
grid on;grid minor;

set(gca,'YScale','log')

legend(name_vec,'Interpreter','latex');

xlabel([],'Interpreter','latex')
ylabel("$\sigma$",'Interpreter','latex')

%% Singular value difference

data_vec = {data_eit,data_mdeit_x,data_mdeit_y,data_mdeit_z,data_mdeit_3};
name_vec = {'eit','mdeit-x','mdeit-y','mdeit-z','mdeit-3'};

for i = 1:numel(data_vec)
    subplot(1,numel(data_vec),i)
    rank_ij = data_vec{i}.rank(1);

    absolute_differences = abs(data_vec{i}.singular_values(1:rank_ij,1)-data_vec{i}.singular_values(1:rank_ij,2:end));
    relative_differences = 100*absolute_differences./abs(data_vec{i}.singular_values(1:rank_ij,1));

    hold on
    for j = 1:size(absolute_differences,2)
        plot(1:rank_ij,relative_differences(:,j),...
            'Marker',markers{i},'MarkerSize',marker_size,'Color',colors(i,:));
    end
    hold off

    grid on;grid minor;

    xlabel([],'Interpreter','latex')
    ylabel("$\frac{|\sigma_1-\sigma_j|}{\sigma_1} (\%)$",'Interpreter','latex')
    title(name_vec{i},'Interpreter','latex')
end
