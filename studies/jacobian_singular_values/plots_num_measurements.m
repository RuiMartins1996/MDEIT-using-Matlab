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

data_folder = './data_num_measurements';

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

%% Plot condition number

min_num_measurements = inf;
max_num_measurements = -inf;

hold on
for i = 1:length(data_vec)
    min_num_measurements = min(min_num_measurements,min(data_vec{i}.num_of_measurements));
    max_num_measurements = max(max_num_measurements,max(data_vec{i}.num_of_measurements));
    
    [num_measurements_sorted,ids] = sort(data_vec{i}.num_of_measurements);
    condition_number_vec = compute_condition_number(data_vec{i});
    condition_number_sorted = condition_number_vec(ids);

    plot(num_measurements_sorted,condition_number_sorted,...
        'Marker',markers{i},'MarkerSize',marker_size,'Color',colors(i,:));
end
hold off
grid on;grid minor;

set(gca,'YScale','log')

xlim([min_num_measurements,max_num_measurements])

legend(name_vec,'Interpreter','latex');

xlabel("$N_{meas}$",'Interpreter','latex')
ylabel("$\kappa$",'Interpreter','latex')

%% Plot rank
figure
min_num_measurements = inf;
max_num_measurements = -inf;

all_num_measurements_sorted = ...
    zeros(length(data_vec),numel(data_vec{1}.num_of_measurements));

all_percentages = ...
    zeros(length(data_vec),numel(data_vec{1}.num_of_measurements));
all_ranks_sorted = ...
    zeros(length(data_vec),numel(data_vec{1}.num_of_measurements));
hold on
for i = 1:length(data_vec)

    min_num_measurements = min(min_num_measurements,min(data_vec{i}.num_of_measurements));
    max_num_measurements = max(max_num_measurements,max(data_vec{i}.num_of_measurements));
    
    [num_measurements_sorted,ids] = sort(data_vec{i}.num_of_measurements);

    all_num_measurements_sorted(i,:) = num_measurements_sorted;
    all_ranks_sorted(i,:) = data_vec{i}.rank(ids);

    if i == 5
        percentages = data_vec{i}.rank(ids)./(num_measurements_sorted)*100;
        all_percentages(i,:) = percentages;
    else
        percentages = data_vec{i}.rank(ids)./(num_measurements_sorted)*100;
        all_percentages(i,:) = percentages;
    end

    % plot(num_measurements_sorted,rank_sorted,...
    %     'Marker',markers{i},'MarkerSize',marker_size,'Color',colors(i,:));
end

subplot(1,2,1)
hold on
b1 = bar(num_measurements_sorted, all_num_measurements_sorted','FaceColor','flat','HandleVisibility','off');
b2 = bar(num_measurements_sorted, all_ranks_sorted' );
hold off

for i = 1:length(data_vec)
    b1(i).CData = 1;
end

for i = 1:length(data_vec)
    b2(i).CData = colors(i,:);
end

hold off
grid on;grid minor;

set(gca,'YScale','log')

ylim([1,max_num_measurements]);

legend(name_vec,'Interpreter','latex');

xlabel("$N_{meas}$",'Interpreter','latex')
ylabel("rank$(J)$",'Interpreter','latex')

subplot(1,2,2)
hold on
b = bar(num_measurements_sorted, all_percentages' );
hold off

for i = 1:length(data_vec)
    b(i).CData = colors(i,:);
end

hold off
grid on;grid minor;

ylim([0,100]);

legend(name_vec,'Interpreter','latex');

xlabel("$N_{meas}$",'Interpreter','latex')
ylabel("$\frac{rank(J)}{N_{meas}} (\%)$",'Interpreter','latex')

