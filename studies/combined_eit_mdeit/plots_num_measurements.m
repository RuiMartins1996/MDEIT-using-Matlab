clc;clear all;close all;

% colors = [228,26,28;... %Colors for figure representation
%     55,126,184;...
%     77,175,74;...
%     152,78,163;...
%     255,127,0;...
%     255,255,51;...
%     166,86,40;...
%     247,129,191;
%     202,178,214;
%     106,61,154]/255;


colors = [228,26,28;
55,126,184;
77,175,74;
152,78,163;
255,127,0;
255,255,51;
166,86,40;
247,129,191]/255;

% Convert the bright yellow to a more darker brightness
c = [255 255 51] / 255;   % normalize
hsv = rgb2hsv(c);

hsv(3) = 0.7 * hsv(3);    % reduce brightness (Value)

c_darker = round(hsv2rgb(hsv) * 255);

colors(6,:) = c_darker/255;

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

data_folder = './data';

file_name_eit = strcat(data_folder,'/singular_values_eit.mat');
file_name_mdeit_x = strcat(data_folder,'/singular_values_mdeit_x.mat');
file_name_mdeit_y = strcat(data_folder,'/singular_values_mdeit_y.mat');
file_name_mdeit_z = strcat(data_folder,'/singular_values_mdeit_z.mat');
file_name_mdeit_3 = strcat(data_folder,'/singular_values_mdeit_3.mat');
file_name_aug = strcat(data_folder,'/singular_values_aug.mat');

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
data_aug = load_data(file_name_aug);

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


data_vec = {data_eit,data_mdeit_x,data_mdeit_y,data_mdeit_z,data_mdeit_3,data_aug};
name_vec = {'eit','mdeit-x','mdeit-y','mdeit-z','mdeit-3','eit+mdeit'};

%% Plot condition number

min_num_measurements = inf;
max_num_measurements = -inf;

number_of_elements = zeros(numel(data_vec),numel(data_vec{1}.num_of_elements));

figure
hold on
set(gca, 'Color', [0.95 0.95 0.95])

for i = 1:length(data_vec)
    number_of_elements(i,:) = data_vec{1}.num_of_elements;

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
box on;
set(gca,'YScale','log')
set(gca,'XScale','log')

xlim([min_num_measurements,max_num_measurements])

legend(name_vec,'Interpreter','latex','Location','northwest');

xlabel('Number of measurements','Interpreter','latex')
ylabel("$\kappa$",'Interpreter','latex')

%% Plot rank
min_num_measurements = inf;
max_num_measurements = -inf;


show_ids = [1 4 5 6]; %skip showing data for mdeit-x and mdeit-y, since its the same as for mdeit-z

all_number_of_sensors_sorted = cell(length(data_vec),1);
all_num_measurements_sorted = cell(length(data_vec),1);
all_percentages = cell(length(data_vec),1);

all_ranks_sorted = cell(length(data_vec),1);

ids_valid = ...
    false(length(data_vec),numel(data_vec{1}.num_of_measurements)); %checks where rank is bigger than number of elements

for i = 1:length(data_vec)
    
    % Which entries are valid <=> the rank is smaller than the number of
    % elements
    ids_valid(i,:) = data_vec{i}.rank<data_vec{i}.num_of_elements;
    
    min_num_measurements = min(...
        min_num_measurements,...
        min(data_vec{i}.num_of_measurements(ids_valid(i,:))));
    max_num_measurements = max(...
        max_num_measurements,...
        max(data_vec{i}.num_of_measurements(ids_valid(i,:))));
     
    num_of_measurements_valid = data_vec{i}.num_of_measurements(ids_valid(i,:));
    valid_ranks = data_vec{i}.rank(ids_valid(i,:));
    
    valid_num_of_sensors = data_vec{i}.num_of_sensors(ids_valid(i,:));
    
    % Sort valid entries
    [num_measurements_sorted,ids] = sort(num_of_measurements_valid);
    
    all_num_measurements_sorted{i} = num_measurements_sorted;
    all_number_of_sensors_sorted{i} = valid_num_of_sensors(ids);
    all_ranks_sorted{i} = valid_ranks(ids);
    
    percentages = valid_ranks(ids)./(num_measurements_sorted)*100;
    all_percentages{i} = percentages;
end


figure
set(gca, 'Color', [0.95 0.95 0.95])

hold on
x = linspace(...
    min(vertcat(all_num_measurements_sorted{:})),...
    max(vertcat(all_num_measurements_sorted{:})));

for i = show_ids
    ids = ids_valid(i,:);
    
    plot(all_num_measurements_sorted{i},all_ranks_sorted{i},'Color',colors(i,:),'Marker',markers{i})
end


name_vec_m = name_vec;
for i = show_ids
    ids = ids_valid(i,:);

    poly{i} = polyfit(log10(all_num_measurements_sorted{i}),log10(all_ranks_sorted{i}),1);
    y = x.^(poly{i}(1))*10^(poly{i}(2));
    
    plot(x,y,'--','Color',colors(i,:))

    name_vec_m{i} = strcat(name_vec{i},' ($m = ',num2str(poly{i}(1)),'$)');
end

% Skip plotting these data points, where the rank is the same as the number
% of elements!

% for i = 1:length(data_vec)
%     plot(all_num_measurements_sorted(i,~ids),all_ranks_sorted(i,~ids),'Color','k','Marker',markers{i})
% end
% hold off

xlabel('Number of measurements','Interpreter','latex')
ylabel('Rank','Interpreter','latex')

set(gca,'YScale','log')
set(gca,'XScale','log')

grid on;grid minor;
box on;

legend(name_vec_m{show_ids},'Interpreter','latex','Location','southeast')


%% Plot of rank percentage w.r.t. number of measurements
figure
set(gca, 'Color', [0.95 0.95 0.95])
hold on
for i = show_ids
    ids = ids_valid(i,:);
    
    plot(all_num_measurements_sorted{i},all_percentages{i},...
        'Color',colors(i,:),'Marker',markers{i},'LineWidth',1)
end

legend(name_vec_m{show_ids},'Interpreter','latex','Location','southeast')

ylim([0,100]);
set(gca,'XScale','log')

grid on;grid minor;
box on;

xlabel('Number of measurements','Interpreter','latex')
ylabel('$\rho (\%)$','Interpreter','latex')

%% Plot of rank w.r.t. number of electrodes

figure
set(gca, 'Color', [0.95 0.95 0.95])
hold on
for i = show_ids
    ids = ids_valid(i,:);
    
    plot(all_number_of_sensors_sorted{i},all_ranks_sorted{i},...
        'Color',colors(i,:),'Marker',markers{i},'LineWidth',1)
end

legend(name_vec_m{show_ids},'Interpreter','latex','Location','southeast')

set(gca,'YScale','log')
set(gca,'XScale','log')

grid on;grid minor;
box on;

xlabel('Number of sensors/electrodes','Interpreter','latex')
ylabel('Rank','Interpreter','latex')

%% Plot of percentages w.r.t. number of electrodes

figure
set(gca, 'Color', [0.95 0.95 0.95])
hold on
for i = show_ids
    ids = ids_valid(i,:);
    
    plot(all_number_of_sensors_sorted{i},all_percentages{i},...
        'Color',colors(i,:),'Marker',markers{i},'LineWidth',1,'MarkerSize',2*(max(show_ids)-i+1))
end

legend(name_vec_m{show_ids},'Interpreter','latex','Location','southeast')

ylim([0,100]);
set(gca,'XScale','log')

grid on;grid minor;
box on;

ylabel('$\rho (\%)$','Interpreter','latex')
xlabel('Number of sensors/electrodes','Interpreter','latex')
