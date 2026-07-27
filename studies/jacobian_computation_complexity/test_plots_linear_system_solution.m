clc;clear all; close all;

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

%%

file_name = strcat(script_folder,'/data/data_linear_system_solvers.mat');

if isfile(file_name)
    loaded = load(file_name, 'T');
    T = loaded.T;  % retrieve the table
else
    error('No data file found named %s\n',file_name);
end


%% Fetch data from table

electrode_count = T.num_electrodes_per_ring.*T.num_rings;
sensor_count = T.num_sensors;

times_pcg = T.times_pcg;
times_ldl = T.times_ldl;
times_backslash = T.times_backslash;
times_pcg_1_e_10 = T.times_pcg_1_e_10;
times_pcg_1_e_12 = T.times_pcg_1_e_12;

std_pcg = T.std_pcg;
std_ldl = T.std_ldl;
std_backslash = T.std_backslash;
std_pcg_1_e_10 = T.std_pcg_1_e_10;
std_pcg_1_e_12 =  T.std_pcg_1_e_12;

n_elems = T.n_elems;


%% PLOTS
figure;
cla;

fit_data_points = 2:numel(n_elems);

[n_elems_sorted,sorted_ids] = sort(n_elems);

times_pcg_sorted = times_pcg(sorted_ids);
times_ldl_sorted = times_ldl(sorted_ids);
times_backslash_sorted = times_backslash(sorted_ids);
times_pcg_1_e_10_sorted = times_pcg_1_e_10(sorted_ids);
times_pcg_1_e_12_sorted = times_pcg_1_e_12(sorted_ids);

std_pcg_sorted = std_pcg(sorted_ids);
std_ldl_sorted = std_ldl(sorted_ids);
std_backslash_sorted = std_backslash(sorted_ids);
std_pcg_1_e_10_sorted = std_pcg_1_e_10(sorted_ids);
std_pcg_1_e_12_sorted = std_pcg_1_e_12(sorted_ids);

hold on
errorbar(n_elems_sorted,times_pcg_sorted,std_pcg_sorted,'o','MarkerSize',5,'Color',colors(3,:))
errorbar(n_elems_sorted,times_ldl_sorted,std_ldl_sorted,'d','MarkerSize',5,'Color',colors(5,:))
errorbar(n_elems_sorted,times_backslash_sorted,std_backslash_sorted,'p','MarkerSize',5,'Color',colors(7,:))
errorbar(n_elems_sorted,times_pcg_1_e_10_sorted,std_pcg_1_e_10_sorted,'h','MarkerSize',5,'Color',colors(9,:))
errorbar(n_elems_sorted,times_pcg_1_e_12_sorted,std_pcg_1_e_12_sorted,'h','MarkerSize',5,'Color',colors(10,:))


hold off

p_pcg = polyfit(...
    log10(n_elems_sorted(fit_data_points)),...
    log10(times_pcg_sorted(fit_data_points)),...
    1);

p_ldl = polyfit(...
    log10(n_elems_sorted(fit_data_points)),...
    log10(times_ldl_sorted(fit_data_points)),...
    1);

p_backslash = polyfit(...
    log10(n_elems_sorted(fit_data_points)),...
    log10(times_backslash_sorted(fit_data_points)),...
    1);

p_pcg_1_e_10 = polyfit(...
    log10(n_elems_sorted(fit_data_points)),...
    log10(times_pcg_1_e_10_sorted(fit_data_points)),...
    1);

p_pcg_1_e_12 = polyfit(...
    log10(n_elems_sorted(fit_data_points)),...
    log10(times_pcg_1_e_12_sorted(fit_data_points)),...
    1);

% We're checking a fit of the type t = 10^(b)*K^m

hold on
x = linspace(min(n_elems(fit_data_points)),max(n_elems(fit_data_points)));

plot(x,10^p_pcg(2)*x.^p_pcg(1),'LineStyle','--','Color',colors(3,:))
plot(x,10^p_ldl(2)*x.^p_ldl(1),'LineStyle','--','Color',colors(5,:))
plot(x,10^p_backslash(2)*x.^p_backslash(1),'LineStyle','--','Color',colors(7,:))
plot(x,10^p_pcg_1_e_10(2)*x.^p_pcg_1_e_10(1),'LineStyle','--','Color',colors(9,:))
plot(x,10^p_pcg_1_e_12(2)*x.^p_pcg_1_e_12(1),'LineStyle','--','Color',colors(10,:))


hold off

% msg1 = strcat('EIT $(m = ',num2str(p_eit(1)),'$)'); 
msg2 = strcat('pcg - tol = 1e-5 $(t \sim K^{',num2str(p_pcg(1),2),'}$)'); 
msg3 = strcat('ldl $(t \sim K^{',num2str(p_ldl(1),2),'}$)'); 
msg4 = strcat('backslash $(t \sim K^{',num2str(p_backslash(1),2),'}$)'); 
msg5 = strcat('pcg - tol = 1e-10 $(t \sim K^{',num2str(p_pcg_1_e_10(1),2),'}$) | r = ',num2str(max(T.rel_residual_pcg_1_e_10))); 
msg6 = strcat('pcg - tol = 1e-12 $(t \sim K^{',num2str(p_pcg_1_e_12(1),2),'}$)'); 

legend({msg2,msg3,msg4,msg5,msg6},'Interpreter','latex','Location','northwest')

box on;
grid on;grid minor;

set(gca,'YScale','log');
set(gca,'XScale','log');

min_x = min(min(n_elems_sorted(fit_data_points)),10^floor(log10(min(n_elems_sorted(fit_data_points)))));
max_x = max(max(n_elems_sorted(fit_data_points)),10^ceil(log10(max(n_elems_sorted(fit_data_points)))));

min_temp = min([times_pcg_sorted(fit_data_points);times_ldl_sorted(fit_data_points);times_backslash_sorted(fit_data_points)]);
min_y = min(...
    [min_temp,...
    10^floor(log10(min_temp))]...
    );
max_temp = max([times_pcg_sorted(fit_data_points);times_ldl_sorted(fit_data_points);times_backslash_sorted(fit_data_points)]);
max_y = max(...
    [max_temp,...
    10^ceil(log10(max_temp))]...
    );
xlabel('$K$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')

axis([min_x,max_x,min_y,max_y])

title('Adjoint linear system computation time log-scale','Interpreter','latex')

