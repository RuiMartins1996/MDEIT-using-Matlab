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

file_name = strcat(script_folder,'/data/data_unfactorized_factorized.mat');

if isfile(file_name)
    loaded = load(file_name, 'T');
    T = loaded.T;  % retrieve the table
else
    error('No data file found named %s\n',file_name);
end


%% Fetch data from table

electrode_count = T.num_electrodes_per_ring.*T.num_rings;
sensor_count = T.num_sensors;

times_eit = T.times_eit;
times_mdeit = T.times_mdeit;
times_mdeit_factorized = T.times_mdeit_fact;
times_mdeit_factorized_fwd = T.times_mdeit_fact_fwd;
std_eit = T.std_eit;
std_mdeit = T.std_mdeit;
std_mdeit_factorized = T.std_mdeit_fact;
std_mdeit_factorized_fwd = T.std_mdeit_fact_fwd;
time_forward_solve_eit = T.time_forward_solve_eit;
time_forward_solve_mdeit = T.time_forward_solve_mdeit;
std_forward_solve_eit = T.std_forward_solve_eit;
std_forward_solve_mdeit = T.std_forward_solve_mdeit;
n_elems = T.n_elems;


%% PLOTS
figure;
cla;

fit_data_points_mdeit = length(n_elems)-10:length(n_elems);
fit_data_points_eit = length(n_elems)-10:length(n_elems);

[n_elems_sorted,sorted_ids] = sort(n_elems);
times_mdeit_factorized_fwd_sorted = times_mdeit_factorized_fwd(sorted_ids);
times_mdeit_factorized_sorted = times_mdeit_factorized(sorted_ids);
times_mdeit_sorted = times_mdeit(sorted_ids);
times_eit_sorted = times_eit(sorted_ids);

std_eit_sorted = std_eit(sorted_ids);
std_mdeit_sorted = std_mdeit(sorted_ids);
std_mdeit_factorized_sorted = std_mdeit_factorized(sorted_ids);
std_mdeit_factorized_fwd_sorted = std_mdeit_factorized_fwd(sorted_ids);

hold on
errorbar(n_elems_sorted(fit_data_points_eit),...
    times_eit_sorted(fit_data_points_eit),std_eit_sorted(fit_data_points_eit),'o','MarkerSize',5,'Color',colors(3,:))
errorbar(n_elems_sorted(fit_data_points_mdeit),...
    times_mdeit_sorted(fit_data_points_mdeit),std_mdeit_sorted(fit_data_points_mdeit),'d','MarkerSize',5,'Color',colors(5,:))
errorbar(n_elems_sorted(fit_data_points_mdeit),...
    times_mdeit_factorized_sorted(fit_data_points_mdeit),std_mdeit_factorized_sorted(fit_data_points_mdeit),'p','MarkerSize',5,'Color',colors(7,:))
errorbar(n_elems_sorted(fit_data_points_mdeit),...
    times_mdeit_factorized_fwd_sorted(fit_data_points_mdeit),std_mdeit_factorized_fwd(fit_data_points_mdeit),'h','MarkerSize',5,'Color',colors(8,:))

hold off

p_mdeit_fact_fwd = polyfit(...
    log10(n_elems_sorted(fit_data_points_mdeit)),...
    log10(times_mdeit_factorized_fwd_sorted(fit_data_points_mdeit)),...
    1);

p_mdeit_fact = polyfit(...
    log10(n_elems_sorted(fit_data_points_mdeit)),...
    log10(times_mdeit_factorized_sorted(fit_data_points_mdeit)),...
    1);

p_mdeit = polyfit(...
    log10(n_elems_sorted(fit_data_points_mdeit)),...
    log10(times_mdeit_sorted(fit_data_points_mdeit)),...
    1);

p_eit = polyfit(...
    log10(n_elems_sorted(fit_data_points_eit)),...
    log10(times_eit_sorted(fit_data_points_eit)),...
    1);

% We're checking a fit of the type t = 10^(b)*K^m

hold on
x_mdeit = linspace(min(n_elems_sorted(fit_data_points_mdeit)),max(n_elems_sorted(fit_data_points_mdeit)));
x_eit = linspace(min(n_elems_sorted(fit_data_points_eit)),max(n_elems_sorted(fit_data_points_eit)));

plot(x_mdeit,10^p_mdeit(2)*x_mdeit.^p_mdeit(1),'LineStyle','--','Color',colors(5,:))
plot(x_mdeit,10^p_mdeit_fact(2)*x_mdeit.^p_mdeit_fact(1),'LineStyle','--','Color',colors(7,:))
plot(x_mdeit,10^p_mdeit_fact_fwd(2)*x_mdeit.^p_mdeit_fact_fwd(1),'LineStyle','--','Color',colors(8,:))
plot(x_eit,10^p_eit(2)*x_eit.^p_eit(1),'LineStyle','--','Color',colors(3,:))


time_to_compute = 3600;
number_of_elements_for_time_to_compute = 10^((log10(time_to_compute) - p_mdeit(2))/p_mdeit(1));

number_of_elements = 3.2e6; %head mesh number of elements
time_to_compute_for_number_of_elements_mdeit = 10^p_mdeit(2)*number_of_elements.^p_mdeit(1);
time_to_compute_for_number_of_elements_fact = 10^p_mdeit_fact(2)*number_of_elements.^p_mdeit_fact(1);
time_to_compute_for_number_of_elements_fact_fwd = 10^p_mdeit_fact_fwd(2)*number_of_elements.^p_mdeit_fact_fwd(1);

fprintf('Expected exec time for head mesh (3.2M)\n')
fprintf('\t regular | t: %2.2f (min)\n',time_to_compute_for_number_of_elements_mdeit/60)
fprintf('\t fact | t: %2.2f (min)\n',time_to_compute_for_number_of_elements_fact/60)
fprintf('\t fact/fwd | t: %2.2f (min)\n',time_to_compute_for_number_of_elements_fact_fwd/60)

hold off

msg1 = strcat('EIT $(t \sim K^{',num2str(p_eit(1),2),'}$)'); 
msg2 = strcat('MDEIT $(t \sim K^{',num2str(p_mdeit(1),2),'}$)'); 
msg3 = strcat('MDEIT - fact $(t \sim K^{',num2str(p_mdeit_fact(1),2),'}$)'); 
msg4 = strcat('MDEIT - fact/fwd $(t \sim K^{',num2str(p_mdeit_fact(1),2),'}$)'); 


legend({msg1,msg2,msg3,msg4},'Interpreter','latex','Location','northwest')

box on;
grid on;grid minor;

set(gca,'YScale','log');
set(gca,'XScale','log');

min_x = min(min(n_elems_sorted(fit_data_points_mdeit)),10^floor(log10(min(n_elems_sorted(fit_data_points_mdeit)))));
max_x = max(max(n_elems_sorted(fit_data_points_mdeit)),10^ceil(log10(max(n_elems_sorted(fit_data_points_mdeit)))));

min_temp = min([times_mdeit_sorted(fit_data_points_mdeit);times_eit_sorted(fit_data_points_eit)]);
min_y = min(...
    [min_temp,...
    10^floor(log10(min_temp))]...
    );

max_temp = max([times_mdeit_sorted(fit_data_points_mdeit);times_eit_sorted(fit_data_points_eit)]);
max_y = max(...
    [max_temp,...
    10^ceil(log10(max_temp))]...
    );
xlabel('$K$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')

axis([min_x,max_x,min_y,max_y])

title('Jacobian computation time log-scale','Interpreter','latex')

