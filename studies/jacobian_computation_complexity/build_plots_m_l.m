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

file_name = strcat(script_folder,'/data/data_m_l.mat');

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
std_eit = T.std_eit;
std_mdeit = T.std_mdeit;
time_forward_solve_eit = T.time_forward_solve_eit;
time_forward_solve_mdeit = T.time_forward_solve_mdeit;
std_forward_solve_eit = T.std_forward_solve_eit;
std_forward_solve_mdeit = T.std_forward_solve_mdeit;
n_elems = T.n_elems;


%% PLOTS
% SHOULD PACK THIS INTO A FUNCTION!!!!!!!!
figure;
cla;

% EIT Jacobian computation time plot
subplot(2,2,1)
hold on
errorbar(electrode_count.^2,times_eit,std_eit,'o','MarkerSize',5,'Color',colors(3,:))

% THE LAST THREE POINTS ARE OUTLIERS! BECAUSE OF EIDORS REMOVING CACHE. FOR
% NOW, REMOVE THEM FROM FIT.

[sorted_time_eit,eit_ids] = sort(times_eit);
sorted_electrode_count = electrode_count(eit_ids);

outlier_id = 1 + find(diff(sorted_time_eit) == max(diff(sorted_time_eit))); %detect jump

p_eit_1 = polyfit(sorted_electrode_count(1:outlier_id-1).^2,sorted_time_eit(1:outlier_id-1),1);
p_eit_2 = polyfit(sorted_electrode_count(outlier_id:end).^2,sorted_time_eit(outlier_id:end),1);
x1 = linspace(min(sorted_electrode_count(1:outlier_id-1).^2),max(sorted_electrode_count(1:outlier_id-1).^2));
x2 = linspace(min(sorted_electrode_count(outlier_id:end).^2),max(sorted_electrode_count(outlier_id:end).^2));

y_plot = polyval(p_eit_1,x1);
y_plot_2 = polyval(p_eit_2,x2);

y_fit_1 = polyval(p_eit_1,sorted_electrode_count(1:outlier_id-1).^2);
y_fit_2 = polyval(p_eit_2,sorted_electrode_count(outlier_id:end).^2);

y1 = sorted_time_eit(1:outlier_id-1);
y2 = sorted_time_eit(outlier_id:end);

SS_res = sum((y1(:) - y_fit_1(:)).^2);
SS_tot = sum((y1(:) - mean(y1(:))).^2);
R2 = 1 - (SS_res / SS_tot);

SS_res_2 = sum((y2(:) - y_fit_2(:)).^2);
SS_tot_2 = sum((y2(:) - mean(y2(:))).^2);
R2_2 = 1 - (SS_res_2 / SS_tot_2);


plot(x1,y_plot,'.','Color',colors(4,:),'LineWidth',2);
plot(x2,y_plot_2,'.','Color',colors(6,:),'LineWidth',2);

hold off
msg1 = strcat('linear fit - $R^2 = ',num2str(R2),'$');
msg2 = strcat('linear fit - $R^2 = ',num2str(R2_2),'$');
legend('EIT',msg1,msg2,'interpreter','latex','location','southeast');

xlabel('$L^2$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;

title('EIT Jacobian computation','Interpreter','latex')

set(gca,'YScale','log');
set(gca,'XScale','log');
axis square


% MDEIT Jacobian computation time plot
subplot(2,2,2)
hold on
errorbar(electrode_count.*sensor_count,times_mdeit,std_mdeit,'d','MarkerSize',5,'Color',colors(5,:))

p_mdeit = polyfit(electrode_count.*sensor_count,times_mdeit,1);
x1 = linspace(min(electrode_count.*sensor_count),max(electrode_count.*sensor_count));

y_plot = polyval(p_mdeit,x1);
y_fit_1 = polyval(p_mdeit,electrode_count.*sensor_count);
y1 = times_mdeit;

SS_res = sum((y1(:) - y_fit_1(:)).^2);
SS_tot = sum((y1(:) - mean(y1(:))).^2);

R2 = 1 - (SS_res / SS_tot);

plot(x1,y_plot,'.','Color',colors(4,:),'LineWidth',2);
hold off
msg1 = strcat('linear fit - $R^2 = ',num2str(R2),'$');
legend('$1$-axis MDEIT',msg1,'interpreter','latex','location','southeast');

set(gca,'YScale','log');
set(gca,'XScale','log');

title('$1$-axis MDEIT Jacobian computation','Interpreter','latex')

xlabel('$M \times L$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;
axis square


% EIT forward computation time plot
subplot(2,2,3)
hold on
errorbar(electrode_count,time_forward_solve_eit,std_forward_solve_eit,'o','MarkerSize',5,'Color',colors(3,:))

p_eit_1 = polyfit(electrode_count,time_forward_solve_eit,1);
x1 = linspace(min(electrode_count),max(electrode_count));

y_plot = polyval(p_eit_1,x1);
y_fit_1 = polyval(p_eit_1,electrode_count);
y1 = time_forward_solve_eit;

SS_res = sum((y1(:) - y_fit_1(:)).^2);
SS_tot = sum((y1(:) - mean(y1(:))).^2);

R2 = 1 - (SS_res / SS_tot);

plot(x1,y_plot,'.','Color',colors(4,:),'LineWidth',2);
hold off
msg1 = strcat('linear fit - $R^2 = ',num2str(R2),'$');
legend('EIT',msg1,'interpreter','latex','location','southeast');

xlabel('$L$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;

title('EIT forward solve','Interpreter','latex')

set(gca,'YScale','log');
set(gca,'XScale','log');
axis square

% MDEIT forward computation time plot
subplot(2,2,4)
hold on
% errorbar(sensor_count,time_forward_solve_mdeit,std_forward_solve_mdeit,'o','MarkerSize',5,'Color',colors(1,:))

p_mdeit = polyfit(sensor_count,time_forward_solve_mdeit,1);
x1 = linspace(min(sensor_count),max(sensor_count));

y_fit_1 = polyval(p_mdeit,sensor_count);
y1 = time_forward_solve_mdeit;

% Remove outliers
residual = y1(:)-y_fit_1(:);
outliers = abs(residual) > 3*std(residual);

sensor_count_clean = sensor_count(~outliers);
y_clean = y1(~outliers);
std_clean = std_forward_solve_mdeit(~outliers);
p_mdeit_clean = polyfit(sensor_count_clean,y_clean,1);

y_plot = polyval(p_mdeit_clean,x1);
y_fit_1 = polyval(p_mdeit_clean,sensor_count_clean);

SS_res = sum((y_clean(:) - y_fit_1(:)).^2);
SS_tot = sum((y_clean(:) - mean(y_clean(:))).^2);

R2 = 1 - (SS_res / SS_tot);

errorbar(sensor_count_clean,y_clean,std_clean,'o','MarkerSize',5,'Color',colors(5,:))
plot(x1,y_plot,'.','Color',colors(4,:),'LineWidth',2);

% plot(sensor_count(out_id),time_forward_solve_mdeit(out_id),'k*')
hold off
msg1 = strcat('linear fit - $R^2 = ',num2str(R2),'$');
legend('$1$-axis MDEIT',msg1,'interpreter','latex','location','southeast');

xlabel('$M$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;

title('MDEIT forward solve','Interpreter','latex')
set(gca,'YScale','log');
set(gca,'XScale','log');
axis square



%% 4 figures from the 2x2 grid above, separated
figure('Name','EIT Jacobian Computation')
hold on
errorbar(electrode_count.^2,times_eit,std_eit,'o','MarkerSize',5,'Color',colors(3,:))

% THE LAST THREE POINTS ARE OUTLIERS! BECAUSE OF EIDORS REMOVING CACHE. FOR
% NOW, REMOVE THEM FROM FIT.

[sorted_time_eit,eit_ids] = sort(times_eit);
sorted_electrode_count = electrode_count(eit_ids);

outlier_id = 1 + find(diff(sorted_time_eit) == max(diff(sorted_time_eit))); %detect jump

p_eit_1 = polyfit(sorted_electrode_count(1:outlier_id-1).^2,sorted_time_eit(1:outlier_id-1),1);
p_eit_2 = polyfit(sorted_electrode_count(outlier_id:end).^2,sorted_time_eit(outlier_id:end),1);
x1 = linspace(min(sorted_electrode_count(1:outlier_id-1).^2),max(sorted_electrode_count(1:outlier_id-1).^2));
x2 = linspace(min(sorted_electrode_count(outlier_id:end).^2),max(sorted_electrode_count(outlier_id:end).^2));

y_plot = polyval(p_eit_1,x1);
y_plot_2 = polyval(p_eit_2,x2);

y_fit_1 = polyval(p_eit_1,sorted_electrode_count(1:outlier_id-1).^2);
y_fit_2 = polyval(p_eit_2,sorted_electrode_count(outlier_id:end).^2);

y1 = sorted_time_eit(1:outlier_id-1);
y2 = sorted_time_eit(outlier_id:end);

SS_res = sum((y1(:) - y_fit_1(:)).^2);
SS_tot = sum((y1(:) - mean(y1(:))).^2);
R2 = 1 - (SS_res / SS_tot);

SS_res_2 = sum((y2(:) - y_fit_2(:)).^2);
SS_tot_2 = sum((y2(:) - mean(y2(:))).^2);
R2_2 = 1 - (SS_res_2 / SS_tot_2);


plot(x1,y_plot,'.','Color',colors(4,:),'LineWidth',2);
plot(x2,y_plot_2,'.','Color',colors(6,:),'LineWidth',2);

hold off
msg1 = strcat('linear fit - $R^2 = ',num2str(R2),'$');
msg2 = strcat('linear fit - $R^2 = ',num2str(R2_2),'$');
legend('EIT',msg1,msg2,'interpreter','latex','location','southeast');

xlabel('$L^2$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;
axis square

title('EIT Jacobian computation','Interpreter','latex')

% set(gca,'YScale','log');
% set(gca,'XScale','log');

figure('Name','1-axis MDEIT Jacobian Computation')
hold on
errorbar(electrode_count.*sensor_count,times_mdeit,std_mdeit,'d','MarkerSize',5,'Color',colors(5,:))

p_mdeit = polyfit(electrode_count.*sensor_count,times_mdeit,1);
x1 = linspace(min(electrode_count.*sensor_count),max(electrode_count.*sensor_count));

y_plot = polyval(p_mdeit,x1);
y_fit_1 = polyval(p_mdeit,electrode_count.*sensor_count);
y1 = times_mdeit;

SS_res = sum((y1(:) - y_fit_1(:)).^2);
SS_tot = sum((y1(:) - mean(y1(:))).^2);

R2 = 1 - (SS_res / SS_tot);

plot(x1,y_plot,'.','Color',colors(4,:),'LineWidth',2);
hold off
msg1 = strcat(' linear fit - $R^2 = ',num2str(R2),'$');
legend('$1$-axis MDEIT',msg1,'interpreter','latex','location','southeast');

xlabel('$M \times L$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')

title('$1$-axis MDEIT Jacobian computation','Interpreter','latex')

% set(gca,'YScale','log');
% set(gca,'XScale','log');

box on;
grid on;grid minor;
axis square



figure('Name','EIT forward solve')
hold on
errorbar(electrode_count,time_forward_solve_eit,std_forward_solve_eit,'o','MarkerSize',5,'Color',colors(3,:))

p_eit_1 = polyfit(electrode_count,time_forward_solve_eit,1);
x1 = linspace(min(electrode_count),max(electrode_count));

y_plot = polyval(p_eit_1,x1);
y_fit_1 = polyval(p_eit_1,electrode_count);
y1 = time_forward_solve_eit;

SS_res = sum((y1(:) - y_fit_1(:)).^2);
SS_tot = sum((y1(:) - mean(y1(:))).^2);

R2 = 1 - (SS_res / SS_tot);

plot(x1,y_plot,'.','Color',colors(4,:),'LineWidth',2);
hold off
msg1 = strcat(' linear fit - $R^2 = ',num2str(R2),'$');
legend('EIT',msg1,'interpreter','latex','location','southeast');

xlabel('$L$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;

title('EIT forward solve','Interpreter','latex')

set(gca,'YScale','log');
set(gca,'XScale','log');
axis square



figure('Name','MDEIT forward solve')
hold on
% errorbar(sensor_count,time_forward_solve_mdeit,std_forward_solve_mdeit,'o','MarkerSize',5,'Color',colors(1,:))

p_mdeit = polyfit(sensor_count,time_forward_solve_mdeit,1);
x1 = linspace(min(sensor_count),max(sensor_count));

y_fit_1 = polyval(p_mdeit,sensor_count);
y1 = time_forward_solve_mdeit;

% Remove outliers
residual = y1(:)-y_fit_1(:);
outliers = abs(residual) > 3*std(residual);

sensor_count_clean = sensor_count(~outliers);
y_clean = y1(~outliers);
std_clean = std_forward_solve_mdeit(~outliers);
p_mdeit_clean = polyfit(sensor_count_clean,y_clean,1);

y_plot = polyval(p_mdeit_clean,x1);
y_fit_1 = polyval(p_mdeit_clean,sensor_count_clean);

SS_res = sum((y_clean(:) - y_fit_1(:)).^2);
SS_tot = sum((y_clean(:) - mean(y_clean(:))).^2);

R2 = 1 - (SS_res / SS_tot);

errorbar(sensor_count_clean,y_clean,std_clean,'o','MarkerSize',5,'Color',colors(5,:))
plot(x1,y_plot,'.','Color',colors(4,:),'LineWidth',2);

% plot(sensor_count(out_id),time_forward_solve_mdeit(out_id),'k*')
hold off
msg1 = strcat('linear fit - $R^2 = ',num2str(R2),'$');
legend('$1$-axis MDEIT',msg1,'interpreter','latex','location','southeast');

xlabel('$M$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')
box on;
grid on;grid minor;

title('MDEIT forward solve','Interpreter','latex')
set(gca,'YScale','log');
set(gca,'XScale','log');
axis square

%% Execution time of the forward and Jacobian solves for EIT and MDEIT on the same graphs

figure('Name','Jacobian Computation')

hold on
errorbar(electrode_count,times_eit,std_eit,'o','MarkerSize',5,'Color',colors(3,:))
errorbar(electrode_count,times_mdeit,std_mdeit,'d','MarkerSize',5,'Color',colors(5,:))

hold off
% Legend, title and axis parameters
legend('EIT','$1$-axis MDEIT','interpreter','latex','location','northwest');

title('Jacobian computation','Interpreter','latex')

xlabel('$L$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')

box on;
grid on;grid minor;
% set(gca,'YScale','log');
% set(gca,'XScale','log');
axis square

figure('Name','Forward Solve')

hold on
errorbar(electrode_count,time_forward_solve_eit,std_forward_solve_eit,'o','MarkerSize',5,'Color',colors(3,:))
errorbar(electrode_count,time_forward_solve_mdeit,std_forward_solve_mdeit,'d','MarkerSize',5,'Color',colors(5,:))

hold off
% Legend, title and axis parameters
legend('EIT','$1$-axis MDEIT','interpreter','latex','location','southeast');

title('Forward computation','Interpreter','latex')

xlabel('$L$','Interpreter','latex');
ylabel('$t(s)$','Interpreter','latex')

box on;
grid on;grid minor;
set(gca,'YScale','log');
set(gca,'XScale','log');
axis square