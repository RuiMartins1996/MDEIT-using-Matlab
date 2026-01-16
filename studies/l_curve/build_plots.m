clc; clear all; close all;
%% Setup EIDORS

% Get the full path of the current script
fullpath = mfilename('fullpath');
% Extract just the folder
script_folder = fileparts(fullpath);
cd(script_folder)

eidors_folder = setupEidors("..\..\");
clc;

run("globalParameters.m")


%% Fetch data
cd(script_folder)
file_name = 'data/data_E_4_R_2_M_8_noiseless';

var = load(file_name);

l_curve_struct = var.l_curve_struct;
model_struct = var.model_struct;

lambda_vec = l_curve_struct.lambda_vec;

residual_norms_mdeit1 = l_curve_struct.s_mdeit1.residual_norms;
x_norms_mdeit1 = l_curve_struct.s_mdeit1.x_norms;
optimal_id_mdeit1 = l_curve_struct.s_mdeit1.optimal_id;

residual_norms_mdeit3 = l_curve_struct.s_mdeit3.residual_norms;
x_norms_mdeit3 = l_curve_struct.s_mdeit3.x_norms;
optimal_id_mdeit3 = l_curve_struct.s_mdeit3.optimal_id;


%% Plot results (1-axis MDEIT)
figure
subplot(1,2,1);
hold on
plot(residual_norms_mdeit1,x_norms_mdeit1,'b.-','MarkerSize',10);

xlabel('$|| r(\sigma_{\lambda}) ||_2$','Interpreter','latex');
ylabel('$|| \sigma_{\lambda} ||_2$','Interpreter','latex');

plot(residual_norms_mdeit1(optimal_id_mdeit1),x_norms_mdeit1(optimal_id_mdeit1),'.','MarkerSize',20,'Color',colors(1,:))
str = strcat('$\lambda = ',num2str(lambda_vec(optimal_id_mdeit1)),'$');
text(1.05*residual_norms_mdeit1(optimal_id_mdeit1),x_norms_mdeit1(optimal_id_mdeit1),str,'Interpreter','latex')
box on
grid on;grid minor
set(gca,'YScale','log')
set(gca,'XScale','log')

title('L-curve ($1$-axis MDEIT)','Interpreter','latex')

subplot(1,2,2);
plot(lambda_vec,residual_norms_mdeit1,'b.-','MarkerSize',10);
ylabel('$|| r(\sigma_{\lambda}) ||_2$','Interpreter','latex');
xlabel('$\lambda$','Interpreter','latex');

box on
grid on;grid minor
set(gca,'YScale','log')
set(gca,'XScale','log')

%% Plot results (3-axis MDEIT)

figure
hold on
plot(residual_norms_mdeit3,x_norms_mdeit3,'b.-','MarkerSize',10);

xlabel('$|| r(\sigma_{\lambda}) ||_2$','Interpreter','latex');
ylabel('$|| \sigma_{\lambda} ||_2$','Interpreter','latex');

plot(residual_norms_mdeit3(optimal_id_mdeit3),x_norms_mdeit3(optimal_id_mdeit3),'.','MarkerSize',20,'Color',colors(1,:))
str = strcat('$\lambda = ',num2str(lambda_vec(optimal_id_mdeit3)),'$');
text(1.05*residual_norms_mdeit3(optimal_id_mdeit3),x_norms_mdeit3(optimal_id_mdeit3),str,'Interpreter','latex')
box on
grid on;grid minor
set(gca,'YScale','log')
set(gca,'XScale','log')

title('L-curve ($1$-axis MDEIT)','Interpreter','latex')