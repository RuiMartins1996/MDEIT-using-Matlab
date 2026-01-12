clc; clear all; close all;
%% Prepare workspace

% Get the full path of the current script
fullpath = mfilename('fullpath');

% Extract just the folder
script_folder = fileparts(fullpath);

cd(script_folder);

cd("..\..\");

addpath(genpath("functions"));
addpath(genpath("libraries"));

run("globalParameters.m")

%Setup EIDORS
eidors_folder = setupEidors(cd);

cd(script_folder);

clc;
%% Fetch data
file_name = 'data/data_E_4_R_4_M_16_mode_mdeit3_SNR_100';

var = load(file_name);

model_struct = var.model_struct;
gcv_struct = var.gcv_struct;
loocv_struct = var.loocv_struct;

% Number of lines in the 1-axis Jacobian matrix
n = numel(model_struct.imdl.fwd_model.sensors)*numel(model_struct.imdl.fwd_model.stimulation);

if numel(fieldnames(model_struct))~=0
    imdl = model_struct.imdl;
    
    imgi = model_struct.imgi;
    imgi = mk_image_mdeit(imgi);
    
    imgh = model_struct.imgh;
    imgh = mk_image_mdeit(imgh);
else
    warning('model_struct is empty')
end

if numel(fieldnames(gcv_struct))~=0
    V_mu = gcv_struct.V_mu;
    mu_min = gcv_struct.mu_min;
    gcv_dsigma = gcv_struct.conductivity_optimal_reg_par;
    mu_vector = gcv_struct.mu_vector;

    optimal_id = gcv_struct.optimal_id;
else
    warning('gcv_struct is empty');
end


%% Plot ground truth and optimal mu reconstruction
figure('Position',[200,200,1000,500])
if numel(fieldnames(model_struct))~=0
    subplot(1,2,1)
    show_fem(imgi);
    eidors_colourbar(imgi);
    title('Ground truth','interpreter','latex');
    box on

    subplot(1,2,2)
    img_output = mk_image_mdeit(imdl.fwd_model,1);
    img_output.elem_data =  gcv_dsigma;
    show_fem(img_output);
    eidors_colourbar(img_output);
    title_str = strcat('Reconstruction $\lambda =',num2str(mu_min),'$');
    title(title_str,'interpreter','latex')
    box on
end

%% Plot the GCV cost function
if numel(fieldnames(gcv_struct))~=0
    figure;
    hold on
    
    range = max(V_mu)-min(V_mu);

    plot(mu_vector,V_mu,'b.');
    plot(mu_vector(optimal_id),V_mu(optimal_id),'r.','MarkerSize',10);

    str = strcat('$n\lambda = ',sprintf('%.3g',mu_min),'$');
    text(mu_vector(optimal_id),V_mu(optimal_id)+range/100,str,'interpreter','latex')


    set(gca,'YScale','log');
    set(gca,'XScale','log');
    title('GCV cost function','interpreter','latex');
    
    xlabel('$\lambda$','Interpreter','latex')
    ylabel('$V(\lambda)$','Interpreter','latex')
    grid on; grid minor;
    box on
    legend('$V(\lambda)$','interpreter','latex')
end