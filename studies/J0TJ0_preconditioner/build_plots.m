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

%% Fetch data

data_file_name=  'data/data_file.mat';

if isfile(data_file_name)
    s = load(data_file_name);
    s = s.data_struct;
else
    error('File not found');
end

%%
fields = fieldnames(s);
for k = 1:numel(fields)
    assignin('base', fields{k}, s.(fields{k}));
end

%% Time w.r.t tolerances
figure('Position',[100,100,1000,500]);
for i = 1:numel(lambdas)
    subplot(1,numel(lambdas),i)

    hold on
    errorbar(tolerances,time_lm_no_pre(i,:),std_lm_no_pre(i,:),'b.-');
    errorbar(tolerances,time_lm(i,:),std_lm(i,:),'r.-');
    errorbar(tolerances,time_gn(i,:),std_gn(i,:),'g.-');
    errorbar(tolerances,time_gn_pre(i,:),std_gn_pre(i,:),'k.-');
    hold off
    
    set(gca,'YScale','log')
    legend('LM - No Pre','LM - Pre','GN - No Pre','GN - Pre')
    set(gca,'XScale','log');
    grid on;grid minor;
    box on;
    xlabel('tol','Interpreter','latex');
    ylabel('time(s)','Interpreter','latex')
    title_str = strcat('$\lambda = ',num2str(lambdas(i)),'$');
    title(title_str,'Interpreter','latex')
end

%% 

figure('Position',[100,100,1000,500]);
for i = 1:numel(lambdas)
        subplot(1,numel(lambdas),i)

    hold on
    plot(residual_lm_no_pre(i,:),time_lm_no_pre(i,:),'b.')
    plot(residual_lm(i,:),time_lm(i,:),'r * ')

    hold off

    legend({'LM - No Pre','LM - Pre'})
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    grid on;grid minor;
    box on;
    ylabel('time(s)','Interpreter','latex')
    xlabel('error','Interpreter','latex')
    title_str = strcat('$\lambda = ',num2str(lambdas(i)),'$');
    title(title_str,'Interpreter','latex')
end

%% Plots of number of PCG iterations per LM iteration
figure('Position',[100,100,1000,500]);

names = {'LM - No Pre','LM - Pre','GN - No Pre','GN - Pre'};
markers = {'p','s','o','p','h'};

all_plots = {pcg_iterations_gn,pcg_iterations_gn_pre};

id = 0;
for i = 1:numel(lambdas)
    id = id+1;
    subplot(numel(lambdas),2,id)
    ct = 1;

    % this_plot = pcg_iterations_lm_no_pre;
    this_plot = all_plots{i};
    hold on

    for j = numel(tolerances):-1:1
        marker = markers{exit_flags_lm_no_pre(j)};
        color = colors(j,:);
        plot(this_plot{i,j},'Color',color,'LineWidth',1,'Marker',marker)
        legendStr{ct} = sprintf('tol = %.1d',tolerances(j));
        ct = ct+1;
    end

    hold off
    legend(legendStr)

    title_str = strcat('LM - No Pre $\lambda = ',num2str(lambdas(i)),'$');
    title(title_str,'Interpreter','latex')    

    box on;
    grid on

    ylabel('number of PCG iterations','Interpreter','latex')
    xlabel('number of LM iterations','Interpreter','latex')
    
    id = id+1;
    subplot(numel(lambdas),2,id)
    ct = 1;
    legendStr = {};

    this_plot = pcg_iterations_lm;
    hold on

    for j = numel(tolerances):-1:1
        marker = markers{exit_flags_lm(j)};
        color = colors(j,:);
        plot(this_plot{i,j},'Color',color,'LineWidth',1,'Marker',marker)
        legendStr{ct} = sprintf('tol = %.1d',tolerances(j));
        ct = ct+1;
    end

    hold off
    legend(legendStr)

    title_str = strcat('LM - Pre $\lambda = ',num2str(lambdas(i)),'$');
    title(title_str,'Interpreter','latex')  

    box on;
    grid on

    ylabel('number of PCG iterations','Interpreter','latex')
    xlabel('number of LM iterations','Interpreter','latex')
end

%% Statistics
i = 1;
for j = 1:numel(tolerances)

    a1 = pcg_iterations_lm_no_pre{i,j};
    a2 = pcg_iterations_lm{i,j};


    fprintf('_____________________________\n');
    fprintf('tol = %.2g\n',tolerances(j));
    fprintf('Total number of PCG iterations LM: %i\n',sum(a1));
    fprintf('Total number of PCG iterations LM - Pre: %i\n',sum(a2));
end

for i = 1:numel(tolerances)
    
    a1 = residual_lm_no_pre(i);
    a2 = residual_lm(i);

    fprintf('_____________________________\n');
    fprintf('tol = %.2g\n',tolerances(i));
    fprintf('Residual LM: %i\n',a1);
    fprintf('Residual LM - Pre: %i\n',a2);
end

for i = 1:numel(tolerances)
    
    a1 = time_lm_no_pre(i);
    a2 = time_lm(i);

    b1 = std_lm_no_pre(i);
    b2 = std_lm(i);

    fprintf('_____________________________\n');
    fprintf('tol = %.2g\n',tolerances(i));
    fprintf('Time LM: %.2f +/- %.2f\n',a1,b1);
    fprintf('Time LM - Pre: %.2f +/- %.2f\n',a2,b2);
end
