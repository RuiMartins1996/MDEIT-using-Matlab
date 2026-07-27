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

data_file_name=  'data/data_E_4_R_4_M_16.mat';

data_cylindrical_file_name =  'data/data_cylindrical_E_4_R_4_M_16.mat';

if isfile(data_file_name)
    s = load(data_file_name);
    s = s.s;
else
    error('File not found');
end


if isfile(data_cylindrical_file_name)
    s_c = load(data_cylindrical_file_name);
    s_c = s_c.s;
else
    error('File not found');
end

r_min = s.r_min;
r_max = s.r_max;
n_steps = s.n_steps;
condition_number_array = s.condition_number_array;

r_min_c = s_c.r_min;
r_max_c = s_c.r_max;
n_steps_c = s_c.n_steps;
condition_number_array_c = s_c.condition_number_array;

%% Plot
markerSize = 5;

min_y = min([condition_number_array(:);condition_number_array_c(:)]);
max_y = max([condition_number_array(:);condition_number_array_c(:)]);

figure
subplot(1,2,1)
hold on
for d = 1:3
    plot(linspace(r_min,r_max,n_steps),condition_number_array(:,d),'Color',colors(d,:),'Marker','o','LineStyle','-','MarkerSize',markerSize);
end

legendStr = {};
legendWords = {'x','y','z'};

for d = 1:3
    legendStr = [legendStr,sprintf('$1$-axis MDEIT: axis %s',legendWords{d})];
end
legend(legendStr,'Interpreter','latex','Location','best')


xlabel('$R$','Interpreter','latex')
ylabel('$\kappa$','Interpreter','latex')

set(gca,'YScale','log')

title('Cylindrical Magnetometer Configuration')

ylim([min_y,max_y])
xlim([r_min,r_max])
grid on;grid minor;box on;

subplot(1,2,2)
hold on
for d = 1:3
    plot(linspace(r_min_c,r_max_c,n_steps_c),condition_number_array_c(:,d),'Color',colors(d,:),'Marker','o','LineStyle','-','MarkerSize',markerSize);
end

legendStr = {};
legendWords = {'r','\theta','z'};
for d = 1:3
    legendStr = [legendStr,sprintf('$1$-axis MDEIT: axis $%s$',legendWords{d})];
end
legend(legendStr,'Interpreter','latex','Location','best')


xlabel('$R$','Interpreter','latex')
ylabel('$\kappa$','Interpreter','latex')

set(gca,'YScale','log')

title('Cylindrical Magnetometer Configuration')

ylim([min_y,max_y])
xlim([r_min_c,r_max_c])
grid on;grid minor;box on;


