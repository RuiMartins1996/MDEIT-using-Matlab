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
s = load('data/data.mat');

r_min = s.r_min;
r_max = s.r_max;
n_steps = s.n_steps;
condition_number_array = s.condition_number_array;

%% Plot
markerSize = 5;

figure
hold on
for d = 1:3
    plot(linspace(r_min,r_max,n_steps),condition_number_array(:,d),'Color',colors(d,:),'Marker','o','LineStyle','-','MarkerSize',markerSize);
end

legendStr = {};
for d = 1:3
    legendStr = [legendStr,sprintf('$1$-axis MDEIT: axis %s',num2str(d))];
end
legend(legendStr,'Interpreter','latex','Location','best')


xlabel('$R$','Interpreter','latex')
ylabel('$\kappa$','Interpreter','latex')

set(gca,'YScale','log')

title('Cylindrical Magnetometer Configuration')

xlim([r_min,r_max])
grid on;grid minor;box on;

