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

data_circular_file_name =  'data/data_circular_E_8_R_4_M_48.mat';
data_square_file_name =  'data/data_square_E_8_R_4_M_48.mat';
data_triangle_file_name =  'data/data_triangle_E_8_R_4_M_48.mat';

if isfile(data_circular_file_name)
    var = load(data_circular_file_name);
    s_circular = var.s_circular;
else
    error('File not found');
end

if isfile(data_square_file_name)
    var = load(data_square_file_name);
    s_square = var.s_square;
else
    error('File not found');
end

if isfile(data_triangle_file_name)
    var = load(data_triangle_file_name);
    s_triangle = var.s_triangle;
else
    error('File not found');
end

%% Plot
markerSize = 5;

min_y = min([s_circular.condition_number_array(:);s_square.condition_number_array(:);s_triangle.condition_number_array(:)]);
max_y = max([s_circular.condition_number_array(:);s_square.condition_number_array(:);s_triangle.condition_number_array(:)]);

all_s = [s_circular,s_square,s_triangle];
all_titles={'Circular','Square','Triangle'};

figure

for i = 1:numel(all_s)

    s = all_s(i);

    subplot(1,3,i)
    hold on
    for d = 1:3
        plot(s.sensor_radius,s.condition_number_array(:,d),'Color',colors(d,:),'Marker','o','LineStyle','-','MarkerSize',markerSize);
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

    ylim([min_y,max_y])
    xlim([min(s.sensor_radius),max(s.sensor_radius)])
    grid on;grid minor;box on;

    title(all_titles{i})

end

