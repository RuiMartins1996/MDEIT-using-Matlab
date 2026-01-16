
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


% Get the full path of the current script
fullpath = mfilename('fullpath');
% Extract just the folder
script_folder = fileparts(fullpath);
cd(script_folder);

model_folder = prepare_workspace(script_folder);
%% Load data

data_folder = 'data/anomaly_position_0';
% data_folder = 'data';
files = dir(fullfile(data_folder, '*.mat'));

SNR_vector = zeros(numel(files),1);

for k = 1:numel(files)

    file_name = fullfile(data_folder, files(k).name);
    data = load(file_name);

    token = regexp(file_name, 'SNR_(\d+)\.mat', 'tokens');
    SNR_vector(k) = str2num(token{1}{1});

    all_data{k} = data;
end

[SNR_vector,ids] = sort(SNR_vector);
for i =1:numel(all_data)
    all_data_sorted{i} = all_data{ids(i)};
end

all_data = all_data_sorted;

%% Compute GREIT parameters
levels =[inf,inf,0];

% Anomaly position/radius
xyzr = [0.5;0;0;0.3125];

imgi = all_data{1}.imgi;
params_real = eval_GREIT_fig_merit(imgi, xyzr);


for i = 1:numel(all_data)

    imgr_eit = all_data{i}.img_output_eit;
    imgr_mdeit = all_data{i}.img_output_mdeit_1;

    img = all_data{i}.img_output_eit;
    img.elem_data = img.elem_data./all_data{i}.sigma_std_deviation_eit;
    imgr_eit_corrected = img;

    img = all_data{i}.img_output_mdeit_1;
    img.elem_data = img.elem_data./all_data{i}.sigma_std_deviation_mdeit;
    imgr_mdeit_corrected = img;

    % show_slices(imgr_mdeit, levels);
    % imgr_mdeit.calc_colours.npoints = 128;
    % imgr_mdeit.calc_slices.levels=levels;
    %
    % rimg= calc_slices( imgr_mdeit, levels );
    % r_img = mk_mosaic(rimg, 0);
    % c_img = calc_colours( r_img, imgr_mdeit);
    % out_img= reshape(c_img, size(r_img,1), size(r_img,2) ,[]);
    % image(out_img)

    params_eit = eval_GREIT_fig_merit(imgr_eit, xyzr);
    params_mdeit = eval_GREIT_fig_merit(imgr_mdeit, xyzr);

    params_eit_corrected = eval_GREIT_fig_merit(imgr_eit_corrected, xyzr);
    params_mdeit_corrected = eval_GREIT_fig_merit(imgr_mdeit_corrected, xyzr);

    p_eit{i} = params_eit;
    p_mdeit{i} = params_mdeit;
    p_eit_corrected{i} = params_eit_corrected;
    p_mdeit_corrected{i} = params_mdeit_corrected;

    % p_names = {'AR','PE','RES','SD','RNG','AR'};
    % fprintf('MDEIT: \n');
    % for j=1:numel(p_names)
    %     fprintf('Image %s : %3.3g | Real %s : %3.3g \n',...
    %         p_names{j},params_mdeit(j,:), ...
    %         p_names{j},params_real(j,:));
    % end
    %
    % fprintf('EIT: \n');
    % for j=1:numel(p_names)
    %     fprintf('Image %s : %3.3g | Real %s : %3.3g \n',...
    %         p_names{j},params_eit(j,:), ...
    %         p_names{j},params_real(j,:));
    % end

end

figure
p_names = {'AR','PE','RES','SD','RNG','AR'};
choose = [2,3,4,5];
ct = 1;
for j = choose

    subplot(numel(choose),1,ct);

    for i = 1:numel(all_data)
        p1(i) = p_eit{i}(j);
        p2(i) = p_mdeit{i}(j);
        p3(i) = p_eit_corrected{i}(j);
        p4(i) = p_mdeit_corrected{i}(j);
    end

    hold on
    plot(SNR_vector,p1);
    plot(SNR_vector,p2);
    plot(SNR_vector,p3);
    plot(SNR_vector,p4);
    plot(SNR_vector,params_real(j)*ones(size(p1)))
    hold off

    legend('EIT','MDEIT','EIT NBC','MDEIT NBC','Real')
    ylabel(p_names{j},'Interpreter','latex');
    xlabel('SNR','Interpreter','latex')

    ct = ct+1;
end

% voxel_size = 0.05;
% voxel_img = voxelize(imgr,voxel_size);
% mesh(voxel_img.uVoxel)
% greitEIT = computeGreitFigures(voxel_img ,background_conductivity);

%% Plot reconstructions

figure
show_fem_transparent_edges(all_data{1}.imgi)
box on
title('Ground-Truth','Interpreter','latex')

figure('Position',[10,10,1550,800],'Name','Reconstruction with 1-axis MDEIT');

% MDEIT reconstruction
for i = 1:length(SNR_vector)
    subplot(4,length(SNR_vector),i)
    box on
    show_fem_transparent_edges(all_data{i}.img_output_mdeit_1)

    title_str = sprintf('$1$-MDEIT, SNR = %.1f, $l = %.1g$', ...
        SNR_vector(i), ...
        all_data{i}.imdl_mdeit_1.hyperparameter.value);

    % title_str = strcat(...
    %     ['$1$-MDEIT, ' ...
    %     'SNR = '],num2str(SNR_vector(i)), ...
    %     ', $\lambda = $',num2str(all_data{i}.imdl_mdeit_1.hyperparameter.value));
    title(title_str,'Interpreter','latex')
end

% MDEIT noise correction
for i = 1:length(SNR_vector)
    subplot(4,length(SNR_vector),1*length(SNR_vector)+i)
    box on
    img = all_data{i}.img_output_mdeit_1;
    img.elem_data = img.elem_data./all_data{i}.sigma_std_deviation_mdeit;
    show_fem_transparent_edges(img)
    % title_str = strcat('$1$-MDEIT, SNR = ',num2str(SNR_vector(i)));
    % title(title_str,'Interpreter','latex')
end

% EIT reconstruction
subplot(4,length(SNR_vector),1)
for i = 1:length(SNR_vector)
    subplot(4,length(SNR_vector),2*length(SNR_vector)+i)
    box on
    show_fem_transparent_edges(all_data{i}.img_output_eit)

    title_str = sprintf('EIT, SNR = %.1f, $l = %.1g$', ...
        SNR_vector(i), ...
        all_data{i}.imdl_eit.hyperparameter.value);

    % title_str = strcat('EIT, SNR = ',num2str(SNR_vector(i)));
    title(title_str,'Interpreter','latex')
end

% EIT noise correction
subplot(4,length(SNR_vector),1)
for i = 1:length(SNR_vector)
    subplot(4,length(SNR_vector),3*length(SNR_vector)+i)
    box on

    img = all_data{i}.img_output_eit;
    img.elem_data = img.elem_data./all_data{i}.sigma_std_deviation_eit;

    show_fem_transparent_edges(img)
    % title_str = strcat('EIT, SNR = ',num2str(SNR_vector(i)));
    % title(title_str,'Interpreter','latex')
end

%%
% cla
% show_fem(img)
% 
% electrodes = img.fwd_model.electrode;
% nodes = img.fwd_model.nodes;
% 
% elec_ids = 1:numel(electrodes);
% 
% % Given:
% w = 0.1;           % width
% h = 0.3;             % height
% 
% minX = min(nodes(:,1))-h;
% maxX = max(nodes(:,1))+h;
% 
% minY = min(nodes(:,2))-h;
% maxY = max(nodes(:,2))+h;
% 
% 
% hold on
% for i = 1:numel(elec_ids)
%     elec_nodes = electrodes(elec_ids(i)).nodes;
% 
%     elec_center = mean(nodes(elec_nodes,:));
% 
%     c = elec_center;
% 
%     % Normalize orientation vector
%     n = elec_center / norm(elec_center);
% 
%     % Compute perpendicular vector
%     t = [-n(2), n(1)];   % rotate n by +90 degrees      
% 
%     h = sin(elec_center(1))/5;
% 
%     % Half dimensions
%     w2 = w/2;
%     h2 = h;
% 
%     % Compute rectangle corners relative to bottom center p
%     % bottom left
%     bl = c - w2*t;
%     % bottom right
%     br = c + w2*t;
%     % top right
%     tr = br + h*n;
%     % top left
%     tl = bl + h*n;
% 
%     % Coordinates
%     X = [bl(1) br(1) tr(1) tl(1)];
%     Y = [bl(2) br(2) tr(2) tl(2)];
% 
%     if h>0
%         color = 'red';
%     else
%         color = 'blue';
%     end
%     % Draw the rectangle
%     fill(X, Y, color, 'FaceAlpha', 1);  % red filled
% 
%     for j = 1:length(elec_nodes)-1
%         node1 = nodes(elec_nodes(j),:);
%         node2 = nodes(elec_nodes(j+1),:);
% 
%         % plot(node1(1),node1(2),'r.','MarkerSize',10)
%         plot([node1(1) node2(1)],[node1(2),node2(2)],'g','LineWidth',3);
%     end
% end
% hold off
% axis([minX maxX minY maxY])
% axis off
%% FUNCTIONS
function show_fem_transparent_edges(img)

hh = show_fem(img);                % draw the model (hh may be a handle or array)
% find the patch objects that actually draw the elements and remove their edges
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    % sometimes hh is an axes handle or figure; search the axes too:
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);

end

