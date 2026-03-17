clc; clear all; close all;


%%

% This script takes the figure and outputs another figure with only the
% second row of data in a different view, corresponding to the
% noise-corrected reconstruction

name = 'snr_10';
ext  = '.fig';

fullpath = fullfile(pwd, [name ext]);

uiopen(fullpath,1);


for i = 1:6
    ax1 = subplot(4,6,i);
    ax2 = subplot(4,6,6+i);
    view(ax2,180,0)
    ax3 = subplot(4,6,12+i);
    ax4 = subplot(4,6,18+i);

    delete(ax1);
    delete(ax3);
    delete(ax4);
end

% Get all axes in figure
axAll = findall(gcf,'Type','axes');

% Remove legends and colorbars
axAll = axAll(~arrayfun(@(a) ...
    isa(a,'matlab.graphics.illustration.Legend') || ...
    isa(a,'matlab.graphics.illustration.ColorBar'), axAll));

% Keep only valid handles
axAll = axAll(isvalid(axAll));

% Sort left-to-right based on X position
[~, idx1] = sort(arrayfun(@(a) a.Position(1), axAll));
axAll = axAll(idx1);

newFig = figure('Position',[200 200 1000 500]);


% Recreate as 1x6
for i = 1:numel(axAll)
    newAx = subplot(1,6,i);
    
    colormap(jet)     % or jet, turbo, hot, etc.

    copyobj(axAll(i).Children, newAx);
    % Copy axis properties
    newAx.XLim = axAll(i).XLim;
    newAx.YLim = axAll(i).YLim;
    newAx.ZLim = axAll(i).ZLim;
    newAx.View = axAll(i).View;
    newAx.Title.String = axAll(i).Title.String;

    box on;
    axis([-1 1 -1 1 0 3])

    daspect([1 1 1])      % <-- enforce equal unit scaling
    pbaspect([1 1 1])     % optional: make box visually cubic

    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
end

set(gcf,'Position',[200 200 1000 300])

%%

clc; clear all; close all;
% This script takes the figure and outputs another figure with the first 
% two rows of data in a different view, corresponding to the reconstructed
% conductivity and the noise-corrected reconstruction

name = 'snr_20';
ext  = '.fig';

fullpath = fullfile(pwd, [name ext]);

uiopen(fullpath,1);


for i = 1:6
    ax1 = subplot(4,6,i);
    view(ax1,180,0)
    ax2 = subplot(4,6,6+i);
    view(ax2,180,0)
    ax3 = subplot(4,6,12+i);
    ax4 = subplot(4,6,18+i);

    delete(ax3);
    delete(ax4);
end

% Get all axes in figure
axAll = findall(gcf,'Type','axes');

% Remove legends and colorbars
axAll = axAll(~arrayfun(@(a) ...
    isa(a,'matlab.graphics.illustration.Legend') || ...
    isa(a,'matlab.graphics.illustration.ColorBar'), axAll));

% Keep only valid handles
axAll = axAll(isvalid(axAll));

axAll_row_1 = axAll(1:6);
axAll_row_2 = axAll(7:12);

% Sort left-to-right based on X position
[~, idx1] = sort(arrayfun(@(a) a.Position(1), axAll_row_1));
axAll_row_1 = axAll_row_1(idx1);

[~, idx2] = sort(arrayfun(@(a) a.Position(1), axAll_row_2));
axAll_row_2 = axAll_row_2(idx2);

newFig = figure('Position',[200 200 1000 500]);


% Recreate as 1x6
for i = 1:numel(axAll)
    newAx = subplot(2,6,i);
    
    colormap(jet)     % or jet, turbo, hot, etc.

    if i<7
        ax = axAll_row_1(i);
    else
        ax = axAll_row_2(i-6);
    end

    copyobj(ax.Children, newAx);
    % Copy axis properties
    newAx.XLim = ax.XLim;
    newAx.YLim = ax.YLim;
    newAx.ZLim = ax.ZLim;
    newAx.View = ax.View;
    newAx.Title.String = ax.Title.String;

    box on;
    axis([-1 1 -1 1 0 3])

    daspect([1 1 1])      % <-- enforce equal unit scaling
    pbaspect([1 1 1])     % optional: make box visually cubic

    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
end

set(gcf,'Position',[200 200 1000 300])


%%

clc; clear all; close all;
% This script takes the figure and outputs another figure with the second
% and the fourth rows of datain a different view, corresponding to the
% noise-corrected reconstruction and to the horizontal slice of that
% reconstruction


name = 'snr_1_side_anomaly';
ext  = '.fig';

fullpath = fullfile(pwd, [name ext]);

uiopen(fullpath,1);

% Determine Layout
fig = gcf;
ax  = findall(fig,'Type','axes','-not','Tag','Colorbar');

pos = vertcat(ax.Position);   % [x y width height]

% Estimate grid structure
xCenters = pos(:,1) + pos(:,3)/2;
yCenters = pos(:,2) + pos(:,4)/2;

cols = numel(unique(round(xCenters,2)));
rows = numel(unique(round(yCenters,2)));

for i = 1:cols
    ax1 = subplot(rows,cols,i);
    
    ax2 = subplot(rows,cols,cols+i);
    view(ax2,180,0)

    ax3 = subplot(rows,cols,2*cols+i);
    
    ax4 = subplot(rows,cols,3*cols+i);

    delete(ax1);
    delete(ax3);
end

% Get all axes in figure
axAll = findall(gcf,'Type','axes');

% Remove legends and colorbars
axAll = axAll(~arrayfun(@(a) ...
    isa(a,'matlab.graphics.illustration.Legend') || ...
    isa(a,'matlab.graphics.illustration.ColorBar'), axAll));

% Keep only valid handles
axAll = axAll(isvalid(axAll));

axAll_row_1 = axAll(1:cols);
axAll_row_2 = axAll(cols+1:2*cols);

% Sort left-to-right based on X position
[~, idx1] = sort(arrayfun(@(a) a.Position(1), axAll_row_1));
axAll_row_1 = axAll_row_1(idx1);

[~, idx2] = sort(arrayfun(@(a) a.Position(1), axAll_row_2));
axAll_row_2 = axAll_row_2(idx2);

newFig = figure('Position',[200 200 1300 500]);


% Recreate as 1xcols
for i = 1:numel(axAll)
    newAx = subplot(2,cols,i);
    
    colormap('jet')     % or jet, turbo, hot, etc.

    if i<cols+1
        ax = axAll_row_1(i);
    else
        ax = axAll_row_2(i-cols);
    end

    copyobj(ax.Children, newAx);
    % Copy axis properties
    newAx.XLim = ax.XLim;
    newAx.YLim = ax.YLim;
    newAx.ZLim = ax.ZLim;
    newAx.View = ax.View;
    newAx.Title.String = ax.Title.String;
    newAx.Title.Interpreter = 'latex';

    box on;
    axis([-1 1 -1 1 0 3])

    daspect([1 1 1])      % <-- enforce equal unit scaling
    pbaspect([1 1 1])     % optional: make box visually cubic

    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
end

set(gcf,'Position',[200 200 1300 300])