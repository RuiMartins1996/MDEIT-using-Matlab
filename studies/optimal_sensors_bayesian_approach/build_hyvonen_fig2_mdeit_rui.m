%% Figure builder for the MDEIT analogue of Hyvonen et al. (2014) Fig. 5.1
%
% Loads data/hyvonen_fig2_mdeit_results.mat (produced by
% hyvonen_fig2_mdeit.m) and draws the 2x2 panel figure, positionally
% identical to the paper:
%
%   [ prior variance      ]     [ brute force, tr(Gamma)     ] 
%   [ even theta + post var]    [ fminunc, tr(Gamma)     ]

% Split from the compute script so the figure can be restyled without a
% rerun of the (expensive) brute-force sweep.

fullpath = mfilename('fullpath');
script_folder = fileparts(fullpath);
cd(script_folder);

grandparent_folder = fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));

data_path = fullfile(script_folder,'data','hyvonen_fig2_mdeit_results.mat');
assert(isfile(data_path),'Run hyvonen_fig2_mdeit first: %s not found',data_path);
S = load(data_path,'results');
R = S.results;

nodes = R.fmdl_recon.nodes;
elems = R.fmdl_recon.elems;
electrode_nodes = R.fmdl_recon.electrode_nodes;
rs = R.rs;
model_radius = R.model_parameters.radius;

clim_max = max(R.prior_variance);

figures_folder = fullfile(script_folder,'figures');
if ~exist(figures_folder,'dir'), mkdir(figures_folder); end

fig = figure('Name','Hyvonen Fig. 2 (MDEIT)','Position',[100 75 1200 650]);
t = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

pct = @(phi_init,phi) 100*(phi_init-phi)/phi_init;
nats = @(phi_init,phi) (phi_init-phi)/2;

%% Tile 1 (top-left): prior variance, electrodes only
ax1 = nexttile(t);
draw_panel(ax1,nodes,elems,electrode_nodes,R.prior_variance,clim_max,[],rs,model_radius);
title(ax1,'Prior variance (5.2)');

%% Tile 2: brute-force optimum, A-opt (tr(Gamma))
% Titles use the PURE (alpha_rep=0) tr(Gamma) value, not the
% repulsion-penalized objective fminunc/brute-force actually minimized
% (results.phi_bf_polished etc.) -- see hyvonen_fig2_mdeit.m header.
ax2 = nexttile(t);
draw_panel(ax2,nodes,elems,electrode_nodes,R.post_var_bf_polished{1}{1},clim_max,...
    R.theta_bf_polished{1}{1},rs,model_radius);
title(ax2,sprintf('Brute force, tr(\\Gamma): \\phi=%.3e\n(%.1f%% reduction)',...
    R.phi_pure_bf_polished(1),pct(R.phi_pure_even(1),R.phi_pure_bf_polished(1))));

%% Tile 3 (bottom-left): initial guess theta_even + its posterior variance
ax4 = nexttile(t);
draw_panel(ax4,nodes,elems,electrode_nodes,R.post_var_even,clim_max,R.theta_even,rs,model_radius);
title(ax4,sprintf('Even spacing: \\phi_A=%.3e',R.phi_pure_even(1)));

%% Tile 4: fminunc optimum, A-opt (tr(Gamma))
ax3 = nexttile(t);
draw_panel(ax3,nodes,elems,electrode_nodes,R.post_var_fminunc{1}{1},clim_max,...
    R.theta_fminunc{1}{1},rs,model_radius);
title(ax3,sprintf('fminunc, tr(\\Gamma): \\phi=%.3e\n(%.1f%% reduction)',...
    R.phi_pure_fminunc(1),pct(R.phi_pure_even(1),R.phi_pure_fminunc(1))));

%% Colormap and save
colormap(fig,hot);
cb = colorbar(ax4);
cb.Layout.Tile = 'east';
for ax = [ax1 ax2 ax3 ax4]
    clim(ax,[0 clim_max]);
end

sgtitle(t,sprintf(['MDEIT analogue of Hyvonen, Seppanen and Staboulis (2014) Fig. 5.1 -- ' ...
    '4 electrodes (fixed), 4 magnetometers (optimized azimuths)\n' ...
    'Magnetometer markers are an in-plane projection (sensors sit at height z_s above the z=0 current sheet)']),...
    'Interpreter','none');

savefig(fig,fullfile(figures_folder,'hyvonen_fig2_mdeit.fig'));
saveas(fig,fullfile(figures_folder,'hyvonen_fig2_mdeit.png'));
fprintf('Saved figures/hyvonen_fig2_mdeit.{fig,png}\n');

%% Optional: cost landscape (sorted phi_grid with fminunc / BF-polished marked)
figL = figure('Name','Cost landscape','Position',[100 100 1000 450]);
opt_modes = R.opt_modes{1};
for iom = 1:numel(opt_modes)
    subplot(1,numel(opt_modes),iom);
    phis = sort(R.phi_grid(:,iom));
    plot(phis,'k-','LineWidth',1.2); hold on;
    yline(R.phi_fminunc(iom),'r--','LineWidth',1.5,'DisplayName','fminunc');
    yline(R.phi_bf_polished(iom),'b--','LineWidth',1.5,'DisplayName','BF polished');
    xlabel('grid combination (sorted)'); ylabel('\phi');
    title(opt_modes{iom});
    legend('grid \phi','fminunc','BF polished','Location','best');
    grid on;
end
sgtitle(figL,'Cost landscape: sorted brute-force grid values vs. optimized points');
saveas(figL,fullfile(figures_folder,'hyvonen_fig2_mdeit_landscape.png'));
fprintf('Saved figures/hyvonen_fig2_mdeit_landscape.png\n');

%% ======================== LOCAL FUNCTIONS ========================

function draw_panel(ax,nodes,elems,electrode_nodes,post_var,clim_max,theta,rs,model_radius)
axes(ax); %#ok<LAXES>
patch(ax,'Faces',elems,'Vertices',nodes,'FaceVertexCData',post_var,...
    'FaceColor','flat','EdgeColor','none');
axis(ax,'equal','off');
clim(ax,[0 clim_max]);
hold(ax,'on');

draw_electrodes(ax,nodes,electrode_nodes,model_radius);

if ~isempty(theta)
    draw_sensors(ax,theta,rs);
end

lim = 1.25*rs;
xlim(ax,[-lim lim]); ylim(ax,[-lim lim]);
hold(ax,'off');
end

%% Electrodes: thick black arcs on the boundary, identical in every panel.
% Electrode 1 (stimulation(1)'s injecting electrode under
% mk_stim_patterns(...,'adjacent')) gets a longer radial "cord" marking it
% as the current-feeding electrode, as in the paper's figure.
function draw_electrodes(ax,nodes,electrode_nodes,model_radius)
for e = 1:numel(electrode_nodes)
    en = electrode_nodes{e};
    xy = nodes(en,:);
    ang = atan2(xy(:,2),xy(:,1));
    [~,ord] = sort(ang);
    xy = xy(ord,:);
    plot(ax,xy(:,1),xy(:,2),'k-','LineWidth',4);

    if e == 1
        mean_ang = atan2(mean(xy(:,2)),mean(xy(:,1)));
        cord = [0.15*model_radius*cos(mean_ang), 0.15*model_radius*sin(mean_ang)];
        outer = mean(xy,1) + cord;
        plot(ax,[mean(xy(:,1)) outer(1)],[mean(xy(:,2)) outer(2)],'k-','LineWidth',4);
    end
end
end

%% Magnetometers: markers at the in-plane projection (rs*cos, rs*sin), per-panel theta.
function draw_sensors(ax,theta,rs)
theta = theta(:);
xs = rs*cos(theta);
ys = rs*sin(theta);
plot(ax,xs,ys,'s','MarkerSize',10,'MarkerFaceColor',[0.2 0.4 0.9],...
    'MarkerEdgeColor','k','LineWidth',1.2);
end
