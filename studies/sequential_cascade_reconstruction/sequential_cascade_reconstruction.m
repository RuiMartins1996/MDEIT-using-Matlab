clc; clear all; close all;
%SEQUENTIAL_CASCADE_RECONSTRUCTION
%
% Approach 3 of joint_eit_mdeit_strategies.pdf: sequential (cascaded)
% EIT + MDEIT reconstruction with cross-modality priors.
%
% Rather than stacking the EIT and MDEIT residuals into a single augmented
% least-squares problem (which the rank analysis shows adds no information
% and which re-imports EIT's vertical elongation), this study uses one
% modality to shape the *regularizer* of the other:
%
%   Stage 1  : single-modality Tikhonov reconstructions (EIT and 3-axis MDEIT)
%   Consensus: a soft spatial weight map W built from the AGREEMENT of both
%              images; the EIT image contributes only its z-max-projected
%              transverse support, MDEIT caps the vertical extent
%   Stage 2  : re-solve MDEIT with penalty lambda * ||W dsigma||^2, iterated
%              (IRLS) so it approximates an edge-preserving / sparsity prior
%
% The mechanism directly attacks the vertical-elongation artifact and
% degrades gracefully to MDEIT-only when EIT is uninformative.
%
% The heavy lifting lives in three reusable helpers in this folder:
%   run_consensus_cascade.m      - Algorithm 1 (the whole cascade)
%   build_consensus_weight_map.m - the consensus soft-weight map + z-projection
%   solve_weighted_tikhonov.m    - weighted Tikhonov solve with GCV lambda
%   cascade_metrics.m            - localization error + elongation proxy
%
% Reported per SNR: EIT-only, MDEIT-only, the cascade, the consensus map, a
% homogeneous null test, held-out EIT misfit, localization error and the
% vertical/transverse elongation ratio.

%% Prepare workspace
fullpath = mfilename('fullpath');
script_folder = fileparts(fullpath);
cd(script_folder);

grandparent_folder = fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));

model_folder = prepare_workspace(script_folder);
data_folder  = fullfile(script_folder,'data');
clc;
rng(1)

%% Characteristic scales (SI) -- same convention as the other studies
z0 = 0.0058;    % (Ohm m^2) contact impedance (CEM article, 58 Ohm cm^2)
l0 = 40e-3;     % (m) tank radius
I0 = 2.4e-3;    % (A) injected current magnitude

V0 = z0*I0/(l0^2);   % (V)   voltage scale
sigma0 = l0/z0;      % (S/m) conductivity scale
J0 = I0/(l0^2);      % current-density scale

background_conductivity = 3.28e-1/sigma0;

%% Experiment configuration (edit here)
SNR_list      = [100 20 1];      % SNRs to sweep (dB), per the experimental plan
anomaly_case  = 'side';          % 'side' (EIT localizes well) or 'center'
contrast      = 2.0;             % anomaly conductivity / background (conductive sphere)

% Mesh sizes (characteristic units).  Coarser = faster; refine for accuracy.
maxsz_data           = 0.20;     % data (forward) mesh
maxsz_reconstruction = 0.30;     % reconstruction (inverse) mesh

% Cascade hyperparameters (Algorithm 1 defaults)
lambda_vector = logspace(-12,3,30);
cascade_T       = 3;             % IRLS iterations
cascade_epsilon = 0.05;          % soft-weight floor
cascade_p       = 0.5;           % soft-weight exponent (p=1/2 -> l1-type prior)
cascade_rule    = 'geomean';     % consensus rule: 'geomean' or 'min'

%% Build model parameters (handoff geometry: R=40mm, H=120mm)
model_parameters = create_kai_3d_model_parameters(l0, z0, sigma0, I0);

model_parameters.radius = 40e-3/l0;    % = 1
model_parameters.height = 120e-3/l0;   % = 3
model_parameters.numOfRings              = 4;
model_parameters.numOfElectrodesPerRing  = 8;   % 4 rings x 8 = 32 electrodes
model_parameters.numOfSensors            = 32;  % 8 per ring x 4 rings
model_parameters.sensorRadius            = 140e-3/l0;  % external magnetometer ring
model_parameters.electrodeRadius         = 8e-3/l0;
model_parameters.electrodeContactImpedance = 0.0058/z0;
model_parameters.mu0 = 1;              % compute H in units of I0/l0
model_parameters.maxsz = maxsz_data;

anomaly_radius = 12.5e-3/l0;           % = 0.3125
switch lower(anomaly_case)
    case 'center'
        anomaly_position = [0, 0, 60e-3/l0];
    case 'side'
        anomaly_position = [20e-3/l0, 0, 60e-3/l0];
    otherwise
        error('anomaly_case must be ''center'' or ''side''');
end
anomaly_conductivity = contrast * background_conductivity;

% A material carves the anomaly geometry into the data mesh
model_parameters.material = struct( ...
    'type','spherical', 'name','sphere_anomaly', ...
    'radius',anomaly_radius, 'position',anomaly_position);
model_parameters.anomaly = struct( ...
    'type','spherical', 'conductivity',anomaly_conductivity, ...
    'radius',anomaly_radius, 'position',anomaly_position);

%% Build data (forward) and reconstruction (inverse) forward models
[~, fmdls_data] = mk_mdeit_model(model_parameters, model_folder);
fmdl_data = fmdls_data{1};

mp_rec = model_parameters;
mp_rec.material = struct();            % homogeneous reconstruction mesh
mp_rec.maxsz    = maxsz_reconstruction;
[~, fmdls_rec] = mk_mdeit_model(mp_rec, model_folder);
fmdl_rec = fmdls_rec{1};

%% Skip-2 stimulation pattern (2.4 mA), applied to both meshes
current_amplitude = 2.4e-3/I0;         % = 1
inj = [0 3]; meas = [0 3];             % skip-2
stim = mk_stim_patterns(numel(fmdl_data.electrode), 1, inj, meas, {}, current_amplitude);
fmdl_data.stimulation = stim;
fmdl_rec.stimulation  = stim;

n_elem = size(fmdl_rec.elems,1);

%% Ground-truth data (fine mesh): homogeneous vs anomaly
imgh = mk_image_mdeit(fmdl_data, background_conductivity);
imgi = add_material_properties(imgh, [background_conductivity, anomaly_conductivity]);

[datah, uh] = fwd_solve_mdeit(imgh);
[datai, ui] = fwd_solve_mdeit(imgi);

mdeit3vec = @(d)[d.Bx(:); d.By(:); d.Bz(:)];   % 3-axis stack, matches [Jx;Jy;Jz]
dB_real = mdeit3vec(datai) - mdeit3vec(datah); % MDEIT difference signal
du_real = ui.meas          - uh.meas;          % EIT   difference signal

%% Reconstruction-mesh Jacobians (linearized at the background)
% 3-axis MDEIT Jacobian: obtained through inv_solve_mdeit, which returns it.
imdl_M = eidors_obj('inv_model','cascade_mdeit3');
imdl_M.fwd_model      = fmdl_rec;
imdl_M.jacobian_bkgnd = struct('value',background_conductivity);
imdl_M.recon_mode     = 'mdeit3';
imdl_M.recon_type     = 'difference';
imdl_M.solver         = 'gn';
imdl_M.RtR_prior      = @(x,Jx) x;
imdl_M.hyperparameter = struct('value',1e-6);
imdl_M.verbose        = false;

fprintf('Assembling 3-axis MDEIT Jacobian ...\n');
tmpM = inv_solve_mdeit(imdl_M, mdeit3vec(datai), mdeit3vec(datah));
J_M  = tmpM.jacobian;

% EIT Jacobian via EIDORS
fprintf('Assembling EIT Jacobian ...\n');
img_eit_b = mk_image(fmdl_rec, background_conductivity);
img_eit_b.fwd_model.stimulation = stim;
J_E = calc_jacobian(img_eit_b);

%% Element centroids (for the z-max projection and metrics)
E = fmdl_rec.elems; N = fmdl_rec.nodes;
centroids = zeros(n_elem, size(N,2));
for k = 1:size(E,2)
    centroids = centroids + N(E(:,k),:);
end
centroids = centroids / size(E,2);
if size(centroids,2) == 2                 % 2D safety: pad a z-column
    centroids(:,3) = 0;
end

proj_cell = 0.5 * anomaly_radius;         % transverse cell for z-max projection

%% Sweep SNR
results = struct('SNRdb',{},'out',{},'out_null',{}, ...
                 'metrics_eit',{},'metrics_mdeit',{},'metrics_cascade',{}, ...
                 'metrics_null',{},'heldout_eit_misfit',{});

for si = 1:numel(SNR_list)
    SNRdb = SNR_list(si);
    fprintf('\n==================== SNR = %g dB ====================\n', SNRdb);

    rng(1)   % reproducible noise across SNRs

    % Noise levels (max-amplitude convention, matches the other studies)
    nlM = max(abs(dB_real - mean(dB_real))) / 10^(SNRdb/20);
    nlE = max(abs(du_real - mean(du_real))) / 10^(SNRdb/20);

    dyM = dB_real + nlM * randn(size(dB_real));
    dyE = du_real + nlE * randn(size(du_real));

    % --- Cascade options ---
    copts = struct( ...
        'sigma_M', nlM, 'sigma_E', nlE, ...
        'lambda_vector', lambda_vector, ...
        'T', cascade_T, 'epsilon', cascade_epsilon, 'p', cascade_p, ...
        'consensus', cascade_rule, 'use_projection', true, ...
        'proj_cell', proj_cell, 'verbose', true);

    % --- Run the cascade ---
    out = run_consensus_cascade(J_M, dyM, J_E, dyE, centroids, copts);

    % --- Null test: homogeneous tank (noise only, no anomaly) ---
    dyM0 = nlM * randn(size(dB_real));
    dyE0 = nlE * randn(size(du_real));
    copts_null = copts; copts_null.verbose = false;
    out_null = run_consensus_cascade(J_M, dyM0, J_E, dyE0, centroids, copts_null);

    % --- Metrics ---
    mE = cascade_metrics(out.dsE,     centroids, anomaly_position);
    mM = cascade_metrics(out.dsM,     centroids, anomaly_position);
    mC = cascade_metrics(out.dsigma,  centroids, anomaly_position);
    m0 = cascade_metrics(out_null.dsigma, centroids, anomaly_position);

    % Held-out modality check: EIT is NOT in the Stage-2 data term.
    heldout = norm(J_E*out.dsigma - dyE) / norm(dyE);

    fprintf('\n--- Metrics (SNR = %g dB) ---\n', SNRdb);
    fprintf('%-14s | loc.err  | elong (rmsZ/rmsXY)\n','method');
    fprintf('%-14s | %7.4f  | %6.3f\n','EIT-only',    mE.localization_error, mE.elongation);
    fprintf('%-14s | %7.4f  | %6.3f\n','MDEIT-only',  mM.localization_error, mM.elongation);
    fprintf('%-14s | %7.4f  | %6.3f\n','Cascade',     mC.localization_error, mC.elongation);
    fprintf('Held-out EIT misfit (rel.): %.3f   (null-test peak |dsigma|: %.2e)\n', ...
        heldout, max(abs(out_null.dsigma)));

    % --- Store ---
    results(si).SNRdb              = SNRdb;
    results(si).out                = out;
    results(si).out_null           = out_null;
    results(si).metrics_eit        = mE;
    results(si).metrics_mdeit      = mM;
    results(si).metrics_cascade    = mC;
    results(si).metrics_null       = m0;
    results(si).heldout_eit_misfit = heldout;

    % --- Figure: reconstructions for this SNR ---
    figure('Name',sprintf('Cascade reconstruction, SNR=%g',SNRdb), ...
           'Position',[100 80 1400 700]);

    subplot(2,3,1); show_on_mesh(imgi);                              title('Ground truth','Interpreter','latex');
    subplot(2,3,2); show_on_mesh(mk_dimg(fmdl_rec, out.dsE));        title('EIT-only (Stage 1a)','Interpreter','latex');
    subplot(2,3,3); show_on_mesh(mk_dimg(fmdl_rec, out.dsM));        title('MDEIT-only (Stage 1b)','Interpreter','latex');
    subplot(2,3,4); show_on_mesh(mk_dimg(fmdl_rec, out.s));          title('Consensus map $s$','Interpreter','latex');
    subplot(2,3,5); show_on_mesh(mk_dimg(fmdl_rec, out.dsigma));     title('Cascade (Stage 2)','Interpreter','latex');
    subplot(2,3,6); show_on_mesh(mk_dimg(fmdl_rec, out_null.dsigma));title('Null test (no anomaly)','Interpreter','latex');
    sgtitle(sprintf('Sequential EIT+MDEIT cascade, %s anomaly, SNR = %g dB', anomaly_case, SNRdb));
    pause(1e-3);
end

%% Summary metrics across SNR
loc_mdeit = arrayfun(@(r) r.metrics_mdeit.localization_error,   results);
loc_casc  = arrayfun(@(r) r.metrics_cascade.localization_error, results);
elo_mdeit = arrayfun(@(r) r.metrics_mdeit.elongation,           results);
elo_casc  = arrayfun(@(r) r.metrics_cascade.elongation,         results);
snrs      = [results.SNRdb];

figure('Name','Cascade vs MDEIT-only: metrics','Position',[200 200 1000 400]);
subplot(1,2,1); hold on; box on; grid on; grid minor;
plot(snrs, loc_mdeit, 'o-'); plot(snrs, loc_casc, 's-');
set(gca,'XScale','log'); xlabel('SNR (dB)','Interpreter','latex');
ylabel('localization error','Interpreter','latex');
legend('MDEIT-only','Cascade','Location','best'); title('Localization error','Interpreter','latex');

subplot(1,2,2); hold on; box on; grid on; grid minor;
plot(snrs, elo_mdeit, 'o-'); plot(snrs, elo_casc, 's-');
yline(1,'k--');
set(gca,'XScale','log'); xlabel('SNR (dB)','Interpreter','latex');
ylabel('elongation $r_{z}/r_{xy}$','Interpreter','latex');
legend('MDEIT-only','Cascade','Location','best'); title('Vertical elongation','Interpreter','latex');

%% Save
save_name = fullfile(data_folder, sprintf('cascade_%s_anomaly.mat', anomaly_case));
save(save_name, 'results','anomaly_case','anomaly_position','anomaly_radius', ...
     'contrast','background_conductivity','SNR_list','lambda_vector', ...
     'cascade_T','cascade_epsilon','cascade_p','cascade_rule','fmdl_rec','imgi','-v7.3');
fprintf('\nSaved results to %s\n', save_name);


%% ===================== local plotting helpers =====================
function img = mk_dimg(fmdl, v)
%MK_DIMG  Wrap an element-wise vector into an EIDORS image on fmdl.
img = mk_image(fmdl, 0);
img.elem_data = v(:);
end

function show_on_mesh(img)
%SHOW_ON_MESH  show_fem with faint element edges and a colourbar.
hh = show_fem(img);
patches = findobj(hh, 'Type', 'Patch');
if isempty(patches)
    patches = findobj(gca, 'Type', 'Patch');
end
set(patches, 'EdgeAlpha', 0.1);
try, eidors_colourbar(img); catch, end
box on;
end
