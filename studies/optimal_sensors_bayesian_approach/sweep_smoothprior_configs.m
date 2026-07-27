%% Sweep driver: how large can the MDEIT A-optimality reduction be made?
%
% Companion to example_anomaly_circle_2d_smoothprior.m. Answers "is a ~50%
% A-opt cost reduction (like the EIT/Hyvonen case) reachable for MDEIT,
% and under what conditions?"
%
% IMPORTANT: each config is run in its OWN MATLAB process. The MDEIT
% script shares a workspace with whatever calls it and reuses common
% variable names (e.g. it assigns 'c' at the reconstruction step), so
% looping over configs inside one workspace silently corrupts the driver's
% own loop state. Each child process writes data/sweep_<tag>.mat; this
% driver only aggregates.
%
% Levers swept (all motivated by the ceiling analysis in
% studies/hyvonen2014_electrode_optimization/analysis_aopt_reduction_ceiling.m):
%   lambda_fac   - prior correlation length / ROI radius. Larger => lower
%                  effective prior rank => higher information ceiling.
%   d_target     - how many data modes sit above the noise floor.
%   n_sensors    - more sensors = more data modes AND more freedom to
%                  cluster them near the anomaly.
%   anom_offset  - anomaly distance from centre (x tank radius). Larger =>
%                  steeper Biot-Savart contrast between a sensor beside
%                  the anomaly and one on the far side, which is what lets
%                  an OPTIMIZED array beat an EVEN one.
%   n_electrodes - fewer current patterns puts more of the information
%                  burden on sensor placement.
%
% Usage:  matlab -batch "sweep_smoothprior_configs"
% (spawns child MATLAB processes; expect a few minutes per config)

script_folder = fileparts(mfilename('fullpath'));
cd(script_folder);

matlab_exe = fullfile(matlabroot,'bin','matlab.exe');

% tag, n_elec, n_sens, d_target, anom_offset, lambda_fac, bg_std
%
% Round-1 findings (kept for the record):
%   - MORE SENSORS HURTS this metric: 4->8->12 sensors gave 32.7->24.0->22.6%
%     even though the ceiling rose 87.8->95.1->97.7%. Reason: the metric is
%     optimized-vs-EVEN, and a dense even ring already has a sensor near the
%     anomaly, so there is little left for placement to win. Same effect the
%     study's own handoff noted for electrodes ("fewer electrodes helps").
%   - Smoother prior helps a lot: lambda 1x->2x ROI radius gave 20.4->32.7%.
%   - Moving the anomaly nearer the wall alone did not help (20.4->18.8%).
% Round 2 therefore pushes: few sensors, low-rank prior, and a SMALLER
% background prior std (the background trace is an irreducible floor in phi
% that dilutes the percentage -- see the note in the summary).
%   - Round 2: lambda 2x->4x ROI radius gave 32.7->39.0%. Lowering the
%     background prior std HURT (39.0->27.4%) -- the background trace is
%     not the binding constraint. 3 sensors gave 28.1%, i.e. fewer than 4
%     is also worse: 4 sensors is near a sweet spot for this geometry.
% Round 3 therefore pushes only the one lever that is monotonically
% helping (prior correlation length) plus the noise level d_target.
cfgs = {
    'r3_lam8',      2, 4, 4, 0.75,  8.0, 0.03
    'r3_lam16',     2, 4, 4, 0.75, 16.0, 0.03
    'r3_lam8_d2',   2, 4, 2, 0.75,  8.0, 0.03
    'r3_lam8_d1',   2, 4, 1, 0.75,  8.0, 0.03
    'r3_lam16_d2',  2, 4, 2, 0.75, 16.0, 0.03
    };

run_children = true;   % set false to only re-aggregate existing .mat files

if run_children
    for k = 1:size(cfgs,1)
        tag = cfgs{k,1};
        fprintf('\n########## CONFIG %i/%i: %s ##########\n', k, size(cfgs,1), tag);
        cmd = sprintf(['"%s" -batch "cd(''%s''); sw_tag=''%s''; sw_n_electrodes=%d; ' ...
            'sw_n_sensors=%d; sw_d_target=%d; sw_anom_offset=%g; sw_lambda_fac=%g; sw_bg_std=%g; ' ...
            'example_anomaly_circle_2d_smoothprior"'], ...
            matlab_exe, script_folder, tag, cfgs{k,2}, cfgs{k,3}, cfgs{k,4}, cfgs{k,5}, cfgs{k,6}, cfgs{k,7});
        logf = fullfile(script_folder,sprintf('sweep_%s.log',tag));
        status = system([cmd ' > "' logf '" 2>&1']);
        if status ~= 0
            fprintf('  config "%s" FAILED (see %s)\n', tag, logf);
        end
    end
end

%% Aggregate
fprintf('\n\n============================== SWEEP SUMMARY ==============================\n');
fprintf('%-12s %6s %6s %6s %8s %9s %10s %11s %10s\n', ...
    'config','n_ele','n_sen','d_tgt','effrank','ceiling%','phi_even','phi_opt','reduction%');
for k = 1:size(cfgs,1)
    tag = cfgs{k,1};
    f = fullfile(script_folder,'data',['sweep_' tag '.mat']);
    if ~exist(f,'file')
        fprintf('%-12s  (no result)\n', tag); continue
    end
    S = load(f); r = S.sweep_row;
    fprintf('%-12s %6d %6d %6d %8.2f %8.1f%% %10.4e %11.4e %9.1f%%\n', ...
        r.tag, r.n_elec, r.n_sens, r.d_target, r.eff_rank, r.ceiling, ...
        r.phi_even, r.phi_opt, r.reduction);
end
fprintf('===========================================================================\n');
