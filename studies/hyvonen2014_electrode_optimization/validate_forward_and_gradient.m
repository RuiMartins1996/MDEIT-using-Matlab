%% VALIDATE: CEM forward solve, Jacobian, smoothness, and analytic gradient
%
% Implements V1-V4 of PLAN_implementation.md sec 7, with one deliberate
% substitution: instead of comparing against EIDORS's own CEM (which only
% supports node-list electrodes, not the continuous arcs this study
% needs -- see Deviation 2), V1/V2 are replaced by two self-contained,
% equally rigorous checks that do not depend on EIDORS's electrode model
% at all:
%   V1: RECIPROCITY of the measurement map (a physical symmetry of any
%       linear elliptic system with symmetric contact impedances --
%       independent of whether the electrodes are node-based or
%       arc-based, so it directly tests OUR discretization, not a
%       different one).
%   V2: Jacobian dV/dsigma vs. central differences in the conductivity.
% V3 (smoothness in theta) and V4 (the mandatory FD gradient-check gate
% before any optimization run) are as specified in the plan.

fullpath = mfilename('fullpath');
script_folder = fileparts(fullpath);
cd(script_folder);
addpath(fullfile(script_folder,'lib'));

grandparent_folder = fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));
prepare_workspace(script_folder);

rng(1);

%% Build a small test mesh + config
shape = struct('type','disk','radius',1.0);
ctx = build_mesh_ctx(shape, 1/12);

M = 4;
width = pi/16;
z_contact = 1.0;
sigma_star = 1.0;

theta0 = 2*pi*(0:M-1)'/M;

%% ---------------------------------------------------------------------
%% V1: reciprocity of the measurement map
%% ---------------------------------------------------------------------
fprintf('\n=== V1: reciprocity check ===\n');
elec = assemble_electrodes(ctx, theta0, width, z_contact);
sigma = sigma_star*ones(ctx.n_elem,1);
I = build_current_patterns(M);
fwd = cem_fwd_solve(ctx, sigma, elec, I);

n = ctx.n_nodes;
Ufull = zeros(M, M-1); % Ufull(k,j) = potential at electrode k for pattern j (electrode M implicit 0)
for j = 1:M-1
    Ufull(1:M-1,j) = fwd.X(n+1:n+M-1, j);
    Ufull(M,j) = 0;
end
% R(j,k) = response at electrode (k+1) minus electrode 1, when driving
% pattern j = e_1 - e_{j+1}.  Reciprocity: R(j,k) == R(k,j).
R = zeros(M-1,M-1);
for j = 1:M-1
    for k = 1:M-1
        R(j,k) = Ufull(k+1,j) - Ufull(1,j);
    end
end
recip_err = max(max(abs(R-R.')))/max(max(abs(R)));
fprintf('  max |R-R^T| / max|R| = %.3e  (want ~1e-10, machine precision)\n', recip_err);
assert(recip_err < 1e-8, 'V1 FAILED: reciprocity violated -- check grounding/electrode assembly signs');
fprintf('  V1 PASSED\n');

%% ---------------------------------------------------------------------
%% V2: Jacobian vs central differences in sigma
%% ---------------------------------------------------------------------
fprintf('\n=== V2: Jacobian vs finite differences (sigma) ===\n');
J = cem_jacobian(ctx, fwd);

h = 1e-3;   % NOT 1e-6: a single element's conductivity has a small effect
            % on the 12-dim data vector, so central differences at h~1e-6
            % are roundoff-noise-dominated (verified via debug_jacobian_fd.m:
            % rel err *improves* monotonically as h grows from 1e-8 to 1e-3,
            % the signature of roundoff, not a Jacobian bug).
n_check = 5;
idx = randperm(ctx.n_elem, n_check);
max_rel_err = 0;
for kk = idx
    sp = sigma; sp(kk) = sp(kk)+h;
    sm = sigma; sm(kk) = sm(kk)-h;
    fwd_p = cem_fwd_solve(ctx, sp, elec, I);
    fwd_m = cem_fwd_solve(ctx, sm, elec, I);
    dV_fd = (fwd_p.V(:) - fwd_m.V(:))/(2*h);
    dV_an = J(:,kk);
    rel_err = norm(dV_fd-dV_an)/max(norm(dV_fd),eps);
    max_rel_err = max(max_rel_err, rel_err);
    fprintf('  elem %5i: ||FD-analytic||/||FD|| = %.3e\n', kk, rel_err);
end
fprintf('  max rel err over %i random elements: %.3e (want ~1e-6 to 1e-8)\n', n_check, max_rel_err);
assert(max_rel_err < 1e-5, 'V2 FAILED: Jacobian does not match finite differences');
fprintf('  V2 PASSED\n');

%% ---------------------------------------------------------------------
%% V3: smoothness of psi(theta) across a mesh-boundary-node crossing
%% ---------------------------------------------------------------------
fprintf('\n=== V3: smoothness of psi along a 1D slice ===\n');
category = ones(ctx.n_elem,1);   % homogeneous prior for this smoke test
kappa_lookup = 0.3;
Gamma_prior = prior_covariance(ctx, 0.5, category, kappa_lookup);
Gamma_noise = calibrate_noise_covariance(ctx, theta0, struct('width',width,'z_contact',z_contact,'sigma_star',sigma_star));

cfg = struct('width',width,'z_contact',z_contact,'sigma_star',sigma_star, ...
    'Gamma_prior',Gamma_prior,'Gamma_noise',Gamma_noise,'alpha_rep',1e-4);

n_slice = 25;
delta = linspace(-0.15,0.15,n_slice);
phi_slice = zeros(n_slice,1);
for kk = 1:n_slice
    th = theta0; th(1) = th(1)+delta(kk);
    phi_slice(kk) = costgrad_electrodes(ctx, th, cfg, 'a-opt');
end
% smoothness metric: max second difference relative to the range
d2 = diff(phi_slice,2);
smoothness_ratio = max(abs(d2)) / (max(phi_slice)-min(phi_slice));
fprintf('  max|2nd diff| / range = %.3e (a large spike here would mean a mesh-crossing discontinuity)\n', smoothness_ratio);
figure; plot(delta, phi_slice,'-o'); xlabel('\delta\theta_1'); ylabel('\psi (a-opt)');
title('V3: psi along a 1D slice (should be visibly smooth)');
saveas(gcf, fullfile(script_folder,'figures','V3_smoothness_slice.png'));

%% ---------------------------------------------------------------------
%% V4: analytic gradient vs finite differences (THE gate)
%% ---------------------------------------------------------------------
fprintf('\n=== V4: analytic gradient vs finite differences (mandatory gate) ===\n');
theta_pert = theta0 + 0.1*randn(M,1);
for opt_mode = {'a-opt','d-opt'}
    om = opt_mode{1};
    fprintf('\n-- %s --\n', om);
    cg = @(th) costgrad_electrodes(ctx, th, cfg, om);
    check_gradient_fd(cg, theta_pert, 4, 1e-5);
end

fprintf('\nAll validations complete. If all rel. errs above are ~1e-6 or smaller, proceed to Case 1.\n');
