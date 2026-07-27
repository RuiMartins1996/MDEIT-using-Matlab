function Gamma_noise = calibrate_noise_covariance(ctx, theta_init, cfg)
%CALIBRATE_NOISE_COVARIANCE  Gamma_noise = (1e-3 * max_{k,l}|U_k-U_l|)^2 * I
%at the INITIAL electrode configuration (paper sec 5), computed once and
%held fixed thereafter (PLAN_implementation.md sec 5 -- do not recompute
%this inside the objective, that would change what is being optimized).
%
%   Gamma_noise = calibrate_noise_covariance(ctx, theta_init, cfg)

M = numel(theta_init);
elec = assemble_electrodes(ctx, theta_init, cfg.width, cfg.z_contact);
sigma = cfg.sigma_star*ones(ctx.n_elem,1);
I = build_current_patterns(M);
fwd = cem_fwd_solve(ctx, sigma, elec, I);

Uall = fwd.V(:);   % mean-free electrode potentials for all patterns
spread = max(Uall) - min(Uall);
noise_std = 1e-3*spread;

n_data = M*(M-1);
Gamma_noise = (noise_std^2)*speye(n_data);
end
