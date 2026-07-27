function [post_var, phi_a, phi_d] = posterior_variance_diag(ctx, theta_minus, cfg)
%POSTERIOR_VARIANCE_DIAG  diag(Gamma_*(theta)) via the same Woodbury
%data-space objects as oed_criterion.m (never forms the n_elem x n_elem
%posterior covariance). Also returns both phi_a=trace(Gamma_*) and
%phi_d=logdet(Gamma_*) in one pass (P,Y,S,Ls are mode-independent).
%
%   [post_var, phi_a, phi_d] = posterior_variance_diag(ctx, theta_minus, cfg)

M = numel(theta_minus);
elec = assemble_electrodes(ctx, theta_minus, cfg.width, cfg.z_contact);
sigma = cfg.sigma_star*ones(ctx.n_elem,1);
I = build_current_patterns(M);
fwd = cem_fwd_solve(ctx, sigma, elec, I);
J = cem_jacobian(ctx, fwd);

[phi_a, ~, extra] = oed_criterion(J, cfg.Gamma_prior, cfg.Gamma_noise, 'a-opt');
post_var = diag(cfg.Gamma_prior) - sum(extra.P.*extra.Y,1).';

if nargout > 2
    logdet_S = 2*sum(log(diag(extra.Ls)));
    n_elem = ctx.n_elem;
    jitter_p = 1e-12*trace(cfg.Gamma_prior)/n_elem;
    Lp = chol(full(cfg.Gamma_prior) + jitter_p*eye(n_elem),'lower');
    logdet_Gpr = 2*sum(log(diag(Lp)));
    logdet_Gn = sum(log(diag(cfg.Gamma_noise)));
    phi_d = logdet_Gpr + logdet_Gn - logdet_S;
end
end
