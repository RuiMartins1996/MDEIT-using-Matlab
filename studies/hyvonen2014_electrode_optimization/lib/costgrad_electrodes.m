function [phi, grad, diagnostics] = costgrad_electrodes(ctx, theta_minus, cfg, opt_mode)
%COSTGRAD_ELECTRODES  Cost psi(theta_minus) and analytic gradient, ready
%for fminunc('SpecifyObjectiveGradient',true) or Algorithm 1
%(PLAN_implementation.md sec 4, 6). This is the single entry point tying
%together the CEM forward/adjoint solves, the Jacobian and its
%theta-derivative, the Woodbury OED criterion, and the repulsion term.
%
%   [phi, grad, diagnostics] = costgrad_electrodes(ctx, theta_minus, cfg, opt_mode)
%
% cfg fields (see lib/hyvonen_config.m):
%   width, z_contact, sigma_star, Gamma_prior, Gamma_noise, alpha_rep
% opt_mode: 'a-opt' or 'd-opt'
%
% grad is [M x 1], d(psi)/d(theta_minus).

M = numel(theta_minus);
n_elem = ctx.n_elem;

elec = assemble_electrodes(ctx, theta_minus, cfg.width, cfg.z_contact);

sigma = cfg.sigma_star*ones(n_elem,1);
I = build_current_patterns(M);

fwd = cem_fwd_solve(ctx, sigma, elec, I);
J = cem_jacobian(ctx, fwd);

[phi_data, W] = oed_criterion(J, cfg.Gamma_prior, cfg.Gamma_noise, opt_mode);

dJ = cem_jacobian_dtheta(ctx, fwd, elec);
grad_data = zeros(M,1);
for m = 1:M
    grad_data(m) = -2*sum(sum(W.*dJ{m}));
end

[phi_rep, grad_rep] = repulsion_term(elec, cfg.alpha_rep);

phi  = phi_data + phi_rep;
grad = grad_data + grad_rep;

if nargout > 2
    diagnostics = struct('phi_data',phi_data,'phi_rep',phi_rep, ...
        'J',J,'W',W,'elec',elec,'fwd',fwd);
end
end
