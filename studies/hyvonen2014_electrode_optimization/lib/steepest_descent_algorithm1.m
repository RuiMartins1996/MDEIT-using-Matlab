function [theta_opt, phi_opt, info] = steepest_descent_algorithm1(costgrad, gaps_fn, theta0, opts)
%STEEPEST_DESCENT_ALGORITHM1  Hyvonen et al. 2014, Algorithm 1: normalized
%steepest descent with a bounded 1-D line search (PLAN_implementation.md
%sec 6 / paper sec 4, eq before (4.3) and after).
%
%   [theta_opt, phi_opt, info] = steepest_descent_algorithm1(costgrad, gaps_fn, theta0, opts)
%
% costgrad(theta) -> [phi, grad]           (grad = d(psi)/d(theta_minus))
% gaps_fn(theta)  -> [M x 1] electrode gaps g_m(theta) (cheap: geometry
%                    only, no FEM solve -- assemble_electrodes without a
%                    forward solve)
% opts fields (all optional, defaults shown):
%   max_iters      = 100
%   grad_tol       = 1e-6      stop if |grad| <= grad_tol
%   rel_tol        = 1e-8      stop if relative decrease in phi <= rel_tol
%   min_gap_frac   = 0.05      the line search never lets any gap fall
%                              below min_gap_frac * (smallest gap at theta0)
%   verbose        = true
%
% theta_new = theta - t_min * grad/|grad|,  t_min = argmin_{t in [0,t_max]} psi(theta - t*grad/|grad|)
% t_max is found by bisection so every gap stays above the floor (this is
% the paper's "Delta subset [0,infty) chosen so the gaps remain
% positive" -- delta is not further specified in the paper, so we enforce
% a strictly-positive floor rather than positivity in the limit, which
% would make t_max ill-defined/unstable numerically).

if nargin < 4, opts = struct(); end
if ~isfield(opts,'max_iters'), opts.max_iters = 100; end
if ~isfield(opts,'grad_tol'), opts.grad_tol = 1e-6; end
if ~isfield(opts,'rel_tol'), opts.rel_tol = 1e-8; end
if ~isfield(opts,'min_gap_frac'), opts.min_gap_frac = 0.05; end
if ~isfield(opts,'verbose'), opts.verbose = true; end

theta = theta0(:);
g0 = gaps_fn(theta);
gap_floor = opts.min_gap_frac * min(g0);

phi_hist = zeros(opts.max_iters+1,1);
[phi, grad] = costgrad(theta);
phi_hist(1) = phi;

info = struct('converged',false,'n_iters',0,'phi_hist',[]);

for it = 1:opts.max_iters
    gnorm = norm(grad);
    if gnorm < opts.grad_tol
        if opts.verbose, fprintf('  [Algorithm 1] converged: |grad|=%.3e < tol at iter %i\n',gnorm,it-1); end
        info.converged = true;
        break
    end
    d = grad/gnorm;   % descent direction step is theta - t*d

    t_max = find_t_max(gaps_fn, theta, d, gap_floor);
    if t_max <= 0
        if opts.verbose, fprintf('  [Algorithm 1] stalled: t_max<=0 at iter %i\n',it); end
        break
    end

    line_fun = @(t) costgrad(theta - t*d);
    t_star = fminbnd(line_fun, 0, t_max, optimset('TolX',1e-10*max(t_max,eps)));

    theta_new = theta - t_star*d;
    phi_new = costgrad(theta_new);

    if opts.verbose
        fprintf('  [Algorithm 1] iter %3i: psi=%.8e  |grad|=%.3e  t*=%.3e  t_max=%.3e\n',...
            it, phi_new, gnorm, t_star, t_max);
    end

    rel_dec = (phi - phi_new)/max(abs(phi),eps);
    theta = theta_new;
    [phi, grad] = costgrad(theta);
    phi_hist(it+1) = phi;

    if rel_dec < opts.rel_tol && rel_dec >= 0
        if opts.verbose, fprintf('  [Algorithm 1] converged: rel. decrease %.3e < tol at iter %i\n',rel_dec,it); end
        info.converged = true;
        info.n_iters = it;
        break
    end
    info.n_iters = it;
end

theta_opt = theta;
phi_opt = phi;
info.phi_hist = phi_hist(1:info.n_iters+1);
end

function t_max = find_t_max(gaps_fn, theta, d, gap_floor)
% Largest t>=0 such that min(gaps_fn(theta-t*d)) >= gap_floor, found by
% bracket-then-bisect (gaps_fn is cheap: no FEM solve).
t_hi = 0.01;
ok = @(t) min(gaps_fn(theta - t*d)) >= gap_floor;
if ~ok(0)
    t_max = 0; return
end
while ok(t_hi) && t_hi < 4*pi
    t_hi = t_hi*2;
end
if ok(t_hi)
    t_max = t_hi; return
end
lo = 0; hi = t_hi;
for iter = 1:60
    mid = 0.5*(lo+hi);
    if ok(mid), lo = mid; else, hi = mid; end
end
t_max = lo;
end
