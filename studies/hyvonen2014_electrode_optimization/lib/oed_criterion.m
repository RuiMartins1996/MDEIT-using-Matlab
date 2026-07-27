function [phi, W, extra] = oed_criterion(J, Gamma_prior, Gamma_noise, opt_mode)
%OED_CRITERION  A-/D-optimality cost and gradient weight, via the
%data-space Woodbury form (PLAN_implementation.md sec 4.1). Never forms
%the n_elem x n_elem posterior covariance.
%
%   [phi, W, extra] = oed_criterion(J, Gamma_prior, Gamma_noise, opt_mode)
%
% J           : [n_data x n_elem] linearized measurement Jacobian
% Gamma_prior : [n_elem x n_elem] prior covariance (dense or sparse; for
%               Hyvonen et al.'s Gaussian smoothness prior (5.1) it is
%               generally DENSE within each kappa-block, unlike the
%               diagonal prior used in studies/optimal_sensors_bayesian_approach)
% Gamma_noise : [n_data x n_data] noise covariance (diagonal here)
% opt_mode    : 'a-opt' -> phi = trace(Gamma_*),      W = Gn^-1 J Gamma_*^2
%               'd-opt' -> phi = logdet(Gamma_*),      W = Gn^-1 J Gamma_*
%               (both minimized; see PLAN sec 4.1/sec 5 for the sign
%               derivation -- D-optimality here MINIMIZES logdet(Gamma_*),
%               which is exactly the paper's eq (3.12)/(4.3) criterion.)
%
% Derivation (Woodbury): with H = J'*Gn^-1*J + Gpr^-1, Gamma_* = H^-1,
%   P = J*Gpr,  S = Gn + P*J' = Gn + J*Gpr*J'  (SPD, n_data x n_data)
%   Y = S \ P                                   (= Gn^-1 J Gamma_* algebraically)
%   trace(Gamma_*) = trace(Gpr) - sum(sum(P.*Y))
%   logdet(Gamma_*) = logdet(Gpr) + logdet(Gn) - logdet(S)
%   dphi/dtheta = -2<W,dJ/dtheta>_F  with
%     W_D = Y
%     W_A = Yd - (Yd*J')*Y,   Yd := Y*Gpr

n_data = size(J,1);

P  = J*Gamma_prior;                          % [n_data x n_elem]
S  = Gamma_noise + P*J.';                    % [n_data x n_data]
S  = (S+S.')/2;
jitter = 1e-12*trace(S)/n_data;
Ls = chol(S + jitter*speye(n_data),'lower');
Y  = Ls.' \ (Ls \ P);                        % S^-1 P   [n_data x n_elem]

switch opt_mode
    case 'a-opt'
        phi = trace(Gamma_prior) - sum(sum(P.*Y));
        Yd  = Y*Gamma_prior;
        W   = Yd - (Yd*J.')*Y;
    case 'd-opt'
        logdet_S = 2*sum(log(diag(Ls)));
        n_elem = size(Gamma_prior,1);
        jitter_p = 1e-12*trace(Gamma_prior)/n_elem;
        Lp = chol(full(Gamma_prior) + jitter_p*eye(n_elem),'lower');
        logdet_Gpr = 2*sum(log(diag(Lp)));
        logdet_Gn  = sum(log(diag(Gamma_noise)));
        phi = logdet_Gpr + logdet_Gn - logdet_S;
        W   = Y;
    otherwise
        error('oed_criterion:unknownMode','opt_mode "%s" not recognized',opt_mode);
end

extra = struct('P',P,'S',S,'Ls',Ls,'Y',Y);
end
