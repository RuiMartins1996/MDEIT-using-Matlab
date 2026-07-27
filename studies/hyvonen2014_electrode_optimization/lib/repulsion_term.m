function [phi_rep, grad_rep] = repulsion_term(elec, alpha)
%REPULSION_TERM  psi_rep = alpha * sum_m 1/g_m  (paper eq (4.3), the term
%that keeps electrodes from colliding) and its exact gradient w.r.t. the
%M design variables theta_minus (PLAN_implementation.md sec 4.3).
%
%   [phi_rep, grad_rep] = repulsion_term(elec, alpha)
%
% elec = assemble_electrodes(...) output. g_m = elec(m).gap_next (arc
% length from electrode m's right end to electrode m+1's left end,
% cyclic). Derivative (derived in PLAN sec 4.3, re-derived and simplified
% here for our exact endpoint convention):
%   dg_k/dtheta_minus(k)     = -speed_plus(k)*dtp_dtm(k)
%   dg_{k-1}/dtheta_minus(k) = +speed_minus(k)
%   => dpsi_rep/dtheta_minus(k) = alpha*[ speed_plus(k)*dtp_dtm(k)/g(k)^2
%                                        - speed_minus(k)/g(k-1)^2 ]

M = numel(elec);
g = [elec.gap_next]';

phi_rep = alpha*sum(1./g);

grad_rep = zeros(M,1);
for k = 1:M
    kprev = mod(k-2,M)+1;
    grad_rep(k) = alpha*( elec(k).speed_plus*elec(k).dtp_dtm/g(k)^2 ...
                        - elec(k).speed_minus/g(kprev)^2 );
end
end
