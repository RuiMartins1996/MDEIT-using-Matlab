function dA = cem_dA_dthetaminus(ctx, elec, m)
%CEM_DA_DTHETAMINUS  dA_red/dtheta_minus(m) -- exact, via the fundamental
%theorem of calculus applied to the electrode arc integrals
%(PLAN_implementation.md sec 3.3, combined with (4.2) for the fixed-width
%chain rule dtheta_plus/dtheta_minus).
%
%   dA = cem_dA_dthetaminus(ctx, elec, m)   % sparse (n+M-1)x(n+M-1)
%
% Only electrode m's blocks are touched, and only at the (<=4) mesh nodes
% adjacent to its two arc endpoints -- this is a handful of nonzeros
% regardless of mesh size, exactly the "point evaluations at the
% electrode edges" structure the paper's continuum formula (2.20)/(4.1)
% predicts (see PLAN sec 3.3 for the derivation).

n = ctx.n_nodes; M = numel(elec);
e = elec(m);
zc = e.z;

dBm = e.dtp_dtm*e.speed_plus*(e.phi_plus*e.phi_plus.') ...
     - e.speed_minus*(e.phi_minus*e.phi_minus.');
dcm = e.dtp_dtm*e.speed_plus*e.phi_plus - e.speed_minus*e.phi_minus;
dElen = e.dtp_dtm*e.speed_plus - e.speed_minus;

dKB = dBm/zc;   % [n x n], contributes to the (u,u) block for ALL m (incl. grounded M)

if m < M
    dC = -dcm/zc;          % [n x 1], the (u,U_m) cross column
    dD = dElen/zc;          % scalar, the (U_m,U_m) diagonal entry
    Cblock = sparse(n, M-1);
    Cblock(:,m) = dC;
    Ddiag = zeros(M-1,1);
    Ddiag(m) = dD;
    D = spdiags(Ddiag,0,M-1,M-1);
else
    Cblock = sparse(n, M-1);
    D = sparse(M-1,M-1);
end

dA = [dKB, Cblock; Cblock.', D];
dA = (dA+dA.')/2;
end
